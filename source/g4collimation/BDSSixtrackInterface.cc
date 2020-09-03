#include "BDSBunchSixTrackLink.hh"
#include "BDSException.hh"
#include "BDSIMLink.hh"
#include "BDSIonDefinition.hh"
#include "BDSParticleCoordsFull.hh"
#include "BDSParticleDefinition.hh"
#include "BDSPhysicsUtilities.hh"

#include "G4Electron.hh"
#include "G4GenericIon.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Types.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <G4GenericIon.hh>

std::string CleanFortranString(char* str, size_t count);

bool debugBDS = false;

BDSIMLink* bds = nullptr;
BDSBunchSixTrackLink* stp = nullptr;

extern "C"
void g4_collimation_init(double* referenceEk,
			 int*    seed,
			 double* relativeEnergyCut,
			 double* absoluteEnergyCut,
			 double* /*rcut*/,
			 double* /*rangecut_mm*/,
			 double* /*v0*/,
			 char*   /*PhysicsSelect*/,
			 bool*   /*g4_debug*/,
			 bool*   /*g4_keep_stable*/,
			 bool*   /*DoEnergyDeposition*/)
{
  stp = new BDSBunchSixTrackLink();
  bds = new BDSIMLink(stp);

  std::string seedStr = std::to_string(*seed);
  std::vector<std::string> arguments = {"--file=lhccrystals.gmad","--file=lhccrystals.gmad", "--vis_debug", "--batch",
                                        "--seed="+seedStr, "--outfile=output_"+seedStr};
  //std::vector<std::string> arguments = {"--file=lhccrystals.gmad","--batch"};
  std::vector<char*> argv;
  for (const auto& arg : arguments)
    {argv.push_back((char*)arg.data());}
  argv.push_back(nullptr);

  // absolute energy cut is in GeV from sixtrack
  double relEKCut = *relativeEnergyCut;
  if (relEKCut < 1e-6) // defaults to 0 which means 0eV cut which is bad
    {relEKCut = 1.0;}
  // referenceEk is in MeV, but absoluteEnergyCut is in GeV incoming here
  double minimumEK = std::min(relEKCut * (*referenceEk), (*absoluteEnergyCut)*1000);
  G4cout << "Minimum kinetic energy " << minimumEK << " MeV" << G4endl;
  try
    {bds->Initialise(argv.size() - 1, argv.data(), true, minimumEK / 1000.0, false);} // minimumEk in GeV
  catch (const std::exception& e)
    {std::cout << e.what() << std::endl; exit(1);}
}

extern "C"
void g4_add_collimator(char*   name,
		       char*   material,
		       double* lengthIn,
		       double* apertureIn,
		       double* /*rotationIn*/,
		       double* /*xOffsetIn*/,
		       double* /*yOffsetIn*/,
		       int*    side,
		       bool*   isACrystalIn,
		       double* crystalAngleIn)
{
  //  keep 48 value in sync with mNameLen in common_modules.f90
  std::string collimatorName = CleanFortranString(name, 48);
  std::string materialName   = CleanFortranString(material, 4);

  bool buildLeft  = *side == 0 || *side == 1;
  bool buildRight = *side == 0 || *side == 2;
  double length   = *lengthIn   * CLHEP::m;
  double aperture = *apertureIn * CLHEP::m;
  double crystalAngle = *crystalAngleIn * 1e-3;
  std::cout << "Aperture " << aperture << std::endl;
  std::cout << "Crystal angle " << crystalAngle << std::endl;
  bds->AddLinkCollimatorJaw(collimatorName,
			 materialName,
			 length,
			 0.5*aperture,
			 0.5*aperture,
			 0, // rotation done in sixtrack on particles so actually not needed
			 0, // x and y offsets are done on sixtrack side so actually not needed in the interface at all
			 0,
			 buildLeft,
			 buildRight,
			 *isACrystalIn,
			 crystalAngle);
}

extern "C"
void g4_terminate()
{
  //delete bds;
  bds = nullptr;
}

extern "C"
void g4_set_collimator(char* name)
{
  //  keep 48 value in sync with mNameLen in common_modules.f90
  std::string collimatorName = CleanFortranString(name, 48);
  bds->SelectLinkElement(collimatorName);
}

extern "C"
void g4_add_particle(double*  xIn,
		     double*  yIn,
		     double*  xpIn,
		     double*  ypIn,
		     double*  totalEnergyIn,
		     int32_t* pdgIDIn,
		     int16_t* nzz,
		     int16_t* naa,
		     int16_t* nqq,
		     double*  /*massIn*/,
		     double*  /*sigmv*/,
		     double*  /*sx*/,
		     double*  /*sy*/,
		     double*  /*sz*/)
{
  G4double totalEnergy = (*totalEnergyIn) * CLHEP::GeV;
  G4double xp          = (*xpIn);
  G4double yp          = (*ypIn);
  G4double zp = 1;
  try
    {zp = BDSBunch::CalculateZp(xp,yp,1);}
  catch (const BDSException& e)
    {
      G4cout << e.what() << "Particle not tracked" << G4endl;
      return;
    }
  BDSParticleCoordsFull coords = BDSParticleCoordsFull((*xIn) * CLHEP::m,
						       (*yIn) * CLHEP::m,
						       0,
						       xp,
						       yp,
						       zp,
						       0,
						       0,
						       totalEnergy,
						       1);

  long long int pdgID = (long long int)(*pdgIDIn);
  
  // Wrap in our class that calculates momentum and kinetic energy.
  // Requires that one of E, Ek, P be non-zero (only one).
  BDSParticleDefinition* particleDefinition = nullptr;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  try
    {
      // PDG for ions = 10LZZZAAAI
      if (pdgID > 1000000000) // is an ion
	{
	  G4GenericIon::GenericIonDefinition(); // construct general ion particle
	  G4int a = (G4int)(*naa);
	  G4int z = (G4int)(*nzz);
	  G4double q = (G4double)(*nqq) * CLHEP::eplus;
	  BDSIonDefinition* ionDef = new BDSIonDefinition(a,z,q);
	  
	  G4IonTable* ionTable = particleTable->GetIonTable();
	  G4double mass   = ionTable->GetIonMass(ionDef->Z(), ionDef->A());
	  mass += ionDef->NElectrons()*G4Electron::Definition()->GetPDGMass();
	  G4double charge = ionDef->Charge(); // correct even if overridden
	  particleDefinition = new BDSParticleDefinition("ion", mass, charge,
							 totalEnergy, 0, 0, 1, ionDef, pdgID);
	}
      else
	{
	  G4ParticleDefinition* particleDef = nullptr;
	  particleDef = particleTable->FindParticle(pdgID);
	  if (!particleDef)
	    {
	      G4cerr << "Particle \"" << std::to_string(pdgID) << "\" not found." << G4endl;
	      return;
	    }
	  particleDefinition = new BDSParticleDefinition(particleDef, totalEnergy, 0, 0, 1, nullptr);
	}
    }
  catch (const BDSException& e)
    {// if we throw an exception the object is invalid for the delete on the next loop
      particleDefinition = nullptr; // reset back to nullptr for safe delete
      return;
    }
  
  stp->AddParticle(particleDefinition, coords);
}

extern "C"
void g4_collimate()
{
  bds->ClearSamplerHits();
  bds->BeamOn((G4int)stp->Size());
}

extern "C"
void g4_collimate_return(int*     j,
			 double*  x,
			 double*  y,
			 double*  xp,
			 double*  yp,
			 double*  e,
			 int32_t* pdgid,
			 double*  m,
			 int16_t* z,
			 int16_t* a,
			 int16_t* q,
			 double*  sigmv,
			 int*     /*part_hit*/,
			 int*     /*part_abs*/,
			 double*  /*part_impact*/,
			 double*  /*part_indiv*/,
			 double*  /*part_linteract*/,
			 double*  sx,
			 double*  sy,
			 double*  sz)
{
  //if (*j == 0)
  //  {return;}
  /*
    part_hit(j), part_abs(j), part_impact(j), part_indiv(j),
    & part_linteract(j))
    !++  Output information:
    !++
    !++  PART_HIT(MAX_NPART)     Hit flag for last hit (10000*element# + turn#)
    !++  PART_ABS(MAX_NPART)     Abs flag (10000*element# + turn#)
    !++  PART_IMPACT(MAX_NPART)  Impact parameter (0 for inner face)
    !++  PART_INDIV(MAX_NPART)   Divergence of impacting particles
  */
  
  // here the units have been converted back to GeV and m (in the tracking action)
  const BDSHitsCollectionSamplerLink* hits = bds->SamplerHits();
  G4int ind = *j;
  if (ind > (G4int)hits->entries())
    {return;}
  if (debugBDS)
    {G4cout << "Returning particle with j and index " << *j << " " << ind << G4endl;}

  BDSHitSamplerLink* hit = (*hits)[ind];
  const BDSParticleCoordsFull& coords = hit->coords;
  *x  = coords.x / CLHEP::m;
  *y  = coords.y / CLHEP::m;
  // remember, sixtrack xp, yp are p_x / p_total
  *xp = coords.xp;
  *yp = coords.yp;
  *e  = coords.totalEnergy / CLHEP::GeV;
  *pdgid  = (int32_t)hit->pdgID;
  *z = (int16_t)hit->Z;
  *a = (int16_t)hit->A;
  *q = (int16_t)hit->charge;
  
  // nucm is in MeV on both sides
  *m = hit->mass;
  
  // spin
  *sx = 0.0;
  *sy = 0.0;
  *sz = 1.0;
  
  // time, must be converted for using with sigmv
  *sigmv  = coords.T / CLHEP::s;
}

std::string CleanFortranString(char* str, size_t count)
{
  // fortran string nightmares
  std::string whitespace(" \t\n\r");
  std::string sstring(str,count);
  
  size_t lpos = sstring.find_last_not_of(whitespace);
  if (lpos != std::string::npos)
    {sstring.erase(lpos+1);}

  std::transform(sstring.begin(), sstring.end(), sstring.begin(), ::toupper);
  // fortran string happy place
  return sstring;
}

extern "C"
void g4_get_particle_count(int* g4_npart)
{
  int count = 1;
  auto hits = bds->SamplerHits();
  if (hits)
    {count = hits->entries();}
  *g4_npart = count;
  if (debugBDS)
    {std::cout << "Returning " << count << " -> " << bds->NPrimariesToReturn() << " primaries and " << bds->NSecondariesToReturn() << " secondaries" << std::endl;}
}

extern "C"
void g4_collimation_clear()
{
  bds->ClearSamplerHits();
  stp->ClearParticles();
}

extern "C"
void g4_keep_id(int* /*id*/)
{
  //keep_ids.insert(*id);
}

extern "C"
void g4_add_edep(char* /*name_in*/,
		 int* /*type*/,
		 double* /*xstep*/,
		 double* /*ystep*/,
		 double* /*zstep*/,
		 double* /*xmax*/,
		 double* /*ymax*/,
		 double* /*zmax*/,
		 double* /*xbigstep*/,
		 double* /*ybigstep*/,
		 double* /*zbigstep*/,
		 double* /*xbigmax*/,
		 double* /*ybigmax*/,
		 double* /*zbigmax*/)
{;}

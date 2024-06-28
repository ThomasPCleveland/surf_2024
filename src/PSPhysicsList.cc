#include "PSPhysicsList.hh"
#include "G4CMPPhysics.hh"
#include "PSPhysics.hh"



#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"

#include "G4HadronInelasticQBBC.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4ChargeExchangePhysics.hh"
#include "G4IonPhysicsXS.hh"
#include "G4IonElasticPhysics.hh"

#include "G4OpticalPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
//#include "CDMSHadronicBiasing.hh"
// Constructor and destructor

//register necessary physics for simulation.
PSPhysicsList::PSPhysicsList(G4int verbose) : G4VModularPhysicsList() {
    SetVerboseLevel(verbose);
    if (verbose) G4cout << "G4CMPPhysicsList::constructor" << G4endl;

    defaultCutValue = DBL_MIN;	// 100*mm;

    RegisterPhysics(new G4CMPPhysics);		// Phonons and charge-carriers



    G4int ver = 0;

    // QBBC Physics
    RegisterPhysics(new G4EmStandardPhysics(ver));

// physics->SetCuts();

  // // Synchroton Radiation & GN Physics
  // physics->RegisterPhysics(new G4EmExtraPhysics(ver));

  // // Decays
  // physics->RegisterPhysics(new G4DecayPhysics(ver));

  // // Hadron Physics
  // physics->RegisterPhysics(new G4HadronElasticPhysicsXS(ver));

  // physics->RegisterPhysics(new G4StoppingPhysics(ver));

  // physics->RegisterPhysics(new G4IonPhysicsXS(ver));

  // physics->RegisterPhysics(new G4IonElasticPhysics(ver));

  // physics->RegisterPhysics(new G4HadronInelasticQBBC(ver));

  // // optical physics
  // physics->RegisterPhysics(new G4OpticalPhysics(ver));

  // // Neutron tracking cut
  // physics->RegisterPhysics(new G4NeutronTrackingCut(ver));

  // physics->RegisterPhysics(new G4ParallelWorldPhysics(ver));
  // physics->RegisterPhysics(new G4RadioactiveDecayPhysics(ver));
  //  physics->RegisterPhysics(new CDMSHadronicBiasing(ver));



    //// Synchroton Radiation & GN Physics
    //RegisterPhysics(new G4EmExtraPhysics(ver));

    //// Decays
    //RegisterPhysics(new G4DecayPhysics(ver));

    //// Hadron Physics
    //RegisterPhysics(new G4HadronElasticPhysicsXS(ver));

    //RegisterPhysics(new G4StoppingPhysics(ver));

    //RegisterPhysics(new G4IonPhysicsXS(ver));

    //RegisterPhysics(new G4IonElasticPhysics(ver));

    //RegisterPhysics(new G4HadronInelasticQBBC(ver));

    // optical physics
    RegisterPhysics(new G4OpticalPhysics(ver));

    //// Neutron tracking cut
    //RegisterPhysics(new G4NeutronTrackingCut(ver));

    ///*RegisterPhysics(new G4ParallelWorldPhysics(ver));
    //RegisterPhysics(new G4RadioactiveDecayPhysics(ver));*/
    ////RegisterPhysics(new CDMSHadronicBiasing(ver));



    ////Custom Physics class I made for g4cmpSecondaryProduction process.
    ////RegisterPhysics(new PSPhysics());
}

PSPhysicsList::~PSPhysicsList() { ; }

// These values are used as the default production thresholds
// for the world volume.

void PSPhysicsList::SetCuts() {
    SetCutsWithDefault();
}
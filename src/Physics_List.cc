#include "Physics_List.hh"
#include "Custom_Physics.hh"

#include "G4CMPPhysics.hh" // from G4CMP
// #include "G4CMPKaplanQP.hh"

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
// #include "CDMSHadronicBiasing.hh"

Physics_List::Physics_List(G4int verbosity) : G4VModularPhysicsList()
{
    SetVerboseLevel(verbosity);
    if (verbosity)
        G4cout << "G4CMPPhysicsList::constructor" << G4endl;
    defaultCutValue = DBL_MIN; // 100*mm;

    verbosity = 1;
    RegisterPhysics(new G4CMPPhysics);
    RegisterPhysics(new G4EmStandardPhysics(verbosity));
    // RegisterPhysics(new G4OpticalPhysics(verbosity));
    SetCuts();

    // RegisterPhysics(new G4EmExtraPhysics(verbosity));
    // RegisterPhysics(new G4DecayPhysics(verbosity));
    // RegisterPhysics(new G4HadronElasticPhysicsXS(verbosity));
    // RegisterPhysics(new G4StoppingPhysics(verbosity));
    // RegisterPhysics(new G4IonPhysicsXS(verbosity));
    // RegisterPhysics(new G4IonElasticPhysics(verbosity));
    // RegisterPhysics(new G4HadronInelasticQBBC(verbosity));
    // RegisterPhysics(new G4OpticalPhysics(verbosity));
    // RegisterPhysics(new G4NeutronTrackingCut(verbosity));
    // RegisterPhysics(new G4ParallelWorldPhysics(verbosity));
    // RegisterPhysics(new G4RadioactiveDecayPhysics(verbosity));
    // RegisterPhysics(new CDMSHadronicBiasing(verbosity));

    ////Custom Physics class I made for g4cmpSecondaryProduction process.
    ////RegisterPhysics(new PSPhysics());
}

Physics_List::~Physics_List() { ; }

void Physics_List::SetCuts()
{
    SetCutsWithDefault();
}
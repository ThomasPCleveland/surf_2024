#include "Custom_Physics.hh"

//#include "G4CMPPhysics.hh"
#include "G4CMPConfigManager.hh"
#include "G4CMPDriftBoundaryProcess.hh"
#include "G4CMPDriftElectron.hh"
#include "G4CMPDriftHole.hh"
#include "G4CMPDriftRecombinationProcess.hh"
#include "G4CMPDriftTrappingProcess.hh"
#include "G4CMPDriftTrapIonization.hh"
#include "G4CMPInterValleyScattering.hh"
#include "G4CMPLukeScattering.hh"
#include "G4CMPPhononBoundaryProcess.hh"
#include "G4CMPSecondaryProduction.hh"
#include "G4CMPTimeStepper.hh"
#include "G4CMPTrackLimiter.hh"
#include "G4GenericIon.hh"
#include "G4ParticleTable.hh"
#include "G4PhononDownconversion.hh"
#include "G4PhononLong.hh"
#include "G4PhononScattering.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4ProcessManager.hh"


#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"
#include "G4Neutron.hh"
#include "G4OpAbsorption.hh"
#include "G4StepLimiter.hh"

#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpWLS.hh"
#include "G4OpBoundaryProcess.hh"

//Note: this whole class is not used.

// Constructor sets global verbosity

Custom_Physics::Custom_Physics(const G4String& name)
    : G4VPhysicsConstructor(name) {
    SetVerboseLevel(G4CMPConfigManager::GetVerboseLevel());
}

// Create phonon and charage carrier particles for later use

void Custom_Physics::ConstructParticle() {
    //G4Gamma::Gamma();
    //G4Neutron::Definition();
    //G4OpticalPhoton::OpticalPhotonDefinition();
}

// Add physics processes to appropriate particles

void Custom_Physics::ConstructProcess() {
    // Only make processes once; will be deleted when physics list goes away
    //G4VProcess* maker = new G4CMPSecondaryProduction;
    //G4VProcess* stepLimiter = new G4StepLimiter;

    // Set process verbosity to match physics list, for diagnostics
    //if (1 > 0) {
        //maker->SetVerboseLevel(1);
    //}

    //G4ParticleDefinition* particle = 0;	// Reusable buffer for convenience
    //particle = G4Gamma::Gamma();
    //RegisterProcess(maker, particle);
    //RegisterProcess(stepLimiter, particle);

    //particle = G4Neutron::Neutron();
    //RegisterProcess(maker, particle);
    //RegisterProcess(stepLimiter, particle);

    //particle = G4OpticalPhoton::OpticalPhoton();
    //RegisterProcess(maker, particle);
    //RegisterProcess(stepLimiter, particle);
    /*G4VProcess* fAbsorptionProcess = new G4OpAbsorption();
    RegisterProcess(fAbsorptionProcess, particle);*/

    /*RegisterProcess(new G4OpAbsorption, particle);
    RegisterProcess(new G4OpRayleigh, particle);
    RegisterProcess(new G4OpMieHG, particle);
    RegisterProcess(new G4OpWLS, particle);
    RegisterProcess(new G4OpBoundaryProcess, particle);*/

    // Add processes only to locally known particles
}
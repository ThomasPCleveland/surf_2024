#include "PSPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4Geantino.hh"
#include "G4ParticleGun.hh"
#include "G4RandomDirection.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhononLong.hh"
#include "G4SystemOfUnits.hh"
#include "G4Neutron.hh"
#include "G4Gamma.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4OpticalPhoton.hh"

using namespace std;

PSPrimaryGeneratorAction::PSPrimaryGeneratorAction() {
  //All of these settings can be set in the command line interface either
  //through "ApplyCommand" in PhononSimulation.cc or by running the commands
  //in an interactive ui session.



  //number of particles to be created for each event.
  G4int n_particle = 1;
  //intantiate particle gun to create particles for each event
  fParticleGun  = new G4ParticleGun(n_particle);   

  //// default particle kinematics ("geantino" triggers random phonon choice)
  //fParticleGun->SetParticleDefinition(G4Geantino::Definition());
  //fParticleGun->SetParticlePosition(G4ThreeVector(0.0,0.0,0.0));
  //fParticleGun->SetParticleEnergy(15*eV);
  //fParticleGun->SetParticleMomentumDirection(G4RandomDirection());

  //what type of particle is it - neutron or photon?
  //G4ParticleDefinition* particle = G4Neutron::Neutron();
  //gamma in geant = photon (G4OpticalPhoton only has refractive and reflectiv properties of visible light
  // i think).
  G4ParticleDefinition* particle = G4Gamma::Gamma();

  //what polarization does the photon have?
  fParticleGun->SetParticlePolarization(G4ThreeVector(1., 0., 0.));
  //sets definition of particle
  fParticleGun->SetParticleDefinition(particle);

  //what momentum direction should the particle have?
  //fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.0, 0.0, 0.0));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, -1.0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));

  //where is the particle created?
  //fParticleGun->SetParticlePosition(G4ThreeVector(-17000.0 * um, 27000.0 * um, 0.0 * um));
  //fParticleGun->SetParticlePosition(G4ThreeVector(410.0 * um, 6750.0 * um, 0.0));
  //fParticleGun->SetParticlePosition(G4ThreeVector(410.0 * um, 6750.0 * um, -499.0 * um));
  fParticleGun->SetParticlePosition(G4ThreeVector());

  //what energy does it have?
  fParticleGun->SetParticleEnergy(2.1 * eV);
  //fParticleGun->SetParticleEnergy(10 * eV);
  //fParticleGun->SetParticleEnergy(20 * eV);
  //fParticleGun->SetParticleEnergy(100 * eV);
  //fParticleGun->SetParticleEnergy(1 * keV);
  //fParticleGun->SetParticleEnergy(1 * eV);
  //fParticleGun->SetParticleEnergy(1* GeV); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


PSPrimaryGeneratorAction::~PSPrimaryGeneratorAction() {
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//I think this function is run each time a new event happens, allows 
// you to have different particle properties for each event.
void PSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  if (fParticleGun->GetParticleDefinition() == G4Geantino::Definition()) {
    G4double selector = G4UniformRand();
    if (selector<0.53539) {
      fParticleGun->SetParticleDefinition(G4PhononTransSlow::Definition()); 
    } else if (selector<0.90217) {
      fParticleGun->SetParticleDefinition(G4PhononTransFast::Definition());
    } else {
      fParticleGun->SetParticleDefinition(G4PhononLong::Definition());
    }
  }

  //uncomment this to make primary particle travel in random direction.
  fParticleGun->SetParticleMomentumDirection(G4RandomDirection());
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



#include "Primary_Action.hh"

#include "G4Event.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4GeneralParticleSource.hh"
#include "G4RandomDirection.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4PhononLong.hh"
#include "G4SystemOfUnits.hh"

using namespace std;

Primary_Action::Primary_Action()
{
  source = new G4GeneralParticleSource();
  source->SetParticlePolarization(G4ThreeVector(1., 0., 0.));
  source->SetParticlePosition(G4ThreeVector(0, 1, 0));
}

Primary_Action::~Primary_Action()
{
  delete source;
}

void Primary_Action::GeneratePrimaries(G4Event *anEvent)
{
  source->GeneratePrimaryVertex(anEvent);
}
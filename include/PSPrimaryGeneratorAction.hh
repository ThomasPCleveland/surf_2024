#ifndef PSPrimaryGeneratorAction_h
#define PSPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

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


class G4ParticleGun;
class G4Event;

class PSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    PSPrimaryGeneratorAction();
    virtual ~PSPrimaryGeneratorAction();

public:
    virtual void GeneratePrimaries(G4Event*);

private:
    G4ParticleGun* fParticleGun;

};


#endif



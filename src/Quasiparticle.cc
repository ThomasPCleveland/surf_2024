#include "globals.hh"
#include "Quasiparticle.hh"
#include "G4CMPConfigManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <numeric>

Quasiparticle::Quasiparticle()
{
    energy = 0.0;
    time = 0;
    qpLifetime = 438 * ns;
    atGap = false;
    isZeroEnergy = false;
    id = 0;
}

Quasiparticle::Quasiparticle(G4double e)
{
    energy = e;
    time = 0;
    qpLifetime = 438 * ns;
    atGap = false;
    isZeroEnergy = false;
    id = 0;
}

Quasiparticle::Quasiparticle(G4double e, G4double t)
{
    energy = e;
    time = t;
    qpLifetime = 438 * ns;
    atGap = false;
    isZeroEnergy = false;
    id = 0;
}

Quasiparticle::Quasiparticle(G4double e, G4double t, G4double l)
{
    energy = e;
    time = t;
    qpLifetime = l;
    atGap = false;
    isZeroEnergy = false;
    id = 0;
}

G4double Quasiparticle::GetEnergy()
{
    return energy;
}

G4double Quasiparticle::GetTime()
{
    return time;
}

G4double Quasiparticle::GetLifeTime()
{
    return qpLifetime;
}

bool Quasiparticle::GetAtGap()
{
    return atGap;
}

bool Quasiparticle::GetIsZeroEnergy()
{
    return isZeroEnergy;
}

int Quasiparticle::GetId()
{
    return id;
}

void Quasiparticle::SetEnergy(G4double e)
{
    energy = e;
}

void Quasiparticle::SetTime(G4double t)
{
    time = t;
}

void Quasiparticle::SetLifeTime(G4double l)
{
    qpLifetime = l;
}

void Quasiparticle::SetAtGap(bool ag)
{
    atGap = ag;
}

void Quasiparticle::SetIsZeroEnergy(bool ze)
{
    isZeroEnergy = ze;
}

void Quasiparticle::SetId(int qpId)
{
    id = qpId;
}

//increase time by default QP lifetime
void Quasiparticle::IncrementTime()
{
    time += qpLifetime;
}

//increase QP's time by t.
void Quasiparticle::IncrementTime(G4double t)
{
    time += t;
}

Quasiparticle::~Quasiparticle()
{
    
}
#include "globals.hh"
#include "PSQP.hh"
#include "G4CMPConfigManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <numeric>

PSQP::PSQP()
{
    energy = 0.0;
    time = 0;
    qpLifetime = 438 * ns;
    atGap = false;
    isZeroEnergy = false;
    id = 0;
}

PSQP::PSQP(G4double e)
{
    energy = e;
    time = 0;
    qpLifetime = 438 * ns;
    atGap = false;
    isZeroEnergy = false;
    id = 0;
}

PSQP::PSQP(G4double e, G4double t)
{
    energy = e;
    time = t;
    qpLifetime = 438 * ns;
    atGap = false;
    isZeroEnergy = false;
    id = 0;
}

PSQP::PSQP(G4double e, G4double t, G4double l)
{
    energy = e;
    time = t;
    qpLifetime = l;
    atGap = false;
    isZeroEnergy = false;
    id = 0;
}

G4double PSQP::GetEnergy()
{
    return energy;
}

G4double PSQP::GetTime()
{
    return time;
}

G4double PSQP::GetLifeTime()
{
    return qpLifetime;
}

bool PSQP::GetAtGap()
{
    return atGap;
}

bool PSQP::GetIsZeroEnergy()
{
    return isZeroEnergy;
}

int PSQP::GetId()
{
    return id;
}

void PSQP::SetEnergy(G4double e)
{
    energy = e;
}

void PSQP::SetTime(G4double t)
{
    time = t;
}

void PSQP::SetLifeTime(G4double l)
{
    qpLifetime = l;
}

void PSQP::SetAtGap(bool ag)
{
    atGap = ag;
}

void PSQP::SetIsZeroEnergy(bool ze)
{
    isZeroEnergy = ze;
}

void PSQP::SetId(int qpId)
{
    id = qpId;
}

//increase time by default QP lifetime
void PSQP::IncrementTime()
{
    time += qpLifetime;
}

//increase QP's time by t.
void PSQP::IncrementTime(G4double t)
{
    time += t;
}

PSQP::~PSQP()
{
    
}
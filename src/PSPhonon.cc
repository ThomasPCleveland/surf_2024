#include "globals.hh"
#include "PSPhonon.hh"
#include "G4CMPConfigManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <numeric>
#include <math.h>

PSPhonon::PSPhonon(G4double e)
{
    energy = e;
    time = 0;
    status = true;
    RandomizeDirection();
}

PSPhonon::PSPhonon(G4double e, G4double t)
{
    energy = e;
    time = t;
    status = true;
    RandomizeDirection();
}

PSPhonon::PSPhonon(G4double e, G4double t, G4double z)
{
    energy = e;
    time = t;
    zmin = z;
    status = true;
    RandomizeDirection();
}

G4double PSPhonon::GetEnergy()
{
    return energy;
}

G4double PSPhonon::GetTime()
{
    return time;
}

bool PSPhonon::GetStatus()
{
    return status;
}

bool PSPhonon::GetStartTime() {
    return startTime;
}

G4ThreeVector PSPhonon::GetDirection()
{
    return direction;
}

G4double PSPhonon::GetZMin()
{
    return zmin;
}

void PSPhonon::SetEnergy(G4double e)
{
    energy = e;
}

void PSPhonon::SetTime(G4double t)
{
    time = t;
}

void PSPhonon::SetStatus(bool stat)
{
    status = stat;
}

void PSPhonon::SetStartTime(G4double t) {
    startTime = t;
}

void PSPhonon::SetDirection(G4ThreeVector dir)
{
    direction = dir;
}

void PSPhonon::SetZMin(G4double z)
{
    zmin = z;
}

//generates a random direction for phonon such that equal areas on the unit
//sphere have equal probability for containing direction of phonon chosen.

//uses cylindrical equal-area projection.
void PSPhonon::RandomizeDirection()
{
    double z;
    if (G4UniformRand() < 0.5)
    {
        //vector is up
        z = G4UniformRand() * (1 - zmin) + zmin;
    }
    else
    {
        //vector is down
        z = G4UniformRand() * (-zmin + 1) - zmin;
    }
    double theta = 2.0 * pi * G4UniformRand();
    double r = sqrt(1 - z * z);
    double x = r * cos(theta);
    double y = r * sin(theta);
    direction = G4ThreeVector(x, y, z);
}

//simply reflects z component of direction since the surface normal is assumed to be (0, 0, 1)
//might have to reprogram if a surface normal is different from this for one of the superconducting
//surfaces in the future.
void PSPhonon::ReflectDirection()
{
    double x = direction.x();
    double y = direction.y();
    double z = direction.z();
    direction = G4ThreeVector(x, y, -z);
}

PSPhonon::~PSPhonon()
{
    
}
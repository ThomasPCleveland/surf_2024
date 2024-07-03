#ifndef PSPhonon_hh
#define PSPhonon_hh 1

#include "G4Types.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include <vector>
#include <string>
#include <fstream>

class Phonon
{
public:
    Phonon(G4double e);
    Phonon(G4double e, G4double t);
    Phonon(G4double e, G4double t, G4double z);
    G4double GetEnergy();
    G4double GetTime();
    bool GetStatus();
    bool GetStartTime();
    G4ThreeVector GetDirection();
    G4double GetZMin();
    void SetEnergy(G4double e);
    void SetTime(G4double t);
    void SetStatus(bool stat);
    void SetStartTime(G4double t);
    void SetDirection(G4ThreeVector dir);
    void SetZMin(G4double z);
    void RandomizeDirection();
    void ReflectDirection();
    virtual ~Phonon();
private:
    G4double energy; //its energy
    G4double time; //its time
    bool status; //will it create QPs?
    G4double startTime; //birth time of phonon
    G4ThreeVector direction; //its direction
    G4double zmin; //min value for its z component. This is used to avoid numerical
                   // overflow error in Kaplan_QP.cc since a z value of 0 will cause
                   // a divide-by-zero error and small values can also result in an extremely
                   // high value for the distance the phonon travels, causing an overflow.
};

#endif

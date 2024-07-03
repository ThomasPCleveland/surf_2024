#ifndef PSQP_hh
#define PSQP_hh 1

#include "G4Types.hh"
#include "G4SystemOfUnits.hh"
#include <vector>
#include <string>
#include <fstream>

class Quasiparticle
{
public:
    Quasiparticle();
    Quasiparticle(G4double e);
    Quasiparticle(G4double e, G4double t);
    Quasiparticle(G4double e, G4double t, G4double l);
    G4double GetEnergy();
    G4double GetTime();
    G4double GetLifeTime();
    bool GetAtGap();
    bool GetIsZeroEnergy();
    int GetId();
    void SetEnergy(G4double e);
    void SetTime(G4double t);
    void SetLifeTime(G4double l);
    void SetAtGap(bool ag);
    void SetIsZeroEnergy(bool ze);
    void SetId(int qpId);
    void IncrementTime();
    void IncrementTime(G4double t);
    virtual ~Quasiparticle();
private:
    G4double energy; //its energy
    G4double time; //its time
    G4double qpLifetime; //its default lifetime (not used)
    bool atGap; //is its energy equal to delta?
    bool isZeroEnergy; //is its energy equal to 0
    int id; //its id.
};

#endif

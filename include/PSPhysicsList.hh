#ifndef PSPhysicsList_h
#define PSPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class PSPhysicsList : public G4VModularPhysicsList {
public:
    PSPhysicsList(G4int verbose = 0);
    virtual ~PSPhysicsList();

    // Use default copy/move semantics
    PSPhysicsList(const PSPhysicsList&) = default;
    PSPhysicsList(PSPhysicsList&&) = default;
    PSPhysicsList& operator=(const PSPhysicsList&) = default;
    PSPhysicsList& operator=(PSPhysicsList&&) = default;

public:
    void SetCuts();
};

#endif
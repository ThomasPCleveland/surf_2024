#ifndef PSPhysicsList_h
#define PSPhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class Physics_List : public G4VModularPhysicsList {
public:
    Physics_List(G4int verbose = 0);
    virtual ~Physics_List();

    // Use default copy/move semantics
    Physics_List(const Physics_List&) = default;
    Physics_List(Physics_List&&) = default;
    Physics_List& operator=(const Physics_List&) = default;
    Physics_List& operator=(Physics_List&&) = default;

public:
    void SetCuts();
};

#endif
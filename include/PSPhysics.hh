#ifndef PSPhysics_hh
#define PSPhysics_hh 1

#include "G4VPhysicsConstructor.hh"


class PSPhysics : public G4VPhysicsConstructor {
public:
    //PSPhysics(const G4String& name = "G4CMPPhysics");
    PSPhysics(const G4String& name = "PSPhysics");
    virtual ~PSPhysics() { ; }

public:
    virtual void ConstructParticle();	// Creates phonons and "drifters"
    virtual void ConstructProcess();	// Adds processes to physics list

protected:
    void AddSecondaryProduction();	// All charged particles make e/h, phn

private:
    PSPhysics(const PSPhysics& rhs);		// Copying is forbidden
    PSPhysics& operator=(const PSPhysics& rhs);
};

#endif	/* G4CMPPhysics_hh */
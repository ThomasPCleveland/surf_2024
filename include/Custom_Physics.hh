#ifndef PSPhysics_hh
#define PSPhysics_hh 1

#include "G4VPhysicsConstructor.hh"


class Custom_Physics : public G4VPhysicsConstructor {
public:
    //PSPhysics(const G4String& name = "G4CMPPhysics");
    Custom_Physics(const G4String& name = "PSPhysics");
    virtual ~Custom_Physics() { ; }

public:
    virtual void ConstructParticle();	// Creates phonons and "drifters"
    virtual void ConstructProcess();	// Adds processes to physics list

protected:
    void AddSecondaryProduction();	// All charged particles make e/h, phn

private:
    Custom_Physics(const Custom_Physics& rhs);		// Copying is forbidden
    Custom_Physics& operator=(const Custom_Physics& rhs);
};

#endif	/* G4CMPPhysics_hh */
#ifndef PSActionInitialization_hh
#define PSActionInitialization_hh 1

#include "G4VUserActionInitialization.hh"
#include "PSSteppingAction.hh"

class PSActionInitialization : public G4VUserActionInitialization {
public:
  PSActionInitialization() {;}
  virtual ~PSActionInitialization() {;}
  virtual void Build() const;
  PSSteppingAction* steppingAction;
};

#endif	/* PSActionInitialization_hh */
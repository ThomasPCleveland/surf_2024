#ifndef PSEventAction_h
#define PSEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "PSSteppingAction.hh"

class PSEventAction : public G4UserEventAction
{
public:
    PSEventAction(PSSteppingAction* stepAct);
    virtual ~PSEventAction();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

private:
    PSSteppingAction* steppingAction;
};
#endif
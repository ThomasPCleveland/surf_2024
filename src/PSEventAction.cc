#include "PSEventAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4VHit.hh"

PSEventAction::PSEventAction(PSSteppingAction* stepAct) : G4UserEventAction()
{
    steppingAction = stepAct;
}

PSEventAction::~PSEventAction()
{}

//code inside this function will execute at beginning of each event
void PSEventAction::BeginOfEventAction(const G4Event* /*anEvent*/)
{}

//code inside this function will execute at end of each event.
//tells steping action to record data and outputs a message saying
//the event has ended in the command line interface.
void PSEventAction::EndOfEventAction(const G4Event* anEvent)
{
    steppingAction->NewEvent();
    G4cout << "End of Event " << ((int)anEvent->GetEventID() + 1) << "\n";
}
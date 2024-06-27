#include "PSActionInitialization.hh"
#include "PSPrimaryGeneratorAction.hh"
#include "G4CMPStackingAction.hh"
#include "G4UserRunAction.hh"
#include "PSEventAction.hh"

void PSActionInitialization::Build() const {
  //tells geant to use PSPrimaryGeneratorAction for primary generator class
  SetUserAction(new PSPrimaryGeneratorAction);
  //not sure what this does, but its probably necessary
  SetUserAction(new G4CMPStackingAction);
  //create new event action, will run what ever code you put in its functions
  //at the beginning or end of each event
  PSEventAction* eventAction = new PSEventAction(steppingAction);
  //need to set it as user action, or else these methods will not be run
  //because geant does not know which event action class object to use
  SetUserAction(eventAction);
  //set steppingAction as stepping action class.
  SetUserAction(steppingAction);
} 

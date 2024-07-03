#include "Config_Messenger.hh"
#include "Config_Manager.hh"

#include "G4UIcmdWithAString.hh"

Config_Messenger::Config_Messenger(Config_Manager *mgr)
    : G4UImessenger("/g4cmp/", "User configuration for G4CMP phonon example"),
      theManager(mgr), hitsCmd(0)
{
  hitsCmd = CreateCommand<G4UIcmdWithAString>("HitsFile",
                                              "Set filename for output of phonon hit locations");
}

Config_Messenger::~Config_Messenger()
{
  delete hitsCmd;
  hitsCmd = 0;
}

// Parse user input and add to configuration
void Config_Messenger::SetNewValue(G4UIcommand *cmd, G4String value)
{
  if (cmd == hitsCmd)
    theManager->SetHitOutput(value);
}

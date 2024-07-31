#include "Config_Manager.hh"
#include "Config_Messenger.hh"
#include "G4RunManager.hh"
#include <stdlib.h>

Config_Manager *Config_Manager::singleton = 0;

Config_Manager *Config_Manager::Instance()
{
  if (!singleton) // if DNE
    singleton = new Config_Manager;
  return singleton;
}

Config_Manager::Config_Manager()
    : Hit_file(getenv("G4CMP_HIT_FILE") ? getenv("G4CMP_HIT_FILE") : "deposits.csv"),
      Primary_file("primaries.csv"),
      messenger(new Config_Messenger(this)) { ; }

Config_Manager::~Config_Manager()
{
  delete messenger;
  messenger = 0;
}

// Trigger rebuild of geometry if parameters change

void Config_Manager::UpdateGeometry()
{
  G4RunManager::GetRunManager()->ReinitializeGeometry(true);
}

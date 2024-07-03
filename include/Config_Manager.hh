/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef RISQTutorialConfigManager_hh
#define RISQTutorialConfigManager_hh 1

// $Id$
// File:  RISQTutorialConfigManager.hh
//
// Description:	Singleton container class for user configuration of G4CMP
//		phonon example. Looks for environment variables	at
//		initialization to set default values; active values may be
//		changed via macro commands (see Config_Messenger).
//
// 20170816  M. Kelsey -- Extract hit filename from G4CMPConfigManager.

#include "globals.hh"

class Config_Messenger;


class Config_Manager {
public:
  ~Config_Manager();	// Must be public for end-of-job cleanup
  static Config_Manager* Instance();   // Only needed by static accessors

  // Access current values
  static const G4String& GetHitOutput()  { return Instance()->Hit_file; }
  static const G4String& GetPrimaryOutput()  { return Instance()->Primary_file; }

  // Change values (e.g., via Messenger)
  static void SetHitOutput(const G4String& name)
    { Instance()->Hit_file=name; UpdateGeometry(); }

  // Change values (e.g., via Messenger)
  static void SetPrimaryOutput(const G4String& name)
    { Instance()->Hit_file=name; UpdateGeometry(); }

  
  static void UpdateGeometry();

private:
  Config_Manager();		// Singleton: only constructed on request
  Config_Manager(const Config_Manager&) = delete;
  Config_Manager(Config_Manager&&) = delete;
  Config_Manager& operator=(const Config_Manager&) = delete;
  Config_Manager& operator=(Config_Manager&&) = delete;

  static Config_Manager* singleton;

private:
  G4String Hit_file;	// Output file of e/h hits ($G4CMP_HIT_FILE)
  G4String Primary_file;	// Output file of primaries

  Config_Messenger* messenger;
};

#endif	/* RISQTutorialConfigManager_hh */

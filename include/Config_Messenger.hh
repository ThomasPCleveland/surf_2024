/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef RISQTutorialConfigMessenger_hh
#define RISQTutorialConfigMessenger_hh 1

// $Id$
// File:  Config_Messenger.hh
//
// Description:	Macro command defitions to set user configuration in
//		RISQTutorialConfigManager.
//
// 20170816  Michael Kelsey

#include "G4UImessenger.hh"

class Config_Manager;
class G4UIcmdWithAString;
class G4UIcommand;


class Config_Messenger : public G4UImessenger {
public:
  Config_Messenger(Config_Manager* theData);
  virtual ~Config_Messenger();

  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  Config_Manager* theManager;
  G4UIcmdWithAString* hitsCmd;

private:
  Config_Messenger(const Config_Messenger&);	// Copying is forbidden
  Config_Messenger& operator=(const Config_Messenger&);
};

#endif /* RISQTutorialConfigMessenger_hh */

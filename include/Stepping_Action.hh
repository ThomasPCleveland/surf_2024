/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef Stepping_Action
#define Stepping_Action_hh 1

#include "G4UserSteppingAction.hh"

#include <fstream>

class G4Step;

class Stepping_Action : public G4UserSteppingAction
{
public:

  Stepping_Action();
  virtual ~Stepping_Action();
  // virtual void UserSteppingAction(const G4Step* step);
  void ExportStepInformation( const G4Step * step );
  
private:

  //Step info output file
  std::ofstream fOutputFile;
  
  
  
  
};

#endif


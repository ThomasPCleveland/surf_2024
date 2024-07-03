/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: e58a61fedbb99b167e16dafebc9c8664ae0c7b94 $

#ifndef Action_Initialization_hh
#define Action_Initialization_hh 1

#include "G4VUserActionInitialization.hh"

class Action_Initialization : public G4VUserActionInitialization {
public:
  Action_Initialization() {;}
  virtual ~Action_Initialization() {;}
  virtual void Build() const;
};

#endif

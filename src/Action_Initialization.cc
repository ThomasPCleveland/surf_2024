/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

// $Id: 539f524339ae53ad098a07cfa3bebd07784d23dd $

#include "Action_Initialization.hh"
#include "Primary_Action.hh"
#include "Stepping_Action.hh"

#include "G4CMPStackingAction.hh"

void Action_Initialization::Build() const {
  SetUserAction(new Primary_Action);
  SetUserAction(new G4CMPStackingAction);
  SetUserAction(new Stepping_Action);
} 

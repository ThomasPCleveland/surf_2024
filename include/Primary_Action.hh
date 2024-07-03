/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/include/Primary_Action.hh
/// \brief Definition of the Primary_Action class
//
// $Id: ecbf57649dfaeb88e0fac25491bf8fb68c9308ec $
//

#ifndef RISQTutorialPrimaryGeneratorAction_h
#define RISQTutorialPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;

class Primary_Action : public G4VUserPrimaryGeneratorAction
{
public:
  Primary_Action();
  virtual ~Primary_Action();

public:
  virtual void GeneratePrimaries(G4Event *);

private:
  G4GeneralParticleSource *source;
};

#endif

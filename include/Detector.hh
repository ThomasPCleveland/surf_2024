/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/include/PhononDetectorConstruction.hh
/// \brief Definition of the RISQTutorialDetectorConstruction class
//
// $Id: 4c06153e9ea08f2a90b22c53e5c39bde4b847c07 $
//

#ifndef Detector_hh
#define Detector_hh 1

#include "Stepping_Action.hh"
#include "Electrode.hh"
#include "G4NistManager.hh"
#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4CMPVElectrodePattern.hh"
#include "G4LatticeManager.hh"
#include "G4CMPPhononElectrode.hh"

class G4Material;
class G4VPhysicalVolume;
class G4CMPSurfaceProperty;
class G4CMPElectrodeSensitivity;

class Detector : public G4VUserDetectorConstruction
{
public:
  Detector();
  virtual ~Detector();

public:
  virtual G4VPhysicalVolume *Construct();
  // Stepping_Action *steppingAction;

private:
  void DefineMaterials();
  void SetupGeometry();

private:
  G4NistManager *nist;
  G4Material *fSilicon;
  G4Material *fAluminum;
  G4Material *kidMat;
  G4Material *capMat;
  G4Material *feedMat;
  G4VPhysicalVolume *physWorld;
  G4CMPSurfaceProperty *indSurf;
  G4CMPSurfaceProperty *capSurf;
  G4CMPSurfaceProperty *feedSurf;
  G4CMPSurfaceProperty *worldSurf;
  G4CMPElectrodeSensitivity* fSuperconductorSensitivity;
  G4CMPPhononElectrode *kidElectrode;
  G4CMPPhononElectrode *capElectrode;
  G4CMPPhononElectrode *feedElectrode;
  G4bool fConstructed;
  G4bool fIfField;

public:
  inline void Field(G4bool bl) { fIfField = bl; }
};

#endif

/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/include/PhononDetectorConstruction.hh
/// \brief Definition of the RISQTutorialDetectorConstruction class
//
// $Id: 4c06153e9ea08f2a90b22c53e5c39bde4b847c07 $
//

#ifndef RISQTutorialDetectorConstruction_h
#define RISQTutorialDetectorConstruction_h 1


#include "PSSteppingAction.hh"
#include "PSElectrode.hh"
#include "G4NistManager.hh"
#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4CMPVElectrodePattern.hh"
#include "G4LatticeManager.hh"

class G4Material;
class G4VPhysicalVolume;
class G4CMPSurfaceProperty;
class G4CMPElectrodeSensitivity;

class RISQTutorialDetectorConstruction : public G4VUserDetectorConstruction {
public:
  RISQTutorialDetectorConstruction();
  virtual ~RISQTutorialDetectorConstruction();
  
public:
  virtual G4VPhysicalVolume* Construct();
  PSSteppingAction* steppingAction;
  
private:
  void DefineMaterials();
  void SetupGeometry();
  
private:
    G4NistManager* nist;
  G4Material* fSilicon;
  G4Material* fAluminum;
  G4Material* kidMat;
  G4Material* capMat;
  G4Material* feedMat;
  G4VPhysicalVolume* physWorld;
  G4CMPSurfaceProperty* indSurf;
  G4CMPSurfaceProperty* capSurf;
  G4CMPSurfaceProperty* feedSurf;
  G4CMPSurfaceProperty* worldSurf;
  PSElectrode* kidElectrode;
  PSElectrode* capElectrode;
  PSElectrode* feedElectrode;
  G4bool fConstructed;
  G4bool fIfField;
  
public:
  inline void Field(G4bool bl) { fIfField = bl; }
};


#endif


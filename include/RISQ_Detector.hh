/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/include/PhononDetectorConstruction.hh
/// \brief Definition of the RISQ_Detector class
//
// $Id: 4c06153e9ea08f2a90b22c53e5c39bde4b847c07 $
//

#ifndef RISQ_Detector_h
#define RISQ_Detector_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Material;
class G4VPhysicalVolume;
class G4CMPSurfaceProperty;
class G4CMPElectrodeSensitivity;

class RISQ_Detector : public G4VUserDetectorConstruction {
public:
  RISQ_Detector();
  virtual ~RISQ_Detector();
  
public:
  virtual G4VPhysicalVolume* Construct();
  
private:
  void DefineMaterials();
  void SetupGeometry();
  void AttachPhononSensor(G4CMPSurfaceProperty * surfProp);

  
private:
  G4Material *fSilicon;
  G4Material *fAluminum;
  G4Material *kidMat;
  G4Material *capMat;
  G4Material *feedMat;
  G4VPhysicalVolume *fWorldPhys;
  G4CMPSurfaceProperty *indSurf;
  G4CMPSurfaceProperty *capSurf;
  G4CMPSurfaceProperty *feedSurf;
  G4CMPSurfaceProperty *worldSurf;
  G4CMPElectrodeSensitivity* fSuperconductorSensitivity;
  G4bool fConstructed;
  G4bool fIfField;
  
public:
  inline void Field(G4bool bl) { fIfField = bl; }
};


#endif


/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

/// \file exoticphysics/phonon/src/PhononDetectorConstruction.cc \brief
/// Implementation of the PhononDetectorConstruction class
//
// $Id: a2016d29cc7d1e75482bfc623a533d20b60390da $
//
// 20140321  Drop passing placement transform to G4LatticePhysical
// 20211207  Replace G4Logical*Surface with G4CMP-specific versions.
// 20220809  [ For M. Hui ] -- Add frequency dependent surface properties.

#include "Single_KID.hh"
#include "RISQ_Sensitivity.hh"
// #include "RISQTutorialQubitHousing.hh"
// #include "RISQTutorialPad.hh"
// #include "RISQTutorialTransmissionLine.hh"
// #include "RISQTutorialStraightFluxLine.hh"
// #include "RISQTutorialCornerFluxLine.hh"
// #include "RISQTutorialResonatorAssembly.hh"
#include "G4CMPPhononElectrode.hh"
#include "G4CMPElectrodeSensitivity.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4FieldManager.hh"
#include "G4GeometryManager.hh"
#include "G4LatticeLogical.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPLogicalBorderSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4Sphere.hh"
#include "G4SystemOfUnits.hh"
#include "G4TransportationManager.hh"
#include "G4Tubs.hh"
#include "G4UniformMagField.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

Single_KID::Single_KID()
    : fLiquidHelium(0), fGermanium(0), fAluminum(0), fTungsten(0), fSilicon(0), fWorldPhys(0),
      fSuperconductorSensitivity(0), fNiobium(0), indSurfProp(0), capSurfProp(0), feedSurfProp(0), wallSurfProp(0), fConstructed(false) { ; } //, fIfField(true) {;}

Single_KID::~Single_KID()
{
  delete indSurfProp;
  delete capSurfProp;
  delete feedSurfProp;
  delete wallSurfProp;
}

G4VPhysicalVolume *Single_KID::Construct()
{
  if (fConstructed)
  {
    if (!G4RunManager::IfGeometryHasBeenDestroyed())
    {
      // Run manager hasn't cleaned volume stores. This code shouldn't execute
      G4GeometryManager::GetInstance()->OpenGeometry();
      G4PhysicalVolumeStore::GetInstance()->Clean();
      G4LogicalVolumeStore::GetInstance()->Clean();
      G4SolidStore::GetInstance()->Clean();
    }
    // Have to completely remove all lattices to avoid warning on reconstruction
    G4LatticeManager::GetLatticeManager()->Reset();
    // Clear all LogicalSurfaces
    // NOTE: No need to redefine the G4CMPSurfaceProperties
    G4CMPLogicalBorderSurface::CleanSurfaceTable();
  }

  DefineMaterials();
  SetupGeometry();
  fConstructed = true;

  return fWorldPhys;
}

void Single_KID::DefineMaterials()
{
  G4NistManager *nistManager = G4NistManager::Instance();

  fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected
  fGermanium = nistManager->FindOrBuildMaterial("G4_Ge");
  fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  fTungsten = nistManager->FindOrBuildMaterial("G4_W");
  fNiobium = nistManager->FindOrBuildMaterial("G4_Nb");
}

void Single_KID::SetupGeometry()
{
  G4VSolid *solid_world = new G4Box("World", 16. * cm, 16. * cm, 16. * cm);
  G4LogicalVolume *log_world = new G4LogicalVolume(solid_world, fLiquidHelium, "World");
  // worldLogical->SetUserLimits(new G4UserLimits(10*mm, DBL_MAX, DBL_MAX, 0, 0));
  log_world->SetVisAttributes(G4VisAttributes::Invisible);
  fWorldPhys = new G4PVPlacement(0,
                                 G4ThreeVector(),
                                 log_world,
                                 "World",
                                 0,
                                 false,
                                 0);

  bool checkOverlaps = true;

  // create geometry for outer box
  G4double siliconThickness = 0.05 * cm / 2.0;
  G4double siliconWidth = 2.0 * cm / 2.0;

  // old name is QubitChip_solid for future reference

  G4Box *solid_siliconChip = new G4Box("fSiliconSolid",
                                       siliconWidth,
                                       siliconWidth,
                                       siliconThickness);

  // Now attribute a physical material to the chip
  G4LogicalVolume *log_siliconChip = new G4LogicalVolume(solid_siliconChip,
                                                         fSilicon,
                                                         "SiliconChip_log");

  // Now, create a physical volume and G4PVPlacement for storing as the final output
  G4ThreeVector siliconChipTranslate(0, 0, 0);
  G4VPhysicalVolume *phys_siliconChip = new G4PVPlacement(0,
                                                          siliconChipTranslate,
                                                          log_siliconChip,
                                                          "SiliconChip",
                                                          log_world,
                                                          false,
                                                          0,
                                                          checkOverlaps);

  G4VisAttributes *siliconChipVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  siliconChipVisAtt->SetVisibility(true);
  log_siliconChip->SetVisAttributes(siliconChipVisAtt);

  // Set up the G4CMP silicon lattice information using the G4LatticeManager
  //  G4LatticeManager gives physics processes access to lattices by volume
  G4LatticeManager *LM = G4LatticeManager::GetLatticeManager();
  G4LatticeLogical *log_siliconLattice = LM->LoadLattice(fSilicon, "Si");

  // G4LatticePhysical assigns G4LatticeLogical a physical orientation
  G4LatticePhysical *phys_siliconLattice = new G4LatticePhysical(log_siliconLattice);
  phys_siliconLattice->SetMillerOrientation(1, 0, 0);
  LM->RegisterLattice(phys_siliconChip, phys_siliconLattice);

  // Set up border surfaces
  G4CMPLogicalBorderSurface *border_siliconChip_world = new G4CMPLogicalBorderSurface("border_siliconChip_world", phys_siliconChip, fWorldPhys, fSiVacuumInterface);

  // KIDs and feed
  // Add capacitor and inductor for KID. They are made of aluminum.
  float indLength = 0.72 * mm / 2.0f;
  float capLength = 0.43 * mm / 2.0f;
  float kidWidth = 1.13 * mm / 2.0f;
  float kidThickness = 0.03 * um / 2.0f;
  float feedThickness = 0.07 * um;
  float feedWidth = 0.2 * mm / 2.0f;
  float feedIndGap = 30 * um;
  float capIndSpacing = 80 * um;

  // Inductor.
  G4VSolid *fIndSolid = new G4Box("Ind", kidWidth, indLength, kidThickness);
  G4LogicalVolume *fIndLogical = new G4LogicalVolume(fIndSolid, fAluminum, "fIndLogical");
  G4VPhysicalVolume *IndPhys = new G4PVPlacement(0, G4ThreeVector(0, indLength + feedWidth / 2.0 + feedIndGap, siliconThickness + kidThickness),
                                                 fIndLogical, "fIndLogical", log_world, false, 0, checkOverlaps);

  // Capacitor.
  // capLength + 2.0 * indLength + capIndSpacing
  G4VSolid *fCapSolid = new G4Box("Cap", kidWidth, capLength, kidThickness);
  G4LogicalVolume *fCapLogical = new G4LogicalVolume(fCapSolid, fAluminum, "fCapLogical");
  G4VPhysicalVolume *CapPhys = new G4PVPlacement(0, G4ThreeVector(0, capLength + feedWidth / 2.0 + feedIndGap + (2 * indLength) + capIndSpacing, siliconThickness + kidThickness),
                                                 fCapLogical, "fCapLogical", log_world, false, 0, checkOverlaps);

  // Feedline.
  G4VSolid *fFeedSolid = new G4Box("Feed", siliconWidth, feedWidth, feedThickness);
  G4LogicalVolume *fFeedLogical = new G4LogicalVolume(fFeedSolid, fNiobium, "fFeedLogical");
  G4VPhysicalVolume *FeedPhys = new G4PVPlacement(0, G4ThreeVector(0, -feedWidth / 2.0, siliconThickness + feedThickness),
                                                  fFeedLogical, "fFeedLogical", log_world, false, 0, checkOverlaps);

  // Now we establish a sensitivity object
  G4SDManager *SDman = G4SDManager::GetSDMpointer();
  if (!fSuperconductorSensitivity)
    fSuperconductorSensitivity = new RISQ_Sensitivity("PhononElectrode");
  SDman->AddNewDetector(fSuperconductorSensitivity);
  log_siliconChip->SetSensitiveDetector(fSuperconductorSensitivity);

  // First, define border surface properties that can be referenced later
  const G4double GHz = 1e9 * hertz;

  // the following coefficients and cutoff values are not well-motivated
  // the code below is used only to demonstrate how to set these values.
  const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 0};  // Turn this off temporarily
  const std::vector<G4double> diffCoeffs = {1, 0, 0, 0, 0, 0}; // Explicitly make this 1 for now
  const std::vector<G4double> specCoeffs = {0, 0, 0, 0, 0, 0}; // Turn this off temporarily
  const G4double anhCutoff = 520., reflCutoff = 350.;          // Units external

  // These are just the definitions of the interface TYPES, not the interfaces themselves. These must be called in a set of loops
  // below, and invoke these surface definitions.
  if (!fConstructed)
  {
    const G4double GHz = 1e9 * hertz;

    // Copied from Erik
    float pAbsProb = 0.1;
    // drop as well (0.01)
    indSurfProp = new G4CMPSurfaceProperty("IndSurf", 0.0, 1.0, 0.0, 0.0, pAbsProb, 1.0, 1.0, 0.0);
    capSurfProp = new G4CMPSurfaceProperty("CapSurf", 0.0, 1.0, 0.0, 0.0, pAbsProb, 1.0, 1.0, 0.0);
    feedSurfProp = new G4CMPSurfaceProperty("FeedSurf", 0.0, 1.0, 0.0, 0.0, pAbsProb, 1.0, 1.0, 0.0);

    // Add sensors to inductor, capacitor, feedline
    AttachPhononSensor(indSurfProp);
    AttachPhononSensor(capSurfProp);
    AttachPhononSensor(feedSurfProp);

    /*
    // Wall properties: Copied from example (says it's not "well motivated")
    //the following coefficients and cutoff values are not well-motivated
    //the code below is used only to demonstrate how to set these values.
    const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
    const std::vector<G4double> diffCoeffs =
      {5.88e-2, 7.83e-4, -2.47e-6, 1.71e-8, -2.98e-11};
    const std::vector<G4double> specCoeffs =
      {0,928, -2.03e-4, -3.21e-6, 3.1e-9, 2.9e-13};
    const G4double anhCutoff = 0., reflCutoff = 350.;   // Units external
    */

    //.001 probability (prob much lower) 1/10 or 1/100 of metal prob

    wallSurfProp = new G4CMPSurfaceProperty("WallSurf", 0.0, 1.0, 0.0, 0.0, 0.0001, 1.0, 1.0, 0.0);
    // wallSurfProp->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs, diffCoeffs, specCoeffs, GHz, GHz,GHz);
  }

  //
  // Separate surfaces for sensors vs. bare sidewall
  //
  new G4CMPLogicalBorderSurface("detInd", phys_siliconChip, IndPhys, indSurfProp);
  new G4CMPLogicalBorderSurface("detCap", phys_siliconChip, CapPhys, capSurfProp);
  new G4CMPLogicalBorderSurface("detFeed", phys_siliconChip, FeedPhys, feedSurfProp);
  new G4CMPLogicalBorderSurface("detWall", phys_siliconChip, fWorldPhys, wallSurfProp);
}

// Set up a phonon sensor for this surface property object. I'm pretty sure that this
// phonon sensor doesn't get stapled to individual geometrical objects, but rather gets
// stapled to a surface property, but I'm not sure... have to ask mKelsey
void Single_KID::AttachPhononSensor(G4CMPSurfaceProperty *surfProp)
{

  float kidThickness = 0.03 * um;

  // If no surface, don't do anything
  if (!surfProp)
    return;
  // Properties must be added to existing surface-property table
  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();

  // Check what kind of object we're putting a sensor on
  G4String name = surfProp->GetName();
  float filmAbsProb = 0.5; // From phonon example code

  if (name == "IndSurf")
  {
    sensorProp->AddConstProperty("gapEnergy", 0.000361 / 2.0 * eV); // Aluminum 0.00036/2
    sensorProp->AddConstProperty("phononLifetime", 0.242 * ns);
    sensorProp->AddConstProperty("phononLifetimeSlope", 0.29); // why did Erik set it to pi??
    sensorProp->AddConstProperty("vSound", 3260 * m / s);
    sensorProp->AddConstProperty("vSubSound", 8860 * m / s);
    sensorProp->AddConstProperty("filmThickness", kidThickness);
    sensorProp->AddConstProperty("lowQPLimit", 3);
    sensorProp->AddConstProperty("subgapAbsorption", 0);
    sensorProp->AddConstProperty("filmAbsorption", filmAbsProb);
    sensorProp->AddConstProperty("qpLifetime", 2350); // matches close to true qp lifetime
  }
  else if (name == "CapSurf")
  {
    sensorProp->AddConstProperty("gapEnergy", 0.000361 / 2.0 * eV); // Aluminum
    // sensorProp->AddConstProperty("phononLifetime", 0.00417 * ns);
    sensorProp->AddConstProperty("phononLifetime", 0.242 * ns);
    sensorProp->AddConstProperty("phononLifetimeSlope", 0.29);
    sensorProp->AddConstProperty("vSound", 3260 * m / s);
    sensorProp->AddConstProperty("vSubSound", 8860 * m / s);
    sensorProp->AddConstProperty("filmThickness", kidThickness);
    sensorProp->AddConstProperty("lowQPLimit", 3);
    sensorProp->AddConstProperty("subgapAbsorption", 0);
    sensorProp->AddConstProperty("filmAbsorption", filmAbsProb);
    sensorProp->AddConstProperty("qpLifetime", 2350);
  }
  else if (name == "FeedSurf")
  {
    sensorProp->AddConstProperty("gapEnergy", 0.00279 / 2.0 * eV); // original 0.00279 / 2.0
    sensorProp->AddConstProperty("phononLifetime", 0.242 * ns);
    sensorProp->AddConstProperty("phononLifetimeSlope", 3.14159);
    sensorProp->AddConstProperty("vSound", 3480 * m / s);
    sensorProp->AddConstProperty("vSubSound", 8860 * m / s);
    sensorProp->AddConstProperty("filmThickness", kidThickness);
    sensorProp->AddConstProperty("lowQPLimit", 3);
    sensorProp->AddConstProperty("subgapAbsorption", 0);
    sensorProp->AddConstProperty("filmAbsorption", filmAbsProb);
    sensorProp->AddConstProperty("qpLifetime", 0.149);
  }

  surfProp->SetPhononElectrode(new G4CMPPhononElectrode);
}
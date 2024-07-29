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

#include "RISQ_Detector.hh"
#include "RISQ_Sensitivity.hh"

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
// RISQ_Detector::RISQ_Detector()
//     : fLiquidHelium(0), fGermanium(0), fAluminum(0), fTungsten(0),
//       fWorldPhys(0),
//       fSuperconductorSensitivity(0), fConstructed(false) { ; } //, fIfField(true) {;}

RISQ_Detector::RISQ_Detector()
    : /*fLiquidHelium(0), fGermanium(0), fAluminum(0), fTungsten(0),
      fWorldPhys(0), topSurfProp(0), botSurfProp(0), wallSurfProp(0),
      electrodeSensitivity(0), fConstructed(false), fIfField(true)*/
      fSilicon(0), fAluminum(0), fSuperconductorSensitivity(0), fConstructed(false), fIfField(true)
{
}

RISQ_Detector::~RISQ_Detector() { ; }

G4VPhysicalVolume *RISQ_Detector::Construct()
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

void RISQ_Detector::DefineMaterials()
{
  G4NistManager *nistManager = G4NistManager::Instance();

  // fLiquidHelium = nistManager->FindOrBuildMaterial("G4_AIR"); // to be corrected
  // fGermanium = nistManager->FindOrBuildMaterial("G4_Ge");
  // fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
  // fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  // fTungsten = nistManager->FindOrBuildMaterial("G4_W");
  // fNiobium = nistManager->FindOrBuildMaterial("G4_Nb");

  fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
  fAluminum = nistManager->FindOrBuildMaterial("G4_Al");
  kidMat = nistManager->FindOrBuildMaterial("G4_Al");
  feedMat = nistManager->FindOrBuildMaterial("G4_Nb");
}

void RISQ_Detector::SetupGeometry()
{
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008 * g / mole;
  G4double density = 1.e-25 * g / cm3;
  G4double temperature = 2.73 * kelvin;
  G4double pressure = 3.e-18 * pascal;
  G4Material *vacuum =
      new G4Material("interGalactic", atomicNumber,
                     massOfMole, density, kStateGas,
                     temperature, pressure);

  G4bool checkOverlaps = true;

  // world
  G4double worldSizeXZ = 12 * cm;
  G4double worldSizeY = 10 * cm;
  G4Material *worldMat = vacuum;
  G4String worldName = "World";

  G4Box *solidWorld =
      new G4Box(worldName,                                               // its name
                0.5 * worldSizeXZ, 0.5 * worldSizeY, 0.5 * worldSizeXZ); // its size

  G4LogicalVolume *logicWorld =
      new G4LogicalVolume(solidWorld, // its solid
                          worldMat,   // its material
                          worldName); // its name

  fWorldPhys =
      new G4PVPlacement(0,               // no rotation
                        G4ThreeVector(), // at (0,0,0)
                        logicWorld,      // its logical volume
                        worldName,       // its name
                        0,               // its mother  volume
                        false,           // no boolean operation
                        0,               // copy number
                        checkOverlaps);  // overlaps checking

  double sThick = 1 * mm / 2.0;
  // double sThick = 1 * mm;
  //  Silicon cylinder - this is the volume in which we will propagate phonons
  G4VSolid *solidTarget = new G4Tubs("TargetS", 0. * cm, 38.5 * mm,
                                     sThick, 0. * deg, 360. * deg);
  G4LogicalVolume *logicTarget =
      new G4LogicalVolume(solidTarget, fSilicon, "TargetL");
  /*G4VPhysicalVolume* SiPhys =
      new G4PVPlacement(new G4RotationMatrix(180 * deg, 90.0 * deg, 0.0 * deg), G4ThreeVector(), logicTarget, "TargetP",
          logicWorld, false, 0, checkOverlaps);*/
  G4VPhysicalVolume *SiPhys =
      new G4PVPlacement(0, G4ThreeVector(), logicTarget, "fSiliconPhysical",
                        logicWorld, true, 0);

  // double sThick = 1 * mm / 2.0;
  ////double sThick = 1 * mm;
  //// Silicon cylinder - this is the volume in which we will propagate phonons
  // G4VSolid* solidTarget = new G4Box("TargetS", 11 * mm, 11 * mm,
  //     sThick);
  // G4LogicalVolume* logicTarget =
  //     new G4LogicalVolume(solidTarget, fSilicon, "TargetL");
  ///*G4VPhysicalVolume* SiPhys =
  //    new G4PVPlacement(new G4RotationMatrix(180 * deg, 90.0 * deg, 0.0 * deg), G4ThreeVector(), logicTarget, "TargetP",
  //        logicWorld, false, 0, checkOverlaps);*/
  // G4VPhysicalVolume* SiPhys =
  //    new G4PVPlacement(0, G4ThreeVector(), logicTarget, "fGermaniumPhysical",
  //        logicWorld, true, 0);

  ////Hit Location
  // G4Colour hitCol(0.0, 1.0, 0.0);
  // G4VisAttributes* hitAtt = new G4VisAttributes(hitCol);
  // G4VSolid* solidHit = new G4Box("HitS", 0.25 * mm, 0.25 * mm,
  //     0.01 * mm);
  // G4LogicalVolume* logicHit =
  //     new G4LogicalVolume(solidHit, fSilicon, "HitL");
  // logicHit->SetVisAttributes(hitAtt);
  // G4VPhysicalVolume* HitPhys =
  //     new G4PVPlacement(0, G4ThreeVector(0.0 * um, 0.0 * um, 2000 * um), logicHit, "HitPhysical",
  //         logicWorld, true, 0);

  // Silicon lattice information
  //  G4LatticeManager gives physics processes access to lattices by volume
  G4LatticeManager *LM = G4LatticeManager::GetLatticeManager();
  G4LatticeLogical *SiLogical = LM->LoadLattice(fSilicon, "Si");
  // G4LatticePhysical assigns G4LatticeLogical a physical orientation
  G4LatticePhysical *SiPhysical = new G4LatticePhysical(SiLogical);
  SiPhysical->SetMillerOrientation(1, 0, 0);
  // SiPhysical->SetMillerOrientation(0, 1, 0);
  // SiPhysical->SetMillerOrientation(0, 0, 1);
  // SiPhysical->SetMillerOrientation(1, 1, 0);
  // SiPhysical->SetMillerOrientation(1, 0, 1);
  // SiPhysical->SetMillerOrientation(0, 1, 1);
  // SiPhysical->SetMillerOrientation(1, 1, 1);
  LM->RegisterLattice(SiPhys, SiPhysical);

  // KIDs
  const double kidPos[] = {
      -10617,
      -28176,
      -3438,
      -28176,
      3740,
      -28176,
      10919,
      -28176,
      -10617,
      -22654,
      -3438,
      -22654,
      3740,
      -22654,
      10919,
      -22654,
      -24975,
      -17131,
      -17796,
      -17131,
      -10617,
      -17131,
      -3438,
      -17131,
      3740,
      -17131,
      10919,
      -17131,
      18098,
      -17131,
      25277,
      -17131,
      -32154,
      -11609,
      -24975,
      -11609,
      -17796,
      -11609,
      -10617,
      -11609,
      -3438,
      -11609,
      3740,
      -11609,
      10919,
      -11609,
      18098,
      -11609,
      25277,
      -11609,
      32456,
      -11609,
      -32154,
      -6087,
      -24975,
      -6087,
      -17796,
      -6087,
      -10617,
      -6087,
      -3438,
      -6087,
      3740,
      -6087,
      10919,
      -6087,
      18098,
      -6087,
      25277,
      -6087,
      32456,
      -6087,
      -32154,
      -565,
      -24975,
      -565,
      -17796,
      -565,
      -10617,
      -565,
      -3438,
      -565,
      3740,
      -565,
      10919,
      -565,
      18098,
      -565,
      25277,
      -565,
      32456,
      -565,
      -32154,
      4957,
      -24975,
      4957,
      -17796,
      4957,
      -10617,
      4957,
      -3438,
      4957,
      3740,
      4957,
      10919,
      4957,
      18098,
      4957,
      25277,
      4957,
      32456,
      4957,
      -17796,
      10479,
      -10617,
      10479,
      -3438,
      10479,
      3740,
      10479,
      10919,
      10479,
      18098,
      10479,
      -17796,
      16001,
      -10617,
      16001,
      -3438,
      16001,
      3740,
      16001,
      10919,
      16001,
      18098,
      16001,
      -17796,
      21524,
      -10617,
      21524,
      -3438,
      21524,
      3740,
      21524,
      10919,
      21524,
      18098,
      21524,
      -17796,
      27046,
      -10617,
      27046,
      -3438,
      27046,
      3740,
      27046,
      10919,
      27046,
      18098,
      27046,
  };

  const double kidLength[] = {1256.409, 1252.532, 1175.089, 1171.582,
                              1260.304, 1248.673, 1178.612, 1168.091,
                              1337.891, 1333.631, 1264.217, 1244.832, 1182.151, 1164.615, 1107.823, 1104.613,
                              1372.717, 1342.171, 1329.391, 1268.149, 1241.008, 1185.706, 1161.155, 1111.048, 1101.416, 1079.424,
                              1368.29, 1346.472, 1325.172, 1272.099, 1237.202, 1189.276, 1157.71, 1114.286, 1098.234, 1082.525,
                              1363.884, 1350.794, 1320.972, 1276.067, 1233.413, 1192.863, 1154.28, 1117.539, 1095.065, 1085.64,
                              1359.5, 1355.136, 1316.793, 1280.054, 1229.642, 1196.466, 1150.866, 1120.806, 1091.909, 1088.768,
                              1312.633, 1284.06, 1225.888, 1200.086, 1147.466, 1124.088,
                              1308.493, 1288.084, 1222.152, 1203.722, 1144.082, 1127.384,
                              1304.373, 1292.128, 1218.432, 1207.374, 1140.713, 1130.694,
                              1300.272, 1296.19, 1214.729, 1211.043, 1137.358, 1134.019};

  const int kidIds[] = {76, 77, 78, 79, 72, 73, 74, 75, 64, 65, 66, 67, 68, 69, 70, 71, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 18, 19, 20, 21, 22, 23, 12, 13, 14, 15, 16, 17, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5};
  const int numKids = (sizeof(kidPos) / (2 * sizeof(double)));
  G4Box *kidBoxes[2 * numKids];
  G4LogicalVolume *kidVolumes[2 * numKids];
  G4PVPlacement *kidPlacements[2 * numKids];
  int index;
  float xpos;
  float ypos;
  float indLength;
  float capLength = 0.44 * mm / 2.0f;
  float kidWidth = 1.13 * mm / 2.0f;
  float capIndSpacing = 80 * um;
  float thickness = 0.03 * um / 2.0f;
  G4Colour indCol(0.0, 0.0, 1.0);
  G4Colour capCol(1.0, 0.0, 0.0);
  G4VisAttributes *indAtt = new G4VisAttributes(indCol);
  G4VisAttributes *capAtt = new G4VisAttributes(capCol);
  std::string name;
  for (int i = 0; i < numKids; i++)
  {
    name = "kid_ind" + std::to_string(kidIds[i]);
    xpos = kidPos[2 * i];
    ypos = kidPos[2 * i + 1];
    indLength = kidLength[i] / 2 * um;
    // Inductor
    kidBoxes[2 * i] = new G4Box(name, indLength, kidWidth, thickness);
    kidVolumes[2 * i] = new G4LogicalVolume(kidBoxes[2 * i], kidMat, name);
    kidVolumes[2 * i]->SetVisAttributes(indAtt);
    kidPlacements[2 * i] = new G4PVPlacement(0, G4ThreeVector(xpos * um + indLength, ypos * um + kidWidth, sThick + thickness),
                                             kidVolumes[2 * i], name, logicWorld, false, 0, checkOverlaps);
    // Capacitor
    name = "kid_cap" + std::to_string(kidIds[i]);
    kidBoxes[2 * i + 1] = new G4Box(name, capLength, kidWidth, thickness);
    kidVolumes[2 * i + 1] = new G4LogicalVolume(kidBoxes[2 * i + 1], kidMat, name);
    kidVolumes[2 * i + 1]->SetVisAttributes(capAtt);
    kidPlacements[2 * i + 1] = new G4PVPlacement(0, G4ThreeVector(xpos * um + capLength + 2.0 * indLength + capIndSpacing, ypos * um + kidWidth, sThick + thickness),
                                                 kidVolumes[2 * i + 1], name, logicWorld, false, 0, checkOverlaps);
  }

  // Feedline
  const double feedPos[] = {
      -32426,
      -14931,
      241.5,
      26969,
      -32185,
      11796,
      6937,
      241.5,
      -25247,
      -22175,
      241.5,
      34212,
      -25006,
      -22175,
      6937,
      241.5,
      -18068,
      -22175,
      241.5,
      52290,
      -17827,
      29874,
      6937,
      241.5,
      -10889,
      -31476,
      241.5,
      61591,
      -10648,
      -31476,
      6937,
      241.5,
      -3710,
      -31476,
      241.5,
      64952,
      -3469,
      33234,
      6937,
      241.5,
      3469,
      -31476,
      241.5,
      64952,
      3710,
      -31476,
      6937,
      241.5,
      10648,
      -31476,
      241.5,
      61591,
      10889,
      29874,
      6937,
      241.5,
      17827,
      -22175,
      241.5,
      52290,
      18068,
      -22175,
      6937,
      241.5,
      25006,
      -22175,
      241.5,
      34212,
      25247,
      11796,
      6937,
      241.5,
      32185,
      -14931,
      241.5,
      26969,
      -9138,
      -35270,
      241.5,
      733.725,
      8896,
      -35270,
      241.5,
      733.725,
      -9983,
      -36070,
      1931.613,
      800,
      8051,
      -36070,
      1931.613,
      800,
  };

  const double feedPos1[] = {
      -29961.722, -18715.491, 241.5, 8804.148, -58.009,
      -23692.592, -26032.247, 241.5, 10646.936, -42.306,
      -14430.254, -32017.76, 241.5, 11684.231, -24.281,
      29961.722, -18715.491, 241.5, 8804.148, 58.009,
      23692.592, -26032.247, 241.5, 10646.936, 42.306,
      14430.254, -32017.76, 241.5, 11684.231, 24.281};

  const int numFeed = sizeof(feedPos) / (4 * sizeof(double));
  const int numFeed1 = sizeof(feedPos1) / (5 * sizeof(double));
  const int totalFeed = numFeed + numFeed1 + 3 * numKids;
  float rot;
  float width;
  float height;
  G4RotationMatrix *rotMat;
  G4Box *feedBoxes[totalFeed];
  G4LogicalVolume *feedVolumes[totalFeed];
  G4PVPlacement *feedPlacements[totalFeed];
  G4Colour feedCol(1.0, 0.0, 1.0);
  G4VisAttributes *feedAtt = new G4VisAttributes(feedCol);
  for (int i = 0; i < numFeed; i++)
  {
    name = "feed" + std::to_string(i);
    xpos = feedPos[4 * i + 0] * um;
    ypos = feedPos[4 * i + 1] * um;
    width = feedPos[4 * i + 2] * um / 2;
    height = feedPos[4 * i + 3] * um / 2;
    feedBoxes[i] = new G4Box(name, width, height, thickness);
    feedVolumes[i] = new G4LogicalVolume(feedBoxes[i], feedMat, name);
    feedVolumes[i]->SetVisAttributes(feedAtt);
    feedPlacements[i] = new G4PVPlacement(0, G4ThreeVector(xpos + width, ypos + height, sThick + thickness), feedVolumes[i], name, logicWorld, false, 0, checkOverlaps);
  }
  for (int i = 0; i < numFeed1; i++)
  {
    name = "feed" + std::to_string(i + numFeed);
    xpos = feedPos1[5 * i + 0] * um;
    ypos = feedPos1[5 * i + 1] * um;
    width = feedPos1[5 * i + 2] * um / 2;
    height = feedPos1[5 * i + 3] * um / 2;
    rot = feedPos1[5 * i + 4] * deg;
    index = i + numFeed;
    feedBoxes[index] = new G4Box(name, width, height, thickness);
    feedVolumes[index] = new G4LogicalVolume(feedBoxes[index], feedMat, name);
    feedVolumes[index]->SetVisAttributes(feedAtt);
    rotMat = new G4RotationMatrix;
    rotMat->rotateZ(90 * deg - rot);
    feedPlacements[index] = new G4PVPlacement(rotMat, G4ThreeVector(xpos, ypos, sThick + thickness), feedVolumes[index], name, logicWorld, false, 0, checkOverlaps);
  }

  for (int i = 0; i < numKids; i++)
  {
    indLength = kidLength[i] / 2 * um;
    xpos = kidPos[2 * i] * um + indLength + capLength + 80 * um;
    ypos = kidPos[2 * i + 1] * um - 70 * um;

    // lower kid feed
    index = 3 * i + numFeed + numFeed1;
    name = "feed_kid_lower" + std::to_string(i);
    feedBoxes[index] = new G4Box(name, indLength + capLength + 110 * um, 40 * um, thickness);
    feedVolumes[index] = new G4LogicalVolume(feedBoxes[index], feedMat, name);
    feedVolumes[index]->SetVisAttributes(feedAtt);
    feedPlacements[index] = new G4PVPlacement(0, G4ThreeVector(xpos + um, ypos, sThick + thickness), feedVolumes[index], name, logicWorld, false, 0, checkOverlaps);

    // upper kid feed
    index = 3 * i + numFeed + numFeed1 + 1;
    name = "feed_kid_upper" + std::to_string(i);
    feedBoxes[index] = new G4Box(name, indLength + capLength + 110 * um, 40 * um, thickness);
    feedVolumes[index] = new G4LogicalVolume(feedBoxes[index], feedMat, name);
    feedVolumes[index]->SetVisAttributes(feedAtt);
    feedPlacements[index] = new G4PVPlacement(0, G4ThreeVector(xpos + um, ypos + 2 * kidWidth + 140 * um, sThick + thickness), feedVolumes[index], name, logicWorld, false, 0, checkOverlaps);

    // right kid feed
    xpos = kidPos[2 * i] * um + 2 * indLength + 2 * capLength + 150 * um;
    ypos = kidPos[2 * i + 1] * um + kidWidth;
    index = 3 * i + numFeed + numFeed1 + 2;
    name = "feed_kid_right" + std::to_string(i);
    feedBoxes[index] = new G4Box(name, 40 * um, kidWidth + 29 * um, thickness);
    feedVolumes[index] = new G4LogicalVolume(feedBoxes[index], feedMat, name);
    feedVolumes[index]->SetVisAttributes(feedAtt);
    feedPlacements[index] = new G4PVPlacement(0, G4ThreeVector(xpos, ypos, sThick + thickness), feedVolumes[index], name, logicWorld, false, 0, checkOverlaps);
  }

  // double absProb = 0.65;
  double absProb = 0.1;
  double pAbsProb = 0.1; // these have only been discriminated because that's how the original code was
  // double absProb = 0.5;
  // double absProb = 0.075;

  // 0.0, 1.0, 0.0, 0.0, pAbsProb, 1.0, 1.0, 0.0
  indSurf = new G4CMPSurfaceProperty("Surf", 0.0, 1.0, 0.0, 0.0, pAbsProb, 1.0, 1.0, 0.0);
  capSurf = new G4CMPSurfaceProperty("CapSurf", 0.0, 1.0, 0.0, 0.0, pAbsProb, 1.0, 1.0, 0.0);
  feedSurf = new G4CMPSurfaceProperty("FeedSurf", 0.0, 1.0, 0.0, 0.0, pAbsProb, 1.0, 1.0, 0.0);
  worldSurf = new G4CMPSurfaceProperty("worldSurf", 0.0, 1.0, 0.0, 0.0, 0.01, 1.0, 0.0, 0.0);

  /*
  //the following coefficients and cutoff values are not well-motivated
  //the code below is used only to demonstrate how to set these values.
  const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 1.51e-14};
  const std::vector<G4double> diffCoeffs =
    {5.88e-2, 7.83e-4, -2.47e-6, 1.71e-8, -2.98e-11};
  const std::vector<G4double> specCoeffs =
    {0,928, -2.03e-4, -3.21e-6, 3.1e-9, 2.9e-13};

  const G4double anhCutoff = 520., reflCutoff = 350.;   // Units external

  indSurf->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                   diffCoeffs, specCoeffs, GHz, GHz, GHz);

  capSurf->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                   diffCoeffs, specCoeffs, GHz, GHz, GHz);

  feedSurf->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                   diffCoeffs, specCoeffs, GHz, GHz, GHz);

  worldSurf->AddScatteringProperties(anhCutoff, reflCutoff, anhCoeffs,
                   diffCoeffs, specCoeffs, GHz, GHz, GHz);
  */

  /*indSurf = new G4CMPSurfaceProperty("Surf", 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
  capSurf = new G4CMPSurfaceProperty("CapSurf", 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
  feedSurf = new G4CMPSurfaceProperty("FeedSurf", 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);*/

  // worldSurf = new G4CMPSurfaceProperty("worldSurf", 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
  // worldSurf = new G4CMPSurfaceProperty("worldSurf", 0.0, 1.0, 0.0, 0.0, 0.1, 1.0, 0.0, 0.0);
  // worldSurf = new G4CMPSurfaceProperty("worldSurf", 0.0, 1.0, 0.0, 0.0, 0.01, 1.0, 1.0, 0.0);
  // worldSurf = new G4CMPSurfaceProperty("worldSurf", 0.0, 1.0, 0.0, 0.0, 0.008, 1.0, 0.0, 0.0);
  // worldSurf = new G4CMPSurfaceProperty("worldSurf", 0.0, 1.0, 0.0, 0.0, 0.005, 1.0, 0.0, 0.0);
  // worldSurf = new G4CMPSurfaceProperty("worldSurf", 0.0, 1.0, 0.0, 0.0, 0.001, 1.0, 0.0, 0.0);
  // worldSurf = new G4CMPSurfaceProperty("worldSurf", 0.0, 1.0, 0.0, 0.0, 0.0001, 1.0, 0.0, 0.0);
  // worldSurf = new G4CMPSurfaceProperty("worldSurf", 0.0, 1.0, 0.0, 0.0, 0.00001, 1.0, 1.0, 0.0);

  // add in properties
  G4MaterialPropertiesTable *indMatTab = indSurf->GetPhononMaterialPropertiesTablePointer();

  /*indMatTab->AddConstProperty("gapEnergy", 0.000361 / 2.0 * eV);
  indMatTab->AddConstProperty("phononLifetime", 0.242 * ns);
  indMatTab->AddConstProperty("phononLifetimeSlope", 3.14159);
  indMatTab->AddConstProperty("vSound", 6320 * m / s);
  indMatTab->AddConstProperty("filmThickness", 30 * nm);
  indMatTab->AddConstProperty("lowQPLimit", 3);
  indMatTab->AddConstProperty("subgapAbsorption", 0.01);
  indMatTab->AddConstProperty("filmAbsorption", 0.65);*/

  indMatTab->AddConstProperty("gapEnergy", 0.000361 / 2.0 * eV);
  indMatTab->AddConstProperty("phononLifetime", 0.242 * ns);
  indMatTab->AddConstProperty("phononLifetimeSlope", 3.14159);
  indMatTab->AddConstProperty("vSound", 6320 * m / s);
  indMatTab->AddConstProperty("vSubSound", 8860 * m / s);
  indMatTab->AddConstProperty("filmThickness", 30 * nm);
  // indMatTab->AddConstProperty("filmThickness", 3000 * nm);
  // indMatTab->AddConstProperty("filmThickness", 1 * mm);
  // indMatTab->AddConstProperty("filmThickness", 25 * nm);
  indMatTab->AddConstProperty("lowQPLimit", 3);
  // indMatTab->AddConstProperty("lowQPLimit", 2);
  // indMatTab->AddConstProperty("subgapAbsorption", 0.01);
  // indMatTab->AddConstProperty("subgapAbsorption", 0.1);
  indMatTab->AddConstProperty("subgapAbsorption", 0.);
  // indMatTab->AddConstProperty("subgapAbsorption", 1.);
  // indMatTab->AddConstProperty("filmAbsorption", 1);
  indMatTab->AddConstProperty("filmAbsorption", absProb);
  // indMatTab->AddConstProperty("filmAbsorption", 0.1);
  // indMatTab->AddConstProperty("qpLifetime", 438.);
  // indMatTab->AddConstProperty("qpLifetime", 12900);
  indMatTab->AddConstProperty("qpLifetime", 2350); // matches close to true qp lifetime
  // indMatTab->AddConstProperty("qpLifetime", 4000);

  G4MaterialPropertiesTable *capMatTab = capSurf->GetPhononMaterialPropertiesTablePointer();
  capMatTab->AddConstProperty("gapEnergy", 0.000361 / 2.0 * eV);
  // capMatTab->AddConstProperty("phononLifetime", 0.00417 * ns);
  capMatTab->AddConstProperty("phononLifetime", 0.242 * ns);
  capMatTab->AddConstProperty("phononLifetimeSlope", 3.14159);
  capMatTab->AddConstProperty("vSound", 3480 * m / s);
  capMatTab->AddConstProperty("vSubSound", 8860 * m / s);
  capMatTab->AddConstProperty("filmThickness", 30 * nm);
  capMatTab->AddConstProperty("lowQPLimit", 3);
  capMatTab->AddConstProperty("subgapAbsorption", 0.);
  capMatTab->AddConstProperty("filmAbsorption", absProb);
  capMatTab->AddConstProperty("qpLifetime", 0.149);

  G4MaterialPropertiesTable *feedMatTab = feedSurf->GetPhononMaterialPropertiesTablePointer();
  feedMatTab->AddConstProperty("gapEnergy", 0.00279 / 2.0 * eV);
  // feedMatTab->AddConstProperty("phononLifetime", 0.00417 * ns);
  feedMatTab->AddConstProperty("phononLifetime", 0.242 * ns);
  feedMatTab->AddConstProperty("phononLifetimeSlope", 3.14159);
  feedMatTab->AddConstProperty("vSound", 3480 * m / s);
  feedMatTab->AddConstProperty("vSubSound", 8860 * m / s);
  // feedMatTab->AddConstProperty("filmThickness", 30 * nm);
  feedMatTab->AddConstProperty("filmThickness", 300 * nm);
  feedMatTab->AddConstProperty("lowQPLimit", 3);
  feedMatTab->AddConstProperty("subgapAbsorption", 0.);
  // feedMatTab->AddConstProperty("filmAbsorption", 1);
  feedMatTab->AddConstProperty("filmAbsorption", absProb);
  // feedMatTab->AddConstProperty("filmAbsorption", 0.1);
  // feedMatTab->AddConstProperty("qpLifetime", 438.);
  feedMatTab->AddConstProperty("qpLifetime", 0.149);

  /*feedMatTab->AddConstProperty("gapEnergy", 0.00279 / 2.0 * eV);
  feedMatTab->AddConstProperty("phononLifetime", 0.00417 * ns);
  feedMatTab->AddConstProperty("phononLifetimeSlope", 3.14159);
  feedMatTab->AddConstProperty("vSound", 3480 * m / s);
  feedMatTab->AddConstProperty("filmThickness", 30 * nm);
  feedMatTab->AddConstProperty("lowQPLimit", 3);
  feedMatTab->AddConstProperty("subgapAbsorption", 0.01);
  feedMatTab->AddConstProperty("filmAbsorption", 0.65);*/

  /*feedMatTab->AddConstProperty("gapEnergy", 0.000361 / 2.0 * eV);
  feedMatTab->AddConstProperty("phononLifetime", 0.242 * ns);
  feedMatTab->AddConstProperty("phononLifetimeSlope", 3.14159);
  feedMatTab->AddConstProperty("vSound", 6320 * m / s);
  feedMatTab->AddConstProperty("filmThickness", 30 * nm);
  feedMatTab->AddConstProperty("lowQPLimit", 3);
  feedMatTab->AddConstProperty("subgapAbsorption", 0.01);
  feedMatTab->AddConstProperty("filmAbsorption", 0.65);*/

  // removed to replace with RISQ-implemented code
  // kidElectrode = new G4CMPPhononElectrode();
  // capElectrode = new G4CMPPhononElectrode();
  // feedElectrode = new G4CMPPhononElectrode();

  // // kidElectrode = new Electrode(indMatTab);
  // // kidElectrode->lattice = SiPhysical;
  // // steppingAction->kidElectrode = kidElectrode;
  // indSurf->SetPhononElectrode(kidElectrode);

  // // capElectrode = new PSElectrode(capMatTab);
  // // capElectrode = new Electrode(indMatTab);
  // // capElectrode->lattice = SiPhysical;
  // // steppingAction->capElectrode = capElectrode;
  // capSurf->SetPhononElectrode(capElectrode);

  // // feedElectrode = new Electrode(feedMatTab);
  // // feedElectrode->lattice = SiPhysical;
  // // steppingAction->feedElectrode = feedElectrode;
  // feedSurf->SetPhononElectrode(feedElectrode);

  AttachPhononSensor(indSurf);
  AttachPhononSensor(capSurf);
  AttachPhononSensor(feedSurf);
  AttachPhononSensor(worldSurf);  // seems to be needed for detection of all particles

  for (int i = 0; i < numKids; i++)
  {
    new G4CMPLogicalBorderSurface("KidIndSurf" + std::to_string(i), SiPhys, kidPlacements[2 * i], indSurf);
    new G4CMPLogicalBorderSurface("KidCapSurf" + std::to_string(i), SiPhys, kidPlacements[2 * i + 1], capSurf); // CAUTION:  this was originally `indSurf`.  I changed it to `capSurf` as it seemed accidental
    // new G4LogicalBorderSurface("KidCapSurf" + std::to_string(i), SiPhys, kidPlacements[2 * i + 1], capSurf);
  }

  for (int i = 0; i < totalFeed; i++)
  {
    new G4CMPLogicalBorderSurface("FeedSurf" + std::to_string(i), SiPhys, feedPlacements[i], feedSurf);
  }

  new G4CMPLogicalBorderSurface("WorldSurface", SiPhys, fWorldPhys, worldSurf);
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

  G4SDManager *SDman = G4SDManager::GetSDMpointer();
  if (!fSuperconductorSensitivity)
    fSuperconductorSensitivity = new RISQ_Sensitivity("PhononElectrode");
  SDman->AddNewDetector(fSuperconductorSensitivity);
  logicTarget->SetSensitiveDetector(fSuperconductorSensitivity);

  // begin RISQ stuff

  // // the following coefficients and cutoff values are not well-motivated
  // // the code below is used only to demonstrate how to set these values.
  // const std::vector<G4double> anhCoeffs = {0, 0, 0, 0, 0, 0};  // Turn this off temporarily
  // const std::vector<G4double> diffCoeffs = {1, 0, 0, 0, 0, 0}; // Explicitly make this 1 for now
  // const std::vector<G4double> specCoeffs = {0, 0, 0, 0, 0, 0}; // Turn this off temporarily
  // const G4double anhCutoff = 520., reflCutoff = 350.;          // Units external
}

// Set up a phonon sensor for this surface property object. I'm pretty sure that this
// phonon sensor doesn't get stapled to individual geometrical objects, but rather gets
// stapled to a surface property, but I'm not sure... have to ask mKelsey
void RISQ_Detector::AttachPhononSensor(G4CMPSurfaceProperty *surfProp)
{
  // If no surface, don't do anything
  if (!surfProp)
    return;

  // Specify properties of the niobium sensors
  auto sensorProp = surfProp->GetPhononMaterialPropertiesTablePointer();
  sensorProp->AddConstProperty("filmAbsorption", 0.0);                  // NOT WELL MOTIVATED - probably parametrize and put on slider?
  sensorProp->AddConstProperty("filmThickness", 90. * CLHEP::nm);       // Accurate for our thin film.
  sensorProp->AddConstProperty("gapEnergy", 1.6e-3 * CLHEP::eV);        // Reasonably motivated. Actually, looks like Novotny and Meincke are quoting 2Delta, and this is delta. Nuss and Goossen mention that Nb has a delta value closer to this.
  sensorProp->AddConstProperty("lowQPLimit", 3.);                       // NOT WELL MOTIVATED YET -- Dunno how to inform this...
  sensorProp->AddConstProperty("phononLifetime", 4.17 * CLHEP::ps);     // Kaplan paper says 242ps for Al, same table says 4.17ps for characteristic time for Nb.
  sensorProp->AddConstProperty("phononLifetimeSlope", 0.29);            // Based on guessing from Kaplan paper, I think this is material-agnostic?
  sensorProp->AddConstProperty("vSound", 3.480 * CLHEP::km / CLHEP::s); // True for room temperature, probably good to 10%ish - should follow up
  sensorProp->AddConstProperty("subgapAbsorption", 0.0);                // Assuming that since we're mostly sensitive to quasiparticle density, phonon "heat" here isn't something that we're sensitive to? Unsure how to select this.

  //  sensorProp->AddConstProperty("gapEnergy",3.0e-3*CLHEP::eV);      //Reasonably motivated. Novotny and Meincke, 1975 (2.8-3.14 meV)
  //  sensorProp->AddConstProperty("phononLifetime",242.*ps);      //Kaplan paper says 242ps for Al, same table says 4.17ps for characteristic time for Nb.

  surfProp->SetPhononElectrode(new G4CMPPhononElectrode);
}

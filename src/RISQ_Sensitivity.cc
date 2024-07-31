/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#include "Config_Manager.hh"

#include "RISQ_Sensitivity.hh"
#include "G4CMPElectrodeHit.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Navigator.hh"
#include <fstream>
#include <iostream>

RISQ_Sensitivity::RISQ_Sensitivity(G4String name) : G4CMPElectrodeSensitivity(name), primaryFileName(""), hitFileName("")
{
  SetHitOutputFile(Config_Manager::GetHitOutput());
  SetPrimaryOutputFile(Config_Manager::GetPrimaryOutput());
}

RISQ_Sensitivity::~RISQ_Sensitivity()
{
  // Close file and check: primaries
  if (primaryOutput.is_open())
    primaryOutput.close();
  if (!primaryOutput.good())
  {
    G4cerr << "Error closing primary output file, " << primaryFileName << ".\n"
           << "Expect bad things like loss of data.";
  }

  // Close file and check: hits
  if (hitOutput.is_open())
    hitOutput.close();
  if (!hitOutput.good())
  {
    G4cerr << "Error closing hit output file, " << hitFileName << ".\n"
           << "Expect bad things like loss of data.";
  }
}

void RISQ_Sensitivity::EndOfEvent(G4HCofThisEvent *HCE)
{
  G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection);
  auto *hitCol = static_cast<G4CMPElectrodeHitsCollection *>(HCE->GetHC(HCID));
  std::vector<G4CMPElectrodeHit *> *hitVec = hitCol->GetVector();

  G4RunManager *runMan = G4RunManager::GetRunManager();

  // Do primary output writing to file
  if (primaryOutput.good())
  {
    primaryOutput << runMan->GetCurrentRun()->GetRunID() << ","                                                                     // run
                  << runMan->GetCurrentEvent()->GetEventID() << ","                                                                 // event
                  << runMan->GetCurrentEvent()->GetPrimaryVertex()->GetPrimary()->GetParticleDefinition()->GetParticleName() << "," // particle
                  << runMan->GetCurrentEvent()->GetPrimaryVertex()->GetPrimary()->GetTotalEnergy() / eV << ","                      // e_i
                  << runMan->GetCurrentEvent()->GetPrimaryVertex()->GetX0() / mm << ","                                             // x_i
                  << runMan->GetCurrentEvent()->GetPrimaryVertex()->GetY0() / mm << ","                                             // y_i
                  << runMan->GetCurrentEvent()->GetPrimaryVertex()->GetZ0() / mm << ","                                             // z_i
                  << runMan->GetCurrentEvent()->GetPrimaryVertex()->GetT0() / ns << "\n";                                           // t_i
  }

  // Do hit output writing to file
  if (hitOutput.good())
  {
    for (G4CMPElectrodeHit *hit : *hitVec)
    {
      hitOutput << runMan->GetCurrentRun()->GetRunID() << ','     // run
                << runMan->GetCurrentEvent()->GetEventID() << ',' // event
                << hit->GetTrackID() << ','                       // track
                << hit->GetParticleName() << ','                  // type
                << hit->GetStartEnergy() / eV << ','              // e_i
                << hit->GetStartPosition().getX() / mm << ','     // x_i
                << hit->GetStartPosition().getY() / mm << ','     // y_i
                << hit->GetStartPosition().getZ() / mm << ','     // z_i
                << hit->GetStartTime() / ns << ','                // t_i
                << hit->GetEnergyDeposit() / eV << ','            // e_dep
                // << hit->GetWeight() << ','
                << hit->GetFinalPosition().getX() / mm << ',' // x_f
                << hit->GetFinalPosition().getY() / mm << ',' // y_f
                << hit->GetFinalPosition().getZ() / mm << ',' // z_f
                << hit->GetFinalTime() / ns << '\n';          // t_f
    }
  }
}

void RISQ_Sensitivity::SetHitOutputFile(const G4String &fn)
{
  if (hitFileName != fn)
  {
    if (hitOutput.is_open())
      hitOutput.close();
    hitFileName = fn;
    hitOutput.open(hitFileName, std::ios_base::trunc);
    if (!hitOutput.good())
    {
      G4ExceptionDescription msg;
      msg << "Error opening hit output file " << hitFileName;
      G4Exception("RISQ_Sensitivity::SetHitOutputFile", "PhonSense003",
                  FatalException, msg);
      hitOutput.close();
    }
    else
    {
      hitOutput << "run,event,track,type," // note:  removed `weight`
                << "e_i,x_i,y_i,z_i,t_i,"
                << "e_dep,x_f,y_f,z_f,t_f\n";
    }
  }

  // need to write to `volume_hits.log` first line with header information
}

void RISQ_Sensitivity::SetPrimaryOutputFile(const G4String &fn)
{
  if (primaryFileName != fn)
  {
    if (primaryOutput.is_open())
      primaryOutput.close();
    primaryFileName = fn;
    primaryOutput.open(primaryFileName, std::ios_base::trunc);
    if (!primaryOutput.good())
    {
      G4ExceptionDescription msg;
      msg << "Error opening output file " << primaryFileName;
      G4Exception("RISQ_Sensitivity::SetPrimaryOutputFile", "PhonSense003",
                  FatalException, msg);
      primaryOutput.close();
    }
    else
    {
      primaryOutput << "run,event,type,"
                    << "e_i,x_i,y_i,z_i,t_i\n";
    }
  }
}

G4bool RISQ_Sensitivity::IsHit(const G4Step *step,
                               const G4TouchableHistory *) const
{

  // Establish track/step information
  G4RunManager *runMan = G4RunManager::GetRunManager();
  const G4Track *track = step->GetTrack();
  const G4StepPoint *postStepPoint = step->GetPostStepPoint();
  const G4ParticleDefinition *particle = track->GetDefinition();

  //-------------------------------------------------------------------
  // Set criterion for what counts as a "hit" that should be recorded.
  bool selectTargetVolumes = false;

  // Option one: a phonon that is stopped and killed at a boundary with a
  // nonzero energy deposition.
  G4bool correctParticle = particle == G4PhononLong::Definition() ||
                           particle == G4PhononTransFast::Definition() ||
                           particle == G4PhononTransSlow::Definition();

  G4TrackStatus status = step->GetTrack()->GetTrackStatus();
  G4StepStatus stepStatus = postStepPoint->GetStepStatus();
  G4double energyDeposit = step->GetNonIonizingEnergyDeposit();

  G4bool correctStatus = status == fStopAndKill &&
                         stepStatus == fGeomBoundary &&
                         energyDeposit > 0;

  G4String volumeName = postStepPoint->GetPhysicalVolume()->GetName();
  G4bool landedOnTargetSurface = (volumeName.find("kid_ind") != std::string::npos);

  if (correctParticle && status == fStopAndKill && stepStatus == fGeomBoundary)
  {
    std::ofstream hitsLog;
    hitsLog.open("hits.csv", std::ios_base::app);                                            // header is hardcoded in main.cc
    hitsLog << runMan->GetCurrentRun()->GetRunID() << ','                                    // run
            << runMan->GetCurrentEvent()->GetEventID() << ','                                // event
            << track->GetTrackID() << ','                                                    // track
            << track->GetParticleDefinition()->GetParticleName() << ','                      // type
            << volumeName << ','                                                             // volume
            << track->GetVertexKineticEnergy() / eV << ','                                   // e_i
            << track->GetVertexPosition().getX() / mm << ','                                 // x_i
            << track->GetVertexPosition().getY() / mm << ','                                 // y_i
            << track->GetVertexPosition().getZ() / mm << ','                                 // z_i
            << step->GetPostStepPoint()->GetGlobalTime() - track->GetLocalTime() / ns << ',' // t_i
            << step->GetNonIonizingEnergyDeposit() / eV << ','                               // e_dep
            << step->GetPostStepPoint()->GetPosition().getX() / mm << ','                    // x_f
            << step->GetPostStepPoint()->GetPosition().getY() / mm << ','                    // y_f
            << step->GetPostStepPoint()->GetPosition().getZ() / mm << ','                    // z_f
            << step->GetPostStepPoint()->GetGlobalTime() / ns << '\n';                       // t_f

    hitsLog.close();
  }

  // Now select which critera matter:
  // Option one: a phonon that is stopped and killed at a boundary with a
  // nonzero energy deposition.
  if (!selectTargetVolumes)
  {
    return correctParticle && correctStatus;
  }

  // Option two: a phonon that satisfies all of the above things, but also landed in a specific
  // volume name. Here, we're looking for a volume that contains the words "shuntConductor", which
  // in this tutorial's geometry is one of the qubit crosses. (Can also just put this info in
  // the output file and sort through this in analysis, but this helps us minimize output filesize.)
  else
  {
    return correctParticle && correctStatus && landedOnTargetSurface;
  }
}

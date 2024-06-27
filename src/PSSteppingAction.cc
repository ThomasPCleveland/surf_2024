#include "PSSteppingAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4StateManager.hh"
#include "G4PhononLong.hh"
#include "G4PhononTransFast.hh"
#include "G4PhononTransSlow.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"

PSSteppingAction::PSSteppingAction()
    : G4UserSteppingAction()
{
    //stores data of all geant4/g4cmp particles.
    stepData = new std::string();
    //Sometimes, the electron created by the photon does not conserve energy
    //and creates ~61meV too much energy in the form of secondary phonons.
    //When this happens, endEvent can be set to change to true or false.
    //if it set to false on line ~180, events that do not conserve energy will not be discarded,
    //meaning their data will be recorded in the files directly below. If set to true,
    //then the data will not be recorded and all particles will be killed after the electron
    //creates the secondaries with too much energy.
    endEvent = false;
    //if any of the data variables exceed this size in bytes, they will be written to the data files
    //and cleared to avoid extremely large data variables (>1MB).
    maxData = 1000000;

    //these files store the cumulative data for all events
    remove("StepData.txt");
    stepFile.open("StepData.txt", std::ios_base::app);
    remove("ElectrodeData.txt");
    electrodeFile.open("ElectrodeData.txt", std::ios_base::app);
    remove("ReflData.txt");
    reflFile.open("ReflData.txt", std::ios_base::app);

    //these files only store the data for the current event that is being run.
    //after the event ends, the data may or may not be appended to the cumulative
    //files depending on the value of endEvent.
    stepFileTemp.open("StepDataTemp.txt", std::ios_base::app);
    electrodeFileTemp.open("ElectrodeDataTemp.txt", std::ios_base::app);
    reflFileTemp.open("ReflDataTemp.txt", std::ios_base::app);
}

PSSteppingAction::~PSSteppingAction()
{
    //remove temporary data files which only store data for 1 event.
    remove("StepDataTemp.txt");
    remove("ElectrodeDataTemp.txt");
    remove("ReflDataTemp.txt");
}

void PSSteppingAction::WriteData()
{
    //write data to temporary files
    stepFileTemp << *(stepData);
    electrodeFileTemp << *(feedElectrode->GetData());
    electrodeFileTemp << *(capElectrode->GetData());
    electrodeFileTemp << *(kidElectrode->GetData());

    reflFileTemp << *(feedElectrode->reflData);
    reflFileTemp << *(capElectrode->reflData);
    reflFileTemp << *(kidElectrode->reflData);

    //clear data variables afterwards.
    stepData->clear();
    feedElectrode->ClearData();
    capElectrode->ClearData();
    kidElectrode->ClearData();

    feedElectrode->reflData->clear();
    capElectrode->reflData->clear();
    kidElectrode->reflData->clear();
}

void PSSteppingAction::NewEvent()
{
    WriteData();

    //copy data from temp files to the cumulative files at end of each event.

    //close files to avoid reading and writing at same time.
    stepFileTemp.close();
    electrodeFileTemp.close();
    reflFileTemp.close();

    //allows us to read from temp files
    std::ifstream stepTemp("StepDataTemp.txt");
    std::ifstream electrodeTemp("ElectrodeDataTemp.txt");
    std::ifstream reflTemp("ReflDataTemp.txt");

    //loops through each line in temp file and appends it to cumulative file
    if (!endEvent)
    {
        std::string line;
        while (std::getline(stepTemp, line))
        {
            stepFile << line << "\n";
        }
        while (std::getline(electrodeTemp, line))
        {
            electrodeFile << line << "\n";
        }
        while (std::getline(reflTemp, line))
        {
            reflFile << line << "\n";
        }

        //denotes where the data for each event begins and ends
        stepFile << "New Event\n";
        electrodeFile << "New Event\n";
    }
    //no longer need to read from files.
    stepTemp.close();
    electrodeTemp.close();
    reflTemp.close();

    //trunc is used to clear the data from the previous event
    stepFileTemp.open("StepDataTemp.txt", std::ios_base::trunc);
    electrodeFileTemp.open("ElectrodeDataTemp.txt", std::ios_base::trunc);
    reflFileTemp.open("ReflDataTemp.txt", std::ios_base::trunc);

    //set back to false so next event is not ended if the previous one was.
    endEvent = false;
}

void PSSteppingAction::UserSteppingAction(const G4Step* step)
{
    //*
    // Get the track associated with this step
    G4Track* track = step->GetTrack();
    G4int trackID = track->GetTrackID();
    G4int parentId = step->GetTrack()->GetParentID();

    // See if we are looking at a phonon. (If not, continue)
    const G4ParticleDefinition* particle = track->GetDefinition();
    G4bool isPhonon = true;
    G4String particleName = "";
    if (particle == G4PhononLong::Definition()) particleName = "phononL";
    else if (particle == G4PhononTransFast::Definition()) particleName = "phononTF";
    else if (particle == G4PhononTransSlow::Definition()) particleName = "phononTS";
    else isPhonon = false;

    //if (the particle is a phonon and it is its first/last step) OR (the particle is not a phonon), then its step data is recorded
    G4bool particleBirth = step->GetTrack()->GetCurrentStepNumber() == 1;
    G4bool particleDeath = step->GetTrack()->GetTrackStatus() == fStopAndKill;


    // Save the information to StepData.txt if phonon
    if (isPhonon && (particleBirth || particleDeath)) {

        // Find time, energy, and volume, at birth or at death
        G4double kineticEnergy = 0;
        G4double globalTime = 0;
        G4String volumeName = "";
        G4String processName = "";
        G4double depositedEnergy = step->GetNonIonizingEnergyDeposit();
        const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
        if (particleBirth) {
            // Access the pre-step point of the step to get KE
            const G4StepPoint* preStepPoint = step->GetPreStepPoint();
            kineticEnergy = preStepPoint->GetKineticEnergy()/eV;
            globalTime = preStepPoint->GetGlobalTime();
            volumeName = preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
        } else {
            // Access the post-step point of the step to get KE
            const G4StepPoint* postStepPoint = step->GetPostStepPoint();
            kineticEnergy = postStepPoint->GetKineticEnergy()/eV;
            globalTime = postStepPoint->GetGlobalTime();
            volumeName = postStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
        }
        if (process) processName = process->GetProcessName();


        if (particleBirth) {
            stepData->append(std::to_string(trackID) + ",");
            stepData->append(std::to_string(globalTime) + ",");
            stepData->append(std::to_string(kineticEnergy) + ",");
            stepData->append(volumeName + ",");
        } 
        if (particleDeath) {
            stepData->append(std::to_string(globalTime) + ",");
            stepData->append(std::to_string(kineticEnergy) + ",");
            stepData->append(volumeName + ",");
            stepData->append(std::to_string(parentId) + ",");
            stepData->append(particleName + ",");
            stepData->append(std::to_string(depositedEnergy) + "\n");
        }

    } 
    //*/


    /*
    const G4Track* track = step->GetTrack();
    const G4ParticleDefinition* particle = track->GetDefinition();
    //Check that it is a phonon
    G4bool correctParticle = particle == G4PhononLong::Definition() ||
        particle == G4PhononTransFast::Definition() ||
        particle == G4PhononTransSlow::Definition();

    //particles after this time are killed
    //G4double maxTime = 5000000
    G4double maxTime = 1000000;
    //how much energy did the particle lose/deposit this step?
    G4double deposit = 1000.0 * step->GetNonIonizingEnergyDeposit() / eV;
    //what is the particle's global time?
    G4double time = step->GetTrack()->GetGlobalTime();
    //difference in time between prestep point and poststep point
    G4double deltaTime = step->GetDeltaTime();
    //id of particle
    G4int particleId = track->GetTrackID();
    //type of particle
    std::string particleName = step->GetTrack()->GetDefinition()->GetParticleName();
    //id of parent particle (the particle that created it)
    G4int parentId = step->GetTrack()->GetParentID();
    //KE of particle
    G4double kineticEnergy = 1000.0 * step->GetTrack()->GetKineticEnergy() / eV;
    //change in energy from prestep point to poststep point.
    G4double deltaEnergy = 1000.0 * step->GetDeltaEnergy() / eV;
    //name of volume the particle is in
    std::string volName;
    if (step->GetPostStepPoint()->GetTouchableHandle()->GetVolume())
    {
        volName = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
    }
    else
    {
        //if the volume does not exist, it is assumed that it is in the world since it escaped from all other geometries.
        volName = "World";
        G4cout << "Volume not found!\n";
    }
    //get total secondary energy of particle on this step.
    const G4TrackVector* tv = step->GetSecondary();
    G4double secEnergy = 0.0;
    for (int i = 0; i < (int)tv->size(); i++)
    {
        G4Track const* secTrack = tv->operator[](i);
        secEnergy += 1000.0 * secTrack->GetKineticEnergy() / eV;
    }
    //if the particle exists after the max time, or the event should be ended due to non-conservation of energy,
    //the particle is killed.
    if(time > maxTime || endEvent)
    {
        step->GetTrack()->SetTrackStatus(fStopAndKill);
    }
    //record the data at the prestep point of first step of every particle
    if (step->GetTrack()->GetCurrentStepNumber() == 1)
    {
        stepData->append(std::to_string(particleId) + ",");
        stepData->append(std::to_string(kineticEnergy - deltaEnergy) + ",");
        stepData->append(std::to_string(time - deltaTime) + ",");
        stepData->append(particleName + ",");
        stepData->append(std::to_string(0.0) + ",");
        stepData->append(volName + ",");
        stepData->append(std::to_string(parentId) + "\n");
    }
    //if (the particle is a phonon and it is its first/last step) OR (the particle is not a phonon), then its step data is recorded
    if (!correctParticle || step->GetTrack()->GetCurrentStepNumber() == 1 || step->GetTrack()->GetTrackStatus() == fStopAndKill)
    {
        stepData->append(std::to_string(particleId) + ",");
        stepData->append(std::to_string(kineticEnergy) + ",");
        stepData->append(std::to_string(time) + ",");
        stepData->append(particleName + ",");
        stepData->append(std::to_string(deposit) + ",");
        stepData->append(volName + "\n");
    } 

    //is the particle killed and an electron?
    if (step->GetTrack()->GetTrackStatus() == fStopAndKill && particle == G4Electron::Definition())
    {
        //is the particle conserving energy?
        if (secEnergy > kineticEnergy - deltaEnergy)
        {
            //set to false to keep data for events that do not conserve energy.
            //I recommend this because only a small amount of energy is added to simulation (~61meV)
            //and this non-conservation bug happens to a lot of events (so a large fraction would be
            //discarded if set to true, though the events would take less time to run because particles
            //are killed after this happens).

            //set to true to discard data for events that do not conserve energy
            endEvent = false;
            //endEvent = true;
            G4cout << "Electron not conserving energy. Energy: " << kineticEnergy - deltaEnergy << "meV, Secondary Energy: " << secEnergy << "meV,";
            G4cout << "Difference: " << secEnergy - (kineticEnergy - deltaEnergy) << "meV\n";
        }
    }
    //*/
    
    //if the size of any of the data variables are greater than the maximum data size, they will be written to the data file
    //and cleared.
    if (stepData->length() > maxData || feedElectrode->GetDataLength() > maxData || capElectrode->GetDataLength() > maxData || kidElectrode->GetDataLength() > maxData)
    {
        WriteData();
    }
}
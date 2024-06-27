#ifndef PSSteppingAction_h
#define PSSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "PSElectrode.hh"
#include "globals.hh"

#include <vector>
#include <iostream>
#include <fstream>

class PSSteppingAction : public G4UserSteppingAction
{
public:
    PSSteppingAction();
    virtual ~PSSteppingAction();

    // method from the base class
    void WriteData();
    void NewEvent();
    virtual void UserSteppingAction(const G4Step*);
    PSElectrode* kidElectrode;
    PSElectrode* capElectrode;
    PSElectrode* feedElectrode;

private:
    mutable std::ofstream stepFile;
    std::ofstream electrodeFile;
    std::ofstream reflFile;

    mutable std::ofstream stepFileTemp;
    std::ofstream electrodeFileTemp;
    std::ofstream reflFileTemp;

    mutable std::string* stepData;
    mutable bool endEvent;
    size_t maxData;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

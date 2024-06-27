#ifndef PSElectrode_hh
#define PSElectrode_hh 1
#include "G4CMPVElectrodePattern.hh"
#include "PSKaplanQP.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include <string>
#include <fstream>
#include <vector>

class G4ParticleChange;
class G4Step;
class G4Track;

class PSElectrode : public G4CMPVElectrodePattern {
public:
    PSElectrode(G4MaterialPropertiesTable* prop);
    virtual ~PSElectrode();
    // Implement cloning function along with default copiers
    PSElectrode(const PSElectrode&) = default;
    PSElectrode(PSElectrode&&) = default;
    PSElectrode& operator=(const PSElectrode&) = default;
    PSElectrode& operator=(PSElectrode&&) = default;
    virtual G4CMPVElectrodePattern* Clone() const {
        return new PSElectrode(*this);
    }

    virtual G4bool IsNearElectrode(const G4Step&) const;
    virtual void AbsorbAtElectrode(const G4Track&,
        const G4Step&,
        G4ParticleChange&) const;
    std::string* GetData() { return elecData; };
    size_t GetDataLength() { return elecData->length(); };
    void ClearData() { elecData->clear(); };
    G4LatticePhysical* lattice;
    mutable std::string* reflData;
    mutable std::string* escapeData; // stores energies of phonons that escape the films
    
private:
    G4ThreeVector SpecularReflection(G4ThreeVector primaryPhononDir, G4ThreeVector surfNorm) const;
    G4MaterialPropertiesTable* matTab;
    PSKaplanQP* kaplanQP;
    mutable std::string* elecData;

};

#endif

#ifndef PSKaplanQP_hh
#define PSKaplanQP_hh 1

#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include <fstream>
#include <vector>

class G4MaterialPropertiesTable;

#include "PSPhonon.hh"
#include "PSQP.hh"

class PSKaplanQP {
public:
    PSKaplanQP(G4MaterialPropertiesTable* prop, G4int vb = 0);
    virtual ~PSKaplanQP();

    // Turn on diagnostic messages
    void SetVerboseLevel(G4int vb) { verboseLevel = vb; }
    G4int GetVerboseLevel() const { return verboseLevel; }

    // Configure thin film (QET, metalization, etc.) for phonon absorption
    void SetFilmProperties(G4MaterialPropertiesTable* prop);

    // Do absorption on sensor/metalization film
    // Returns absorbed energy, fills list of re-emitted phonons
    G4double AbsorbPhonon(PSPhonon phonon,
        std::vector<PSPhonon>& reflectedEnergies) const;

    std::string* GetQPData() { return qpData; };
    std::string* GetEnergyGainData() { return energyGainData; };
    size_t GetQPDataLength() { return qpData->length(); };

protected:
    // Compute the probability of a phonon reentering the crystal without breaking
    // any Cooper pairs.
    G4double CalcEscapeProbability(PSPhonon phonon,
        G4double thicknessFrac) const;

    // Model the phonons (phonEnergies) breaking Cooper pairs into quasiparticles
    // (qpEnergies).
    G4double CalcQPEnergies(std::vector<PSPhonon>& phonEnergies,
        std::vector<PSQP>& qpEnergies) const;

    // Model the quasiparticles (qpEnergies) emitting phonons (phonEnergies) in
    // the superconductor.
    G4double CalcPhononEnergies(std::vector<PSPhonon>& phonEnergies,
        std::vector<PSQP>& qpEnergies) const;

    // Calculate energies of phonon tracks that have reentered the crystal.
    void CalcReflectedPhononEnergies(std::vector<PSPhonon>& phonEnergies,
        std::vector<PSPhonon>& reflectedEnergies) const;

    // Compute probability of absorbing phonon below Cooper-pair breaking
    // NOTE:  Function should ONLY be called for energy < 2.*gapEnergy
    G4double CalcSubgapAbsorption(PSPhonon phonon,
        std::vector<PSPhonon>& keepEnergies) const;

    // Handle absorption of quasiparticle energies below Cooper-pair breaking
    // If qpEnergy < 3*Delta, radiate a phonon, absorb bandgap minimum
    G4double CalcQPAbsorption(PSQP qp,
        std::vector<PSPhonon>& phonEnergies,
        std::vector<PSQP>& qpEnergies) const;

    // Compute quasiparticle energy distribution from broken Cooper pair.
    PSQP QPEnergyRand(PSPhonon phonon) const;
    G4double QPEnergyPDF(G4double E, G4double x) const;

    // Compute phonon energy distribution from quasiparticle in superconductor.
    PSPhonon PhononEnergyRand(PSQP qp) const;
    G4double PhononEnergyPDF(G4double E, G4double x) const;

    // Encapsulate below-bandgap logic
    G4bool IsSubgap(PSPhonon phonon) const { return phonon.GetEnergy() < 2. * gapEnergy; }

    //G4double CalcQPLifetime();
private:


    G4int verboseLevel;		// For diagnostic messages

    G4MaterialPropertiesTable* filmProperties;
    G4double filmAbsorption;
    G4double filmThickness;	// Quantities extracted from properties table
    G4double gapEnergy;		// Bandgap energy (delta)
    G4double lowQPLimit;		// Minimum X*delta to keep as a quasiparticle
    G4double subgapAbsorption;	// Probability to absorb energy below bandgap
    G4double phononLifetime;	// Lifetime of phonons in film at 2*delta
    G4double phononLifetimeSlope;	// Energy dependence of phonon lifetime
    G4double vSound;		// Speed of sound in film
    G4double vSubSound;		// Speed of sound in substrate
    G4double zmin;		// minimum value for z component of phonons.
    G4double thetaCrit;		//Critical angle for total internal reflection of phonons.
    G4double tirProb;		//this variable represents the maximum value for the z component
                            // of a phonon that will undergo total internal reflection.
                            // since the z components of the phonons are distributed uniformly
                            // from -1 to 1 when created from QP recombination, this variable
                            // also represents the probability that a random phonon will be 
                            // subject to tir, hence the name. If the way the phonon directions
                            // are chosen is changed, it might no longer represent this prob.
                            //this value can be set manually in PSKaplanQP.cc so that it is not
                            // calculated using snell's law. This provides a helpful knob to turn
                            // in the simulation
    G4double qpLifetime;    //average qp lifetime in film
    G4ThreeVector norm;		//Normal of film-substrate boundary;
    mutable int qpId;       //stores the id of the last QP created

    mutable std::string* qpData;            //stores all QP data (time, energy, id) of every QP.
    mutable std::string* energyGainData;    //stores what energy is lost and gained by getting 
                                            //a quiescent delta QP from film or losing a delta QP
                                            //into the film.

};

#endif	/* PSKaplanQP_hh */

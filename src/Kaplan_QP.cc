#include "globals.hh"
#include "Kaplan_QP.hh"
#include "G4CMPConfigManager.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <math.h>
#include <numeric>
#include <algorithm>


// Class constructor and destructor

Kaplan_QP::Kaplan_QP(G4MaterialPropertiesTable* prop, G4int vb)
    : verboseLevel(vb), filmProperties(0), filmThickness(0.), gapEnergy(0.),
    lowQPLimit(3.), subgapAbsorption(0.), phononLifetime(0.),
    phononLifetimeSlope(0.), vSound(0.), vSubSound(0.), zmin(0.00001), thetaCrit(pi), 
    tirProb(0.9), qpLifetime(438.), norm(0, 0, 1), qpId(1) {
    SetFilmProperties(prop);
}

Kaplan_QP::~Kaplan_QP() {
#ifdef G4CMP_DEBUG
    if (output.is_open()) output.close();
#endif
}

// Configure thin film (QET, metalization, etc.) for phonon absorption

void Kaplan_QP::SetFilmProperties(G4MaterialPropertiesTable* prop) {
    if (!prop) {
        G4Exception("Kaplan_QP::SetFilmProperties()", "G4CMP001",
            RunMustBeAborted, "Null MaterialPropertiesTable vector.");
    }

    // Check that the MaterialPropertiesTable has everything we need. If it came
    // from a G4CMPSurfaceProperty, then it will be fine.
    if (!(prop->ConstPropertyExists("gapEnergy") &&
        prop->ConstPropertyExists("phononLifetime") &&
        prop->ConstPropertyExists("phononLifetimeSlope") &&
        prop->ConstPropertyExists("vSound") &&
        prop->ConstPropertyExists("filmThickness"))) {
        G4Exception("Kaplan_QP::SetFilmProperties()", "G4CMP002",
            RunMustBeAborted,
            "Insufficient info in MaterialPropertiesTable.");
    }

    // Extract values from table here for convenience in functions
    if (filmProperties != prop) {
        filmThickness = prop->GetConstProperty("filmThickness");
        gapEnergy = prop->GetConstProperty("gapEnergy");
        phononLifetime = prop->GetConstProperty("phononLifetime");
        phononLifetimeSlope = prop->GetConstProperty("phononLifetimeSlope");
        vSound = prop->GetConstProperty("vSound");
        vSubSound = prop->GetConstProperty("vSubSound");

        lowQPLimit = (prop->ConstPropertyExists("lowQPLimit")
            ? prop->GetConstProperty("lowQPLimit") : 3.);

        subgapAbsorption = (prop->ConstPropertyExists("subgapAbsorption")
            ? prop->GetConstProperty("subgapAbsorption") : 0.);

        qpLifetime = (prop->ConstPropertyExists("qpLifetime")
            ? prop->GetConstProperty("qpLifetime") : 438.);
        filmAbsorption = prop->GetConstProperty("filmAbsorption");

        filmProperties = prop;
    }

    //calculate critical angle
    if (vSound < vSubSound)
    {
        //use this line of code to calculate critical angle from snell's law.
        //thetaCrit = asin(vSound / vSubSound);
        
        //use this line of code to adjust the probability of TIR using tirProb.
        //the value for tirProb is set in the constructor of this class.
        thetaCrit = asin(sqrt(1 - tirProb * tirProb));
    }
    //thetaCrit = pi;
}

//commented by Mike Kelsey. This is somewhat misleading because the physics has been altered
// from the original G4CMPKaplanQP.cc where this code (and this comment) was copied from.
// for example, there are no longer any "energy deposits" (this only applies to TES-like devices
// which this code was originally used for).

// This is the main function for the Kaplan quasiparticle downconversion
// process. Based on the energy of the incoming phonon and the properties
// of the superconductor, we return the total energy deposited as well
// as fill a vector of energies that correspond to newly created phonons
// that are emitted back into the crystal.

G4double Kaplan_QP::
AbsorbPhonon(Phonon primary, std::vector<Phonon>& reflectedEnergies) const {
    //initialize qpData and energyGainData.
    //"|" is used as a delimiter.
    qpData = new std::string("|");
    energyGainData = new std::string("|");
    //get energy of primary phonon
    G4double energy = primary.GetEnergy();
    //set minimum z component to zmin
    primary.SetZMin(zmin);

    //errors/warnings
    if (!filmProperties) {
        G4Exception("Kaplan_QP::AbsorbPhonon()", "G4CMP001",
            RunMustBeAborted, "Null MaterialPropertiesTable vector.");
    }

    if (verboseLevel)
        G4cout << "Kaplan_QP::AbsorbPhonon " << primary.GetEnergy() << G4endl;

    if (reflectedEnergies.size() > 0) {
        G4Exception("Kaplan_QP::AbsorbPhonon", "G4CMP007", JustWarning,
            "Passed a nonempty reflectedEnergies vector.");
        // FIXME: Should we discard previous contents?
    }

#ifdef G4CMP_DEBUG
    if (!output.good()) {
        output.open("kaplanqp_stats");
        if (output.good()) {
            output << "Incident Energy [eV],Absorbed Energy [eV],"
                << "Reflected Energy [eV],Reflected Phonons" << std::endl;
        }
    }
#endif

    // For the phonon to not break a Cooper pair, it must go 2*thickness vertically
    G4double frac = 2.0;

    // If phonon is does not create QPs, reflect it back into crystal.

    //phonon is subgap
    if (IsSubgap(primary)) {
        primary.ReflectDirection();
        //not important, just ends function and adds phonon to list of reflected phonons.
        return CalcSubgapAbsorption(primary, reflectedEnergies);
    }
    else if (G4UniformRand() <= CalcEscapeProbability(primary, frac)) { //phonon is not subgap, does not create QPs
        if (verboseLevel > 1) G4cout << " Not absorbed." << G4endl;
        G4double angle = abs(primary.GetDirection().angle(norm));
        if (angle > thetaCrit) //is it subject to total internal reflection? If this is true, then it does create QPs.
        {
            primary.SetStatus(true);
        }
        else
        {
            if (G4UniformRand() > filmAbsorption) //does it reflect due to transmission coefficient?
            {
                primary.SetStatus(false);
            }
            else //escapes back into crystal at this point
            {
                primary.ReflectDirection();
                reflectedEnergies.push_back(primary);
                return 0.;
            }
        }
    }
    // Phonon goes into superconductor and gets partitioned into
    // quasiparticles and new phonons
    G4double EDep = 0.;
    std::vector<Quasiparticle> qpEnergies;
    std::vector<Phonon> phonEnergies{ primary };
    while (qpEnergies.size() > 0 || phonEnergies.size() > 0) {

        //loops through all phonons and determines whether they create QPs and the energies
        // of the QPs they create.
        if (phonEnergies.size() > 0) {
            // Partition the phonons' energies into quasi-particles according to
            // a PDF defined in CalcQPEnergies().
            // NOTE: Both energy vectors mutate.
            EDep += CalcQPEnergies(phonEnergies, qpEnergies);
        }

        //record qp data
        for (Quasiparticle qp : qpEnergies)
        {
            qpData->append(std::to_string(qp.GetId()) + ",");
            qpData->append(std::to_string(1000.0 * qp.GetEnergy() / eV) + ",");
            qpData->append(std::to_string(qp.GetTime()) + ",");
        }

        //loops through all QPs and finds the energy(s) of the phonons the QPs create either
        //through scattering or recombination
        if (qpEnergies.size() > 0) {
            // Quasiparticles can also excite phonons.
            // NOTE: Both energy vectors mutate.

            EDep += CalcPhononEnergies(phonEnergies, qpEnergies);
        }

        //loops through all phonons and determines which ones escape out of crystal and which
        //ones will create QPs the next time this whole "while" loop is run.
        if (phonEnergies.size() > 0) {
            // Some phonons will escape back into the crystal.
            // NOTE: Both energy vectors mutate.
            CalcReflectedPhononEnergies(phonEnergies, reflectedEnergies);
        }
    }

    //after this while loop ends, there are no more QPs or phonons in the film that
    // were ultimately created from the incident phonon.

    // Sanity check -- Reflected + Absorbed should equal input
    /*G4double ERefl = std::accumulate(reflectedEnergies.begin(),
        reflectedEnergies.end(), 0.);*/

    G4double ERefl = 0.0;
    for (Phonon phonon : reflectedEnergies)
    {
        ERefl += phonon.GetEnergy();
    }
    //EDep = energy - ERefl;
    if (verboseLevel > 1) {
        G4cout << " Reflected " << ERefl << " (" << reflectedEnergies.size()
            << ")\n Absorbed " << EDep << G4endl;
    }

    if (fabs(energy - ERefl - EDep) / energy > 1e-5) {
        G4cerr << "WARNING Kaplan_QP missing " << (energy - ERefl - EDep) / eV
            << " eV" << G4endl;
    }

#ifdef G4CMP_DEBUG
    if (output.good()) {
        output << energy / eV << "," << EDep / eV << "," << ERefl / eV << ","
            << reflectedEnergies.size() << std::endl;
    }
#endif

    return EDep;
}


// Compute the probability of phonon reaching the film-crystal interface without breaking
// any Cooper pairs.

G4double Kaplan_QP::CalcEscapeProbability(Phonon phonon,
    G4double thicknessFrac) const {
    G4double energy = phonon.GetEnergy();
    if (verboseLevel > 1) {
        G4cout << "Kaplan_QP::CalcEscapeProbability E " << energy
            << " thickFrac " << thicknessFrac << G4endl;
    }

    // Compute energy-dependent mean free path for phonons in film
    if (gapEnergy <= 0.) return 1.;

    //calculate mfp of phonon given the velocity of sound in substrate, phonon lifetime,
    //its energy, and the phonon lifetime slope.
    G4double mfp = vSound * phononLifetime /
        (1. + phononLifetimeSlope * (energy / gapEnergy - 2.));

    /*if (verboseLevel > 2) {
        G4cout << " mfp " << mfp << " returning "
            << std::exp(-2. * thicknessFrac * filmThickness / mfp) << G4endl;
    }*/

    //this line assumes the phonons are traveling perfectly vertically, which is not true
    //return std::exp(-2. * thicknessFrac * filmThickness / mfp);

    //this line of code correctly (i think) accounts for slanted direction of phonons when 
    //calculating the distance the phonon must travel to get to the interface, and therefore,
    //the probability that it will create QPs.
    return std::exp(-2. * thicknessFrac * filmThickness / (mfp * fabs(phonon.GetDirection().z())));
}


// Model the phonons (phonEnergies) breaking Cooper pairs into quasiparticles
// (qpEnergies).

G4double
Kaplan_QP::CalcQPEnergies(std::vector<Phonon>& phonEnergies,
    std::vector<Quasiparticle>& qpEnergies) const {
    if (verboseLevel > 1) {
        G4cout << "Kaplan_QP::CalcQPEnergies QPcut " << lowQPLimit * gapEnergy
            << G4endl;
    }

    // Phonons above the bandgap give all of its energy to the qp pair it breaks.
    G4double EDep = 0.;
    std::vector<Phonon> newPhonEnergies;

    //loop through every phonon
    for (Phonon phonon : phonEnergies) {
        const G4double& E = phonon.GetEnergy();

        //if the phonon is subgap or its status is false, it does not create QPs.
        if (IsSubgap(phonon) || !phonon.GetStatus()) {
            if (verboseLevel > 2) G4cout << " Skipping phononE " << E << G4endl;
            newPhonEnergies.push_back(phonon);
            continue;
        }
        
        //energy of first QP created is calculated using QPEnergyRand
        Quasiparticle qp = QPEnergyRand(phonon);
        qp.SetId(qpId);
        //increment QP id
        qpId++;
        G4double qpE = qp.GetEnergy();
        if (verboseLevel > 2) G4cout << " phononE " << E << " qpE " << qpE << G4endl;

        //EDep += CalcQPAbsorption(qpE, newPhonEnergies, qpEnergies);
        //EDep += CalcQPAbsorption(Quasiparticle(E - qpE), newPhonEnergies, qpEnergies);
        qpEnergies.push_back(qp);

        //energy of second QP gets remainder of phonon energy (Eph - Eqp)
        Quasiparticle qp1(E - qp.GetEnergy(), qp.GetTime(), qpLifetime);
        qp1.SetId(qpId);
        //increment QP id
        qpId++;
        qpEnergies.push_back(qp1);
    }	// for (E: ...)

    if (verboseLevel > 1)
        G4cout << " replacing phonEnergies, returning EDep " << EDep << G4endl;

    //update list of phonons since some were killed when creating QPs.
    phonEnergies.swap(newPhonEnergies);
    return EDep;
}


// Model the quasiparticles (qpEnergies) emitting phonons (phonEnergies) in
// the superconductor.

G4double
Kaplan_QP::CalcPhononEnergies(std::vector<Phonon>& phonEnergies,
    std::vector<Quasiparticle>& qpEnergies) const {
    if (verboseLevel > 1) {
        G4cout << "Kaplan_QP::CalcPhononEnergies 2*gap " << 2. * gapEnergy
            << " QPcut " << lowQPLimit * gapEnergy << G4endl;
    }

    G4double EDep = 0.;
    G4double lifetime;
    std::vector<Quasiparticle> newQPEnergies;

    G4double E;
    //loop through every QP.
    for (Quasiparticle qp : qpEnergies) {
        //qp.IncrementTime();
        //const G4double& E = qp.GetEnergy();

        //if the QP has 0 energy because it recombined, it will be deleted and not create any phonons.
        if (qp.GetIsZeroEnergy())
        {
            //qpIndex++;
            continue;
        }
        E = qp.GetEnergy();

        //if QP is above 3 delta (default lowQPLimit value is 3), then it will emit a phonon by scattering.
        //the phonon energy is calculated using PhononEnergyRand.

        //notice the increment time. Currently it is set to 0 since the scattering lifetime assumed to be
        //negligible. However, this could be an erroneous assumption, so a different value for the scattering
        //lifetime. The lifetime could also be calculated in terms of the QP's energy or other variables.
        if (E >= lowQPLimit * gapEnergy)
        {
            qp.IncrementTime(0.0);
            //qp.IncrementTime(qpLifetime / (3 * pow((E / eV) / gapEnergy, 3.0)));
            //qp.IncrementTime(0.2 * qpLifetime);
            //qp.IncrementTime(qpLifetime);
            Phonon phon = PhononEnergyRand(qp);
            phonEnergies.push_back(phon);
            qp.SetEnergy(E - phon.GetEnergy());
            newQPEnergies.push_back(qp);
        }
        //if its energy is above delta but not above 3 delta, it will emit the remainder of its energy above
        // the gap into a phonon (E - gapEnergy). Again, the scattering lifetime set to 0.
        else if (E > gapEnergy)
            //else if (!qp.GetAtGap())
        {
            qp.IncrementTime(0.0);
            //qp.IncrementTime(qpLifetime);
            Phonon phon = Phonon(E - gapEnergy, qp.GetTime(), zmin);
            phonEnergies.push_back(phon);
            qp.SetEnergy(gapEnergy);
            qp.SetAtGap(true);
            newQPEnergies.push_back(qp);
        }
        else //qp is at delta
        {
            //lifetime = qpLifetime;

            //random recombination lifetime is chosen from exponential decay probablity function with
            // average of qpLifetime.
            lifetime = -log(G4UniformRand()) * qpLifetime;

            //if even id, the QP recombines with quiescent QP in film.
            // therefore, it creates a 2 delta phonon and its energy is
            // set to 0.
            //The delta energy gained from the quiescent QP in the film (and the time)
            // is also logged so the total energy in the analysis stays constant.
            if (qp.GetId() % 2 == 0)
            {
                qp.IncrementTime(lifetime);
                Phonon phon = Phonon(2 * gapEnergy, qp.GetTime(), zmin);

                qp.SetEnergy(0.0);
                qp.SetIsZeroEnergy(true);
                phonEnergies.push_back(phon);
                newQPEnergies.push_back(qp);

                energyGainData->append(std::to_string(qp.GetTime()) + ",");
                energyGainData->append(std::to_string(1) + ",");

                EDep -= gapEnergy;
            }
            //if odd id, the QP is lost into film and becomes a quiescent QP. Again,
            // the energy loss and the time when the energy is lost to the film is recorded.
            // Because exactly half of delta QPs create 2 delta phonons and other half are
            // deleted, energy in the simulation and the QP density in the film is conserved.
            else
            {
                qp.IncrementTime(lifetime);
                //qp.SetEnergy(0.0);
                qp.SetIsZeroEnergy(true);
                newQPEnergies.push_back(qp);

                energyGainData->append(std::to_string(qp.GetTime()) + ",");
                energyGainData->append(std::to_string(-1) + ",");

                EDep += gapEnergy;
            }
        }
    }
    if (verboseLevel > 1)
        G4cout << " replacing qpEnergies, returning EDep " << EDep << G4endl;

    //update list of QPs because some were killed due to having 0 energy.
    qpEnergies.swap(newQPEnergies);
    return EDep;
}


// Calculate energies of phonon tracks that have reentered the crystal.

void Kaplan_QP::
CalcReflectedPhononEnergies(std::vector<Phonon>& phonEnergies,
    std::vector<Phonon>& reflectedEnergies) const {
    if (verboseLevel > 1)
        G4cout << "Kaplan_QP::CalcReflectedPhononEnergies " << G4endl;

    std::vector<Phonon> newPhonEnergies;

    //loop through all phonons
    for (Phonon phonon : phonEnergies) {
        //energy of phonon
        const G4double& E = phonon.GetEnergy();
        if (verboseLevel > 2) G4cout << " phononE " << E << G4endl;

        // Phonons below the bandgap are unconditionally reflected
        if (IsSubgap(phonon)) {
            if (verboseLevel > 2) G4cout << " phononE got reflected" << G4endl;
            phonon.ReflectDirection();
            reflectedEnergies.push_back(phonon);
            continue;
        }

        // frac = 1.5 for phonons headed away from the subst. 0.5 for toward.
        // This assumes that, on average, the phonons are spawned at the center
        // of the superconductor. This could be a false assumption, almost certianly
        // is for Nb since phonon mfp is much less than film thickness in current detector.
        G4double frac;
        if (phonon.GetDirection().z() > 0.)//heading up = 1.5 * film thickness vertically to get to interface
        {
            frac = 1.5;
        }
        //heading downwards = 0.5 film thickness vertically since all phonons created through recombination start
        //in the center of the crystal.
        else
        {
            frac = 0.5;
        }
        //if the phonon did not create QPs the previous loop, it must have gotten to the interface and been reflected
        //there to still be in crystal and not create QPs (status = false). Thus, it must pass through twice the film
        //thickness in the vericle direction to get to the interface again.
        if (!phonon.GetStatus())
        {
            frac = 2.0;
        }

        //does the phonon create QPs?
        //if this "if statement" executes, then the phonon does not create QPs before getting to the interface
        if (G4UniformRand() < CalcEscapeProbability(phonon, frac)) {
            G4double angle = abs(phonon.GetDirection().angle(norm));
            //does it get reflected due to tir? (yes-> will inevitably create QPs since the top and bottom of the film
            // are assumed to be perfectly smooth and parallel. Therefore, the phonon will bounce around in inside of film
            // until it creates QPs. This could be a false assumption.
            if (angle > thetaCrit)
            {
                phonon.SetStatus(true);
                newPhonEnergies.push_back(phonon);
            }
            else
            {
                //does it get reflected due to film absorption (yes->stays in crystal, but does not create QPs next loop,
                // no->escapes).
                if (G4UniformRand() > filmAbsorption)
                {
                    phonon.SetStatus(false);
                    newPhonEnergies.push_back(phonon);
                }
                else //phonon escapes into crystal.
                {
                    if (phonon.GetDirection().z() > 0.)
                    {
                        phonon.ReflectDirection();
                    }
                    reflectedEnergies.push_back(phonon);
                }
            }
        }
        else {
            //if (verboseLevel > 2) G4cout << " phononE stays in film" << G4endl;

            //if the phonon does create QPs, it status is set to true.
            phonon.SetStatus(true);
            newPhonEnergies.push_back(phonon);
        }
    }	// for (E: ...)

    phonEnergies.swap(newPhonEnergies);
}


// Compute probability of absorbing phonon below 2 * gap energy.
// set to 0%, subgap phonons will always escape with no energy deposition
// or QP creation.

G4double
Kaplan_QP::CalcSubgapAbsorption(Phonon phonon,
    std::vector<Phonon>& keepEnergies) const {
    G4double energy = phonon.GetEnergy();
    if (G4UniformRand() < subgapAbsorption) {
        if (verboseLevel > 2)
            G4cout << " Deposit phonon " << energy << " as heat" << G4endl;

        return energy;
    }
    else {
        if (verboseLevel > 2)
            G4cout << " Record phonon " << energy << " for processing" << G4endl;

        keepEnergies.push_back(phonon);
        return 0.;
    }
}


// Handle absorption of quasiparticle energies below Cooper-pair breaking
// If qpEnergy < 3*Delta, radiate a phonon, absorb bandgap minimum


//Not used anymore
G4double
Kaplan_QP::CalcQPAbsorption(Quasiparticle qp,
    std::vector<Phonon>& phonEnergies,
    std::vector<Quasiparticle>& qpEnergies) const {
    G4double EDep = 0.;		// Energy lost by this QP into the film

    G4double qpE = qp.GetEnergy();
    if (qpE >= lowQPLimit * gapEnergy) {
        if (verboseLevel > 2) G4cout << " Storing qpE in qpEnergies" << G4endl;
        qpEnergies.push_back(qp);
    }
    else if (qpE > gapEnergy) {
        if (verboseLevel > 2) G4cout << " Reducing qpE to gapEnergy" << G4endl;
        EDep += CalcSubgapAbsorption(Phonon(qpE - gapEnergy), phonEnergies);
        EDep += gapEnergy;
    }
    else {
        EDep += qpE;
    }

    return EDep;
}

// Compute quasiparticle energy distribution from broken Cooper pair.

Quasiparticle Kaplan_QP::QPEnergyRand(Phonon phonon) const {
    // PDF is not integrable, so we can't do an inverse transform sampling.
    // Instead, we'll do a rejection method.
    //
    // PDF(E') = (E'*(Energy - E') + gapEnergy*gapEnergy)
    //           /
    //           sqrt((E'*E' - gapEnergy*gapEnergy) *
    //                ((Energy - E')*(Energy - E') - gapEnergy*gapEnergy));
    // The shape of the PDF is like a U, so the max values are at the endpoints:
    // E' = gapEnergy and E' = Energy - gapEnergy

    // Add buffer so first/last bins don't give zero denominator in pdfSum

    G4double Energy = phonon.GetEnergy();

    const G4double BUFF = 1000.;
    G4double xmin = gapEnergy + (Energy - 2. * gapEnergy) / BUFF;
    G4double xmax = gapEnergy + (Energy - 2. * gapEnergy) * (BUFF - 1.) / BUFF;
    G4double ymax = QPEnergyPDF(Energy, xmin);

    G4double xtest = 0., ytest = ymax;
    do {
        ytest = G4UniformRand() * ymax;
        xtest = G4UniformRand() * (xmax - xmin) + xmin;
    } while (ytest > QPEnergyPDF(Energy, xtest));

    return Quasiparticle(xtest, phonon.GetTime(), qpLifetime);
}

G4double Kaplan_QP::QPEnergyPDF(G4double E, G4double x) const {
    const G4double gapsq = gapEnergy * gapEnergy;
    return ((x * (E - x) + gapsq) / sqrt((x * x - gapsq) * ((E - x) * (E - x) - gapsq)));
}


// Compute phonon energy distribution from quasiparticle in superconductor.

Phonon Kaplan_QP::PhononEnergyRand(Quasiparticle qp) const {
    // PDF is not integrable, so we can't do an inverse transform sampling.
    // Instead, we'll do a rejection method.
    //
    // PDF(E') = (E'*(Energy-E')*(Energy-E') * (E'-gapEnergy*gapEnergy/Energy))
    //           /
    //           sqrt((E'*E' - gapEnergy*gapEnergy);

    // Add buffer so first bin doesn't give zero denominator in pdfSum

    G4double Energy = qp.GetEnergy();

    const G4double BUFF = 1000.;
    G4double xmin = gapEnergy + gapEnergy / BUFF;
    G4double xmax = Energy;
    G4double ymax = PhononEnergyPDF(Energy, xmin);

    G4double xtest = 0., ytest = ymax;
    do {
        ytest = G4UniformRand() * ymax;
        xtest = G4UniformRand() * (xmax - xmin) + xmin;
    } while (ytest > PhononEnergyPDF(Energy, xtest));

    return Phonon(Energy - xtest, qp.GetTime(), zmin);
}

G4double Kaplan_QP::PhononEnergyPDF(G4double E, G4double x) const {
    const G4double gapsq = gapEnergy * gapEnergy;
    return (x * (E - x) * (E - x) * (x - gapsq / E) / sqrt(x * x - gapsq));
}

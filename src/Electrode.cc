#include "Electrode.hh"
#include "G4AffineTransform.hh"
#include "G4FieldManager.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4ParticleChange.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "Randomize.hh"

//#ifdef CDMS_USE_G4CMP
#include "G4CMPGeometryUtils.hh"
#include "G4CMPPhononTrackInfo.hh"
#include "G4CMPSecondaryUtils.hh"
#include "G4CMPSurfaceProperty.hh"
#include "G4CMPTrackUtils.hh"
#include "G4CMPUtils.hh"
#include "G4LatticeManager.hh"
#include "G4LatticePhysical.hh"
#include "G4CMPKaplanQP.hh"
//#endif

#include "G4SystemOfUnits.hh"

#include <stdio.h>
#include <stdlib.h>


#include "Phonon.hh"
// For phonons, electrode hits only happen at crystal-aluminum boundary
// See CDMSZipConstruction::AddLatticeSurfaces()

Electrode::Electrode(G4MaterialPropertiesTable* prop)
{
    //initialize material properties table, verboseLevel, Kaplan_QP,
    //and strings to store data from PSKQP and reflected phonons
    matTab = prop;
    //leave at 0 unless you want this class to output information.
    //in general, higher verbose = more information outputed from that class.
    verboseLevel = 0;
    //does all QP and phonon physics inside film.
    kaplanQP = new Kaplan_QP(matTab);
    //contains data about each QP created in film that this electrode is
    //attatched to.
    elecData = new std::string();

    //not super necessary for analysis, could remove all places where it is 
    //used if you do not want simulation to output ReflData.txt.
    //could also comment out where it is written to file in Stepping_Action.cc
    reflData = new std::string();
}

Electrode::~Electrode()
{
    delete kaplanQP;
}

G4bool Electrode::IsNearElectrode(const G4Step& /*step*/) const {
    //filmAbsorption - prob. that a phonon which hit the interface from crystal side
    // gets transmitted into film.
    //if this function returns true, the phonon goes into film.
    //if false, it gets reflected.
    return G4UniformRand() < GetMaterialProperty("filmAbsorption");
    //return false;
    //return true;
}

// Phonon gets killed upon breaking a Cooper pair. New phonons may
// be created and emitted back into the crystal

void Electrode::
AbsorbAtElectrode(const G4Track& track, const G4Step& step, G4ParticleChange& particleChange) const
{
    // Incident phonon is killed, all energy transferred out
    particleChange.ProposeTrackStatus(fStopAndKill);
    particleChange.ProposeEnergy(0.);

    // Transfer phonon energy into superconducting film
    //energy of initial phonon
    G4double phonEnergy = GetKineticEnergy(track);
    //wave vector of initial photon
    G4ThreeVector k = GetLocalWaveVector(track);
    //normalize
    G4ThreeVector kunit = k.unit();
    //if the phonon goes in and creates no QPs, it should be reflected out of the crystal specularly
    //since Snell's law applies to both sides of the interface.
    //doSpecular keeps track of whether it should be reflected specularly or not. Set to true at first,
    //and changed to false if QPs are created.
    G4bool doSpecular = true;

    //used to store data about all the phonons which escape from the film and are child phonons
    //of the original incident phonon.
    std::vector<Phonon> reflPhonons;
    //get the time when the phonon went into the film.
    G4double time = track.GetGlobalTime();

    //primary phonon - abstract representation of the incident geant4 particle phonon which went into the film.
    //only keeps track of its energy and direction.
    //the time of the primary phonon is set to 0 so that the times of all phonons/QPS created in the film are 
    //equal to the difference in time between the phonon hitting the film and them being created.
    //when the secondaries are created the original global time of the phonon is added to their times to give
    //them the correct global times. The reason this is done is to prevent large time values from being assigned
    //to QPs since the true global time of the phonon could be up to millions of nanoseconds.
    Phonon primary = Phonon(phonEnergy, 0);

    //k unit will be the direction of the the primary phonon.
    //*****NOTE: Snell's law is not accounted for here, should be implemented by transforming kunit vector
    //and using that new vector as the primary's direction.*****
    primary.SetDirection(kunit);
    //give kaplanQP access to the primary phonon so it can carry out physics inside simulation
    G4double EDep = kaplanQP->AbsorbPhonon(primary, reflPhonons);
    particleChange.ProposeNonIonizingEnergyDeposit(EDep);


    //get name of volume that phonon went into
    //need to use if statement to check if volume exists to avoid possibility of segmentation fault.
    std::string volName;
    if (step.GetPostStepPoint()->GetTouchableHandle()->GetVolume())
    {
        volName = step.GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
    }
    //Write QP data
    //incudes the time of the phonon hit, the volume, and all the QP data and QP gain data.
    G4double minEnergy = 100;
    if (kaplanQP->GetQPDataLength() > 1)
    {
        //elecData += std::to_string(time) + "," + vol_name + "," + std::to_string(sqp->qpLifetime) + sqp->data + "\n";
        elecData->append(std::to_string(time) + "," + volName + "," + *kaplanQP->GetQPData() + *kaplanQP->GetEnergyGainData() + "\n");
        //phonon created QPs since the length of the QP data is greater than 1, so doSpecular is set to false
        doSpecular = false;
    }

    //get normal of surface that phonon went into
    //should always be (0, 0, 1)
    CLHEP::Hep3Vector surfNorm = G4CMP::GetSurfaceNormal(step);
    //CLHEP::Hep3Vector surfNorm = G4ThreeVector(0, 0, 1);
    
    //how many secondaries should the initial phonon have by creating QPs which recombine to create secondary phonons
    //that escape?

    //all secondaries are phonons
    particleChange.SetNumberOfSecondaries(reflPhonons.size());

    //loops through all relflected phonons calculated in kaplanQP.
    //for each of these, it gets its energy and magnitude of k vector
    //the polarization of the outgoing phonons is chosen randomly. This might be a false assumption.
    //could include polarization of phonons as part of physics in Kaplan_QP, not currently accounted for.

    //Reflection models for outgoing phonons
    
    //1. all diffuse - the directions of every phonon to come from the film is generated randomly using
    // the DMM model (calculated with "LambertReflection"). This is incorrect because phonons that go
    // into the film and create no QPs should be refleced according to specular reflection (AMM).
    // also, the directions of the abstract PSPhonons inside the film are not used to calculate the directions
    // of the actual geant4 phonons exiting the film. Instead, the directions are generated completely randomly
    // with LambertReflection.

    //2. partial specular - same as 1. but initial phonon is reflected specularly (and therefore correctly).

    //1 vs 2 has practically no effect on output.

    //*****NOTE: what should really happens is that the directions of the PSPhonons are refracted according
    //to Snell's law and those new directions are used as the reflected phonon directions.*****
    G4double ERefl = 0.0;
    for (Phonon reflPhonon : reflPhonons) {
        G4double E = reflPhonon.GetEnergy();
        //comments put here by Mike Kelsey
        // TODO: Map E to K properly.
        // TODO: Is this the correct way to assign k magnitudes
        G4double kmag = k.mag() * E / phonEnergy;
        G4int pol = ChoosePhononPolarization();

        //all diffuse - what is used by default in supersim

        G4ThreeVector reflectedKDir = k.unit();
        do {
            reflectedKDir = G4CMP::LambertReflection(surfNorm);
        } while (!G4CMP::PhononVelocityIsInward(theLattice, pol,
            kmag * reflectedKDir, surfNorm));

        //partial specular (commented out for now)
        
        /*G4ThreeVector reflectedKDir = k.unit();
        if (doSpecular)
        {
            if(!G4CMP::PhononVelocityIsInward(theLattice, pol,
                kmag * reflectedKDir, surfNorm))
                reflectedKDir = SpecularReflection(kunit, surfNorm);
        }
        else
        {
            do {
                reflectedKDir = G4CMP::LambertReflection(surfNorm);
            } while (!G4CMP::PhononVelocityIsInward(theLattice, pol,
                kmag * reflectedKDir, surfNorm));
        }*/

        //set direction of reflected phonon
        reflPhonon.SetDirection(reflectedKDir);

        //ensures phonon is heading away from film, back into substrate so it will not be immediately reabsorbed.
        if (!G4CMP::PhononVelocityIsInward(theLattice, pol, kmag * reflPhonon.GetDirection(), surfNorm))
        {
            reflPhonon.SetDirection(SpecularReflection(reflPhonon.GetDirection(), surfNorm));
        }
        
        //record direction of every phonon and whether or not it is heading away from film.
        reflData->append(std::to_string(reflPhonon.GetDirection().x()) + ",");
        reflData->append(std::to_string(reflPhonon.GetDirection().y()) + ",");
        reflData->append(std::to_string(reflPhonon.GetDirection().z()) + ",");
        reflData->append(std::to_string(G4CMP::PhononVelocityIsInward(theLattice, pol, kmag * reflPhonon.GetDirection(), surfNorm)) + ",");
        reflData->append(std::to_string(reflPhonon.GetEnergy() / eV) + ",");
        reflData->append(std::to_string(reflPhonon.GetTime()) + "\n");


        //create actual geant4 secondary reflected phonon
        particleChange.AddSecondary(G4CMP::CreatePhonon(GetCurrentTouchable(),
            pol, kmag * reflPhonon.GetDirection(),
            E, time + reflPhonon.GetTime(),
            track.GetPosition()));
    }

    // Sanity check: secondaries' energy should equal assigned E
    if (verboseLevel > 1) {
        G4double Esum = 0.;
        for (G4int i = 0; i < particleChange.GetNumberOfSecondaries(); i++) {
            Esum += particleChange.GetSecondary(i)->GetKineticEnergy();
            G4cout << " secondary " << i << " GetKineticEnergy() "
                << particleChange.GetSecondary(i)->GetKineticEnergy()
                << G4endl;
        }

        if (fabs(phonEnergy - EDep - Esum) > 1e-6) {
            G4cerr << "ERROR: Energy non-conservation: Ekin-EDep " << phonEnergy - EDep
                << " vs. sum of secondaries " << Esum << G4endl;
        }
    }
    
}

//does specular reflection for phonon hitting a surface with normal surfNorm
G4ThreeVector Electrode::SpecularReflection(G4ThreeVector primaryPhononDir, G4ThreeVector surfNorm) const
{
    return primaryPhononDir - 2 * (primaryPhononDir.dot(surfNorm)) * surfNorm;
}
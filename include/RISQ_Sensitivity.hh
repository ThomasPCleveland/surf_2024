/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/

#ifndef RISQTutorialSensitivity_h
#define RISQTutorialSensitivity_h 1

#include "G4CMPElectrodeSensitivity.hh"

class RISQ_Sensitivity final : public G4CMPElectrodeSensitivity {
public:
  RISQ_Sensitivity(G4String name);
  virtual ~RISQ_Sensitivity();
  // No copies
  RISQ_Sensitivity(const RISQ_Sensitivity&) = delete;
  RISQ_Sensitivity& operator=(const RISQ_Sensitivity&) = delete;
  /* Move is disabled for now because old versions of GCC can't move ofstream
  // Move OK
  RISQ_Sensitivity(RISQ_Sensitivity&&);
  RISQ_Sensitivity& operator=(RISQ_Sensitivity&&);
  */
  RISQ_Sensitivity(RISQ_Sensitivity&&) = delete;
  RISQ_Sensitivity& operator=(RISQ_Sensitivity&&) = delete;

  virtual void EndOfEvent(G4HCofThisEvent*);

  void SetHitOutputFile(const G4String& fn);
  void SetPrimaryOutputFile(const G4String& fn);

protected:
  virtual G4bool IsHit(const G4Step*, const G4TouchableHistory*) const;

private:
  std::ofstream primaryOutput;
  std::ofstream hitOutput;
  G4String primaryFileName;
  G4String hitFileName;
  
};

#endif

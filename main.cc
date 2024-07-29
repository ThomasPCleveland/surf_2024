
#include "Custom_Physics.hh"
#include "Physics_List.hh"
#include "Config_Manager.hh"
#include "Detector.hh"
#include "Single_KID.hh"
#include "RISQ_Detector.hh"
#include "Action_Initialization.hh"

#include "G4RunManager.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4CMPPhysicsList.hh"
#include "G4CMPPhysics.hh"
#include "G4CMPConfigManager.hh"
#include "FTFP_BERT.hh"
#include "G4CMPPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4ChargeExchangePhysics.hh"
#include "G4IonPhysicsXS.hh"
#include "G4IonElasticPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

// #include "CDMSHadronicBiasing.hh"

int main(int argc, char **argv)
{
  {
    std::remove("volume_hits.log");
    std::ofstream hitsLog;
    hitsLog.open("volume_hits.log", std::ios_base::app);
    hitsLog << "run event track particle energy deposit x y z t_i t_f volume\n";
    hitsLog.close();
  }

  G4RunManager *runManager = new G4RunManager;
  // Detector *detector = new Detector();
  // Single_KID *detector = new Single_KID();
  RISQ_Detector *detector = new RISQ_Detector();
  runManager->SetUserInitialization(detector);
  Physics_List *physics = new Physics_List();
  runManager->SetUserInitialization(physics);
  runManager->SetUserInitialization(new Action_Initialization);

  // prepare commands
  G4CMPConfigManager::Instance();
  Config_Manager::Instance();

  // setup UX
  G4VisManager *visManager = new G4VisExecutive;
  visManager->Initialize();
  G4UImanager *UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/control/macroPath ../macros");
  // UImanager->ApplyCommand("/control/execute prep"); // doesn't work, needs a delay?

  if (argc == 1) // if interactive
  {
    G4UIExecutive *ui = new G4UIExecutive(argc, argv);
    ui->SessionStart();
    delete ui;
  }
  else // otherwise console
  {
    G4String command_headless = G4String("/control/execute ");
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command_headless + fileName);
  }

  delete visManager;
  delete runManager;
  return 0;
}

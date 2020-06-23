//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file exampleDRAGON.cc
/// \brief Main program of the DRAGON example


#include "DRAGONDetectorConstruction.hh"
#include "DRAGONActionInitialization.hh"
#include "QGSP_BERT.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "InelasticPhysics.hh"
//#include "BinaryReactionPhysics.hh"
#include "json/json.h"

#include "G4ProcessManager.hh"
#include "G4ionIonisation.hh"
#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"
#include "G4VProcess.hh"
#include "G4HadronicProcessStore.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "Randomize.hh"
#include "time.h"

int main(int argc,char** argv)
{
    G4HadronicProcessStore::Instance()->SetVerbose(0);
    G4double recoilMass = 4.;
    G4double recoilCharge = 2.;

    //choose the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
    //set random seed with system time
//  G4long seed = time(NULL);
//  if(argc<2) {seed+=473879*atoi(argv[2]); G4cout << "SNORLAX!!" << argv[2] << G4endl;}
//  CLHEP::HepRandom::setTheSeed(seed);

    //CREATE AND READ JSON CONFIG
    Json::Value config;
    std::string configFileName = argv[1];
    std::ifstream configStream(configFileName.c_str());
    configStream >> config;
    configStream.close();

    //Parse JSON
    G4double gasPressure = config["gasPressure"].asDouble(); //Torr
    G4double gasTemperature = config["gasTemperature"].asDouble(); //K
    G4int processNumber = config["processNumber"].asInt();
    std::map<std::string,double> reactionParams;
    G4String macroName = config["macroName"].asString();
    G4bool isInteractive = config["interactive"].asBool();
    reactionParams["qValue"] = config["qValue"].asDouble();
    reactionParams["lightProductCharge"] = config["lightProduct"][0].asDouble();
    reactionParams["lightProductMass"] = config["lightProduct"][1].asDouble();
    reactionParams["heavyProductCharge"] = config["heavyProduct"][0].asDouble();
    reactionParams["heavyProductMass"] = config["heavyProduct"][1].asDouble();
    reactionParams["targetCharge"] = config["target"][0].asDouble();
    reactionParams["targetMass"] = config["target"][1].asDouble();
    reactionParams["cmEnergy"] = config["cmEnergy"].asDouble();

    DRAGONDetectorConstruction* detector = new DRAGONDetectorConstruction();
    detector->SetGasPressure(gasPressure);

    G4long seed = time(NULL);
    seed+=473879*processNumber;
    CLHEP::HepRandom::setTheSeed(seed);

#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
#else
    G4RunManager* runManager = new G4RunManager;
#endif

    // Mandatory user initialization classes

    runManager->SetUserInitialization(new DRAGONDetectorConstruction);

    // //Changes: new lines for the decay//
    //    G4ProcessManager* pManager = G4GenericIon::GenericIon()->GetProcessManager();
    //    //EM Processes for generic ions
    //    pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
    //    G4cout << "hello 3" << G4endl;
    //    pManager->AddProcess(new G4ionIonisation,-1,2,2);
    //    //6He Decay Physics
    //    G4String the6HennProcessName = "GeneralPhase3_2nDecay";
    //    GeneralPhase3_2nDecay* the6HeProcess = new GeneralPhase3_2nDecay(the6HennProcessName);
    //    G4cout << "hello 3" << G4endl;
    //    the6HeProcess->SetParentNucleus(6,2);                        //6He
    //    G4cout << "hello 3.5" << G4endl;
    //    pManager->AddDiscreteProcess(the6HeProcess);
    //    G4cout << "hello 4" << G4endl;

    G4VModularPhysicsList* physicsList = new QGSP_BERT;

    //G4RadioactiveDecayPhysics* radDecayPhysics = new G4RadioactiveDecayPhysics;
    //pysicsList->RegisterPhysics(radDecayPhysics);
/*    ElasticScatteringPhysics* elasticScatteringPhysics = new ElasticScatteringPhysics;
    elasticScatteringPhysics->SetCMScatteringEnergy(cmEnergy);
    elasticScatteringPhysics->SetRecoilMassCharge(recoilMass,recoilCharge);
    physicsList->RegisterPhysics(elasticScatteringPhysics);
    runManager->SetUserInitialization(physicsList);
  */


//  Physics Lists -- Comment out first three lines to disable reactions.
    InelasticPhysics* reactionPhysics = new InelasticPhysics;
    reactionPhysics->SetReactionParams(reactionParams);
    physicsList->RegisterPhysics(reactionPhysics);
    runManager->SetUserInitialization(physicsList);         //Need this even if the physics lists are commented out.


    // User action initialization

    DRAGONActionInitialization* actionInitialization = new DRAGONActionInitialization;
    actionInitialization->SetRecoilMassCharge(recoilMass,recoilCharge);
    runManager->SetUserInitialization(actionInitialization);


    // Initialize Geant4 kernel
    runManager->Initialize();


#ifdef G4VIS_USE
    // Visualization manager construction
    G4VisManager* visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();
#endif

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (!isInteractive) {
        // execute an argument macro file if exist
        G4String command = "/control/execute ";
        G4String fileName = macroName;
        UImanager->ApplyCommand(command+fileName);
    }
    else {
        // start interactive session
#ifdef G4UI_USE
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
	UImanager->ApplyCommand("/control/execute init_vis.mac");
#else
	UImanager->ApplyCommand("/control/execute init.mac");
#endif
	if (ui->IsGUI())
	  UImanager->ApplyCommand("/control/execute gui.mac");
        ui->SessionStart();
        delete ui;
#endif
    }

    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !

#ifdef G4VIS_USE
    delete visManager;
#endif
    delete runManager;

    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

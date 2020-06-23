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
/// \file DRAGONActionInitialization.cc
/// \brief Implementation of the DRAGONActionInitialization class

#include "DRAGONActionInitialization.hh"
#include "DRAGONPrimaryGeneratorAction.hh"
#include "DRAGONEventAction.hh"
#include "DRAGONRunAction.hh"
#include "DRAGONTrackingAction.hh"

DRAGONActionInitialization::DRAGONActionInitialization() :
        G4VUserActionInitialization() {

}

DRAGONActionInitialization::~DRAGONActionInitialization() {
}

void DRAGONActionInitialization::BuildForMaster() const {
  SetUserAction(new DRAGONRunAction);
}

void DRAGONActionInitialization::Build() const {
  SetUserAction(new DRAGONPrimaryGeneratorAction);
  SetUserAction(new DRAGONRunAction);
//  SetUserAction(new DRAGONEventAction);
  DRAGONTrackingAction *trackingAction = new DRAGONTrackingAction;
  trackingAction->SetRecoilMassCharge(fRecoilMass,fRecoilCharge);
  SetUserAction(trackingAction);
  DRAGONEventAction* eventAction = new DRAGONEventAction(new DRAGONRunAction);
  SetUserAction(eventAction);
}

void DRAGONActionInitialization::SetRecoilMassCharge(G4double mass, G4double charge) {
  fRecoilMass = mass;
  fRecoilCharge = charge;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

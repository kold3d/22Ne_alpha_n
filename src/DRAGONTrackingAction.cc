//
// Created by Sriteja Upadhyayula on 2020-06-18.
//

#include "DRAGONTrackingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"
//#include "DRAGONTrackingInformation.hh"
#include "G4Step.hh"
#include <vector>

DRAGONTrackingAction::DRAGONTrackingAction() {
    fName = "Inelastic";
}

DRAGONTrackingAction::~DRAGONTrackingAction() {
}

void DRAGONTrackingAction::PreUserTrackingAction(const G4Track* track) {
    const G4VProcess* creatorProcess = track->GetCreatorProcess();
    if(!creatorProcess) return;
    if(creatorProcess->GetProcessName()!=fName) return;
    if(track->GetParticleDefinition()->GetAtomicNumber() != G4int(fRecoilCharge) ||
       track->GetParticleDefinition()->GetAtomicMass() != G4int(fRecoilMass)) return;
    G4ThreeVector VertexPosition = track->GetPosition();
    std::vector<double> pos_double(3);
    pos_double[0]=VertexPosition.x()/mm;
    pos_double[1]=VertexPosition.y()/mm;
    pos_double[2]=VertexPosition.z()/mm;

//    DRAGONTrackingInformation* info = (DRAGONTrackingInformation*) track->GetUserInformation();
//
//    DRAGONAnalysis* analysis = DRAGONAnalysis::Instance();
    // analysis->SetCMEnergy(info->GetCMEnergy());
    // analysis->SetCMTheta(info->GetCMAlphaTheta());
    // analysis->SetLabTheta(info->GetLabAlphaTheta());
    // analysis->SetCMPhi(info->GetCMAlphaPhi());
//    analysis->SetVertexPosition(pos_double);
}

void DRAGONTrackingAction::SetRecoilMassCharge(G4double mass, G4double charge) {
    fRecoilMass = mass;
    fRecoilCharge = charge;
}

//
// Created by Sriteja Upadhyayula on 2020-06-18.
//

#ifndef DRAGONTrackingAction_h
#define DRAGONTrackingAction_h

#include "G4String.hh"
#include "G4UserTrackingAction.hh"
#include "G4ThreeVector.hh"

class DRAGONTrackingAction : public G4UserTrackingAction {
public:
    DRAGONTrackingAction();
    ~DRAGONTrackingAction();

    void PreUserTrackingAction(const G4Track*);

    void SetRecoilMassCharge(G4double,G4double);

private:
    G4String fName;
    G4double fRecoilMass;
    G4double fRecoilCharge;
};

#endif

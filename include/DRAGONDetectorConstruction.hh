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
/// \file DRAGONDetectorConstruction.hh
/// \brief Definition of the DRAGONDetectorConstruction class

#ifndef DRAGONDetectorConstruction_h
#define DRAGONDetectorConstruction_h

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4GenericMessenger.hh"

class G4Material;

class DRAGONDetectorConstruction : public G4VUserDetectorConstruction {

public:
    DRAGONDetectorConstruction();
    virtual ~DRAGONDetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
//    virtual void ConstructSDandField();

    void SetGasPressure(G4double pressure){
        fPressureInTorr = pressure;
        G4cout << "Set: Gas pressure: " << fPressureInTorr << G4endl;
    };
    void SetTemperature(G4double);
    G4Material* GetGasMaterial(){return fGasMaterial;};
    G4LogicalVolume* GetWorldVolume() const {return fWorldLogical;}
private:
    double fPressureInTorr;
    double fTemperature;

    G4double fDriftGap;

    G4LogicalVolume* fVacuumLogical;
    G4LogicalVolume* fWorldLogical;

    G4Material* fGasMaterial;

    void ConstructMaterials();
    void SetAttributes();
};

#endif


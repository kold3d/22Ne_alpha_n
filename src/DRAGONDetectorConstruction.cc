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
/// \file DRAGONDetectorConstruction.cc
/// \brief Implementation of the DRAGONDetectorConstruction class

#include "DRAGONDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"

DRAGONDetectorConstruction::DRAGONDetectorConstruction() :
        G4VUserDetectorConstruction(),
        fPressureInTorr(5),
        fTemperature(293.15),
        fDriftGap(6.0) {
}

DRAGONDetectorConstruction::~DRAGONDetectorConstruction() {
    G4cout << "Pressure: " << fPressureInTorr << G4endl;
}

G4VPhysicalVolume* DRAGONDetectorConstruction::Construct() {
//    ConstructMaterials();

    G4NistManager* nistMgr = G4NistManager::Instance();
    //Load pointers to materials
    G4Material* he = nistMgr->FindOrBuildMaterial("G4_He");
    G4Material* co2 = nistMgr->FindOrBuildMaterial("G4_CARBON_DIOXIDE");

    // Specific Gas constant for Helium
    double gasConstant = 2077.1*joule*kelvin/kg;
    //Gas mixtures
    double volFracHe = 1.00;
    double volFracCO2 = 0.00;

    double molarMassHe = 4.0026020;
    double molarMassCO2 = 44.01;

    double massFracHe = volFracHe*molarMassHe/(volFracHe*molarMassHe+volFracCO2*molarMassCO2);
    double massFracCO2 = volFracCO2*molarMassCO2/(volFracHe*molarMassHe+volFracCO2*molarMassCO2);
//    double gasDensity = 1.603e-05*fPressureInTorr*(volFracHe*molarMassHe+volFracCO2*molarMassCO2)/fTemperature*g/cm3;
//    double gasDensity = 2.20865e-06*fPressureInTorr*(volFracHe*molarMassHe+volFracCO2*molarMassCO2)/fTemperature*g/cm3;
//    double gasDensity = 1.313732567606E-6*g/cm3;      // From LISE++
    double gasDensity = 1.0947771396717e-6*g/cm3;      // from internets
    fGasMaterial = new G4Material("He_CO2",gasDensity,2);
    fGasMaterial->AddMaterial(he,massFracHe);
    fGasMaterial->AddMaterial(co2,massFracCO2);

    G4Material* vacuum =
            new G4Material("Vacuum",      //Name as String
                           1,             //Atomic Number,  in this case we use 1 for hydrogen
                           1.008*g/mole,  //Mass per Mole "Atomic Weight"  1.008*g/mole for Hydoren
                           1.e-25*g/cm3,  //Density of Vaccuum  *Cant be Zero, Must be small insted
                           kStateGas,     //kStateGas for Gas
                           2.73*kelvin,   //Temperatuer for ga
                           1.e-25*g/cm3); //Pressure for Vaccum

    //Overlaps flag
    G4bool checkOverlaps = false;

    //Create a vacuum filled world
    G4VSolid* vacuumSolid = new G4Box("vacuumBox", 1.*m,1.*m,1.*m);
    fVacuumLogical = new G4LogicalVolume(vacuumSolid, vacuum, "vacuumLogical");
    G4VPhysicalVolume* vacuumPhysical = new G4PVPlacement(0,G4ThreeVector(),fVacuumLogical,"vacuumPhysical",0,false,0,checkOverlaps);

    //Create He filled world
    G4VSolid* worldSolid
            = new G4Box("worldBox",.2*m,.2*m,.123/2*m);
    fWorldLogical
            = new G4LogicalVolume(worldSolid,fGasMaterial,"worldLogical");
    G4VPhysicalVolume* worldPhysical
            = new G4PVPlacement(0,G4ThreeVector(0.,0.,.123/2.*m),fWorldLogical,"worldPhysical",fVacuumLogical,
//            = new G4PVPlacement(0,G4ThreeVector(0.,0,.0*m),fWorldLogical,"worldPhysical",0,
                                false,0,checkOverlaps);

    G4double maxstep = 0.0001*mm;
    G4UserLimits* userLimits = new G4UserLimits(maxstep);
    fWorldLogical->SetUserLimits(userLimits);

    return vacuumPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

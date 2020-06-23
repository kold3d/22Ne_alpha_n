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
/// \file DRAGONPrimaryGeneratorAction.cc
/// \brief Implementation of the DRAGONPrimaryGeneratorAction class

#include "DRAGONPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4GenericMessenger.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

DRAGONPrimaryGeneratorAction::DRAGONPrimaryGeneratorAction() :
        G4VUserPrimaryGeneratorAction() , fAlpha(0) {


    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    fAlpha = particleTable->FindParticle(particleName="alpha");

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,1));

    DefineCommands();
}

DRAGONPrimaryGeneratorAction::~DRAGONPrimaryGeneratorAction() {
    delete fParticleGun;
    delete fMessenger;
}

void DRAGONPrimaryGeneratorAction:: GeneratePrimaries(G4Event* g4Event) {
    fParticleGun->GeneratePrimaryVertex(g4Event);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,-1.0*mm));
    G4double kinE = 4.7*MeV;
//    G4double sigma = kinE*.00025/2.36;
        G4double sigma = 0;
    kinE = G4RandGauss::shoot(kinE,sigma);

    fParticleGun->SetParticleEnergy(kinE);
//    DRAGONAnalysis* analysis = DRAGONAnalysis::Instance();
//    analysis->SetGunEnergy(kinE);

}

void DRAGONPrimaryGeneratorAction::DefineCommands() {
    fMessenger = new G4GenericMessenger(this,
                                        "/DRAGON/generator/",
                                        "Primary generator control");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


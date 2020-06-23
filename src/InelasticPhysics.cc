//
// Created by Sriteja Upadhyayula on 2020-06-17.
//

#include "InelasticPhysics.hh"

#include "InelasticProcess.hh"
#include "G4GenericIon.hh"
#include "globals.hh"
#include "G4PhysicsListHelper.hh"

//#include "GeneralPhase3_2nDecay.hh"
#include "G4ionIonisation.hh"
#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(InelasticPhysics);

InelasticPhysics::InelasticPhysics(G4int)
        :  G4VPhysicsConstructor("InelasticPhysics")/*, fCMScatteringEnergy(0.)*/ {
}

InelasticPhysics::InelasticPhysics(const G4String& name)
        :  G4VPhysicsConstructor(name) {
}

InelasticPhysics::~InelasticPhysics()
{
}

void InelasticPhysics::ConstructParticle()
{
//  G4GenericIon::GenericIon();
    G4IonConstructor pIons;
    pIons.ConstructParticle();
    G4BaryonConstructor pBaryons;
    pBaryons.ConstructParticle();

//    G4Deuteron::DeuteronDefinition();
//    G4Triton::TritonDefinition();
//    G4Alpha::AlphaDefinition();
//    G4He3::He3Definition();

    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
}

void InelasticPhysics::ConstructProcess()
{
    InelasticProcess* reactionProcess = new InelasticProcess();
    reactionProcess->ParseParams(fReactionParams);
    // reactionProcess->SetCMScatteringEnergy(fCMScatteringEnergy);
    G4PhysicsListHelper::GetPhysicsListHelper()->
            RegisterProcess(reactionProcess, G4GenericIon::GenericIon());

    /*G4hMultipleScattering* multiple = new G4hMultipleScattering;
    G4ionIonisation* ion = new G4ionIonisation;
    G4String the6HennProcessName = "GeneralPhase3_2nDecay";
    GeneralPhase3_2nDecay* the6HeProcess = new GeneralPhase3_2nDecay(the6HennProcessName);
    the6HeProcess->SetParentNucleus(6,2);
    G4PhysicsListHelper::GetPhysicsListHelper()->
      RegisterProcess(multiple, G4GenericIon::GenericIon());
    G4PhysicsListHelper::GetPhysicsListHelper()->
      RegisterProcess(ion, G4GenericIon::GenericIon());
    G4PhysicsListHelper::GetPhysicsListHelper()->
      RegisterProcess(the6HeProcess, G4GenericIon::GenericIon());
    G4cout << "Registering proceses" << G4endl;*/
}




// void InelasticPhysics::SetCMScatteringEnergy(G4double energy) {
//   fCMScatteringEnergy = energy;
// }
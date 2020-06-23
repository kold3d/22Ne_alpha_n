//
// Created by Sriteja Upadhyayula on 2020-06-17.
//

#ifndef InelasticProcess_h
#define InelasticProcess_h

#include "G4VDiscreteProcess.hh"
#include "G4LogicalVolume.hh"

class InelasticProcess : public G4VDiscreteProcess {
public:
    InelasticProcess(const G4String& name = "Inelastic");
    ~InelasticProcess();

    G4double GetMeanFreePath(const G4Track&,G4double,
                             G4ForceCondition*);
    G4VParticleChange* PostStepDoIt(const G4Track&,const G4Step&);

    void StartTracking(G4Track*);
    void SetQValue(G4double qValue) {
        fQValue = qValue;
        std::cout << "SET: Q Value -- " << fQValue << std::endl;
    };
    void SetTarget(G4double charge, G4double mass) {
        fTargetMass = mass;
        fTargetCharge = charge;
        std::cout << "SET: Target -- " << fTargetCharge << ' ' << fTargetMass << std::endl;
    };
    void SetLightProduct(G4double charge, G4double mass) {
        fLightProductMass = mass;
        fLightProductCharge = charge;
        std::cout << "SET: Light Product -- " << fLightProductCharge << ' ' << fLightProductMass << std::endl;
    };
    void SetHeavyProduct(G4double charge, G4double mass) {
        fHeavyProductMass = mass;
        fHeavyProductCharge = charge;
        std::cout << "SET: Heavy Product -- " << fHeavyProductCharge << ' ' << fHeavyProductMass << std::endl;
    };
    G4double GetLightProductMass() {
        return fLightProductMass;
    }
    G4double GetLightProductCharge() {
        return fLightProductCharge;
    }
    void ParseParams(std::map<std::string,double>&);

    G4double DNdp3(G4double Mass1, G4double Mass2, G4double Mass3, G4double Eng, G4double P3);

    void ThreeBodyPhaseSpace(G4double EngDecay);

    void CMFrameToLabFrame();

    // void SetCMScatteringEnergy(G4double);

private:

    G4double Absolute(G4double Num);      // Gives absolute value of number
    G4LorentzVector LorentzBoost(G4LorentzVector Mom, G4ThreeVector Beta);
    G4double VectModulus(G4ThreeVector Mom);
    G4double ScalarProduct(G4LorentzVector V1, G4ThreeVector V2);

    G4double fScatteringEnergy;
    G4double fQValue;
    G4double fExEnergy;
    G4double fTargetMass;
    G4double fTargetCharge;
    G4double fLightProductMass;
    G4double fLightProductCharge;
    G4double fHeavyProductMass;
    G4double fHeavyProductCharge;
    G4double fCMScatteringEnergy;
    G4double fCMTheta;
    G4double fDecayE;
    G4bool fScatter;
    G4double fenergy_cm_rounded;

    G4double E_rel;

    G4double Pi;

    G4double Mass_Particle; // Stores MeV/c2 Mass of Parent Nucleus
    G4int A_Part;        // Stores Mass of Parnet Nucleus (int)
    G4int Z_Part;        // Stores Charge of Parnet Nucleus
    G4int A_Rec;         // Stores Mass of Daughter Nucleus
    G4int Z_Rec;         // Stores Charge of Daughter Nucleus

    G4double Mass_Neutron1;
    G4double Mass_Neutron2;
    G4double Mass_Recoil;

    G4LorentzVector P_Part;  // original momentum of particle
    G4LorentzVector P_Neutron1;
    G4LorentzVector P_Neutron2;
    G4LorentzVector P_Recoil;

    G4LorentzVector P_light;  // momentum of the light recoil (in this case alpha)

};

#endif

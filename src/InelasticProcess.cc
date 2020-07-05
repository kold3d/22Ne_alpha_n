//
// Created by Sriteja Upadhyayula on 2020-06-17.
//

#include "InelasticProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "DRAGONDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include <cmath>
#include <math.h>
#include "G4Ions.hh"

#include "globals.hh"
#include "G4PhysicalConstants.hh"

#include "G4VDiscreteProcess.hh"
#include "G4VParticleChange.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

#include "G4Material.hh"
#include "G4UnitsTable.hh"

#include <iostream>
#include <fstream>

InelasticProcess::InelasticProcess(const G4String& processName)
        : G4VDiscreteProcess(processName,fHadronic), fScatteringEnergy(1.e6), fCMScatteringEnergy(0.) {
    SetProcessSubType(111);
//    fQValue = -0.478;
}

InelasticProcess::~InelasticProcess() {
}

G4double InelasticProcess::GetMeanFreePath(const G4Track& aTrack,
                                           G4double /*previousStepSize*/,
                                           G4ForceCondition* condition) {

    G4double energy = aTrack.GetKineticEnergy()/MeV;

    const DRAGONDetectorConstruction* detectorConstruction
            = static_cast<const DRAGONDetectorConstruction*>
            (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    G4LogicalVolume* world = detectorConstruction->GetWorldVolume();
    G4LogicalVolume* currentVolume = aTrack.GetStep()->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();


    G4double mfp = (energy>fScatteringEnergy ||
                    aTrack.GetTrackID()>1) ? DBL_MAX : 0;

    *condition = NotForced;
    G4cout << "energy: " << energy << '\t' << '\t' << "QValue: " << fQValue <<  '\t' << "fScatteringE: " << fScatteringEnergy << '\t' << "MFP: " << mfp << G4endl;
    return mfp;
}

G4VParticleChange* InelasticProcess::PostStepDoIt( const G4Track& aTrack,
                                                   const G4Step& aStep) {


    const DRAGONDetectorConstruction* detectorConstruction
            = static_cast<const DRAGONDetectorConstruction*>
            (G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    G4LogicalVolume* world = detectorConstruction->GetWorldVolume();


    G4StepPoint* postStepPoint = aStep.GetPostStepPoint();
    G4StepPoint* preStepPoint = aStep.GetPreStepPoint();
    if (postStepPoint->GetStepStatus()==fGeomBoundary) {
        return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }

    if (preStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume()!= world &&
        postStepPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume()!=world) {
        return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }
    aParticleChange.Initialize(aTrack);

    G4double energy = aTrack.GetKineticEnergy()/MeV;
    G4double totalEnergy = energy+fQValue;
    G4ParticleDefinition* projectile = aTrack.GetDefinition();
    G4double projectileMass = projectile->GetAtomicMass();


  G4DynamicParticle* lightRec = new G4DynamicParticle;                      //Creating a G4DynamicParticle to extract 4-Momentum
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* light;
//    if(particleTable->GetIonTable()->FindIon(fLightProductCharge,fLightProductMass,0.0))
//        light = particleTable->GetIonTable()->FindIon(fLightProductCharge,fLightProductMass,0.0);
//    else light = particleTable->GetIonTable()->GetIon(fLightProductCharge,fLightProductMass,0.0);
    light = particleTable->FindParticle("neutron");
    lightRec->SetDefinition(light);                                           //Filling in the definitions of the G4Dynamicparticle

    fenergy_cm_rounded = roundf(energy*fTargetMass/(projectileMass+fTargetMass)*100.)/100.;
//    if (fenergy_cm_rounded < 2. || fenergy_cm_rounded > 8.) {
//        return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
//    }
    fDecayE = (fCMScatteringEnergy>0.) ? fCMScatteringEnergy + fQValue : roundf(energy*fTargetMass/(projectileMass+fTargetMass)*100.)/100. + fQValue;

    G4double randAngle = 180.*G4UniformRand();           //Included for the file 0-180 degrees.
        randAngle = (G4double)floor(randAngle*10+0.5)/10.;
            fCMTheta = (fQValue ==0) ? 180. - randAngle : randAngle;
    G4double aAngleLightCM = 2*3.14159*G4UniformRand();                     // Isotropic distribution
    G4double pAngleLightCM = fCMTheta*3.14159/180.;
  fExEnergy = fabs(fQValue)*MeV;
//    fExEnergy = 0.;

//  printf("TotalE: %f, ExEnergy: %f\n", totalEnergy, fExEnergy);

    G4DynamicParticle* heavyRec = new G4DynamicParticle;
    G4ParticleDefinition* heavy;
    if(particleTable->GetIonTable()->FindIon(fHeavyProductCharge,fHeavyProductMass,fExEnergy))
        heavy = particleTable->GetIonTable()->FindIon(fHeavyProductCharge,fHeavyProductMass,fExEnergy);
    else heavy = particleTable->GetIonTable()->GetIon(fHeavyProductCharge,fHeavyProductMass,fExEnergy);
    heavyRec->SetDefinition(heavy);
//    heavyRec->SetKineticEnergy(heavyEnergyLab*MeV);
//    heavyRec->SetMomentumDirection(heavyLab.unit());

//    Kinematics
    //Rotate so that angles are relative to the beam direction
    G4ThreeVector momentumDirection = aTrack.GetMomentumDirection();
    G4ThreeVector v = G4ThreeVector(0.,0.,1.).cross(momentumDirection);
    G4double rotAngle = acos(momentumDirection.z());
    G4ThreeVector dir = G4ThreeVector(sin(pAngleLightCM)*sin(aAngleLightCM),sin(pAngleLightCM)*cos(aAngleLightCM),cos(pAngleLightCM));
    if (v.getR()>0) dir.rotate(v,rotAngle);

    if(fCMScatteringEnergy<0.) return &aParticleChange; //sub-threshold
    G4double LightProductMass = fLightProductMass * 931.5*MeV;
    G4double HeavyProductMass = fHeavyProductMass * 931.5*MeV;
    G4double pLight = sqrt(2.*LightProductMass*fDecayE*(1.*HeavyProductMass/(LightProductMass+HeavyProductMass)));      //E_1 = m2/(m1+m2)* E_t
    G4double pHeavy = sqrt(2.*HeavyProductMass*fDecayE*(1.*LightProductMass/(LightProductMass+HeavyProductMass)));      //E_2 = m1/(m1+m2)* E_t
    G4ThreeVector pNewLight = pLight*dir;
    G4ThreeVector pNewHeavy = -pNewLight;
    G4ThreeVector pN = aTrack.GetMomentum();


    pNewLight += pN * (1. * LightProductMass/(LightProductMass+HeavyProductMass));
    pNewHeavy += pN * (1. * HeavyProductMass/(LightProductMass+HeavyProductMass));
    G4cout << "incomingPLab: " << pN << '\t' << "outgoingPLab: " << pNewLight << '\t' << pNewHeavy << "kLightLab: " << pow(pNewLight.getR(),2.)/(2.*LightProductMass)/MeV << "kHeavyLab: " << pow(pNewHeavy.getR(),2.)/(2.*HeavyProductMass)/MeV << G4endl;

    std::ofstream myfile;
    myfile.open("22Ne_a_n_simOutput.txt", std::ios::app);
    myfile << "incomingPLab: " << pN << '\t' << "outgoingPLab neutron: " << pNewLight << '\t' << "outputPLab 25Mg: " << pNewHeavy << G4endl;
    myfile.close();

    lightRec->SetMomentum(pNewLight);
    heavyRec->SetMomentum(pNewHeavy);

    G4double totalpLight = pNewLight.getR();
    G4double totalpHeavy = pNewHeavy.getR();

    lightRec->SetKineticEnergy(pow(totalpLight,2.)/(2.*LightProductMass));
    heavyRec->SetKineticEnergy(pow(totalpHeavy,2.)/(2.*HeavyProductMass));
    G4double pAngleLightLab = acos(pNewLight.getZ()/pNewLight.getR());

/*    G4Track* sec1 = new G4Track(new G4DynamicParticle(light,lightLab.unit(), lightEnergyLab*MeV),
                                aTrack.GetGlobalTime(),
                                aTrack.GetPosition());
    sec1->SetUserInformation(new DRAGONTrackingInformation(energy*fTargetMass/(projectileMass+fTargetMass),pAngleLightCM,pAngleLightLab,
                                                       aAngleLightCM,aTrack.GetPosition()));*/

    G4Track* sec1 = new G4Track(lightRec,
                                aTrack.GetGlobalTime(),
                                aTrack.GetPosition());
    G4cout <<"postition: " << sec1->GetPosition() << G4endl;
    G4Track* sec2 = new G4Track(heavyRec,
                                aTrack.GetGlobalTime(),
                                aTrack.GetPosition());

//    DRAGONAnalysis* analysis = DRAGONAnalysis::Instance();
//    analysis->SetCMEnergy(energy*fTargetMass/(projectileMass+fTargetMass));
//    analysis->SetCMTheta(pAngleLightCM);
//    analysis->SetLabTheta(pAngleLightLab);
//    analysis->SetCMPhi(aAngleLightCM);
//    G4cout << "fCMScatteringEnergy: " << fCMScatteringEnergy << '\t' << "decayE: " << fDecayE << '\t' << "cmERounded: " << fenergy_cm_rounded << '\t' << "SetCMEnergy: " << energy*fTargetMass/(projectileMass+fTargetMass) << G4endl;
//    G4cout << fDecayE << '\t' << fenergy_cm_rounded << '\t' << energy*fTargetMass/(projectileMass+fTargetMass) << G4endl;


    /*(printf("IncomingZ: %d incomingM: %d recoilZ: %d recoilM: %d \n Lab Energy Incoming: %f Light Product Energy: %f Heavy Product Energy: %f \n Light Product Angle: %f Heavy Product Angle %f\n",
           projectile->GetAtomicNumber(),projectile->GetAtomicMass(),heavy->GetAtomicNumber(),heavy->GetAtomicMass(),
           energy,lightEnergyLab,heavyEnergyLab,pAngleLightLab,pAngleHeavyLab);
*/

    aParticleChange.SetNumberOfSecondaries(2);            //change to 4 if adding neutron tracks
    aParticleChange.AddSecondary(sec1);                  //Scattered alpha track

    // Kill parent track
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    return &aParticleChange;
}

void InelasticProcess::StartTracking(G4Track* track) {
    G4VProcess::StartTracking(track); // Apply base class actions
/*    fCMScatteringEnergy = fCMScatteringEnergy*MeV;
    G4bool goodCS = false;
    while (!goodCS) {
        G4bool goodCME = false;
        G4double randE;
        while (!goodCME) {
            randE = 6.*G4UniformRand() + 2.;
            randE = (G4double)floor(randE*100+0.5)/100.;
//            Biasing with effective target thickness
            G4double cmETT = getCMETT(randE);
            G4double randAttackTT = G4UniformRand();
            if (randAttackTT < cmETT) {
                goodCME = true;
            }
        }
        G4double randAngle = 180.*G4UniformRand();           //Included for the file 0-180 degrees.
        randAngle = (G4double)floor(randAngle*10+0.5)/10.;
        G4double csCDF = getCSNorm(randE, randAngle);
        G4double randAttack = G4UniformRand();
        if (randAttack < csCDF) {
            fScatteringEnergy = (fCMScatteringEnergy>0.) ? (fCMScatteringEnergy*(4.+6.)/4.) : ((4.+6.)/4.)*randE/MeV;
            fDecayE = (fCMScatteringEnergy>0.) ? fCMScatteringEnergy + fQValue : randE + fQValue;
            fCMTheta = randAngle;
            fScatter = true;
        } else {
     //       randE = 25.;
     //       fScatteringEnergy = (fCMScatteringEnergy>0.) ? (fCMScatteringEnergy*(4.+6.)/4.) : ((4.+6.)/4.)*randE/MeV;
     //       fDecayE = (fCMScatteringEnergy>0.) ? fCMScatteringEnergy + fQValue : randE + fQValue;
     //       fCMTheta = randAngle;
            fScatter = false;
        }
        goodCS = true;
//       G4cout << "fScatterStartTracking " << fScatter << G4endl;
    }*/
    fScatteringEnergy = (fCMScatteringEnergy>0.) ? (fCMScatteringEnergy*(4.+22.)/4.) : track->GetKineticEnergy()*G4UniformRand()/MeV;
//     G4cout << "fCMScattering Energy: " << fCMScatteringEnergy << '\t' << fLightProductCharge << " " << fLightProductMass << '\t' << fHeavyProductCharge << " " << fHeavyProductMass << G4endl;
}

void InelasticProcess::ParseParams(std::map<std::string,double> &params) {
    double lightProductMass = -1;
    double lightProductCharge = -1;
    double heavyProductMass = -1;
    double heavyProductCharge = -1;
    double targetMass = -1;
    double targetCharge = -1;
    // double cmEnergy = -1;
    for(std::map<std::string,double>::const_iterator it = params.begin();
        it!= params.end();it++) {
        if(it->first == "qValue"  ) {
            SetQValue(it->second);
        } else if(it->first == "lightProductMass" ) {
            lightProductMass = it->second;
        } else if(it->first == "lightProductCharge" ) {
            lightProductCharge = it->second;
        } else if(it->first == "heavyProductMass" ) {
            heavyProductMass = it->second;
        } else if(it->first == "heavyProductCharge" ) {
            heavyProductCharge = it->second;
        } else if(it->first == "targetMass" ) {
            targetMass = it->second;
        } else if(it->first == "targetCharge" ) {
            targetCharge = it->second;
        } else if(it->first == "cmEnergy") {
            fCMScatteringEnergy = it->second;
        }
    }
    if(lightProductCharge>=0 && lightProductMass>=0) {
        SetLightProduct(lightProductCharge,lightProductMass);
    }
    if(heavyProductCharge>0 && heavyProductMass>0) {
        SetHeavyProduct(heavyProductCharge,heavyProductMass);
    }
    if(targetCharge>0 && targetMass>0) {
        SetTarget(targetCharge,targetMass);
    }
}

// void InelasticProcess::SetCMScatteringEnergy(G4double energy) {
//   fCMScatteringEnergy = energy;
// }

G4double InelasticProcess::Absolute(G4double Num)
{
    if(Num < 0.)
    {Num *= -1.;}

    return Num;
}

G4double InelasticProcess::DNdp3(G4double Mass1, G4double Mass2, G4double Mass3, G4double Eng, G4double P3)
{
    // Momentum Probability Distribution Function

    G4double Eng3 = sqrt(pow(Mass3,2)+pow(P3,2));

    G4double Fact1 = pow(P3,2)/Eng3;

    G4double Fact2 = (pow(Eng,2)+pow(Mass3,2)-2.*Eng*Eng3-pow((Mass1+Mass2),2))*
                     (pow(Eng,2)+pow(Mass3,2)-2.*Eng*Eng3-pow((Mass1-Mass2),2));

    G4double Denom = pow(Eng,2)+pow(Mass3,2)-2.*Eng*Eng3;

    G4double dNdp3 = Fact1*sqrt(Absolute(Fact2))/Denom;

    return dNdp3;
}

/*void InelasticProcess::ThreeBodyPhaseSpace(G4double EnergyDecay)
{
    const G4int NDiv = 200;

    G4double Mass1 = Mass_Neutron1;
    G4double Mass2 = Mass_Neutron2;
    G4double Mass3 = Mass_Recoil;

    // Four - Vectors for Phase-Space Generation
    G4LorentzVector P1;
    G4LorentzVector P2;
    G4LorentzVector P3;
    G4ThreeVector Beta;

    G4double Energy = Mass1+Mass2+Mass3+EnergyDecay;

    G4double P3M = sqrt(Absolute((pow(Energy,2)-pow((Mass1+Mass2+Mass3),2))*
                                 (pow(Energy,2)-pow((Mass1+Mass2-Mass3),2))))/(2.*Energy);

    // Generate Distribution Function

    G4double F_Old = 0.;
    G4double P_Old = 0.;
    G4double StepSize = P3M/(static_cast<G4double>(NDiv));

    G4double P_Max = P3M - StepSize;
    G4double F_Max = DNdp3(Mass1,Mass2,Mass3,Energy,P_Max);

    while(F_Max > F_Old)
    {
        P_Old = P_Max;
        F_Old = F_Max;
        P_Max = P_Old-StepSize;
        F_Max = DNdp3(Mass1,Mass2,Mass3,Energy,P_Max);
    }
    F_Max = DNdp3(Mass1,Mass2,Mass3,Energy,(P_Max+StepSize/2.));

    G4double P3_Mod;
    G4double Func3;

    do{
        P3_Mod = P3M*G4UniformRand();
        Func3 = F_Max*G4UniformRand();
    }
    while(Func3 > DNdp3(Mass1,Mass2,Mass3,Energy,P3_Mod));

    // Generate an isotropic theta, phi
    G4double theta = acos(1.-2.*G4UniformRand());
    G4double phi = 2.*M_PI*G4UniformRand();
//    G4cout << "phi" << phi << '\t' << M_PI << '\t' << "theta" << theta <<G4endl;


    // Create Momentum vectors for Particle3

    P3[0] = P3_Mod*sin(theta)*cos(phi);
    P3[1] = P3_Mod*sin(theta)*sin(phi);
    P3[2] = P3_Mod*cos(theta);
    P3[3] = sqrt((pow(Mass3,2))+(pow(P3_Mod,2)));

    // Continue with a two-body phase space in (12) center of mass

    G4double Mass12 = sqrt(pow((Energy-P3[3]),2) - pow(P3_Mod,2));

    G4double P1_Mod = sqrt(Absolute((pow(Mass12,2)-pow((Mass1+Mass2),2))*
                                    (pow(Mass12,2)-pow((Mass1-Mass2),2))))/(2.*Mass12);

//    DRAGONAnalysis* analysis = DRAGONAnalysis::Instance();
//    analysis->SetLabPhi(phi);
//    analysis->SetLabTheta(theta);
//    P1[0] = P1_Mod*sin(theta)*cos(phi);
//    P1[1] = P1_Mod*sin(theta)*sin(phi);
//    P1[2] = P1_Mod*cos(theta);
//    P1[3] = sqrt((pow(Mass1,2))+(pow(P1_Mod,2)));
//
//    for(G4int i=0;i<3;i++)
//    {P2[i] = - P1[i];}
//    P2[3] = sqrt((pow(Mass2,2))+(pow(P1_Mod,2)));
//
//    // Now put V1 and V2 into frame of V3
//
//    for(G4int i=0;i<3;i++)
//    { Beta[i] = -P3[i]/(Energy-P3[3]); }
//
//    G4LorentzVector P1_Final;
//    G4LorentzVector P2_Final;
//    G4LorentzVector P3_Final;
//
//    P1_Final = LorentzBoost(P1,Beta);
//    P2_Final = LorentzBoost(P2,Beta);
//    P3_Final = P3;
//    // Do nothing with P3_Final (P3)
//
//    // Update Particle Vectors
//    P_Neutron1 = P1_Final;
//    P_Neutron2 = P2_Final;
//    P_Recoil = P3_Final;

    // End of Routine
}*/

G4LorentzVector InelasticProcess::LorentzBoost(G4LorentzVector Mom, G4ThreeVector Beta)
{
    // Boost routine originally written in Fortran by F.Miguel Marques
    // "Adds" a Momentum vector in the rest frame to the Beta
    // thereby giving a Lorentz Transformation of two 4-vectors
    // Translated to C++/GEANT4 by : BTR

    G4double Norm_Beta = VectModulus(Beta);

    if(Norm_Beta == 1.)
    {
        G4cout << " Error! Vector |Beta| is at Speed of Light! - Norm_Beta = " << Norm_Beta << G4endl;
        exit(0);
    }
    else if(Norm_Beta == 0.)
    {
        G4cout << "No Boost, |Beta| = 0" << G4endl;
        return Mom;
    }

    G4double gamma = 1./sqrt(1.-pow(Norm_Beta,2));    // Calculates gamma factor
    // G4cout << "gamma = " << gamma << G4endl;

    // Now add the Lorentz Vectors by Scalar Dot Product
    // and Matrix Transform

    G4double theScalarAmp = ScalarProduct(Mom,Beta);
    G4double P_Parallel = theScalarAmp/Norm_Beta;
    G4double E_total = Mom[3];

    G4ThreeVector P_Perp;

    for(G4int i=0; i<3; i++)
    {P_Perp[i] = Mom[i]-(P_Parallel*Beta[i]/Norm_Beta);}

    P_Parallel = gamma*(P_Parallel+(Norm_Beta*E_total));
    Mom[3] = gamma*(E_total+theScalarAmp);
    // Update Output Momentum Vector

    for(G4int i=0; i<3; i++)
    {Mom[i] = P_Perp[i]+(P_Parallel*Beta[i]/Norm_Beta);}

    return Mom;

    // End of function Boost!
}

G4double InelasticProcess::VectModulus(G4ThreeVector Mom)
{
    // Gives Modulus (normal) of a 3-vector
    G4double V_Modulus = sqrt(pow(Mom[0],2)+pow(Mom[1],2)+pow(Mom[2],2));
    return V_Modulus;
}

G4double InelasticProcess::ScalarProduct(G4LorentzVector V1, G4ThreeVector V2)
{
    // Gives Scalar Product of a Lorentz Vector and a ThreeVector
    G4double theScalarProduct = V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2];
    return theScalarProduct;
}

void InelasticProcess::CMFrameToLabFrame()
{
    // Takes CM distribition in BacktoBack() and Boosts to Lab Frame
    // Follows routine written in Fortran for Bzb written by: JL Lecouey
    // Translated and updated for C++/GEANT by : BTR

    //G4cout << " The Particle Momentum Direction (P_Part) is : " << P_Part << G4endl;

    G4ThreeVector Beta_Part;    // Turn Beam Momentum into normalized velocity vector
    for(G4int i=0;i<3;i++)
    {Beta_Part[i] = P_Part[i]/P_Part[3];}

// Now add Particle Lorentz Vector to Momentum Vectors for the neutrons and recoil

    P_Neutron1 = LorentzBoost(P_Neutron1, Beta_Part);
    P_Neutron2 = LorentzBoost(P_Neutron2, Beta_Part);
    P_Recoil = LorentzBoost(P_Recoil, Beta_Part);
}

double efficency (double cmEnergy) {
    double a = 2.74519603e+01;
    double b = 9.85632981e+00;
    double c = 2.73735136e-05;
    double eff = a/pow(cmEnergy,b)+c;
    return eff;
}
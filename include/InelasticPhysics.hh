//
// Created by Sriteja Upadhyayula on 2020-06-17.
//

#ifndef InelasticPhysics_h
#define InelasticPhysics_h

#include "G4VPhysicsConstructor.hh"
#include <string>
#include <map>

class InelasticPhysics : public G4VPhysicsConstructor {
public:
    InelasticPhysics(G4int verbose =1);
    InelasticPhysics(const G4String& name);
    // void SetCMScatteringEnergy(G4double);
    virtual ~InelasticPhysics();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

    void SetReactionParams(std::map<std::string,double> params) {
        fReactionParams=params;
    };

private:
    std::map<std::string,double> fReactionParams;
    // G4double fCMScatteringEnergy;
};

#endif

/* 
 * File:   NewPhysicsSTU.h
 * Author: mishima
 */

#ifndef NEWPHYSICSSTU_H
#define	NEWPHYSICSSTU_H

#include "StandardModel.h"


class NewPhysicsSTU : public StandardModel {
public:
    static const int NSTUvars = 3;
    static const std::string STUvars[NSTUvars];
    
    /**
     * @brief NewPhysicsSTU constructor
     */
    NewPhysicsSTU();

    virtual std::string ModelName() const {
        return "NewPhysicsSTU";
    }

    virtual bool Update(const std::map<std::string, double>& DPars);
    virtual bool Init(const std::map<std::string, double>& DPars);    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    
    ////////////////////////////////////////////////////////////////////////     

    /**
     * @return Oblique parameter S
     */
    double obliqueS() const {
        return myObliqueS;
    }

    /**
     * @return Oblique parameter T
     */
    double obliqueT() const {
        return myObliqueT;
    }

    /**
     * @return Oblique parameter U
     */
    double obliqueU() const {
        return myObliqueU;
    }
    
    
    ////////////////////////////////////////////////////////////////////////     
    
protected:    
    virtual void SetParameter(const std::string name, const double& value);
    double myObliqueS, myObliqueT, myObliqueU;

    ////////////////////////////////////////////////////////////////////////     
    
private:

};

#endif	/* NEWPHYSICSSTU_H */


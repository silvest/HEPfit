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
    double getObliqueS() const {
        return obliqueS;
    }

    //void setObliqueS(double obliqueS) {
    //    this->obliqueS = obliqueS;
    //}

    /**
     * @return Oblique parameter T
     */
    double getObliqueT() const {
        return obliqueT;
    }

    //void setObliqueT(double obliqueT) {
    //    this->obliqueT = obliqueT;
    //}

    /**
     * @return Oblique parameter U
     */
    double getObliqueU() const {
        return obliqueU;
    }

    //void setObliqueU(double obliqueU) {
    //    this->obliqueU = obliqueU;
    //}
    
    
    ////////////////////////////////////////////////////////////////////////     
    
protected:    
    virtual void SetParameter(const std::string name, const double& value);
    double obliqueS, obliqueT, obliqueU;


    ////////////////////////////////////////////////////////////////////////     
    
private:

};

#endif	/* NEWPHYSICSSTU_H */


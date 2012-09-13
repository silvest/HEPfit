/* 
 * File:   NewPhysicsSTU.h
 * Author: mishima
 */

#ifndef NEWPHYSICSSTU_H
#define	NEWPHYSICSSTU_H

#include "StandardModel.h"


class NewPhysicsSTU : public StandardModel {
public:
    static const int NSTUvars = 7;
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

    /**
     * @return Oblique parameter V
     */
    double obliqueV() const {
        return myObliqueV;
    }

    /**
     * @return Oblique parameter W
     */
    double obliqueW() const {
        return myObliqueW;
    }

    /**
     * @return Oblique parameter X
     */
    double obliqueX() const {
        return myObliqueX;
    }

    /**
     * @return Oblique parameter Y
     */
    double obliqueY() const {
        return myObliqueY;
    }
    
    
    ////////////////////////////////////////////////////////////////////////     
    
protected:    
    virtual void SetParameter(const std::string name, const double& value);
    double myObliqueS, myObliqueT, myObliqueU;
    double myObliqueV, myObliqueW, myObliqueX, myObliqueY;
    

    ////////////////////////////////////////////////////////////////////////     
    
private:

};

#endif	/* NEWPHYSICSSTU_H */


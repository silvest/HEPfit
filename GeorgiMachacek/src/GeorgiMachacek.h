/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GEORGIMACHACEK_H
#define	GEORGIMACHACEK_H

#include <StandardModel.h>

class GMcache; //forward reference to GMcache class

/**
 * @class GeorgiMachacek
 * @brief The Georgi-Machacek model
 */
class GeorgiMachacek: public StandardModel {
public:

    static const int NGMvars = 8;
    //The parameters of the Higgs potential for Georgi Machacek model according to 1502.01275v1
    //We choose the physical basis: tantheta,alpha,mHh,mH3,mH5,M1,M2
    static const std::string GMvars[NGMvars];
    
    /**
     * @brief GeorgiMachacek constructor
     */
    GeorgiMachacek();
    
    /**
     * @brief GeorgiMachacek destructor
     */
    ~GeorgiMachacek();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
//    virtual bool setFlag(const std::string name, const bool value);

    GMcache* getMyGMCache() const
    {
        return myGMcache;
    }

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    /**
     *
     * @return \f$\log(\tan \theta_H)$\f
     */
    double getLogttheta() const {
        return logttheta;
    }

    /**
     *
     * @return \f$\tan \theta_H$\f
     */
    double getTantheta() const {
        return tantheta;
    }

    /**
     *
     * @return \f$\sin \theta_H$\f
     */
    double getSintheta() const {
        return sintheta;
    }

    /**
     *
     * @return \f$\cos \theta_H$\f
     */
    double getCostheta() const {
        return costheta;
    }

    /**
     *
     * @return \f$\alpha$\f
     */
    double getalpha() const {
        return alpha;
    }

   /**
     *
     * @return \f$\sin(\theta_H-\alpha)$\f
     */
    double getsin_tma() const {
        return sin_tma;
    }

    double computeCosa() const;

    double computeSina() const;

    double computelambda1() const;

    double computelambda2() const;

    double computelambda3() const;

    double computelambda4() const;

    double computelambda5() const;

    double computemu2sq() const;

    double computemu3sq() const;

    /**
     *
     * @return heavy singlet mass
     */
    double getMHh() const {
        return mHh;
    }

    /**
     *
     * @return triplet mass
     */
    double getMH3() const {
        return mH3;
    }

    /**
     *
     * @return quintuplet mass
     */
    double getMH5() const {
        return mH5;
    }

    /**
     *
     * @return massive parameter of the scalar potential \f$M_1$\f
     */
    double getM1() const {
        return M1;
    }

    /**
     *
     * @return massive parameter of the scalar potential \f$M_2$\f
     */
    double getM2() const {
        return M2;
    }

    /**
     *
     * @return Georgi-Machacek scale
     */
    double getQ_GM() const {
        return Q_GM;
    }


protected: 
    
    virtual void setParameter(const std::string, const double&);

private:

    GMcache* myGMcache;

    double logttheta, tantheta, sintheta, costheta, alpha, sin_tma, mHh, mH3, mH5, M1, M2, Q_GM;
    double sign(const double x) const;

};

#endif	/* GEORGIMACHACEK_H */

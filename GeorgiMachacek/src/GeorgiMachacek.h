/* 
 * Copyright (C) 2015 SusyFit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GEORGIMACHACEK_H
#define	GEORGIMACHACEK_H

#include "StandardModel.h"
#include "GMMatching.h"

class GMcache; //forward reference to GMcache class

/**
 * @class GeorgiMachacek
 * @ingroup GeorgiMachacek
 * @brief A base class for the Georgi-Machacek model.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */
class GeorgiMachacek: public StandardModel {
public:

    static const int NGMvars = 8;
    //The parameters of the Higgs potential for Georgi Machacek model according to 1511.00865v1
    //We choose the physical basis: tanbeta, alpha, mHh=mH1, mA=mH3, mH5, Mu1, Mu2
    static std::string GMvars[NGMvars];
    
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

    /**
     * @brief A get method to access the member reference of type StandardModelMatching.
     * @return a reference to a StandardModelMatching object
     */
    virtual GMMatching& getMatching() const
    {
        return GMM.getObj();
    }

    GMcache* getMyGMCache() const
    {
        return myGMcache;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Flags

//    virtual bool setFlagStr(const std::string name, const std::string value);
    virtual bool setFlag(const std::string, const bool);

    ///////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////

    /**
     *
     * @return \f$\log(\tan \beta)$\f
     */
    double getlogtb() const {
        return logtb;
    }

    /**
     *
     * @return \f$\tan \beta$\f
     */
    double gettanb() const {
        return tanb;
    }

    /**
     *
     * @return \f$\sin \beta$\f
     */
    double getsinb() const {
        return sinb;
    }

    /**
     *
     * @return @f$\cos \beta@f$
     */
    double getcosb() const {
        return cosb;
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
     * @return @f$\beta-\alpha@f$
     */
    double getbma() const {
        return atan(tanb)-alpha;
    }

    /**
     *
     * @return @f$\sin(\beta-\alpha)@f$
     */
    double getsin_ba() const {
        return sin_ba;
    }

    /**
     *
     * @return @f$\cos \alpha@f$
     */
    double getcosa() const{
        return cos(alpha);
    }

    /**
     *
     * @return @f$\sin \alpha@f$
     */
    double getsina() const{
        return sin(alpha);
    }

    /**
     *
     * @return squared mass of the lighter singlet Higgs
     */
    double getmHl2() const {
        if(flag_use_sq_masses) {
            if(mHhsq < mHl2) {
                return mHhsq;
            }
            else
            {
                return mHl2;
            }
        }
        else
        {
            if(mHh*mHh < mHl2) {
                return mHh*mHh;
            }
            else
            {
                return mHl2;
            }
        }
    }

    /**
     *
     * @return squared mass of the heavier singlet Higgs
     */
    double getmHh2() const {
        if(flag_use_sq_masses) {
            if(mHhsq < 0.) {
                return 0.;
            }
            else if(mHhsq < mHl2) {
                return mHl2;
            }
            else
            {
                return mHhsq;
            }
        }
        else
        {
            if(mHh*mHh < mHl2) {
                return mHl2;
            }
            else
            {
                return mHh*mHh;
            }
        }
    }

    /**
     *
     * @return mass of the heavier singlet Higgs
     */
    double getmHh() const {
        if(flag_use_sq_masses) {
            if(mHhsq < 0.) {
                return 0.;
            }
            else if(mHhsq < mHl2) {
                return sqrt(mHl2);
            }
            else
            {
                return sqrt(mHhsq);
            }
        }
        else
        {
            if(mHh*mHh < mHl2) {
                return sqrt(mHl2);
            }
            else
            {
                return mHh;
            }
        }
    }

    /**
     *
     * @return squared mass of the singlet Higgs input
     */
    double getinputmHh2() const {
        if(flag_use_sq_masses) {
            if(mHhsq < 0.) {
                return 0.;
            }
                return mHhsq;
        }
        else
        {
                return mHh*mHh;
        }
    }

    /**
     *
     * @return squared triplet mass
     */
    double getmAsq() const {
        if(flag_use_sq_masses) {
            return mAsq;
        }
        else
        {
            return mA*mA;
        }
    }

    /**
     *
     * @return triplet mass
     */
    double getmA() const {
        if(flag_use_sq_masses) {
            if(mA < 0.) {
                return 0.;
            }
            else
            {
                return sqrt(mAsq);
            }
        }
        else
        {
                return mA;
        }
    }

    /**
     *
     * @return squared quintet mass
     */
    double getmH5sq() const {
        if(flag_use_sq_masses) {
            return mH5sq;
        }
        else
        {
            return mH5*mH5;
        }
    }

    /**
     *
     * @return quintet mass
     */
    double getmH5() const {
        if(flag_use_sq_masses) {
            if(mH5 < 0.) {
                return 0.;
            }
            else
            {
                return sqrt(mH5sq);
            }
        }
        else
        {
                return mH5;
        }
    }

    /**
     *
     * @return massive parameter of the scalar potential \f$\mu_1$\f
     */
    double getMu1() const {
        return Mu1;
    }

    /**
     *
     * @return massive parameter of the scalar potential \f$\mu_2$\f
     */
    double getMu2() const {
        return Mu2;
    }

    /**
     *
     * @return Georgi-Machacek scale
     */
    double getQ_GM() const {
        return Q_GM;
    }

//    double computelambda1() const;
//
//    double computelambda2() const;
//
//    double computelambda3() const;
//
//    double computelambda4() const;
//
//    double computelambda5() const;
//
//    double computemu2sq() const;
//
//    double computemu3sq() const;


protected: 
    
    virtual void setParameter(const std::string, const double&);

    mutable Matching<GMMatching,GeorgiMachacek> GMM; ///< An object of type Matching.

private:

    GMcache* myGMcache;

    double logtb, tanb, sinb, cosb, alpha, sin_ba, mHh, mA, mH5, mHhsq, mAsq, mH5sq, Mu1, Mu2, Q_GM;
    double mHl2;
    bool flag_use_sq_masses;
//    double sign(const double x) const;

};

#endif	/* GEORGIMACHACEK_H */

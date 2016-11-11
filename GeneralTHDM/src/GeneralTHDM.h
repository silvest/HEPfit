/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDM_H
#define	GENERALTHDM_H

#include "StandardModel.h"
//#include "GeneralTHDMMatching.h"

/**
 * @addtogroup GeneralTHDM
 * @brief A module for general Two-Higgs-Doublet models.
 * @{
 */


class GeneralTHDMcache; //forward reference to GeneralTHDMcache class

/**
 * @class GeneralTHDM
 * @ingroup GeneralTHDM
 * @brief A base class for generic Two-Higgs-Doublet models.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 *
 * 
 * @anchor GeneralTHDMParameters
 * <h3>%Model parameters</h3>
 */
class GeneralTHDM: public StandardModel {
public:

    static const int NGeneralTHDMvars = 12;
    static const std::string GeneralTHDMvars[NGeneralTHDMvars];

    /**
     * @brief GeneralTHDM constructor
     */
    GeneralTHDM();

    /**
     * @brief GeneralTHDM destructor
     */
    ~GeneralTHDM();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);

//    /**
//     * @brief A get method to access the member reference of type StandardModelMatching.
//     * @return a reference to a StandardModelMatching object
//     */
//    virtual GeneralTHDMMatching& getMatching() const
//    {
//        return GTHDMM.getObj();
//    }

//    virtual void setMatching(GeneralTHDMMatching& GTHDMMr)
//    {
//        GTHDMM.setObj(GTHDMMr);
//    }


    ///////////////////////////////////////////////////////////////////////////
    // Flags

//    virtual bool setFlagStr(const std::string name, const std::string value);
//    virtual bool setFlag(const std::string, const bool);

    GeneralTHDMcache* getMyGTHDMCache() const
    {
        return myGTHDMcache;
    }

    ///////////////////////////////////////////////////////////////////////////

    /**
     * @return @f$\log_{10}(\tan \beta)@f$
     */
    double getlogtb() const {
        return logtb;
    }

    /**
     * @return @f$\tan \beta@f$
     */
    double gettanb() const {
        return tanb;
    }

    /**
     * @return @f$\sin \beta@f$
     */
    double getsinb() const {
        return sinb;
    }

    /**
     * @return @f$\cos \beta@f$
     */
    double getcosb() const {
        return cosb;
    }

    /**
     * @return mass squared of the charged Higgs
     */
    double getmHp2() const {
        return mHp2;
    }

    /**
     * @return mass of the charged Higgs
     */
    double getmHp() const {
    if(mHp2 < 0.) {
            throw std::runtime_error("error in GeneralTHDM: mHp2 < 0!");
    }
        return sqrt(mHp2);
    }

    /**
     * @return rotation angle between the first and second neutral mass eigenstate
     */
    double getalpha1() const {
        return alpha1;
    }

    /**
     * @return cos of the rotation angle between the first and second neutral mass eigenstate
     */
    double getcosalpha1() const {
        return cosalpha1;
    }

    /**
     * @return sin of the rotation angle between the first and second neutral mass eigenstate
     */
    double getsinalpha1() const {
        return sinalpha1;
    }

    /**
     * @return rotation angle between the first and third neutral mass eigenstate
     */
    double getalpha2() const {
        return alpha2;
    }

    /**
     * @return cos of the rotation angle between the first and third neutral mass eigenstate
     */
    double getcosalpha2() const {
        return cosalpha2;
    }

    /**
     * @return sin of the rotation angle between the first and third neutral mass eigenstate
     */
    double getsinalpha2() const {
        return sinalpha2;
    }

    /**
     * @return rotation angle between the second and third neutral mass eigenstate
     */
    double getalpha3() const {
        return alpha3;
    }

    /**
     * @return cos of the rotation angle between the second and third neutral mass eigenstate
     */
    double getcosalpha3() const {
        return cosalpha3;
    }

    /**
     * @return sin of the rotation angle between the second and third neutral mass eigenstate
     */
    double getsinalpha3() const {
        return sinalpha3;
    }

    /**
     * @return Real part of the Higgs potential parameter @f$m_{12}^2@f$
     */
    double getRem12_2() const {
        return Rem12_2;
    }
    
    /**
     * @return Imaginary part of the Higgs potential parameter @f$lambda_{5}@f$
     */
    double getImlambda5() const {
        return Imlambda5;
    }
    
    /**
     * @return Real part of the Higgs potential parameter @f$lambda_{6}@f$
     */
    double getRelambda6() const {
        return Relambda6;
    }

    /**
     * @return Imaginary part of the Higgs potential parameter @f$lambda_{6}@f$
     */
    double getImlambda6() const {
        return Imlambda6;
    }

    /**
     * @return Real part of the Higgs potential parameter @f$lambda_{7}@f$
     */
    double getRelambda7() const {
        return Relambda7;
    }

    /**
     * @return Imaginary part of the Higgs potential parameter @f$lambda_{7}@f$
     */
    double getImlambda7() const {
        return Imlambda7;
    }

    /**
     * @return General THDM scale
     */
    double getQ_GTHDM() const {
        return Q_GTHDM;
    }
    
    /**
     * 
     * @return @f$M^2@f$ of @cite Kanemura:2015ska
     */
    double getM2() const {
            return Rem12_2/sinb/cosb;
    }

protected:

    virtual void setParameter(const std::string, const double&);

//    /**
//     * @brief A method to check if the model type name in string form is valid.
//     * @param[in] THDM model type name
//     * @return a boolean that is true if the model type name is valid
//     */
//    bool checkmodelType(const std::string modeltype) const
//    {
//        if (modeltype.compare("type1") == 0
//                || modeltype.compare("type2") == 0
//                || modeltype.compare("typeX") == 0
//                || modeltype.compare("typeY") == 0)
//            return true;
//        else
//            return false;
//    }

//    mutable Matching<THDMMatching,THDM> THDMM; ///< An object of type Matching.

private:

    GeneralTHDMcache* myGTHDMcache;

    double logtb, tanb, sinb, cosb, mHp2, alpha1, cosalpha1, sinalpha1, alpha2, cosalpha2, sinalpha2,
            alpha3, cosalpha3, sinalpha3, Rem12_2, Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7, Q_GTHDM;
};

/**
 * @}
 */

#endif	/* GENERALTHDM_H */

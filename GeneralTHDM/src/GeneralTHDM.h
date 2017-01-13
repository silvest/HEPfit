/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDM_H
#define	GENERALTHDM_H

#include "StandardModel.h"
#include "GeneralTHDMMatching.h"

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

    static const int NGeneralTHDMvars = 66;
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

    /**
     * @brief A get method to access the member reference of type StandardModelMatching.
     * @return a reference to a StandardModelMatching object
     */
    virtual GeneralTHDMMatching& getMatching() const
    {
        return GTHDMM.getObj();
    }

//    virtual void setMatching(GeneralTHDMMatching& GTHDMMr)
//    {
//        GTHDMM.setObj(GTHDMMr);
//    }


    ///////////////////////////////////////////////////////////////////////////
    // Flags

//    virtual bool setFlagStr(const std::string name, const std::string value);
    virtual bool setFlag(const std::string, const bool);

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

    double getYtu_11r() const {
        return Ytu_11r;
    }

    double getYtu_11i() const {
        return Ytu_11i;
    }

    double getYtu_12r() const {
        return Ytu_12r;
    }

    double getYtu_12i() const {
        return Ytu_12i;
    }

    double getYtu_13r() const {
        return Ytu_13r;
    }

    double getYtu_13i() const {
        return Ytu_13i;
    }

    double getYtu_21r() const {
        return Ytu_21r;
    }

    double getYtu_21i() const {
        return Ytu_21i;
    }

    double getYtu_22r() const {
        return Ytu_22r;
    }

    double getYtu_22i() const {
        return Ytu_22i;
    }

    double getYtu_23r() const {
        return Ytu_23r;
    }

    double getYtu_23i() const {
        return Ytu_23i;
    }

    double getYtu_31r() const {
        return Ytu_31r;
    }

    double getYtu_31i() const {
        return Ytu_31i;
    }

    double getYtu_32r() const {
        return Ytu_32r;
    }

    double getYtu_32i() const {
        return Ytu_32i;
    }

    double getYtu_33r() const {
        return Ytu_33r;
    }

    double getYtu_33i() const {
        return Ytu_33i;
    }

    double getYtd_11r() const {
        return Ytd_11r;
    }

    double getYtd_11i() const {
        return Ytd_11i;
    }

    double getYtd_12r() const {
        return Ytd_12r;
    }

    double getYtd_12i() const {
        return Ytd_12i;
    }

    double getYtd_13r() const {
        return Ytd_13r;
    }

    double getYtd_13i() const {
        return Ytd_13i;
    }

    double getYtd_21r() const {
        return Ytd_21r;
    }

    double getYtd_21i() const {
        return Ytd_21i;
    }

    double getYtd_22r() const {
        return Ytd_22r;
    }

    double getYtd_22i() const {
        return Ytd_22i;
    }

    double getYtd_23r() const {
        return Ytd_23r;
    }

    double getYtd_23i() const {
        return Ytd_23i;
    }

    double getYtd_31r() const {
        return Ytd_31r;
    }

    double getYtd_31i() const {
        return Ytd_31i;
    }

    double getYtd_32r() const {
        return Ytd_32r;
    }

    double getYtd_32i() const {
        return Ytd_32i;
    }

    double getYtd_33r() const {
        return Ytd_33r;
    }

    double getYtd_33i() const {
        return Ytd_33i;
    }

    double getYtl_11r() const {
        return Ytl_11r;
    }

    double getYtl_11i() const {
        return Ytl_11i;
    }

    double getYtl_12r() const {
        return Ytl_12r;
    }

    double getYtl_12i() const {
        return Ytl_12i;
    }

    double getYtl_13r() const {
        return Ytl_13r;
    }

    double getYtl_13i() const {
        return Ytl_13i;
    }

    double getYtl_21r() const {
        return Ytl_21r;
    }

    double getYtl_21i() const {
        return Ytl_21i;
    }

    double getYtl_22r() const {
        return Ytl_22r;
    }

    double getYtl_22i() const {
        return Ytl_22i;
    }

    double getYtl_23r() const {
        return Ytl_23r;
    }

    double getYtl_23i() const {
        return Ytl_23i;
    }

    double getYtl_31r() const {
        return Ytl_31r;
    }

    double getYtl_31i() const {
        return Ytl_31i;
    }

    double getYtl_32r() const {
        return Ytl_32r;
    }

    double getYtl_32i() const {
        return Ytl_32i;
    }

    double getYtl_33r() const {
        return Ytl_33r;
    }

    double getYtl_33i() const {
        return Ytl_33i;
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

    mutable Matching<GeneralTHDMMatching,GeneralTHDM> GTHDMM; ///< An object of type Matching.

private:

    GeneralTHDMcache* myGTHDMcache;

    double logtb, tanb, sinb, cosb, mHp2, alpha1, cosalpha1, sinalpha1, alpha2, cosalpha2, sinalpha2,
            alpha3, cosalpha3, sinalpha3, Rem12_2, Imlambda5, Relambda6, Imlambda6, Relambda7, Imlambda7, 
            Ytu_11r, Ytu_11i, Ytu_12r, Ytu_12i, Ytu_13r, Ytu_13i, 
            Ytu_21r, Ytu_21i, Ytu_22r, Ytu_22i, Ytu_23r, Ytu_23i, 
            Ytu_31r, Ytu_31i, Ytu_32r, Ytu_32i, Ytu_33r, Ytu_33i, 
            Ytd_11r, Ytd_11i, Ytd_12r, Ytd_12i, Ytd_13r, Ytd_13i, 
            Ytd_21r, Ytd_21i, Ytd_22r, Ytd_22i, Ytd_23r, Ytd_23i, 
            Ytd_31r, Ytd_31i, Ytd_32r, Ytd_32i, Ytd_33r, Ytd_33i, 
            Ytl_11r, Ytl_11i, Ytl_12r, Ytl_12i, Ytl_13r, Ytl_13i, 
            Ytl_21r, Ytl_21i, Ytl_22r, Ytl_22i, Ytl_23r, Ytl_23i, 
            Ytl_31r, Ytl_31i, Ytl_32r, Ytl_32i, Ytl_33r, Ytl_33i, 
            Q_GTHDM;
    bool flag_ATHDM;
};

/**
 * @}
 */

#endif	/* GENERALTHDM_H */

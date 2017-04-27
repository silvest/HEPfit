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

    /**
     * 
     * @return Nu matrix element 11
     */
    gslpp::complex getNu_11() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_11r + i*Nu_11i);
    }
    
    /**
     * 
     * @return Nu matrix element 12
     */
    gslpp::complex getNu_12() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_12r + i*Nu_12i);
    }
    
    /**
     * 
     * @return Nu matrix element 13
     */
    gslpp::complex getNu_13() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_13r + i*Nu_13i);
    }
    
    /**
     * 
     * @return Nu matrix element 21
     */
    gslpp::complex getNu_21() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_21r + i*Nu_21i);
    }
    
    /**
     * 
     * @return Nu matrix element 22
     */
    gslpp::complex getNu_22() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_22r + i*Nu_22i);
    }
    
    /**
     * 
     * @return Nu matrix element 23
     */
    gslpp::complex getNu_23() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_23r + i*Nu_23i);
    }
    
    /**
     * 
     * @return Nu matrix element 31
     */
    gslpp::complex getNu_31() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_31r + i*Nu_31i);
    }
    
    /**
     * 
     * @return Nu matrix element 32
     */
    gslpp::complex getNu_32() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_32r + i*Nu_32i);
    }
    
    /**
     * 
     * @return Nu matrix element 33
     */
    gslpp::complex getNu_33() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nu_33r + i*Nu_33i);
    }
    
    /**
     * 
     * @return Nd matrix element 11
     */
    gslpp::complex getNd_11() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_11r + i*Nd_11i);
    }
    
    /**
     * 
     * @return Nd matrix element 12
     */
    gslpp::complex getNd_12() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_12r + i*Nd_12i);
    }
    
    /**
     * 
     * @return Nd matrix element 13
     */
    gslpp::complex getNd_13() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_13r + i*Nd_13i);
    }
    
    /**
     * 
     * @return Nd matrix element 21
     */
    gslpp::complex getNd_21() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_21r + i*Nd_21i);
    }
    
    /**
     * 
     * @return Nd matrix element 22
     */
    gslpp::complex getNd_22() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_22r + i*Nd_22i);
    }
    
    /**
     * 
     * @return Nd matrix element 23
     */
    gslpp::complex getNd_23() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_23r + i*Nd_23i);
    }
    
    /**
     * 
     * @return Nd matrix element 31
     */
    gslpp::complex getNd_31() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_31r + i*Nd_31i);
    }
    
    /**
     * 
     * @return Nd matrix element 32
     */
    gslpp::complex getNd_32() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_32r + i*Nd_32i);
    }
    
    /**
     * 
     * @return Nd matrix element 33
     */
    gslpp::complex getNd_33() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nd_33r + i*Nd_33i);
    }
    
    /**
     * 
     * @return Nl matrix element 11
     */
    gslpp::complex getNl_11() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_11r + i*Nl_11i);
    }
    
    /**
     * 
     * @return Nl matrix element 12
     */
    gslpp::complex getNl_12() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_12r + i*Nl_12i);
    }
    
    /**
     * 
     * @return Nl matrix element 13
     */
    gslpp::complex getNl_13() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_13r + i*Nl_13i);
    }
    
    /**
     * 
     * @return Nl matrix element 21
     */
    gslpp::complex getNl_21() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_21r + i*Nl_21i);
    }
    
    /**
     * 
     * @return Nl matrix element 22
     */
    gslpp::complex getNl_22() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_22r + i*Nl_22i);
    }
    
    /**
     * 
     * @return Nl matrix element 23
     */
    gslpp::complex getNl_23() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_23r + i*Nl_23i);
    }
    
    /**
     * 
     * @return Nl matrix element 31
     */
    gslpp::complex getNl_31() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_31r + i*Nl_31i);
    }
    
    /**
     * 
     * @return Nl matrix element 32
     */
    gslpp::complex getNl_32() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_32r + i*Nl_32i);
    }
    
    /**
     * 
     * @return Nl matrix element 33
     */
    gslpp::complex getNl_33() const {
        gslpp::complex i = gslpp::complex::i();
        return (Nl_33r + i*Nl_33i);
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
    
    /**
     *
     * @return Flag to switch on the Aligned Two Higgs Doublet Model
     */
    bool getATHDMflag() const {
        return flag_ATHDM;
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
            Nu_11r, Nu_11i, Nu_12r, Nu_12i, Nu_13r, Nu_13i, 
            Nu_21r, Nu_21i, Nu_22r, Nu_22i, Nu_23r, Nu_23i, 
            Nu_31r, Nu_31i, Nu_32r, Nu_32i, Nu_33r, Nu_33i, 
            Nd_11r, Nd_11i, Nd_12r, Nd_12i, Nd_13r, Nd_13i, 
            Nd_21r, Nd_21i, Nd_22r, Nd_22i, Nd_23r, Nd_23i, 
            Nd_31r, Nd_31i, Nd_32r, Nd_32i, Nd_33r, Nd_33i, 
            Nl_11r, Nl_11i, Nl_12r, Nl_12i, Nl_13r, Nl_13i, 
            Nl_21r, Nl_21i, Nl_22r, Nl_22i, Nl_23r, Nl_23i, 
            Nl_31r, Nl_31i, Nl_32r, Nl_32i, Nl_33r, Nl_33i, 
            Q_GTHDM;
    bool flag_ATHDM;
};

/**
 * @}
 */

#endif	/* GENERALTHDM_H */

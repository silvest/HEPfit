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
#include "NPbase.h"

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
class GeneralTHDM: public NPbase {
public:

    static const int NGeneralTHDMvars = 68;
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
     * @brief A get method to retrieve the matching scale @f$\mu_W@f$ around
     * the weak scale.
     * @return @f$\mu_W@f$ in GeV
     */
    double getMuw() const
    {
        return muw;
    }


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

    virtual bool setFlagStr(const std::string name, const std::string value);
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
     * @return Second neutral mass eigenvalue
     */
    double getmH2sq() const {
    if(mH2sq < 0.) {
            throw std::runtime_error("error in GeneralTHDM: mH2sq < 0!");
    }
        return mH2sq;
    }

    /**
     * @return Third neutral mass eigenvalue
     */
    double getmH3sq() const {
    if(mH3sq < 0.) {
            throw std::runtime_error("error in GeneralTHDM: mH3sq < 0!");
    }
        return mH3sq;
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
        if(flag_CPconservation==true) {
            return 0.0;
        }
        else {
            return alpha2;
        }
    }

    /**
     * @return cos of the rotation angle between the first and third neutral mass eigenstate
     */
    double getcosalpha2() const {
        if(flag_CPconservation==true) {
            return 1.0;
        }
        else {
            return cosalpha2;
        }
    }

    /**
     * @return sin of the rotation angle between the first and third neutral mass eigenstate
     */
    double getsinalpha2() const {
        if(flag_CPconservation==true) {
            return 0.0;
        }
        else {
            return sinalpha2;
        }
    }

    /**
     * @return rotation angle between the second and third neutral mass eigenstate
     */
    double getalpha3() const {
        if(flag_CPconservation==true) {
            return 0.0;
        }
        else {
            return alpha3;
        }
    }

    /**
     * @return cos of the rotation angle between the second and third neutral mass eigenstate
     */
    double getcosalpha3() const {
        if(flag_CPconservation==true) {
            return 1.0;
        }
        else {
            return cosalpha3;
        }
    }

    /**
     * @return sin of the rotation angle between the second and third neutral mass eigenstate
     */
    double getsinalpha3() const {
        if(flag_CPconservation==true) {
            return 0.0;
        }
        else {
            return sinalpha3;
        }
    }

    /**
     * @return Real part of the Higgs potential parameter @f$lambda_{5}@f$
     */
    double getRelambda5() const {
        return Relambda5;
    }
    
    /**
     * @return Imaginary part of the Higgs potential parameter @f$lambda_{5}@f$
     */
    double getImlambda5() const {
        if(flag_CPconservation==true) {
            return 0.0;
        }
        else {
            return Imlambda5;
        }
    }
    
    /**
     * @return Real part of the Higgs potential parameter @f$lambda_{6}@f$
     */
    double getRelambda6() const {
        return Relambda6;
    }

    /**
     * @return Real part of the Higgs potential parameter @f$lambda_{7}@f$
     */
    double getRelambda7() const {
        return Relambda7;
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
     * @return Minimal R' value
     */
    double getRpepsGTHDM() const {
        return RpepsGTHDM;
    }

    /**
     *
     * @return Minimal NLO unitarity check scale
     */
    double getNLOuniscaleGTHDM() const {
        return NLOuniscaleGTHDM;
    }

    /**
     *
     * @return Flag to switch on the Aligned Two Higgs Doublet Model
     */
    bool getATHDMflag() const {
        return flag_ATHDM;
    }

    /**
     *
     * @return Flag to switch on CP conservation
     */
    bool getCPconservationflag() const {
        return flag_CPconservation;
    }

    /**
     *
     * @return Switch for NLO RGE and approximate NLO RGE
     */
    std::string getRGEorderflag() const {
        return flag_RGEorder;
    }

    virtual double muggH(const double sqrt_s) const;
    virtual double muVBF(const double sqrt_s) const;
    virtual double mueeWBF(const double sqrt_s) const;
    virtual double muWH(const double sqrt_s) const;
    virtual double muZH(const double sqrt_s) const;
    virtual double mueeZH(const double sqrt_s) const;
    virtual double muVH(const double sqrt_s) const;
    virtual double muVBFpVH(const double sqrt_s) const;
    virtual double muttH(const double sqrt_s) const;
    virtual double GammaTotal() const;
    virtual double BrHggRatio() const;
    virtual double BrHWWRatio() const;
    virtual double BrHZZRatio() const;
    virtual double BrHZgaRatio() const;
    virtual double BrHgagaRatio() const;
    virtual double BrHmumuRatio() const;
    virtual double BrHtautauRatio() const;
    virtual double BrHccRatio() const;
    virtual double BrHbbRatio() const;
    virtual double muggHgaga(const double sqrt_s) const;
    virtual double muVBFHgaga(const double sqrt_s) const;
    virtual double muVHgaga(const double sqrt_s) const;
    virtual double muttHgaga(const double sqrt_s) const;
    virtual double muggHZZ(const double sqrt_s) const;
    virtual double muVBFHZZ(const double sqrt_s) const;
    virtual double muVHZZ(const double sqrt_s) const;
    virtual double muttHZZ(const double sqrt_s) const;
    virtual double muggHWW(const double sqrt_s) const;
    virtual double muVBFHWW(const double sqrt_s) const;
    virtual double muVHWW(const double sqrt_s) const;
    virtual double muttHWW(const double sqrt_s) const;
    virtual double muggHtautau(const double sqrt_s) const;
    virtual double muVBFHtautau(const double sqrt_s) const;
    virtual double muVHtautau(const double sqrt_s) const;
    virtual double muttHtautau(const double sqrt_s) const;
    virtual double muggHbb(const double sqrt_s) const;
    virtual double muVBFHbb(const double sqrt_s) const;
    virtual double muVHbb(const double sqrt_s) const;
    virtual double muttHbb(const double sqrt_s) const;
    virtual double muppHmumu(const double sqrt_s) const;
    virtual double muppHZga(const double sqrt_s) const;
    virtual double computeGammaTotalRatio() const;

    /*
     * @brief Method overriding the NPbase Mw with the tree-level W mass
     */
    virtual double Mw() const;

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

    /**
     * @brief A method to check if the RGE order name in string form is valid.
     * @param[in] GeneralTHDM RGE order
     * @return a boolean that is true if the RGE order string is valid
     */
    bool checkRGEorder(const std::string RGEorder) const
    {
        if (RGEorder.compare("LO") == 0)
//                || RGEorder.compare("approxNLO") == 0)
            return true;
        else
            return false;
    }

    mutable Matching<GeneralTHDMMatching,GeneralTHDM> GTHDMM; ///< An object of type Matching.

private:

    GeneralTHDMcache* myGTHDMcache;

    double logtb, tanb, sinb, cosb, mHp2, mH2sq, mH3sq, alpha1, cosalpha1, sinalpha1, alpha2, cosalpha2, sinalpha2,
            alpha3, cosalpha3, sinalpha3, Relambda5, Imlambda5, Relambda6, Relambda7,
            Nu_11r, Nu_11i, Nu_12r, Nu_12i, Nu_13r, Nu_13i, 
            Nu_21r, Nu_21i, Nu_22r, Nu_22i, Nu_23r, Nu_23i, 
            Nu_31r, Nu_31i, Nu_32r, Nu_32i, Nu_33r, Nu_33i, 
            Nd_11r, Nd_11i, Nd_12r, Nd_12i, Nd_13r, Nd_13i, 
            Nd_21r, Nd_21i, Nd_22r, Nd_22i, Nd_23r, Nd_23i, 
            Nd_31r, Nd_31i, Nd_32r, Nd_32i, Nd_33r, Nd_33i, 
            Nl_11r, Nl_11i, Nl_12r, Nl_12i, Nl_13r, Nl_13i, 
            Nl_21r, Nl_21i, Nl_22r, Nl_22i, Nl_23r, Nl_23i, 
            Nl_31r, Nl_31i, Nl_32r, Nl_32i, Nl_33r, Nl_33i, 
            Q_GTHDM, RpepsGTHDM, NLOuniscaleGTHDM;
    bool flag_ATHDM, flag_CPconservation;
    std::string flag_RGEorder;
};

/**
 * @}
 */

#endif	/* GENERALTHDM_H */

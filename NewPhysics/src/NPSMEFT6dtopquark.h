/*
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef NPSMEFT6DTOPQUARK_H
#define NPSMEFT6DTOPQUARK_H


#include "gslpp.h"
#include "ThObservable.h"
#include "NPbase.h"
#include <stdexcept>

class NPSMEFT6dtopquark : public NPbase {
public:
    
    
/**
*　@brief The number of the model parameters in %NPSMEFT6dtopquarkVars. 
*/
static const int NNPSMEFT6dtopquarkVars = 67;
   
/*
@brief A string array containing the labels of the model parameters in NPSMEFT6dtopquark 
*/
    static const std::string NPSMEFT6dtopquarkVars[NNPSMEFT6dtopquarkVars];
    
    
   NPSMEFT6dtopquark();
   
   
       /**
     * @brief The parameter \f$\varepsilon_1\f$.
     * @return the SM value (FlagEpsilon1SM=true) or the SM plus new physics
     * value (FlagEpsilon1SM=false) of \f$\varepsilon_1\f$
     */
    //virtual double myC_phit() const;
    
    /**
     *
     * @return @f$\log_{10}(\tan \beta)@f$
     */
    
    double getNPSMEFT6dtopquark_C_phit() const {
        return C_phit;
    }
    
    double getNPSMEFT6dtopquark_C_phiQ3() const {
        return C_phiQ3;
    }
    
    double getNPSMEFT6dtopquark_C_phiQ1() const {
        return C_phiQ1;
    }
    
    double getNPSMEFT6dtopquark_C_phiQm() const {
        return C_phiQm;
    }
    
    double getNPSMEFT6dtopquark_C_tW() const {
        return C_tW;
    }
    
    double getNPSMEFT6dtopquark_C_tZ() const {
        return C_tZ;
    }
    
    double getNPSMEFT6dtopquark_C_tB() const {
        return C_tB;
    }
    
    double getNPSMEFT6dtopquark_C_tphi() const {
        return C_tphi;
    }
    
    double getNPSMEFT6dtopquark_C_phib() const {
        return C_phib;
    }
    
    double getNPSMEFT6dtopquark_C_bW() const {
        return C_bW;
    }
    
    double getNPSMEFT6dtopquark_C_bB() const {
        return C_bB;
    }
    
    double getNPSMEFT6dtopquark_C_bZ() const {
        return C_bZ;
    }
    
    double getNPSMEFT6dtopquark_C_phitb() const {
        return C_phitb;
    }
    
    double getNPSMEFT6dtopquark_C_tG() const {
        return C_tG;
    }
    
    double getNPSMEFT6dtopquark_C_ed() const {
        return C_ed;
    }
    
    double getNPSMEFT6dtopquark_C_eq() const {
        return C_eq;
    }
    
    double getNPSMEFT6dtopquark_C_ld() const {
        return C_ld;
    }
    
    double getNPSMEFT6dtopquark_C_lqP() const {
        return C_lqP;
    }
    
    double getNPSMEFT6dtopquark_C_eu() const {
        return C_eu;
    }
    
    double getNPSMEFT6dtopquark_C_lu() const {
        return C_lu;
    }
    
    double getNPSMEFT6dtopquark_C_lqM() const {
        return C_lqM;
    }
    
    
    double getNPSMEFT6dtopquark_C_tu8() const {
        return C_tu8;
    }

    double getNPSMEFT6dtopquark_C_td8() const {
        return C_td8;
    }
    double getNPSMEFT6dtopquark_C_Qq18() const {
        return C_Qq18;
    }
    double getNPSMEFT6dtopquark_C_tq8() const {
        return C_tq8;
    }
    double getNPSMEFT6dtopquark_C_Qq38() const {
        return C_Qq38;
    }
    double getNPSMEFT6dtopquark_C_Qu8() const {
        return C_Qu8;
    }
    double getNPSMEFT6dtopquark_C_Qd8() const {
        return C_Qd8;
    }
    
    
    
     double getNPSMEFT6dtopquark_SM_tAq_inc() const {
        return SM_tAq_inc;
    }
     
    double getNPSMEFT6dtopquark_SM_ttZ_bin_0_40() const {
        return SM_ttZ_bin_0_40;
    }
    
    double getNPSMEFT6dtopquark_SM_ttZ_bin_40_70() const {
        return SM_ttZ_bin_40_70;
    }
    
    double getNPSMEFT6dtopquark_SM_ttZ_bin_70_110() const {
        return SM_ttZ_bin_70_110;
    }
    
    double getNPSMEFT6dtopquark_SM_ttZ_bin_110_160() const {
        return SM_ttZ_bin_110_160;
    }
    
    double getNPSMEFT6dtopquark_SM_ttZ_bin_160_220() const {
        return SM_ttZ_bin_160_220;
    }
    
    double getNPSMEFT6dtopquark_SM_ttZ_bin_220_290() const {
        return SM_ttZ_bin_220_290;
    }
    
    double getNPSMEFT6dtopquark_SM_ttZ_bin_290_400() const {
        return SM_ttZ_bin_290_400;
    }
        
    
    
    
    
    
    
    
    
    double getNPSMEFT6dtopquark_SM_ttA_bin_20_25() const {
        return SM_ttA_bin_20_25;
    }
    double getNPSMEFT6dtopquark_SM_ttA_bin_25_30() const {
        return SM_ttA_bin_25_30;
    }
    double getNPSMEFT6dtopquark_SM_ttA_bin_30_35() const {
        return SM_ttA_bin_30_35;
    }
    double getNPSMEFT6dtopquark_SM_ttA_bin_35_40() const {
        return SM_ttA_bin_35_40;
    }
    double getNPSMEFT6dtopquark_SM_ttA_bin_40_47() const {
        return SM_ttA_bin_40_47;
    }
    double getNPSMEFT6dtopquark_SM_ttA_bin_47_55() const {
        return SM_ttA_bin_47_55;
    }
    double getNPSMEFT6dtopquark_SM_ttA_bin_55_70() const {
        return SM_ttA_bin_55_70;
    }
    double getNPSMEFT6dtopquark_SM_ttA_bin_70_85() const {
        return SM_ttA_bin_70_85;
    }
    double getNPSMEFT6dtopquark_SM_ttA_bin_85_132() const {
        return SM_ttA_bin_85_132;
    }
    double getNPSMEFT6dtopquark_SM_ttA_bin_132_180() const {
        return SM_ttA_bin_132_180;
    }
    double getNPSMEFT6dtopquark_SM_ttA_bin_180_300() const {
        return SM_ttA_bin_180_300;
    }
    
    
    
    
    
    
    double getNPSMEFT6dtopquark_SM_tZQ_inc() const {
        return SM_tZQ_inc;
    }
    
    double getNPSMEFT6dtopquark_SM_ttA_inc() const {
        return SM_ttA_inc;
    }
   
    double getNPSMEFT6dtopquark_SM_ttZ_inc() const {
        return SM_ttZ_inc;
    }
    
    double getNPSMEFT6dtopquark_SM_ttH_inc() const {
        return SM_ttH_inc;
    }
    
    double getNPSMEFT6dtopquark_SM_ttW_inc() const {
        return SM_ttW_inc;
    }
    
     double getNPSMEFT6dtopquark_SM_tW_inc() const {
        return SM_tW_inc;
    }
      double getNPSMEFT6dtopquark_SM_tW_inc_8TeV() const {
        return SM_tW_inc_8TeV;
    }
        double getNPSMEFT6dtopquark_SM_sigmatchannel13() const{
        return SM_sigmatchannel13;
    }
        
        double getNPSMEFT6dtopquark_SM_sigmatchannel8() const{
        return SM_sigmatchannel8;
    }
        
        double getNPSMEFT6dtopquark_SM_sigmaschannel8() const{
        return SM_sigmaschannel8;
    }
        
        double getNPSMEFT6dtopquark_SM_sigmaschannelTev() const{
        return SM_sigmaschannelTev;
    }
        
        double getNPSMEFT6dtopquark_SM_tH_tchan_value() const{
        return SM_tH_tchan_value;
    }
      
    double getNPSMEFT6dtopquark_flag_Quadratic() const {
        return flag_Quadratic;
    }
 
    double getNPSMEFT6dtopquark_flag_LHC_WG_Basis() const {
        return flag_LHC_WG_Basis;
    }
    
    double getNPSMEFT6dtopquark_SM_ttbar_LHC13_value() const{
        return SM_ttbar_LHC13;
    }
	
    double getNPSMEFT6dtopquark_SM_ttbar_LHC8_value() const{
        return SM_ttbar_LHC8;
    }
	
    double getNPSMEFT6dtopquark_SM_ttbar_Tev_value() const{
        return SM_ttbar_Tev;
    }
    
     double getNPSMEFT6dtopquark_F0_SM() const{
        return F0_SM;
    }
      double getNPSMEFT6dtopquark_FL_SM() const{
        return FL_SM;
    }
      
    double getNPSMEFT6dtopquark_Rb_SM() const{
        return Rb_SM;
    }
     
    double getNPSMEFT6dtopquark_AFBLR_SM() const{
        return AFBLR_SM;
    }
     
    double getNPSMEFT6dtopquark_ttWqEM_SM() const{
        return ttWqEM_SM;
    }



    /////////////////////////////////////////////////////////////////
    protected:
    virtual void setParameter(const std::string name, const double& value);

    double C_phit; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phiQ3; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phiQ1; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phiQm; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tW; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tZ; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tB; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tphi; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phib; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_bW; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_bB; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_bZ; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phitb; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tG; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_ed; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_eq; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_ld; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_lqP; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_eu; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_lu; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_lqM; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tu8; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_td8; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_Qq18; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tq8; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_Qq38; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_Qu8; //< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_Qd8; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double SM_tAq_inc; ///< The SM value for tgammaq
    double SM_ttZ_bin_0_40; ///< The SM value for the differential cross section ttZ in the bin from 0 to 40 GeV.
    double SM_ttZ_bin_40_70; ///< The SM value for the differential cross section ttZ in the bin from 40 to 70 GeV.
    double SM_ttZ_bin_70_110; ///< The SM value for the differential cross section ttZ in the bin from 70 to 110 GeV.
    double SM_ttZ_bin_110_160; ///< The SM value for the differential cross section ttZ in the bin from 110 to 160 GeV.
    double SM_ttZ_bin_160_220; ///< The SM value for the differential cross section ttZ in the bin from 160 to 220 GeV.
    double SM_ttZ_bin_220_290; ///< The SM value for the differential cross section ttZ in the bin from 220 to 290 GeV.
    double SM_ttZ_bin_290_400; ///< The SM value for the differential cross section ttZ in the bin from 290 to 400 GeV.

    double SM_ttA_bin_20_25; ///< The SM value for the differential cross section ttA in the bin from 20 to 25 GeV.
    double SM_ttA_bin_25_30; ///< The SM value for the differential cross section ttA in the bin from 25 to 30 GeV.
    double SM_ttA_bin_30_35; ///< The SM value for the differential cross section ttA in the bin from 30 to 35 GeV.
    double SM_ttA_bin_35_40; ///< The SM value for the differential cross section ttA in the bin from 35 to 40 GeV.
    double SM_ttA_bin_40_47; ///< The SM value for the differential cross section ttA in the bin from 40 to 47 GeV.
    double SM_ttA_bin_47_55; ///< The SM value for the differential cross section ttA in the bin from 47 to 55 GeV.
    double SM_ttA_bin_55_70; ///< The SM value for the differential cross section ttA in the bin from 55 to 70 GeV.
    double SM_ttA_bin_70_85; ///< The SM value for the differential cross section ttA in the bin from 70 to 85 GeV.
    double SM_ttA_bin_85_132; ///< The SM value for the differential cross section ttA in the bin from 85 to 132 GeV.
    double SM_ttA_bin_132_180; ///< The SM value for the differential cross section ttA in the bin from 132 to 180 GeV.
    double SM_ttA_bin_180_300; ///< The SM value for the differential cross section ttA in the bin from 180 to 300 GeV.



    double SM_tZQ_inc; ///< The SM value for the inclusive cross section tZQ.
    double SM_ttA_inc; ///< The SM value for the inclusive cross section ttA.
    double SM_ttZ_inc; ///< The SM value for the inclusive cross section ttZ.
    double SM_ttH_inc; ///< The SM value for the inclusive cross section ttH.
    double SM_ttW_inc; ///< The SM value for the inclusive cross section ttW.
    double SM_tW_inc; ///< The SM value for the inclusive cross section tW.
    double SM_tW_inc_8TeV; ///< The SM value for the inclusive cross section tW.
    double SM_sigmatchannel13; ///< The SM value for the inclusive cross section single top production at 13 TeV in the t-channel.
    double SM_sigmatchannel8; ///< The SM value for the inclusive cross section single top production at 8 TeV in the t-channel.
    double SM_sigmaschannel8; ///< The SM value for the inclusive cross section single top production at 8 TeV in the s-channel.
    double SM_sigmaschannelTev; ///< The SM value for the inclusive cross section single top production at 1.96 TeV in the s-channel.

    double SM_tH_tchan_value; ///< The SM value for the inclusive cross section single top production at 8 TeV in the s-channel.

    
    double SM_ttbar_LHC13;
    double SM_ttbar_LHC8;
    double SM_ttbar_Tev;
    
    double F0_SM;
    double FL_SM;
    
    double Rb_SM;
    double AFBLR_SM;
    
    double ttWqEM_SM;
        
    virtual bool setFlag(const std::string, const bool);

        
////////////////////////////////////////////////////////////////////////
private:

    bool flag_LHC_WG_Basis;
    bool flag_Quadratic;
};


/**
 * @class C_phit
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_phit: public ThObservable {
public:

    /**
     * @brief C_phit constructor.
     */
    C_phit(const StandardModel& SM_i);

    /**
     * @return The value of C_phit
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_phiQ3
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phiQ3
 */
class C_phiQ3: public ThObservable {
public:

    /**
     * @brief C_phiQ3 constructor.
     */
    C_phiQ3(const StandardModel& SM_i);

    /**
     * @return The value of C_phiQ3
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_phiQ1
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phiQ1
 */
class C_phiQ1: public ThObservable {
public:

    /**
     * @brief C_phiQ1 constructor.
     */
    C_phiQ1(const StandardModel& SM_i);

    /**
     * @return The value of C_phiQ1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_phiQm
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_phiQm: public ThObservable {
public:

    /**
     * @brief C_phit constructor.
     */
    C_phiQm(const StandardModel& SM_i);

    /**
     * @return The value of C_phiQm
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_tW
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_tW
 */
class C_tW: public ThObservable {
public:

    /**
     * @brief C_phit constructor.
     */
    C_tW(const StandardModel& SM_i);

    /**
     * @return The value of C_tW
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_tZ
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_tZ: public ThObservable {
public:

    /**
     * @brief C_tZ constructor.
     */
    C_tZ(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_tB
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_tB
 */
class C_tB: public ThObservable {
public:

    /**
     * @brief C_tB constructor.
     */
    C_tB(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_tphi
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_tphi
 */
class C_tphi: public ThObservable {
public:

    /**
     * @brief C_phit constructor.
     */
    C_tphi(const StandardModel& SM_i);

    /**
     * @return The value of C_tphi
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



/**
 * @class C_phib
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_phib: public ThObservable {
public:

    /**
     * @brief C_phib constructor.
     */
    C_phib(const StandardModel& SM_i);

    /**
     * @return The value of C_phib
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_bW
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_bW
 */
class C_bW: public ThObservable {
public:

    /**
     * @brief C_bW constructor.
     */
    C_bW(const StandardModel& SM_i);

    /**
     * @return The value of C_bW
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_bB
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_bB
 */
class C_bB: public ThObservable {
public:

    /**
     * @brief C_bB constructor.
     */
    C_bB(const StandardModel& SM_i);

    /**
     * @return The value of C_bB
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

/**
 * @class C_bZ
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_bZ
 */
class C_bZ: public ThObservable {
public:

    /**
     * @brief C_bZ constructor.
     */
    C_bZ(const StandardModel& SM_i);

    /**
     * @return The value of C_bZ
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_tG
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_tG
 */
class C_tG: public ThObservable {
public:

    /**
     * @brief C_tG constructor.
     */
    C_tG(const StandardModel& SM_i);

    /**
     * @return The value of C_tG
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_phitb
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phitb
 */
class C_phitb: public ThObservable {
public:

    /**
     * @brief C_phitb constructor.
     */
    C_phitb(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_ed
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_ed
 */
class C_ed: public ThObservable {
public:

    /**
     * @brief C_ed constructor.
     */
    C_ed(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_eq
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_eq
 */
class C_eq: public ThObservable {
public:

    /**
     * @brief C_eq constructor.
     */
    C_eq(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_ld
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_ld
 */
class C_ld: public ThObservable {
public:

    /**
     * @brief C_ld constructor.
     */
    C_ld(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_lqP
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_lqP
 */
class C_lqP: public ThObservable {
public:

    /**
     * @brief C_lqP constructor.
     */
    C_lqP(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_eu
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_eu
 */
class C_eu: public ThObservable {
public:

    /**
     * @brief C_eu constructor.
     */
    C_eu(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_lu
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_lu
 */
class C_lu: public ThObservable {
public:

    /**
     * @brief C_lu constructor.
     */
    C_lu(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_lqM
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_lqM
 */
class C_lqM: public ThObservable {
public:

    /**
     * @brief C_lqP constructor.
     */
    C_lqM(const StandardModel& SM_i);

    /**
     * @return The value of mH1
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


/**
 * @class C_tu8
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_tu8
 */
class C_tu8: public ThObservable {
public:

    /**
     * @brief C_tu8 constructor.
     */
    C_tu8(const StandardModel& SM_i);

    /**
     * @return The value of C_tu8
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

/**
 * @class C_td8
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_td8: public ThObservable {
public:

    /**
     * @brief C_td8 constructor.
     */
    C_td8(const StandardModel& SM_i);

    /**
     * @return The value of C_td8
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

/**
 * @class C_Qq18
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_Qq18: public ThObservable {
public:

    /**
     * @brief C_Qq18 constructor.
     */
    C_Qq18(const StandardModel& SM_i);

    /**
     * @return The value of C_Qq18
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

/**
 * @class C_tq8
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_tq8: public ThObservable {
public:

    /**
     * @brief C_tq8 constructor.
     */
    C_tq8(const StandardModel& SM_i);

    /**
     * @return The value of C_tq8
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

/**
 * @class C_Qq38
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_Qq38: public ThObservable {
public:

    /**
     * @brief C_Qq38 constructor.
     */
    C_Qq38(const StandardModel& SM_i);

    /**
     * @return The value of C_Qq38
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

/**
 * @class C_Qu8
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_Qu8: public ThObservable {
public:

    /**
     * @brief C_Qu8 constructor.
     */
    C_Qu8(const StandardModel& SM_i);

    /**
     * @return The value of C_Qu8
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

/**
 * @class C_Qd8
 * @ingroup NPSMEFT6dtopquark 
 * @brief Wilson coefficient C_phit
 */
class C_Qd8: public ThObservable {
public:

    /**
     * @brief C_Qd8 constructor.
     */
    C_Qd8(const StandardModel& SM_i);

    /**
     * @return The value of C_Qd8
     */
    double computeThValue();

private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};




//Observables from LEP

class Rb_NPSMEFT6dtopquark : public ThObservable {
public:   

    Rb_NPSMEFT6dtopquark(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class AFBLR : public ThObservable {
public:   

    AFBLR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


//Observables Tevatron



class sigmattbarTev : public ThObservable {
public:   

    sigmattbarTev(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};




class sigmaschannelTev : public ThObservable {
public:   

    sigmaschannelTev(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};





//Observables from LHC



class F0 : public ThObservable {
public:   

    F0(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class FL : public ThObservable {
public:   

    FL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};




class sigmattbarLHC13 : public ThObservable {
public:   

    sigmattbarLHC13(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



class sigmattbarLHC8 : public ThObservable {
public:   

    sigmattbarLHC8(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};




class sigmattZ : public ThObservable {
public:   

    sigmattZ(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigmattA : public ThObservable {
public:   

    sigmattA(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class sigmattH : public ThObservable {
public:   

    sigmattH(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class sigmattW : public ThObservable {
public:   

    sigmattW(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class sigmatchannel13 : public ThObservable {
public:   

    sigmatchannel13(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigmatchannel8 : public ThObservable {
public:   

    sigmatchannel8(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigmaschannel8 : public ThObservable {
public:   

    sigmaschannel8(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigmatW : public ThObservable {
public:   

    sigmatW(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class sigmatW_8TeV : public ThObservable {
public:   

    sigmatW_8TeV(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class sigmatqZ : public ThObservable {
public:   

    sigmatqZ(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


//sigma tAq

class sigmatAq : public ThObservable {
public:   

    sigmatAq(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class tH_tchan : public ThObservable {
public:   

    tH_tchan(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class ttHSUM : public ThObservable {
public:   

    ttHSUM(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class ttWqEM : public ThObservable {
public:   

    ttWqEM(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class ttWqSUM : public ThObservable {
public:   

    ttWqSUM(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

// ttZ differential cross section for different bins



class ttZ_bin_0_40 : public ThObservable {
public:   

    ttZ_bin_0_40(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class ttZ_bin_40_70 : public ThObservable {
public:   

    ttZ_bin_40_70(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



class ttZ_bin_70_110 : public ThObservable {
public:   

    ttZ_bin_70_110(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



class ttZ_bin_110_160 : public ThObservable {
public:   

    ttZ_bin_110_160(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class ttZ_bin_160_220 : public ThObservable {
public:   

    ttZ_bin_160_220(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class ttZ_bin_220_290 : public ThObservable {
public:   

    ttZ_bin_220_290(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class ttZ_bin_290_400 : public ThObservable {
public:   

    ttZ_bin_290_400(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};






// ttA differential cross section for different bins


class ttA_bin_20_25 : public ThObservable {
public:   

    ttA_bin_20_25(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class ttA_bin_25_30 : public ThObservable {
public:   

    ttA_bin_25_30(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class ttA_bin_30_35 : public ThObservable {
public:   

    ttA_bin_30_35(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class ttA_bin_35_40 : public ThObservable {
public:   

    ttA_bin_35_40(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class ttA_bin_40_47 : public ThObservable {
public:   

    ttA_bin_40_47(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class ttA_bin_47_55 : public ThObservable {
public:   

    ttA_bin_47_55(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class ttA_bin_55_70 : public ThObservable {
public:   

    ttA_bin_55_70(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class ttA_bin_70_85 : public ThObservable {
public:   

    ttA_bin_70_85(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class ttA_bin_85_132 : public ThObservable {
public:   

    ttA_bin_85_132(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class ttA_bin_132_180 : public ThObservable {
public:   

    ttA_bin_132_180(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class ttA_bin_180_300 : public ThObservable {
public:   

    ttA_bin_180_300(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};






//Prospects of Linear Collider at 250 GeV
//250 bb observables


class sigma_250_bb_eLpR : public ThObservable {
public:   

    sigma_250_bb_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class a_250_bb_eLpR : public ThObservable {
public:   

    a_250_bb_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigma_250_bb_eRpL : public ThObservable {
public:   

    sigma_250_bb_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class a_250_bb_eRpL : public ThObservable {
public:   

    a_250_bb_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

//Prospects of Linear Colders at 500 GeV
//500 bb observables


class sigma_500_bb_eLpR : public ThObservable {
public:   

    sigma_500_bb_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class a_500_bb_eLpR : public ThObservable {
public:   

    a_500_bb_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigma_500_bb_eRpL : public ThObservable {
public:   

    sigma_500_bb_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class a_500_bb_eRpL : public ThObservable {
public:   

    a_500_bb_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};




//500 tt observables
//Should not be used in combination with 500 tt optimal observables (they are redundant)

class sigma_500_tt_eLpR : public ThObservable {
public:   

    sigma_500_tt_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class a_500_tt_eLpR : public ThObservable {
public:   

    a_500_tt_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigma_500_tt_eRpL : public ThObservable {
public:   

    sigma_500_tt_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class a_500_tt_eRpL : public ThObservable {
public:   

    a_500_tt_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class pt_500_tt_eLpR : public ThObservable {
public:   

    pt_500_tt_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class pt_500_tt_eRpL : public ThObservable {
public:   

    pt_500_tt_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


// Optimal observables for 500 GeV linear collider

class op1 : public ThObservable {
public:   

    op1(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class op2 : public ThObservable {
public:   

    op2(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



class op3 : public ThObservable {
public:   

    op3(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class op4 : public ThObservable {
public:   

    op4(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};








//Prospects of Linear Colders at 1000 GeV
//1000 bb observables


class sigma_1000_bb_eLpR : public ThObservable {
public:   

    sigma_1000_bb_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class a_1000_bb_eLpR : public ThObservable {
public:   

    a_1000_bb_eLpR(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigma_1000_bb_eRpL : public ThObservable {
public:   

    sigma_1000_bb_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class a_1000_bb_eRpL : public ThObservable {
public:   

    a_1000_bb_eRpL(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



// Optimal observables for 1000 GeV linear collider


class op_1000_1 : public ThObservable {
public:   

    op_1000_1(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class op_1000_2 : public ThObservable {
public:   

    op_1000_2(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class op_1000_3 : public ThObservable {
public:   

    op_1000_3(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class op_1000_4 : public ThObservable {
public:   

    op_1000_4(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class op_1000_5 : public ThObservable {
public:   

    op_1000_5(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class op_1000_6 : public ThObservable {
public:   

    op_1000_6(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class op_1000_7 : public ThObservable {
public:   

    op_1000_7(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class op_1000_8 : public ThObservable {
public:   

    op_1000_8(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


// Relation with other basis

class gLt : public ThObservable {
public:   

    gLt(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class gLb : public ThObservable {
public:   

    gLb(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};

class gRt : public ThObservable {
public:   

    gRt(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class gRb : public ThObservable {
public:   

    gRb(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


#endif /* NPSMEFT6DTOPQUARK_H */

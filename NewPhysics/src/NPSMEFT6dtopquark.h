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
//static const int NNPSMEFT6dtopquarkVars = 104;
static const int NNPSMEFT6dtopquarkVars = 35 + 7;
  

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
    
    
    
    
    
    
    
    
    double getNPSMEFT6dtopquark_C_Qd1() const {
        return C_Qd1;
    }
    double getNPSMEFT6dtopquark_C_Qu1() const {
        return C_Qu1;
    }
    double getNPSMEFT6dtopquark_C_td1() const {
        return C_td1;
    }
    double getNPSMEFT6dtopquark_C_tu1() const {
        return C_tu1;
    }
    double getNPSMEFT6dtopquark_C_tq1() const {
        return C_tq1;
    }
    double getNPSMEFT6dtopquark_C_Qq11() const {
        return C_Qq11;
    }
    double getNPSMEFT6dtopquark_C_Qq31() const {
        return C_Qq31;
    }
    
    
    
    double getNPSMEFT6dtopquark_C_phiG() const {
        return C_phiG;
    }
    double getNPSMEFT6dtopquark_C_phiGtil() const {
        return C_phiGtil;
    }
    double getNPSMEFT6dtopquark_C_phiW() const {
        return C_phiW;
    }
    double getNPSMEFT6dtopquark_C_phiWtil() const {
        return C_phiWtil;
    }
    double getNPSMEFT6dtopquark_C_tphiIm() const {
        return C_tphiIm;
    }
    double getNPSMEFT6dtopquark_C_tGIm() const {
        return C_tGIm;
    }
    double getNPSMEFT6dtopquark_C_tWIm() const {
        return C_tWIm;
    }

    
    
    double getNPSMEFT6dtopquark_flag_Quadratic() const {
        return flag_Quadratic;
    }
 
    double getNPSMEFT6dtopquark_flag_LHC_WG_Basis() const {
        return flag_LHC_WG_Basis;
    }
    
    
    /*
    
    
    
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
    
    
    
    
    double getNPSMEFT6dtopquark_SM_tt_bin_250_400() const {
        return SM_tt_bin_250_400;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_400_480() const {
        return SM_tt_bin_400_480;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_480_560() const {
        return SM_tt_bin_480_560;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_560_640() const {
        return SM_tt_bin_560_640;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_640_720() const {
        return SM_tt_bin_640_720;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_720_800() const {
        return SM_tt_bin_720_800;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_800_900() const {
        return SM_tt_bin_800_900;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_900_1000() const {
        return SM_tt_bin_900_1000;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_1000_1150() const {
        return SM_tt_bin_1000_1150;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_1150_1300() const {
        return SM_tt_bin_1150_1300;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_1300_1500() const {
        return SM_tt_bin_1300_1500;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_1500_1700() const {
        return SM_tt_bin_1500_1700;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_1700_2000() const {
        return SM_tt_bin_1700_2000;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_2000_2300() const {
        return SM_tt_bin_2000_2300;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_2300_2600() const {
        return SM_tt_bin_2300_2600;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_2600_3000() const {
        return SM_tt_bin_2600_3000;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_3000_3500() const {
        return SM_tt_bin_3000_3500;
    }
        
    double getNPSMEFT6dtopquark_SM_tt_bin_3500_4000() const {
        return SM_tt_bin_3500_4000;
    }
        
    
    
    
    
    
    
    double getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_0_500() const {
        return SM_Charge_Asymmetry_bin_tt_0_500;
    }
    
    double getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_500_750() const {
        return SM_Charge_Asymmetry_bin_tt_500_750;
    }
    
    double getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_750_1000() const {
        return SM_Charge_Asymmetry_bin_tt_750_1000;
    }
    
    double getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_1000_1500() const {
        return SM_Charge_Asymmetry_bin_tt_1000_1500;
    }
    
    double getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_1500_2000() const {
        return SM_Charge_Asymmetry_bin_tt_1500_2000;
    }
    
    double getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_2000_2500() const {
        return SM_Charge_Asymmetry_bin_tt_2000_2500;
    }
    
    double getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_2500_3000() const {
        return SM_Charge_Asymmetry_bin_tt_2500_3000;
    }
    
    
    
    
    
    double getNPSMEFT6dtopquark_SM_ttll_bin_100_120() const {
        return SM_ttll_bin_100_120;
    }
    double getNPSMEFT6dtopquark_SM_ttll_bin_120_140() const {
        return SM_ttll_bin_120_140;
    }
    double getNPSMEFT6dtopquark_SM_ttll_bin_140_180() const {
        return SM_ttll_bin_140_180;
    }
    double getNPSMEFT6dtopquark_SM_ttll_bin_180_500() const {
        return SM_ttll_bin_180_500;
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
    
    double getNPSMEFT6dtopquark_SM_Asymmetry_leptonic_charge_rapidity_ttW() const {
        return SM_Asymmetry_leptonic_charge_rapidity_ttW;
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
    */



    /////////////////////////////////////////////////////////////////
    protected:
    virtual void setParameter(const std::string name, const double& value);

    
    
    
    double C_phit = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phiQ3 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phiQ1 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phiQm = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tW = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tZ = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tB = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tphi = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phib = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_bW = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_bB = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_bZ = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phitb = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tG = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_ed = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_eq = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_ld = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_lqP = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_eu = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_lu = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_lqM = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tu8 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_td8 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_Qq18 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tq8 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_Qq38 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_Qu8 = 0; //< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_Qd8 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    
    double C_Qd1 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_Qu1 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_td1 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tu1 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tq1 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_Qq11 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_Qq31 = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    
    
    double C_phiG = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phiGtil = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phiW = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_phiWtil = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tphiIm = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tGIm = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    double C_tWIm = 0; ///< The dimension-6 operator coefficient \f$C_{G}\f$.
    
    
    /*
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

    
    
    double SM_tt_bin_250_400; ///< The SM value for the differential cross section tt in the bin from 250 to 400 GeV.
    double SM_tt_bin_400_480; ///< The SM value for the differential cross section tt in the bin from  400 to 480 GeV.
    double SM_tt_bin_480_560; ///< The SM value for the differential cross section tt in the bin from  480 to 560 GeV.
    double SM_tt_bin_560_640; ///< The SM value for the differential cross section tt in the bin from  560 to 640 GeV.
    double SM_tt_bin_640_720; ///< The SM value for the differential cross section tt in the bin from  640 to 720 GeV.
    double SM_tt_bin_720_800; ///< The SM value for the differential cross section tt in the bin from  720 to 800 GeV.
    double SM_tt_bin_800_900; ///< The SM value for the differential cross section tt in the bin from  800 to 900 GeV.
    double SM_tt_bin_900_1000; ///< The SM value for the differential cross section tt in the bin from  900 to 1000 GeV.
    double SM_tt_bin_1000_1150; ///< The SM value for the differential cross section tt in the bin from  1000 to 1150 GeV.
    double SM_tt_bin_1150_1300; ///< The SM value for the differential cross section tt in the bin from  1150 to 1300 GeV.
    double SM_tt_bin_1300_1500; ///< The SM value for the differential cross section tt in the bin from  1300 to 1500 GeV.
    double SM_tt_bin_1500_1700; ///< The SM value for the differential cross section tt in the bin from  1500 to 1700 GeV.
    double SM_tt_bin_1700_2000; ///< The SM value for the differential cross section tt in the bin from  1700 to 2000 GeV.
    double SM_tt_bin_2000_2300; ///< The SM value for the differential cross section tt in the bin from  2000 to 2300 GeV.
    double SM_tt_bin_2300_2600; ///< The SM value for the differential cross section tt in the bin from  2300 to 2600 GeV.
    double SM_tt_bin_2600_3000; ///< The SM value for the differential cross section tt in the bin from  2600 to 3000 GeV.
    double SM_tt_bin_3000_3500; ///< The SM value for the differential cross section tt in the bin from  3000 to 3500 GeV.
    double SM_tt_bin_3500_4000; ///< The SM value for the differential cross section tt in the bin from  3500 to 4000 GeV.
    
    double SM_Charge_Asymmetry_bin_tt_0_500;     ///< The SM value for the charged asymmetry of tt differential in the bin from  0 to 500 GeV.
    double SM_Charge_Asymmetry_bin_tt_500_750;   ///< The SM value for the charged asymmetry of tt differential in the bin from  500 to 750 GeV.
    double SM_Charge_Asymmetry_bin_tt_750_1000;  ///< The SM value for the charged asymmetry of tt differential in the bin from  750 to 1000 GeV.
    double SM_Charge_Asymmetry_bin_tt_1000_1500; ///< The SM value for the charged asymmetry of tt differential in the bin from  1000 to 1500 GeV.
    double SM_Charge_Asymmetry_bin_tt_1500_2000; ///< The SM value for the charged asymmetry of tt differential in the bin from  1500 to 2000 GeV.
    double SM_Charge_Asymmetry_bin_tt_2000_2500; ///< The SM value for the charged asymmetry of tt differential in the bin from  2000 to 2500 GeV.
    double SM_Charge_Asymmetry_bin_tt_2500_3000; ///< The SM value for the charged asymmetry of tt differential in the bin from  2500 to 3000 GeV.
    
    double SM_ttll_bin_100_120; ///< The SM value for the differential cross section tt in the bin from 100 to 120 GeV.
    double SM_ttll_bin_120_140; ///< The SM value for the differential cross section tt in the bin from 120 to 140 GeV.
    double SM_ttll_bin_140_180; ///< The SM value for the differential cross section tt in the bin from 140 to 180 GeV.
    double SM_ttll_bin_180_500; ///< The SM value for the differential cross section tt in the bin from 180 to 500 GeV.
            
            

    double SM_tZQ_inc; ///< The SM value for the inclusive cross section tZQ.
    double SM_ttA_inc; ///< The SM value for the inclusive cross section ttA.
    double SM_ttZ_inc; ///< The SM value for the inclusive cross section ttZ.
    double SM_ttH_inc; ///< The SM value for the inclusive cross section ttH.
    double SM_ttW_inc; ///< The SM value for the inclusive cross section ttW.
    double SM_Asymmetry_leptonic_charge_rapidity_ttW; ///< The SM value for the inclusive cross section SM_Asymmetry_leptonic_charge_rapidity_ttW.
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
    */
    
    virtual bool setFlag(const std::string, const bool);

        
////////////////////////////////////////////////////////////////////////
private:

    bool flag_LHC_WG_Basis;
    bool flag_Quadratic;
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






class FB_asymmetry_Tevatron_tt_diff_mtt_top_basis_LO : public ThObservable {
public:   

    FB_asymmetry_Tevatron_tt_diff_mtt_top_basis_LO(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    
    double b_min = 0.;
    double b_max = 0.;
        
        
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


class Asymmetry_leptonic_charge_rapidity_ttW : public ThObservable {
public:   

    Asymmetry_leptonic_charge_rapidity_ttW(const StandardModel& SM_i);
    
   
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




// ttH differential cross section for different bins

/**
    * @class sigma_ttH_diff
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_ttH_diff_NLO_ATLAS_220700092: public ThObservable {
    public:

    /**
     * @brief sigma_ttH_diff_NLO constructor.
     */
    sigma_ttH_diff_NLO_ATLAS_220700092(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
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

/**
    * @class sigma_ttz_diff
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_ttz_diff_NLO_ATLAS_210312603: public ThObservable {
    public:

    /**
     * @brief sigma_ttz_diff_LO constructor.
     */
    sigma_ttz_diff_NLO_ATLAS_210312603(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
        
        
    };


    
    /**
    * @class sigma_ttz_diff
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_ttz_diff_NLO_ATLAS_231204450: public ThObservable {
    public:

    /**
     * @brief sigma_ttz_diff_NLO constructor.
     */
    sigma_ttz_diff_NLO_ATLAS_231204450(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
        
        
    };




// ttA differential cross section for different bins

    
/**
    * @class sigma_tta_diff
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tta_diff_NLO_ATLAS_emu_200706946: public ThObservable {
    public:

    /**
     * @brief sigma_tta_diff_NLO_ATLAS_emu_200706946 constructor.
     */
    sigma_tta_diff_NLO_ATLAS_emu_200706946(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    /**
    * @class sigma_tta_diff
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tta_diff_NLO_CMS_dilepton_220107301: public ThObservable {
    public:

    /**
     * @brief sigma_tta_diff_NLO_CMS_dilepton_220107301 constructor.
     */
    sigma_tta_diff_NLO_CMS_dilepton_220107301(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    

// tt differential cross section for different bins


    
    /**
    * @class sigma_ttbar_diff
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tt_diff_NLO: public ThObservable {
    public:

    /**
     * @brief sigma_tt_diff_NLO constructor.
     */
    sigma_tt_diff_NLO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    




// Charge Asymmetry ttbar for different bins

/**
    * @class charge_asymmetry_tt_diff_mtt_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class charge_asymmetry_tt_diff_mtt_NLO: public ThObservable {
    public:

    /**
     * @brief charge_asymmetry_tt_diff_mtt_LO constructor.
     */
    charge_asymmetry_tt_diff_mtt_NLO(const StandardModel& SM_i);

    /**
     * @return The value of charge_asymmetry_tt_diff_mtt_LO
     */
    double computeThValue();
    
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
// ttll

    /**
    * @class sigma_ttll_diff
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_ttll_diff_LO: public ThObservable {
    public:

    /**
     * @brief sigma_tt_diff_NLO constructor.
     */
    sigma_ttll_diff_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
//Prospects of muon collider
  
    
class sigma_mumu_VBF_3TeV_tt : public ThObservable {
public:   

    sigma_mumu_VBF_3TeV_tt(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};




class sigma_mumu_3TeV_ttH : public ThObservable {
public:   

    sigma_mumu_3TeV_ttH(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigma_mumu_3TeV_bb : public ThObservable {
public:   

    sigma_mumu_3TeV_bb(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



class sigma_mumu_VBF_10TeV_tt : public ThObservable {
public:   

    sigma_mumu_VBF_10TeV_tt(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};




class sigma_mumu_10TeV_ttH : public ThObservable {
public:   

    sigma_mumu_10TeV_ttH(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigma_mumu_10TeV_bb : public ThObservable {
public:   

    sigma_mumu_10TeV_bb(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



class sigma_mumu_VBF_30TeV_tt : public ThObservable {
public:   

    sigma_mumu_VBF_30TeV_tt(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};


class sigma_mumu_30TeV_ttH : public ThObservable {
public:   

    sigma_mumu_30TeV_ttH(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



class sigma_mumu_30TeV_bb : public ThObservable {
public:   

    sigma_mumu_30TeV_bb(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



//Prospects of Linear Collider at 250 GeV
//250 bb observables


class sigma_250_bb_eP_P30_eM_M80 : public ThObservable {
public:   

    sigma_250_bb_eP_P30_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



        
class sigma_250_bb_eP_M30_eM_P80 : public ThObservable {
public:   

    sigma_250_bb_eP_M30_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};
        
      



class a_250_bb_eP_P30_eM_M80 : public ThObservable {
public:   

    a_250_bb_eP_P30_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



        
class a_250_bb_eP_M30_eM_P80 : public ThObservable {
public:   

    a_250_bb_eP_M30_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};





//500 bb observables
     

class sigma_500_bb_eP_P30_eM_M80 : public ThObservable {
public:   

    sigma_500_bb_eP_P30_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 





class sigma_500_bb_eP_M30_eM_P80 : public ThObservable {
public:   

    sigma_500_bb_eP_M30_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 
        
        
        





class a_500_bb_eP_P30_eM_M80 : public ThObservable {
public:   

    a_500_bb_eP_P30_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 





class a_500_bb_eP_M30_eM_P80 : public ThObservable {
public:   

    a_500_bb_eP_M30_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};










class sigma_500_ttH_eP_P30_eM_M80 : public ThObservable {
public:   

    sigma_500_ttH_eP_P30_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 





class sigma_500_ttH_eP_M30_eM_P80 : public ThObservable {
public:   

    sigma_500_ttH_eP_M30_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};








//1000 bb observables
     

class sigma_1000_bb_eP_P30_eM_M80 : public ThObservable {
public:   

    sigma_1000_bb_eP_P30_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 





class sigma_1000_bb_eP_M30_eM_P80 : public ThObservable {
public:   

    sigma_1000_bb_eP_M30_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 








class a_1000_bb_eP_P30_eM_M80 : public ThObservable {
public:   

    a_1000_bb_eP_P30_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 





class a_1000_bb_eP_M30_eM_P80 : public ThObservable {
public:   

    a_1000_bb_eP_M30_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 







class sigma_1000_ttH_eP_P30_eM_M80 : public ThObservable {
public:   

    sigma_1000_ttH_eP_P30_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 





class sigma_1000_ttH_eP_M30_eM_P80 : public ThObservable {
public:   

    sigma_1000_ttH_eP_M30_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 


//380 bb observables


class sigma_380_bb_eP_0_eM_P80 : public ThObservable {
public:   

    sigma_380_bb_eP_0_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 




class sigma_380_bb_eP_0_eM_M80 : public ThObservable {
public:   

    sigma_380_bb_eP_0_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 





class a_380_bb_eP_0_eM_P80 : public ThObservable {
public:   

    a_380_bb_eP_0_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 




class a_380_bb_eP_0_eM_M80 : public ThObservable {
public:   

    a_380_bb_eP_0_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 




//1400 bb observables

class sigma_1400_bb_eP_0_eM_P80 : public ThObservable {
public:   

    sigma_1400_bb_eP_0_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 




class sigma_1400_bb_eP_0_eM_M80 : public ThObservable {
public:   

    sigma_1400_bb_eP_0_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 





class a_1400_bb_eP_0_eM_P80 : public ThObservable {
public:   

    a_1400_bb_eP_0_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 




class a_1400_bb_eP_0_eM_M80 : public ThObservable {
public:   

    a_1400_bb_eP_0_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 







class sigma_1500_ttH_eP_0_eM_P80 : public ThObservable {
public:   

    sigma_1500_ttH_eP_0_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 




class sigma_1500_ttH_eP_0_eM_M80 : public ThObservable {
public:   

    sigma_1500_ttH_eP_0_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 





//3000 bb observables

class sigma_3000_bb_eP_0_eM_P80 : public ThObservable {
public:   

    sigma_3000_bb_eP_0_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 




class sigma_3000_bb_eP_0_eM_M80 : public ThObservable {
public:   

    sigma_3000_bb_eP_0_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 





class a_3000_bb_eP_0_eM_P80 : public ThObservable {
public:   

    a_3000_bb_eP_0_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 




class a_3000_bb_eP_0_eM_M80 : public ThObservable {
public:   

    a_3000_bb_eP_0_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 





class sigma_3000_ttH_eP_0_eM_P80 : public ThObservable {
public:   

    sigma_3000_ttH_eP_0_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 




class sigma_3000_ttH_eP_0_eM_M80 : public ThObservable {
public:   

    sigma_3000_ttH_eP_0_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 


/////// Prospects of Linear Colders at 240 GeV //////////////////////////////////////////////////////////


class sigma_240_bb : public ThObservable {
public:   

    sigma_240_bb(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 



class a_240_bb : public ThObservable {
public:   

    a_240_bb(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 



/////// Prospects of Linear Colders at 360 GeV //////////////////////////////////////////////////////////


class sigma_360_bb : public ThObservable {
public:   

    sigma_360_bb(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



class a_360_bb : public ThObservable {
public:   

    a_360_bb(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};




/////// Prospects of Linear Colders at Z-Pole  //////////////////////////////////////////////////////////


class sigma_Z_pole_bb : public ThObservable {
public:   

    sigma_Z_pole_bb(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 




class sigma_Z_pole_bb_eP_0_eM_P80 : public ThObservable {
public:   

    sigma_Z_pole_bb_eP_0_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 



class sigma_Z_pole_bb_eP_0_eM_M80 : public ThObservable {
public:   

    sigma_Z_pole_bb_eP_0_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 





class sigma_Z_pole_bb_eP_M30_eM_P80 : public ThObservable {
public:   

    sigma_Z_pole_bb_eP_M30_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 








class sigma_Z_pole_bb_eP_P30_eM_M80 : public ThObservable {
public:   

    sigma_Z_pole_bb_eP_P30_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 




//////// asymmetry //////////////



class a_Z_pole_bb : public ThObservable {
public:   

    a_Z_pole_bb(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 




class a_Z_pole_bb_eP_0_eM_P80 : public ThObservable {
public:   

    a_Z_pole_bb_eP_0_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 



class a_Z_pole_bb_eP_0_eM_M80 : public ThObservable {
public:   

    a_Z_pole_bb_eP_0_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 





class a_Z_pole_bb_eP_M30_eM_P80 : public ThObservable {
public:   

    a_Z_pole_bb_eP_M30_eM_P80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 








class a_Z_pole_bb_eP_P30_eM_M80 : public ThObservable {
public:   

    a_Z_pole_bb_eP_P30_eM_M80(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
}; 






    /**
    * @class opt_obs_ilc_500_M30_P80
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for ILC
     * with polarisation for the (e+, e-)
     * -30% and 80% at 500GeV and assuming
     * a luminosity of 4iab and an 
     * efficiency*acceptance of 10%
    */
    class opt_obs_ilc_500_M30_P80: public ThObservable {
    public:

    /**
     * @brief opt_obs_ilc_500_M30_P80 constructor.
     */
    opt_obs_ilc_500_M30_P80(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_ilc_500_M30_P80
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    /**
    * @class opt_obs_ilc_500_P30_M80
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for ILC
     * with polarisation for the (e+, e-)
     * 30% and -80% at 500GeV and assuming
     * a luminosity of 4iab and an 
     * efficiency*acceptance of 10%
    */
    class opt_obs_ilc_500_P30_M80: public ThObservable {
    public:

    /**
     * @brief opt_obs_ilc_500_P30_M80 constructor.
     */
    opt_obs_ilc_500_P30_M80(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_ilc_500_P30_M80
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    /**
    * @class opt_obs_ilc_1000_M30_P80
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for ILC
     * with polarisation for the (e+, e-)
     * -30% and 80% at 500GeV and assuming
     * a luminosity of 8iab and an 
     * efficiency*acceptance of 6%
    */
    class opt_obs_ilc_1000_M30_P80: public ThObservable {
    public:

    /**
     * @brief opt_obs_ilc_1000_M30_P80 constructor.
     */
    opt_obs_ilc_1000_M30_P80(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_ilc_1000_M30_P80
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    /**
    * @class opt_obs_ilc_1000_P30_M80
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for ILC
     * with polarisation for the (e+, e-)
     * 30% and -80% at 1000GeV and assuming
     * a luminosity of 8iab and an 
     * efficiency*acceptance of 6%
    */
    class opt_obs_ilc_1000_P30_M80: public ThObservable {
    public:

    /**
     * @brief opt_obs_ilc_1000_P30_M80 constructor.
     */
    opt_obs_ilc_1000_P30_M80(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_ilc_1000_P30_M80
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    /**
    * @class opt_obs_clic_380_0_M80
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for CLIC
     * with polarisation for the (e+, e-)
     * 0% and -80% at 380GeV and assuming
     * a luminosity of 1iab and an 
     * efficiency*acceptance of 10%
    */
    class opt_obs_clic_380_0_M80: public ThObservable {
    public:

    /**
     * @brief opt_obs_clic_380_0_M80 constructor.
     */
    opt_obs_clic_380_0_M80(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_clic_380_0_M80
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    /**
    * @class opt_obs_clic_380_0_P80
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for CLIC
     * with polarisation for the (e+, e-)
     * 0% and +80% at 380GeV and assuming
     * a luminosity of 1iab and an 
     * efficiency*acceptance of 10%
    */
    class opt_obs_clic_380_0_P80: public ThObservable {
    public:

    /**
     * @brief opt_obs_clic_380_0_P80 constructor.
     */
    opt_obs_clic_380_0_P80(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_clic_380_0_P80
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    
    /**
    * @class opt_obs_clic_1500_0_M80
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for CLIC
     * with polarisation for the (e+, e-)
     * 0% and -80% at 1500GeV and assuming
     * a luminosity of 2.5iab and an 
     * efficiency*acceptance of 6%
    */
    class opt_obs_clic_1500_0_M80: public ThObservable {
    public:

    /**
     * @brief opt_obs_clic_1500_0_M80 constructor.
     */
    opt_obs_clic_1500_0_M80(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_clic_1500_0_M80
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    /**
    * @class opt_obs_clic_1500_0_P80
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for CLIC
     * with polarisation for the (e+, e-)
     * 0% and +80% at 1500GeV and assuming
     * a luminosity of 2.5iab and an 
     * efficiency*acceptance of 6%
    */
    class opt_obs_clic_1500_0_P80: public ThObservable {
    public:

    /**
     * @brief opt_obs_clic_1500_0_P80 constructor.
     */
    opt_obs_clic_1500_0_P80(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_clic_1500_0_P80
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    
    
    
    /**
    * @class opt_obs_clic_3000_0_M80
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for CLIC
     * with polarisation for the (e+, e-)
     * 0% and -80% at 3000GeV and assuming
     * a luminosity of 5iab and an 
     * efficiency*acceptance of 5%
    */
    class opt_obs_clic_3000_0_M80: public ThObservable {
    public:

    /**
     * @brief opt_obs_clic_3000_0_M80 constructor.
     */
    opt_obs_clic_3000_0_M80(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_clic_3000_0_M80
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    /**
    * @class opt_obs_clic_3000_0_P80
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for CLIC
     * with polarisation for the (e+, e-)
     * 0% and +80% at 3000GeV and assuming
     * a luminosity of 5iab and an 
     * efficiency*acceptance of 5%
    */
    class opt_obs_clic_3000_0_P80: public ThObservable {
    public:

    /**
     * @brief opt_obs_clic_3000_0_P80 constructor.
     */
    opt_obs_clic_3000_0_P80(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_clic_3000_0_P80
     */
    double computeThValue();
    

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    
    
    /**
    * @class opt_obs_fcc_350
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for FCCee
     * at 350GeV and assuming
     * a luminosity of 0.2iab and an 
     * efficiency*acceptance of 10%
    */
    class opt_obs_fcc_350: public ThObservable {
    public:

    /**
     * @brief opt_obs_fcc_350 constructor.
     */
    opt_obs_fcc_350(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_fcc_350
     */
    double computeThValue();
    
    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };


    
    /**
    * @class opt_obs_fcc_365
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for FCCee
     * at 350GeV and assuming
     * a luminosity of 1.5iab and an 
     * efficiency*acceptance of 10%
    */
    class opt_obs_fcc_365: public ThObservable {
    public:

    /**
     * @brief opt_obs_fcc_365 constructor.
     */
    opt_obs_fcc_365(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_fcc_365
     */
    double computeThValue();
    
    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    /**
    * @class opt_obs_cepc_350
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for CEPC
     * at 350GeV and assuming
     * a luminosity of 0.2iab and an 
     * efficiency*acceptance of 10%
    */
    class opt_obs_cepc_350: public ThObservable {
    public:

    /**
     * @brief opt_obs_cepc_350 constructor.
     */
    opt_obs_cepc_350(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_cepc_350
     */
    double computeThValue();
    
    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    /**
    * @class opt_obs_cepc_360
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for CEPC
     * at 360GeV and assuming
     * a luminosity of 1iab and an 
     * efficiency*acceptance of 10%
    */
    class opt_obs_cepc_360: public ThObservable {
    public:

    /**
     * @brief opt_obs_cepc_360 constructor.
     */
    opt_obs_cepc_360(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_cepc_360
     */
    double computeThValue();
    
    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
        
    /**
    * @class opt_obs_muon_3TeV
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for a
     * Muon Collider at 3TeV and assuming
     * a luminosity of 1iab and an 
     * efficiency*acceptance of 5%
    */
    class opt_obs_muon_3TeV: public ThObservable {
    public:

    /**
     * @brief opt_obs_muon_3TeV constructor.
     */
    opt_obs_muon_3TeV(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_muon_3TeV
     */
    double computeThValue();
    
    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    /**
    * @class opt_obs_muon_10TeV
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for a
     * Muon Collider at 10TeV and assuming
     * a luminosity of 10iab and an 
     * efficiency*acceptance of 2.5%
    */
    class opt_obs_muon_10TeV: public ThObservable {
    public:

    /**
     * @brief opt_obs_muon_10TeV constructor.
     */
    opt_obs_muon_10TeV(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_muon_10TeV
     */
    double computeThValue();
    
    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    /**
    * @class opt_obs_muon_30TeV
    * @ingroup TopQuarkObservables 
    * @brief Optimal observables for a
     * Muon Collider at 30TeV and assuming
     * a luminosity of 90iab and an 
     * efficiency*acceptance of 1.%
    */
    class opt_obs_muon_30TeV: public ThObservable {
    public:

    /**
     * @brief opt_obs_muon_30TeV constructor.
     */
    opt_obs_muon_30TeV(const StandardModel& SM_i);

    /**
     * @return The value of opt_obs_muon_30TeV
     */
    double computeThValue();
    
    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    

    
    
    /*
    * @class b4_ttH_LO
    * @ingroup TopQuarkObservables 
    * @brief b4 observable for ttH
    */
    class b4_ttH_LO: public ThObservable {
    public:

    /**
     * @brief b4_ttH_LO constructor.
     */
    b4_ttH_LO(const StandardModel& SM_i);

    /**
     * @return The value of b4_ttH_LO
     */
    double computeThValue();
    
    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    /*
    * @class Asymmetry_Dazi_ord_ttH
    * @ingroup TopQuarkObservables 
    * @brief Asymmetry of the variation
     * of the azimuthal angle of the 
     * top pairs in the tth process
     * ordered by the pseudorapidty
    */
    class Asymmetry_Dazi_ord_ttH: public ThObservable {
    public:

    /**
     * @brief Asymmetry_Dazi_ord_ttH constructor.
     */
    Asymmetry_Dazi_ord_ttH(const StandardModel& SM_i);

    /**
     * @return The value of Asymmetry_Dazi_ord_ttH
     */
    double computeThValue();
    
    private:

        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    /*
    * @class Asymmetry_Dazi_ord_ttH
    * @ingroup TopQuarkObservables 
    * @brief Asymmetry of the variation
     * of the azimuthal angle of the 
     * final electron pairs in the 
     * tth process ordered by the 
     * pseudorapidty
    */
    class Asymmetry_Dazi_ord_ttH_ee: public ThObservable {
    public:

    /**
     * @brief Asymmetry_Dazi_ord_ttH_ee constructor.
     */
    Asymmetry_Dazi_ord_ttH_ee(const StandardModel& SM_i);

    /**
     * @return The value of Asymmetry_Dazi_ord_ttH_ee
     */
    double computeThValue();
    
    private:

        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    /*
    * @class sigma_ttH_diff_LO_mtth
    * @ingroup TopQuarkObservables 
    * @brief ttH differential cross
     * section w.r.t. the invariant
     * mass of tth
    */
    class sigma_ttH_diff_LO_mtth: public ThObservable {
    public:

    /**
     * @brief sigma_ttH_diff_LO_mtth constructor.
     */
    sigma_ttH_diff_LO_mtth(const StandardModel& SM_i);

    /**
     * @return The value of sigma_ttH_diff_LO_mtth
     */
    double computeThValue();
    
    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    /*
    * @class Asymmetry_trip_prod_pt_pe_pp_ttH
    * @ingroup TopQuarkObservables 
    * @brief Asymmetry of the triple product
     * p_t\cdot(p_{e^-}\times p_{e^+})
     * in the tth process 
    */
    class Asymmetry_trip_prod_pt_pe_pp_ttH: public ThObservable {
    public:

    /**
     * @brief Asymmetry_trip_prod_pt_pe_pp_ttH constructor.
     */
    Asymmetry_trip_prod_pt_pe_pp_ttH(const StandardModel& SM_i);

    /**
     * @return The value of Asymmetry_trip_prod_pt_pe_pp_ttH
     */
    double computeThValue();
    
    private:
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    /*
    * @class Asymmetry_sign_trip_prod_pe_pp_ttH
    * @ingroup TopQuarkObservables 
    * @brief Asymmetry of the scalar product
     * p_{e^-}\cdot p_{e^+} in the tth process 
     * scaled by the sign of the triple product
     * p_t\cdot(p_{e^-}\times p_{e^+})
    */
    class Asymmetry_sign_trip_prod_pe_pp_ttH: public ThObservable {
    public:

    /**
     * @brief Asymmetry_sign_trip_prod_pe_pp_ttH constructor.
     */
    Asymmetry_sign_trip_prod_pe_pp_ttH(const StandardModel& SM_i);

    /**
     * @return The value of Asymmetry_sign_trip_prod_pe_pp_ttH
     */
    double computeThValue();
    
    private:
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    
    
    
    
    
    /*
    * @class Asymmetry_cos_je_tHj
    * @ingroup TopQuarkObservables 
    * @brief Asymmetry of the cosine of the 
     * polarisation angle defined with respect
     * to the axis "je" in the thj process.
    */
    class Asymmetry_cos_je_tHj: public ThObservable {
    public:

    /**
     * @brief Asymmetry_cos_je_tHj constructor.
     */
    Asymmetry_cos_je_tHj(const StandardModel& SM_i);

    /**
     * @return The value of Asymmetry_cos_je_tHj
     */
    double computeThValue();
    
    private:
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    /*
    * @class Asymmetry_cos_se_tHj
    * @ingroup TopQuarkObservables 
    * @brief Asymmetry of the cosine of the 
     * polarisation angle defined with respect
     * to the axis "se" in the thj process.
    */
    class Asymmetry_cos_se_tHj: public ThObservable {
    public:

    /**
     * @brief Asymmetry_cos_se_tHj constructor.
     */
    Asymmetry_cos_se_tHj(const StandardModel& SM_i);

    /**
     * @return The value of Asymmetry_cos_se_tHj
     */
    double computeThValue();
    
    private:
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    /*
    * @class Asymmetry_cos_ye_tHj
    * @ingroup TopQuarkObservables 
    * @brief Asymmetry of the cosine of the 
     * polarisation angle defined with respect
     * to the axis "ye" in the thj process.
    */
    class Asymmetry_cos_ye_tHj: public ThObservable {
    public:

    /**
     * @brief Asymmetry_cos_ye_tHj constructor.
     */
    Asymmetry_cos_ye_tHj(const StandardModel& SM_i);

    /**
     * @return The value of Asymmetry_cos_ye_tHj
     */
    double computeThValue();
    
    private:
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    
    
    /*
    * @class Asymmetry_trip_prod_ph_pt_pj_tHj
    * @ingroup TopQuarkObservables 
    * @brief Asymmetry of the triple product
     * ph\cross(pt\cdot pj) in the thj process.
    */
    class Asymmetry_trip_prod_ph_pt_pj_tHj: public ThObservable {
    public:

    /**
     * @brief Asymmetry_trip_prod_ph_pt_pj_tHj constructor.
     */
    Asymmetry_trip_prod_ph_pt_pj_tHj(const StandardModel& SM_i);

    /**
     * @return The value of Asymmetry_trip_prod_ph_pt_pj_tHj
     */
    double computeThValue();
    
    private:
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    
    
    
    /*
    * @class sigma_tHj_diff_LO_Del_R_th
    * @ingroup TopQuarkObservables 
    * @brief Differential cross section
     * with respect to the \Delta R(th)
     *  in the thj process.
    */
    class sigma_tHj_diff_LO_Del_R_th: public ThObservable {
    public:

    /**
     * @brief sigma_tHj_diff_LO_Del_R_th constructor.
     */
    sigma_tHj_diff_LO_Del_R_th(const StandardModel& SM_i);

    /**
     * @return The value of sigma_tHj_diff_LO_Del_R_th
     */
    double computeThValue();
    
    private:
        
        double b_min = 0.;
        double b_max = 0.;
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    
    
    
    /*
    * @class sigma_tHj_diff_LO_mth
    * @ingroup TopQuarkObservables 
    * @brief Differential cross section
     * with respect to the invariant
     * mass of the Higgs and top
     *  in the thj process.
    */
    class sigma_tHj_diff_LO_mth: public ThObservable {
    public:

    /**
     * @brief sigma_tHj_diff_LO_mth constructor.
     */
    sigma_tHj_diff_LO_mth(const StandardModel& SM_i);

    /**
     * @return The value of sigma_tHj_diff_LO_mth
     */
    double computeThValue();
    
    private:
        
        double b_min = 0.;
        double b_max = 0.;
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
    

    
    
    /*
    * @class sigma_tHj_diff_LO_trip_prod_z_pt_pj
    * @ingroup TopQuarkObservables 
    * @brief Differential cross section with
     * respect to the triple product of the 
     * 3-momenta of the proton beam the 
     * top-quark and the jet in the thj process.
    */
    class sigma_tHj_diff_LO_trip_prod_z_pt_pj: public ThObservable {
    public:

    /**
     * @brief sigma_tHj_diff_LO_trip_prod_z_pt_pj constructor.
     */
    sigma_tHj_diff_LO_trip_prod_z_pt_pj(const StandardModel& SM_i);

    /**
     * @return The value of sigma_tHj_diff_LO_trip_prod_z_pt_pj
     */
    double computeThValue();
    
    private:
        
        double b_min = 0.;
        double b_max = 0.;
        
        const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;

    };
        
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// OLD //////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////        
//////////////////////////////////////////////////////////////////////////////////////////////////////        
        
//250 bb observables OLD



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



class op_eigen_ttll_1 : public ThObservable {
public:   

    op_eigen_ttll_1(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



class op_eigen_ttll_2 : public ThObservable {
public:   

    op_eigen_ttll_2(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



class op_eigen_ttll_3 : public ThObservable {
public:   

    op_eigen_ttll_3(const StandardModel& SM_i);
    
   
    double computeThValue();
    
private:
    const NPSMEFT6dtopquark& myNPSMEFT6dtopquark;
};



class op_eigen_ttll_4 : public ThObservable {
public:   

    op_eigen_ttll_4(const StandardModel& SM_i);
    
   
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

/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   TopQuarkObservables.h
 * Author: silvest
 *
 * Created on 21 settembre 2023, 16.45
 */

#ifndef TOPQUARKOBSERVABLES_H
#define TOPQUARKOBSERVABLES_H

#include "NPSMEFTd6General.h"

class TopQuarkObservables {
public:
    
    /**
    *　@brief The number of the model parameters in %TopQuarkObservablesVars. 
    */
    static const int NTopQuarkObservablesVars = 104;
   
    /*
    @brief A string array containing the labels of the model parameters in  TopQuarkObservables
    */
    static const std::string TopQuarkObservablesVars[NTopQuarkObservablesVars];
    
    
    
    TopQuarkObservables(const NPSMEFTd6General& NP_in);
    TopQuarkObservables(const TopQuarkObservables& orig);
    virtual ~TopQuarkObservables();
    
    double GetAFBLR_SM() const
    {
        return AFBLR_SM;
    }

    void SetAFBLR_SM(double AFBLR_SM)
    {
        this->AFBLR_SM = AFBLR_SM;
    }

    double GetF0_SM() const
    {
        return F0_SM;
    }

    void SetF0_SM(double F0_SM)
    {
        this->F0_SM = F0_SM;
    }

    double GetFL_SM() const
    {
        return FL_SM;
    }

    void SetFL_SM(double FL_SM)
    {
        this->FL_SM = FL_SM;
    }

    double GetRb_SM() const
    {
        return Rb_SM;
    }

    void SetRb_SM(double Rb_SM)
    {
        this->Rb_SM = Rb_SM;
    }

    double GetSM_Asymmetry_leptonic_charge_rapidity_ttW() const
    {
        return SM_Asymmetry_leptonic_charge_rapidity_ttW;
    }

    void SetSM_Asymmetry_leptonic_charge_rapidity_ttW(double SM_Asymmetry_leptonic_charge_rapidity_ttW)
    {
        this->SM_Asymmetry_leptonic_charge_rapidity_ttW = SM_Asymmetry_leptonic_charge_rapidity_ttW;
    }

    double GetSM_Charge_Asymmetry_bin_tt_0_500() const
    {
        return SM_Charge_Asymmetry_bin_tt_0_500;
    }

    void SetSM_Charge_Asymmetry_bin_tt_0_500(double SM_Charge_Asymmetry_bin_tt_0_500)
    {
        this->SM_Charge_Asymmetry_bin_tt_0_500 = SM_Charge_Asymmetry_bin_tt_0_500;
    }

    double GetSM_Charge_Asymmetry_bin_tt_1000_1500() const
    {
        return SM_Charge_Asymmetry_bin_tt_1000_1500;
    }

    void SetSM_Charge_Asymmetry_bin_tt_1000_1500(double SM_Charge_Asymmetry_bin_tt_1000_1500)
    {
        this->SM_Charge_Asymmetry_bin_tt_1000_1500 = SM_Charge_Asymmetry_bin_tt_1000_1500;
    }

    double GetSM_Charge_Asymmetry_bin_tt_1500_2000() const
    {
        return SM_Charge_Asymmetry_bin_tt_1500_2000;
    }

    void SetSM_Charge_Asymmetry_bin_tt_1500_2000(double SM_Charge_Asymmetry_bin_tt_1500_2000)
    {
        this->SM_Charge_Asymmetry_bin_tt_1500_2000 = SM_Charge_Asymmetry_bin_tt_1500_2000;
    }

    double GetSM_Charge_Asymmetry_bin_tt_2000_2500() const
    {
        return SM_Charge_Asymmetry_bin_tt_2000_2500;
    }

    void SetSM_Charge_Asymmetry_bin_tt_2000_2500(double SM_Charge_Asymmetry_bin_tt_2000_2500)
    {
        this->SM_Charge_Asymmetry_bin_tt_2000_2500 = SM_Charge_Asymmetry_bin_tt_2000_2500;
    }

    double GetSM_Charge_Asymmetry_bin_tt_2500_3000() const
    {
        return SM_Charge_Asymmetry_bin_tt_2500_3000;
    }

    void SetSM_Charge_Asymmetry_bin_tt_2500_3000(double SM_Charge_Asymmetry_bin_tt_2500_3000)
    {
        this->SM_Charge_Asymmetry_bin_tt_2500_3000 = SM_Charge_Asymmetry_bin_tt_2500_3000;
    }

    double GetSM_Charge_Asymmetry_bin_tt_500_750() const
    {
        return SM_Charge_Asymmetry_bin_tt_500_750;
    }

    void SetSM_Charge_Asymmetry_bin_tt_500_750(double SM_Charge_Asymmetry_bin_tt_500_750)
    {
        this->SM_Charge_Asymmetry_bin_tt_500_750 = SM_Charge_Asymmetry_bin_tt_500_750;
    }

    double GetSM_Charge_Asymmetry_bin_tt_750_1000() const
    {
        return SM_Charge_Asymmetry_bin_tt_750_1000;
    }

    void SetSM_Charge_Asymmetry_bin_tt_750_1000(double SM_Charge_Asymmetry_bin_tt_750_1000)
    {
        this->SM_Charge_Asymmetry_bin_tt_750_1000 = SM_Charge_Asymmetry_bin_tt_750_1000;
    }

    double GetSM_sigmaschannel8() const
    {
        return SM_sigmaschannel8;
    }

    void SetSM_sigmaschannel8(double SM_sigmaschannel8)
    {
        this->SM_sigmaschannel8 = SM_sigmaschannel8;
    }

    double GetSM_sigmaschannelTev() const
    {
        return SM_sigmaschannelTev;
    }

    void SetSM_sigmaschannelTev(double SM_sigmaschannelTev)
    {
        this->SM_sigmaschannelTev = SM_sigmaschannelTev;
    }

    double GetSM_sigmatchannel13() const
    {
        return SM_sigmatchannel13;
    }

    void SetSM_sigmatchannel13(double SM_sigmatchannel13)
    {
        this->SM_sigmatchannel13 = SM_sigmatchannel13;
    }

    double GetSM_sigmatchannel8() const
    {
        return SM_sigmatchannel8;
    }

    void SetSM_sigmatchannel8(double SM_sigmatchannel8)
    {
        this->SM_sigmatchannel8 = SM_sigmatchannel8;
    }

    double GetSM_tAq_inc() const
    {
        return SM_tAq_inc;
    }

    void SetSM_tAq_inc(double SM_tAq_inc)
    {
        this->SM_tAq_inc = SM_tAq_inc;
    }

    double GetSM_tH_tchan_value() const
    {
        return SM_tH_tchan_value;
    }

    void SetSM_tH_tchan_value(double SM_tH_tchan_value)
    {
        this->SM_tH_tchan_value = SM_tH_tchan_value;
    }

    double GetSM_tW_inc() const
    {
        return SM_tW_inc;
    }

    void SetSM_tW_inc(double SM_tW_inc)
    {
        this->SM_tW_inc = SM_tW_inc;
    }

    double GetSM_tW_inc_8TeV() const
    {
        return SM_tW_inc_8TeV;
    }

    void SetSM_tW_inc_8TeV(double SM_tW_inc_8TeV)
    {
        this->SM_tW_inc_8TeV = SM_tW_inc_8TeV;
    }

    double GetSM_tZQ_inc() const
    {
        return SM_tZQ_inc;
    }

    void SetSM_tZQ_inc(double SM_tZQ_inc)
    {
        this->SM_tZQ_inc = SM_tZQ_inc;
    }

    double GetSM_ttA_bin_132_180() const
    {
        return SM_ttA_bin_132_180;
    }

    void SetSM_ttA_bin_132_180(double SM_ttA_bin_132_180)
    {
        this->SM_ttA_bin_132_180 = SM_ttA_bin_132_180;
    }

    double GetSM_ttA_bin_180_300() const
    {
        return SM_ttA_bin_180_300;
    }

    void SetSM_ttA_bin_180_300(double SM_ttA_bin_180_300)
    {
        this->SM_ttA_bin_180_300 = SM_ttA_bin_180_300;
    }

    double GetSM_ttA_bin_20_25() const
    {
        return SM_ttA_bin_20_25;
    }

    void SetSM_ttA_bin_20_25(double SM_ttA_bin_20_25)
    {
        this->SM_ttA_bin_20_25 = SM_ttA_bin_20_25;
    }

    double GetSM_ttA_bin_25_30() const
    {
        return SM_ttA_bin_25_30;
    }

    void SetSM_ttA_bin_25_30(double SM_ttA_bin_25_30)
    {
        this->SM_ttA_bin_25_30 = SM_ttA_bin_25_30;
    }

    double GetSM_ttA_bin_30_35() const
    {
        return SM_ttA_bin_30_35;
    }

    void SetSM_ttA_bin_30_35(double SM_ttA_bin_30_35)
    {
        this->SM_ttA_bin_30_35 = SM_ttA_bin_30_35;
    }

    double GetSM_ttA_bin_35_40() const
    {
        return SM_ttA_bin_35_40;
    }

    void SetSM_ttA_bin_35_40(double SM_ttA_bin_35_40)
    {
        this->SM_ttA_bin_35_40 = SM_ttA_bin_35_40;
    }

    double GetSM_ttA_bin_40_47() const
    {
        return SM_ttA_bin_40_47;
    }

    void SetSM_ttA_bin_40_47(double SM_ttA_bin_40_47)
    {
        this->SM_ttA_bin_40_47 = SM_ttA_bin_40_47;
    }

    double GetSM_ttA_bin_47_55() const
    {
        return SM_ttA_bin_47_55;
    }

    void SetSM_ttA_bin_47_55(double SM_ttA_bin_47_55)
    {
        this->SM_ttA_bin_47_55 = SM_ttA_bin_47_55;
    }

    double GetSM_ttA_bin_55_70() const
    {
        return SM_ttA_bin_55_70;
    }

    void SetSM_ttA_bin_55_70(double SM_ttA_bin_55_70)
    {
        this->SM_ttA_bin_55_70 = SM_ttA_bin_55_70;
    }

    double GetSM_ttA_bin_70_85() const
    {
        return SM_ttA_bin_70_85;
    }

    void SetSM_ttA_bin_70_85(double SM_ttA_bin_70_85)
    {
        this->SM_ttA_bin_70_85 = SM_ttA_bin_70_85;
    }

    double GetSM_ttA_bin_85_132() const
    {
        return SM_ttA_bin_85_132;
    }

    void SetSM_ttA_bin_85_132(double SM_ttA_bin_85_132)
    {
        this->SM_ttA_bin_85_132 = SM_ttA_bin_85_132;
    }

    double GetSM_ttA_inc() const
    {
        return SM_ttA_inc;
    }

    void SetSM_ttA_inc(double SM_ttA_inc)
    {
        this->SM_ttA_inc = SM_ttA_inc;
    }

    double GetSM_ttH_inc() const
    {
        return SM_ttH_inc;
    }

    void SetSM_ttH_inc(double SM_ttH_inc)
    {
        this->SM_ttH_inc = SM_ttH_inc;
    }

    double GetSM_ttW_inc() const
    {
        return SM_ttW_inc;
    }

    void SetSM_ttW_inc(double SM_ttW_inc)
    {
        this->SM_ttW_inc = SM_ttW_inc;
    }

    double GetSM_ttZ_bin_0_40() const
    {
        return SM_ttZ_bin_0_40;
    }

    void SetSM_ttZ_bin_0_40(double SM_ttZ_bin_0_40)
    {
        this->SM_ttZ_bin_0_40 = SM_ttZ_bin_0_40;
    }

    double GetSM_ttZ_bin_110_160() const
    {
        return SM_ttZ_bin_110_160;
    }

    void SetSM_ttZ_bin_110_160(double SM_ttZ_bin_110_160)
    {
        this->SM_ttZ_bin_110_160 = SM_ttZ_bin_110_160;
    }

    double GetSM_ttZ_bin_160_220() const
    {
        return SM_ttZ_bin_160_220;
    }

    void SetSM_ttZ_bin_160_220(double SM_ttZ_bin_160_220)
    {
        this->SM_ttZ_bin_160_220 = SM_ttZ_bin_160_220;
    }

    double GetSM_ttZ_bin_220_290() const
    {
        return SM_ttZ_bin_220_290;
    }

    void SetSM_ttZ_bin_220_290(double SM_ttZ_bin_220_290)
    {
        this->SM_ttZ_bin_220_290 = SM_ttZ_bin_220_290;
    }

    double GetSM_ttZ_bin_290_400() const
    {
        return SM_ttZ_bin_290_400;
    }

    void SetSM_ttZ_bin_290_400(double SM_ttZ_bin_290_400)
    {
        this->SM_ttZ_bin_290_400 = SM_ttZ_bin_290_400;
    }

    double GetSM_ttZ_bin_40_70() const
    {
        return SM_ttZ_bin_40_70;
    }

    void SetSM_ttZ_bin_40_70(double SM_ttZ_bin_40_70)
    {
        this->SM_ttZ_bin_40_70 = SM_ttZ_bin_40_70;
    }

    double GetSM_ttZ_bin_70_110() const
    {
        return SM_ttZ_bin_70_110;
    }

    void SetSM_ttZ_bin_70_110(double SM_ttZ_bin_70_110)
    {
        this->SM_ttZ_bin_70_110 = SM_ttZ_bin_70_110;
    }

    double GetSM_ttZ_inc() const
    {
        return SM_ttZ_inc;
    }

    void SetSM_ttZ_inc(double SM_ttZ_inc)
    {
        this->SM_ttZ_inc = SM_ttZ_inc;
    }

    double GetSM_tt_bin_1000_1150() const
    {
        return SM_tt_bin_1000_1150;
    }

    void SetSM_tt_bin_1000_1150(double SM_tt_bin_1000_1150)
    {
        this->SM_tt_bin_1000_1150 = SM_tt_bin_1000_1150;
    }

    double GetSM_tt_bin_1150_1300() const
    {
        return SM_tt_bin_1150_1300;
    }

    void SetSM_tt_bin_1150_1300(double SM_tt_bin_1150_1300)
    {
        this->SM_tt_bin_1150_1300 = SM_tt_bin_1150_1300;
    }

    double GetSM_tt_bin_1300_1500() const
    {
        return SM_tt_bin_1300_1500;
    }

    void SetSM_tt_bin_1300_1500(double SM_tt_bin_1300_1500)
    {
        this->SM_tt_bin_1300_1500 = SM_tt_bin_1300_1500;
    }

    double GetSM_tt_bin_1500_1700() const
    {
        return SM_tt_bin_1500_1700;
    }

    void SetSM_tt_bin_1500_1700(double SM_tt_bin_1500_1700)
    {
        this->SM_tt_bin_1500_1700 = SM_tt_bin_1500_1700;
    }

    double GetSM_tt_bin_1700_2000() const
    {
        return SM_tt_bin_1700_2000;
    }

    void SetSM_tt_bin_1700_2000(double SM_tt_bin_1700_2000)
    {
        this->SM_tt_bin_1700_2000 = SM_tt_bin_1700_2000;
    }

    double GetSM_tt_bin_2000_2300() const
    {
        return SM_tt_bin_2000_2300;
    }

    void SetSM_tt_bin_2000_2300(double SM_tt_bin_2000_2300)
    {
        this->SM_tt_bin_2000_2300 = SM_tt_bin_2000_2300;
    }

    double GetSM_tt_bin_2300_2600() const
    {
        return SM_tt_bin_2300_2600;
    }

    void SetSM_tt_bin_2300_2600(double SM_tt_bin_2300_2600)
    {
        this->SM_tt_bin_2300_2600 = SM_tt_bin_2300_2600;
    }

    double GetSM_tt_bin_250_400() const
    {
        return SM_tt_bin_250_400;
    }

    void SetSM_tt_bin_250_400(double SM_tt_bin_250_400)
    {
        this->SM_tt_bin_250_400 = SM_tt_bin_250_400;
    }

    double GetSM_tt_bin_2600_3000() const
    {
        return SM_tt_bin_2600_3000;
    }

    void SetSM_tt_bin_2600_3000(double SM_tt_bin_2600_3000)
    {
        this->SM_tt_bin_2600_3000 = SM_tt_bin_2600_3000;
    }

    double GetSM_tt_bin_3000_3500() const
    {
        return SM_tt_bin_3000_3500;
    }

    void SetSM_tt_bin_3000_3500(double SM_tt_bin_3000_3500)
    {
        this->SM_tt_bin_3000_3500 = SM_tt_bin_3000_3500;
    }

    double GetSM_tt_bin_3500_4000() const
    {
        return SM_tt_bin_3500_4000;
    }

    void SetSM_tt_bin_3500_4000(double SM_tt_bin_3500_4000)
    {
        this->SM_tt_bin_3500_4000 = SM_tt_bin_3500_4000;
    }

    double GetSM_tt_bin_400_480() const
    {
        return SM_tt_bin_400_480;
    }

    void SetSM_tt_bin_400_480(double SM_tt_bin_400_480)
    {
        this->SM_tt_bin_400_480 = SM_tt_bin_400_480;
    }

    double GetSM_tt_bin_480_560() const
    {
        return SM_tt_bin_480_560;
    }

    void SetSM_tt_bin_480_560(double SM_tt_bin_480_560)
    {
        this->SM_tt_bin_480_560 = SM_tt_bin_480_560;
    }

    double GetSM_tt_bin_560_640() const
    {
        return SM_tt_bin_560_640;
    }

    void SetSM_tt_bin_560_640(double SM_tt_bin_560_640)
    {
        this->SM_tt_bin_560_640 = SM_tt_bin_560_640;
    }

    double GetSM_tt_bin_640_720() const
    {
        return SM_tt_bin_640_720;
    }

    void SetSM_tt_bin_640_720(double SM_tt_bin_640_720)
    {
        this->SM_tt_bin_640_720 = SM_tt_bin_640_720;
    }

    double GetSM_tt_bin_720_800() const
    {
        return SM_tt_bin_720_800;
    }

    void SetSM_tt_bin_720_800(double SM_tt_bin_720_800)
    {
        this->SM_tt_bin_720_800 = SM_tt_bin_720_800;
    }

    double GetSM_tt_bin_800_900() const
    {
        return SM_tt_bin_800_900;
    }

    void SetSM_tt_bin_800_900(double SM_tt_bin_800_900)
    {
        this->SM_tt_bin_800_900 = SM_tt_bin_800_900;
    }

    double GetSM_tt_bin_900_1000() const
    {
        return SM_tt_bin_900_1000;
    }

    void SetSM_tt_bin_900_1000(double SM_tt_bin_900_1000)
    {
        this->SM_tt_bin_900_1000 = SM_tt_bin_900_1000;
    }

    double GetSM_ttbar_LHC13() const
    {
        return SM_ttbar_LHC13;
    }

    void SetSM_ttbar_LHC13(double SM_ttbar_LHC13)
    {
        this->SM_ttbar_LHC13 = SM_ttbar_LHC13;
    }

    double GetSM_ttbar_LHC8() const
    {
        return SM_ttbar_LHC8;
    }

    void SetSM_ttbar_LHC8(double SM_ttbar_LHC8)
    {
        this->SM_ttbar_LHC8 = SM_ttbar_LHC8;
    }

    double GetSM_ttbar_Tev() const
    {
        return SM_ttbar_Tev;
    }

    void SetSM_ttbar_Tev(double SM_ttbar_Tev)
    {
        this->SM_ttbar_Tev = SM_ttbar_Tev;
    }

    double GetSM_ttll_bin_100_120() const
    {
        return SM_ttll_bin_100_120;
    }

    void SetSM_ttll_bin_100_120(double SM_ttll_bin_100_120)
    {
        this->SM_ttll_bin_100_120 = SM_ttll_bin_100_120;
    }

    double GetSM_ttll_bin_120_140() const
    {
        return SM_ttll_bin_120_140;
    }

    void SetSM_ttll_bin_120_140(double SM_ttll_bin_120_140)
    {
        this->SM_ttll_bin_120_140 = SM_ttll_bin_120_140;
    }

    double GetSM_ttll_bin_140_180() const
    {
        return SM_ttll_bin_140_180;
    }

    void SetSM_ttll_bin_140_180(double SM_ttll_bin_140_180)
    {
        this->SM_ttll_bin_140_180 = SM_ttll_bin_140_180;
    }

    double GetSM_ttll_bin_180_500() const
    {
        return SM_ttll_bin_180_500;
    }

    void SetSM_ttll_bin_180_500(double SM_ttll_bin_180_500)
    {
        this->SM_ttll_bin_180_500 = SM_ttll_bin_180_500;
    }

    double GetTtWqEM_SM() const
    {
        return ttWqEM_SM;
    }

    void SetTtWqEM_SM(double ttWqEM_SM)
    {
        this->ttWqEM_SM = ttWqEM_SM;
    }

protected:
        
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

private:
    
    const NPSMEFTd6General& NP;

};

#endif /* TOPQUARKOBSERVABLES_H */


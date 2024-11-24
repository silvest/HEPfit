/*
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


/* Effective dimension-6 operators describing four-fermion operators involved in
* tt production at hadron colliders. Such basis is described in 1807.02121
*/

#include "NPSMEFT6dtopquark.h"
#include <limits>
#include <gsl/gsl_sf.h>
#include "std_make_vector.h"


const std::string NPSMEFT6dtopquark::NPSMEFT6dtopquarkVars[NNPSMEFT6dtopquarkVars]
        = {"C_phit","C_phiQ3","C_phiQ1","C_phiQm","C_tW",
        "C_tZ","C_tB","C_tphi","C_phib","C_bW","C_bB","C_bZ",
        "C_phitb","C_tG","C_ed","C_eq","C_ld","C_lqP","C_eu",
        "C_lu","C_lqM","C_tu8","C_td8", "C_Qq18","C_tq8",
        "C_Qq38","C_Qu8","C_Qd8","C_Qd1","C_Qu1","C_td1",
        "C_tu1","C_tq1","C_Qq11","C_Qq31"
        //Operators for CPV study in tth and thj
        ,"C_phiG","C_phiGtil","C_phiW","C_phiWtil",
        "C_tphiIm","C_tGIm","C_tWIm","C_tBIm","C_tZIm"
        };


NPSMEFT6dtopquark::NPSMEFT6dtopquark()
: NPbase()
{
    
    
    ModelParamMap.insert(std::make_pair("C_phit", std::cref(C_phit)));
    ModelParamMap.insert(std::make_pair("C_phiQ3", std::cref(C_phiQ3)));
    ModelParamMap.insert(std::make_pair("C_phiQ1", std::cref(C_phiQ1)));
    ModelParamMap.insert(std::make_pair("C_phiQm", std::cref(C_phiQm)));
    ModelParamMap.insert(std::make_pair("C_tW", std::cref(C_tW)));
    ModelParamMap.insert(std::make_pair("C_tZ", std::cref(C_tZ)));
    ModelParamMap.insert(std::make_pair("C_tB", std::cref(C_tB)));
    ModelParamMap.insert(std::make_pair("C_tphi", std::cref(C_tphi)));
    ModelParamMap.insert(std::make_pair("C_phib", std::cref(C_phib)));
    ModelParamMap.insert(std::make_pair("C_bW", std::cref(C_bW)));
    ModelParamMap.insert(std::make_pair("C_bB", std::cref(C_bB)));
    ModelParamMap.insert(std::make_pair("C_bZ", std::cref(C_bZ)));
    ModelParamMap.insert(std::make_pair("C_phitb", std::cref(C_phitb)));
    ModelParamMap.insert(std::make_pair("C_tG", std::cref(C_tG)));
    
    ModelParamMap.insert(std::make_pair("C_ed", std::cref(C_ed)));
    ModelParamMap.insert(std::make_pair("C_eq", std::cref(C_eq)));
    ModelParamMap.insert(std::make_pair("C_ld", std::cref(C_ld)));
    ModelParamMap.insert(std::make_pair("C_lqP", std::cref(C_lqP)));
    ModelParamMap.insert(std::make_pair("C_eu", std::cref(C_eu)));
    ModelParamMap.insert(std::make_pair("C_lu", std::cref(C_lu)));
    ModelParamMap.insert(std::make_pair("C_lqM", std::cref(C_lqM)));
    
    ModelParamMap.insert(std::make_pair("C_tu8", std::cref(C_tu8)));
    ModelParamMap.insert(std::make_pair("C_td8", std::cref(C_td8)));
    ModelParamMap.insert(std::make_pair("C_Qq18", std::cref(C_Qq18)));
    ModelParamMap.insert(std::make_pair("C_tq8", std::cref(C_tq8)));
    ModelParamMap.insert(std::make_pair("C_Qq38", std::cref(C_Qq38)));
    ModelParamMap.insert(std::make_pair("C_Qu8", std::cref(C_Qu8)));
    ModelParamMap.insert(std::make_pair("C_Qd8", std::cref(C_Qd8)));
    ModelParamMap.insert(std::make_pair("C_Qd1", std::cref(C_Qd1)));
    ModelParamMap.insert(std::make_pair("C_Qu1", std::cref(C_Qu1)));
    ModelParamMap.insert(std::make_pair("C_td1", std::cref(C_td1)));
    ModelParamMap.insert(std::make_pair("C_tu1", std::cref(C_tu1)));
    ModelParamMap.insert(std::make_pair("C_tq1", std::cref(C_tq1)));
    ModelParamMap.insert(std::make_pair("C_Qq11", std::cref(C_Qq11)));
    ModelParamMap.insert(std::make_pair("C_Qq31", std::cref(C_Qq31)));
    
    ModelParamMap.insert(std::make_pair("C_phiG", std::cref(C_phiG)));
    ModelParamMap.insert(std::make_pair("C_phiGtil", std::cref(C_phiGtil)));
    ModelParamMap.insert(std::make_pair("C_phiW", std::cref(C_phiW)));
    ModelParamMap.insert(std::make_pair("C_phiWtil", std::cref(C_phiWtil)));
    ModelParamMap.insert(std::make_pair("C_tphiIm", std::cref(C_tphiIm)));
    ModelParamMap.insert(std::make_pair("C_tGIm", std::cref(C_tGIm)));
    ModelParamMap.insert(std::make_pair("C_tWIm", std::cref(C_tWIm)));
    ModelParamMap.insert(std::make_pair("C_tBIm", std::cref(C_tBIm)));
    ModelParamMap.insert(std::make_pair("C_tZIm", std::cref(C_tZIm)));

    
    flag_LHC_WG_Basis = false;
    flag_Quadratic = false;
}



void NPSMEFT6dtopquark::setParameter(const std::string name, const double& value)
{
    if (name.compare("C_phit") == 0)
        C_phit = value;
    else if (name.compare("C_phiQ3") == 0)
        C_phiQ3 = value;
    else if (name.compare("C_phiQ1") == 0)
        C_phiQ1 = value;
    else if (name.compare("C_phiQm") == 0)
        C_phiQm = value;
    else if (name.compare("C_tW") == 0)
        C_tW = value;
    else if (name.compare("C_tZ") == 0)
        C_tZ = value;
    else if (name.compare("C_tB") == 0)
        C_tB = value;
    else if (name.compare("C_tphi") == 0)
        C_tphi = value;
    else if (name.compare("C_phib") == 0)
        C_phib = value;
    else if (name.compare("C_bW") == 0)
        C_bW = value;
    else if (name.compare("C_bB") == 0)
        C_bB = value;
    else if (name.compare("C_bZ") == 0)
        C_bZ = value;
    else if (name.compare("C_phitb") == 0)
        C_phitb = value;
    else if (name.compare("C_tG") == 0)
        C_tG = value;
    else if (name.compare("C_ed") == 0)
        C_ed = value;
    else if (name.compare("C_eq") == 0)
        C_eq = value;
    else if (name.compare("C_ld") == 0)
        C_ld = value;
    else if (name.compare("C_lqP") == 0)
        C_lqP = value;
    else if (name.compare("C_eu") == 0)
        C_eu = value;
    else if (name.compare("C_lu") == 0)
        C_lu = value;
    else if (name.compare("C_lqM") == 0)
        C_lqM = value;
    else if (name.compare("C_tu8") == 0)
        C_tu8 = value;
    else if (name.compare("C_td8") == 0)
        C_td8 = value;
    else if (name.compare("C_Qq18") == 0)
        C_Qq18 = value;
    else if (name.compare("C_tq8") == 0)
        C_tq8 = value;
    else if (name.compare("C_Qq38") == 0)
        C_Qq38 = value;
    else if (name.compare("C_Qu8") == 0)
        C_Qu8 = value;
    else if (name.compare("C_Qd8") == 0)
        C_Qd8 = value;
    else if (name.compare("C_Qd1") == 0)
        C_Qd1 = value;
    else if (name.compare("C_Qu1") == 0)
        C_Qu1 = value;
    else if (name.compare("C_td1") == 0)
        C_td1 = value;
    else if (name.compare("C_tu1") == 0)
        C_tu1 = value;
    else if (name.compare("C_tq1") == 0)
        C_tq1 = value;
    else if (name.compare("C_Qq11") == 0)
        C_Qq11 = value;
    else if (name.compare("C_Qq31") == 0)
        C_Qq31 = value;
    
    

    
    else if (name.compare("C_phiG") == 0)
        C_phiG = value;
    else if (name.compare("C_phiGtil") == 0)
        C_phiGtil = value;
    else if (name.compare("C_phiW") == 0)
        C_phiW = value;
    else if (name.compare("C_phiWtil") == 0)
        C_phiWtil = value;
    else if (name.compare("C_tphiIm") == 0)
        C_tphiIm = value;
    else if (name.compare("C_tGIm") == 0)
        C_tGIm = value;
    else if (name.compare("C_tWIm") == 0)
        C_tWIm = value;
    else if (name.compare("C_tBIm") == 0)
        C_tBIm = value;
    else if (name.compare("C_tZIm") == 0)
        C_tZIm = value;
    
    else
        NPbase::setParameter(name, value);
}




// Flag to swith on/off the quadratic terms. flag = true means quadratic terms
//are switched on

bool NPSMEFT6dtopquark::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if(name.compare("Quadraticflag") == 0) {
        std::cout<<"Quadraticflag = "<< value <<std::endl;
        flag_Quadratic = value;
        res = true;
    }
    else if(name.compare("flag_LHC_WG_Basis") == 0) {
        std::cout<<"flag_LHC_WG_Basis = "<< value <<std::endl;
        flag_LHC_WG_Basis = value;
        res = true;
    }
    else
        res = StandardModel::setFlag(name,value);

    return(res);
}








/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////



////// ttbar Forward Backward asymmetry at Tevatron ///////////////////////////////////


FB_asymmetry_Tevatron_tt_diff_mtt_NPSMEFT6dtopquark::FB_asymmetry_Tevatron_tt_diff_mtt_NPSMEFT6dtopquark(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() 
            << "SM_FB_asymmetry_bin_mtt_350_450" << "SM_FB_asymmetry_deno_bin_mtt_350_450"
            << "SM_FB_asymmetry_bin_mtt_450_550" << "SM_FB_asymmetry_deno_bin_mtt_450_550"
            << "SM_FB_asymmetry_bin_mtt_550_650" << "SM_FB_asymmetry_deno_bin_mtt_550_650"
            << "SM_FB_asymmetry_bin_mtt_650_750" << "SM_FB_asymmetry_deno_bin_mtt_650_750");
            
    
}


double FB_asymmetry_Tevatron_tt_diff_mtt_NPSMEFT6dtopquark::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
        
    
    
    
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    
    
    
    if(b_min == 350 && b_max == 450){
        
        
        double SM_FB_asymmetry_bin_mtt_350_450 = SM.getOptionalParameter("SM_FB_asymmetry_bin_mtt_350_450");
        double SM_FB_asymmetry_deno_bin_mtt_350_450 = SM.getOptionalParameter("SM_FB_asymmetry_deno_bin_mtt_350_450");
        double SM_FB_asymmetry_num_bin_mtt_350_450 = SM_FB_asymmetry_bin_mtt_350_450*SM_FB_asymmetry_deno_bin_mtt_350_450;
        
        double SM_sigma_pos_bin_mtt_350_450 =0.5*(SM_FB_asymmetry_num_bin_mtt_350_450+SM_FB_asymmetry_deno_bin_mtt_350_450);
        double SM_sigma_neg_bin_mtt_350_450 =0.5*(SM_FB_asymmetry_deno_bin_mtt_350_450-SM_FB_asymmetry_num_bin_mtt_350_450);

                
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double sigma_pos_bin_mtt_350_450_madgraph = 1.708900;
            double sigma_pos_bin_mtt_350_450_NP = 0.014331*C_Qq11+0.020387*C_Qq31+0.057020*C_Qq38+
            0.008985*C_Qu1+0.075826*C_Qu8+-0.005241*C_td1+0.013312*C_td8+0.011673*C_tq1+0.075137*C_tu8+
            -0.000215*C_Qd1+0.015275*C_Qd8+0.089885*C_Qq18+0.461893*C_tG+0.088284*C_tq8+0.011721*C_tu1;
            double sigma_neg_bin_mtt_350_450_madgraph = 1.711870;
            double sigma_neg_bin_mtt_350_450_NP = 0.015137*C_Qq11+0.020418*C_Qq31+0.061615*C_Qq38+
            0.005954*C_Qu1+0.072189*C_Qu8+0.002100*C_td1+0.015387*C_td8+0.011721*C_tq1+0.072523*C_tu8+
            -0.000797*C_Qd1+0.013454*C_Qd8+0.088964*C_Qq18+0.457535*C_tG+0.090719*C_tq8+0.008783*C_tu1;
            
            double sigma_pos_bin_mtt_350_450_NP_Corrected = SM_sigma_pos_bin_mtt_350_450*sigma_pos_bin_mtt_350_450_NP/sigma_pos_bin_mtt_350_450_madgraph;
            double sigma_neg_bin_mtt_350_450_NP_Corrected = SM_sigma_neg_bin_mtt_350_450*sigma_neg_bin_mtt_350_450_NP/sigma_neg_bin_mtt_350_450_madgraph;
            
            double NP_FB_asymmetry_num_bin_mtt_350_450 = sigma_pos_bin_mtt_350_450_NP_Corrected - sigma_neg_bin_mtt_350_450_NP_Corrected;
            double NP_FB_asymmetry_deno_bin_mtt_350_450 = sigma_pos_bin_mtt_350_450_NP_Corrected + sigma_neg_bin_mtt_350_450_NP_Corrected;;
            
            //return  SM_FB_asymmetry_bin_mtt_350_450*(1+(NP_FB_asymmetry_num_bin_mtt_350_450-NP_FB_asymmetry_deno_bin_mtt_350_450)/SM_FB_asymmetry_deno_bin_mtt_350_450);            

            return SM_FB_asymmetry_bin_mtt_350_450 + (NP_FB_asymmetry_num_bin_mtt_350_450
                    - SM_FB_asymmetry_bin_mtt_350_450*NP_FB_asymmetry_deno_bin_mtt_350_450
                    )/SM_FB_asymmetry_deno_bin_mtt_350_450;
            
            
        }
        
    } else if(b_min == 450 && b_max == 550){
        
        
        double SM_FB_asymmetry_bin_mtt_450_550 = SM.getOptionalParameter("SM_FB_asymmetry_bin_mtt_450_550");
        double SM_FB_asymmetry_deno_bin_mtt_450_550 = SM.getOptionalParameter("SM_FB_asymmetry_deno_bin_mtt_450_550");
        double SM_FB_asymmetry_num_bin_mtt_450_550 = SM_FB_asymmetry_bin_mtt_450_550*SM_FB_asymmetry_deno_bin_mtt_450_550;
        
        double SM_sigma_pos_bin_mtt_450_550 =0.5*(SM_FB_asymmetry_num_bin_mtt_450_550+SM_FB_asymmetry_deno_bin_mtt_450_550);
        double SM_sigma_neg_bin_mtt_450_550 =0.5*(SM_FB_asymmetry_deno_bin_mtt_450_550-SM_FB_asymmetry_num_bin_mtt_450_550);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double sigma_pos_bin_mtt_450_550_madgraph = 0.700156;
            double sigma_pos_bin_mtt_450_550_NP = 0.000293*C_Qd1+0.006654*C_Qd8+0.056437*C_Qq18+
            0.202154*C_tG+0.057218*C_tq8+0.005943*C_tu1+0.012571*C_Qq11+0.014036*C_Qq31+0.041744*C_Qq38+
            0.004409*C_Qu1+0.048109*C_Qu8+-0.000686*C_td1+0.007816*C_td8+0.005690*C_tq1+0.048751*C_tu8;
            double sigma_neg_bin_mtt_450_550_madgraph = 0.697721;
            double sigma_neg_bin_mtt_450_550_NP = +-0.000320*C_Qd1+0.008928*C_Qd8+0.056860*C_Qq18+
            0.200120*C_tG+0.055877*C_tq8+0.008468*C_tu1+0.010415*C_Qq11+0.015990*C_Qq31+0.039785*C_Qq38+
            0.004021*C_Qu1+0.049076*C_Qu8+0.000080*C_td1+0.007790*C_td8+0.007160*C_tq1+0.048573*C_tu8;
            
            double sigma_pos_bin_mtt_450_550_NP_Corrected = SM_sigma_pos_bin_mtt_450_550*sigma_pos_bin_mtt_450_550_NP/sigma_pos_bin_mtt_450_550_madgraph;
            double sigma_neg_bin_mtt_450_550_NP_Corrected = SM_sigma_neg_bin_mtt_450_550*sigma_neg_bin_mtt_450_550_NP/sigma_neg_bin_mtt_450_550_madgraph;
            
            double NP_FB_asymmetry_num_bin_mtt_450_550 = sigma_pos_bin_mtt_450_550_NP_Corrected - sigma_neg_bin_mtt_450_550_NP_Corrected;
            double NP_FB_asymmetry_deno_bin_mtt_450_550 = sigma_pos_bin_mtt_450_550_NP_Corrected + sigma_neg_bin_mtt_450_550_NP_Corrected;;
            
            //return  SM_FB_asymmetry_bin_mtt_450_550*(1+(NP_FB_asymmetry_num_bin_mtt_450_550-NP_FB_asymmetry_deno_bin_mtt_450_550)/SM_FB_asymmetry_deno_bin_mtt_450_550);            
        
            
            return SM_FB_asymmetry_bin_mtt_450_550 + (NP_FB_asymmetry_num_bin_mtt_450_550
                    - SM_FB_asymmetry_bin_mtt_450_550*NP_FB_asymmetry_deno_bin_mtt_450_550
                    )/SM_FB_asymmetry_deno_bin_mtt_450_550;
            
            
        }
        
    } else if(b_min == 550 && b_max == 650){
        
        
        double SM_FB_asymmetry_bin_mtt_550_650 = SM.getOptionalParameter("SM_FB_asymmetry_bin_mtt_550_650");
        double SM_FB_asymmetry_deno_bin_mtt_550_650 = SM.getOptionalParameter("SM_FB_asymmetry_deno_bin_mtt_550_650");
        double SM_FB_asymmetry_num_bin_mtt_550_650 = SM_FB_asymmetry_bin_mtt_550_650*SM_FB_asymmetry_deno_bin_mtt_550_650;
        
        double SM_sigma_pos_bin_mtt_550_650 =0.5*(SM_FB_asymmetry_num_bin_mtt_550_650+SM_FB_asymmetry_deno_bin_mtt_550_650);
        double SM_sigma_neg_bin_mtt_550_650 =0.5*(SM_FB_asymmetry_deno_bin_mtt_550_650-SM_FB_asymmetry_num_bin_mtt_550_650);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double sigma_pos_bin_mtt_550_650_madgraph = 0.225585;
            double sigma_pos_bin_mtt_550_650_NP = 0.005195*C_Qq11+0.007938*C_Qq31+0.020820*C_Qq38+
            0.001724*C_Qu1+0.024030*C_Qu8+-0.000224*C_td1+0.003242*C_td8+0.002375*C_tq1+0.023554*C_tu8+
            0.000035*C_Qd1+0.003262*C_Qd8+0.026734*C_Qq18+0.068150*C_tG+0.027216*C_tq8+0.003859*C_tu1;
            double sigma_neg_bin_mtt_550_650_madgraph = 0.226014;
            double sigma_neg_bin_mtt_550_650_NP = 0.007480*C_Qq11+0.008103*C_Qq31+0.020939*C_Qq38+
            0.001998*C_Qu1+0.023779*C_Qu8+-0.000302*C_td1+0.003158*C_td8+0.002697*C_tq1+0.024310*C_tu8+
            -0.000427*C_Qd1+0.003031*C_Qd8+0.027516*C_Qq18+0.069044*C_tG+0.026921*C_tq8+0.003865*C_tu1;
            
            double sigma_pos_bin_mtt_550_650_NP_Corrected = SM_sigma_pos_bin_mtt_550_650*sigma_pos_bin_mtt_550_650_NP/sigma_pos_bin_mtt_550_650_madgraph;
            double sigma_neg_bin_mtt_550_650_NP_Corrected = SM_sigma_neg_bin_mtt_550_650*sigma_neg_bin_mtt_550_650_NP/sigma_neg_bin_mtt_550_650_madgraph;
            
            double NP_FB_asymmetry_num_bin_mtt_550_650 = sigma_pos_bin_mtt_550_650_NP_Corrected - sigma_neg_bin_mtt_550_650_NP_Corrected;
            double NP_FB_asymmetry_deno_bin_mtt_550_650 = sigma_pos_bin_mtt_550_650_NP_Corrected + sigma_neg_bin_mtt_550_650_NP_Corrected;;
            
            //return  SM_FB_asymmetry_bin_mtt_550_650*(1+(NP_FB_asymmetry_num_bin_mtt_550_650-NP_FB_asymmetry_deno_bin_mtt_550_650)/SM_FB_asymmetry_deno_bin_mtt_550_650);            
            
            return SM_FB_asymmetry_bin_mtt_550_650 + (NP_FB_asymmetry_num_bin_mtt_550_650
                    - SM_FB_asymmetry_bin_mtt_550_650*NP_FB_asymmetry_deno_bin_mtt_550_650
                    )/SM_FB_asymmetry_deno_bin_mtt_550_650;
            
        }
        
    } else if(b_min == 650 && b_max == 750){
        
        
        double SM_FB_asymmetry_bin_mtt_650_750 = SM.getOptionalParameter("SM_FB_asymmetry_bin_mtt_650_750");
        double SM_FB_asymmetry_deno_bin_mtt_650_750 = SM.getOptionalParameter("SM_FB_asymmetry_deno_bin_mtt_650_750");
        double SM_FB_asymmetry_num_bin_mtt_650_750 = SM_FB_asymmetry_bin_mtt_650_750*SM_FB_asymmetry_deno_bin_mtt_650_750;
        
        double SM_sigma_pos_bin_mtt_650_750 =0.5*(SM_FB_asymmetry_num_bin_mtt_650_750+SM_FB_asymmetry_deno_bin_mtt_650_750);
        double SM_sigma_neg_bin_mtt_650_750 =0.5*(SM_FB_asymmetry_deno_bin_mtt_650_750-SM_FB_asymmetry_num_bin_mtt_650_750);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double sigma_pos_bin_mtt_650_750_madgraph = 0.104697;
            double sigma_pos_bin_mtt_650_750_NP = +-0.000060*C_Qd1+0.001870*C_Qd8+0.020042*C_Qq18+
            0.033534*C_tG+0.020162*C_tq8+0.002968*C_tu1+0.005506*C_Qq11+0.006207*C_Qq31+0.016592*C_Qq38+
            0.001106*C_Qu1+0.018423*C_Qu8+-0.000255*C_td1+0.001826*C_td8+0.001752*C_tq1+0.018350*C_tu8;
            double sigma_neg_bin_mtt_650_750_madgraph = 0.104541;
            double sigma_neg_bin_mtt_650_750_NP = -0.000081*C_Qd1+0.001783*C_Qd8+0.020151*C_Qq18+
            0.033203*C_tG+0.020170*C_tq8+0.003204*C_tu1+0.005148*C_Qq11+0.006367*C_Qq31+0.016418*C_Qq38+
            0.001303*C_Qu1+0.018244*C_Qu8+-0.000097*C_td1+0.001814*C_td8+0.001599*C_tq1+0.018284*C_tu8;
            
            double sigma_pos_bin_mtt_650_750_NP_Corrected = SM_sigma_pos_bin_mtt_650_750*sigma_pos_bin_mtt_650_750_NP/sigma_pos_bin_mtt_650_750_madgraph;
            double sigma_neg_bin_mtt_650_750_NP_Corrected = SM_sigma_neg_bin_mtt_650_750*sigma_neg_bin_mtt_650_750_NP/sigma_neg_bin_mtt_650_750_madgraph;
            
            double NP_FB_asymmetry_num_bin_mtt_650_750 = sigma_pos_bin_mtt_650_750_NP_Corrected - sigma_neg_bin_mtt_650_750_NP_Corrected;
            double NP_FB_asymmetry_deno_bin_mtt_650_750 = sigma_pos_bin_mtt_650_750_NP_Corrected + sigma_neg_bin_mtt_650_750_NP_Corrected;;
            
            //return  SM_FB_asymmetry_bin_mtt_650_750*(1+(NP_FB_asymmetry_num_bin_mtt_650_750-NP_FB_asymmetry_deno_bin_mtt_650_750)/SM_FB_asymmetry_deno_bin_mtt_650_750);            
        
            return SM_FB_asymmetry_bin_mtt_650_750 + (NP_FB_asymmetry_num_bin_mtt_650_750
                    - SM_FB_asymmetry_bin_mtt_650_750*NP_FB_asymmetry_deno_bin_mtt_650_750
                    )/SM_FB_asymmetry_deno_bin_mtt_650_750;
            
            
        }
        
    } else {
        throw std::runtime_error("\nERROR: Please specify a correct binning range for FB_asymmetry_Tevatron_tt_diff_mtt_LO.\n");
    }
    
    
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Observables from LEP ////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

Rb_NPSMEFT6dtopquark::Rb_NPSMEFT6dtopquark(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "Rb_SM" );

}

double Rb_NPSMEFT6dtopquark::computeThValue()
{
    double Rb_SM = SM.getOptionalParameter("Rb_SM");
    double lep_bb_madgraph = 0.22;
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return Rb_SM + (0.023*(C_phiQ3/0.998)+0.023*((C_phiQm+C_phiQ3)/0.998)-0.005*(C_phib/(0.998))+0.0018*C_bW*(-1)*C_bW*(-1)+0.002*C_bW*(-1)*((0.881533*C_bW-C_bZ)/(0.472123))*(-1))*(Rb_SM/lep_bb_madgraph);
            }
            else{
                return Rb_SM + (0.023*(C_phiQ3/0.998)+0.023*((C_phiQm+C_phiQ3)/0.998)-0.005*(C_phib/(0.998)))*(Rb_SM/lep_bb_madgraph);
            }
    }
    else{
            if(flag_Quadratic){
                return Rb_SM + (0.023*(C_phiQ3/0.998)+0.023*((C_phiQ1)/0.998)-0.005*(C_phib/(0.998))+0.0018*C_bW*(-1)*C_bW*(-1)+0.002*C_bW*(-1)*(C_bB)*(-1))*(Rb_SM/lep_bb_madgraph);
            }
            else{
                return Rb_SM + (0.023*(C_phiQ3/0.998)+0.023*((C_phiQ1)/0.998)-0.005*(C_phib/(0.998)))*(Rb_SM/lep_bb_madgraph);
            }
    }
}


AFBLR::AFBLR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "AFBLR_SM" );
}

double AFBLR::computeThValue()
{
    double AFBLR_SM = SM.getOptionalParameter("AFBLR_SM");
    double lep_alr_madgraph = 0.66;
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return AFBLR_SM + (0.008*(C_phiQ3/0.998)+0.008*((C_phiQm+C_phiQ3)/0.998)+0.034*(C_phib/(0.998))+0.0056*C_bW*(-1)*C_bW*(-1)-0.002*(C_phib/(0.998))*((0.881533*C_bW-C_bZ)/(0.472123))*(-1)+0.0015*((0.881533*C_bW-C_bZ)/(0.472123))*(-1)*((0.881533*C_bW-C_bZ)/(0.472123))*(-1)+0.0076*C_bW*(-1)*((0.881533*C_bW-C_bZ)/(0.472123))*(-1))*(AFBLR_SM/lep_alr_madgraph);
            }
            else{
                return AFBLR_SM + (0.008*(C_phiQ3/0.998)+0.008*((C_phiQm+C_phiQ3)/0.998)+0.034*(C_phib/(0.998)))*(AFBLR_SM/lep_alr_madgraph);
            }    
    }
    else{
            if(flag_Quadratic){
                return AFBLR_SM + (0.008*(C_phiQ3/0.998)+0.008*((C_phiQ1)/0.998)+0.034*(C_phib/(0.998))+0.0056*C_bW*(-1)*C_bW*(-1)-0.002*(C_phib/(0.998))*(C_bB)*(-1)+0.0015*(C_bB)*(-1)*(C_bB)*(-1)+0.0076*C_bW*(-1)*(C_bB)*(-1))*(AFBLR_SM/lep_alr_madgraph);
            }
            else{
                return AFBLR_SM + (0.008*(C_phiQ3/0.998)+0.008*((C_phiQ1)/0.998)+0.034*(C_phib/(0.998)))*(AFBLR_SM/lep_alr_madgraph);
            }    
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////// Observables from Tevatron ////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


//// ttbar xsection

sigmattbarTev::sigmattbarTev(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_ttbar_tev" );
}

double sigmattbarTev::computeThValue()
{
    
    double SM_ttbar_tev = SM.getOptionalParameter("SM_ttbar_tev");
    double   xttbarTev_madgraph_LO = 4948.8;//fb
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double mC_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double mC_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();


    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return SM_ttbar_tev + (0.303*mC_tu8+ 0.051*mC_td8 +1.378*C_tG + 0.022*mC_tu8*mC_tu8+ 0.004*mC_td8*mC_td8+0.123*C_tG*C_tG + 0.381130*C_Qq18 + 0.380930*C_tq8 + 0.271890*C_Qq38 + 0.326390*C_Qu8 + 0.054694*C_Qd8)*(1000*SM_ttbar_tev/xttbarTev_madgraph_LO);
            }
            else{
                return SM_ttbar_tev + (0.303*mC_tu8+ 0.051*mC_td8 +1.378*C_tG + 0.381130*C_Qq18 + 0.380930*C_tq8 + 0.271890*C_Qq38 + 0.326390*C_Qu8 + 0.054694*C_Qd8)*(1000*SM_ttbar_tev/xttbarTev_madgraph_LO);
            }    
    }   
    else{
            if(flag_Quadratic){
                return SM_ttbar_tev + (0.303*mC_tu8+ 0.051*mC_td8 +1.378*C_tG + 0.022*mC_tu8*mC_tu8+ 0.004*mC_td8*mC_td8+0.123*C_tG*C_tG)*(1000*SM_ttbar_tev/xttbarTev_madgraph_LO);
            }
            else{
                return SM_ttbar_tev + (0.303*mC_tu8+ 0.051*mC_td8 +1.378*C_tG )*(1000*SM_ttbar_tev/xttbarTev_madgraph_LO);
            }    
    }
}


////// single production s-channel

sigmaschannelTev::sigmaschannelTev(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigmaschannelTev" );
}

double sigmaschannelTev::computeThValue()
    {
    double SM_sigmaschannelTev = SM.getOptionalParameter("SM_sigmaschannelTev");
    double xschan_madgraph_LO = 659;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_phitb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();


    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_sigmaschannelTev + (0.30395*C_tW+0.07974*C_phiQ3-0.00023*C_tG
                        +0.04789*C_tW*C_tW+0.00268*C_phiQ3*C_phiQ3+0.00758*C_tG*C_tG+0.02019*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532))+0.0006*(-C_phitb/(0.998))*(-C_phitb/(0.998)))*(1000*SM_sigmaschannelTev/xschan_madgraph_LO) ;      
            }
            else{
                return SM_sigmaschannelTev + (0.30395*C_tW+0.07974*C_phiQ3-0.00023*C_tG)*(1000*SM_sigmaschannelTev/xschan_madgraph_LO) ;      
            }      
    }
    
    else{
            if(flag_Quadratic){
                return  SM_sigmaschannelTev + (0.30395*C_tW+0.07974*C_phiQ3-0.00023*C_tG
                        +0.04789*C_tW*C_tW+0.00268*C_phiQ3*C_phiQ3+0.00758*C_tG*C_tG+0.02019*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532))+0.0006*(-C_phitb/(0.998))*(-C_phitb/(0.998)))*(1000*SM_sigmaschannelTev/xschan_madgraph_LO) ;      
            }
            else{
                return SM_sigmaschannelTev + (0.30395*C_tW+0.07974*C_phiQ3-0.00023*C_tG)*(1000*SM_sigmaschannelTev/xschan_madgraph_LO) ;      
            }      
    }
}




////// Forward-Backward asymmetry ttbar

FB_asymmetry_Tevatron_tt_diff_mtt_top_basis_LO::FB_asymmetry_Tevatron_tt_diff_mtt_top_basis_LO(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() 
            << "SM_FB_asymmetry_bin_mtt_350_450" << "SM_FB_asymmetry_deno_bin_mtt_350_450"
            << "SM_FB_asymmetry_bin_mtt_450_550" << "SM_FB_asymmetry_deno_bin_mtt_450_550"
            << "SM_FB_asymmetry_bin_mtt_550_650" << "SM_FB_asymmetry_deno_bin_mtt_550_650"
            << "SM_FB_asymmetry_bin_mtt_650_750" << "SM_FB_asymmetry_deno_bin_mtt_650_750");
            
    
}


double FB_asymmetry_Tevatron_tt_diff_mtt_top_basis_LO::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    bool   flag_Quadratic= false; //Needs to be properly defined
        
    
    
    if(b_min == 350 && b_max == 450){
        
        
        double SM_FB_asymmetry_bin_mtt_350_450 = SM.getOptionalParameter("SM_FB_asymmetry_bin_mtt_350_450");
        double SM_FB_asymmetry_deno_bin_mtt_350_450 = SM.getOptionalParameter("SM_FB_asymmetry_deno_bin_mtt_350_450");
        double SM_FB_asymmetry_num_bin_mtt_350_450 = SM_FB_asymmetry_bin_mtt_350_450*SM_FB_asymmetry_deno_bin_mtt_350_450;
        
        double SM_sigma_pos_bin_mtt_350_450 =0.5*(SM_FB_asymmetry_num_bin_mtt_350_450+SM_FB_asymmetry_deno_bin_mtt_350_450);
        double SM_sigma_neg_bin_mtt_350_450 =0.5*(SM_FB_asymmetry_deno_bin_mtt_350_450-SM_FB_asymmetry_num_bin_mtt_350_450);

                
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double sigma_pos_bin_mtt_350_450_madgraph = 1.675650;
            double sigma_pos_bin_mtt_350_450_NP = 0.;
            double sigma_neg_bin_mtt_350_450_madgraph = 1.674510;
            double sigma_neg_bin_mtt_350_450_NP = 0.;
            
            double sigma_pos_bin_mtt_350_450_NP_Corrected = SM_sigma_pos_bin_mtt_350_450*sigma_pos_bin_mtt_350_450_NP/sigma_pos_bin_mtt_350_450_madgraph;
            double sigma_neg_bin_mtt_350_450_NP_Corrected = SM_sigma_neg_bin_mtt_350_450*sigma_neg_bin_mtt_350_450_NP/sigma_neg_bin_mtt_350_450_madgraph;
            
            double NP_FB_asymmetry_num_bin_mtt_350_450 = sigma_pos_bin_mtt_350_450_NP_Corrected - sigma_neg_bin_mtt_350_450_NP_Corrected;
            double NP_FB_asymmetry_deno_bin_mtt_350_450 = sigma_pos_bin_mtt_350_450_NP_Corrected + sigma_neg_bin_mtt_350_450_NP_Corrected;;
            
            return  SM_FB_asymmetry_bin_mtt_350_450*(1+(NP_FB_asymmetry_num_bin_mtt_350_450-NP_FB_asymmetry_deno_bin_mtt_350_450)/SM_FB_asymmetry_deno_bin_mtt_350_450);            
        
        }
        
    } else if(b_min == 450 && b_max == 550){
        
        
        double SM_FB_asymmetry_bin_mtt_450_550 = SM.getOptionalParameter("SM_FB_asymmetry_bin_mtt_450_550");
        double SM_FB_asymmetry_deno_bin_mtt_450_550 = SM.getOptionalParameter("SM_FB_asymmetry_deno_bin_mtt_450_550");
        double SM_FB_asymmetry_num_bin_mtt_450_550 = SM_FB_asymmetry_bin_mtt_450_550*SM_FB_asymmetry_deno_bin_mtt_450_550;
        
        double SM_sigma_pos_bin_mtt_450_550 =0.5*(SM_FB_asymmetry_num_bin_mtt_450_550+SM_FB_asymmetry_deno_bin_mtt_450_550);
        double SM_sigma_neg_bin_mtt_450_550 =0.5*(SM_FB_asymmetry_deno_bin_mtt_450_550-SM_FB_asymmetry_num_bin_mtt_450_550);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double sigma_pos_bin_mtt_450_550_madgraph = 0.683803;
            double sigma_pos_bin_mtt_450_550_NP = 0.;
            double sigma_neg_bin_mtt_450_550_madgraph = 0.685429;
            double sigma_neg_bin_mtt_450_550_NP = 0.;
            
            double sigma_pos_bin_mtt_450_550_NP_Corrected = SM_sigma_pos_bin_mtt_450_550*sigma_pos_bin_mtt_450_550_NP/sigma_pos_bin_mtt_450_550_madgraph;
            double sigma_neg_bin_mtt_450_550_NP_Corrected = SM_sigma_neg_bin_mtt_450_550*sigma_neg_bin_mtt_450_550_NP/sigma_neg_bin_mtt_450_550_madgraph;
            
            double NP_FB_asymmetry_num_bin_mtt_450_550 = sigma_pos_bin_mtt_450_550_NP_Corrected - sigma_neg_bin_mtt_450_550_NP_Corrected;
            double NP_FB_asymmetry_deno_bin_mtt_450_550 = sigma_pos_bin_mtt_450_550_NP_Corrected + sigma_neg_bin_mtt_450_550_NP_Corrected;;
            
            return  SM_FB_asymmetry_bin_mtt_450_550*(1+(NP_FB_asymmetry_num_bin_mtt_450_550-NP_FB_asymmetry_deno_bin_mtt_450_550)/SM_FB_asymmetry_deno_bin_mtt_450_550);            
        
        }
        
    } else if(b_min == 550 && b_max == 650){
        
        
        double SM_FB_asymmetry_bin_mtt_550_650 = SM.getOptionalParameter("SM_FB_asymmetry_bin_mtt_550_650");
        double SM_FB_asymmetry_deno_bin_mtt_550_650 = SM.getOptionalParameter("SM_FB_asymmetry_deno_bin_mtt_550_650");
        double SM_FB_asymmetry_num_bin_mtt_550_650 = SM_FB_asymmetry_bin_mtt_550_650*SM_FB_asymmetry_deno_bin_mtt_550_650;
        
        double SM_sigma_pos_bin_mtt_550_650 =0.5*(SM_FB_asymmetry_num_bin_mtt_550_650+SM_FB_asymmetry_deno_bin_mtt_550_650);
        double SM_sigma_neg_bin_mtt_550_650 =0.5*(SM_FB_asymmetry_deno_bin_mtt_550_650-SM_FB_asymmetry_num_bin_mtt_550_650);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double sigma_pos_bin_mtt_550_650_madgraph = 0.216436;
            double sigma_pos_bin_mtt_550_650_NP = 0.;
            double sigma_neg_bin_mtt_550_650_madgraph = 0.216949;
            double sigma_neg_bin_mtt_550_650_NP = 0.;
            
            double sigma_pos_bin_mtt_550_650_NP_Corrected = SM_sigma_pos_bin_mtt_550_650*sigma_pos_bin_mtt_550_650_NP/sigma_pos_bin_mtt_550_650_madgraph;
            double sigma_neg_bin_mtt_550_650_NP_Corrected = SM_sigma_neg_bin_mtt_550_650*sigma_neg_bin_mtt_550_650_NP/sigma_neg_bin_mtt_550_650_madgraph;
            
            double NP_FB_asymmetry_num_bin_mtt_550_650 = sigma_pos_bin_mtt_550_650_NP_Corrected - sigma_neg_bin_mtt_550_650_NP_Corrected;
            double NP_FB_asymmetry_deno_bin_mtt_550_650 = sigma_pos_bin_mtt_550_650_NP_Corrected + sigma_neg_bin_mtt_550_650_NP_Corrected;;
            
            return  SM_FB_asymmetry_bin_mtt_550_650*(1+(NP_FB_asymmetry_num_bin_mtt_550_650-NP_FB_asymmetry_deno_bin_mtt_550_650)/SM_FB_asymmetry_deno_bin_mtt_550_650);            
        
        }
        
    } else if(b_min == 650 && b_max == 750){
        
        
        double SM_FB_asymmetry_bin_mtt_650_750 = SM.getOptionalParameter("SM_FB_asymmetry_bin_mtt_650_750");
        double SM_FB_asymmetry_deno_bin_mtt_650_750 = SM.getOptionalParameter("SM_FB_asymmetry_deno_bin_mtt_650_750");
        double SM_FB_asymmetry_num_bin_mtt_650_750 = SM_FB_asymmetry_bin_mtt_650_750*SM_FB_asymmetry_deno_bin_mtt_650_750;
        
        double SM_sigma_pos_bin_mtt_650_750 =0.5*(SM_FB_asymmetry_num_bin_mtt_650_750+SM_FB_asymmetry_deno_bin_mtt_650_750);
        double SM_sigma_neg_bin_mtt_650_750 =0.5*(SM_FB_asymmetry_deno_bin_mtt_650_750-SM_FB_asymmetry_num_bin_mtt_650_750);

        
        if(flag_Quadratic){
        
            return  0.;
            
        }
        else{
            
            
            
            double sigma_pos_bin_mtt_650_750_madgraph = 0.097727;
            double sigma_pos_bin_mtt_650_750_NP = 0.;
            double sigma_neg_bin_mtt_650_750_madgraph = 0.097766;
            double sigma_neg_bin_mtt_650_750_NP = 0.;
            
            double sigma_pos_bin_mtt_650_750_NP_Corrected = SM_sigma_pos_bin_mtt_650_750*sigma_pos_bin_mtt_650_750_NP/sigma_pos_bin_mtt_650_750_madgraph;
            double sigma_neg_bin_mtt_650_750_NP_Corrected = SM_sigma_neg_bin_mtt_650_750*sigma_neg_bin_mtt_650_750_NP/sigma_neg_bin_mtt_650_750_madgraph;
            
            double NP_FB_asymmetry_num_bin_mtt_650_750 = sigma_pos_bin_mtt_650_750_NP_Corrected - sigma_neg_bin_mtt_650_750_NP_Corrected;
            double NP_FB_asymmetry_deno_bin_mtt_650_750 = sigma_pos_bin_mtt_650_750_NP_Corrected + sigma_neg_bin_mtt_650_750_NP_Corrected;;
            
            return  SM_FB_asymmetry_bin_mtt_650_750*(1+(NP_FB_asymmetry_num_bin_mtt_650_750-NP_FB_asymmetry_deno_bin_mtt_650_750)/SM_FB_asymmetry_deno_bin_mtt_650_750);            
        
        }
        
    } else {
        throw std::runtime_error("\nERROR: Please specify a correct binning range for FB_asymmetry_Tevatron_tt_diff_mtt_LO.\n");
    }
    
    
}







/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Observables from LHC ////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////



/////// Helicities  /////////////////////////////////////////////////////////////////////////////////////


F0::F0(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "F0_SM" );
}

double F0::computeThValue()
{
    double F0_SM = SM.getOptionalParameter("F0_SM");
    double F0_madgraph = 0.70381;
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  F0_SM + (-0.0569231*C_tW + 0.0022956*C_tW*C_tW -
                        0.00215306*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532))
                        )*(F0_SM/F0_madgraph);
            }
            else{
                return F0_SM + (-0.0569231*(C_tW) )*(F0_SM/F0_madgraph);
            }
    }
    
    else{
            if(flag_Quadratic){
                return  F0_SM + (-0.0569231*C_tW + 0.0022956*C_tW*C_tW -
                        0.00215306*((C_bW/(0.999*0.6532)))*((C_bW/(0.999*0.6532)))
                        )*(F0_SM/F0_madgraph);
            }
            else{
                return F0_SM + (-0.0569231*(C_tW) )*(F0_SM/F0_madgraph);
            }
        }
}


FL::FL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "FL_SM" );
}

double FL::computeThValue()
{
    double FL_SM = SM.getOptionalParameter("FL_SM");
    double FL_madgraph = 0.29619;
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_phitb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  FL_SM + ( 0.0569227*C_tW - 0.0022956*C_tW*C_tW 
                        - 0.0011115*((C_bW/(0.999*0.6532)))*((C_bW/(0.999*0.6532)))
                        - 0.000105115*((-C_phitb/(0.998)))*((-C_phitb/(0.998)))
                        )*(FL_SM/FL_madgraph);
            }
            else{
                return FL_SM + (0.0569227*(C_tW) )*(FL_SM/FL_madgraph);
            }
    }
    
    else{
            if(flag_Quadratic){
                return  FL_SM + ( 0.0569227*C_tW - 0.0022956*C_tW*C_tW 
                        - 0.0011115*((C_bW/(0.999*0.6532)))*((C_bW/(0.999*0.6532)))
                        - 0.000105115*((-C_phitb/(0.998)))*((-C_phitb/(0.998)))
                        )*(FL_SM/FL_madgraph);
            }
            else{
                return FL_SM + (0.0569227*(C_tW) )*(FL_SM/FL_madgraph);
            }
    }
    
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////




/////// Inclusive cross sections ////////////////////////////////////////////////////////////////////////

sigmattbarLHC13::sigmattbarLHC13(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_ttbar_LHC13" );
}

double sigmattbarLHC13::computeThValue()
{
    double SM_ttbar_LHC13 = SM.getOptionalParameter("SM_ttbar_LHC13");
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double   xttbarLHC13_madgraph_LO = 453100;//fb
    double   xttbarLHC13_madgraph_NLO = 661300;//fb
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double mC_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double mC_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    


    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return SM_ttbar_LHC13+ (3.950*mC_tu8+2.719*mC_td8+1.007*mC_tu8*mC_tu8+0.558*mC_td8*mC_td8 + 7.100600*C_Qq18 + 7.107300*C_tq8 + 1.438000*C_Qq38 + 4.261800*C_Qu8 + 2.904600*C_Qd8)*(1000*SM_ttbar_LHC13/xttbarLHC13_madgraph_LO)+(213.4*C_tG+31.2*C_tG*C_tG)*(SM_ttbar_LHC13/xttbarLHC13_madgraph_NLO);
            }
            else{
                return SM_ttbar_LHC13+ (3.950*mC_tu8+2.719*mC_td8 + 7.100600*C_Qq18 + 7.107300*C_tq8 + 1.438000*C_Qq38 + 4.261800*C_Qu8 + 2.904600*C_Qd8)*(SM_ttbar_LHC13/xttbarLHC13_madgraph_LO)+(213.4*C_tG)*(1000*SM_ttbar_LHC13/xttbarLHC13_madgraph_NLO);
            }    
    }
    
    else{
            if(flag_Quadratic){
                return SM_ttbar_LHC13+ (3.950*mC_tu8+2.719*mC_td8+1.007*mC_tu8*mC_tu8+0.558*mC_td8*mC_td8)*(1000*SM_ttbar_LHC13/xttbarLHC13_madgraph_LO)+(213.4*C_tG+31.2*C_tG*C_tG)*(SM_ttbar_LHC13/xttbarLHC13_madgraph_NLO);
            }
            else{
                return SM_ttbar_LHC13+ (3.950*mC_tu8+2.719*mC_td8)*(SM_ttbar_LHC13/xttbarLHC13_madgraph_LO)+(213.4*C_tG)*(1000*SM_ttbar_LHC13/xttbarLHC13_madgraph_NLO);
            }    

    }
    
      
}



sigmattbarLHC8::sigmattbarLHC8(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_ttbar_LHC8" );
}

double sigmattbarLHC8::computeThValue()
{
    
    double SM_ttbar_LHC8 = SM.getOptionalParameter("SM_ttbar_LHC8");
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double   xttbarLHC8_madgraph_LO = 138910;//fb
    double   xttbarLHC8_madgraph_NLO = 202500;//fb
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double mC_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double mC_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();

    

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return SM_ttbar_LHC8 + (1.762*mC_tu8+ 1.064*mC_td8 + 0.211*mC_tu8*mC_tu8+ 0.136*mC_td8*mC_td8 + 2.964900*C_Qq18 + 2.963700*C_tq8 + 0.694830*C_Qq38 + 1.828200*C_Qu8 + 1.152000*C_Qd8)*(1000*SM_ttbar_LHC8/xttbarLHC8_madgraph_LO)+ (64.0*C_tG+11.2*C_tG*C_tG)*(SM_ttbar_LHC8/xttbarLHC8_madgraph_NLO);
            }
            else{
                return SM_ttbar_LHC8 + (1.762*mC_tu8+ 1.064*mC_td8 + 2.964900*C_Qq18 + 2.963700*C_tq8 + 0.694830*C_Qq38 + 1.828200*C_Qu8 + 1.152000*C_Qd8)*(SM_ttbar_LHC8/xttbarLHC8_madgraph_LO)+ (64.0*C_tG)*(1000*SM_ttbar_LHC8/xttbarLHC8_madgraph_NLO);
            }    
    }
    
    else{
            if(flag_Quadratic){
                return SM_ttbar_LHC8 + (1.762*mC_tu8+ 1.064*mC_td8 + 0.211*mC_tu8*mC_tu8+ 0.136*mC_td8*mC_td8)*(1000*SM_ttbar_LHC8/xttbarLHC8_madgraph_LO)+ (64.0*C_tG+11.2*C_tG*C_tG)*(SM_ttbar_LHC8/xttbarLHC8_madgraph_NLO);
            }
            else{
                return SM_ttbar_LHC8 + (1.762*mC_tu8+ 1.064*mC_td8 )*(SM_ttbar_LHC8/xttbarLHC8_madgraph_LO)+ (64.0*C_tG)*(1000*SM_ttbar_LHC8/xttbarLHC8_madgraph_NLO);
            }    
    }
      
}



sigmattZ::sigmattZ(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_ttZ_inc" );
}

double sigmattZ::computeThValue()
{
    double SM_ttZ_inc = SM.getOptionalParameter("SM_ttZ_inc");
    double xttzNLO_madgraph = 740;//fb
    double xttzLO_madgraph = 682.;//fb
    //double xttzLO_madgraph = 510.23;//fb
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    
    
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic) {
                return  SM_ttZ_inc+ (-3.25*C_tZ+ 288.50*C_tG + 48.59*C_phit -76.50*C_phiQm 
                        + 62.33*C_tZ*C_tZ+ 291.80*C_tG*C_tG +3.20*C_phit*C_phit
                        +3.75*C_phiQm*C_phiQm-0.0016*C_phit*C_phiQm)*(SM_ttZ_inc/(xttzNLO_madgraph)) 
                        +(-0.017852*C_Qq38+0.031201*C_Qq38*C_Qq38+0.011527*C_tu8+0.005967*C_tu8*C_tu8+0.067759*C_Qq18+0.031474*C_Qq18*C_Qq18+0.009777*C_td8+0.003648*C_td8*C_td8+0.015271*C_Qd8+0.007609*C_Qd8*C_Qd8+0.017474*C_Qu8+0.011569*C_Qu8*C_Qu8+0.055888*C_tq8+0.021441*C_tq8*C_tq8)*(1000*SM_ttZ_inc/xttzLO_madgraph);

            }
            else{
                return  SM_ttZ_inc + (-3.25*C_tZ+288.50*C_tG+48.59*C_phit
                        -76.50*C_phiQm)*(SM_ttZ_inc/(xttzNLO_madgraph))+
                        (-0.017852*C_Qq38+0.011527*C_tu8+0.067759*C_Qq18+0.009777*C_td8+0.015271*C_Qd8+0.017474*C_Qu8+0.055888*C_tq8+0.002146*C_phiQ3)*(1000*SM_ttZ_inc/xttzLO_madgraph);
            }
    }
    else{
        if(flag_Quadratic){
                return  SM_ttZ_inc+ (-3.25*(-0.472123*C_tB + 0.881533*C_tW)+ 288.50*C_tG + 48.59*C_phit -76.50*(C_phiQ1-C_phiQ3) 
                        + 62.33*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW)+ 291.80*C_tG*C_tG +3.20*C_phit*C_phit
                        +3.75*(C_phiQ1-C_phiQ3)*(C_phiQ1-C_phiQ3)-0.0016*C_phit*(C_phiQ1-C_phiQ3))*(SM_ttZ_inc/(xttzNLO_madgraph)) 
                        +(-0.017852*C_Qq38+0.031201*C_Qq38*C_Qq38+0.011527*C_tu8+0.005967*C_tu8*C_tu8+0.067759*C_Qq18+0.031474*C_Qq18*C_Qq18+0.009777*C_td8+0.003648*C_td8*C_td8+0.015271*C_Qd8+0.007609*C_Qd8*C_Qd8+0.017474*C_Qu8+0.011569*C_Qu8*C_Qu8+0.055888*C_tq8+0.021441*C_tq8*C_tq8)*(1000*SM_ttZ_inc/xttzLO_madgraph);
            }
            else{
                return  SM_ttZ_inc + (-3.25*(-0.472123*C_tB + 0.881533*C_tW)+288.50*C_tG+48.59*C_phit
                        -76.50*(C_phiQ1-C_phiQ3))*(SM_ttZ_inc/(xttzNLO_madgraph))+
                        (-0.017852*C_Qq38+0.011527*C_tu8+0.067759*C_Qq18+0.009777*C_td8+0.015271*C_Qd8+0.017474*C_Qu8+0.055888*C_tq8+0.002146*C_phiQ3)*(1000*SM_ttZ_inc/xttzLO_madgraph);
            }
    }
}


sigmattA::sigmattA(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_ttA_inc" );
}

double sigmattA::computeThValue()
{   
    double SM_ttA_inc = SM.getOptionalParameter("SM_ttA_inc");
    double xtta_madgraph_NLO = 1986;//fb
    double xtta_madgraph_LO = 2511;//fb
    //double xtta_madgraph_LO = 1317.6;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttA_inc + (-88.4173*C_tZ + 110.514*C_tW + 578*C_tG  + 291.96*C_tZ*C_tZ + 380.115*C_tW*C_tW +152.0*C_tG*C_tG - 660.969*C_tZ*C_tW)*(SM_ttA_inc/xtta_madgraph_NLO)+
                        (0.029197*C_Qq38+0.039433*C_Qq38*C_Qq38+0.066549*C_tu8+0.023649*C_tu8*C_tu8+0.103103*C_Qq18+0.039474*C_Qq18*C_Qq18+0.036565*C_td8+0.014590*C_td8*C_td8+0.019167*C_Qd8+0.007880*C_Qd8*C_Qd8+0.119846*C_Qu8+0.049185*C_Qu8*C_Qu8+0.139730*C_tq8+0.058076*C_tq8*C_tq8)*(1000*SM_ttA_inc/xtta_madgraph_LO);

            }
            else{
                return  SM_ttA_inc + (-88.4173*C_tZ + 110.514*C_tW + 578*C_tG  )*(SM_ttA_inc/xtta_madgraph_NLO)+
                        (0.029197*C_Qq38+0.066549*C_tu8+0.103103*C_Qq18+0.036565*C_td8+0.019167*C_Qd8+0.119846*C_Qu8+0.139730*C_tq8)*(1000*SM_ttA_inc/xtta_madgraph_LO); 
            }
    }    
    
    else{
        if(flag_Quadratic){
            return  SM_ttA_inc + (-88.4173*(-0.472123*C_tB + 0.881533*C_tW) + 110.514*C_tW + 578*C_tG  + 291.96*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 380.115*C_tW*C_tW +152.0*C_tG*C_tG - 660.969*(-0.472123*C_tB + 0.881533*C_tW)*C_tW)*(SM_ttA_inc/xtta_madgraph_NLO)+
                        (0.029197*C_Qq38+0.039433*C_Qq38*C_Qq38+0.066549*C_tu8+0.023649*C_tu8*C_tu8+0.103103*C_Qq18+0.039474*C_Qq18*C_Qq18+0.036565*C_td8+0.014590*C_td8*C_td8+0.019167*C_Qd8+0.007880*C_Qd8*C_Qd8+0.119846*C_Qu8+0.049185*C_Qu8*C_Qu8+0.139730*C_tq8+0.058076*C_tq8*C_tq8)*(1000*SM_ttA_inc/xtta_madgraph_LO);
        }
        else{
            return  SM_ttA_inc + (-88.4173*(-0.472123*C_tB + 0.881533*C_tW) + 110.514*C_tW + 578*C_tG  )*(SM_ttA_inc/xtta_madgraph_NLO)+
                        (0.029197*C_Qq38+0.066549*C_tu8+0.103103*C_Qq18+0.036565*C_td8+0.019167*C_Qd8+0.119846*C_Qu8+0.139730*C_tq8)*(1000*SM_ttA_inc/xtta_madgraph_LO); 
        }
    }
        
}


sigmattW::sigmattW(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_ttW_inc" );
}

double sigmattW::computeThValue()
{
    double SM_ttW_inc = SM.getOptionalParameter("SM_ttW_inc");
    double NLOxttw_madgraph = 543.3;//fb
    double LOxttw_madgraph = 407.2;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();



    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttW_inc + (158.3*C_tG  + 23.4*C_tG*C_tG)*(SM_ttW_inc/NLOxttw_madgraph)+
                        (119.3*C_Qq18 + 42.0*C_Qq18*C_Qq18 - 43.7*C_Qq38 + 95.8*C_Qq38*C_Qq38 + 119.4*C_tq8 + 42.1*C_tq8*C_tq8 - 31.9*C_Qq18*C_Qq38 + 7.1*C_Qq18*C_tq8 - 3.4*C_Qq38*C_tq8  )*(SM_ttW_inc/LOxttw_madgraph);
            }
            else{
                return  (SM_ttW_inc + (158.3*C_tG)*(SM_ttW_inc/NLOxttw_madgraph))+
                        (119.3*C_Qq18 - 43.7*C_Qq38 + 119.4*C_tq8 )*(SM_ttW_inc/LOxttw_madgraph);
            }      
    }
    
    else{
            if(flag_Quadratic){
                return  (SM_ttW_inc + (158.3*C_tG  + 23.4*C_tG*C_tG)*(SM_ttW_inc/NLOxttw_madgraph))+
                        (119.3*C_Qq18 + 42.0*C_Qq18*C_Qq18 - 43.7*C_Qq38 + 95.8*C_Qq38*C_Qq38 + 119.4*C_tq8 + 42.1*C_tq8*C_tq8 - 31.9*C_Qq18*C_Qq38 + 7.1*C_Qq18*C_tq8 - 3.4*C_Qq38*C_tq8  )*(SM_ttW_inc/LOxttw_madgraph);
            }
            else{
                return  (SM_ttW_inc + (158.3*C_tG)*(SM_ttW_inc/NLOxttw_madgraph))+
                        (119.3*C_Qq18 - 43.7*C_Qq38 + 119.4*C_tq8 )*(SM_ttW_inc/LOxttw_madgraph);
            }      
    }
}




Asymmetry_leptonic_charge_rapidity_ttW::Asymmetry_leptonic_charge_rapidity_ttW(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_Asymmetry_leptonic_charge_rapidity_ttW"
            << "SM_ttW_inc");
    
}

double Asymmetry_leptonic_charge_rapidity_ttW::computeThValue()
{
    double SM_Asymmetry_leptonic_charge_rapidity_ttW = SM.getOptionalParameter("SM_Asymmetry_leptonic_charge_rapidity_ttW");
    double SM_ttW_inc = SM.getOptionalParameter("SM_ttW_inc");
    double SM_numerator = SM_Asymmetry_leptonic_charge_rapidity_ttW*SM_ttW_inc;
    double LOxttw_madgraph = 407.2;//fb
    double numerator_madgraph =-74.59904 ;//fb
    //double Asymmetry_leptonic_charge_rapidity_ttW_madgraph=numerator_madgraph/LOxttw_madgraph;
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    //double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double xsection;
    double numerator;

    //We remove C_tG from the xsection and from the numerator. This must be checked later though!!!!!

    xsection= SM_ttW_inc + (119.3*C_Qq18 + 42.0*C_Qq18*C_Qq18 - 43.7*C_Qq38 + 95.8*C_Qq38*C_Qq38 + 119.4*C_tq8 + 42.1*C_tq8*C_tq8 - 31.9*C_Qq18*C_Qq38 + 7.1*C_Qq18*C_tq8 - 3.4*C_Qq38*C_tq8 )*(SM_ttW_inc/LOxttw_madgraph);
    
    numerator= SM_numerator + (9.7*C_Qq18 + 12.2*C_Qq18*C_Qq18 - 7.0*C_Qq38 + 23.5*C_Qq38*C_Qq38 - 40.4*C_tq8 - 16.4*C_tq8*C_tq8 - 12.8*C_Qq18*C_Qq38 - 1.6*C_Qq18*C_tq8 + 0.5*C_Qq38*C_tq8 )*(SM_numerator/numerator_madgraph);


    //std::cout << "SM_Asymmetry_leptonic_charge_rapidity_ttW = " << SM_Asymmetry_leptonic_charge_rapidity_ttW << std::endl;
    //std::cout << " "   << std::endl;
    //std::cout << "SM_numerator = " << SM_numerator << std::endl;
    //std::cout << " "   << std::endl;
    
    //std::cout << "SM_ttW_inc = " << SM_ttW_inc << std::endl;
    //std::cout << "xsection = " << xsection << std::endl;
    //std::cout << " "   << std::endl;
    //std::cout << "SM_numerator = " << SM_numerator << std::endl;
    //std::cout << "numerator = " << numerator << std::endl;
    //std::cout << " "   << std::endl;
    
    //std::cout << "numerator/xsection = " << numerator/xsection << std::endl;
    
    
    
    return numerator/xsection;
            
}




sigmatchannel13::sigmatchannel13(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigmatchannel13" );
}

double sigmatchannel13::computeThValue()
{
    double SM_sigmatchannel13 = SM.getOptionalParameter("SM_sigmatchannel13");
    double NLOxtq_madgraph_top = 211233;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_phitb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_sigmatchannel13 + (6.24167*C_tW + 25.4583*C_phiQ3 -0.9*C_tG + 3.24583*C_tW*C_tW 
                        + 0.845833*C_phiQ3*C_phiQ3 + 2.2*C_tG*C_tG+ 1.4*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532)) + 0.17*(-C_phitb/(0.998))*(-C_phitb/(0.998))+ 0.65625*C_tW*C_phiQ3+ 0.12*(-C_phitb/(0.998))*(C_bW/(0.999*0.6532)) )*(1000*SM_sigmatchannel13/NLOxtq_madgraph_top) ;
            }
            else{
                return  SM_sigmatchannel13 + (6.24167*C_tW + 25.4583*C_phiQ3 -0.9*C_tG )*(1000*SM_sigmatchannel13/NLOxtq_madgraph_top) ;     
            }      
    }
    else{
            if(flag_Quadratic){
                return  SM_sigmatchannel13 + (6.24167*C_tW + 25.4583*C_phiQ3 -0.9*C_tG + 3.24583*C_tW*C_tW 
                        + 0.845833*C_phiQ3*C_phiQ3 + 2.2*C_tG*C_tG+ 1.4*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532)) + 0.17*(-C_phitb/(0.998))*(-C_phitb/(0.998))+ 0.65625*C_tW*C_phiQ3+ 0.12*(-C_phitb/(0.998))*(C_bW/(0.999*0.6532)) )*(1000*SM_sigmatchannel13/NLOxtq_madgraph_top) ;
            }
            else{
                return  SM_sigmatchannel13 + (6.24167*C_tW + 25.4583*C_phiQ3 -0.9*C_tG )*(1000*SM_sigmatchannel13/NLOxtq_madgraph_top) ;     
            }      
    }
    
}


sigmatchannel8::sigmatchannel8(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigmatchannel8" );
}

double sigmatchannel8::computeThValue()
    {
    double  SM_sigmatchannel8 = SM.getOptionalParameter("SM_sigmatchannel8");
    double NLOxtq_madgraph_top = 82600;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_phitb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_sigmatchannel8 + (2.465*C_tW+9.92667*C_phiQ3 +1.27*C_tW*C_tW + 0.2275*C_phiQ3*C_phiQ3 + 0.46*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532)) + 0.075*(-C_phitb/(0.998))*(-C_phitb/(0.998)) + 0.163125*C_tW*C_phiQ3 +0.022*(C_bW/(0.999*0.6532))*(-C_phitb/(0.998)) )*(1000*SM_sigmatchannel8/NLOxtq_madgraph_top) ;      
            }
            else{
                    return  SM_sigmatchannel8 + (2.465*C_tW+9.92667*C_phiQ3 )*(1000*SM_sigmatchannel8/NLOxtq_madgraph_top) ;      
            }      
    }
    
    else{
            if(flag_Quadratic){
                return  SM_sigmatchannel8 + (2.465*C_tW+9.92667*C_phiQ3 +1.27*C_tW*C_tW + 0.2275*C_phiQ3*C_phiQ3 + 0.46*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532)) + 0.075*(-C_phitb/(0.998))*(-C_phitb/(0.998)) + 0.163125*C_tW*C_phiQ3 +0.022*(C_bW/(0.999*0.6532))*(-C_phitb/(0.998)) )*(1000*SM_sigmatchannel8/NLOxtq_madgraph_top) ;      
            }
            else{
                    return  SM_sigmatchannel8 + (2.465*C_tW+9.92667*C_phiQ3 )*(1000*SM_sigmatchannel8/NLOxtq_madgraph_top) ;      
            }      
    }
    
}


sigmaschannel8::sigmaschannel8(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigmaschannel8" );
}

double sigmaschannel8::computeThValue()
{
    double  SM_sigmaschannel8 = SM.getOptionalParameter("SM_sigmaschannel8");
    double NLOxtq_madgraph_top = 4843;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_phitb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_sigmaschannel8 + (2276*C_tW + 593.5*C_phiQ3 -4.0*C_tG+ 421.593*C_tW*C_tW 
                        + 21.75*C_phiQ3*C_phiQ3 + 88.0*C_tG*C_tG+ 178.8*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532)) + 4.2*(-C_phitb/(0.998))*(-C_phitb/(0.998)) 
                        + 140.375*C_tW*C_phiQ3 )*(SM_sigmaschannel8/NLOxtq_madgraph_top) ;      
            }
            else{
                return SM_sigmaschannel8 + (2276*C_tW + 593.5*C_phiQ3 -4.0*C_tG)*(SM_sigmaschannel8/NLOxtq_madgraph_top) ;     
            }      
    }
    
    else{
            if(flag_Quadratic){
                return  SM_sigmaschannel8 + (2276*C_tW + 593.5*C_phiQ3 -4.0*C_tG+ 421.593*C_tW*C_tW 
                        + 21.75*C_phiQ3*C_phiQ3 + 88.0*C_tG*C_tG+ 178.8*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532)) + 4.2*(-C_phitb/(0.998))*(-C_phitb/(0.998)) 
                        + 140.375*C_tW*C_phiQ3 )*(SM_sigmaschannel8/NLOxtq_madgraph_top) ;      
            }
            else{
                return SM_sigmaschannel8 + (2276*C_tW + 593.5*C_phiQ3 -4.0*C_tG)*(SM_sigmaschannel8/NLOxtq_madgraph_top) ;     
            }      
    }
    
}


sigmatW::sigmatW(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_tW_inc" );
}

double sigmatW::computeThValue()
{
    double SM_tW_inc = SM.getOptionalParameter("SM_tW_inc");
    double NLOxgbtw_madgraph = 56188.5;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_tW_inc + (-4.84212*C_tW + 6.79739*C_phiQ3+ 4.2*C_tG+ 1.3335*C_tW*C_tW + 0.207038*C_phiQ3*C_phiQ3 + 2*C_tG*C_tG -0.292105*C_tW*C_phiQ3)*(1000*SM_tW_inc/NLOxgbtw_madgraph) ;
            }
            else{
                return  SM_tW_inc + (-4.84212*C_tW + 6.79739*C_phiQ3+ 4.2*C_tG)*(1000*SM_tW_inc/NLOxgbtw_madgraph) ;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  SM_tW_inc + (-4.84212*C_tW + 6.79739*C_phiQ3+ 4.2*C_tG+ 1.3335*C_tW*C_tW + 0.207038*C_phiQ3*C_phiQ3 + 2*C_tG*C_tG -0.292105*C_tW*C_phiQ3)*(1000*SM_tW_inc/NLOxgbtw_madgraph) ;
            }
            else{
                return  SM_tW_inc + (-4.84212*C_tW + 6.79739*C_phiQ3+ 4.2*C_tG)*(1000*SM_tW_inc/NLOxgbtw_madgraph) ;
            }
    }
}



sigmatW_8TeV::sigmatW_8TeV(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_tW_inc_8TeV" );
}

double sigmatW_8TeV::computeThValue()
{
    double  SM_tW_inc_8TeV = SM.getOptionalParameter("SM_tW_inc_8TeV");
    double NLOxgbtw_madgraph_8TeV = 17.3969;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_phitb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_tW_inc_8TeV + (-1.98*C_tW + 0.49*C_tW*C_tW + 2.72*C_phiQ3 + 0.06*C_phiQ3*C_phiQ3
                        + 1.30178*C_tG + 0.129086*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532)) + 0.0143561*(-C_phitb/(0.998))*(-C_phitb/(0.998)) -0.0175843*(C_bW/(0.999*0.6532))*(-C_phitb/(0.998)))*(SM_tW_inc_8TeV/NLOxgbtw_madgraph_8TeV) ;
            }
            else{
                return  SM_tW_inc_8TeV + (-1.98*C_tW  + 2.72*C_phiQ3 
                        + 1.30178*C_tG)*(SM_tW_inc_8TeV/NLOxgbtw_madgraph_8TeV);
            }
    }
    
    else{
            if(flag_Quadratic){
                return  SM_tW_inc_8TeV + (-1.98*C_tW + 0.49*C_tW*C_tW + 2.72*C_phiQ3 + 0.06*C_phiQ3*C_phiQ3
                        + 1.30178*C_tG + 0.129086*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532)) + 0.0143561*(-C_phitb/(0.998))*(-C_phitb/(0.998)) -0.0175843*(C_bW/(0.999*0.6532))*(-C_phitb/(0.998)))*(SM_tW_inc_8TeV/NLOxgbtw_madgraph_8TeV) ;
            }
            else{
                return  SM_tW_inc_8TeV + (-1.98*C_tW  + 2.72*C_phiQ3 
                        + 1.30178*C_tG)*(SM_tW_inc_8TeV/NLOxgbtw_madgraph_8TeV);
            }
    }
    
}



sigmatqZ::sigmatqZ(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_tZQ_inc" );
}

double sigmatqZ::computeThValue()
{
    double  SM_tZQ_inc = SM.getOptionalParameter("SM_tZQ_inc");
    double xztq_NLO_madgraph = 806.99;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_phitb =  myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return SM_tZQ_inc +(-5.9*C_tZ + 22.2*C_tW + 176.075*C_phiQ3 + 
                        17.9006*C_phiQm + 4.66561*C_phit - 3.9*C_tG
                        + 15.3*C_tZ*C_tZ + 80.0*C_tW*C_tW + 25.439*C_phiQ3*C_phiQ3
                        + 1.23048*C_phiQm*C_phiQm + 0.492986*C_phit*C_phit + 5.4*C_tG*C_tG
                        + 48.7*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532))+3.4*(-C_phitb/(0.998))*(-C_phitb/(0.998))
                        + 4.19554*C_phiQ3*C_phiQm -0.531018*C_phiQ3*C_phit -0.818108*C_phiQm*C_phit-1.1*(C_bW/(0.999*0.6532))*(-C_phitb/(0.998))
                        )*(SM_tZQ_inc/xztq_NLO_madgraph);
            }
            else{
                return SM_tZQ_inc +(-5.9*C_tZ + 22.2*C_tW + 176.075*C_phiQ3 + 
                        17.9006*C_phiQm + 4.66561*C_phit - 3.9*C_tG)*(SM_tZQ_inc/xztq_NLO_madgraph);
            }
    }
    
    else{
            if(flag_Quadratic){
                return SM_tZQ_inc +(-5.9*(-0.472123*C_tB + 0.881533*C_tW) + 22.2*C_tW + 176.075*C_phiQ3 + 
                        17.9006*(C_phiQ1-C_phiQ3) + 4.66561*C_phit - 3.9*C_tG
                        + 15.3*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 80.0*C_tW*C_tW + 25.439*C_phiQ3*C_phiQ3
                        + 1.23048*(C_phiQ1-C_phiQ3)*(C_phiQ1-C_phiQ3) + 0.492986*C_phit*C_phit + 5.4*C_tG*C_tG
                        + 48.7*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532))+3.4*(-C_phitb/(0.998))*(-C_phitb/(0.998))
                        + 4.19554*C_phiQ3*(C_phiQ1-C_phiQ3) -0.531018*C_phiQ3*C_phit -0.818108*(C_phiQ1-C_phiQ3)*C_phit-1.1*(C_bW/(0.999*0.6532))*(-C_phitb/(0.998))
                        )*(SM_tZQ_inc/xztq_NLO_madgraph);
            }
            else{
                return SM_tZQ_inc +(-5.9*(-0.472123*C_tB + 0.881533*C_tW) + 22.2*C_tW + 176.075*C_phiQ3 + 
                        17.9006*(C_phiQ1-C_phiQ3) + 4.66561*C_phit - 3.9*C_tG)*(SM_tZQ_inc/xztq_NLO_madgraph);
            }
    }
}



sigmatAq::sigmatAq(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_tAq_inc" );
}

double sigmatAq::computeThValue()
{
    double  SM_tAq_inc = SM.getOptionalParameter("SM_tAq_inc");
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_phitb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double xtAqNLO_madgraph = 869.5;

    
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_tAq_inc + (-2.1957*C_tZ + 45.587*C_tZ*C_tZ + 22.9802*C_tW + 31.2686*C_tW*C_tW 
                        + 105.377*C_phiQ3  
                        + 0.792383*C_phiQ3*C_phiQ3
                        + 4.6*C_tG -6.4*C_tG*C_tG -1.18213*C_phiQm*C_phiQ3 + 0.0499714*C_phiQ3*C_phit 
                        + 0.215421*C_phiQm*C_phit -73.4469*C_tZ*C_tW
                        +13.1*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532))
                        - 0.8*(-C_phitb/(0.998))*(-C_phitb/(0.998)))*(SM_tAq_inc/(xtAqNLO_madgraph));
            }
            else{
                return  SM_tAq_inc + (-2.1957*C_tZ  + 22.9802*C_tW  + 0.539446*C_phiQm  
                        + 105.377*C_phiQ3   -0.197342*C_phit
                        )*(SM_tAq_inc/(xtAqNLO_madgraph));
            }
    }
    
    else{
            if(flag_Quadratic){
                return  SM_tAq_inc + (-2.1957*(-0.472123*C_tB + 0.881533*C_tW) + 45.587*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 22.9802*C_tW + 31.2686*C_tW*C_tW 
                        + 105.377*C_phiQ3  
                        + 0.792383*C_phiQ3*C_phiQ3
                        + 4.6*C_tG -6.4*C_tG*C_tG -1.18213*(C_phiQ1-C_phiQ3)*C_phiQ3 + 0.0499714*C_phiQ3*C_phit 
                        + 0.215421*(C_phiQ1-C_phiQ3)*C_phit -73.4469*(-0.472123*C_tB + 0.881533*C_tW)*C_tW
                        +13.1*(C_bW/(0.999*0.6532))*(C_bW/(0.999*0.6532))
                        - 0.8*(-C_phitb/(0.998))*(-C_phitb/(0.998)))*(SM_tAq_inc/(xtAqNLO_madgraph));
            }
            else{
                return  SM_tAq_inc + (-2.1957*(-0.472123*C_tB + 0.881533*C_tW)  + 22.9802*C_tW  + 0.539446*(C_phiQ1-C_phiQ3)  
                        + 105.377*C_phiQ3   -0.197342*C_phit
                        )*(SM_tAq_inc/(xtAqNLO_madgraph));
            }
    }
    
}

tH_tchan::tH_tchan(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_tH_tchan_value" );
}

double tH_tchan::computeThValue()
{
    double  SM_tH_tchan_value = SM.getOptionalParameter("SM_tH_tchan_value");
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double NLO_madgraph = 71.18;//fb
    double upper_lim = 8*SM_tH_tchan_value;
    
    //We introduce the ratio Th/Exp to take into account the upper limit
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                double  cross_sect =  SM_tH_tchan_value + ((21.1*C_tW-8.3*C_phiQ3-0.1*C_tG-0.5*C_tphi+56.1*C_tW*C_tW+13.4*C_phiQ3*C_phiQ3-1.3*C_tG*C_tG+0.8*C_tphi*C_tphi))*(SM_tH_tchan_value/NLO_madgraph);
                return  cross_sect/upper_lim;
            }
            else{
                double  cross_sect = SM_tH_tchan_value + ((21.1*C_tW-8.3*C_phiQ3-0.1*C_tG-0.5*C_tphi))*(SM_tH_tchan_value/NLO_madgraph);
                return  cross_sect/upper_lim;
            }
    }
    
    else{
            if(flag_Quadratic){
                double  cross_sect =  SM_tH_tchan_value + ((21.1*C_tW-8.3*C_phiQ3-0.1*C_tG-0.5*C_tphi+56.1*C_tW*C_tW+13.4*C_phiQ3*C_phiQ3-1.3*C_tG*C_tG+0.8*C_tphi*C_tphi))*(SM_tH_tchan_value/NLO_madgraph);
                return  cross_sect/upper_lim;
            }
            else{
                double  cross_sect = SM_tH_tchan_value + ((21.1*C_tW-8.3*C_phiQ3-0.1*C_tG-0.5*C_tphi))*(SM_tH_tchan_value/NLO_madgraph);
                return  cross_sect/upper_lim;
            }
    }
}




sigmattH::sigmattH(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_ttH_inc" );
}

double sigmattH::computeThValue()
{
    
    
    double  SM_ttH_inc = SM.getOptionalParameter("SM_ttH_inc");
    double xtth_madgraph_NLO = 459.7;//fb
//    double xtth_madgraph_LO = 354.57;//fb //with dinamical scale
    double xtth_madgraph_LO = 459;//fb //with fixed scale
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    //double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    //double C_phitb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    double mC_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double mC_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttH_inc + (-59.3*C_tphi+ 458.9*C_tG +1.9*C_tphi*C_tphi+ 598.6*C_tG*C_tG)*(SM_ttH_inc/xtth_madgraph_NLO)+
                        (0.0203*mC_tu8+0.0126*mC_td8+0.0088*mC_tu8*mC_tu8+0.0046*mC_td8*mC_td8 + 0.038028*C_Qq18 + 0.03802*C_tq8 + 0.009163*C_Qq38 + 0.023543*C_Qu8 + 0.014638*C_Qd8)*(1000*SM_ttH_inc/xtth_madgraph_LO) ;
            }
            else{
                return  SM_ttH_inc + (-59.3*C_tphi+ 458.9*C_tG )*(SM_ttH_inc/xtth_madgraph_NLO)+
                        (0.0203*mC_tu8+0.0126*mC_td8 + 0.038028*C_Qq18 + 0.03802*C_tq8 + 0.009163*C_Qq38 + 0.023543*C_Qu8 + 0.014638*C_Qd8)*(1000*SM_ttH_inc/xtth_madgraph_LO);
            }
    }
    
    else{
            if(flag_Quadratic){
                return  SM_ttH_inc + (-59.3*C_tphi+ 458.9*C_tG +1.9*C_tphi*C_tphi+ 598.6*C_tG*C_tG)*(SM_ttH_inc/xtth_madgraph_NLO)+
                        (0.0203*mC_tu8+0.0126*mC_td8+0.0088*mC_tu8*mC_tu8+0.0046*mC_td8*mC_td8)*(1000*SM_ttH_inc/xtth_madgraph_LO) ;
            }
            else{
                return  SM_ttH_inc + (-59.3*C_tphi+ 458.9*C_tG )*(SM_ttH_inc/xtth_madgraph_NLO)+
                        (0.0203*mC_tu8+0.0126*mC_td8)*(1000*SM_ttH_inc/xtth_madgraph_LO);
            }
    }
}


ttHSUM::ttHSUM(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_ttH_inc" << "SM_tH_tchan_value" );
}

double ttHSUM::computeThValue()
{
    
    double SM_ttH_inc = SM.getOptionalParameter("SM_ttH_inc");
    double SM_tH_tchan_value = SM.getOptionalParameter("SM_tH_tchan_value");
    
    double xtth_madgraph_NLO = 459.7;//fb
//    double xtth_madgraph_LO = 354.57;//fb //with dinamical scale
    double xtth_madgraph_LO = 459;//fb //with fixed scale
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_phitb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    double mC_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double mC_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    


    
    double xtH_madgraph_NLO = 71.18;//fb
    double xtH_madgraph_LO = 66.6;//fb

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                double  x_sec_ttH = SM_ttH_inc + (-59.3*C_tphi+ 458.9*C_tG +1.9*C_tphi*C_tphi+ 598.6*C_tG*C_tG)*(SM_ttH_inc/xtth_madgraph_NLO)+
                        (0.0203*mC_tu8+0.0126*mC_td8+0.0088*mC_tu8*mC_tu8+0.0046*mC_td8*mC_td8 + 0.038028*C_Qq18 + 0.03802*C_tq8 + 0.009163*C_Qq38 + 0.023543*C_Qu8 + 0.014638*C_Qd8)*(1000*SM_ttH_inc/xtth_madgraph_LO) ;

                double  x_sec_tH =  SM_tH_tchan_value + ((21.1*C_tW-8.3*C_phiQ3-0.1*C_tG-0.5*C_tphi+56.1*C_tW*C_tW+13.4*C_phiQ3*C_phiQ3-1.3*C_tG*C_tG+0.8*C_tphi*C_tphi))*(SM_tH_tchan_value/xtH_madgraph_NLO)+
                (1.9*C_bW*C_bW+0.06*C_phitb*C_phitb)*(SM_tH_tchan_value/xtH_madgraph_LO);
                
                return  x_sec_ttH+x_sec_tH;
            }
            else{
                double  x_sec_tH = SM_tH_tchan_value + ((21.1*C_tW-8.3*C_phiQ3-0.1*C_tG-0.5*C_tphi))*(SM_tH_tchan_value/xtH_madgraph_NLO);
                
                double  x_sec_ttH = SM_ttH_inc + (-59.3*C_tphi+ 458.9*C_tG )*(SM_ttH_inc/xtth_madgraph_NLO)+
                        (0.0203*mC_tu8+0.0126*mC_td8 + 0.038028*C_Qq18 + 0.03802*C_tq8 + 0.009163*C_Qq38 + 0.023543*C_Qu8 + 0.014638*C_Qd8)*(1000*SM_ttH_inc/xtth_madgraph_LO);
                
                return  x_sec_tH+x_sec_ttH;
            }
    }
    
    else{
            if(flag_Quadratic){
                
                double  x_sec_tH  =  SM_tH_tchan_value + ((21.1*C_tW-8.3*C_phiQ3-0.1*C_tG-0.5*C_tphi+56.1*C_tW*C_tW+13.4*C_phiQ3*C_phiQ3-1.3*C_tG*C_tG+0.8*C_tphi*C_tphi))*(SM_tH_tchan_value/xtH_madgraph_NLO)+
                        (1.9*C_bW*C_bW+0.06*C_phitb*C_phitb)*(SM_tH_tchan_value/xtH_madgraph_LO);
                
                double  x_sec_ttH  =  SM_ttH_inc + (-59.3*C_tphi+ 458.9*C_tG +1.9*C_tphi*C_tphi+ 598.6*C_tG*C_tG)*(SM_ttH_inc/xtth_madgraph_NLO)+
                        (0.0203*mC_tu8+0.0126*mC_td8+0.0088*mC_tu8*mC_tu8+0.0046*mC_td8*mC_td8)*(1000*SM_ttH_inc/xtth_madgraph_LO) ;
                
                return  x_sec_tH+x_sec_ttH;
            }
            else{
                double  x_sec_tH  = SM_tH_tchan_value + ((21.1*C_tW-8.3*C_phiQ3-0.1*C_tG-0.5*C_tphi))*(SM_tH_tchan_value/xtH_madgraph_NLO);
                
                double  x_sec_ttH  = SM_ttH_inc + (-59.3*C_tphi+ 458.9*C_tG )*(SM_ttH_inc/xtth_madgraph_NLO)+
                        (0.0203*mC_tu8+0.0126*mC_td8)*(1000*SM_ttH_inc/xtth_madgraph_LO);

                return  x_sec_tH+x_sec_ttH;
            }
    }
}





sigma_ttH_diff_NLO_ATLAS_220700092::sigma_ttH_diff_NLO_ATLAS_220700092(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_ttH_bin_0_60_ATLAS_220700092" << "SM_sigma_ttH_bin_60_120_ATLAS_220700092"
            << "SM_sigma_ttH_bin_120_200_ATLAS_220700092" << "SM_sigma_ttH_bin_200_300_ATLAS_220700092" << "SM_sigma_ttH_bin_300_450_ATLAS_220700092" 
            << "SM_sigma_ttH_bin_450_inf_ATLAS_220700092");

    
}

double sigma_ttH_diff_NLO_ATLAS_220700092::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    //double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    //double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    //double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    //double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    
   
    if(b_min == 0 && b_max == 60){
        
        double SM_ttH_bin_0_60 = SM.getOptionalParameter("SM_sigma_ttH_bin_0_60_ATLAS_220700092");
        double ttH_bin_0_60_madgraph_NLO = 0.122104;//pb
 
        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttH_bin_0_60 + (0.003954*C_tq8+0.000476*C_tq8*C_tq8+-0.000291*C_tu1+0.002843*C_tu1*C_tu1+
                        -0.000041*C_Qu1+0.002855*C_Qu1*C_Qu1+0.002340*C_tu8+0.000365*C_tu8*C_tu8+0.000001*C_Qq31+
                        0.004657*C_Qq31*C_Qq31+0.000787*C_Qq38+0.000650*C_Qq38*C_Qq38+0.001542*C_td8+
                        0.000155*C_td8*C_td8+0.112300*C_tG+0.071001*C_tG*C_tG+-0.015672*C_tphi+0.000377*C_tphi*C_tphi+
                        -0.000194*C_td1+0.001737*C_td1*C_td1+0.001686*C_Qd8+0.000047*C_Qd8*C_Qd8+0.000146*C_tq1+
                        0.004604*C_tq1*C_tq1+0.000147*C_Qd1+0.001638*C_Qd1*C_Qd1+0.003818*C_Qq18+0.000672*C_Qq18*C_Qq18+
                        -0.000349*C_Qq11+0.004600*C_Qq11*C_Qq11+0.002331*C_Qu8+0.000259*C_Qu8*C_Qu8)*(SM_ttH_bin_0_60/ttH_bin_0_60_madgraph_NLO);

            }
            else{
                return  SM_ttH_bin_0_60 + (0.003954*C_tq8+-0.000291*C_tu1+-0.000041*C_Qu1+0.002340*C_tu8+
                        0.000001*C_Qq31+0.000787*C_Qq38+0.001542*C_td8+0.112300*C_tG+-0.015672*C_tphi+
                        -0.000194*C_td1+0.001686*C_Qd8+0.000146*C_tq1+0.000147*C_Qd1+0.003818*C_Qq18+
                        -0.000349*C_Qq11+0.002331*C_Qu8)*(SM_ttH_bin_0_60/ttH_bin_0_60_madgraph_NLO);
            }
        }
        else{
            throw std::runtime_error("sigma_ttH_diff_NLO_ATLAS_231204450 is only available in the LHC top-quark WC basis");
            
        }
        
    } else if(b_min == 60 && b_max == 120){
        
        double SM_ttH_bin_60_120 = SM.getOptionalParameter("SM_sigma_ttH_bin_60_120_ATLAS_220700092");
        
        double ttH_bin_60_120_madgraph_NLO = 0.182098;//pb
        
 

        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttH_bin_60_120  + (0.006999*C_tq8+0.001290*C_tq8*C_tq8+-0.000220*C_tu1+0.006516*C_tu1*C_tu1+
                        -0.000219*C_Qu1+0.006475*C_Qu1*C_Qu1+0.004522*C_tu8+0.000887*C_tu8*C_tu8+0.000060*C_Qq31+
                        0.010350*C_Qq31*C_Qq31+0.001469*C_Qq38+0.001545*C_Qq38*C_Qq38+0.002741*C_td8+
                        0.000622*C_td8*C_td8+0.179050*C_tG+0.142520*C_tG*C_tG+-0.023288*C_tphi+0.000643*C_tphi*C_tphi+
                        -0.000075*C_td1+0.003796*C_td1*C_td1+0.003172*C_Qd8+0.000454*C_Qd8*C_Qd8+-0.000406*C_tq1+
                        0.010404*C_tq1*C_tq1+0.000287*C_Qd1+0.003741*C_Qd1*C_Qd1+0.007293*C_Qq18+0.001702*C_Qq18*C_Qq18+
                        -0.000690*C_Qq11+0.010380*C_Qq11*C_Qq11+0.004294*C_Qu8+0.000740*C_Qu8*C_Qu8)*(SM_ttH_bin_60_120/ttH_bin_60_120_madgraph_NLO);
            }
            else{
                return SM_ttH_bin_60_120  + (0.006999*C_tq8+-0.000220*C_tu1+-0.000219*C_Qu1+0.004522*C_tu8+
                        0.000060*C_Qq31+0.001469*C_Qq38+0.002741*C_td8+0.179050*C_tG+-0.023288*C_tphi+
                        -0.000075*C_td1+0.003172*C_Qd8+-0.000406*C_tq1+0.000287*C_Qd1+0.007293*C_Qq18+
                        -0.000690*C_Qq11+0.004294*C_Qu8)*(SM_ttH_bin_60_120/ttH_bin_60_120_madgraph_NLO);
            }
        }
    
        else{
            throw std::runtime_error("sigma_ttH_diff_NLO_ATLAS_231204450 is only available in the LHC top-quark WC basis");
        }
        
    } else if(b_min == 120 && b_max == 200){
        
        double SM_ttH_bin_120_200 = SM.getOptionalParameter("SM_sigma_ttH_bin_120_200_ATLAS_220700092");        
        
        double ttH_bin_120_200_madgraph_NLO = 0.127644 ;//pb
      
 

        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttH_bin_120_200 + (0.007289*C_tq8+0.001984*C_tq8*C_tq8+0.000077*C_tu1+0.009008*C_tu1*C_tu1+
                        -0.000507*C_Qu1+0.009023*C_Qu1*C_Qu1+0.004404*C_tu8+0.001658*C_tu8*C_tu8+0.000139*C_Qq31+
                        0.014176*C_Qq31*C_Qq31+0.001568*C_Qq38+0.002337*C_Qq38*C_Qq38+0.002972*C_td8+
                        0.001020*C_td8*C_td8+0.137374*C_tG+0.169359*C_tG*C_tG+-0.016649*C_tphi+0.000766*C_tphi*C_tphi+
                        -0.000406*C_td1+0.005443*C_td1*C_td1+0.003068*C_Qd8+0.000855*C_Qd8*C_Qd8+-0.000509*C_tq1+
                        0.014199*C_tq1*C_tq1+-0.000002*C_Qd1+0.005280*C_Qd1*C_Qd1+0.007130*C_Qq18+0.002483*C_Qq18*C_Qq18+
                        -0.000146*C_Qq11+0.014233*C_Qq11*C_Qq11+0.004402*C_Qu8+0.001360*C_Qu8*C_Qu8)*(SM_ttH_bin_120_200/ttH_bin_120_200_madgraph_NLO) ;
            }
            else{
                return  SM_ttH_bin_120_200 + (0.007289*C_tq8+0.000077*C_tu1+-0.000507*C_Qu1+0.004404*C_tu8+
                        0.000139*C_Qq31+0.001568*C_Qq38+0.002972*C_td8+0.137374*C_tG+-0.016649*C_tphi+
                        -0.000406*C_td1+0.003068*C_Qd8+-0.000509*C_tq1+-0.000002*C_Qd1+0.007130*C_Qq18+
                        -0.000146*C_Qq11+0.004402*C_Qu8)*(SM_ttH_bin_120_200/ttH_bin_120_200_madgraph_NLO);
            }
        }
    
        else{
            throw std::runtime_error("sigma_ttH_diff_NLO_ATLAS_231204450 is only available in the LHC top-quark WC basis");
        }
        
    } else if(b_min == 200 && b_max == 300){
        
        double SM_ttH_bin_200_300 = SM.getOptionalParameter("SM_sigma_ttH_bin_200_300_ATLAS_220700092");
        
        double ttH_bin_200_300_madgraph_NLO = 0.053671;//pb
 


        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttH_bin_200_300 + (0.004487*C_tq8+0.001602*C_tq8*C_tq8+-0.000015*C_tu1+0.008477*C_tu1*C_tu1+
                        -0.000260*C_Qu1+0.008442*C_Qu1*C_Qu1+0.002946*C_tu8+0.001315*C_tu8*C_tu8+0.000075*C_Qq31+
                        0.013171*C_Qq31*C_Qq31+0.001182*C_Qq38+0.001998*C_Qq38*C_Qq38+0.001907*C_td8+
                        0.000713*C_td8*C_td8+0.057738*C_tG+0.134324*C_tG*C_tG+-0.006926*C_tphi+0.000175*C_tphi*C_tphi+
                        -0.000128*C_td1+0.004725*C_td1*C_td1+0.001854*C_Qd8+0.000559*C_Qd8*C_Qd8+-0.000546*C_tq1+
                        0.013203*C_tq1*C_tq1+-0.000176*C_Qd1+0.004707*C_Qd1*C_Qd1+0.004815*C_Qq18+0.002123*C_Qq18*C_Qq18+
                        0.000140*C_Qq11+0.013192*C_Qq11*C_Qq11+0.002754*C_Qu8+0.001011*C_Qu8*C_Qu8)*(SM_ttH_bin_200_300/ttH_bin_200_300_madgraph_NLO);
            }
            else{
                return  SM_ttH_bin_200_300 + (0.004487*C_tq8+-0.000015*C_tu1+-0.000260*C_Qu1+0.002946*C_tu8+
                        0.000075*C_Qq31+0.001182*C_Qq38+0.001907*C_td8+0.057738*C_tG+-0.006926*C_tphi+
                        -0.000128*C_td1+0.001854*C_Qd8+-0.000546*C_tq1+-0.000176*C_Qd1+0.004815*C_Qq18+
                        0.000140*C_Qq11+0.002754*C_Qu8)*(SM_ttH_bin_200_300/ttH_bin_200_300_madgraph_NLO);
            }
        }
   
        else{
            throw std::runtime_error("sigma_ttH_diff_NLO_ATLAS_231204450 is only available in the LHC top-quark WC basis");
        }
        
    } else if(b_min == 300 && b_max == 450){
        
        double SM_ttH_bin_300_450 = SM.getOptionalParameter("SM_sigma_ttH_bin_300_450_ATLAS_220700092");
        
        double ttH_bin_300_450_madgraph_NLO = 0.019497;//pb

        
        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttH_bin_300_450  + (0.002479*C_tq8+0.001491*C_tq8*C_tq8+0.000105*C_tu1+0.008123*C_tu1*C_tu1+
                        -0.000323*C_Qu1+0.008109*C_Qu1*C_Qu1+0.001760*C_tu8+0.001309*C_tu8*C_tu8+0.000102*C_Qq31+
                        0.012322*C_Qq31*C_Qq31+0.000914331*C_Qq38+0.002613*C_Qq38*C_Qq38+0.001086*C_td8+
                        0.000703*C_td8*C_td8+0.023817*C_tG+0.1089*C_tG*C_tG+-0.002531*C_tphi+0.000083*C_tphi*C_tphi+
                        0.000010*C_td1+0.004281*C_td1*C_td1+0.000971*C_Qd8+0.000525*C_Qd8*C_Qd8+-0.000531*C_tq1+
                        0.012363*C_tq1*C_tq1+-0.000205*C_Qd1+0.004275*C_Qd1*C_Qd1+0.002854*C_Qq18+0.002093*C_Qq18*C_Qq18+
                        0.000142*C_Qq11+0.012312*C_Qq11*C_Qq11+0.001566*C_Qu8+0.000973*C_Qu8*C_Qu8)*(SM_ttH_bin_300_450/ttH_bin_300_450_madgraph_NLO);
            }
            else{
                return  SM_ttH_bin_300_450  + (0.002479*C_tq8+0.000105*C_tu1+-0.000323*C_Qu1+0.001760*C_tu8+
                        0.000102*C_Qq31+0.000914331*C_Qq38+0.001086*C_td8+0.023817*C_tG+-0.002531*C_tphi+
                        0.000010*C_td1+0.000971*C_Qd8+-0.000531*C_tq1+-0.000205*C_Qd1+0.002854*C_Qq18+
                        0.000142*C_Qq11+0.001566*C_Qu8)*(SM_ttH_bin_300_450/ttH_bin_300_450_madgraph_NLO);
            }
        }
    
        else{
            throw std::runtime_error("sigma_ttH_diff_NLO_ATLAS_231204450 is only available in the LHC top-quark WC basis");
        }
        
    } else if(b_min == 450){
        
        double SM_ttH_bin_450_inf = SM.getOptionalParameter("SM_sigma_ttH_bin_450_inf_ATLAS_220700092");
        
        double ttH_bin_450_inf_madgraph_NLO = 0.005472;//pb
 

        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttH_bin_450_inf + (0.001209*C_tq8+0.001536*C_tq8*C_tq8+0.000217*C_tu1+0.010956*C_tu1*C_tu1+
                        -0.000432*C_Qu1+0.010925*C_Qu1*C_Qu1+0.001005*C_tu8+0.001709*C_tu8*C_tu8+-0.000109*C_Qq31+
                        0.015616*C_Qq31*C_Qq31+0.000683838*C_Qq38+0.004062*C_Qq38*C_Qq38+0.000611*C_td8+
                        0.000753*C_td8*C_td8+0.0063*C_tG+0.0913*C_tG*C_tG+-0.000571182*C_tphi+0.000024*C_tphi*C_tphi+
                        0.000124*C_td1+0.004728*C_td1*C_td1+0.000577*C_Qd8+0.000506*C_Qd8*C_Qd8+-0.000386*C_tq1+
                        0.015709*C_tq1*C_tq1+-0.000108*C_Qd1+0.004683*C_Qd1*C_Qd1+0.001527*C_Qq18+0.002670*C_Qq18*C_Qq18+
                        0.000222*C_Qq11+0.015621*C_Qq11*C_Qq11+0.000795*C_Qu8+0.001018*C_Qu8*C_Qu8)*(SM_ttH_bin_450_inf/ttH_bin_450_inf_madgraph_NLO);
            }
            else{
                return  SM_ttH_bin_450_inf + (0.001209*C_tq8+0.000217*C_tu1+-0.000432*C_Qu1+0.001005*C_tu8+
                        -0.000109*C_Qq31+0.000683838*C_Qq38+0.000611*C_td8+0.0063*C_tG+-0.000571182*C_tphi+
                        0.000124*C_td1+0.000577*C_Qd8+-0.000386*C_tq1+-0.000108*C_Qd1+0.001527*C_Qq18+0.000222*C_Qq11+
                        0.000795*C_Qu8)*(SM_ttH_bin_450_inf/ttH_bin_450_inf_madgraph_NLO);
            }
        }
        else{
            throw std::runtime_error("sigma_ttH_diff_NLO_ATLAS_231204450 is only available in the LHC top-quark WC basis");
        }
        
    } else throw std::runtime_error("wrong bin choice in sigma_ttH_diff_NLO_ATLAS_231204450");
   
}


ttWqEM::ttWqEM(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "ttWqEM_SM" );
}

double ttWqEM::computeThValue()
{
    double ttWqEM_SM = SM.getOptionalParameter("ttWqEM_SM");
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double ttWqEM_madgraph = 47.53;
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  ttWqEM_SM + ((-0.4*C_tZ + 1.1*C_tZ*C_tZ + 12.0*C_tW + 27.3*C_tW*C_tW + 3.6*C_phiQ3 
                        + 3.0*C_phiQ3*C_phiQ3 -3.6*C_phiQm + 0.8*C_phiQm*C_phiQm 
                        + 1.8*C_phit + 0.8*C_phit*C_phit + 4.9*C_tG + 1.1*C_tG*C_tG ))*(ttWqEM_SM/ttWqEM_madgraph);
            }
            else{
                return  ttWqEM_SM + ((-0.4*C_tZ + 12.0*C_tW + 3.6*C_phiQ3 -3.6*C_phiQm 
                        + 1.8*C_phit + 4.9*C_tG ))*(ttWqEM_SM/ttWqEM_madgraph);
            }
    }
    else{
            if(flag_Quadratic){
                return  ttWqEM_SM + ((-0.4*(-0.472123*C_tB + 0.881533*C_tW) + 1.1*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 12.0*C_tW + 27.3*C_tW*C_tW + 3.6*C_phiQ3 
                        + 3.0*C_phiQ3*C_phiQ3 -3.6*(C_phiQ1-C_phiQ3) + 0.8*(C_phiQ1-C_phiQ3)*(C_phiQ1-C_phiQ3) 
                        + 1.8*C_phit + 0.8*C_phit*C_phit + 4.9*C_tG + 1.1*C_tG*C_tG ))*(ttWqEM_SM/ttWqEM_madgraph);
            }
            else{
                return  ttWqEM_SM + ((-0.4*(-0.472123*C_tB + 0.881533*C_tW) + 12.0*C_tW + 3.6*C_phiQ3 -3.6*(C_phiQ1-C_phiQ3) 
                        + 1.8*C_phit + 4.9*C_tG  ))*(ttWqEM_SM/ttWqEM_madgraph);
            }
    }
}



ttWqSUM::ttWqSUM(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "ttWqEM_SM" << "SM_ttW_inc");
}

double ttWqSUM::computeThValue()
{
    double ttWqEM_SM = SM.getOptionalParameter("ttWqEM_SM");
    double SM_ttW_inc = SM.getOptionalParameter("SM_ttW_inc");

    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();


    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double ttWqEM_madgraph = 47.53;//fb
    
    
    double NLOxttw_madgraph = 543.3;//fb
    double LOxttw_madgraph = 361.2;//fb

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
            
                double  ttW = (SM_ttW_inc + (158.3*C_tG + 23.4*C_tG*C_tG)*(SM_ttW_inc/NLOxttw_madgraph) + (0.1193*C_Qq18 + 0.1194*C_tq8)*(SM_ttW_inc/LOxttw_madgraph));
  
            
                double  ttWqEM = ttWqEM_SM + ((-0.4*C_tZ + 1.1*C_tZ*C_tZ + 12.0*C_tW + 27.3*C_tW*C_tW + 3.6*C_phiQ3 
                        + 3.0*C_phiQ3*C_phiQ3 -3.6*C_phiQm + 0.8*C_phiQm*C_phiQm 
                        + 1.8*C_phit + 0.8*C_phit*C_phit + 4.9*C_tG + 1.1*C_tG*C_tG ))*(ttWqEM_SM/ttWqEM_madgraph);
                
                return  ttW + ttWqEM;
            }
            else{
                
                double ttW =  (SM_ttW_inc + (158.3*C_tG)*(SM_ttW_inc/NLOxttw_madgraph) + (0.1193*C_Qq18 + 0.1194*C_tq8)*(SM_ttW_inc/LOxttw_madgraph));
            
                double ttWqEM =  ttWqEM_SM + ((-0.4*C_tZ + 12.0*C_tW + 3.6*C_phiQ3 -3.6*C_phiQm 
                        + 1.8*C_phit + 4.9*C_tG  ))*(ttWqEM_SM/ttWqEM_madgraph);
                
                return ttW+ ttWqEM;
            }
    }
    
    else{
            if(flag_Quadratic){
            
                double  ttW = (SM_ttW_inc + (158.3*C_tG + 23.4*C_tG*C_tG)*(SM_ttW_inc/NLOxttw_madgraph));
  
            
                double  ttWqEM = ttWqEM_SM + ((-0.4*(-0.472123*C_tB + 0.881533*C_tW) + 1.1*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 12.0*C_tW + 27.3*C_tW*C_tW + 3.6*C_phiQ3 
                        + 3.0*C_phiQ3*C_phiQ3 -3.6*(C_phiQ1-C_phiQ3) + 0.8*(C_phiQ1-C_phiQ3)*(C_phiQ1-C_phiQ3) 
                        + 1.8*C_phit + 0.8*C_phit*C_phit + 4.9*C_tG + 1.1*C_tG*C_tG ))*(ttWqEM_SM/ttWqEM_madgraph);
                
                return  ttW + ttWqEM;
            }
            
            else{
                double ttW =  (SM_ttW_inc + (158.3*C_tG )*(SM_ttW_inc/NLOxttw_madgraph));
            
                double ttWqEM =  ttWqEM_SM + ((-0.4*(-0.472123*C_tB + 0.881533*C_tW) + 12.0*C_tW + 3.6*C_phiQ3 -3.6*(C_phiQ1-C_phiQ3) 
                        + 1.8*C_phit + 4.9*C_tG  ))*(ttWqEM_SM/ttWqEM_madgraph);
                
                return ttW+ ttWqEM;
            }
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////




/////// ttZ differential cross section for different bins ///////////////////////////////////////////////


sigma_ttz_diff_NLO_ATLAS_210312603::sigma_ttz_diff_NLO_ATLAS_210312603(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_ttz_bin_0_40_ATLAS_210312603" << "SM_sigma_ttz_bin_40_70_ATLAS_210312603"
            << "SM_sigma_ttz_bin_70_110_ATLAS_210312603" << "SM_sigma_ttz_bin_110_160_ATLAS_210312603" << "SM_sigma_ttz_bin_160_220_ATLAS_210312603" 
            << "SM_sigma_ttz_bin_220_290_ATLAS_210312603" << "SM_sigma_ttz_bin_290_400_ATLAS_210312603");

    
}

double sigma_ttz_diff_NLO_ATLAS_210312603::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    
   
    if(b_min == 0 && b_max == 40){
        
        double SM_ttZ_bin_0_40 = SM.getOptionalParameter("SM_sigma_ttz_bin_0_40_ATLAS_210312603");
        double ttZ_bin_0_40_madgraph_NLO = 2.0765;//fb
        double ttZ_bin_0_40_madgraph_LO = 1.67745;//fb
 
        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_0_40 + ((- C_tZ*0.0240625+ 0.66*C_tG +  0.14212675*C_phit -0.13142025*C_phiQm) + (C_tZ*C_tZ*0.0475938+ 0.103*C_tG*C_tG + C_phit*C_phit*0.015851325  + C_phiQm*C_phiQm*0.034422 + C_phit*C_phiQm*0.02448735))*(SM_ttZ_bin_0_40/ttZ_bin_0_40_madgraph_NLO)
                 + (0.0130463*C_tu8 + 0.01893339329*C_td8 + 0.0068088*C_tu8*C_tu8 + 0.00660006299*C_td8*C_td8 + 0.174928*C_Qq18 + 0.151665*C_tq8 -0.072533*C_Qq38 + 0.025343*C_Qu8 + 0.031515*C_Qd8)*(SM_ttZ_bin_0_40/ttZ_bin_0_40_madgraph_LO);

            }
            else{
                return  SM_ttZ_bin_0_40 + ((- C_tZ*0.0240625+ 0.66*C_tG +  0.14212675*C_phit -0.13142025*C_phiQm))*(SM_ttZ_bin_0_40/ttZ_bin_0_40_madgraph_NLO)
                        + (0.0130463*C_tu8 + 0.01893339329*C_td8 + 0.174928*C_Qq18 + 0.151665*C_tq8 -0.072533*C_Qq38 + 0.025343*C_Qu8 + 0.031515*C_Qd8)*(SM_ttZ_bin_0_40/ttZ_bin_0_40_madgraph_LO);
            }
        }
        else{
            if(flag_Quadratic){
                return  SM_ttZ_bin_0_40 + ((- (-0.472123*C_tB + 0.881533*C_tW)*0.0240625+ 0.66*C_tG +  0.14212675*C_phit -0.13142025*(C_phiQ1-C_phiQ3)) + ((-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW)*0.0475938+ 0.103*C_tG*C_tG + C_phit*C_phit*0.015851325  + (C_phiQ1-C_phiQ3)*(C_phiQ1-C_phiQ3)*0.034422 + C_phit*(C_phiQ1-C_phiQ3)*0.02448735))*(SM_ttZ_bin_0_40/ttZ_bin_0_40_madgraph_NLO)
                        + (+ 0.0130463*C_tu8 + 0.01893339329*C_td8 + 0.0068088*C_tu8*C_tu8 + 0.00660006299*C_td8*C_td8)*(SM_ttZ_bin_0_40/ttZ_bin_0_40_madgraph_LO);
            }
            else{
                return  SM_ttZ_bin_0_40 + ((- (-0.472123*C_tB + 0.881533*C_tW)*0.0240625+ 0.66*C_tG +  0.14212675*C_phit -0.13142025*(C_phiQ1-C_phiQ3) ))*(SM_ttZ_bin_0_40/ttZ_bin_0_40_madgraph_NLO)
                        + (0.0130463*C_tu8 + 0.01893339329*C_td8)*(SM_ttZ_bin_0_40/ttZ_bin_0_40_madgraph_LO);
            }
        }
        
    } else if(b_min == 40 && b_max == 70){
        
        double SM_ttZ_bin_40_70 = SM.getOptionalParameter("SM_sigma_ttz_bin_40_70_ATLAS_210312603");
        
        double ttZ_bin_40_70_madgraph_NLO = 4.15667;//fb
        double ttZ_bin_40_70_madgraph_LO =3.42935;//fb
 

        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_40_70  + ((- C_tZ*0.02333333333 + 1.383333333*C_tG + C_phit*0.2299266667 - C_phiQm*0.4160166667)+ (C_tZ*C_tZ*0.1450000+ 0.3666666667*C_tG*C_tG + C_phit*C_phit*0.0113882  + C_phiQm*C_phiQm*0.04426666667 - C_phit*C_phiQm*0.006941766667))*(SM_ttZ_bin_40_70/ttZ_bin_40_70_madgraph_NLO)
                        + ( 0.0304314*C_tu8 + 0.04316951263*C_td8 + 0.0072957*C_tu8*C_tu8 + 0.007787554675*C_td8*C_td8 + 0.336267*C_Qq18 + 0.288633*C_tq8 -0.127963*C_Qq38 + 0.057187*C_Qu8 + 0.064633*C_Qd8)*(SM_ttZ_bin_40_70/ttZ_bin_40_70_madgraph_LO);
            }
            else{
                return SM_ttZ_bin_40_70  + ((- C_tZ*0.02333333333 + 1.383333333*C_tG + C_phit*0.2299266667 - C_phiQm*0.4160166667 ))*(SM_ttZ_bin_40_70/ttZ_bin_40_70_madgraph_NLO)
                        + (0.0304314*C_tu8 + 0.04316951263*C_td8 + 0.336267*C_Qq18 + 0.288633*C_tq8 -0.127963*C_Qq38 + 0.057187*C_Qu8 + 0.064633*C_Qd8)*(SM_ttZ_bin_40_70/ttZ_bin_40_70_madgraph_LO);
            }
        }
    
        else{
            if(flag_Quadratic){
                return  SM_ttZ_bin_40_70  + ((- (-0.472123*C_tB + 0.881533*C_tW)*0.02333333333 + 1.383333333*C_tG + C_phit*0.2299266667 - (C_phiQ1-C_phiQ3)*0.4160166667)+ ((-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW)*0.1450000+ 0.3666666667*C_tG*C_tG + C_phit*C_phit*0.0113882  + (C_phiQ1-C_phiQ3)*(C_phiQ1-C_phiQ3)*0.04426666667 - C_phit*(C_phiQ1-C_phiQ3)*0.006941766667))*(SM_ttZ_bin_40_70/ttZ_bin_40_70_madgraph_NLO)
                        + (0.0304314*C_tu8 + 0.04316951263*C_td8 + 0.0072957*C_tu8*C_tu8 + 0.007787554675*C_td8*C_td8)*(SM_ttZ_bin_40_70/ttZ_bin_40_70_madgraph_LO);
            }
            else{
                return SM_ttZ_bin_40_70  + ((- (-0.472123*C_tB + 0.881533*C_tW)*0.02333333333 + 1.383333333*C_tG + C_phit*0.2299266667 - (C_phiQ1-C_phiQ3)*0.4160166667))*(SM_ttZ_bin_40_70/ttZ_bin_40_70_madgraph_NLO)
                        + (0.0304314*C_tu8 + 0.04316951263*C_td8)*(SM_ttZ_bin_40_70/ttZ_bin_40_70_madgraph_LO);
            }
        }
        
    } else if(b_min == 70 && b_max == 110){
        
        double SM_ttZ_bin_70_110 = SM.getOptionalParameter("SM_sigma_ttz_bin_70_110_ATLAS_210312603");        
        
        double ttZ_bin_70_110_madgraph_NLO = 4.1175;//fb
        double ttZ_bin_70_110_madgraph_LO = 3.43425;//fb
 

        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_70_110 + ((- C_tZ*0.02875+ 1.49*C_tG + C_phit*0.2634875 - C_phiQm*0.4159025) + (C_tZ*C_tZ*0.1900000+ 0.48*C_tG*C_tG +C_phit*C_phit*0.0133862  + C_phiQm*C_phiQm*0.024605025 - C_phit*C_phiQm*0.00894125 ))*(SM_ttZ_bin_70_110/ttZ_bin_70_110_madgraph_NLO)
                        + (0.0424788*C_tu8 + 0.04089323825*C_td8 + 0.0163593*C_tu8*C_tu8 + 0.01208651739*C_td8*C_td8 + 0.317150*C_Qq18 + 0.265025*C_tq8 -0.105523*C_Qq38 + 0.065980*C_Qu8 + 0.066175*C_Qd8)*(SM_ttZ_bin_70_110/ttZ_bin_70_110_madgraph_LO) ;
            }
            else{
                return  SM_ttZ_bin_70_110 + ((- C_tZ*0.02875+ 1.49*C_tG + C_phit*0.2634875 - C_phiQm*0.4159025))*(SM_ttZ_bin_70_110/ttZ_bin_70_110_madgraph_NLO)
                        + (0.0424788*C_tu8 + 0.04089323825*C_td8 + 0.317150*C_Qq18 + 0.265025*C_tq8 -0.105523*C_Qq38 + 0.065980*C_Qu8 + 0.066175*C_Qd8)*(SM_ttZ_bin_70_110/ttZ_bin_70_110_madgraph_LO) ;
            }
        }
    
        else{
            if(flag_Quadratic){
                return  SM_ttZ_bin_70_110 + ((- (-0.472123*C_tB + 0.881533*C_tW)*0.02875+ 1.49*C_tG + C_phit*0.2634875 - (C_phiQ1-C_phiQ3)*0.4159025) + ((-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW)*0.1900000+ 0.48*C_tG*C_tG +C_phit*C_phit*0.0133862  + (C_phiQ1-C_phiQ3)*(C_phiQ1-C_phiQ3)*0.024605025 - C_phit*(C_phiQ1-C_phiQ3)*0.00894125))*(SM_ttZ_bin_70_110/ttZ_bin_70_110_madgraph_NLO)
                        + (0.0424788*C_tu8 + 0.04089323825*C_td8 + 0.0163593*C_tu8*C_tu8 + 0.01208651739*C_td8*C_td8)*(SM_ttZ_bin_70_110/ttZ_bin_70_110_madgraph_LO);
            }
            else{
                return  SM_ttZ_bin_70_110 + ((- (-0.472123*C_tB + 0.881533*C_tW)*0.02875+ 1.49*C_tG + C_phit*0.2634875 - (C_phiQ1-C_phiQ3)*0.4159025))*(SM_ttZ_bin_70_110/ttZ_bin_70_110_madgraph_NLO)
                        + (0.0424788*C_tu8 + 0.04089323825*C_td8)*(SM_ttZ_bin_70_110/ttZ_bin_70_110_madgraph_LO);
            }
        }
        
    } else if(b_min == 110 && b_max == 160){
        
        double SM_ttZ_bin_110_160 = SM.getOptionalParameter("SM_sigma_ttz_bin_110_160_ATLAS_210312603");
        
        double ttZ_bin_110_160_madgraph_NLO = 3.008;//fb
        double ttZ_bin_110_160_madgraph_LO = 2.61775;//fb
 


        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_110_160 + ((-C_tZ*0.006+ 1.06*C_tG + C_phit*0.218294 - C_phiQm*0.326156) + (C_tZ*C_tZ*0.2215000+ 0.64*C_tG*C_tG+  C_phit*C_phit*0.0113604  + C_phiQm*C_phiQm*0.0240982 - C_phit*C_phiQm*0.00787972 ))*(SM_ttZ_bin_110_160/ttZ_bin_110_160_madgraph_NLO)
                        + (0.0382339*C_tu8 + 0.03493132846*C_td8 + 0.0130857*C_tu8*C_tu8 + 0.01034808693*C_td8*C_td8 + 0.230740*C_Qq18 + 0.188158*C_tq8 -0.062276*C_Qq38 + 0.058940*C_Qu8 + 0.053080*C_Qd8)*(SM_ttZ_bin_110_160/ttZ_bin_110_160_madgraph_LO);
            }
            else{
                return  SM_ttZ_bin_110_160 + ((-C_tZ*0.006+ 1.06*C_tG + C_phit*0.218294 - C_phiQm*0.326156))*(SM_ttZ_bin_110_160/ttZ_bin_110_160_madgraph_NLO)
                        + (0.0382339*C_tu8 + 0.03493132846*C_td8 + 0.230740*C_Qq18 + 0.188158*C_tq8 -0.062276*C_Qq38 + 0.058940*C_Qu8 + 0.053080*C_Qd8)*(SM_ttZ_bin_110_160/ttZ_bin_110_160_madgraph_LO);
            }
        }
   
        else{
            if(flag_Quadratic){
                return  SM_ttZ_bin_110_160 + ((-(-0.472123*C_tB + 0.881533*C_tW)*0.006+ 1.06*C_tG + C_phit*0.218294 - (C_phiQ1-C_phiQ3)*0.326156) + ((-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW)*0.2215000+ 0.64*C_tG*C_tG+  C_phit*C_phit*0.0113604  + (C_phiQ1-C_phiQ3)*(C_phiQ1-C_phiQ3)*0.0240982 - C_phit*(C_phiQ1-C_phiQ3)*0.00787972))*(SM_ttZ_bin_110_160/ttZ_bin_110_160_madgraph_NLO)
                        + (0.0382339*C_tu8 + 0.03493132846*C_td8 +  0.0130857*C_tu8*C_tu8 + 0.01034808693*C_td8*C_td8)*(SM_ttZ_bin_110_160/ttZ_bin_110_160_madgraph_LO);
            }
            else{
                return  SM_ttZ_bin_110_160 + ((-(-0.472123*C_tB + 0.881533*C_tW)*0.006+ 1.06*C_tG + C_phit*0.218294 - (C_phiQ1-C_phiQ3)*0.326156))*(SM_ttZ_bin_110_160/ttZ_bin_110_160_madgraph_NLO)
                        + (0.0382339*C_tu8 + 0.03493132846*C_td8)*(SM_ttZ_bin_110_160/ttZ_bin_110_160_madgraph_LO);
            }
        }
        
    } else if(b_min == 160 && b_max == 220){
        
        double SM_ttZ_bin_160_220 = SM.getOptionalParameter("SM_sigma_ttz_bin_160_220_ATLAS_210312603");
        
        double ttZ_bin_160_220_madgraph_NLO = 1.76667;//fb
        double ttZ_bin_160_220_madgraph_LO = 1.60498;//fb
        
        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_160_220  + ((-C_tZ*0.004166666667+ 0.6475*C_tG + C_phit*0.1374968333 - C_phiQm*0.2022933333) + (C_tZ*C_tZ*0.1941667+ 0.545*C_tG*C_tG + C_phit*C_phit*0.005696683333  + C_phiQm*C_phiQm*0.01418198333 - C_phit*C_phiQm*0.005393633333))*(SM_ttZ_bin_160_220/ttZ_bin_160_220_madgraph_NLO)
                        + (0.0296376*C_tu8 + 0.02403566565*C_td8 + 0.0141296*C_tu8*C_tu8 + 0.009765608033*C_td8*C_td8 + 0.145108*C_Qq18 + 0.117050*C_tq8 -0.029508*C_Qq38 + 0.044727*C_Qu8 + 0.036490*C_Qd8)*(SM_ttZ_bin_160_220/ttZ_bin_160_220_madgraph_LO);
            }
            else{
                return  SM_ttZ_bin_160_220  + ((-C_tZ*0.004166666667+ 0.6475*C_tG + C_phit*0.1374968333 - C_phiQm*0.2022933333))*(SM_ttZ_bin_160_220/ttZ_bin_160_220_madgraph_NLO)
                        + (0.0296376*C_tu8 + 0.02403566565*C_td8 + 0.145108*C_Qq18 + 0.117050*C_tq8 -0.029508*C_Qq38 + 0.044727*C_Qu8 + 0.036490*C_Qd8)*(SM_ttZ_bin_160_220/ttZ_bin_160_220_madgraph_LO);
            }
        }
    
        else{
            if(flag_Quadratic){
                return  SM_ttZ_bin_160_220  + ((-(-0.472123*C_tB + 0.881533*C_tW)*0.004166666667+ 0.6475*C_tG + C_phit*0.1374968333 - (C_phiQ1-C_phiQ3)*0.2022933333) + ((-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW)*0.1941667+ 0.545*C_tG*C_tG + C_phit*C_phit*0.005696683333  + (C_phiQ1-C_phiQ3)*(C_phiQ1-C_phiQ3)*0.01418198333 - C_phit*(C_phiQ1-C_phiQ3)*0.005393633333))*(SM_ttZ_bin_160_220/ttZ_bin_160_220_madgraph_NLO)
                        + (0.0296376*C_tu8 + 0.02403566565*C_td8 + 0.0141296*C_tu8*C_tu8 + 0.009765608033*C_td8*C_td8)*(SM_ttZ_bin_160_220/ttZ_bin_160_220_madgraph_LO);
            }
            else{
                return  SM_ttZ_bin_160_220  + ((-(-0.472123*C_tB + 0.881533*C_tW)*0.004166666667+ 0.6475*C_tG + C_phit*0.1374968333 - (C_phiQ1-C_phiQ3)*0.2022933333))*(SM_ttZ_bin_160_220/ttZ_bin_160_220_madgraph_NLO)
                        + (0.0296376*C_tu8 + 0.02403566565*C_td8)*(SM_ttZ_bin_160_220/ttZ_bin_160_220_madgraph_LO);
            }
        }
        
    } else if(b_min == 220 && b_max == 290){
        
        double SM_ttZ_bin_220_290 = SM.getOptionalParameter("SM_sigma_ttz_bin_220_290_ATLAS_210312603");
        
        double ttZ_bin_220_290_madgraph_NLO = 0.884;//fb
        double ttZ_bin_220_290_madgraph_LO = 0.83745;//fb
 

        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_220_290 + ((C_tZ*0.0002142857143+ 0.325*C_tG + C_phit*0.07184985714 - C_phiQm*0.1012151429) + (C_tZ*C_tZ*0.1330357+ 0.3917142857*C_tG*C_tG + C_phit*C_phit*0.003090857143  + C_phiQm*C_phiQm*0.007307557143 - C_phit*C_phiQm*0.004132214286))*(SM_ttZ_bin_220_290/ttZ_bin_220_290_madgraph_NLO)
                        + (0.0212577*C_tu8 + 0.01425832419*C_td8 + 0.0108656*C_tu8*C_tu8 + 0.00656670422*C_td8*C_td8 + 0.085773*C_Qq18 + 0.068533*C_tq8 -0.011413*C_Qq38 + 0.030494*C_Qu8 + 0.022837*C_Qd8)*(SM_ttZ_bin_220_290/ttZ_bin_220_290_madgraph_LO);
            }
            else{
                return  SM_ttZ_bin_220_290 + ((C_tZ*0.0002142857143+ 0.325*C_tG + C_phit*0.07184985714 - C_phiQm*0.1012151429))*(SM_ttZ_bin_220_290/ttZ_bin_220_290_madgraph_NLO)
                        + (0.0212577*C_tu8 + 0.01425832419*C_td8 + 0.085773*C_Qq18 + 0.068533*C_tq8 -0.011413*C_Qq38 + 0.030494*C_Qu8 + 0.022837*C_Qd8)*(SM_ttZ_bin_220_290/ttZ_bin_220_290_madgraph_LO);
            }
        }
        else{
            if(flag_Quadratic){
                return  SM_ttZ_bin_220_290 + (((-0.472123*C_tB + 0.881533*C_tW)*0.0002142857143+ 0.325*C_tG + C_phit*0.07184985714 - (C_phiQ1-C_phiQ3)*0.1012151429) + ((-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW)*0.1330357+ 0.3917142857*C_tG*C_tG + C_phit*C_phit*0.003090857143  + (C_phiQ1-C_phiQ3)*(C_phiQ1-C_phiQ3)*0.007307557143 - C_phit*(C_phiQ1-C_phiQ3)*0.004132214286))*(SM_ttZ_bin_220_290/ttZ_bin_220_290_madgraph_NLO)
                        + (0.0212577*C_tu8 + 0.01425832419*C_td8 + 0.0108656*C_tu8*C_tu8 + 0.00656670422*C_td8*C_td8)*(SM_ttZ_bin_220_290/ttZ_bin_220_290_madgraph_LO);
            }
            else{
                return  SM_ttZ_bin_220_290 + (((-0.472123*C_tB + 0.881533*C_tW)*0.0002142857143+ 0.325*C_tG + C_phit*0.07184985714 - (C_phiQ1-C_phiQ3)*0.1012151429))*(SM_ttZ_bin_220_290/ttZ_bin_220_290_madgraph_NLO)
                        + (0.0212577*C_tu8 + 0.01425832419*C_td8)*(SM_ttZ_bin_220_290/ttZ_bin_220_290_madgraph_LO);
            }
        }
        
    } else if(b_min == 290 && b_max == 400){
        
        double SM_ttZ_bin_290_400 = SM.getOptionalParameter("SM_sigma_ttz_bin_290_400_ATLAS_210312603");
        
        double ttZ_bin_290_400_madgraph_NLO = 0.338091;//fb
        double ttZ_bin_290_400_madgraph_LO = 0.33948;//fb
 

        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_290_400 + ((-C_tZ*0.001431818182+ 0.1359090909*C_tG + C_phit*0.02729918182 - C_phiQm*0.03635190909)+ (C_tZ*C_tZ*0.0778750+0.2721818182*C_tG*C_tG + C_phit*C_phit*0.001279009091 + C_phiQm*C_phiQm*0.002421036364- C_phit*C_phiQm*0.001086909091))*(SM_ttZ_bin_290_400/ttZ_bin_290_400_madgraph_NLO)
                        + (0.0129139*C_tu8 + 0.008451176896*C_td8 + 0.0084430*C_tu8*C_tu8 + 0.005011647401*C_td8*C_td8 + 0.044595*C_Qq18 + 0.035566*C_tq8 -0.002726*C_Qq38 + 0.018004*C_Qu8 + 0.012432*C_Qd8)*(SM_ttZ_bin_290_400/ttZ_bin_290_400_madgraph_LO);
            } 
            else{
                return  SM_ttZ_bin_290_400 + ((-C_tZ*0.001431818182+ 0.1359090909*C_tG + C_phit*0.02729918182 - C_phiQm*0.03635190909))*(SM_ttZ_bin_290_400/ttZ_bin_290_400_madgraph_NLO)
                        + (0.0129139*C_tu8 + 0.008451176896*C_td8 + 0.044595*C_Qq18 + 0.035566*C_tq8 -0.002726*C_Qq38 + 0.018004*C_Qu8 + 0.012432*C_Qd8)*(SM_ttZ_bin_290_400/ttZ_bin_290_400_madgraph_LO);
            }
        }
        else{
            if(flag_Quadratic){
                return  SM_ttZ_bin_290_400 + ((-(-0.472123*C_tB + 0.881533*C_tW)*0.001431818182+ 0.1359090909*C_tG + C_phit*0.02729918182 - (C_phiQ1-C_phiQ3)*0.03635190909)+ ((-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW)*0.0778750+0.2721818182*C_tG*C_tG + C_phit*C_phit*0.001279009091 + (C_phiQ1-C_phiQ3)*(C_phiQ1-C_phiQ3)*0.002421036364- C_phit*(C_phiQ1-C_phiQ3)*0.001086909091))*(SM_ttZ_bin_290_400/ttZ_bin_290_400_madgraph_NLO)
                        + (0.0129139*C_tu8 + 0.008451176896*C_td8 + 0.0084430*C_tu8*C_tu8 + 0.005011647401*C_td8*C_td8)*(SM_ttZ_bin_290_400/ttZ_bin_290_400_madgraph_LO);
            } 
            else{
                return  SM_ttZ_bin_290_400 + ((-(-0.472123*C_tB + 0.881533*C_tW)*0.001431818182+ 0.1359090909*C_tG + C_phit*0.02729918182 - (C_phiQ1-C_phiQ3)*0.03635190909))*(SM_ttZ_bin_290_400/ttZ_bin_290_400_madgraph_NLO)
                        + (0.0129139*C_tu8 + 0.008451176896*C_td8)*(SM_ttZ_bin_290_400/ttZ_bin_290_400_madgraph_LO);
            }
        }
        
    }
    else throw std::runtime_error("wrong bin choice in sigma_ttz_diff_LO_ATLAS_210312603");
   
}






sigma_ttz_diff_NLO_ATLAS_231204450::sigma_ttz_diff_NLO_ATLAS_231204450(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_ttz_bin_0_60_ATLAS_231204450" << "SM_sigma_ttz_bin_60_100_ATLAS_231204450"
            << "SM_sigma_ttz_bin_100_140_ATLAS_231204450" << "SM_sigma_ttz_bin_140_180_ATLAS_231204450" << "SM_sigma_ttz_bin_180_230_ATLAS_231204450" 
            << "SM_sigma_ttz_bin_230_280_ATLAS_231204450" << "SM_sigma_ttz_bin_280_350_ATLAS_231204450" << "SM_sigma_ttz_bin_350_1000_ATLAS_231204450");

    
}

double sigma_ttz_diff_NLO_ATLAS_231204450::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    //double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    //double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    
   
    if(b_min == 0 && b_max == 60){
        
        double SM_ttZ_bin_0_60 = SM.getOptionalParameter("SM_sigma_ttz_bin_0_60_ATLAS_231204450");
        double ttZ_bin_0_60_madgraph_NLO = 0.183649;//pb
 
        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_0_60 + (0.011269*C_tq8+0.002461*C_tq8*C_tq8+0.000099*C_tu1+0.002392*C_tu1*C_tu1+-0.000361*C_Qu1+0.004878*C_Qu1*C_Qu1+-0.017563*C_phiQm+0.000756*C_phiQm*C_phiQm+0.001524*C_tu8+0.000451*C_tu8*C_tu8+-0.005258*C_Qq38+0.003777*C_Qq38*C_Qq38+0.009789*C_phit+0.000833*C_phit*C_phit+0.001780*C_td8+0.000306*C_td8*C_td8+0.060984*C_tG+0.015753*C_tG*C_tG+-0.001235*C_tZ+0.005322*C_tZ*C_tZ+-0.000035*C_td1+0.001819*C_td1*C_td1+0.002219*C_Qd8+0.000618*C_Qd8*C_Qd8+-0.001912*C_tq1+0.017956*C_tq1*C_tq1+-0.000454*C_Qd1+0.004078*C_Qd1*C_Qd1+0.013873*C_Qq18+0.004197*C_Qq18*C_Qq18+0.001646*C_Qu8+0.000684*C_Qu8*C_Qu8)*(SM_ttZ_bin_0_60/ttZ_bin_0_60_madgraph_NLO);

            }
            else{
                return  SM_ttZ_bin_0_60 + (0.011269*C_tq8+0.000099*C_tu1+-0.000361*C_Qu1+-0.017563*C_phiQm+0.001524*C_tu8+-0.005258*C_Qq38+0.009789*C_phit+0.001780*C_td8+0.060984*C_tG+-0.001235*C_tZ+-0.000035*C_td1+0.002219*C_Qd8+-0.001912*C_tq1+-0.000454*C_Qd1+0.013873*C_Qq18+0.001646*C_Qu8)*(SM_ttZ_bin_0_60/ttZ_bin_0_60_madgraph_NLO);
            }
        }
        else{
            throw std::runtime_error("sigma_ttz_diff_LO_ATLAS_231204450 is only available in the LHC top-quark WC basis");            
        }
        
    } else if(b_min == 60 && b_max == 100){
        
        double SM_ttZ_bin_60_100 = SM.getOptionalParameter("SM_sigma_ttz_bin_60_100_ATLAS_231204450");
        
        double ttZ_bin_60_100_madgraph_NLO = 0.191991;//pb
        
 

        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_60_100  + (0.010217*C_tq8+0.002213*C_tq8*C_tq8+-0.000137*C_tu1+0.002743*C_tu1*C_tu1+-0.000214*C_Qu1+0.005965*C_Qu1*C_Qu1+-0.019223*C_phiQm+0.000459*C_phiQm*C_phiQm+0.001520*C_tu8+0.000257*C_tu8*C_tu8+-0.004565*C_Qq38+0.003805*C_Qq38*C_Qq38+0.011644*C_phit+0.000647*C_phit*C_phit+0.001569*C_td8+0.000268*C_td8*C_td8+0.065876*C_tG+0.021105*C_tG*C_tG+-0.000950*C_tZ+0.008327*C_tZ*C_tZ+-0.000064*C_td1+0.002031*C_td1*C_td1+0.002114*C_Qd8+0.000501*C_Qd8*C_Qd8+-0.001514*C_tq1+0.017893*C_tq1*C_tq1+-0.000091*C_Qd1+0.004975*C_Qd1*C_Qd1+0.012940*C_Qq18+0.004197*C_Qq18*C_Qq18+0.002223*C_Qu8+0.000797*C_Qu8*C_Qu8)*(SM_ttZ_bin_60_100/ttZ_bin_60_100_madgraph_NLO);
            }
            else{
                return SM_ttZ_bin_60_100  + (0.010217*C_tq8+-0.000137*C_tu1+-0.000214*C_Qu1+-0.019223*C_phiQm+0.001520*C_tu8+-0.004565*C_Qq38+0.011644*C_phit+0.001569*C_td8+0.065876*C_tG+-0.000950*C_tZ+-0.000064*C_td1+0.002114*C_Qd8+-0.001514*C_tq1+-0.000091*C_Qd1+0.012940*C_Qq18+0.002223*C_Qu8)*(SM_ttZ_bin_60_100/ttZ_bin_60_100_madgraph_NLO);
            }
        }
    
        else{
            throw std::runtime_error("sigma_ttz_diff_LO_ATLAS_231204450 is only available in the LHC top-quark WC basis");
        }
        
    } else if(b_min == 100 && b_max == 140){
        
        double SM_ttZ_bin_100_140 = SM.getOptionalParameter("SM_sigma_ttz_bin_100_140_ATLAS_231204450");        
        
        double ttZ_bin_100_140_madgraph_NLO = 0.155546 ;//pb
      
 

        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_100_140 + (0.007390*C_tq8+0.002182*C_tq8*C_tq8+-0.000099*C_tu1+0.003344*C_tu1*C_tu1+-0.000584*C_Qu1+0.006721*C_Qu1*C_Qu1+-0.016592*C_phiQm+0.000727*C_phiQm*C_phiQm+0.001161*C_tu8+0.000590*C_tu8*C_tu8+-0.002943*C_Qq38+0.003662*C_Qq38*C_Qq38+0.010515*C_phit+0.000744*C_phit*C_phit+0.001235*C_td8+0.000617*C_td8*C_td8+0.055205*C_tG+0.024034*C_tG*C_tG+-0.000824*C_tZ+0.010011*C_tZ*C_tZ+-0.000175*C_td1+0.002214*C_td1*C_td1+0.001987*C_Qd8+0.000772*C_Qd8*C_Qd8+-0.001141*C_tq1+0.015456*C_tq1*C_tq1+-0.000606*C_Qd1+0.004785*C_Qd1*C_Qd1+0.009693*C_Qq18+0.003900*C_Qq18*C_Qq18+0.002068*C_Qu8+0.000865*C_Qu8*C_Qu8)*(SM_ttZ_bin_100_140/ttZ_bin_100_140_madgraph_NLO) ;
            }
            else{
                return  SM_ttZ_bin_100_140 + (0.007390*C_tq8+-0.000099*C_tu1+-0.000584*C_Qu1+-0.016592*C_phiQm+0.001161*C_tu8+-0.002943*C_Qq38+0.010515*C_phit+0.001235*C_td8+0.055205*C_tG+-0.000824*C_tZ+-0.000175*C_td1+0.001987*C_Qd8+-0.001141*C_tq1+-0.000606*C_Qd1+0.009693*C_Qq18+0.002068*C_Qu8)*(SM_ttZ_bin_100_140/ttZ_bin_100_140_madgraph_NLO);
            }
        }
    
        else{
            throw std::runtime_error("sigma_ttz_diff_LO_ATLAS_231204450 is only available in the LHC top-quark WC basis");
        }
        
    } else if(b_min == 140 && b_max == 180){
        
        double SM_ttZ_bin_140_180 = SM.getOptionalParameter("SM_sigma_ttz_bin_140_180_ATLAS_231204450");
        
        double ttZ_bin_140_180_madgraph_NLO = 0.111076;//pb
 


        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_140_180 + (0.005149*C_tq8+0.001471*C_tq8*C_tq8+0.000229*C_tu1+0.002882*C_tu1*C_tu1+-0.000283*C_Qu1+0.005958*C_Qu1*C_Qu1+-0.012197*C_phiQm+0.000413*C_phiQm*C_phiQm+0.001216*C_tu8+0.000399*C_tu8*C_tu8+-0.001704*C_Qq38+0.002780*C_Qq38*C_Qq38+0.008411*C_phit+0.000428*C_phit*C_phit+0.000969*C_td8+0.000265*C_td8*C_td8+0.039782*C_tG+0.023535*C_tG*C_tG+-0.000231*C_tZ+0.009815*C_tZ*C_tZ+0.000188*C_td1+0.001906*C_td1*C_td1+0.001649*C_Qd8+0.000513*C_Qd8*C_Qd8+-0.000742*C_tq1+0.011939*C_tq1*C_tq1+-0.000163*C_Qd1+0.004081*C_Qd1*C_Qd1+0.006866*C_Qq18+0.003056*C_Qq18*C_Qq18+0.001595*C_Qu8+0.000664*C_Qu8*C_Qu8)*(SM_ttZ_bin_140_180/ttZ_bin_140_180_madgraph_NLO);
            }
            else{
                return  SM_ttZ_bin_140_180 + (0.005149*C_tq8+0.000229*C_tu1+-0.000283*C_Qu1+-0.012197*C_phiQm+0.001216*C_tu8+-0.001704*C_Qq38+0.008411*C_phit+0.000969*C_td8+0.039782*C_tG+-0.000231*C_tZ+0.000188*C_td1+0.001649*C_Qd8+-0.000742*C_tq1+-0.000163*C_Qd1+0.006866*C_Qq18+0.001595*C_Qu8)*(SM_ttZ_bin_140_180/ttZ_bin_140_180_madgraph_NLO);
            }
        }
   
        else{
            throw std::runtime_error("sigma_ttz_diff_LO_ATLAS_231204450 is only available in the LHC top-quark WC basis");
        }
        
    } else if(b_min == 180 && b_max == 230){
        
        double SM_ttZ_bin_180_230 = SM.getOptionalParameter("SM_sigma_ttz_bin_180_230_ATLAS_231204450");
        
        double ttZ_bin_180_230_madgraph_NLO = 0.088653;//pb
        
        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_180_230  + (0.004128*C_tq8+0.001486*C_tq8*C_tq8+-0.000085*C_tu1+0.003255*C_tu1*C_tu1+-0.000282*C_Qu1+0.006308*C_Qu1*C_Qu1+-0.009886*C_phiQm+0.000290*C_phiQm*C_phiQm+0.001162*C_tu8+0.000569*C_tu8*C_tu8+-0.001204*C_Qq38+0.002711*C_Qq38*C_Qq38+0.007089*C_phit+0.000340*C_phit*C_phit+0.000804*C_td8+0.000341*C_td8*C_td8+0.032375*C_tG+0.027076*C_tG*C_tG+-0.000076*C_tZ+0.010482*C_tZ*C_tZ+0.000104*C_td1+0.001906*C_td1*C_td1+0.001263*C_Qd8+0.000498*C_Qd8*C_Qd8+-0.000831*C_tq1+0.011204*C_tq1*C_tq1+-0.000296*C_Qd1+0.004084*C_Qd1*C_Qd1+0.005648*C_Qq18+0.002898*C_Qq18*C_Qq18+0.001490*C_Qu8+0.000829*C_Qu8*C_Qu8)*(SM_ttZ_bin_180_230/ttZ_bin_180_230_madgraph_NLO);
            }
            else{
                return  SM_ttZ_bin_180_230  + (0.004128*C_tq8+-0.000085*C_tu1+-0.000282*C_Qu1+-0.009886*C_phiQm+0.001162*C_tu8+-0.001204*C_Qq38+0.007089*C_phit+0.000804*C_td8+0.032375*C_tG+-0.000076*C_tZ+0.000104*C_td1+0.001263*C_Qd8+-0.000831*C_tq1+-0.000296*C_Qd1+0.005648*C_Qq18+0.001490*C_Qu8)*(SM_ttZ_bin_180_230/ttZ_bin_180_230_madgraph_NLO);
            }
        }
    
        else{
            throw std::runtime_error("sigma_ttz_diff_LO_ATLAS_231204450 is only available in the LHC top-quark WC basis");
        }
        
    } else if(b_min == 230 && b_max == 280){
        
        double SM_ttZ_bin_230_280 = SM.getOptionalParameter("SM_sigma_ttz_bin_230_280_ATLAS_231204450");
        
        double ttZ_bin_230_280_madgraph_NLO = 0.051815;//fpb
 

        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_230_280 + (0.002612*C_tq8+0.001122*C_tq8*C_tq8+0.000002*C_tu1+0.002810*C_tu1*C_tu1+-0.000086*C_Qu1+0.005276*C_Qu1*C_Qu1+-0.005801*C_phiQm+0.000132*C_phiQm*C_phiQm+0.000799*C_tu8+0.000436*C_tu8*C_tu8+-0.000580*C_Qq38+0.002063*C_Qq38*C_Qq38+0.004124*C_phit+0.000221*C_phit*C_phit+0.000669*C_td8+0.000234*C_td8*C_td8+0.019259*C_tG+0.023624*C_tG*C_tG+-0.000036*C_tZ+0.008174*C_tZ*C_tZ+0.000047*C_td1+0.001629*C_td1*C_td1+0.000868*C_Qd8+0.000364*C_Qd8*C_Qd8+-0.000483*C_tq1+0.008390*C_tq1*C_tq1+-0.000223*C_Qd1+0.003306*C_Qd1*C_Qd1+0.003561*C_Qq18+0.002188*C_Qq18*C_Qq18+0.000920*C_Qu8+0.000649*C_Qu8*C_Qu8)*(SM_ttZ_bin_230_280/ttZ_bin_230_280_madgraph_NLO);
            }
            else{
                return  SM_ttZ_bin_230_280 + (0.002612*C_tq8+0.000002*C_tu1+-0.000086*C_Qu1+-0.005801*C_phiQm+0.000799*C_tu8+-0.000580*C_Qq38+0.004124*C_phit+0.000669*C_td8+0.019259*C_tG+-0.000036*C_tZ+0.000047*C_td1+0.000868*C_Qd8+-0.000483*C_tq1+-0.000223*C_Qd1+0.003561*C_Qq18+0.000920*C_Qu8)*(SM_ttZ_bin_230_280/ttZ_bin_230_280_madgraph_NLO);
            }
        }
        else{
            throw std::runtime_error("sigma_ttz_diff_LO_ATLAS_231204450 is only available in the LHC top-quark WC basis");
        }
        
    } else if(b_min == 280 && b_max == 350){
        
        double SM_ttZ_bin_280_350 = SM.getOptionalParameter("SM_sigma_ttz_bin_280_350_ATLAS_231204450");
        
        double ttZ_bin_280_350_madgraph_NLO = 0.037860;//pb
 

        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_280_350 + (0.002148*C_tq8+0.001112*C_tq8*C_tq8+0.000106*C_tu1+0.003321*C_tu1*C_tu1+
                        -0.000285*C_Qu1+0.005863*C_Qu1*C_Qu1+-0.004249*C_phiQm+0.000122*C_phiQm*C_phiQm+0.000737*C_tu8+
                        0.000571*C_tu8*C_tu8+-0.000383*C_Qq38+0.002095*C_Qq38*C_Qq38+0.003095*C_phit+0.000121*C_phit*C_phit+
                        0.000597*C_td8+0.000346*C_td8*C_td8+0.014211*C_tG+0.026930*C_tG*C_tG+-0.000123*C_tZ+0.007960*C_tZ*C_tZ+
                        0.000170*C_td1+0.001903*C_td1*C_td1+0.000673*C_Qd8+0.000452*C_Qd8*C_Qd8+-0.000567*C_tq1+
                        0.008639*C_tq1*C_tq1+-0.000238*C_Qd1+0.003532*C_Qd1*C_Qd1+0.003025*C_Qq18+0.002312*C_Qq18*C_Qq18+
                        0.000878*C_Qu8+0.000699*C_Qu8*C_Qu8)*(SM_ttZ_bin_280_350/ttZ_bin_280_350_madgraph_NLO);
            } 
            else{
                return  SM_ttZ_bin_280_350 + (0.002148*C_tq8+0.000106*C_tu1+-0.000285*C_Qu1+-0.004249*C_phiQm+
                        0.000737*C_tu8+-0.000383*C_Qq38+0.003095*C_phit+0.000597*C_td8+0.014211*C_tG+
                        -0.000123*C_tZ+0.000170*C_td1+0.000673*C_Qd8+-0.000567*C_tq1+-0.000238*C_Qd1+
                        0.003025*C_Qq18+0.000878*C_Qu8)*(SM_ttZ_bin_280_350/ttZ_bin_280_350_madgraph_NLO);
            }
        }
        else{
            throw std::runtime_error("sigma_ttz_diff_LO_ATLAS_231204450 is only available in the LHC top-quark WC basis");
        }
        
    } else if(b_min == 350 && b_max == 1000){
        
        double SM_ttZ_bin_350_1000 = SM.getOptionalParameter("SM_sigma_ttz_bin_350_1000_ATLAS_231204450");
        
        double ttZ_bin_350_1000_madgraph_NLO = 0.035227;//pb
 

        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttZ_bin_350_1000 + (0.002977*C_tq8+0.002713*C_tq8*C_tq8+0.000208*C_tu1+0.011277*C_tu1*C_tu1+
                        -0.000674*C_Qu1+0.018149*C_Qu1*C_Qu1+-0.003756*C_phiQm+0.000147*C_phiQm*C_phiQm+0.001383*C_tu8+
                        0.001825*C_tu8*C_tu8+-0.000294*C_Qq38+0.005212*C_Qq38*C_Qq38+0.002782*C_phit+0.000149*C_phit*C_phit+
                        0.000900*C_td8+0.000949*C_td8*C_td8+0.014178*C_tG+0.077992*C_tG*C_tG+-0.000043*C_tZ+0.013238*C_tZ*C_tZ+
                        0.000108*C_td1+0.005560*C_td1*C_td1+0.001033*C_Qd8+0.001103*C_Qd8*C_Qd8+-0.000872*C_tq1+
                        0.022731*C_tq1*C_tq1+-0.000480*C_Qd1+0.009184*C_Qd1*C_Qd1+0.004429*C_Qq18+0.005979*C_Qq18*C_Qq18+
                        0.001415*C_Qu8+0.001955*C_Qu8*C_Qu8)*(SM_ttZ_bin_350_1000/ttZ_bin_350_1000_madgraph_NLO);
            } 
            else{
                return  SM_ttZ_bin_350_1000 + (0.002977*C_tq8+0.000208*C_tu1+-0.000674*C_Qu1+-0.003756*C_phiQm+
                        0.001383*C_tu8+-0.000294*C_Qq38+0.002782*C_phit+0.000900*C_td8+0.014178*C_tG+-0.000043*C_tZ+
                        0.000108*C_td1+0.001033*C_Qd8+-0.000872*C_tq1+-0.000480*C_Qd1+0.004429*C_Qq18+
                        0.001415*C_Qu8)*(SM_ttZ_bin_350_1000/ttZ_bin_350_1000_madgraph_NLO);
            }
        }
        else{
            throw std::runtime_error("sigma_ttz_diff_LO_ATLAS_231204450 is only available in the LHC top-quark WC basis");
        }
        
    }
    else throw std::runtime_error("wrong bin choice in sigma_ttz_diff_LO_ATLAS_231204450");
   
}









/////////////////////////////////////////////////////////////////////////////////////////////////////////



/////// ttA differential cross section for different bins ///////////////////////////////////////////////


sigma_tta_diff_NLO_ATLAS_emu_200706946::sigma_tta_diff_NLO_ATLAS_emu_200706946(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tta_bin_20_25_ATLAS_emu" << "SM_sigma_tta_bin_25_30_ATLAS_emu"
            << "SM_sigma_tta_bin_30_35_ATLAS_emu" << "SM_sigma_tta_bin_35_40_ATLAS_emu" << "SM_sigma_tta_bin_40_47_ATLAS_emu" 
            << "SM_sigma_tta_bin_47_55_ATLAS_emu" << "SM_sigma_tta_bin_55_70_ATLAS_emu" << "SM_sigma_tta_bin_70_85_ATLAS_emu"
            << "SM_sigma_tta_bin_85_132_ATLAS_emu" << "SM_sigma_tta_bin_132_180_ATLAS_emu" << "SM_sigma_tta_bin_180_300_ATLAS_emu");

    
}

double sigma_tta_diff_NLO_ATLAS_emu_200706946::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();

    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    

    if(b_min == 20 && b_max == 25){
        
        double SM_ttA_bin_20_25 = SM.getOptionalParameter("SM_sigma_tta_bin_20_25_ATLAS_emu");
        
        double ttA_bin_20_25_madgraph_NLO = 1.38765;//fb
        double ttA_bin_20_25_madgraph_LO = 1.10093;//fb

    
        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttA_bin_20_25 + (-0.073402*C_tZ + 0.004305*C_tZ*C_tZ + 0.035194*C_tW + 0.008250*C_tW*C_tW + 0.003315802469*C_tZ*C_tW + 0.4024691358*C_tG + 0.09382716049*C_tG*C_tG)*(SM_ttA_bin_20_25/ttA_bin_20_25_madgraph_NLO) 
                        + (0.0324910*C_tu8 + 0.0159132*C_td8 + 0.0082019*C_tu8*C_tu8 + 0.0042998*C_td8*C_td8 + 0.060800*C_Qq18 + 0.077047*C_tq8 + 0.024578*C_Qq38 + 0.066356*C_Qu8 + 0.010860*C_Qd8)*(SM_ttA_bin_20_25/ttA_bin_20_25_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_20_25 + ((-0.073402*C_tZ  + 0.035194*C_tW + 0.4024691358*C_tG))*(SM_ttA_bin_20_25/ttA_bin_20_25_madgraph_NLO)
                    + (0.0324910*C_tu8 + 0.0159132*C_td8 + 0.060800*C_Qq18 + 0.077047*C_tq8 + 0.024578*C_Qq38 + 0.066356*C_Qu8 + 0.010860*C_Qd8)*(SM_ttA_bin_20_25/ttA_bin_20_25_madgraph_LO);
            }
        }
        else{
            if(flag_Quadratic){
                return  SM_ttA_bin_20_25 + ((-0.073402*(-0.472123*C_tB + 0.881533*C_tW) + 0.004305*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 0.035194*C_tW + 0.008250*C_tW*C_tW + 0.003315802469*(-0.472123*C_tB + 0.881533*C_tW)*C_tW + 0.4024691358*C_tG + 0.09382716049*C_tG*C_tG ))*(SM_ttA_bin_20_25/ttA_bin_20_25_madgraph_NLO)
                        + (0.0324910*C_tu8 + 0.0159132*C_td8 + 0.0082019*C_tu8*C_tu8 + 0.0042998*C_td8*C_td8)*(SM_ttA_bin_20_25/ttA_bin_20_25_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_20_25 + ((-0.073402*(-0.472123*C_tB + 0.881533*C_tW)  + 0.035194*C_tW + 0.4024691358*C_tG ))*(SM_ttA_bin_20_25/ttA_bin_20_25_madgraph_NLO)
                        + (0.0324910*C_tu8 + 0.0159132*C_td8)*(SM_ttA_bin_20_25/ttA_bin_20_25_madgraph_LO);
            }
        }
    } else if(b_min == 25 && b_max == 30){
        
    
        double SM_ttA_bin_25_30 = SM.getOptionalParameter("SM_sigma_tta_bin_25_30_ATLAS_emu");
        
        double ttA_bin_25_30_madgraph_NLO = 1.11605;//fb
        double ttA_bin_25_30_madgraph_LO = 0.84869;//fb


        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttA_bin_25_30 + (-0.018930*C_tZ + 0.057613*C_tZ*C_tZ + 0.033745*C_tW + 0.067490*C_tW*C_tW - 0.1200617284*C_tZ*C_tW 	+ 0.3209876543*C_tG + 0.04814814815*C_tG*C_tG)*(SM_ttA_bin_25_30/ttA_bin_25_30_madgraph_NLO)
                        + (0.0226520*C_tu8 + 0.0126074*C_td8 + 0.0063043*C_tu8*C_tu8 + 0.0034749*C_td8*C_td8 + 0.045347*C_Qq18 + 0.058257*C_tq8 + 0.017563*C_Qq38 + 0.050267*C_Qu8 + 0.008151*C_Qd8)*(SM_ttA_bin_25_30/ttA_bin_25_30_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_25_30 + (-0.018930*C_tZ  + 0.033745*C_tW + 0.3209876543*C_tG)*(SM_ttA_bin_25_30/ttA_bin_25_30_madgraph_NLO)
                        + (0.0226520*C_tu8 + 0.0126074*C_td8 + 0.045347*C_Qq18 + 0.058257*C_tq8 + 0.017563*C_Qq38 + 0.050267*C_Qu8 + 0.008151*C_Qd8)*(SM_ttA_bin_25_30/ttA_bin_25_30_madgraph_LO);
            }
        }
        else{
            if(flag_Quadratic){
                return  SM_ttA_bin_25_30 + (-0.018930*(-0.472123*C_tB + 0.881533*C_tW) + 0.057613*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 0.033745*C_tW + 0.067490*C_tW*C_tW - 0.1200617284*(-0.472123*C_tB + 0.881533*C_tW)*C_tW 	+ 0.3209876543*C_tG + 0.04814814815*C_tG*C_tG )*(SM_ttA_bin_25_30/ttA_bin_25_30_madgraph_NLO)
                        + (0.0226520*C_tu8 + 0.0126074*C_td8+ 0.0063043*C_tu8*C_tu8 + 0.0034749*C_td8*C_td8)*(SM_ttA_bin_25_30/ttA_bin_25_30_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_25_30 + (-0.018930*(-0.472123*C_tB + 0.881533*C_tW)  + 0.033745*C_tW + 0.3209876543*C_tG)*(SM_ttA_bin_25_30/ttA_bin_25_30_madgraph_NLO)
                        + (0.0226520*C_tu8 + 0.0126074*C_td8)*(SM_ttA_bin_25_30/ttA_bin_25_30_madgraph_LO);
            }
        }
    } else if(b_min == 30 && b_max == 35){
        
    
        double SM_ttA_bin_30_35 = SM.getOptionalParameter("SM_sigma_tta_bin_30_35_ATLAS_emu");
        
        double ttA_bin_30_35_madgraph_NLO = 0.898765;//fb
        double ttA_bin_30_35_madgraph_LO = 0.67594;//fb

    
        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return   SM_ttA_bin_30_35 + (-0.005267*C_tZ + 0.010076*C_tZ*C_tZ - 0.011981*C_tW + 0.007265*C_tW*C_tW - 0.02546444444*C_tZ*C_tW + 0.2443209877*C_tG + 0.04327407407*C_tG*C_tG)*(SM_ttA_bin_30_35/ttA_bin_30_35_madgraph_NLO)
                        + (0.0190330*C_tu8 + 0.0105528*C_td8+ 0.0056393*C_tu8*C_tu8 + 0.0038102*C_td8*C_td8 + 0.035395*C_Qq18 + 0.045853*C_tq8 + 0.013175*C_Qq38 + 0.039633*C_Qu8 + 0.006393*C_Qd8)*(SM_ttA_bin_30_35/ttA_bin_30_35_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_30_35 + (-0.005267*C_tZ  - 0.011981*C_tW  + 0.2443209877*C_tG)*(SM_ttA_bin_30_35/ttA_bin_30_35_madgraph_NLO)
                        + (0.0190330*C_tu8 + 0.0105528*C_td8 + 0.035395*C_Qq18 + 0.045853*C_tq8 + 0.013175*C_Qq38 + 0.039633*C_Qu8 + 0.006393*C_Qd8)*(SM_ttA_bin_30_35/ttA_bin_30_35_madgraph_LO);
            }
        }
    
        else{
            if(flag_Quadratic){
                return   SM_ttA_bin_30_35 + (-0.005267*(-0.472123*C_tB + 0.881533*C_tW) + 0.010076*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) - 0.011981*C_tW + 0.007265*C_tW*C_tW - 0.02546444444*(-0.472123*C_tB + 0.881533*C_tW)*C_tW + 0.2443209877*C_tG + 0.04327407407*C_tG*C_tG)*(SM_ttA_bin_30_35/ttA_bin_30_35_madgraph_NLO)
                        + (0.0190330*C_tu8 + 0.0105528*C_td8+ 0.0056393*C_tu8*C_tu8 + 0.0038102*C_td8*C_td8)*(SM_ttA_bin_30_35/ttA_bin_30_35_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_30_35 + (-0.005267*(-0.472123*C_tB + 0.881533*C_tW)  - 0.011981*C_tW  + 0.2443209877*C_tG)*(SM_ttA_bin_30_35/ttA_bin_30_35_madgraph_NLO)
                        + (0.0190330*C_tu8 + 0.0105528*C_td8)*(SM_ttA_bin_30_35/ttA_bin_30_35_madgraph_LO);
            }
        }
    } else if(b_min == 35 && b_max == 40){
        
        double SM_ttA_bin_35_40 = SM.getOptionalParameter("SM_sigma_tta_bin_35_40_ATLAS_emu");
        
        double ttA_bin_35_40_madgraph_NLO= 0.730864;//fb
        double ttA_bin_35_40_madgraph_LO = 0.55527;//fb

        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttA_bin_35_40 + (0.013992*C_tZ  -0.002881*C_tZ*C_tZ + 0.013169*C_tW + 0.036626*C_tW*C_tW - 0.05092592593*C_tZ*C_tW + 0.222962963*C_tG +  0.04901234568*C_tG*C_tG )*(SM_ttA_bin_35_40/ttA_bin_35_40_madgraph_NLO)
                        + (0.0155991*C_tu8 + 0.0090981*C_td8 + 0.0050532*C_tu8*C_tu8 + 0.0030395*C_td8*C_td8 + 0.028450*C_Qq18 + 0.037207*C_tq8 + 0.010212*C_Qq38 + 0.032156*C_Qu8 + 0.005149*C_Qd8)*(SM_ttA_bin_35_40/ttA_bin_35_40_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_35_40 + (0.013992*C_tZ  + 0.013169*C_tW  + 0.222962963*C_tG )*(SM_ttA_bin_35_40/ttA_bin_35_40_madgraph_NLO)
                        + (0.0155991*C_tu8 + 0.0090981*C_td8 + 0.028450*C_Qq18 + 0.037207*C_tq8 + 0.010212*C_Qq38 + 0.032156*C_Qu8 + 0.005149*C_Qd8)*(SM_ttA_bin_35_40/ttA_bin_35_40_madgraph_LO);
            }
        }
    
        else{
            if(flag_Quadratic){
                return  SM_ttA_bin_35_40 + (0.013992*(-0.472123*C_tB + 0.881533*C_tW)  -0.002881*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 0.013169*C_tW + 0.036626*C_tW*C_tW - 0.05092592593*(-0.472123*C_tB + 0.881533*C_tW)*C_tW + 0.222962963*C_tG +  0.04901234568*C_tG*C_tG)*(SM_ttA_bin_35_40/ttA_bin_35_40_madgraph_NLO)
                        + (0.0155991*C_tu8 + 0.0090981*C_td8 + 0.0050532*C_tu8*C_tu8 + 0.0030395*C_td8*C_td8 )*(SM_ttA_bin_35_40/ttA_bin_35_40_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_35_40 + (0.013992*(-0.472123*C_tB + 0.881533*C_tW)  + 0.013169*C_tW  + 0.222962963*C_tG)*(SM_ttA_bin_35_40/ttA_bin_35_40_madgraph_NLO)
                        + (0.0155991*C_tu8 + 0.0090981*C_td8)*(SM_ttA_bin_35_40/ttA_bin_35_40_madgraph_LO);
            }
        }
    } else if(b_min == 40 && b_max == 47){
        
    
        double SM_ttA_bin_40_47 = SM.getOptionalParameter("SM_sigma_tta_bin_40_47_ATLAS_emu");
        
        double ttA_bin_40_47_madgraph_NLO= 0.592945;//fb
        double ttA_bin_40_47_madgraph_LO = 0.44974;//fb


        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttA_bin_40_47 + (-0.064903*C_tZ + 0.029189*C_tZ*C_tZ + 0.060024*C_tW + 0.037566*C_tW*C_tW - 0.06441798942*C_tZ*C_tW + 0.1672839506*C_tG + 0.02561552028*C_tG*C_tG )*(SM_ttA_bin_40_47/ttA_bin_40_47_madgraph_NLO)
                        +(0.0124768*C_tu8 + 0.0068320*C_td8 + 0.0045958*C_tu8*C_tu8 + 0.0026133*C_td8*C_td8 + 0.022585*C_Qq18 + 0.029873*C_tq8 + 0.007806*C_Qq38 + 0.025840*C_Qu8 + 0.004114*C_Qd8)*(SM_ttA_bin_40_47/ttA_bin_40_47_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_40_47 + (-0.064903*C_tZ + 0.060024*C_tW + 0.1672839506*C_tG)*(SM_ttA_bin_40_47/ttA_bin_40_47_madgraph_NLO)
                        + (0.0124768*C_tu8 + 0.0068320*C_td8 + 0.022585*C_Qq18 + 0.029873*C_tq8 + 0.007806*C_Qq38 + 0.025840*C_Qu8 + 0.004114*C_Qd8)*(SM_ttA_bin_40_47/ttA_bin_40_47_madgraph_LO);
            }
        }
        else{
            if(flag_Quadratic){
                return  SM_ttA_bin_40_47 + (-0.064903*(-0.472123*C_tB + 0.881533*C_tW) + 0.029189*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 0.060024*C_tW + 0.037566*C_tW*C_tW - 0.06441798942*(-0.472123*C_tB + 0.881533*C_tW)*C_tW + 0.1672839506*C_tG + 0.02561552028*C_tG*C_tG)*(SM_ttA_bin_40_47/ttA_bin_40_47_madgraph_NLO)
                        + (0.0124768*C_tu8 + 0.0068320*C_td8  + 0.0045958*C_tu8*C_tu8 + 0.0026133*C_td8*C_td8)*(SM_ttA_bin_40_47/ttA_bin_40_47_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_40_47 + (-0.064903*(-0.472123*C_tB + 0.881533*C_tW) + 0.060024*C_tW + 0.1672839506*C_tG )*(SM_ttA_bin_40_47/ttA_bin_40_47_madgraph_NLO)
                        + (0.0124768*C_tu8 + 0.0068320*C_td8)*(SM_ttA_bin_40_47/ttA_bin_40_47_madgraph_LO);
            }
        }
    } else if(b_min == 47 && b_max == 55){
        
    
        double SM_ttA_bin_47_55 = SM.getOptionalParameter("SM_sigma_tta_bin_47_55_ATLAS_emu");
        
        double ttA_bin_47_55_madgraph_NLO= 0.441049;//fb
        double ttA_bin_47_55_madgraph_LO = 0.35670;//fb


    
        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttA_bin_47_55 + (-0.022374*C_tZ + 0.033832*C_tZ*C_tZ + 0.024549*C_tW + 0.038170*C_tW*C_tW - 0.05701790123*C_tZ*C_tW + 0.1261574074*C_tG + 0.03036419753*C_tG*C_tG)*(SM_ttA_bin_47_55/ttA_bin_47_55_madgraph_NLO)
                        + (0.0097075*C_tu8 + 0.0053889*C_td8 + 0.0034068*C_tu8*C_tu8 + 0.0021192*C_td8*C_td8 + 0.017461*C_Qq18 + 0.023408*C_tq8 + 0.005779*C_Qq38 + 0.020298*C_Qu8 + 0.003200*C_Qd8)*(SM_ttA_bin_47_55/ttA_bin_47_55_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_47_55 + (-0.022374*C_tZ  + 0.024549*C_tW + 0.1261574074*C_tG )*(SM_ttA_bin_47_55/ttA_bin_47_55_madgraph_NLO)
                        + (0.0097075*C_tu8 + 0.0053889*C_td8 + 0.017461*C_Qq18 + 0.023408*C_tq8 + 0.005779*C_Qq38 + 0.020298*C_Qu8 + 0.003200*C_Qd8)*(SM_ttA_bin_47_55/ttA_bin_47_55_madgraph_LO);
            }
        }
        else{
            if(flag_Quadratic){
                return  SM_ttA_bin_47_55 + (-0.022374*(-0.472123*C_tB + 0.881533*C_tW) + 0.033832*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 0.024549*C_tW + 0.038170*C_tW*C_tW - 0.05701790123*(-0.472123*C_tB + 0.881533*C_tW)*C_tW + 0.1261574074*C_tG + 0.03036419753*C_tG*C_tG )*(SM_ttA_bin_47_55/ttA_bin_47_55_madgraph_NLO)
                        + (0.0097075*C_tu8 + 0.0053889*C_td8 + 0.0034068*C_tu8*C_tu8 + 0.0021192*C_td8*C_td8)*(SM_ttA_bin_47_55/ttA_bin_47_55_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_47_55 + (-0.022374*(-0.472123*C_tB + 0.881533*C_tW)  + 0.024549*C_tW + 0.1261574074*C_tG )*(SM_ttA_bin_47_55/ttA_bin_47_55_madgraph_NLO)
                        + (0.0097075*C_tu8 + 0.0053889*C_td8 )*(SM_ttA_bin_47_55/ttA_bin_47_55_madgraph_LO);
            }
        }
    } else if(b_min == 55 && b_max == 70){
        
    
        double SM_ttA_bin_55_70 = SM.getOptionalParameter("SM_sigma_tta_bin_55_70_ATLAS_emu");
        
        double ttA_bin_55_70_madgraph_NLO=0.357037;//fb
        double ttA_bin_55_70_madgraph_LO = 0.26308;//fb


        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return   SM_ttA_bin_55_70 + (-0.014689*C_tZ + 0.030385*C_tZ*C_tZ + 0.018715*C_tW + 0.037912*C_tW*C_tW - 0.07414814815*C_tZ*C_tW + 0.08728395062*C_tG +0.01195555556*C_tG*C_tG)*(SM_ttA_bin_55_70/ttA_bin_55_70_madgraph_NLO)
                        + (0.0063986*C_tu8 + 0.0041057*C_td8 + 0.0024281*C_tu8*C_tu8 + 0.0016332*C_td8*C_td8 + 0.012520*C_Qq18 + 0.017057*C_tq8 + 0.003877*C_Qq38 + 0.014788*C_Qu8 + 0.002307*C_Qd8)*(SM_ttA_bin_55_70/ttA_bin_55_70_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_55_70 + (-0.014689*C_tZ + 0.018715*C_tW + 0.08728395062*C_tG)*(SM_ttA_bin_55_70/ttA_bin_55_70_madgraph_NLO)
                        + (0.0063986*C_tu8 + 0.0041057*C_td8 + 0.012520*C_Qq18 + 0.017057*C_tq8 + 0.003877*C_Qq38 + 0.014788*C_Qu8 + 0.002307*C_Qd8)*(SM_ttA_bin_55_70/ttA_bin_55_70_madgraph_LO);
            }
        }
    
        else{
            if(flag_Quadratic){
                return   SM_ttA_bin_55_70 + (-0.014689*(-0.472123*C_tB + 0.881533*C_tW) + 0.030385*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 0.018715*C_tW + 0.037912*C_tW*C_tW - 0.07414814815*(-0.472123*C_tB + 0.881533*C_tW)*C_tW + 0.08728395062*C_tG +0.01195555556*C_tG*C_tG )*(SM_ttA_bin_55_70/ttA_bin_55_70_madgraph_NLO)
                        + (0.0063986*C_tu8 + 0.0041057*C_td8 + 0.0024281*C_tu8*C_tu8 + 0.0016332*C_td8*C_td8)*(SM_ttA_bin_55_70/ttA_bin_55_70_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_55_70 + (-0.014689*(-0.472123*C_tB + 0.881533*C_tW) + 0.018715*C_tW + 0.08728395062*C_tG)*(SM_ttA_bin_55_70/ttA_bin_55_70_madgraph_NLO)
                        + (0.0063986*C_tu8 + 0.0041057*C_td8 )*(SM_ttA_bin_55_70/ttA_bin_55_70_madgraph_LO);
            }
        }
        
    } else if(b_min == 70 && b_max == 85){
        
    
        double SM_ttA_bin_70_85 = SM.getOptionalParameter("SM_sigma_tta_bin_70_85_ATLAS_emu");
        
        double ttA_bin_70_85_madgraph_NLO = 0.22535;//fb
        double ttA_bin_70_85_madgraph_LO = 0.18502;//fb


        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttA_bin_70_85 + (-0.009959*C_tZ + 0.033236*C_tZ*C_tZ + 0.007004*C_tW + 0.042146*C_tW*C_tW - 0.07424674897*C_tZ*C_tW + 0.06699588477*C_tG +0.0187654321*C_tG*C_tG)*(SM_ttA_bin_70_85/ttA_bin_70_85_madgraph_NLO)
                        + (0.0046269*C_tu8 + 0.0030051*C_td8+ 0.0019743*C_tu8*C_tu8 + 0.0030051*C_td8*C_td8 + 0.008563*C_Qq18 + 0.011882*C_tq8 + 0.002455*C_Qq38 + 0.010339*C_Qu8 + 0.001584*C_Qd8)*(SM_ttA_bin_70_85/ttA_bin_70_85_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_70_85 + (-0.009959*C_tZ + 0.007004*C_tW + 0.06699588477*C_tG)*(SM_ttA_bin_70_85/ttA_bin_70_85_madgraph_NLO)
                     + (0.0046269*C_tu8 + 0.0030051*C_td8 + 0.008563*C_Qq18 + 0.011882*C_tq8 + 0.002455*C_Qq38 + 0.010339*C_Qu8 + 0.001584*C_Qd8)*(SM_ttA_bin_70_85/ttA_bin_70_85_madgraph_LO);
            }
        }
    
        else{
            if(flag_Quadratic){
                return  SM_ttA_bin_70_85 + (-0.009959*(-0.472123*C_tB + 0.881533*C_tW) + 0.033236*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 0.007004*C_tW + 0.042146*C_tW*C_tW - 0.07424674897*(-0.472123*C_tB + 0.881533*C_tW)*C_tW + 0.06699588477*C_tG +0.0187654321*C_tG*C_tG )*(SM_ttA_bin_70_85/ttA_bin_70_85_madgraph_NLO)
                        + (0.0046269*C_tu8 + 0.0030051*C_td8 + 0.0019743*C_tu8*C_tu8 + 0.0030051*C_td8*C_td8)*(SM_ttA_bin_70_85/ttA_bin_70_85_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_70_85 + (-0.009959*(-0.472123*C_tB + 0.881533*C_tW) + 0.007004*C_tW + 0.06699588477*C_tG )*(SM_ttA_bin_70_85/ttA_bin_70_85_madgraph_NLO)
                        + (0.0046269*C_tu8 + 0.0030051*C_td8)*(SM_ttA_bin_70_85/ttA_bin_70_85_madgraph_LO);
            }
        }
    
    } else if(b_min == 85 && b_max == 132){
        
    
        double SM_ttA_bin_85_132 = SM.getOptionalParameter("SM_sigma_tta_bin_85_132_ATLAS_emu");
        
        double ttA_bin_85_132_madgraph_NLO=0.132598;//fb
        double ttA_bin_85_132_madgraph_LO = 0.10484;//fb


        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return SM_ttA_bin_85_132 + (-0.011540*C_tZ + 0.027391*C_tZ*C_tZ + 0.013199*C_tW + 0.036253*C_tW*C_tW - 0.0629676911*C_tZ*C_tW 	+ 0.03607827686*C_tG +0.007558707644*C_tG*C_tG)*(SM_ttA_bin_85_132/ttA_bin_85_132_madgraph_NLO)
                        + (0.0026823*C_tu8 + 0.0016437*C_td8 + 0.0011696*C_tu8*C_tu8 + 0.0008383*C_td8*C_td8 + 0.004705*C_Qq18 + 0.006739*C_tq8 + 0.001183*C_Qq38 + 0.005869*C_Qu8 + 0.000881*C_Qd8)*(SM_ttA_bin_85_132/ttA_bin_85_132_madgraph_LO);
            }
            else{
                return SM_ttA_bin_85_132 + (-0.011540*C_tZ  + 0.013199*C_tW + 0.03607827686*C_tG)*(SM_ttA_bin_85_132/ttA_bin_85_132_madgraph_NLO)
                        + (0.0026823*C_tu8 + 0.0016437*C_td8 + 0.004705*C_Qq18 + 0.006739*C_tq8 + 0.001183*C_Qq38 + 0.005869*C_Qu8 + 0.000881*C_Qd8)*(SM_ttA_bin_85_132/ttA_bin_85_132_madgraph_LO);
            }
        }
    
        else{
            if(flag_Quadratic){
                return SM_ttA_bin_85_132 + (-0.011540*(-0.472123*C_tB + 0.881533*C_tW) + 0.027391*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 0.013199*C_tW + 0.036253*C_tW*C_tW - 0.0629676911*(-0.472123*C_tB + 0.881533*C_tW)*C_tW 	+ 0.03607827686*C_tG +0.007558707644*C_tG*C_tG )*(SM_ttA_bin_85_132/ttA_bin_85_132_madgraph_NLO)
                        + (0.0026823*C_tu8 + 0.0016437*C_td8 + 0.0011696*C_tu8*C_tu8 + 0.0008383*C_td8*C_td8)*(SM_ttA_bin_85_132/ttA_bin_85_132_madgraph_LO);
            }
            else{
                return SM_ttA_bin_85_132 + (-0.011540*(-0.472123*C_tB + 0.881533*C_tW)  + 0.013199*C_tW + 0.03607827686*C_tG)*(SM_ttA_bin_85_132/ttA_bin_85_132_madgraph_NLO)
                        + (0.0026823*C_tu8 + 0.0016437*C_td8)*(SM_ttA_bin_85_132/ttA_bin_85_132_madgraph_LO);
            }
        }
        
    } else if(b_min == 132 && b_max == 180){
        
    
        double SM_ttA_bin_132_180 = SM.getOptionalParameter("SM_sigma_tta_bin_132_180_ATLAS_emu");
        
        double ttA_bin_132_180_madgraph_NLO = 0.0582305;//fb
        double ttA_bin_132_180_madgraph_LO = 0.04850;//fb

    
        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttA_bin_132_180 + (-0.005155*C_tZ + 0.021199*C_tZ*C_tZ + 0.006004*C_tW + 0.029633*C_tW*C_tW - 0.05011851852*C_tZ*C_tW + 0.01400462963*C_tG +0.003517489712*C_tG*C_tG )*(SM_ttA_bin_132_180/ttA_bin_132_180_madgraph_NLO)
                        + (0.0011640*C_tu8 + 0.0008132*C_td8 + 0.0006643*C_tu8*C_tu8 + 0.0004577*C_td8*C_td8 + 0.002180*C_Qq18 + 0.003232*C_tq8 + 0.000454*C_Qq38 + 0.002820*C_Qu8 + 0.000416*C_Qd8)*(SM_ttA_bin_132_180/ttA_bin_132_180_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_132_180 + (-0.005155*C_tZ + 0.006004*C_tW  + 0.01400462963*C_tG)*(SM_ttA_bin_132_180/ttA_bin_132_180_madgraph_NLO)
                        + (0.0011640*C_tu8 + 0.0008132*C_td8 + 0.002180*C_Qq18 + 0.003232*C_tq8 + 0.000454*C_Qq38 + 0.002820*C_Qu8 + 0.000416*C_Qd8)*(SM_ttA_bin_132_180/ttA_bin_132_180_madgraph_LO);
            }
        }
        else{
            if(flag_Quadratic){
                return  SM_ttA_bin_132_180 + (-0.005155*(-0.472123*C_tB + 0.881533*C_tW) + 0.021199*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 0.006004*C_tW + 0.029633*C_tW*C_tW - 0.05011851852*(-0.472123*C_tB + 0.881533*C_tW)*C_tW + 0.01400462963*C_tG +0.003517489712*C_tG*C_tG )*(SM_ttA_bin_132_180/ttA_bin_132_180_madgraph_NLO)
                        + (0.0011640*C_tu8 + 0.0008132*C_td8 + 0.0006643*C_tu8*C_tu8 + 0.0004577*C_td8*C_td8)*(SM_ttA_bin_132_180/ttA_bin_132_180_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_132_180 + (-0.005155*(-0.472123*C_tB + 0.881533*C_tW) + 0.006004*C_tW  + 0.01400462963*C_tG )*(SM_ttA_bin_132_180/ttA_bin_132_180_madgraph_NLO)
                        + (0.0011640*C_tu8 + 0.0008132*C_td8)*(SM_ttA_bin_132_180/ttA_bin_132_180_madgraph_LO);
            }
        }
    } else if(b_min == 180 && b_max == 300){
        
    
        double SM_ttA_bin_180_300 = SM.getOptionalParameter("SM_sigma_tta_bin_180_300_ATLAS_emu");
        
        double ttA_bin_180_300_madgraph_NLO = 0.0256399;//fb
        double ttA_bin_180_300_madgraph_LO = 0.01676;//fb

    
        if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttA_bin_180_300 + (-0.001690*C_tZ + 0.024677*C_tZ*C_tZ + 0.001959*C_tW + 0.031588*C_tW*C_tW - 0.05586666667*C_tZ*C_tW + 0.006676954733*C_tG + 0.002124485597*C_tG*C_tG)*(SM_ttA_bin_180_300/ttA_bin_180_300_madgraph_NLO)
                        + (0.0004632*C_tu8 + 0.0003355*C_td8+ 0.0003386*C_tu8*C_tu8 + 0.0002343*C_td8*C_td8 + 0.000836*C_Qq18 + 0.001286*C_tq8 + 0.000143*C_Qq38 + 0.001123*C_Qu8 + 0.000163*C_Qd8)*(SM_ttA_bin_180_300/ttA_bin_180_300_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_180_300 + (-0.001690*C_tZ  + 0.001959*C_tW + 0.006676954733*C_tG)*(SM_ttA_bin_180_300/ttA_bin_180_300_madgraph_NLO)
                        + (0.0004632*C_tu8 + 0.0003355*C_td8 + 0.000836*C_Qq18 + 0.001286*C_tq8 + 0.000143*C_Qq38 + 0.001123*C_Qu8 + 0.000163*C_Qd8)*(SM_ttA_bin_180_300/ttA_bin_180_300_madgraph_LO);
            }
        }
    
        else{
            if(flag_Quadratic){
                return  SM_ttA_bin_180_300 + (-0.001690*(-0.472123*C_tB + 0.881533*C_tW) + 0.024677*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 0.001959*C_tW + 0.031588*C_tW*C_tW - 0.05586666667*(-0.472123*C_tB + 0.881533*C_tW)*C_tW + 0.006676954733*C_tG + 0.002124485597*C_tG*C_tG )*(SM_ttA_bin_180_300/ttA_bin_180_300_madgraph_NLO)
                        + (0.0004632*C_tu8 + 0.0003355*C_td8 + 0.0003386*C_tu8*C_tu8 + 0.0002343*C_td8*C_td8)*(SM_ttA_bin_180_300/ttA_bin_180_300_madgraph_LO);
            }
            else{
                return  SM_ttA_bin_180_300 + (-0.001690*(-0.472123*C_tB + 0.881533*C_tW)  + 0.001959*C_tW + 0.006676954733*C_tG)*(SM_ttA_bin_180_300/ttA_bin_180_300_madgraph_NLO)
                        + (0.0004632*C_tu8 + 0.0003355*C_td8)*(SM_ttA_bin_180_300/ttA_bin_180_300_madgraph_LO);
            }
        }
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct binning range for sigma_tta_diff_LO_ATLAS_emu.\n");
    }
    
}








//REAL NLO
sigma_tta_diff_NLO_CMS_dilepton_220107301::sigma_tta_diff_NLO_CMS_dilepton_220107301(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tta_bin_20_35_CMS_dilepton" << "SM_sigma_tta_bin_35_50_CMS_dilepton"
            << "SM_sigma_tta_bin_50_70_CMS_dilepton" << "SM_sigma_tta_bin_70_130_CMS_dilepton" 
            << "SM_sigma_tta_bin_130_200_CMS_dilepton" << "SM_sigma_tta_bin_200_300_CMS_dilepton");

    
}

double sigma_tta_diff_NLO_CMS_dilepton_220107301::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    

    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    
    

    if(b_min == 20 && b_max == 35){
        
        double SM_sigma_tta_bin_20_35 = SM.getOptionalParameter("SM_sigma_tta_bin_20_35_CMS_dilepton");
        double sigma_tta_bin_20_35_madgraph = 0.975113; //pb
        
        
        if(flag_Quadratic){
        
            return SM_sigma_tta_bin_20_35 +(0.975113+0.000788*C_Qd1+0.022216*C_Qd1*C_Qd1
                    +0.015615*C_Qd8+0.009971*C_Qd8*C_Qd8+0.003095*C_Qq11+0.074636*C_Qq11*C_Qq11
                    +-0.000392*C_Qq31+0.076939*C_Qq31*C_Qq31+0.025252*C_Qq18+0.029486*C_Qq18*C_Qq18
                    +0.009971*C_Qq38+0.032206*C_Qq38*C_Qq38+0.004873*C_Qu1+0.071282*C_Qu1*C_Qu1
                    +0.027259*C_Qu8+0.022655*C_Qu8*C_Qu8+0.010583*C_td1+0.021511*C_td1*C_td1
                    +0.008297*C_td8+0.022696*C_td8*C_td8+0.293375*C_tG+0.071832*C_tG*C_tG
                    -0.017083*C_tq1+0.069788*C_tq1*C_tq1+0.028813*C_tq8+0.031600*C_tq8*C_tq8
                    -0.000271*C_tu1+0.055283*C_tu1*C_tu1+0.015878*C_tu8+0.028156*C_tu8*C_tu8
                    +0.015359*C_tW+0.021535*C_tW*C_tW+-0.020018*C_tZ+0.030927*C_tZ*C_tZ)
                    *(SM_sigma_tta_bin_20_35/sigma_tta_bin_20_35_madgraph);
            
        }
        else{
            
            
            
            return SM_sigma_tta_bin_20_35 +(0.000788*C_Qd1+0.015615*C_Qd8
                    +0.003095*C_Qq11+-0.000392*C_Qq31+0.025252*C_Qq18+0.009971*C_Qq38
                    +0.004873*C_Qu1+0.027259*C_Qu8+0.010583*C_td1+0.008297*C_td8
                    +0.293375*C_tG+-0.017083*C_tq1+0.028813*C_tq8+-0.000271*C_tu1
                    +0.015878*C_tu8+0.015359*C_tW+-0.020018*C_tZ)
                    *(SM_sigma_tta_bin_20_35/sigma_tta_bin_20_35_madgraph);
        }
    } else if(b_min == 35 && b_max == 50){
        
    
        double SM_sigma_tta_bin_35_50 = SM.getOptionalParameter("SM_sigma_tta_bin_35_50_CMS_dilepton");
        double sigma_tta_bin_35_50_madgraph = 0.523831; //pb
        
        
        if(flag_Quadratic){
        
            return  SM_sigma_tta_bin_35_50 +(-0.002203*C_Qd1+0.002912*C_Qd1*C_Qd1+0.001037*C_Qd8
                    -0.001324*C_Qd8*C_Qd8+-0.001593*C_Qq11+0.021815*C_Qq11*C_Qq11+-0.002413*C_Qq31
                    +0.014750*C_Qq31*C_Qq31+0.023210*C_Qq18+-0.000332*C_Qq18*C_Qq18+0.003673*C_Qq38
                    +0.001051*C_Qq38*C_Qq38+0.000986*C_Qu1+0.033309*C_Qu1*C_Qu1+0.012209*C_Qu8
                    +0.000396*C_Qu8*C_Qu8+0.005764*C_td1+0.007074*C_td1*C_td1+0.002688*C_td8
                    -0.003340*C_td8*C_td8+0.147862*C_tG+0.023437*C_tG*C_tG+-0.004221*C_tq1
                    +0.039220*C_tq1*C_tq1+0.019842*C_tq8+-0.005781*C_tq8*C_tq8+-0.000079*C_tu1
                    +0.010630*C_tu1*C_tu1+0.012863*C_tu8+-0.007084*C_tu8*C_tu8+0.019122*C_tW
                    +0.025787*C_tW*C_tW+-0.021752*C_tZ+0.017695*C_tZ*C_tZ)
                    *(SM_sigma_tta_bin_35_50/sigma_tta_bin_35_50_madgraph);;
            
        }
        else{
            
            return SM_sigma_tta_bin_35_50 +(-0.002203*C_Qd1+0.001037*C_Qd8+-0.001593*C_Qq11
                    -0.002413*C_Qq31+0.023210*C_Qq18+0.003673*C_Qq38+0.000986*C_Qu1+0.012209*C_Qu8
                    +0.005764*C_td1+0.002688*C_td8+0.147862*C_tG+-0.004221*C_tq1+0.019842*C_tq8
                    -0.000079*C_tu1+0.012863*C_tu8+0.019122*C_tW+-0.021752*C_tZ)
                    *(SM_sigma_tta_bin_35_50/sigma_tta_bin_35_50_madgraph);
            
        }
    } else if(b_min == 50 && b_max == 70){
        
    
        double SM_sigma_tta_bin_50_70 = SM.getOptionalParameter("SM_sigma_tta_bin_50_70_CMS_dilepton");
        double sigma_tta_bin_50_70_madgraph = 0.399663; //pb
        
        
        if(flag_Quadratic){
        
            return  SM_sigma_tta_bin_50_70 +(-0.000900*C_Qd1+0.006167*C_Qd1*C_Qd1+0.003927*C_Qd8
                    +0.001326*C_Qd8*C_Qd8+0.001956*C_Qq11+0.021969*C_Qq11*C_Qq11+-0.002388*C_Qq31
                    +0.026569*C_Qq31*C_Qq31+0.009836*C_Qq18+0.005966*C_Qq18*C_Qq18+0.001784*C_Qq38
                    +0.004714*C_Qq38*C_Qq38+-0.008926*C_Qu1+0.026436*C_Qu1*C_Qu1+0.008793*C_Qu8
                    +0.002637*C_Qu8*C_Qu8+-0.004573*C_td1+0.006145*C_td1*C_td1+0.004629*C_td8
                    +0.004192*C_td8*C_td8+0.111319*C_tG+0.023254*C_tG*C_tG+-0.003314*C_tq1
                    +0.030819*C_tq1*C_tq1+0.010340*C_tq8+0.009217*C_tq8*C_tq8+-0.001303*C_tu1
                    +0.018581*C_tu1*C_tu1+0.004583*C_tu8+0.004069*C_tu8*C_tu8+0.024243*C_tW
                    +0.042926*C_tW*C_tW+-0.013215*C_tZ+0.030797*C_tZ*C_tZ)
                    *(SM_sigma_tta_bin_50_70/sigma_tta_bin_50_70_madgraph);;
            
        }
        else{
            
            return SM_sigma_tta_bin_50_70 +(-0.000900*C_Qd1+0.003927*C_Qd8+0.001956*C_Qq11+-0.002388*C_Qq31
                    +0.009836*C_Qq18+0.001784*C_Qq38+-0.008926*C_Qu1+0.008793*C_Qu8+-0.004573*C_td1
                    +0.004629*C_td8+0.111319*C_tG+-0.003314*C_tq1+0.010340*C_tq8+-0.001303*C_tu1
                    +0.004583*C_tu8+0.024243*C_tW+-0.013215*C_tZ)
                    *(SM_sigma_tta_bin_50_70/sigma_tta_bin_50_70_madgraph);
            
        }
    } else if(b_min == 70 && b_max == 130){
        
    
        double SM_sigma_tta_bin_70_130 = SM.getOptionalParameter("SM_sigma_tta_bin_70_130_CMS_dilepton");
        double sigma_tta_bin_70_100_madgraph = 0.330761; //pb
        double sigma_tta_bin_100_130_madgraph =0.182187; //pb
        double sigma_tta_bin_70_130_madgraph =sigma_tta_bin_70_100_madgraph+sigma_tta_bin_100_130_madgraph; //pb
        
        if(flag_Quadratic){
        
            double sigma_tta_bin_70_100_NP = (0.000736*C_Qd1+0.002661*C_Qd1*C_Qd1+0.000607*C_Qd8+0.001692*C_Qd8*C_Qd8
            +-0.000097*C_Qq11+0.020897*C_Qq11*C_Qq11+-0.001074*C_Qq31+0.019638*C_Qq31*C_Qq31+0.007285*C_Qq18
            +0.003571*C_Qq18*C_Qq18+0.001987*C_Qq38+0.001518*C_Qq38*C_Qq38+-0.003159*C_Qu1+0.027137*C_Qu1*C_Qu1
            +0.008729*C_Qu8+0.003259*C_Qu8*C_Qu8+0.002396*C_td1+0.005831*C_td1*C_td1+0.002438*C_td8
            +0.001050*C_td8*C_td8+0.088372*C_tG+0.018412*C_tG*C_tG+-0.002403*C_tq1+0.033166*C_tq1*C_tq1
            +0.010186*C_tq8+0.002150*C_tq8*C_tq8+0.001916*C_tu1+0.010698*C_tu1*C_tu1+0.007129*C_tu8
            +0.002852*C_tu8*C_tu8+0.028690*C_tW+0.055463*C_tW*C_tW+-0.024168*C_tZ+0.044158*C_tZ*C_tZ);
             
            
            double sigma_tta_bin_100_130_NP = (-0.000532*C_Qd1+0.003711*C_Qd1*C_Qd1+-0.000866*C_Qd8+0.000827*C_Qd8*C_Qd8
            +0.000833*C_Qq11+0.014119*C_Qq11*C_Qq11+0.000026*C_Qq31+0.013016*C_Qq31*C_Qq31+0.003946*C_Qq18
            +0.002769*C_Qq18*C_Qq18+0.002527*C_Qq38+0.003527*C_Qq38*C_Qq38+-0.001730*C_Qu1+0.016375*C_Qu1*C_Qu1
            +0.004798*C_Qu8+0.003116*C_Qu8*C_Qu8+-0.000105*C_td1+0.006509*C_td1*C_td1+0.001677*C_td8
            -0.000010*C_td8*C_td8+0.047396*C_tG+0.012131*C_tG*C_tG+-0.001607*C_tq1+0.017941*C_tq1*C_tq1
            +0.003987*C_tq8+0.003210*C_tq8*C_tq8+0.000006*C_tu1+0.009547*C_tu1*C_tu1+0.002691*C_tu8
            +0.001113*C_tu8*C_tu8+0.017670*C_tW+0.054483*C_tW*C_tW+-0.014593*C_tZ+0.042939*C_tZ*C_tZ);
            
            return SM_sigma_tta_bin_70_130 + (sigma_tta_bin_70_100_NP + sigma_tta_bin_100_130_NP)*(SM_sigma_tta_bin_70_130/sigma_tta_bin_70_130_madgraph);
            
            
        }
        else{
            
            double sigma_tta_bin_70_100_NP = (0.000736*C_Qd1+0.000607*C_Qd8+-0.000097*C_Qq11+-0.001074*C_Qq31
            +0.007285*C_Qq18+0.001987*C_Qq38+-0.003159*C_Qu1+0.008729*C_Qu8+0.002396*C_td1+0.002438*C_td8
            +0.088372*C_tG+-0.002403*C_tq1+0.010186*C_tq8+0.001916*C_tu1+0.007129*C_tu8+0.028690*C_tW-0.024168*C_tZ);
             
            
            double sigma_tta_bin_100_130_NP = (-0.000532*C_Qd1+-0.000866*C_Qd8+0.000833*C_Qq11+0.000026*C_Qq31
            +0.003946*C_Qq18+0.002527*C_Qq38+-0.001730*C_Qu1+0.004798*C_Qu8+-0.000105*C_td1+0.001677*C_td8
            +0.047396*C_tG+-0.001607*C_tq1+0.003987*C_tq8+0.000006*C_tu1+0.002691*C_tu8+0.017670*C_tW+-0.014593*C_tZ);
            
            return SM_sigma_tta_bin_70_130 + (sigma_tta_bin_70_100_NP + sigma_tta_bin_100_130_NP)*(SM_sigma_tta_bin_70_130/sigma_tta_bin_70_130_madgraph);
            
        }
        
    } else if(b_min == 130 && b_max == 200){
        
    
        double SM_sigma_tta_bin_130_200 = SM.getOptionalParameter("SM_sigma_tta_bin_130_200_CMS_dilepton");
        double sigma_tta_bin_130_165_madgraph = 0.123480; //pb
        double sigma_tta_bin_165_200_madgraph = 0.071661; //pb
        double sigma_tta_bin_130_200_madgraph = sigma_tta_bin_130_165_madgraph+sigma_tta_bin_165_200_madgraph;
        
        if(flag_Quadratic){
            
            
            double sigma_tta_bin_130_165_NP = (0.123480+-0.000179*C_Qd1+0.001106*C_Qd1*C_Qd1+0.000840*C_Qd8
            +-0.000657*C_Qd8*C_Qd8+0.000217*C_Qq11+0.009123*C_Qq11*C_Qq11+-0.000610*C_Qq31+0.009951*C_Qq31*C_Qq31
            +0.002961*C_Qq18+0.001668*C_Qq18*C_Qq18+0.000933*C_Qq38+0.001121*C_Qq38*C_Qq38+-0.001176*C_Qu1
            +0.012505*C_Qu1*C_Qu1+0.003349*C_Qu8+0.001583*C_Qu8*C_Qu8+0.001208*C_td1+0.003861*C_td1*C_td1
            +0.001677*C_td8+0.001007*C_td8*C_td8+0.033068*C_tG+0.007722*C_tG*C_tG+-0.001873*C_tq1
            +0.014593*C_tq1*C_tq1+0.004056*C_tq8+0.001186*C_tq8*C_tq8+-0.000148*C_tu1+0.005718*C_tu1*C_tu1
            +0.000955*C_tu8+0.001565*C_tu8*C_tu8+0.008850*C_tW+0.054838*C_tW*C_tW+-0.011773*C_tZ+0.042418*C_tZ*C_tZ);

                                
            double sigma_tta_bin_165_200_NP = (0.071661+0.000321*C_Qd1+0.001864*C_Qd1*C_Qd1+0.000523*C_Qd8
            +0.000484*C_Qd8*C_Qd8+-0.000468*C_Qq11+0.006489*C_Qq11*C_Qq11+-0.000323*C_Qq31+0.007170*C_Qq31*C_Qq31
            +0.001760*C_Qq18+0.001505*C_Qq18*C_Qq18+0.000167*C_Qq38+0.001674*C_Qq38*C_Qq38+-0.001324*C_Qu1
            +0.009641*C_Qu1*C_Qu1+0.001312*C_Qu8+0.000612*C_Qu8*C_Qu8+0.000060*C_td1+0.002828*C_td1*C_td1
            +0.000787*C_td8+0.000566*C_td8*C_td8+0.018680*C_tG+0.005598*C_tG*C_tG+-0.001764*C_tq1
            +0.009979*C_tq1*C_tq1+0.001779*C_tq8+0.001441*C_tq8*C_tq8+-0.000282*C_tu1+0.004307*C_tu1*C_tu1
            +0.001148*C_tu8+0.000369*C_tu8*C_tu8+0.008066*C_tW+0.045701*C_tW*C_tW+-0.004856*C_tZ+0.033619*C_tZ*C_tZ);
            
            return SM_sigma_tta_bin_130_200 + (sigma_tta_bin_130_165_NP+sigma_tta_bin_165_200_NP)*(SM_sigma_tta_bin_130_200/sigma_tta_bin_130_200_madgraph);
            
            
        }
        else{
            
            
            double sigma_tta_bin_130_165_NP = (-0.000179*C_Qd1+0.000840*C_Qd8+0.000217*C_Qq11+-0.000610*C_Qq31
            +0.002961*C_Qq18+0.000933*C_Qq38+-0.001176*C_Qu1+0.003349*C_Qu8+0.001208*C_td1+0.001677*C_td8
            +0.033068*C_tG+-0.001873*C_tq1+0.004056*C_tq8+-0.000148*C_tu1+0.000955*C_tu8+0.008850*C_tW+-0.011773*C_tZ);


            double sigma_tta_bin_165_200_NP = (0.000321*C_Qd1+0.000523*C_Qd8+-0.000468*C_Qq11+-0.000323*C_Qq31
            +0.001760*C_Qq18+0.000167*C_Qq38+-0.001324*C_Qu1+0.001312*C_Qu8+0.000060*C_td1+0.000787*C_td8+0.018680*C_tG
            -0.001764*C_tq1+0.001779*C_tq8+-0.000282*C_tu1+0.001148*C_tu8+0.008066*C_tW+-0.004856*C_tZ);
            
            return SM_sigma_tta_bin_130_200 + (sigma_tta_bin_130_165_NP+sigma_tta_bin_165_200_NP)*(SM_sigma_tta_bin_130_200/sigma_tta_bin_130_200_madgraph);
            
        }
    } else if(b_min == 200 && b_max == 300){
        
    
        double SM_sigma_tta_bin_200_300 = SM.getOptionalParameter("SM_sigma_tta_bin_200_300_CMS_dilepton");
        double sigma_tta_bin_200_250_madgraph = 0.057272; //pb
        double sigma_tta_bin_250_300_madgraph = 0.067466; //pb
        double sigma_tta_bin_200_300_madgraph = sigma_tta_bin_200_250_madgraph+sigma_tta_bin_250_300_madgraph; //pb
        
        
        if(flag_Quadratic){
        
            double sigma_tta_bin_200_250_NP = (-0.000047*C_Qd1+0.001342*C_Qd1*C_Qd1+0.000136*C_Qd8+0.000006*C_Qd8*C_Qd8
            +0.000342*C_Qq11+0.006388*C_Qq11*C_Qq11+-0.000150*C_Qq31+0.006218*C_Qq31*C_Qq31+0.001724*C_Qq18
            +0.000750*C_Qq18*C_Qq18+0.000025*C_Qq38+0.000623*C_Qq38*C_Qq38+-0.000462*C_Qu1+0.007979*C_Qu1*C_Qu1
            +0.001728*C_Qu8+0.001000*C_Qu8*C_Qu8+0.000385*C_td1+0.001911*C_td1*C_td1+0.000618*C_td8+0.000138*C_td8*C_td8
            +0.014921*C_tG+0.004175*C_tG*C_tG+-0.000474*C_tq1+0.009899*C_tq1*C_tq1+0.002069*C_tq8+0.000524*C_tq8*C_tq8
            +0.000373*C_tu1+0.003594*C_tu1*C_tu1+0.000824*C_tu8+0.000624*C_tu8*C_tu8+0.003702*C_tW+0.049912*C_tW*C_tW
            +-0.004612*C_tZ+0.038502*C_tZ*C_tZ);
            
            double sigma_tta_bin_250_300_NP = (-0.000159*C_Qd1+0.004318*C_Qd1*C_Qd1+0.000593*C_Qd8+0.000437*C_Qd8*C_Qd8
            +0.000507*C_Qq11+0.019093*C_Qq11*C_Qq11+-0.000174*C_Qq31+0.019040*C_Qq31*C_Qq31+0.002418*C_Qq18
            +0.003576*C_Qq18*C_Qq18+0.000138*C_Qq38+0.003125*C_Qq38*C_Qq38+-0.001614*C_Qu1+0.023126*C_Qu1*C_Qu1
            +0.002333*C_Qu8+0.002169*C_Qu8*C_Qu8+0.000208*C_td1+0.006881*C_td1*C_td1+0.001121*C_td8+0.001259*C_td8*C_td8
            +0.017830*C_tG+0.007512*C_tG*C_tG+-0.001787*C_tq1+0.027531*C_tq1*C_tq1+0.002581*C_tq8+0.002984*C_tq8*C_tq8
            -0.000132*C_tu1+0.012280*C_tu1*C_tu1+0.001367*C_tu8+0.002013*C_tu8*C_tu8+0.003271*C_tW+0.124997*C_tW*C_tW
            -0.002154*C_tZ+0.096300*C_tZ*C_tZ);
            
            return SM_sigma_tta_bin_200_300 + (sigma_tta_bin_200_250_NP + sigma_tta_bin_250_300_NP)*(SM_sigma_tta_bin_200_300/sigma_tta_bin_200_300_madgraph);
            
        }
        else{
            
            double sigma_tta_bin_200_250_NP = (-0.000047*C_Qd1+0.000136*C_Qd8+0.000342*C_Qq11+-0.000150*C_Qq31
            +0.001724*C_Qq18+0.000025*C_Qq38+-0.000462*C_Qu1+0.001728*C_Qu8+0.000385*C_td1+0.000618*C_td8
            +0.014921*C_tG+-0.000474*C_tq1+0.002069*C_tq8+0.000373*C_tu1+0.000824*C_tu8+0.003702*C_tW+-0.004612*C_tZ);
            
            double sigma_tta_bin_250_300_NP = (-0.000159*C_Qd1+0.000593*C_Qd8+0.000507*C_Qq11+-0.000174*C_Qq31
            +0.002418*C_Qq18+0.000138*C_Qq38+-0.001614*C_Qu1+0.002333*C_Qu8+0.000208*C_td1+0.001121*C_td8
            +0.017830*C_tG+-0.001787*C_tq1+0.002581*C_tq8+-0.000132*C_tu1+0.001367*C_tu8+0.003271*C_tW+-0.002154*C_tZ);
            
            return SM_sigma_tta_bin_200_300 + (sigma_tta_bin_200_250_NP + sigma_tta_bin_250_300_NP)*(SM_sigma_tta_bin_200_300/sigma_tta_bin_200_300_madgraph);
            
        }
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct binning range for sigma_tta_diff_LO_CMS_dilepton.\n");
    }

}










sigma_tta_diff_NLO_ATLAS_240309452::sigma_tta_diff_NLO_ATLAS_240309452(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tta_bin_20_40_ATLAS_240309452" 
            << "SM_sigma_tta_bin_40_60_ATLAS_240309452" << "SM_sigma_tta_bin_60_75_ATLAS_240309452" 
            << "SM_sigma_tta_bin_75_95_ATLAS_240309452" << "SM_sigma_tta_bin_95_115_ATLAS_240309452" 
            << "SM_sigma_tta_bin_115_135_ATLAS_240309452" << "SM_sigma_tta_bin_135_160_ATLAS_240309452"
            << "SM_sigma_tta_bin_160_200_ATLAS_240309452" << "SM_sigma_tta_bin_200_260_ATLAS_240309452"
             << "SM_sigma_tta_bin_260_500_ATLAS_240309452");

    
}

double sigma_tta_diff_NLO_ATLAS_240309452::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    
    
    
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    
    

    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    
    

    if(b_min == 20 && b_max == 40){
        
        double SM_sigma_tta_bin_20_40 = SM.getOptionalParameter("SM_sigma_tta_bin_20_40_ATLAS_240309452");
        double sigma_tta_bin_20_40_madgraph = 1177.4069; //fb
        
        
        if(flag_Quadratic){
        
            return SM_sigma_tta_bin_20_40 +(-3.6263*C_phiQ3 + 0.0339*C_phiQ3*C_phiQ3 + -3.6643*C_phiQm +
                    0.0227*C_phiQm*C_phiQm + -3.5919*C_phit + -0.0384*C_phit*C_phit + 12.2447*C_Qd8 + 
                    0.7109*C_Qd8*C_Qd8 + 32.9687*C_Qq18 + 11.2487*C_Qq18*C_Qq18 + 12.8705*C_Qq38 + 
                    10.4471*C_Qq38*C_Qq38 + 32.3790*C_Qu8 + 9.3811*C_Qu8*C_Qu8 + 9.7399*C_td8 + 
                    2.1104*C_td8*C_td8 + 334.4265*C_tG + 61.2791*C_tG*C_tG + -3.3893*C_tphi + 
                    -0.3108*C_tphi*C_tphi + 36.1056*C_tq8 + 10.0010*C_tq8*C_tq8 + 23.6950*C_tu8 + 
                    9.0848*C_tu8*C_tu8 + 34.7636*C_tW + 29.3066*C_tW*C_tW + -31.0740*C_tZ + 29.0832*C_tZ*C_tZ)
                    *(SM_sigma_tta_bin_20_40/sigma_tta_bin_20_40_madgraph);
            
        }
        else{
            
            
            
            return SM_sigma_tta_bin_20_40 +(-3.6263*C_phiQ3 + -3.6643*C_phiQm + -3.5919*C_phit + 12.2447*C_Qd8 +
                    32.9687*C_Qq18 + 12.8705*C_Qq38 + 32.3790*C_Qu8 + 9.7399*C_td8 + 334.4265*C_tG + 
                    -3.3893*C_tphi + 36.1056*C_tq8 + 23.6950*C_tu8 + 34.7636*C_tW + -31.0740*C_tZ)
                    *(SM_sigma_tta_bin_20_40/sigma_tta_bin_20_40_madgraph);
        }
    } else if(b_min == 40 && b_max == 60){
        
    
        double SM_sigma_tta_bin_40_60 = SM.getOptionalParameter("SM_sigma_tta_bin_40_60_ATLAS_240309452");
        double sigma_tta_bin_40_60_madgraph = 529.8128; //fb
        
        
        if(flag_Quadratic){
        
            return  SM_sigma_tta_bin_40_60 +(-1.2328*C_phiQ3 + 1.1934*C_phiQ3*C_phiQ3 + -1.4786*C_phiQm + 1.1073*C_phiQm*C_phiQm +
                    -1.0026*C_phit + 1.4763*C_phit*C_phit + 2.9225*C_Qd8 + 3.3985*C_Qd8*C_Qd8 + 12.2580*C_Qq18 + 5.8612*C_Qq18*C_Qq18 + 
                    3.7126*C_Qq38 + 5.8362*C_Qq38*C_Qq38 + 14.1556*C_Qu8 + 4.1339*C_Qu8*C_Qu8 + 7.2970*C_td8 + 2.3002*C_td8*C_td8 +
                    147.2638*C_tG + 28.7812*C_tG*C_tG + -0.6709*C_tphi + 1.0703*C_tphi*C_tphi + 16.0890*C_tq8 + 7.2375*C_tq8*C_tq8 +
                    10.6510*C_tu8 + 4.9574*C_tu8*C_tu8 + 22.6910*C_tW + 41.7841*C_tW*C_tW + -24.1873*C_tZ + 28.9994*C_tZ*C_tZ)
                    *(SM_sigma_tta_bin_40_60/sigma_tta_bin_40_60_madgraph);;
            
        }
        else{
            
            return SM_sigma_tta_bin_40_60 +(-1.2328*C_phiQ3 + -1.4786*C_phiQm + -1.0026*C_phit + 2.9225*C_Qd8 +
                    12.2580*C_Qq18 + 3.7126*C_Qq38 + 14.1556*C_Qu8 + 7.2970*C_td8 + 147.2638*C_tG + 
                    -0.6709*C_tphi + 16.0890*C_tq8 + 10.6510*C_tu8 + 22.6910*C_tW + -24.1873*C_tZ)
                    *(SM_sigma_tta_bin_40_60/sigma_tta_bin_40_60_madgraph);
            
        }
    } else if(b_min == 60 && b_max == 75){
        
    
        double SM_sigma_tta_bin_60_75 = SM.getOptionalParameter("SM_sigma_tta_bin_60_75_ATLAS_240309452");
        double sigma_tta_bin_60_75_madgraph = 241.6828; //fb
        
        
        if(flag_Quadratic){
        
            return  SM_sigma_tta_bin_60_75 +(-0.4801*C_phiQ3 + -0.3292*C_phiQ3*C_phiQ3 + -0.5938*C_phiQm + -0.4390*C_phiQm*C_phiQm + 
                    -0.6264*C_phit + -0.3886*C_phit*C_phit + 1.2599*C_Qd8 + 0.6517*C_Qd8*C_Qd8 + 6.0605*C_Qq18 + 2.1080*C_Qq18*C_Qq18 +
                    1.7701*C_Qq38 + 2.6314*C_Qq38*C_Qq38 + 4.2427*C_Qu8 + 2.7795*C_Qu8*C_Qu8 + 1.8165*C_td8 + 0.8025*C_td8*C_td8 + 
                    67.5221*C_tG + 14.7931*C_tG*C_tG + -0.6729*C_tphi + -0.4015*C_tphi*C_tphi + 7.7025*C_tq8 + 2.6132*C_tq8*C_tq8 +
                    5.9503*C_tu8 + 1.4553*C_tu8*C_tu8 + 16.0354*C_tW + 30.1549*C_tW*C_tW + -15.4837*C_tZ + 21.8673*C_tZ*C_tZ)
                    *(SM_sigma_tta_bin_60_75/sigma_tta_bin_60_75_madgraph);;
            
        }
        else{
            
            return SM_sigma_tta_bin_60_75 +(-0.4801*C_phiQ3 + -0.5938*C_phiQm + -0.6264*C_phit + 1.2599*C_Qd8 + 
                    6.0605*C_Qq18 + 1.7701*C_Qq38 + 4.2427*C_Qu8 + 1.8165*C_td8 + 67.5221*C_tG + 
                    -0.6729*C_tphi + 7.7025*C_tq8 + 5.9503*C_tu8 + 16.0354*C_tW + -15.4837*C_tZ)
                    *(SM_sigma_tta_bin_60_75/sigma_tta_bin_60_75_madgraph);
            
        }
    } else if(b_min == 75 && b_max == 95){
        
    
        double SM_sigma_tta_bin_75_95 = SM.getOptionalParameter("SM_sigma_tta_bin_75_95_ATLAS_240309452");
        double sigma_tta_bin_75_95_madgraph = 215.0879; //fb
        
        
        if(flag_Quadratic){
        
            return  SM_sigma_tta_bin_75_95 +(-0.6245*C_phiQ3 + -0.1479*C_phiQ3*C_phiQ3 + -0.3324*C_phiQm + 0.0034*C_phiQm*C_phiQm + 
                    -0.4829*C_phit + 0.0294*C_phit*C_phit + 1.6329*C_Qd8 + 0.6331*C_Qd8*C_Qd8 + 4.9233*C_Qq18 + 2.8120*C_Qq18*C_Qq18 +
                    1.1708*C_Qq38 + 3.2268*C_Qq38*C_Qq38 + 7.2407*C_Qu8 + 2.9187*C_Qu8*C_Qu8 + 3.1319*C_td8 + 1.2359*C_td8*C_td8 + 
                    58.3993*C_tG + 12.7909*C_tG*C_tG + -0.1904*C_tphi + -0.2576*C_tphi*C_tphi + 7.2798*C_tq8 + 2.6504*C_tq8*C_tq8 + 
                    3.4404*C_tu8 + 1.8945*C_tu8*C_tu8 + 16.1736*C_tW + 40.0102*C_tW*C_tW + -14.0892*C_tZ + 30.0366*C_tZ*C_tZ)
                    *(SM_sigma_tta_bin_75_95/sigma_tta_bin_75_95_madgraph);;
            
        }
        else{
            
            return SM_sigma_tta_bin_75_95 +(-0.6245*C_phiQ3 + -0.3324*C_phiQm + -0.4829*C_phit + 1.6329*C_Qd8 +
                    4.9233*C_Qq18 + 1.1708*C_Qq38 + 7.2407*C_Qu8 + 3.1319*C_td8 + 58.3993*C_tG + 
                    -0.1904*C_tphi + 7.2798*C_tq8 + 3.4404*C_tu8 + 16.1736*C_tW + -14.0892*C_tZ)
                    *(SM_sigma_tta_bin_75_95/sigma_tta_bin_75_95_madgraph);
            
        }
    } else if(b_min == 95 && b_max == 115){
        
    
        double SM_sigma_tta_bin_95_115 = SM.getOptionalParameter("SM_sigma_tta_bin_95_115_ATLAS_240309452");
        double sigma_tta_bin_95_115_madgraph = 142.5465; //fb
        
        
        if(flag_Quadratic){
        
            return  SM_sigma_tta_bin_95_115 +(-0.0100*C_phiQ3 + 0.9150*C_phiQ3*C_phiQ3 + -0.1497*C_phiQm + 0.8027*C_phiQm*C_phiQm + 
                    -0.0423*C_phit + 0.8970*C_phit*C_phit + 0.4354*C_Qd8 + 0.2846*C_Qd8*C_Qd8 + 3.0706*C_Qq18 + 1.7216*C_Qq18*C_Qq18 +
                    1.3968*C_Qq38 + 2.1580*C_Qq38*C_Qq38 + 3.0244*C_Qu8 + 1.2824*C_Qu8*C_Qu8 + 2.0654*C_td8 + 1.0597*C_td8*C_td8 + 
                    37.9456*C_tG + 8.6076*C_tG*C_tG + -0.3500*C_tphi + 0.6733*C_tphi*C_tphi + 3.6589*C_tq8 + 2.0057*C_tq8*C_tq8 + 
                    2.9744*C_tu8 + 1.0258*C_tu8*C_tu8 + 13.8124*C_tW + 37.1495*C_tW*C_tW + -10.7876*C_tZ + 28.5651*C_tZ*C_tZ)
                    *(SM_sigma_tta_bin_95_115/sigma_tta_bin_95_115_madgraph);;
            
        }
        else{
            
            return SM_sigma_tta_bin_95_115 +(-0.0100*C_phiQ3 + -0.1497*C_phiQm + -0.0423*C_phit + 0.4354*C_Qd8 +
                    3.0706*C_Qq18 + 1.3968*C_Qq38 + 3.0244*C_Qu8 + 2.0654*C_td8 + 37.9456*C_tG + 
                    -0.3500*C_tphi + 3.6589*C_tq8 + 2.9744*C_tu8 + 13.8124*C_tW + -10.7876*C_tZ)
                    *(SM_sigma_tta_bin_95_115/sigma_tta_bin_95_115_madgraph);
            
        }
    } else if(b_min == 115 && b_max == 135){
        
    
        double SM_sigma_tta_bin_115_135 = SM.getOptionalParameter("SM_sigma_tta_bin_115_135_ATLAS_240309452");
        double sigma_tta_bin_115_135_madgraph = 100.2376; //fb
        
        
        if(flag_Quadratic){
        
            return  SM_sigma_tta_bin_115_135 +(0.0982*C_phiQ3 + -0.0192*C_phiQ3*C_phiQ3 + 0.2759*C_phiQm + -0.0991*C_phiQm*C_phiQm + 
                    0.2592*C_phit + -0.0754*C_phit*C_phit + 0.3715*C_Qd8 + 0.4373*C_Qd8*C_Qd8 + 2.3608*C_Qq18 + 1.2985*C_Qq18*C_Qq18 +
                    0.2884*C_Qq38 + 1.0401*C_Qq38*C_Qq38 + 2.8584*C_Qu8 + 1.2571*C_Qu8*C_Qu8 + 1.0173*C_td8 + 0.6677*C_td8*C_td8 +
                    26.8228*C_tG + 6.2131*C_tG*C_tG + 0.2409*C_tphi + -0.0626*C_tphi*C_tphi + 3.0877*C_tq8 + 1.5817*C_tq8*C_tq8 + 
                    1.2391*C_tu8 + 1.0078*C_tu8*C_tu8 + 10.1158*C_tW + 33.9812*C_tW*C_tW + -8.8975*C_tZ + 26.5316*C_tZ*C_tZ)
                    *(SM_sigma_tta_bin_115_135/sigma_tta_bin_115_135_madgraph);;
            
        }
        else{
            
            return SM_sigma_tta_bin_115_135 +(0.0982*C_phiQ3 + 0.2759*C_phiQm + 0.2592*C_phit + 0.3715*C_Qd8 +
                    2.3608*C_Qq18 + 0.2884*C_Qq38 + 2.8584*C_Qu8 + 1.0173*C_td8 + 26.8228*C_tG + 
                    0.2409*C_tphi + 3.0877*C_tq8 + 1.2391*C_tu8 + 10.1158*C_tW + -8.8975*C_tZ)
                    *(SM_sigma_tta_bin_115_135/sigma_tta_bin_115_135_madgraph);
            
        }
    } else if(b_min == 135 && b_max == 160){
        
    
        double SM_sigma_tta_bin_135_160 = SM.getOptionalParameter("SM_sigma_tta_bin_135_160_ATLAS_240309452");
        double sigma_tta_bin_135_160_madgraph = 86.8255; //fb
        
        
        if(flag_Quadratic){
        
            return  SM_sigma_tta_bin_135_160 +(0.0443*C_phiQ3 + -0.1031*C_phiQ3*C_phiQ3 + 0.0548*C_phiQm + -0.0336*C_phiQm*C_phiQm + 
                    0.1104*C_phit + -0.1402*C_phit*C_phit + 0.4663*C_Qd8 + 0.0807*C_Qd8*C_Qd8 + 1.9391*C_Qq18 + 1.3985*C_Qq18*C_Qq18 +
                    0.1629*C_Qq38 + 0.9372*C_Qq38*C_Qq38 + 1.8818*C_Qu8 + 0.8036*C_Qu8*C_Qu8 + 0.8608*C_td8 + 0.2294*C_td8*C_td8 + 
                    22.6273*C_tG + 5.4599*C_tG*C_tG + 0.0082*C_tphi + -0.0811*C_tphi*C_tphi + 2.5293*C_tq8 + 0.9628*C_tq8*C_tq8 +
                    1.3016*C_tu8 + 0.6664*C_tu8*C_tu8 + 7.8119*C_tW + 39.2686*C_tW*C_tW + -7.0002*C_tZ + 30.1223*C_tZ*C_tZ)
                    *(SM_sigma_tta_bin_135_160/sigma_tta_bin_135_160_madgraph);;
            
        }
        else{
            
            return SM_sigma_tta_bin_135_160 +( 0.0443*C_phiQ3 + 0.0548*C_phiQm + 0.1104*C_phit + 0.4663*C_Qd8 +
                    1.9391*C_Qq18 + 0.1629*C_Qq38 + 1.8818*C_Qu8 + 0.8608*C_td8 + 22.6273*C_tG + 
                    0.0082*C_tphi + 2.5293*C_tq8 + 1.3016*C_tu8 + 7.8119*C_tW + -7.0002*C_tZ)
                    *(SM_sigma_tta_bin_135_160/sigma_tta_bin_135_160_madgraph);
            
        }
    } else if(b_min == 160 && b_max == 200){
        
    
        double SM_sigma_tta_bin_160_200 = SM.getOptionalParameter("SM_sigma_tta_bin_160_200_ATLAS_240309452");
        double sigma_tta_bin_160_200_madgraph = 85.9152; //fb
        
        
        if(flag_Quadratic){
        
            return  SM_sigma_tta_bin_160_200 +(-0.1847*C_phiQ3 + -0.2516*C_phiQ3*C_phiQ3 + -0.1580*C_phiQm + -0.2331*C_phiQm*C_phiQm +
                    -0.1703*C_phit + -0.2045*C_phit*C_phit + 0.5987*C_Qd8 + 0.3884*C_Qd8*C_Qd8 + 2.2115*C_Qq18 + 1.5935*C_Qq18*C_Qq18 +
                    0.3104*C_Qq38 + 1.1680*C_Qq38*C_Qq38 + 2.4475*C_Qu8 + 1.0485*C_Qu8*C_Qu8 + 1.3963*C_td8 + 0.4938*C_td8*C_td8 +
                    21.9619*C_tG + 5.7199*C_tG*C_tG + -0.1752*C_tphi + -0.2107*C_tphi*C_tphi + 2.7676*C_tq8 + 1.1408*C_tq8*C_tq8 +
                    1.8702*C_tu8 + 0.9288*C_tu8*C_tu8 + 7.9401*C_tW + 52.7680*C_tW*C_tW + -8.1942*C_tZ + 40.4677*C_tZ*C_tZ)
                    *(SM_sigma_tta_bin_160_200/sigma_tta_bin_160_200_madgraph);;
            
        }
        else{
            
            return SM_sigma_tta_bin_160_200 +(-0.1847*C_phiQ3 + -0.1580*C_phiQm + -0.1703*C_phit + 0.5987*C_Qd8 + 
                    2.2115*C_Qq18 + 0.3104*C_Qq38 + 2.4475*C_Qu8 + 1.3963*C_td8 + 21.9619*C_tG + 
                    -0.1752*C_tphi + 2.7676*C_tq8 + 1.8702*C_tu8 + 7.9401*C_tW + -8.1942*C_tZ)
                    *(SM_sigma_tta_bin_160_200/sigma_tta_bin_160_200_madgraph);
            
        }
    } else if(b_min == 200 && b_max == 260){
        
    
        double SM_sigma_tta_bin_200_260 = SM.getOptionalParameter("SM_sigma_tta_bin_200_260_ATLAS_240309452");
        double sigma_tta_bin_200_260_madgraph = 64.3789; //fb
        
        
        if(flag_Quadratic){
        
            return  SM_sigma_tta_bin_200_260 +(-0.0159*C_phiQ3 + -0.1659*C_phiQ3*C_phiQ3 + -0.0580*C_phiQm + -0.1858*C_phiQm*C_phiQm + 
                    0.0549*C_phit + -0.2053*C_phit*C_phit + 0.3169*C_Qd8 + -0.0456*C_Qd8*C_Qd8 + 1.6921*C_Qq18 + 1.3495*C_Qq18*C_Qq18 +
                    0.0659*C_Qq38 + 1.1830*C_Qq38*C_Qq38 + 1.7007*C_Qu8 + 0.9087*C_Qu8*C_Qu8 + 0.8325*C_td8 + 0.5036*C_td8*C_td8 + 
                    16.9626*C_tG + 5.1974*C_tG*C_tG + -0.0888*C_tphi + -0.1738*C_tphi*C_tphi + 2.1195*C_tq8 + 1.2096*C_tq8*C_tq8 + 
                    0.9658*C_tu8 + 0.7172*C_tu8*C_tu8 + 4.8947*C_tW + 58.4204*C_tW*C_tW + -4.1328*C_tZ + 44.9640*C_tZ*C_tZ)
                    *(SM_sigma_tta_bin_200_260/sigma_tta_bin_200_260_madgraph);;
            
        }
        else{
            
            return SM_sigma_tta_bin_200_260 +(-0.0159*C_phiQ3 + -0.0580*C_phiQm + 0.0549*C_phit + 0.3169*C_Qd8 
                    + 1.6921*C_Qq18 + 0.0659*C_Qq38 + 1.7007*C_Qu8 + 0.8325*C_td8 + 16.9626*C_tG + 
                    -0.0888*C_tphi + 2.1195*C_tq8 + 0.9658*C_tu8 + 4.8947*C_tW + -4.1328*C_tZ)
                    *(SM_sigma_tta_bin_200_260/sigma_tta_bin_200_260_madgraph);
            
        }
    } else if(b_min == 260 && b_max == 500){
        
    
        double SM_sigma_tta_bin_260_500 = SM.getOptionalParameter("SM_sigma_tta_bin_260_500_ATLAS_240309452");
        double sigma_tta_bin_260_500_madgraph = 54.6387; //fb
        
        
        if(flag_Quadratic){
        
            return  SM_sigma_tta_bin_260_500 +(-0.1799*C_phiQ3 + 0.0335*C_phiQ3*C_phiQ3 + -0.1870*C_phiQm + 0.0253*C_phiQm*C_phiQm + 
                    -0.1714*C_phit + 0.0470*C_phit*C_phit + 0.2100*C_Qd8 + 0.3024*C_Qd8*C_Qd8 + 1.7935*C_Qq18 + 2.2280*C_Qq18*C_Qq18 +
                    0.2762*C_Qq38 + 1.9609*C_Qq38*C_Qq38 + 1.6440*C_Qu8 + 1.5191*C_Qu8*C_Qu8 + 0.8866*C_td8 + 0.7711*C_td8*C_td8 + 
                    14.2313*C_tG + 5.8237*C_tG*C_tG + -0.1895*C_tphi + 0.0588*C_tphi*C_tphi + 2.1582*C_tq8 + 1.9149*C_tq8*C_tq8 + 
                    1.1636*C_tu8 + 1.2606*C_tu8*C_tu8 + 2.4841*C_tW + 93.5780*C_tW*C_tW + -1.0919*C_tZ + 72.0953*C_tZ*C_tZ)
                    *(SM_sigma_tta_bin_260_500/sigma_tta_bin_260_500_madgraph);;
            
        }
        else{
            
            return SM_sigma_tta_bin_260_500 +(-0.1799*C_phiQ3 + -0.1870*C_phiQm + -0.1714*C_phit + 0.2100*C_Qd8 +
                    1.7935*C_Qq18 + 0.2762*C_Qq38 + 1.6440*C_Qu8 + 0.8866*C_td8 + 14.2313*C_tG + 
                    -0.1895*C_tphi + 2.1582*C_tq8 + 1.1636*C_tu8 + 2.4841*C_tW + -1.0919*C_tZ)
                    *(SM_sigma_tta_bin_260_500/sigma_tta_bin_260_500_madgraph);
            
        }
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct binning range for sigma_tta_diff_LO_CMS_dilepton.\n");
    }

}





/////////////////////////////////////////////////////////////////////////////////////////////////////////



/////// tt differential cross section for different bins ///////////////////////////////////////////////

////// ttbar differential cross section ///////////////////////////////////

sigma_tt_diff_NLO::sigma_tt_diff_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tt_bin_250_400" << "SM_sigma_tt_bin_400_480" <<
            "SM_sigma_tt_bin_480_560" << "SM_sigma_tt_bin_560_640" << "SM_sigma_tt_bin_640_720" << "SM_sigma_tt_bin_720_800" <<
            "SM_sigma_tt_bin_800_900" << "SM_sigma_tt_bin_900_1000" << "SM_sigma_tt_bin_1000_1150" << "SM_sigma_tt_bin_1150_1300" <<
            "SM_sigma_tt_bin_1300_1500" << "SM_sigma_tt_bin_1500_1700" << "SM_sigma_tt_bin_1700_2000" << "SM_sigma_tt_bin_2000_2300" <<
            "SM_sigma_tt_bin_2300_3500" << "SM_sigma_tt_bin_2300_2600" << "SM_sigma_tt_bin_2600_3000" << "SM_sigma_tt_bin_3000_3500" <<
            "SM_sigma_tt_bin_3500_4000");
                
}

double sigma_tt_diff_NLO::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    
    
    
    
    if(b_min == 250 && b_max == 400){
        
        double SM_tt_bin_250_400 = SM.getOptionalParameter("SM_sigma_tt_bin_250_400");
                
        //double tt_bin_250_400_Madgraph = 105600.0;
        
        double tt_bin_250_400_Madgraph_NLO = 143484.180;
        
        if(flag_Quadratic){ 
            /*return  ( tt_bin_250_400_Madgraph  +32050.0*C_tG  +3932.0*C_tG*C_tG  +380.7*C_Qd8  
                    +11.0*C_Qd8*C_Qd8  +49.58*C_Qd1*C_Qd1  +68.77*C_Qu1*C_Qu1  +523.0*C_Qu8  
                    +15.29*C_Qu8*C_Qu8  +49.6*C_td1*C_td1  +380.1*C_td8  +11.0*C_td8*C_td8  
                    +116.7*C_tq1*C_tq1  +896.3*C_tq8  +25.91*C_tq8*C_tq8  +68.72*C_tu1*C_tu1  
                    +519.8*C_tu8  +15.27*C_tu8*C_tu8  +116.7*C_Qq11*C_Qq11  +116.7*C_Qq31*C_Qq31  
                    +897.8*C_Qq18  +25.9*C_Qq18*C_Qq18  +155.7*C_Qq38  +25.92*C_Qq38*C_Qq38  
                    +40.06*C_tG*C_Qd8  +56.05*C_tG*C_Qu8  +39.69*C_tG*C_td8  +93.64*C_tG*C_tq8  
                    +55.28*C_tG*C_tu8  +94.38*C_tG*C_Qq18  +17.37*C_tG*C_Qq38  +17.52*C_Qd8*C_td8  
                    +78.71*C_Qd1*C_td1  +109.1*C_Qu1*C_tu1  +24.2*C_Qu8*C_tu8  +185.1*C_tq1*C_Qq11  
                    +32.67*C_tq1*C_Qq31  +41.19*C_tq8*C_Qq18  +7.302*C_tq8*C_Qq38  +41.46*C_Qq11*C_Qq31  
                    +9.201*C_Qq18*C_Qq38 )*(SM_tt_bin_250_400/tt_bin_250_400_Madgraph);
             */
            return SM_tt_bin_250_400 + (51048.27*C_tG + 7436.42*C_tG*C_tG + 18.8*C_Qd1 +
                    92.49*C_Qd1*C_Qd1 + 496.41*C_Qd8 + -0.69*C_Qd8*C_Qd8 + -639.89*C_Qq11 +
                    255.45*C_Qq11*C_Qq11 + -371.94*C_Qq31 + 289.27*C_Qq31*C_Qq31 + 944.56*C_Qq18 +
                    50.77*C_Qq18*C_Qq18 + 390.43*C_Qq38 + 72.73*C_Qq38*C_Qq38 + -21.97*C_Qu1 +
                    172.24*C_Qu1*C_Qu1 + 818.49*C_Qu8 + 38.89*C_Qu8*C_Qu8 + -56.13*C_td1 +
                    127.77*C_td1*C_td1 + 214.57*C_td8 + -0.88*C_td8*C_td8 + 37.75*C_tq1 +
                    267.14*C_tq1*C_tq1 + 1090.76*C_tq8 + 62.92*C_tq8*C_tq8 + -438.56*C_tu1 +
                    184.36*C_tu1*C_tu1 + 659.54*C_tu8 + 35.97*C_tu8*C_tu8)*(SM_tt_bin_250_400/tt_bin_250_400_Madgraph_NLO);
            
              } 
        else{ 
            return  SM_tt_bin_250_400 + (51048.27*C_tG + 18.8*C_Qd1 + 496.41*C_Qd8 + 
                    -639.89*C_Qq11 + -371.94*C_Qq31 + 944.56*C_Qq18 + 390.43*C_Qq38 + 
                    -21.97*C_Qu1 + 818.49*C_Qu8 + -56.13*C_td1 + 214.57*C_td8 + 37.75*C_tq1 
                    + 1090.76*C_tq8 + -438.56*C_tu1 + 659.54*C_tu8)*(SM_tt_bin_250_400/tt_bin_250_400_Madgraph_NLO);
        } 
        
    } else if(b_min == 400 && b_max == 480){
        
        double SM_tt_bin_400_480 = SM.getOptionalParameter("SM_sigma_tt_bin_400_480");
        
        //double tt_bin_400_480_Madgraph = 161200.0;

        double tt_bin_400_480_Madgraph_NLO = 197878.6;

        if(flag_Quadratic){ 
            /*
            return  ( tt_bin_400_480_Madgraph  +42690.0*C_tG  +5543.0*C_tG*C_tG  +585.0*C_Qd8  +27.22*C_Qd8*C_Qd8  
                    +122.3*C_Qd1*C_Qd1  +174.3*C_Qu1*C_Qu1  +843.0*C_Qu8  +38.71*C_Qu8*C_Qu8  +122.4*C_td1*C_td1  
                    +591.4*C_td8  +27.2*C_td8*C_td8  +293.1*C_tq1*C_tq1  +1417.0*C_tq8  +65.19*C_tq8*C_tq8  
                    +174.3*C_tu1*C_tu1  +848.3*C_tu8  +38.71*C_tu8*C_tu8  +293.2*C_Qq11*C_Qq11  +293.0*C_Qq31*C_Qq31  
                    +1413.0*C_Qq18  +65.17*C_Qq18*C_Qq18  +274.0*C_Qq38  +65.11*C_Qq38*C_Qq38  +66.64*C_tG*C_Qd8  
                    +94.74*C_tG*C_Qu8  +66.3*C_tG*C_td8  +159.0*C_tG*C_tq8  +95.27*C_tG*C_tu8  +159.3*C_tG*C_Qq18  
                    +30.01*C_tG*C_Qq38  +29.48*C_Qd8*C_td8  +132.7*C_Qd1*C_td1  +188.4*C_Qu1*C_tu1  +41.99*C_Qu8*C_tu8  
                    +318.0*C_tq1*C_Qq11  +59.73*C_tq1*C_Qq31  +70.63*C_tq8*C_Qq18  +13.32*C_tq8*C_Qq38  +110.4*C_Qq11*C_Qq31  
                    +24.47*C_Qq18*C_Qq38 )*(SM_tt_bin_400_480/tt_bin_400_480_Madgraph);
             */
            return SM_tt_bin_400_480 + (63428.12*C_tG + 10154.31*C_tG*C_tG + 49.03*C_Qd1 +
                    267.6*C_Qd1*C_Qd1 + 722.45*C_Qd8 + 47.41*C_Qd8*C_Qd8 + -225.72*C_Qq11 +
                    610.4*C_Qq11*C_Qq11 + 159.37*C_Qq31 + 578.25*C_Qq31*C_Qq31 + 1421.57*C_Qq18 +
                    116.34*C_Qq18*C_Qq18 + 187.29*C_Qq38 + 110.38*C_Qq38*C_Qq38 + -15.83*C_Qu1 +
                    362.64*C_Qu1*C_Qu1 + 1036.59*C_Qu8 + 66.42*C_Qu8*C_Qu8 + -353.63*C_td1 +
                    266.92*C_td1*C_td1 + 485.47*C_td8 + 36.59*C_td8*C_td8 + 129.28*C_tq1 +
                    579.16*C_tq1*C_tq1 + 1844.81*C_tq8 + 87.31*C_tq8*C_tq8 + -254.37*C_tu1 +
                    370.44*C_tu1*C_tu1 + 874.28*C_tu8 + 49.92*C_tu8*C_tu8)*(SM_tt_bin_400_480/tt_bin_400_480_Madgraph_NLO);
        } 
        else{ 
            return  SM_tt_bin_400_480 + (63428.12*C_tG + 49.03*C_Qd1 + 722.45*C_Qd8 + -225.72*C_Qq11 + 159.37*C_Qq31 +
                    1421.57*C_Qq18 + 187.29*C_Qq38 + -15.83*C_Qu1 + 1036.59*C_Qu8 + -353.63*C_td1 +
                    485.47*C_td8 + 129.28*C_tq1 + 1844.81*C_tq8 + -254.37*C_tu1 + 874.28*C_tu8)*(SM_tt_bin_400_480/tt_bin_400_480_Madgraph_NLO);
        } 

        
    } else if(b_min == 480 && b_max == 560){
        
        double SM_tt_bin_480_560 = SM.getOptionalParameter("SM_sigma_tt_bin_480_560");
        
        //double tt_bin_480_560_Madgraph = 102200.0;

        double tt_bin_480_560_Madgraph_NLO = 120034.78;


        if(flag_Quadratic){ 
            /*
            return  ( tt_bin_480_560_Madgraph  +25440.0*C_tG  +3642.0*C_tG*C_tG  +442.4*C_Qd8  +31.84*C_Qd8*C_Qd8  
                    +143.0*C_Qd1*C_Qd1  +209.6*C_Qu1*C_Qu1  +646.7*C_Qu8  +46.54*C_Qu8*C_Qu8  +143.2*C_td1*C_td1  
                    +439.2*C_td8  +31.82*C_td8*C_td8  +349.2*C_tq1*C_tq1  +1076.0*C_tq8  +77.6*C_tq8*C_tq8  
                    +209.6*C_tu1*C_tu1  +642.9*C_tu8  +46.61*C_tu8*C_tu8  +349.5*C_Qq11*C_Qq11  +349.5*C_Qq31*C_Qq31  
                    +1082.0*C_Qq18  +77.63*C_Qq18*C_Qq18  +204.3*C_Qq38  +77.61*C_Qq38*C_Qq38  +53.38*C_tG*C_Qd8  
                    +78.24*C_tG*C_Qu8  +53.22*C_tG*C_td8  +130.0*C_tG*C_tq8  +77.83*C_tG*C_tu8  +130.6*C_tG*C_Qq18  
                    +25.94*C_tG*C_Qq38  +23.54*C_Qd8*C_td8  +106.7*C_Qd1*C_td1  +155.0*C_Qu1*C_tu1  +34.55*C_Qu8*C_tu8  
                    +258.3*C_tq1*C_Qq11  +51.59*C_tq1*C_Qq31  +57.56*C_tq8*C_Qq18  +11.5*C_tq8*C_Qq38  +139.5*C_Qq11*C_Qq31  
                    +31.04*C_Qq18*C_Qq38 )*(SM_tt_bin_480_560/tt_bin_480_560_Madgraph);
            */
            return SM_tt_bin_480_560 + (35765.48*C_tG + 6133.63*C_tG*C_tG + 43.67*C_Qd1 + 242.09*C_Qd1*C_Qd1 +
                    561.49*C_Qd8 + 46.79*C_Qd8*C_Qd8 + -63.85*C_Qq11 + 611.64*C_Qq11*C_Qq11 + 35.69*C_Qq31 +
                    636.75*C_Qq31*C_Qq31 + 1036.88*C_Qq18 + 110.93*C_Qq18*C_Qq18 + 189.86*C_Qq38 +
                    125.12*C_Qq38*C_Qq38 + 144.46*C_Qu1 + 367.25*C_Qu1*C_Qu1 + 625.84*C_Qu8 + 61.91*C_Qu8*C_Qu8 +
                    -73.57*C_td1 + 270.43*C_td1*C_td1 + 452.89*C_td8 + 38.34*C_td8*C_td8 + 107.19*C_tq1 +
                    631.27*C_tq1*C_tq1 + 1078.92*C_tq8 + 111.62*C_tq8*C_tq8 + -17.36*C_tu1 + 376.76*C_tu1*C_tu1 +
                    649.1*C_tu8 + 73.91*C_tu8*C_tu8)*(SM_tt_bin_480_560/tt_bin_480_560_Madgraph_NLO);
        } 
        else{ 
            return  SM_tt_bin_480_560 + (35765.48*C_tG + 43.67*C_Qd1 + 561.49*C_Qd8 + -63.85*C_Qq11 + 35.69*C_Qq31 +
                    1036.88*C_Qq18 + 189.86*C_Qq38 + 144.46*C_Qu1 + 625.84*C_Qu8 + -73.57*C_td1 + 452.89*C_td8 +
                    107.19*C_tq1 + 1078.92*C_tq8 + -17.36*C_tu1 + 649.1*C_tu8)*(SM_tt_bin_480_560/tt_bin_480_560_Madgraph_NLO);
        } 
        
    } else if(b_min == 560 && b_max == 640){
        
        double SM_tt_bin_560_640 = SM.getOptionalParameter("SM_sigma_tt_bin_560_640");
        
        //double tt_bin_560_640_Madgraph = 60110.0;
        
        double tt_bin_560_640_Madgraph_NLO = 68052.18 ;
        
        if(flag_Quadratic){ 
            /*
            return  ( tt_bin_560_640_Madgraph  +14650.0*C_tG  +2324.0*C_tG*C_tG  +315.2*C_Qd8  +33.06*C_Qd8*C_Qd8  
                    +148.7*C_Qd1*C_Qd1  +223.1*C_Qu1*C_Qu1  +472.2*C_Qu8  +49.62*C_Qu8*C_Qu8  +148.7*C_td1*C_td1  
                    +321.7*C_td8  +33.09*C_td8*C_td8  +368.9*C_tq1*C_tq1  +788.7*C_tq8  +81.95*C_tq8*C_tq8  
                    +223.2*C_tu1*C_tu1  +479.7*C_tu8  +49.65*C_tu8*C_tu8  +369.0*C_Qq11*C_Qq11  +368.9*C_Qq31*C_Qq31  
                    +794.9*C_Qq18  +82.0*C_Qq18*C_Qq18  +155.0*C_Qq38  +81.97*C_Qq38*C_Qq38  +40.46*C_tG*C_Qd8  
                    +60.78*C_tG*C_Qu8  +40.29*C_tG*C_td8  +100.2*C_tG*C_tq8  +60.75*C_tG*C_tu8  +100.5*C_tG*C_Qq18  
                    +21.03*C_tG*C_Qq38  +17.8*C_Qd8*C_td8  +80.35*C_Qd1*C_td1  +120.5*C_Qu1*C_tu1  +26.72*C_Qu8*C_tu8  
                    +198.0*C_tq1*C_Qq11  +41.9*C_tq1*C_Qq31  +44.29*C_tq8*C_Qq18  +9.373*C_tq8*C_Qq38  +154.9*C_Qq11*C_Qq31  
                    +34.33*C_Qq18*C_Qq38 )*(SM_tt_bin_560_640/tt_bin_560_640_Madgraph);
            */
            
            return SM_tt_bin_560_640 + (19579.39*C_tG + 3635.9*C_tG*C_tG + 17.7*C_Qd1 + 259.49*C_Qd1*C_Qd1 +
                    349.67*C_Qd8 + 50.72*C_Qd8*C_Qd8 + -75.33*C_Qq11 + 636.76*C_Qq11*C_Qq11 + -33.2*C_Qq31 +
                    613.42*C_Qq31*C_Qq31 + 741.75*C_Qq18 + 116.38*C_Qq18*C_Qq18 + 178.11*C_Qq38 +
                    100.08*C_Qq38*C_Qq38 + 7.0*C_Qu1 + 377.77*C_Qu1*C_Qu1 + 434.43*C_Qu8 +
                    76.59*C_Qu8*C_Qu8 + 40.96*C_td1 + 257.81*C_td1*C_td1 + 293.63*C_td8 +
                    43.49*C_td8*C_td8 + -131.65*C_tq1 + 597.97*C_tq1*C_tq1 + 909.1*C_tq8 +
                    124.18*C_tq8*C_tq8 + 18.95*C_tu1 + 358.67*C_tu1*C_tu1 + 496.91*C_tu8 +
                    70.62*C_tu8*C_tu8)*(SM_tt_bin_560_640/tt_bin_560_640_Madgraph_NLO);
        } 
        else{ 
            return  SM_tt_bin_560_640 + (19579.39*C_tG + 17.7*C_Qd1 + 349.67*C_Qd8 +  -75.33*C_Qq11 +
                    -33.2*C_Qq31 + 741.75*C_Qq18 + 178.11*C_Qq38 + 7.0*C_Qu1 + 434.43*C_Qu8 +
                    40.96*C_td1 + 293.63*C_td8 + -131.65*C_tq1 + 909.1*C_tq8 + 18.95*C_tu1 +
                    496.91*C_tu8)*(SM_tt_bin_560_640/tt_bin_560_640_Madgraph_NLO);
        } 
        
    } else if(b_min == 640 && b_max == 720){
        
        double SM_tt_bin_640_720 = SM.getOptionalParameter("SM_sigma_tt_bin_640_720");
        
        //double tt_bin_640_720_Madgraph = 35680.0;
        
        double tt_bin_640_720_Madgraph_NLO = 38544.57;


        if(flag_Quadratic){ 
            /*
            return  ( tt_bin_640_720_Madgraph  +8622.0*C_tG  +1513.0*C_tG*C_tG  +232.3*C_Qd8  +32.74*C_Qd8*C_Qd8  
                    +147.4*C_Qd1*C_Qd1  +226.3*C_Qu1*C_Qu1  +365.2*C_Qu8  +50.27*C_Qu8*C_Qu8  +147.5*C_td1*C_td1  
                    +227.6*C_td8  +32.77*C_td8*C_td8  +371.2*C_tq1*C_tq1  +581.3*C_tq8  +82.52*C_tq8*C_tq8  
                    +226.3*C_tu1*C_tu1  +366.3*C_tu8  +50.32*C_tu8*C_tu8  +371.1*C_Qq11*C_Qq11  +371.4*C_Qq31*C_Qq31  
                    +596.7*C_Qq18  +82.47*C_Qq18*C_Qq18  +126.3*C_Qq38  +82.48*C_Qq38*C_Qq38  +30.75*C_tG*C_Qd8  
                    +46.7*C_tG*C_Qu8  +30.34*C_tG*C_td8  +76.88*C_tG*C_tq8  +46.77*C_tG*C_tu8  +76.48*C_tG*C_Qq18  
                    +16.6*C_tG*C_Qq38  +13.52*C_Qd8*C_td8  +60.89*C_Qd1*C_td1  +92.53*C_Qu1*C_tu1  +20.68*C_Qu8*C_tu8  
                    +153.0*C_tq1*C_Qq11  +33.47*C_tq1*C_Qq31  +34.06*C_tq8*C_Qq18  +7.344*C_tq8*C_Qq38  
                    +162.6*C_Qq11*C_Qq31  +36.22*C_Qq18*C_Qq38 )*(SM_tt_bin_640_720/tt_bin_640_720_Madgraph);
             */
            return SM_tt_bin_640_720 + (11159.52*C_tG + 2360.1*C_tG*C_tG + 28.06*C_Qd1 + 240.9*C_Qd1*C_Qd1 +
                    235.18*C_Qd8 + 46.44*C_Qd8*C_Qd8 + 58.0*C_Qq11 + 568.17*C_Qq11*C_Qq11 + 43.91*C_Qq31 +
                    560.2*C_Qq31*C_Qq31 + 512.73*C_Qq18 + 101.03*C_Qq18*C_Qq18 + 104.37*C_Qq38 +
                    107.83*C_Qq38*C_Qq38 + 23.63*C_Qu1 + 350.38*C_Qu1*C_Qu1 + 343.54*C_Qu8 + 
                    79.86*C_Qu8*C_Qu8 + 4.39*C_td1 + 227.77*C_td1*C_td1 + 206.7*C_td8 +
                    39.61*C_td8*C_td8 + 38.64*C_tq1 + 574.41*C_tq1*C_tq1 + 527.14*C_tq8 +
                    112.35*C_tq8*C_tq8 + -125.48*C_tu1 + 359.4*C_tu1*C_tu1 + 331.56*C_tu8 +
                    69.15*C_tu8*C_tu8)*(SM_tt_bin_640_720/tt_bin_640_720_Madgraph_NLO);
        } 
        else{ 
            return  SM_tt_bin_640_720 + (11159.52*C_tG + 28.06*C_Qd1 + 235.18*C_Qd8 + 58.0*C_Qq11 +
                    43.91*C_Qq31 + 512.73*C_Qq18 + 104.37*C_Qq38 + 23.63*C_Qu1 + 343.54*C_Qu8 +
                    4.39*C_td1 + 206.7*C_td8 + 38.64*C_tq1 + 527.14*C_tq8 + -125.48*C_tu1 +
                    331.56*C_tu8)*(SM_tt_bin_640_720/tt_bin_640_720_Madgraph_NLO);
        } 
        
    } else if(b_min == 720 && b_max == 800){
        
        double SM_tt_bin_720_800 = SM.getOptionalParameter("SM_sigma_tt_bin_720_800");

        //double tt_bin_720_800_Madgraph = 21500.0;

        double tt_bin_720_800_Madgraph_NLO = 22430.21;

        if(flag_Quadratic){ 
            /*
            return  (tt_bin_720_800_Madgraph  +5176.0*C_tG  +1012.0*C_tG*C_tG  +176.6*C_Qd8  +31.74*C_Qd8*C_Qd8  
                    +142.8*C_Qd1*C_Qd1  +223.7*C_Qu1*C_Qu1  +273.6*C_Qu8  +49.7*C_Qu8*C_Qu8  +142.9*C_td1*C_td1  
                    +172.5*C_td8  +31.78*C_td8*C_td8  +364.0*C_tq1*C_tq1  +442.9*C_tq8  +80.93*C_tq8*C_tq8  
                    +223.6*C_tu1*C_tu1  +269.2*C_tu8  +49.69*C_tu8*C_tu8  +364.2*C_Qq11*C_Qq11  +364.2*C_Qq31*C_Qq31  
                    +447.9*C_Qq18  +80.96*C_Qq18*C_Qq18  +98.19*C_Qq38  +80.97*C_Qq38*C_Qq38  +22.86*C_tG*C_Qd8  
                    +36.51*C_tG*C_Qu8  +23.09*C_tG*C_td8  +59.14*C_tG*C_tq8  +36.42*C_tG*C_tu8  +59.37*C_tG*C_Qq18  
                    +13.48*C_tG*C_Qq38  +10.24*C_Qd8*C_td8  +46.61*C_Qd1*C_td1  +72.78*C_Qu1*C_tu1  +16.2*C_Qu8*C_tu8  
                    +118.8*C_tq1*C_Qq11  +26.93*C_tq1*C_Qq31  +26.32*C_tq8*C_Qq18  +6.01*C_tq8*C_Qq38  +167.1*C_Qq11*C_Qq31  
                    +36.94*C_Qq18*C_Qq38 )*(SM_tt_bin_720_800/tt_bin_720_800_Madgraph);
            */
            return SM_tt_bin_720_800 + (6356.56*C_tG + 1487.06*C_tG*C_tG + -82.97*C_Qd1 +
                    239.34*C_Qd1*C_Qd1 + 218.16*C_Qd8 + 56.15*C_Qd8*C_Qd8 + 53.21*C_Qq11 +
                    496.4*C_Qq11*C_Qq11 + -84.85*C_Qq31 + 504.56*C_Qq31*C_Qq31 + 361.61*C_Qq18 +
                    98.04*C_Qq18*C_Qq18 + 23.39*C_Qq38 + 142.44*C_Qq38*C_Qq38 + -19.16*C_Qu1 +
                    352.96*C_Qu1*C_Qu1 + 218.82*C_Qu8 + 58.84*C_Qu8*C_Qu8 + -15.52*C_td1 +
                    242.46*C_td1*C_td1 + 185.82*C_td8 + 25.72*C_td8*C_td8 + -94.94*C_tq1 +
                    569.29*C_tq1*C_tq1 + 495.25*C_tq8 + 143.9*C_tq8*C_tq8 + -20.18*C_tu1 +
                    326.49*C_tu1*C_tu1 + 203.44*C_tu8 + 98.99*C_tu8*C_tu8)*(SM_tt_bin_720_800/tt_bin_720_800_Madgraph_NLO);
        } 
        else{ 
            return  SM_tt_bin_720_800 + (6356.56*C_tG + -82.97*C_Qd1 + 218.16*C_Qd8 + 53.21*C_Qq11 +
                    -84.85*C_Qq31 + 361.61*C_Qq18 + 23.39*C_Qq38 + -19.16*C_Qu1 + 218.82*C_Qu8 +
                    -15.52*C_td1 + 185.82*C_td8 + -94.94*C_tq1 + 495.25*C_tq8 + -20.18*C_tu1 +
                    203.44*C_tu8)*(SM_tt_bin_720_800/tt_bin_720_800_Madgraph_NLO);
        } 

        
    } else if(b_min == 800 && b_max == 900){
        
        double SM_tt_bin_800_900 = SM.getOptionalParameter("SM_sigma_tt_bin_800_900");
        
        
        //double tt_bin_800_900_Madgraph = 15720.0;
        
        double tt_bin_800_900_Madgraph_NLO = 15975.89;

        
        if(flag_Quadratic){ 
            /*
            return  ( tt_bin_800_900_Madgraph  +3783.0*C_tG  +831.2*C_tG*C_tG  +163.3*C_Qd8  +37.66*C_Qd8*C_Qd8  
                    +169.4*C_Qd1*C_Qd1  +271.2*C_Qu1*C_Qu1  +257.6*C_Qu8  +60.21*C_Qu8*C_Qu8  +169.4*C_td1*C_td1  
                    +155.3*C_td8  +37.65*C_td8*C_td8  +438.5*C_tq1*C_tq1  +416.1*C_tq8  +97.43*C_tq8*C_tq8  
                    +271.0*C_tu1*C_tu1  +252.0*C_tu8  +60.3*C_tu8*C_tu8  +438.3*C_Qq11*C_Qq11  +438.3*C_Qq31*C_Qq31  
                    +413.8*C_Qq18  +97.38*C_Qq18*C_Qq18  +97.07*C_Qq38  +97.4*C_Qq38*C_Qq38  +22.22*C_tG*C_Qd8  
                    +35.3*C_tG*C_Qu8  +22.04*C_tG*C_td8  +56.83*C_tG*C_tq8  +35.06*C_tG*C_tu8  +57.14*C_tG*C_Qq18  
                    +13.71*C_tG*C_Qq38  +9.766*C_Qd8*C_td8  +43.65*C_Qd1*C_td1  +69.77*C_Qu1*C_tu1  +15.54*C_Qu8*C_tu8  
                    +112.0*C_tq1*C_Qq11  +26.27*C_tq1*C_Qq31  +24.97*C_tq8*C_Qq18  +5.882*C_tq8*C_Qq38  
                    +208.3*C_Qq11*C_Qq31  +46.05*C_Qq18*C_Qq38 )*(SM_tt_bin_800_900/tt_bin_800_900_Madgraph);
            */
            return SM_tt_bin_800_900 + (4480.84*C_tG + 1117.83*C_tG*C_tG + -11.57*C_Qd1 + 223.02*C_Qd1*C_Qd1 +
                    155.41*C_Qd8 + 16.57*C_Qd8*C_Qd8 + 11.71*C_Qq11 + 558.35*C_Qq11*C_Qq11 + 8.23*C_Qq31 +
                    545.38*C_Qq31*C_Qq31 + 345.66*C_Qq18 + 122.04*C_Qq18*C_Qq18 + 66.06*C_Qq38 + 83.51*C_Qq38*C_Qq38 +
                    -29.94*C_Qu1 + 388.28*C_Qu1*C_Qu1 + 229.26*C_Qu8 + 107.91*C_Qu8*C_Qu8 + -10.93*C_td1 +
                    188.8*C_td1*C_td1 + 187.13*C_td8 + -0.04*C_td8*C_td8 + -6.36*C_tq1 + 564.88*C_tq1*C_tq1 +
                    374.3*C_tq8 + 113.29*C_tq8*C_tq8 + -24.56*C_tu1 + 375.1*C_tu1*C_tu1 + 235.34*C_tu8 +
                    41.63*C_tu8*C_tu8)*(SM_tt_bin_800_900/tt_bin_800_900_Madgraph_NLO);
        } 
        
        else{ 
            return  SM_tt_bin_800_900 + (4480.84*C_tG + -11.57*C_Qd1 + 155.41*C_Qd8 + 11.71*C_Qq11 +
                    8.23*C_Qq31 + 345.66*C_Qq18 + 66.06*C_Qq38 + -29.94*C_Qu1 + 229.26*C_Qu8 +
                    -10.93*C_td1 + 187.13*C_td8 + -6.36*C_tq1 + 374.3*C_tq8 + -24.56*C_tu1 +
                    235.34*C_tu8)*(SM_tt_bin_800_900/tt_bin_800_900_Madgraph_NLO);
        } 
        
    } else if(b_min == 900 && b_max == 1000){
        
        double SM_tt_bin_900_1000 = SM.getOptionalParameter("SM_sigma_tt_bin_900_1000");
        
        //double tt_bin_900_1000_Madgraph = 9102.0;
        
        double tt_bin_900_1000_Madgraph_NLO = 8715.37;
        
        if(flag_Quadratic){ 
            /*
            return  ( tt_bin_900_1000_Madgraph  +2187.0*C_tG  +537.9*C_tG*C_tG  +117.2*C_Qd8  +35.16*C_Qd8*C_Qd8  
                    +158.2*C_Qd1*C_Qd1  +258.4*C_Qu1*C_Qu1  +195.1*C_Qu8  +57.45*C_Qu8*C_Qu8  +158.2*C_td1*C_td1  
                    +120.8*C_td8  +35.15*C_td8*C_td8  +414.8*C_tq1*C_tq1  +306.6*C_tq8  +92.23*C_tq8*C_tq8  
                    +258.3*C_tu1*C_tu1  +187.9*C_tu8  +57.43*C_tu8*C_tu8  +414.9*C_Qq11*C_Qq11  +414.5*C_Qq31*C_Qq31  
                    +307.8*C_Qq18  +92.26*C_Qq18*C_Qq18  +72.16*C_Qq38  +92.2*C_Qq38*C_Qq38  +16.01*C_tG*C_Qd8  
                    +26.12*C_tG*C_Qu8  +16.05*C_tG*C_td8  +42.41*C_tG*C_tq8  +26.27*C_tG*C_tu8  +41.78*C_tG*C_Qq18  
                    +10.13*C_tG*C_Qq38  +7.134*C_Qd8*C_td8  +32.04*C_Qd1*C_td1  +53.54*C_Qu1*C_tu1  +11.67*C_Qu8*C_tu8  
                    +84.75*C_tq1*C_Qq11  +20.83*C_tq1*C_Qq31  +18.68*C_tq8*C_Qq18  +4.627*C_tq8*C_Qq38  +203.8*C_Qq11*C_Qq31  
                    +45.21*C_Qq18*C_Qq38 )*(SM_tt_bin_900_1000/tt_bin_900_1000_Madgraph);
        */
            return  SM_tt_bin_900_1000 + (2494.7*C_tG + 735.75*C_tG*C_tG + -12.85*C_Qd1 +
                    209.63*C_Qd1*C_Qd1 + 80.87*C_Qd8 + 65.31*C_Qd8*C_Qd8 + -43.19*C_Qq11 +
                    590.18*C_Qq11*C_Qq11 + -47.86*C_Qq31 + 536.79*C_Qq31*C_Qq31 + 238.66*C_Qq18 +
                    102.0*C_Qq18*C_Qq18 + 18.29*C_Qq38 + 100.5*C_Qq38*C_Qq38 + -35.45*C_Qu1 +
                    362.42*C_Qu1*C_Qu1 + 165.59*C_Qu8 + 98.04*C_Qu8*C_Qu8 + -18.47*C_td1 +
                    212.01*C_td1*C_td1 + 68.99*C_td8 + 37.47*C_td8*C_td8 + -18.92*C_tq1 +
                    584.25*C_tq1*C_tq1 + 248.96*C_tq8 + 119.78*C_tq8*C_tq8 + 7.38*C_tu1 +
                    333.01*C_tu1*C_tu1 + 154.19*C_tu8 + 70.58*C_tu8*C_tu8)*(SM_tt_bin_900_1000/tt_bin_900_1000_Madgraph_NLO);

        } 
        else{ 
            return  SM_tt_bin_900_1000 + (2494.7*C_tG + -12.85*C_Qd1 + 80.87*C_Qd8 + -43.19*C_Qq11 +
                    -47.86*C_Qq31 + 238.66*C_Qq18 + 18.29*C_Qq38 + -35.45*C_Qu1 + 165.59*C_Qu8 +
                    -18.47*C_td1 + 68.99*C_td8 + -18.92*C_tq1 + 248.96*C_tq8 + 7.38*C_tu1 +
                    154.19*C_tu8 )*(SM_tt_bin_900_1000/tt_bin_900_1000_Madgraph_NLO);
        } 
        
    } else if(b_min == 1000 && b_max == 1150){
        
        double SM_tt_bin_1000_1150 = SM.getOptionalParameter("SM_sigma_tt_bin_1000_1150");
        
        //double tt_bin_1000_1150_Madgraph = 7170.0;
        
        double tt_bin_1000_1150_Madgraph_NLO = 6674.6;

        
        if(flag_Quadratic){ 
            /*
            return  ( tt_bin_1000_1150_Madgraph  +1717.0*C_tG  +491.3*C_tG*C_tG  +124.2*C_Qd8  +47.87*C_Qd8*C_Qd8  
                    +215.2*C_Qd1*C_Qd1  +360.4*C_Qu1*C_Qu1  +214.7*C_Qu8  +80.15*C_Qu8*C_Qu8  +215.4*C_td1*C_td1  
                    +122.0*C_td8  +47.86*C_td8*C_td8  +573.5*C_tq1*C_tq1  +326.1*C_tq8  +127.5*C_tq8*C_tq8  
                    +360.3*C_tu1*C_tu1  +198.4*C_tu8  +80.1*C_tu8*C_tu8  +573.7*C_Qq11*C_Qq11  +573.7*C_Qq31*C_Qq31  
                    +327.1*C_Qq18  +127.4*C_Qq18*C_Qq18  +94.46*C_Qq38  +127.5*C_Qq38*C_Qq38  +17.18*C_tG*C_Qd8  
                    +28.64*C_tG*C_Qu8  +17.08*C_tG*C_td8  +45.9*C_tG*C_tq8  +28.77*C_tG*C_tu8  +46.22*C_tG*C_Qq18  
                    +11.44*C_tG*C_Qq38  +7.645*C_Qd8*C_td8  +34.45*C_Qd1*C_td1  +57.58*C_Qu1*C_tu1  +12.63*C_Qu8*C_tu8  
                    +91.17*C_tq1*C_Qq11  +23.68*C_tq1*C_Qq31  +20.3*C_tq8*C_Qq18  +5.195*C_tq8*C_Qq38  
                    +294.7*C_Qq11*C_Qq31  +65.52*C_Qq18*C_Qq38 )*(SM_tt_bin_1000_1150/tt_bin_1000_1150_Madgraph);
            */
            return SM_tt_bin_1000_1150 + (1860.42*C_tG + 596.18*C_tG*C_tG + -3.45*C_Qd1 + 266.03*C_Qd1*C_Qd1 +
                   98.52*C_Qd8 + 46.83*C_Qd8*C_Qd8 + -4.02*C_Qq11 + 681.56*C_Qq11*C_Qq11 + 1.41*C_Qq31 +
                    711.4*C_Qq31*C_Qq31 + 266.48*C_Qq18 + 133.38*C_Qq18*C_Qq18 + 72.8*C_Qq38 +
                    142.45*C_Qq38*C_Qq38 + -26.25*C_Qu1 + 444.18*C_Qu1*C_Qu1 + 166.62*C_Qu8 +
                    91.69*C_Qu8*C_Qu8 + -17.86*C_td1 + 262.28*C_td1*C_td1 + 99.56*C_td8 + 46.39*C_td8*C_td8 +
                    -26.04*C_tq1 + 705.1*C_tq1*C_tq1 + 262.36*C_tq8 + 134.7*C_tq8*C_tq8 + 0.17*C_tu1 +
                    430.39*C_tu1*C_tu1 + 179.34*C_tu8 + 85.09*C_tu8*C_tu8)*(SM_tt_bin_1000_1150/tt_bin_1000_1150_Madgraph_NLO);
        } 
        else{ 
            return  SM_tt_bin_1000_1150 + (1860.42*C_tG + -3.45*C_Qd1 + 98.52*C_Qd8 + -4.02*C_Qq11 +
                    1.41*C_Qq31 + 266.48*C_Qq18 + 72.8*C_Qq38 + -26.25*C_Qu1 + 166.62*C_Qu8 + -17.86*C_td1 +
                    99.56*C_td8 + -26.04*C_tq1 + 262.36*C_tq8 + 0.17*C_tu1 + 179.34*C_tu8)*(SM_tt_bin_1000_1150/tt_bin_1000_1150_Madgraph_NLO);
        } 
        
    } else if(b_min == 1150 && b_max == 1300){
        
        double SM_tt_bin_1150_1300 = SM.getOptionalParameter("SM_sigma_tt_bin_1150_1300");
        
        //double tt_bin_1150_1300_Madgraph = 3513.0;
        
        double tt_bin_1150_1300_Madgraph_NLO = 3064.41;


        if(flag_Quadratic){ 
            
            /*
            return  ( tt_bin_1150_1300_Madgraph  +842.0*C_tG  +281.2*C_tG*C_tG  +82.6*C_Qd8  +42.2*C_Qd8*C_Qd8  
                    +189.9*C_Qd1*C_Qd1  +326.4*C_Qu1*C_Qu1  +131.6*C_Qu8  +72.52*C_Qu8*C_Qu8  +189.9*C_td1*C_td1  
                    +80.17*C_td8  +42.18*C_td8*C_td8  +514.8*C_tq1*C_tq1  +223.2*C_tq8  +114.4*C_tq8*C_tq8  
                    +326.5*C_tu1*C_tu1  +132.9*C_tu8  +72.48*C_tu8*C_tu8  +514.7*C_Qq11*C_Qq11  +514.5*C_Qq31*C_Qq31  
                    +216.9*C_Qq18  +114.4*C_Qq18*C_Qq18  +58.18*C_Qq38  +114.4*C_Qq38*C_Qq38  +11.65*C_tG*C_Qd8  
                    +20.08*C_tG*C_Qu8  +11.8*C_tG*C_td8  +31.39*C_tG*C_tq8  +20.16*C_tG*C_tu8  +31.14*C_tG*C_Qq18  
                    +8.065*C_tG*C_Qq38  +5.124*C_Qd8*C_td8  +22.94*C_Qd1*C_td1  +39.76*C_Qu1*C_tu1  +8.918*C_Qu8*C_tu8  
                    +62.47*C_tq1*C_Qq11  +17.32*C_tq1*C_Qq31  +13.79*C_tq8*C_Qq18  +3.564*C_tq8*C_Qq38  +276.9*C_Qq11*C_Qq31  
                    +61.08*C_Qq18*C_Qq38 )*(SM_tt_bin_1150_1300/tt_bin_1150_1300_Madgraph);
            */
            return SM_tt_bin_1150_1300 + (858.66*C_tG + 323.66*C_tG*C_tG + -6.81*C_Qd1 + 214.97*C_Qd1*C_Qd1 +
                    59.67*C_Qd8 + 47.21*C_Qd8*C_Qd8 + -4.64*C_Qq11 + 596.93*C_Qq11*C_Qq11 + -5.53*C_Qq31 +
                    597.72*C_Qq31*C_Qq31 + 166.04*C_Qq18 + 118.48*C_Qq18*C_Qq18 + 43.59*C_Qq38 +
                    121.4*C_Qq38*C_Qq38 + -17.07*C_Qu1 + 375.94*C_Qu1*C_Qu1 + 99.12*C_Qu8 +
                    78.84*C_Qu8*C_Qu8 + 2.94*C_td1 + 214.14*C_td1*C_td1 + 63.01*C_td8 + 45.1*C_td8*C_td8 +
                    -61.85*C_tq1 + 587.77*C_tq1*C_tq1 + 157.72*C_tq8 + 116.96*C_tq8*C_tq8 + -1.08*C_tu1 +
                    375.84*C_tu1*C_tu1 + 100.16*C_tu8 + 79.08*C_tu8*C_tu8)*(SM_tt_bin_1150_1300/tt_bin_1150_1300_Madgraph_NLO);
        } 
        else{ 
            return SM_tt_bin_1150_1300 + (858.66*C_tG + -6.81*C_Qd1 + 59.67*C_Qd8 + -4.64*C_Qq11 +
                    -5.53*C_Qq31 + 166.04*C_Qq18 + 43.59*C_Qq38 + -17.07*C_Qu1 + 99.12*C_Qu8 + 
                    2.94*C_td1 + 63.01*C_td8 + -61.85*C_tq1 + 157.72*C_tq8 + -1.08*C_tu1 +
                    100.16*C_tu8)*(SM_tt_bin_1150_1300/tt_bin_1150_1300_Madgraph_NLO);
        } 
        
    } else if(b_min == 1300 && b_max == 1500){
        
        double SM_tt_bin_1300_1500 = SM.getOptionalParameter("SM_sigma_tt_bin_1300_1500");
        
        //double tt_bin_1300_1500_Madgraph = 2178.0;
        
        double tt_bin_1300_1500_Madgraph_NLO = 1821.32;

        if(flag_Quadratic){ 
            /*
            return  ( tt_bin_1300_1500_Madgraph  +520.9*C_tG  +207.9*C_tG*C_tG  +71.72*C_Qd8  +48.17*C_Qd8*C_Qd8  
                    +216.5*C_Qd1*C_Qd1  +383.3*C_Qu1*C_Qu1  +118.7*C_Qu8  +85.17*C_Qu8*C_Qu8  +216.5*C_td1*C_td1  
                    +70.4*C_td8  +48.12*C_td8*C_td8  +598.4*C_tq1*C_tq1  +193.2*C_tq8  +133.0*C_tq8*C_tq8  
                    +383.4*C_tu1*C_tu1  +129.2*C_tu8  +85.16*C_tu8*C_tu8  +598.5*C_Qq11*C_Qq11  +598.5*C_Qq31*C_Qq31  
                    +188.4*C_Qq18  +133.0*C_Qq18*C_Qq18  +70.99*C_Qq38  +133.0*C_Qq38*C_Qq38  +10.14*C_tG*C_Qd8  
                    +17.99*C_tG*C_Qu8  +10.25*C_tG*C_td8  +28.47*C_tG*C_tq8  +18.19*C_tG*C_tu8  +28.3*C_tG*C_Qq18  
                    +8.221*C_tG*C_Qq38  +4.415*C_Qd8*C_td8  +20.18*C_Qd1*C_td1  +35.76*C_Qu1*C_tu1  +7.912*C_Qu8*C_tu8  
                    +56.09*C_tq1*C_Qq11  +15.09*C_tq1*C_Qq31  +12.33*C_tq8*C_Qq18  +3.508*C_tq8*C_Qq38  +336.9*C_Qq11*C_Qq31  
                    +74.73*C_Qq18*C_Qq38 )*(SM_tt_bin_1300_1500/tt_bin_1300_1500_Madgraph);
            */
            return SM_tt_bin_1300_1500 + (501.38*C_tG + 215.53*C_tG*C_tG + -6.74*C_Qd1 + 232.32*C_Qd1*C_Qd1 +
                    47.09*C_Qd8 + 49.01*C_Qd8*C_Qd8 + 10.05*C_Qq11 + 644.43*C_Qq11*C_Qq11 + 5.8*C_Qq31 +
                    640.35*C_Qq31*C_Qq31 + 138.04*C_Qq18 + 123.73*C_Qq18*C_Qq18 + 42.11*C_Qq38 +
                    126.81*C_Qq38*C_Qq38 + -19.47*C_Qu1 + 417.54*C_Qu1*C_Qu1 + 82.36*C_Qu8 + 
                    84.15*C_Qu8*C_Qu8 + -2.31*C_td1 + 226.87*C_td1*C_td1 + 54.34*C_td8 + 47.96*C_td8*C_td8 +
                    -13.53*C_tq1 + 638.83*C_tq1*C_tq1 + 129.53*C_tq8 + 131.02*C_tq8*C_tq8 + 9.16*C_tu1 + 
                    407.6*C_tu1*C_tu1 + 88.03*C_tu8 + 84.82*C_tu8*C_tu8)*(SM_tt_bin_1300_1500/tt_bin_1300_1500_Madgraph_NLO);
        } 
        else{ 
            return  SM_tt_bin_1300_1500 + (501.38*C_tG + -6.74*C_Qd1 + 47.09*C_Qd8 + 10.05*C_Qq11 +
                    5.8*C_Qq31 + 138.04*C_Qq18 + 42.11*C_Qq38 + -19.47*C_Qu1 + 82.36*C_Qu8 +
                    -2.31*C_td1 + 54.34*C_td8 + -13.53*C_tq1 + 129.53*C_tq8 + 9.16*C_tu1 +
                    88.03*C_tu8)*(SM_tt_bin_1300_1500/tt_bin_1300_1500_Madgraph_NLO);
        } 
        
    } else if(b_min == 1500 && b_max == 1700){
        
        double SM_tt_bin_1500_1700 = SM.getOptionalParameter("SM_sigma_tt_bin_1500_1700");
        
        //double tt_bin_1500_1700_Madgraph = 959.4;
        
        double tt_bin_1500_1700_Madgraph_NLO = 758.95;

        if(flag_Quadratic){ 
            /*
            return  ( tt_bin_1500_1700_Madgraph  +230.7*C_tG  +111.2*C_tG*C_tG  +44.12*C_Qd8  +39.87*C_Qd8*C_Qd8  
                    +179.6*C_Qd1*C_Qd1  +328.2*C_Qu1*C_Qu1  +90.01*C_Qu8  +72.96*C_Qu8*C_Qu8  +179.6*C_td1*C_td1  
                    +43.53*C_td8  +39.92*C_td8*C_td8  +506.6*C_tq1*C_tq1  +125.8*C_tq8  +112.7*C_tq8*C_tq8  
                    +328.1*C_tu1*C_tu1  +76.09*C_tu8  +72.99*C_tu8*C_tu8  +507.3*C_Qq11*C_Qq11  +507.0*C_Qq31*C_Qq31  
                    +117.8*C_Qq18  +112.7*C_Qq18*C_Qq18  +36.54*C_Qq38  +112.7*C_Qq38*C_Qq38  +6.474*C_tG*C_Qd8  
                    +11.35*C_tG*C_Qu8  +6.107*C_tG*C_td8  +17.93*C_tG*C_tq8  +11.38*C_tG*C_tu8  +18.0*C_tG*C_Qq18  
                    +5.291*C_tG*C_Qq38  +2.847*C_Qd8*C_td8  +12.85*C_Qd1*C_td1  +23.24*C_Qu1*C_tu1  +5.081*C_Qu8*C_tu8  
                    +35.61*C_tq1*C_Qq11  +10.65*C_tq1*C_Qq31  +7.878*C_tq8*C_Qq18  +2.311*C_tq8*C_Qq38  +298.0*C_Qq11*C_Qq31  
                    +66.6*C_Qq18*C_Qq38 )*(SM_tt_bin_1500_1700/tt_bin_1500_1700_Madgraph);
            */
            return SM_tt_bin_1500_1700 + (207.41*C_tG + 106.34*C_tG*C_tG + 1.67*C_Qd1 + 176.27*C_Qd1*C_Qd1 +
                    28.43*C_Qd8 + 38.36*C_Qd8*C_Qd8 + 5.02*C_Qq11 + 502.47*C_Qq11*C_Qq11 + -3.57*C_Qq31 +
                    497.51*C_Qq31*C_Qq31 + 82.76*C_Qq18 + 98.99*C_Qq18*C_Qq18 + 23.66*C_Qq38 +
                    105.82*C_Qq38*C_Qq38 + -9.44*C_Qu1 + 323.9*C_Qu1*C_Qu1 + 50.11*C_Qu8 + 70.11*C_Qu8*C_Qu8 +
                    -0.38*C_td1 + 173.06*C_td1*C_td1 + 29.0*C_td8 + 35.79*C_td8*C_td8 + -13.15*C_tq1 +
                    495.38*C_tq1*C_tq1 + 76.73*C_tq8 + 104.3*C_tq8*C_tq8 + 5.37*C_tu1 + 324.83*C_tu1*C_tu1 +
                    54.11*C_tu8 + 65.72*C_tu8*C_tu8)*(SM_tt_bin_1500_1700/tt_bin_1500_1700_Madgraph_NLO);
        } 
        else{ 
            return  SM_tt_bin_1500_1700 + (207.41*C_tG + 1.67*C_Qd1 + 28.43*C_Qd8 + 5.02*C_Qq11 + 
                    -3.57*C_Qq31 + 82.76*C_Qq18 + 23.66*C_Qq38 + -9.44*C_Qu1 + 50.11*C_Qu8 + 
                    -0.38*C_td1 + 29.0*C_td8 + -13.15*C_tq1 + 76.73*C_tq8 + 5.37*C_tu1 + 
                    54.11*C_tu8)*(SM_tt_bin_1500_1700/tt_bin_1500_1700_Madgraph_NLO);
        } 

        
    } else if(b_min == 1700 && b_max == 2000){
        
        double SM_tt_bin_1700_2000 = SM.getOptionalParameter("SM_sigma_tt_bin_1700_2000");

        //double tt_bin_1700_2000_Madgraph = 559.5;
        
        double tt_bin_1700_2000_Madgraph_NLO = 429.85;
        
        
        if(flag_Quadratic){ 
            /*
            return  ( tt_bin_1700_2000_Madgraph  +134.9*C_tG  +82.39*C_tG*C_tG  +37.57*C_Qd8  +46.9*C_Qd8*C_Qd8  
                    +211.1*C_Qd1*C_Qd1  +401.7*C_Qu1*C_Qu1  +71.21*C_Qu8  +89.28*C_Qu8*C_Qu8  +211.1*C_td1*C_td1  
                    +37.54*C_td8  +46.96*C_td8*C_td8  +612.1*C_tq1*C_tq1  +112.0*C_tq8  +136.1*C_tq8*C_tq8  
                    +401.4*C_tu1*C_tu1  +74.73*C_tu8  +89.17*C_tu8*C_tu8  +611.9*C_Qq11*C_Qq11  +612.2*C_Qq31*C_Qq31  
                    +112.6*C_Qq18  +135.9*C_Qq18*C_Qq18  +36.28*C_Qq38  +136.0*C_Qq38*C_Qq38  +5.761*C_tG*C_Qd8  
                    +10.62*C_tG*C_Qu8  +5.431*C_tG*C_td8  +16.26*C_tG*C_tq8  +11.14*C_tG*C_tu8  +16.89*C_tG*C_Qq18  
                    +5.017*C_tG*C_Qq38  +2.532*C_Qd8*C_td8  +10.6*C_Qd1*C_td1  +21.98*C_Qu1*C_tu1  +4.876*C_Qu8*C_tu8  
                    +32.73*C_tq1*C_Qq11  +10.64*C_tq1*C_Qq31  +7.146*C_tq8*C_Qq18  +2.282*C_tq8*C_Qq38  +383.2*C_Qq11*C_Qq31  
                    +85.13*C_Qq18*C_Qq38 )*(SM_tt_bin_1700_2000/tt_bin_1700_2000_Madgraph);
            */
            return SM_tt_bin_1700_2000 + (116.19*C_tG + 71.81*C_tG*C_tG + -2.55*C_Qd1 + 185.11*C_Qd1*C_Qd1 +
                    20.52*C_Qd8 + 38.57*C_Qd8*C_Qd8 + 6.82*C_Qq11 + 545.12*C_Qq11*C_Qq11 + -0.53*C_Qq31 +
                    542.36*C_Qq31*C_Qq31 + 70.98*C_Qq18 + 110.22*C_Qq18*C_Qq18 + 23.56*C_Qq38 +
                    116.32*C_Qq38*C_Qq38 + -4.96*C_Qu1 + 359.33*C_Qu1*C_Qu1 + 43.27*C_Qu8 +
                    75.48*C_Qu8*C_Qu8 + 3.71*C_td1 + 183.79*C_td1*C_td1 + 22.99*C_td8 + 38.91*C_td8*C_td8 +
                    -9.25*C_tq1 + 546.28*C_tq1*C_tq1 + 66.1*C_tq8 + 112.87*C_tq8*C_tq8 + 7.52*C_tu1 +
                    360.94*C_tu1*C_tu1 + 46.63*C_tu8 + 74.67*C_tu8*C_tu8)*(SM_tt_bin_1700_2000/tt_bin_1700_2000_Madgraph_NLO);
        } 
        else{ 
            return  SM_tt_bin_1700_2000 + (116.19*C_tG + -2.55*C_Qd1 + 20.52*C_Qd8 + 6.82*C_Qq11 +
                    -0.53*C_Qq31 + 70.98*C_Qq18 + 23.56*C_Qq38 + -4.96*C_Qu1 + 43.27*C_Qu8 +
                    3.71*C_td1 + 22.99*C_td8 + -9.25*C_tq1 + 66.1*C_tq8 + 7.52*C_tu1 +
                    46.63*C_tu8)*(SM_tt_bin_1700_2000/tt_bin_1700_2000_Madgraph_NLO);
        } 
        
    } else if(b_min == 2000 && b_max == 2300){
        
        double SM_tt_bin_2000_2300 = SM.getOptionalParameter("SM_sigma_tt_bin_2000_2300");
        
        //double tt_bin_2000_2300_Madgraph = 212.5;
        
        double tt_bin_2000_2300_Madgraph_NLO = 141.79;

        
        if(flag_Quadratic){ 
            /*
            return  ( tt_bin_2000_2300_Madgraph  +50.67*C_tG  +37.06*C_tG*C_tG  +15.98*C_Qd8  +34.4*C_Qd8*C_Qd8  
                    +154.9*C_Qd1*C_Qd1  +310.6*C_Qu1*C_Qu1  +45.02*C_Qu8  +68.98*C_Qu8*C_Qu8  +154.9*C_td1*C_td1  
                    +20.86*C_td8  +34.43*C_td8*C_td8  +465.0*C_tq1*C_tq1  +63.54*C_tq8  +103.3*C_tq8*C_tq8  
                    +310.9*C_tu1*C_tu1  +38.3*C_tu8  +69.09*C_tu8*C_tu8  +464.9*C_Qq11*C_Qq11  +465.1*C_Qq31*C_Qq31  
                    +63.34*C_Qq18  +103.3*C_Qq18*C_Qq18  +19.7*C_Qq38  +103.3*C_Qq38*C_Qq38  +3.111*C_tG*C_Qd8  
                    +6.053*C_tG*C_Qu8  +2.987*C_tG*C_td8  +8.916*C_tG*C_tq8  +5.868*C_tG*C_tu8  +8.867*C_tG*C_Qq18  
                    +3.049*C_tG*C_Qq38  +1.349*C_Qd8*C_td8  +6.249*C_Qd1*C_td1  +11.01*C_Qu1*C_tu1  +2.834*C_Qu8*C_tu8  
                    +18.39*C_tq1*C_Qq11  +6.336*C_tq1*C_Qq31  +4.182*C_tq8*C_Qq18  +1.447*C_tq8*C_Qq38  
                    +312.3*C_Qq11*C_Qq31  +69.2*C_Qq18*C_Qq38 )*(SM_tt_bin_2000_2300/tt_bin_2000_2300_Madgraph);
             */
            return SM_tt_bin_2000_2300 + (38.68*C_tG + 29.11*C_tG*C_tG + -1.98*C_Qd1 + 119.61*C_Qd1*C_Qd1 +
                    10.8*C_Qd8 + 25.82*C_Qd8*C_Qd8 + 0.31*C_Qq11 + 367.4*C_Qq11*C_Qq11 + -1.27*C_Qq31 +
                    366.28*C_Qq31*C_Qq31 + 35.82*C_Qq18 + 77.06*C_Qq18*C_Qq18 + 11.57*C_Qq38 +
                    79.11*C_Qq38*C_Qq38 + -4.81*C_Qu1 + 247.55*C_Qu1*C_Qu1 + 20.68*C_Qu8 + 52.83*C_Qu8*C_Qu8 +
                    1.35*C_td1 + 120.15*C_td1*C_td1 + 11.74*C_td8 + 25.86*C_td8*C_td8 + -3.42*C_tq1 +
                    367.15*C_tq1*C_tq1 + 32.29*C_tq8 + 77.87*C_tq8*C_tq8 + 1.84*C_tu1 + 248.49*C_tu1*C_tu1 +
                    23.7*C_tu8 + 52.57*C_tu8*C_tu8)*(SM_tt_bin_2000_2300/tt_bin_2000_2300_Madgraph_NLO);
        } 
        else{ 
            return  SM_tt_bin_2000_2300 + (38.68*C_tG + -1.98*C_Qd1 + 10.8*C_Qd8 + 0.31*C_Qq11 + 
                    -1.27*C_Qq31 + 35.82*C_Qq18 + 11.57*C_Qq38 + -4.81*C_Qu1 + 20.68*C_Qu8 + 
                    1.35*C_td1 + 11.74*C_td8 + -3.42*C_tq1 + 32.29*C_tq8 + 1.84*C_tu1 +
                    23.7*C_tu8)*(SM_tt_bin_2000_2300/tt_bin_2000_2300_Madgraph_NLO);
        } 
        
    } else if(b_min == 2300 && b_max == 2600){
        
        
        double SM_tt_bin_2300_2600 = SM.getOptionalParameter("SM_sigma_tt_bin_2300_2600");
        
        //double tt_bin_2300_2600_Madgraph = 89.87;
        
        double tt_bin_2300_2600_Madgraph_NLO = 51.27;

        
        if(flag_Quadratic){ 
            
            /*
            return  ( tt_bin_2300_2600_Madgraph  +21.46*C_tG  +17.57*C_tG*C_tG  +12.06*C_Qd8  +24.81*C_Qd8*C_Qd8  
                    +111.5*C_Qd1*C_Qd1  +237.8*C_Qu1*C_Qu1  +20.62*C_Qu8  +52.81*C_Qu8*C_Qu8  +111.5*C_td1*C_td1  
                    +9.636*C_td8  +24.82*C_td8*C_td8  +349.0*C_tq1*C_tq1  +35.87*C_tq8  +77.54*C_tq8*C_tq8  
                    +237.6*C_tu1*C_tu1  +22.63*C_tu8  +52.89*C_tu8*C_tu8  +348.9*C_Qq11*C_Qq11  +348.8*C_Qq31*C_Qq31  
                    +37.08*C_Qq18  +77.54*C_Qq18*C_Qq18  +11.79*C_Qq38  +77.54*C_Qq38*C_Qq38  +1.661*C_tG*C_Qd8  
                    +3.426*C_tG*C_Qu8  +1.553*C_tG*C_td8  +5.365*C_tG*C_tq8  +3.309*C_tG*C_tu8  +5.352*C_tG*C_Qq18  
                    +1.921*C_tG*C_Qq38  +0.6809*C_Qd8*C_td8  +3.425*C_Qd1*C_td1  +7.079*C_Qu1*C_tu1  +1.383*C_Qu8*C_tu8
                    +10.11*C_tq1*C_Qq11  +3.626*C_tq1*C_Qq31  +2.383*C_tq8*C_Qq18  +0.7835*C_tq8*C_Qq38  
                    +252.6*C_Qq11*C_Qq31  +56.29*C_Qq18*C_Qq38 )*(SM_tt_bin_2300_2600/tt_bin_2300_2600_Madgraph);
            */
        
            return  SM_tt_bin_2300_2600 + (13.87*C_tG + 12.5*C_tG*C_tG + -0.82*C_Qd1 + 
                    75.62*C_Qd1*C_Qd1 + 5.16*C_Qd8 + 16.41*C_Qd8*C_Qd8 + 4.42*C_Qq11 + 
                    248.08*C_Qq11*C_Qq11 + 1.66*C_Qq31 + 247.66*C_Qq31*C_Qq31 + 19.12*C_Qq18 +
                    52.25*C_Qq18*C_Qq18 + 7.2*C_Qq38 + 54.16*C_Qq38*C_Qq38 + -0.36*C_Qu1 + 
                    170.66*C_Qu1*C_Qu1 + 11.06*C_Qu8 + 36.43*C_Qu8*C_Qu8 + 1.03*C_td1 + 
                    75.8*C_td1*C_td1 + 6.05*C_td8 + 16.5*C_td8*C_td8 + -4.06*C_tq1 + 
                    246.9*C_tq1*C_tq1 + 16.17*C_tq8 + 52.21*C_tq8*C_tq8 + 1.47*C_tu1 + 
                    170.34*C_tu1*C_tu1 + 13.06*C_tu8 + 36.81*C_tu8*C_tu8)*(SM_tt_bin_2300_2600/tt_bin_2300_2600_Madgraph_NLO);
        } 
        else{ 
            return  SM_tt_bin_2300_2600 + (-0.82*C_Qd1 + 5.16*C_Qd8 + 4.42*C_Qq11 + 1.66*C_Qq31 + 19.12*C_Qq18 + 
                    7.2*C_Qq38 + -0.36*C_Qu1 + 11.06*C_Qu8 + 1.03*C_td1 + 6.05*C_td8 + -4.06*C_tq1 + 16.17*C_tq8 + 
                    1.47*C_tu1 + 13.06*C_tu8 )*(SM_tt_bin_2300_2600/tt_bin_2300_2600_Madgraph_NLO);
        } 
        
    }   else if(b_min == 2300 && b_max == 3500){
        
        double SM_tt_bin_2300_3500 = SM.getOptionalParameter("SM_sigma_tt_bin_2300_3500");
        
        double tt_bin_2300_3500_Madgraph_NLO = 82.53;


        if(flag_Quadratic){ 
            return  SM_tt_bin_2300_3500 + (22.25*C_tG + 22.87*C_tG*C_tG + -1.98*C_Qd1 + 186.86*C_Qd1*C_Qd1 +
                    9.73*C_Qd8 + 40.71*C_Qd8*C_Qd8 + 5.73*C_Qq11 + 711.85*C_Qq11*C_Qq11 + 2.07*C_Qq31 + 
                    712.37*C_Qq31*C_Qq31 + 41.11*C_Qq18 + 155.75*C_Qq18*C_Qq18 + 18.11*C_Qq38 + 
                    161.13*C_Qq38*C_Qq38 + -4.69*C_Qu1 + 525.85*C_Qu1*C_Qu1 + 26.14*C_Qu8 + 
                    113.82*C_Qu8*C_Qu8 + 1.15*C_td1 + 186.44*C_td1*C_td1 + 11.37*C_td8 + 41.9*C_td8*C_td8 +
                    -7.8*C_tq1 + 712.56*C_tq1*C_tq1 + 36.69*C_tq8 + 154.32*C_tq8*C_tq8 + 3.35*C_tu1 + 
                    525.51*C_tu1*C_tu1 + 29.48*C_tu8 + 116.13*C_tu8*C_tu8)*(SM_tt_bin_2300_3500/tt_bin_2300_3500_Madgraph_NLO);
        } 
        else{ 
            return  SM_tt_bin_2300_3500 + (22.25*C_tG + -1.98*C_Qd1 + 9.73*C_Qd8 + 5.73*C_Qq11 +
                    2.07*C_Qq31 + 41.11*C_Qq18 + 18.11*C_Qq38 + -4.69*C_Qu1 + 26.14*C_Qu8 + 
                    1.15*C_td1 + 11.37*C_td8 + -7.8*C_tq1 + 36.69*C_tq8 + 3.35*C_tu1 + 
                    29.48*C_tu8)*(SM_tt_bin_2300_3500/tt_bin_2300_3500_Madgraph_NLO);
        } 
        
    }  else if(b_min == 2600 && b_max == 3000){

        
        double SM_tt_bin_2600_3000 = SM.getOptionalParameter("SM_sigma_tt_bin_2600_3000");
        
        //double tt_bin_2600_3000_Madgraph = 39.36;
        
        double tt_bin_2600_3000_Madgraph_NLO = 23.0;


        if(flag_Quadratic){ 
            /*
            return  ( tt_bin_2600_3000_Madgraph  +9.427*C_tG  +10.55*C_tG*C_tG  +7.978*C_Qd8  +22.07*C_Qd8*C_Qd8  
                    +99.44*C_Qd1*C_Qd1  +230.0*C_Qu1*C_Qu1  +16.83*C_Qu8  +51.09*C_Qu8*C_Qu8  +99.35*C_td1*C_td1  
                    +7.156*C_td8  +22.12*C_td8*C_td8  +329.0*C_tq1*C_tq1  +28.47*C_tq8  +73.12*C_tq8*C_tq8  
                    +229.9*C_tu1*C_tu1  +19.59*C_tu8  +51.07*C_tu8*C_tu8  +329.2*C_Qq11*C_Qq11  +329.3*C_Qq31*C_Qq31  
                    +23.28*C_Qq18  +73.19*C_Qq18*C_Qq18  +10.36*C_Qq38  +73.14*C_Qq38*C_Qq38  +1.314*C_tG*C_Qd8  
                    +2.622*C_tG*C_Qu8  +0.9288*C_tG*C_td8  +3.89*C_tG*C_tq8  +2.703*C_tG*C_tu8  +3.908*C_tG*C_Qq18  
                    +1.633*C_tG*C_Qq38  +0.4681*C_Qd8*C_td8  +2.099*C_Qd1*C_td1  +5.291*C_Qu1*C_tu1  +1.25*C_Qu8*C_tu8  
                    +7.524*C_tq1*C_Qq11  +2.615*C_tq1*C_Qq31  +1.637*C_tq8*C_Qq18  +0.7179*C_tq8*C_Qq38  
                    +260.9*C_Qq11*C_Qq31  +57.97*C_Qq18*C_Qq38 )*(SM_tt_bin_2600_3000/tt_bin_2600_3000_Madgraph);
            */
            
            return  SM_tt_bin_2600_3000 + (6.19*C_tG + 6.69*C_tG*C_tG + -0.45*C_Qd1 + 57.77*C_Qd1*C_Qd1 + 
                    2.83*C_Qd8 + 12.58*C_Qd8*C_Qd8 + 2.02*C_Qq11 + 203.9*C_Qq11*C_Qq11 + -0.01*C_Qq31 + 
                    204.16*C_Qq31*C_Qq31 + 12.27*C_Qq18 + 44.12*C_Qq18*C_Qq18 + 5.27*C_Qq38 + 45.54*C_Qq38*C_Qq38 +
                    -1.73*C_Qu1 + 147.23*C_Qu1*C_Qu1 + 7.74*C_Qu8 + 31.62*C_Qu8*C_Qu8 + 0.26*C_td1 + 
                    58.11*C_td1*C_td1 + 3.44*C_td8 + 12.95*C_td8*C_td8 + -2.52*C_tq1 + 204.61*C_tq1*C_tq1 + 
                    11.2*C_tq8 + 43.89*C_tq8*C_tq8 + 1.82*C_tu1 + 147.52*C_tu1*C_tu1 + 8.61*C_tu8 + 
                    31.93*C_tu8*C_tu8)*(SM_tt_bin_2600_3000/tt_bin_2600_3000_Madgraph_NLO);
            
        } 
        else{ 
            return  SM_tt_bin_2600_3000 + (6.19*C_tG + -0.45*C_Qd1 + 2.83*C_Qd8 + 2.02*C_Qq11 + -0.01*C_Qq31 + 
                    12.27*C_Qq18 + 5.27*C_Qq38 + -1.73*C_Qu1 + 7.74*C_Qu8 + 0.26*C_td1 + 3.44*C_td8 + -2.52*C_tq1 +
                    11.2*C_tq8 + 1.82*C_tu1 + 8.61*C_tu8)*(SM_tt_bin_2600_3000/tt_bin_2600_3000_Madgraph_NLO);
        } 
        
    }  else if(b_min == 3000 && b_max == 3500){

        double SM_tt_bin_3000_3500 = SM.getOptionalParameter("SM_sigma_tt_bin_3000_3500");

        //double tt_bin_3000_3500_Madgraph = 13.76;

        double tt_bin_3000_3500_Madgraph_NLO = 7.71;
        
        if(flag_Quadratic){ 
            /*
            return  SM_tt_bin_3000_3500 + (3.303*C_tG  +4.991*C_tG*C_tG  +5.038*C_Qd8  +15.85*C_Qd8*C_Qd8  
                    +71.36*C_Qd1*C_Qd1  +187.2*C_Qu1*C_Qu1  +10.95*C_Qu8  +41.64*C_Qu8*C_Qu8  +71.4*C_td1*C_td1  
                    +4.086*C_td8  +15.9*C_td8*C_td8  +258.6*C_tq1*C_tq1  +14.98*C_tq8  +57.42*C_tq8*C_tq8  
                    +187.2*C_tu1*C_tu1  +10.9*C_tu8  +41.59*C_tu8*C_tu8  +258.4*C_Qq11*C_Qq11  +258.4*C_Qq31*C_Qq31
                    +15.04*C_Qq18  +57.47*C_Qq18*C_Qq18  +8.259*C_Qq38  +57.37*C_Qq38*C_Qq38  +0.6228*C_tG*C_Qd8  
                    +1.473*C_tG*C_Qu8  +0.4583*C_tG*C_td8  +2.507*C_tG*C_tq8  +1.468*C_tG*C_tu8  +2.297*C_tG*C_Qq18 
                    +1.469*C_tG*C_Qq38  +0.2076*C_Qd8*C_td8  +1.409*C_Qd1*C_td1  +3.14*C_Qu1*C_tu1  +0.77*C_Qu8*C_tu8  
                    +4.855*C_tq1*C_Qq11  +2.366*C_tq1*C_Qq31  +0.8384*C_tq8*C_Qq18  +0.5345*C_tq8*C_Qq38  
                    +232.4*C_Qq11*C_Qq31  +51.59*C_Qq18*C_Qq38)*(SM_tt_bin_3000_3500/tt_bin_3000_3500_Madgraph);
            */
            
            return  SM_tt_bin_3000_3500 + (2.07*C_tG + 2.8*C_tG*C_tG + -0.41*C_Qd1 + 33.75*C_Qd1*C_Qd1 + 1.31*C_Qd8 +
                    7.43*C_Qd8*C_Qd8 + 0.24*C_Qq11 + 136.58*C_Qq11*C_Qq11 + 0.26*C_Qq31 + 136.62*C_Qq31*C_Qq31 + 
                    6.03*C_Qq18 + 30.12*C_Qq18*C_Qq18 + 2.89*C_Qq38 + 31.16*C_Qq38*C_Qq38 + -0.47*C_Qu1 + 103.22*C_Qu1*C_Qu1 +
                    3.99*C_Qu8 + 22.39*C_Qu8*C_Qu8 + 0.19*C_td1 + 33.8*C_td1*C_td1 + 1.5*C_td8 + 7.77*C_td8*C_td8 + 
                    -0.7*C_tq1 + 136.99*C_tq1*C_tq1 + 5.4*C_tq8 + 29.6*C_tq8*C_tq8 + 0.85*C_tu1 + 103.21*C_tu1*C_tu1 +
                    4.63*C_tu8 + 22.83*C_tu8*C_tu8)*(SM_tt_bin_3000_3500/tt_bin_3000_3500_Madgraph_NLO);
            
        } 
        else{ 
            return  SM_tt_bin_3000_3500 + (2.07*C_tG + -0.41*C_Qd1 + 1.31*C_Qd8 + 0.24*C_Qq11 + 0.26*C_Qq31 + 6.03*C_Qq18 +
                    2.89*C_Qq38 + -0.47*C_Qu1 + 3.99*C_Qu8 + 0.19*C_td1 + 1.5*C_td8 + -0.7*C_tq1 + 5.4*C_tq8 + 0.85*C_tu1 +
                    4.63*C_tu8)*(SM_tt_bin_3000_3500/tt_bin_3000_3500_Madgraph_NLO);
        } 
        
    }  else if(b_min == 3500 && b_max == 4000){
        ////////////////////////////////
        //NEEDS TO BE UPDATED!!!!!!!!!//
        ////////////////////////////////
        double SM_tt_bin_3500_4000 = SM.getOptionalParameter("SM_sigma_tt_bin_3500_4000");
        
        //double tt_bin_3500_4000_Madgraph = 3.864;
        
        double tt_bin_3500_4000_Madgraph_NLO = 1.87;


        if(flag_Quadratic){ 
            /*
            return  SM_tt_bin_3500_4000 + (+0.9714*C_tG  +1.769*C_tG*C_tG  +1.842*C_Qd8  +8.195*C_Qd8*C_Qd8  
                    +36.82*C_Qd1*C_Qd1  +113.3*C_Qu1*C_Qu1  +4.798*C_Qu8  +25.16*C_Qu8*C_Qu8  +36.88*C_td1*C_td1  
                    +1.512*C_td8  +8.172*C_td8*C_td8  +149.9*C_tq1*C_tq1  +6.922*C_tq8  +33.32*C_tq8*C_tq8  
                    +113.0*C_tu1*C_tu1  +5.41*C_tu8  +25.1*C_tu8*C_tu8  +150.0*C_Qq11*C_Qq11  +149.9*C_Qq31*C_Qq31  
                    +7.362*C_Qq18  +33.33*C_Qq18*C_Qq18  +4.131*C_Qq38  +33.35*C_Qq38*C_Qq38  +0.1559*C_tG*C_Qd8  
                    +0.712*C_tG*C_Qu8  +0.3765*C_tG*C_td8  +1.038*C_tG*C_tq8  +1.042*C_tG*C_tu8  +0.8589*C_tG*C_Qq18  
                    +0.2499*C_tG*C_Qq38  +0.1047*C_Qd8*C_td8  +0.532*C_Qd1*C_td1  +1.65*C_Qu1*C_tu1  +0.3337*C_Qu8*C_tu8
                    +1.465*C_tq1*C_Qq11  +0.9763*C_tq1*C_Qq31  +0.3631*C_tq8*C_Qq18  +0.2485*C_tq8*C_Qq38  
                    +152.7*C_Qq11*C_Qq31  +34.02*C_Qq18*C_Qq38 )*(SM_tt_bin_3500_4000/tt_bin_3500_4000_Madgraph);
            */ 
            return  SM_tt_bin_3500_4000 + (0.5*C_tG + 0.85*C_tG*C_tG + -0.07*C_Qd1 + 13.16*C_Qd1*C_Qd1 + 0.38*C_Qd8 +
                    2.94*C_Qd8*C_Qd8 + 0.46*C_Qq11 + 66.59*C_Qq11*C_Qq11 + 0.45*C_Qq31 + 66.62*C_Qq31*C_Qq31 + 2.18*C_Qq18 +
                    14.91*C_Qq18*C_Qq18 + 1.44*C_Qq38 + 15.49*C_Qq38*C_Qq38 + -0.17*C_Qu1 + 53.24*C_Qu1*C_Qu1 + 1.61*C_Qu8 +
                    11.67*C_Qu8*C_Qu8 + 0.04*C_td1 + 13.17*C_td1*C_td1 + 0.45*C_td8 + 3.13*C_td8*C_td8 + -0.46*C_tq1 +
                    66.38*C_tq1*C_tq1 + 1.97*C_tq8 + 14.52*C_tq8*C_tq8 + 0.26*C_tu1 + 53.41*C_tu1*C_tu1 + 1.8*C_tu8 +
                    12.11*C_tu8*C_tu8)*(SM_tt_bin_3500_4000/tt_bin_3500_4000_Madgraph_NLO);
            
        } 
        else{ 
            return  SM_tt_bin_3500_4000 + (0.5*C_tG + -0.07*C_Qd1 + 0.38*C_Qd8 + 0.46*C_Qq11 + 0.45*C_Qq31 + 2.18*C_Qq18 + 
                    1.44*C_Qq38 + -0.17*C_Qu1 + 1.61*C_Qu8 + 0.04*C_td1 + 0.45*C_td8 + -0.46*C_tq1 + 1.97*C_tq8 + 
                    0.26*C_tu1 + 1.8*C_tu8)*(SM_tt_bin_3500_4000/tt_bin_3500_4000_Madgraph_NLO);
        } 
        
    } else {
        throw std::runtime_error("\nERROR: Please specify a correct binning range for sigma_tt_diff_NLO.\n");
    }

}







/////////////////////////////////////// Charge Asymmetry ttbar ////////////////////////////////////////////////


charge_asymmetry_tt_diff_mtt_NLO::charge_asymmetry_tt_diff_mtt_NLO(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_charge_asymmetry_bin_mtt_0_500" << "SM_charge_asymmetry_deno_bin_mtt_0_500"
            << "SM_charge_asymmetry_bin_mtt_500_750" << "SM_charge_asymmetry_deno_bin_mtt_500_750" << "SM_charge_asymmetry_bin_mtt_750_1000" 
            << "SM_charge_asymmetry_deno_bin_mtt_750_1000" << "SM_charge_asymmetry_bin_mtt_1000_1500" << "SM_charge_asymmetry_deno_bin_mtt_1000_1500"
            << "SM_charge_asymmetry_bin_mtt_1500_3000" << "SM_charge_asymmetry_deno_bin_mtt_1500_3000" << "SM_charge_asymmetry_bin_mtt_1500_2000"
            << "SM_charge_asymmetry_deno_bin_mtt_1500_2000" << "SM_charge_asymmetry_bin_mtt_2000_2500" << "SM_charge_asymmetry_deno_bin_mtt_2000_2500" 
            << "SM_charge_asymmetry_bin_mtt_2500_3000" << "SM_charge_asymmetry_deno_bin_mtt_2500_3000");
    
}

double charge_asymmetry_tt_diff_mtt_NLO::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    bool   flag_Quadratic = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic(); 
    
    //double SM_Charge_Asymmetry_bin_tt_0_500 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_0_500();

    
    if(b_min == 0 && b_max == 500){
        
        double SM_charge_asymmetry_bin_mtt_0_500 = SM.getOptionalParameter("SM_charge_asymmetry_bin_mtt_0_500");
        double SM_charge_asymmetry_deno_bin_mtt_0_500 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_0_500");
        double SM_charge_asymmetry_num_bin_mtt_0_500 = SM_charge_asymmetry_bin_mtt_0_500*SM_charge_asymmetry_deno_bin_mtt_0_500;
        
        double SM_sigma_pos_bin_mtt_0_500 =0.5*(SM_charge_asymmetry_num_bin_mtt_0_500+SM_charge_asymmetry_deno_bin_mtt_0_500);
        double SM_sigma_neg_bin_mtt_0_500 =0.5*(SM_charge_asymmetry_deno_bin_mtt_0_500-SM_charge_asymmetry_num_bin_mtt_0_500);

        
        double sigma_pos_bin_mtt_0_500_madgraph_NLO = 206522.71; //fb
        double sigma_neg_bin_mtt_0_500_madgraph_NLO = 205736.11; //fb
        
        
        if(flag_Quadratic){
        
            
            double sigma_pos_bin_mtt_0_500_NP_quad = 10856.61*C_tG*C_tG + 186.94*C_Qd1*C_Qd1 + 81.37*C_Qd8*C_Qd8 + 655.08*C_Qq11*C_Qq11 +
            712.66*C_Qq31*C_Qq31 + 589.45*C_Qq18*C_Qq18 + 121.35*C_Qq38*C_Qq38 + 198.11*C_Qu1*C_Qu1 + 83.7*C_Qu8*C_Qu8 + 244.73*C_td1*C_td1 +
            -27.74*C_td8*C_td8 + 499.26*C_tq1*C_tq1 + 35.53*C_tq8*C_tq8 + 380.22*C_tu1*C_tu1 + 128.01*C_tu8*C_tu8;
            double sigma_neg_bin_mtt_0_500_NP_quad = 10844.05*C_tG*C_tG + 124.97*C_Qd1*C_Qd1 + -57.71*C_Qd8*C_Qd8 + 415.01*C_Qq11*C_Qq11 + 
            348.52*C_Qq31*C_Qq31 + 575.55*C_Qq18*C_Qq18 + 148.51*C_Qq38*C_Qq38 + 373.47*C_Qu1*C_Qu1 + 31.75*C_Qu8*C_Qu8 + 258.97*C_td1*C_td1 + 
            1.84*C_td8*C_td8 + 602.71*C_tq1*C_tq1 + 100.72*C_tq8*C_tq8 + 275.67*C_tu1*C_tu1 + -42.82*C_tu8*C_tu8;
            
            
            double sigma_pos_bin_mtt_0_500_NP_lin = 69255.55*C_tG + 86.03*C_Qd1 + 631.63*C_Qd8 + -171.25*C_Qq11 + 11.28*C_Qq31 + 977.78*C_Qq18 +
            383.52*C_Qq38 + 191.56*C_Qu1 + 1018.38*C_Qu8 + -191.45*C_td1 + 677.13*C_td8 + 273.38*C_tq1 + 1321.58*C_tq8 + -7.5*C_tu1 + 1119.97*C_tu8;
            double sigma_neg_bin_mtt_0_500_NP_lin = 68845.04*C_tG + 203.17*C_Qd1 + 727.79*C_Qd8 + -445.67*C_Qq11 + -122.59*C_Qq31 + 423.69*C_Qq18 + 
            72.21*C_Qq38 + -65.45*C_Qu1 + 760.75*C_Qu8 + -267.01*C_td1 + 508.03*C_td8 + 176.2*C_tq1 + 1589.47*C_tq8 + -336.04*C_tu1 + 610.15*C_tu8;
            
              
            double sigma_pos_bin_mtt_0_500_NP = sigma_pos_bin_mtt_0_500_NP_lin + sigma_pos_bin_mtt_0_500_NP_quad;
            double sigma_neg_bin_mtt_0_500_NP = sigma_neg_bin_mtt_0_500_NP_lin + sigma_neg_bin_mtt_0_500_NP_quad;
            

            double sigma_pos_bin_mtt_0_500_NP_Corrected = SM_sigma_pos_bin_mtt_0_500*sigma_pos_bin_mtt_0_500_NP/sigma_pos_bin_mtt_0_500_madgraph_NLO;
            double sigma_neg_bin_mtt_0_500_NP_Corrected = SM_sigma_neg_bin_mtt_0_500*sigma_neg_bin_mtt_0_500_NP/sigma_neg_bin_mtt_0_500_madgraph_NLO;
            
            double NP_charge_asymmetry_num_bin_mtt_0_500 = sigma_pos_bin_mtt_0_500_NP_Corrected - sigma_neg_bin_mtt_0_500_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_0_500 = sigma_pos_bin_mtt_0_500_NP_Corrected + sigma_neg_bin_mtt_0_500_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_0_500*(1+(NP_charge_asymmetry_num_bin_mtt_0_500-NP_charge_asymmetry_deno_bin_mtt_0_500)/SM_charge_asymmetry_deno_bin_mtt_0_500);            
            
            
            return SM_charge_asymmetry_bin_mtt_0_500 + (NP_charge_asymmetry_num_bin_mtt_0_500
                    - SM_charge_asymmetry_bin_mtt_0_500*NP_charge_asymmetry_deno_bin_mtt_0_500
                    )/SM_charge_asymmetry_deno_bin_mtt_0_500;
            
        }
        else{
            

            double sigma_pos_bin_mtt_0_500_NP = 69255.55*C_tG + 86.03*C_Qd1 + 631.63*C_Qd8 + -171.25*C_Qq11 + 11.28*C_Qq31 + 977.78*C_Qq18 +
            383.52*C_Qq38 + 191.56*C_Qu1 + 1018.38*C_Qu8 + -191.45*C_td1 + 677.13*C_td8 + 273.38*C_tq1 + 1321.58*C_tq8 + -7.5*C_tu1 + 1119.97*C_tu8;
            double sigma_neg_bin_mtt_0_500_NP = 68845.04*C_tG + 203.17*C_Qd1 + 727.79*C_Qd8 + -445.67*C_Qq11 + -122.59*C_Qq31 + 423.69*C_Qq18 + 
            72.21*C_Qq38 + -65.45*C_Qu1 + 760.75*C_Qu8 + -267.01*C_td1 + 508.03*C_td8 + 176.2*C_tq1 + 1589.47*C_tq8 + -336.04*C_tu1 + 610.15*C_tu8;
            
;
            
            double sigma_pos_bin_mtt_0_500_NP_Corrected = SM_sigma_pos_bin_mtt_0_500*sigma_pos_bin_mtt_0_500_NP/sigma_pos_bin_mtt_0_500_madgraph_NLO;
            double sigma_neg_bin_mtt_0_500_NP_Corrected = SM_sigma_neg_bin_mtt_0_500*sigma_neg_bin_mtt_0_500_NP/sigma_neg_bin_mtt_0_500_madgraph_NLO;
            
            double NP_charge_asymmetry_num_bin_mtt_0_500 = sigma_pos_bin_mtt_0_500_NP_Corrected - sigma_neg_bin_mtt_0_500_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_0_500 = sigma_pos_bin_mtt_0_500_NP_Corrected + sigma_neg_bin_mtt_0_500_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_0_500*(1+(NP_charge_asymmetry_num_bin_mtt_0_500-NP_charge_asymmetry_deno_bin_mtt_0_500)/SM_charge_asymmetry_deno_bin_mtt_0_500);            
        
            
            return SM_charge_asymmetry_bin_mtt_0_500 + (NP_charge_asymmetry_num_bin_mtt_0_500
                    - SM_charge_asymmetry_bin_mtt_0_500*NP_charge_asymmetry_deno_bin_mtt_0_500
                    )/SM_charge_asymmetry_deno_bin_mtt_0_500;
            
            
        }
        
    } else if(b_min == 500 && b_max == 750){
        
        double SM_charge_asymmetry_bin_mtt_500_750 = SM.getOptionalParameter("SM_charge_asymmetry_bin_mtt_500_750");
        double SM_charge_asymmetry_deno_bin_mtt_500_750 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_500_750");
        double SM_charge_asymmetry_num_bin_mtt_500_750 = SM_charge_asymmetry_bin_mtt_500_750*SM_charge_asymmetry_deno_bin_mtt_500_750;
        
        
        double SM_sigma_pos_bin_mtt_500_750 =0.5*(SM_charge_asymmetry_num_bin_mtt_500_750+SM_charge_asymmetry_deno_bin_mtt_500_750);
        double SM_sigma_neg_bin_mtt_500_750 =0.5*(SM_charge_asymmetry_deno_bin_mtt_500_750-SM_charge_asymmetry_num_bin_mtt_500_750);

        
        double sigma_pos_bin_mtt_500_750_madgraph_NLO = 98972.12;
        double sigma_neg_bin_mtt_500_750_madgraph_NLO = 98471.35;
        
        if(flag_Quadratic){
        
            
            double sigma_pos_bin_mtt_500_750_NP_quad = 5493.57*C_tG*C_tG + 315.65*C_Qd1*C_Qd1 + 52.9*C_Qd8*C_Qd8 + 1274.19*C_Qq11*C_Qq11 +
            1247.82*C_Qq31*C_Qq31 + 234.29*C_Qq18*C_Qq18 + 220.5*C_Qq38*C_Qq38 + 375.02*C_Qu1*C_Qu1 + 81.75*C_Qu8*C_Qu8 + 449.86*C_td1*C_td1 +
            67.51*C_td8*C_td8 + 652.68*C_tq1*C_tq1 + 147.74*C_tq8*C_tq8 + 776.43*C_tu1*C_tu1 + 140.43*C_tu8*C_tu8;
            double sigma_neg_bin_mtt_500_750_NP_quad = 5476.01*C_tG*C_tG + 419.52*C_Qd1*C_Qd1 + 77.02*C_Qd8*C_Qd8 + 616.75*C_Qq11*C_Qq11 + 
            617.3*C_Qq31*C_Qq31 + 141.47*C_Qq18*C_Qq18 + 76.31*C_Qq38*C_Qq38 + 751.03*C_Qu1*C_Qu1 + 118.22*C_Qu8*C_Qu8 + 276.88*C_td1*C_td1 +
            24.92*C_td8*C_td8 + 1237.5*C_tq1*C_tq1 + 200.26*C_tq8*C_tq8 + 359.71*C_tu1*C_tu1 + 52.72*C_tu8*C_tu8;
            

            
            double sigma_pos_bin_mtt_500_750_NP_lin = 28985.27*C_tG + 32.74*C_Qd1 + 382.94*C_Qd8 + -108.77*C_Qq11 + -48.35*C_Qq31 + 1456.01*C_Qq18 +
            341.82*C_Qq38 + 115.28*C_Qu1 + 464.61*C_Qu8 + -46.51*C_td1 + 563.4*C_td8 + 103.98*C_tq1 + 886.74*C_tq8 + -18.49*C_tu1 + 904.87*C_tu8;
            double sigma_neg_bin_mtt_500_750_NP_lin = 28804.73*C_tG + 48.33*C_Qd1 + 591.25*C_Qd8 + -189.14*C_Qq11 + -21.89*C_Qq31 + 737.03*C_Qq18 + 
            69.36*C_Qq38 + -81.18*C_Qu1 + 932.13*C_Qu8 + -82.69*C_td1 + 292.91*C_td8 + -53.01*C_tq1 + 1450.11*C_tq8 + -87.3*C_tu1 + 442.24*C_tu8;
            
            
              
            double sigma_pos_bin_mtt_500_750_NP = sigma_pos_bin_mtt_500_750_NP_lin + sigma_pos_bin_mtt_500_750_NP_quad;
            double sigma_neg_bin_mtt_500_750_NP = sigma_neg_bin_mtt_500_750_NP_lin + sigma_neg_bin_mtt_500_750_NP_quad;
            
            
            double sigma_pos_bin_mtt_500_750_NP_Corrected = SM_sigma_pos_bin_mtt_500_750*sigma_pos_bin_mtt_500_750_NP/sigma_pos_bin_mtt_500_750_madgraph_NLO;
            double sigma_neg_bin_mtt_500_750_NP_Corrected = SM_sigma_neg_bin_mtt_500_750*sigma_neg_bin_mtt_500_750_NP/sigma_neg_bin_mtt_500_750_madgraph_NLO;
            
            double NP_charge_asymmetry_num_bin_mtt_500_750 = sigma_pos_bin_mtt_500_750_NP_Corrected - sigma_neg_bin_mtt_500_750_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_500_750 = sigma_pos_bin_mtt_500_750_NP_Corrected + sigma_neg_bin_mtt_500_750_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_500_750*(1+(NP_charge_asymmetry_num_bin_mtt_500_750-NP_charge_asymmetry_deno_bin_mtt_500_750)/SM_charge_asymmetry_deno_bin_mtt_500_750);            
        
            
            return SM_charge_asymmetry_bin_mtt_500_750 + (NP_charge_asymmetry_num_bin_mtt_500_750
                    - SM_charge_asymmetry_bin_mtt_500_750*NP_charge_asymmetry_deno_bin_mtt_500_750
                    )/SM_charge_asymmetry_deno_bin_mtt_500_750;
            
            
        }
        else{
            

            double sigma_pos_bin_mtt_500_750_NP = 28985.27*C_tG + 32.74*C_Qd1 + 382.94*C_Qd8 + -108.77*C_Qq11 + -48.35*C_Qq31 + 1456.01*C_Qq18 +
            341.82*C_Qq38 + 115.28*C_Qu1 + 464.61*C_Qu8 + -46.51*C_td1 + 563.4*C_td8 + 103.98*C_tq1 + 886.74*C_tq8 + -18.49*C_tu1 + 904.87*C_tu8;
            double sigma_neg_bin_mtt_500_750_NP = 28804.73*C_tG + 48.33*C_Qd1 + 591.25*C_Qd8 + -189.14*C_Qq11 + -21.89*C_Qq31 + 737.03*C_Qq18 + 
            69.36*C_Qq38 + -81.18*C_Qu1 + 932.13*C_Qu8 + -82.69*C_td1 + 292.91*C_td8 + -53.01*C_tq1 + 1450.11*C_tq8 + -87.3*C_tu1 + 442.24*C_tu8;
            
            
            double sigma_pos_bin_mtt_500_750_NP_Corrected = SM_sigma_pos_bin_mtt_500_750*sigma_pos_bin_mtt_500_750_NP/sigma_pos_bin_mtt_500_750_madgraph_NLO;
            double sigma_neg_bin_mtt_500_750_NP_Corrected = SM_sigma_neg_bin_mtt_500_750*sigma_neg_bin_mtt_500_750_NP/sigma_neg_bin_mtt_500_750_madgraph_NLO;
            
            double NP_charge_asymmetry_num_bin_mtt_500_750 = sigma_pos_bin_mtt_500_750_NP_Corrected - sigma_neg_bin_mtt_500_750_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_500_750 = sigma_pos_bin_mtt_500_750_NP_Corrected + sigma_neg_bin_mtt_500_750_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_500_750*(1+(NP_charge_asymmetry_num_bin_mtt_500_750-NP_charge_asymmetry_deno_bin_mtt_500_750)/SM_charge_asymmetry_deno_bin_mtt_500_750);            
        
            
            return SM_charge_asymmetry_bin_mtt_500_750 + (NP_charge_asymmetry_num_bin_mtt_500_750
                    - SM_charge_asymmetry_bin_mtt_500_750*NP_charge_asymmetry_deno_bin_mtt_500_750
                    )/SM_charge_asymmetry_deno_bin_mtt_500_750;
            
            
            
        }
        
    }  else if(b_min == 750 && b_max == 1000){
        
        
        double SM_charge_asymmetry_deno_bin_mtt_750_1000 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_750_1000");
        double SM_charge_asymmetry_bin_mtt_750_1000 = SM.getOptionalParameter("SM_charge_asymmetry_bin_mtt_750_1000");
        double SM_charge_asymmetry_num_bin_mtt_750_1000 = SM_charge_asymmetry_bin_mtt_750_1000*SM_charge_asymmetry_deno_bin_mtt_750_1000;
        
        double SM_sigma_pos_bin_mtt_750_1000 =0.5*(SM_charge_asymmetry_num_bin_mtt_750_1000+SM_charge_asymmetry_deno_bin_mtt_750_1000);
        double SM_sigma_neg_bin_mtt_750_1000 =0.5*(SM_charge_asymmetry_deno_bin_mtt_750_1000-SM_charge_asymmetry_num_bin_mtt_750_1000);

        double sigma_pos_bin_mtt_750_1000_madgraph_NLO = 18609.62;
        double sigma_neg_bin_mtt_750_1000_madgraph_NLO = 18459.50;
        
        
        if(flag_Quadratic){
            
            
            double sigma_pos_bin_mtt_750_1000_NP_quad = 1339.65*C_tG*C_tG + 210.5*C_Qd1*C_Qd1 + 43.98*C_Qd8*C_Qd8 + 1017.96*C_Qq11*C_Qq11 + 
            1008.2*C_Qq31*C_Qq31 + 189.09*C_Qq18*C_Qq18 + 201.45*C_Qq38*C_Qq38 + 275.07*C_Qu1*C_Qu1 + 55.95*C_Qu8*C_Qu8 + 378.54*C_td1*C_td1 + 
            66.48*C_td8*C_td8 + 458.34*C_tq1*C_tq1 + 108.63*C_tq8*C_tq8 + 655.25*C_tu1*C_tu1 + 134.08*C_tu8*C_tu8;
            double sigma_neg_bin_mtt_750_1000_NP_quad = 1343.97*C_tG*C_tG + 370.22*C_Qd1*C_Qd1 + 70.13*C_Qd8*C_Qd8 + 481.64*C_Qq11*C_Qq11 + 
            481.3*C_Qq31*C_Qq31 + 94.58*C_Qq18*C_Qq18 + 96.39*C_Qq38*C_Qq38 + 653.11*C_Qu1*C_Qu1 + 137.27*C_Qu8*C_Qu8 + 206.12*C_td1*C_td1 + 
            54.0*C_td8*C_td8 + 1028.88*C_tq1*C_tq1 + 202.31*C_tq8*C_tq8 + 279.84*C_tu1*C_tu1 + 53.02*C_tu8*C_tu8;
            
            
            
        
            double sigma_pos_bin_mtt_750_1000_NP_lin = 5266.33*C_tG + 17.3*C_Qd1 + 124.0*C_Qd8 + 20.27*C_Qq11 + 7.09*C_Qq31 + 594.31*C_Qq18 + 
            162.54*C_Qq38 + 4.65*C_Qu1 + 167.93*C_Qu8 + 3.88*C_td1 + 210.36*C_td8 + 10.89*C_tq1 + 280.34*C_tq8 + 7.1*C_tu1 + 369.73*C_tu8;
            double sigma_neg_bin_mtt_750_1000_NP_lin = 5217.21*C_tG + -13.12*C_Qd1 + 209.09*C_Qd8 + -27.46*C_Qq11 + -7.0*C_Qq31 + 255.72*C_Qq18 +
            51.39*C_Qq38 + -36.2*C_Qu1 + 356.72*C_Qu8 + -8.9*C_td1 + 120.85*C_td8 + -42.63*C_tq1 + 568.51*C_tq8 + -21.45*C_tu1 + 151.18*C_tu8;
            
              
            double sigma_pos_bin_mtt_750_1000_NP = sigma_pos_bin_mtt_750_1000_NP_lin + sigma_pos_bin_mtt_750_1000_NP_quad;
            double sigma_neg_bin_mtt_750_1000_NP = sigma_neg_bin_mtt_750_1000_NP_lin + sigma_neg_bin_mtt_750_1000_NP_quad;
            
            
            double sigma_pos_bin_mtt_750_1000_NP_Corrected = SM_sigma_pos_bin_mtt_750_1000*sigma_pos_bin_mtt_750_1000_NP/sigma_pos_bin_mtt_750_1000_madgraph_NLO;
            double sigma_neg_bin_mtt_750_1000_NP_Corrected = SM_sigma_neg_bin_mtt_750_1000*sigma_neg_bin_mtt_750_1000_NP/sigma_neg_bin_mtt_750_1000_madgraph_NLO;
            
            double NP_charge_asymmetry_num_bin_mtt_750_1000 = sigma_pos_bin_mtt_750_1000_NP_Corrected - sigma_neg_bin_mtt_750_1000_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_750_1000 = sigma_pos_bin_mtt_750_1000_NP_Corrected + sigma_neg_bin_mtt_750_1000_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_750_1000*(1+(NP_charge_asymmetry_num_bin_mtt_750_1000-NP_charge_asymmetry_deno_bin_mtt_750_1000)/SM_charge_asymmetry_deno_bin_mtt_750_1000);            
        
            return SM_charge_asymmetry_bin_mtt_750_1000 + (NP_charge_asymmetry_num_bin_mtt_750_1000
                    - SM_charge_asymmetry_bin_mtt_750_1000*NP_charge_asymmetry_deno_bin_mtt_750_1000
                    )/SM_charge_asymmetry_deno_bin_mtt_750_1000;
            
        }
        else{
            

            double sigma_pos_bin_mtt_750_1000_NP = 5266.33*C_tG + 17.3*C_Qd1 + 124.0*C_Qd8 + 20.27*C_Qq11 + 7.09*C_Qq31 + 594.31*C_Qq18 + 
            162.54*C_Qq38 + 4.65*C_Qu1 + 167.93*C_Qu8 + 3.88*C_td1 + 210.36*C_td8 + 10.89*C_tq1 + 280.34*C_tq8 + 7.1*C_tu1 + 369.73*C_tu8;
            double sigma_neg_bin_mtt_750_1000_NP = 5217.21*C_tG + -13.12*C_Qd1 + 209.09*C_Qd8 + -27.46*C_Qq11 + -7.0*C_Qq31 + 255.72*C_Qq18 +
            51.39*C_Qq38 + -36.2*C_Qu1 + 356.72*C_Qu8 + -8.9*C_td1 + 120.85*C_td8 + -42.63*C_tq1 + 568.51*C_tq8 + -21.45*C_tu1 + 151.18*C_tu8;
            
            
            double sigma_pos_bin_mtt_750_1000_NP_Corrected = SM_sigma_pos_bin_mtt_750_1000*sigma_pos_bin_mtt_750_1000_NP/sigma_pos_bin_mtt_750_1000_madgraph_NLO;
            double sigma_neg_bin_mtt_750_1000_NP_Corrected = SM_sigma_neg_bin_mtt_750_1000*sigma_neg_bin_mtt_750_1000_NP/sigma_neg_bin_mtt_750_1000_madgraph_NLO;
            
            double NP_charge_asymmetry_num_bin_mtt_750_1000 = sigma_pos_bin_mtt_750_1000_NP_Corrected - sigma_neg_bin_mtt_750_1000_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_750_1000 = sigma_pos_bin_mtt_750_1000_NP_Corrected + sigma_neg_bin_mtt_750_1000_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_750_1000*(1+(NP_charge_asymmetry_num_bin_mtt_750_1000-NP_charge_asymmetry_deno_bin_mtt_750_1000)/SM_charge_asymmetry_deno_bin_mtt_750_1000);            
        
            return SM_charge_asymmetry_bin_mtt_750_1000 + (NP_charge_asymmetry_num_bin_mtt_750_1000
                    - SM_charge_asymmetry_bin_mtt_750_1000*NP_charge_asymmetry_deno_bin_mtt_750_1000
                    )/SM_charge_asymmetry_deno_bin_mtt_750_1000;
            
            
            
        }
        
    }  else if(b_min == 1000 && b_max == 1500){
        
        
        double SM_charge_asymmetry_deno_bin_mtt_1000_1500 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_1000_1500");
        double SM_charge_asymmetry_bin_mtt_1000_1500 = SM.getOptionalParameter("SM_charge_asymmetry_bin_mtt_1000_1500");
        double SM_charge_asymmetry_num_bin_mtt_1000_1500 = SM_charge_asymmetry_bin_mtt_1000_1500*SM_charge_asymmetry_deno_bin_mtt_1000_1500;
        
        
        double SM_sigma_pos_bin_mtt_1000_1500 =0.5*(SM_charge_asymmetry_num_bin_mtt_1000_1500+SM_charge_asymmetry_deno_bin_mtt_1000_1500);
        double SM_sigma_neg_bin_mtt_1000_1500 =0.5*(SM_charge_asymmetry_deno_bin_mtt_1000_1500-SM_charge_asymmetry_num_bin_mtt_1000_1500);

        double sigma_pos_bin_mtt_1000_1500_madgraph_NLO = 5702.90;
        double sigma_neg_bin_mtt_1000_1500_madgraph_NLO = 5641.72;
        
        
        if(flag_Quadratic){
            
            double sigma_pos_bin_mtt_1000_1500_NP_quad = 563.15*C_tG*C_tG + 247.32*C_Qd1*C_Qd1 + 51.83*C_Qd8*C_Qd8 + 1345.59*C_Qq11*C_Qq11 + 
            1346.53*C_Qq31*C_Qq31 + 262.79*C_Qq18*C_Qq18 + 273.66*C_Qq38*C_Qq38 + 355.65*C_Qu1*C_Qu1 + 75.51*C_Qu8*C_Qu8 + 471.12*C_td1*C_td1 +
            94.38*C_td8*C_td8 + 595.18*C_tq1*C_tq1 + 127.3*C_tq8*C_tq8 + 880.02*C_tu1*C_tu1 + 182.03*C_tu8*C_tu8;
            double sigma_neg_bin_mtt_1000_1500_NP_quad = 558.55*C_tG*C_tG + 470.26*C_Qd1*C_Qd1 + 97.31*C_Qd8*C_Qd8 + 597.83*C_Qq11*C_Qq11 +
            598.0*C_Qq31*C_Qq31 + 115.55*C_Qq18*C_Qq18 + 121.64*C_Qq38*C_Qq38 + 881.64*C_Qu1*C_Qu1 + 181.08*C_Qu8*C_Qu8 + 246.57*C_td1*C_td1 +
            47.61*C_td8*C_td8 + 1347.14*C_tq1*C_tq1 + 275.14*C_tq8*C_tq8 + 353.28*C_tu1*C_tu1 + 72.25*C_tu8*C_tu8;
            
            
            
            double sigma_pos_bin_mtt_1000_1500_NP_lin = 1591.22*C_tG + 1.74*C_Qd1 + 67.33*C_Qd8 + 25.81*C_Qq11 + 4.65*C_Qq31 + 411.2*C_Qq18 + 
            127.63*C_Qq38 + 1.19*C_Qu1 + 97.73*C_Qu8 + 10.72*C_td1 + 142.96*C_td8 + -2.4*C_tq1 + 167.79*C_tq8 + 22.4*C_tu1 + 269.18*C_tu8;
            double sigma_neg_bin_mtt_1000_1500_NP_lin = 1572.62*C_tG + -19.54*C_Qd1 + 131.94*C_Qd8 + -9.81*C_Qq11 + 1.09*C_Qq31 + 159.58*C_Qq18 +
            22.5*C_Qq38 + -37.72*C_Qu1 + 245.16*C_Qu8 + -5.86*C_td1 + 68.52*C_td8 + -49.98*C_tq1 + 378.21*C_tq8 + -9.48*C_tu1 + 92.32*C_tu8;
            
            
            double sigma_pos_bin_mtt_1000_1500_NP = sigma_pos_bin_mtt_1000_1500_NP_lin + sigma_pos_bin_mtt_1000_1500_NP_quad;
            double sigma_neg_bin_mtt_1000_1500_NP = sigma_neg_bin_mtt_1000_1500_NP_lin + sigma_neg_bin_mtt_1000_1500_NP_quad;
            
            
            
            double sigma_pos_bin_mtt_1000_1500_NP_Corrected = SM_sigma_pos_bin_mtt_1000_1500*sigma_pos_bin_mtt_1000_1500_NP/sigma_pos_bin_mtt_1000_1500_madgraph_NLO;
            double sigma_neg_bin_mtt_1000_1500_NP_Corrected = SM_sigma_neg_bin_mtt_1000_1500*sigma_neg_bin_mtt_1000_1500_NP/sigma_neg_bin_mtt_1000_1500_madgraph_NLO;
            
            double NP_charge_asymmetry_num_bin_mtt_1000_1500 = sigma_pos_bin_mtt_1000_1500_NP_Corrected - sigma_neg_bin_mtt_1000_1500_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_1000_1500 = sigma_pos_bin_mtt_1000_1500_NP_Corrected + sigma_neg_bin_mtt_1000_1500_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_1000_1500*(1+(NP_charge_asymmetry_num_bin_mtt_1000_1500-NP_charge_asymmetry_deno_bin_mtt_1000_1500)/SM_charge_asymmetry_deno_bin_mtt_1000_1500);            
        
            return SM_charge_asymmetry_bin_mtt_1000_1500 + (NP_charge_asymmetry_num_bin_mtt_1000_1500
                    - SM_charge_asymmetry_bin_mtt_1000_1500*NP_charge_asymmetry_deno_bin_mtt_1000_1500
                    )/SM_charge_asymmetry_deno_bin_mtt_1000_1500;
            
            
        }
        else{
            
            
            
            double sigma_pos_bin_mtt_1000_1500_NP = 1591.22*C_tG + 1.74*C_Qd1 + 67.33*C_Qd8 + 25.81*C_Qq11 + 4.65*C_Qq31 + 411.2*C_Qq18 + 
            127.63*C_Qq38 + 1.19*C_Qu1 + 97.73*C_Qu8 + 10.72*C_td1 + 142.96*C_td8 + -2.4*C_tq1 + 167.79*C_tq8 + 22.4*C_tu1 + 269.18*C_tu8;
            double sigma_neg_bin_mtt_1000_1500_NP = 1572.62*C_tG + -19.54*C_Qd1 + 131.94*C_Qd8 + -9.81*C_Qq11 + 1.09*C_Qq31 + 159.58*C_Qq18 +
            22.5*C_Qq38 + -37.72*C_Qu1 + 245.16*C_Qu8 + -5.86*C_td1 + 68.52*C_td8 + -49.98*C_tq1 + 378.21*C_tq8 + -9.48*C_tu1 + 92.32*C_tu8;
            
            
            double sigma_pos_bin_mtt_1000_1500_NP_Corrected = SM_sigma_pos_bin_mtt_1000_1500*sigma_pos_bin_mtt_1000_1500_NP/sigma_pos_bin_mtt_1000_1500_madgraph_NLO;
            double sigma_neg_bin_mtt_1000_1500_NP_Corrected = SM_sigma_neg_bin_mtt_1000_1500*sigma_neg_bin_mtt_1000_1500_NP/sigma_neg_bin_mtt_1000_1500_madgraph_NLO;
            
            double NP_charge_asymmetry_num_bin_mtt_1000_1500 = sigma_pos_bin_mtt_1000_1500_NP_Corrected - sigma_neg_bin_mtt_1000_1500_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_1000_1500 = sigma_pos_bin_mtt_1000_1500_NP_Corrected + sigma_neg_bin_mtt_1000_1500_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_1000_1500*(1+(NP_charge_asymmetry_num_bin_mtt_1000_1500-NP_charge_asymmetry_deno_bin_mtt_1000_1500)/SM_charge_asymmetry_deno_bin_mtt_1000_1500);            
            
            return SM_charge_asymmetry_bin_mtt_1000_1500 + (NP_charge_asymmetry_num_bin_mtt_1000_1500
                    - SM_charge_asymmetry_bin_mtt_1000_1500*NP_charge_asymmetry_deno_bin_mtt_1000_1500
                    )/SM_charge_asymmetry_deno_bin_mtt_1000_1500;
            
            
            
        }
        
    }     else if(b_min == 1500 && b_max == 3000){
        
        double SM_charge_asymmetry_bin_mtt_1500_3000 = SM.getOptionalParameter("SM_charge_asymmetry_bin_mtt_1500_3000");
        double SM_charge_asymmetry_deno_bin_mtt_1500_3000 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_1500_3000");
        double SM_charge_asymmetry_num_bin_mtt_1500_3000 = SM_charge_asymmetry_bin_mtt_1500_3000*SM_charge_asymmetry_deno_bin_mtt_1500_3000;
        
        double SM_sigma_pos_bin_mtt_1500_3000 =0.5*(SM_charge_asymmetry_num_bin_mtt_1500_3000+SM_charge_asymmetry_deno_bin_mtt_1500_3000);
        double SM_sigma_neg_bin_mtt_1500_3000 =0.5*(SM_charge_asymmetry_deno_bin_mtt_1500_3000-SM_charge_asymmetry_num_bin_mtt_1500_3000);
        
        double sigma_pos_bin_mtt_1500_3000_madgraph_NLO = 746.10;
        double sigma_neg_bin_mtt_1500_3000_madgraph_NLO = 735.27;
        
        
        if(flag_Quadratic){
        
            
            double sigma_pos_bin_mtt_1500_3000_NP_quad = 120.81*C_tG*C_tG + 202.85*C_Qd1*C_Qd1 + 42.63*C_Qd8*C_Qd8 + 1495.21*C_Qq11*C_Qq11 + 
            1494.64*C_Qq31*C_Qq31 + 314.09*C_Qq18*C_Qq18 + 321.03*C_Qq38*C_Qq38 + 449.71*C_Qu1*C_Qu1 + 93.31*C_Qu8*C_Qu8 + 471.54*C_td1*C_td1 +
            100.31*C_td8*C_td8 + 650.11*C_tq1*C_tq1 + 135.25*C_tq8*C_tq8 + 1025.45*C_tu1*C_tu1 + 220.28*C_tu8*C_tu8;
            double sigma_neg_bin_mtt_1500_3000_NP_quad = 119.3*C_tG*C_tG + 471.71*C_Qd1*C_Qd1 + 93.84*C_Qd8*C_Qd8 + 650.6*C_Qq11*C_Qq11 + 
            650.82*C_Qq31*C_Qq31 + 132.25*C_Qq18*C_Qq18 + 135.33*C_Qq38*C_Qq38 + 1025.52*C_Qu1*C_Qu1 + 203.34*C_Qu8*C_Qu8 + 202.66*C_td1*C_td1 +
            42.05*C_td8*C_td8 + 1494.96*C_tq1*C_tq1 + 296.52*C_tq8*C_tq8 + 448.95*C_tu1*C_tu1 + 92.58*C_tu8*C_tu8;
            
            
            double sigma_pos_bin_mtt_1500_3000_NP_lin = 204.01*C_tG + -2.08*C_Qd1 + 21.2*C_Qd8 + 23.76*C_Qq11 + 8.49*C_Qq31 + 172.2*C_Qq18 + 
            59.21*C_Qq38 + -3.52*C_Qu1 + 39.69*C_Qu8 + 7.84*C_td1 + 56.64*C_td8 + -5.58*C_tq1 + 61.85*C_tq8 + 16.14*C_tu1 + 115.53*C_tu8;
            double sigma_neg_bin_mtt_1500_3000_NP_lin = 200.5*C_tG + -11.33*C_Qd1 + 47.73*C_Qd8 + -0.59*C_Qq11 + -1.24*C_Qq31 + 63.93*C_Qq18 + 
            18.88*C_Qq38 + -24.09*C_Qu1 + 98.88*C_Qu8 + 0.43*C_td1 + 22.17*C_td8 + -36.32*C_tq1 + 147.64*C_tq8 + -0.96*C_tu1 + 40.92*C_tu8;
            
            double sigma_pos_bin_mtt_1500_3000_NP = sigma_pos_bin_mtt_1500_3000_NP_lin + sigma_pos_bin_mtt_1500_3000_NP_quad;
            double sigma_neg_bin_mtt_1500_3000_NP = sigma_neg_bin_mtt_1500_3000_NP_lin + sigma_neg_bin_mtt_1500_3000_NP_quad;
            
            
            double sigma_pos_bin_mtt_1500_3000_NP_Corrected = SM_sigma_pos_bin_mtt_1500_3000*sigma_pos_bin_mtt_1500_3000_NP/sigma_pos_bin_mtt_1500_3000_madgraph_NLO;
            double sigma_neg_bin_mtt_1500_3000_NP_Corrected = SM_sigma_neg_bin_mtt_1500_3000*sigma_neg_bin_mtt_1500_3000_NP/sigma_neg_bin_mtt_1500_3000_madgraph_NLO;
            
            double NP_charge_asymmetry_num_bin_mtt_1500_3000 = sigma_pos_bin_mtt_1500_3000_NP_Corrected - sigma_neg_bin_mtt_1500_3000_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_1500_3000 = sigma_pos_bin_mtt_1500_3000_NP_Corrected + sigma_neg_bin_mtt_1500_3000_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_1500_3000*(1+(NP_charge_asymmetry_num_bin_mtt_1500_3000-NP_charge_asymmetry_deno_bin_mtt_1500_3000)/SM_charge_asymmetry_deno_bin_mtt_1500_3000);            
            
            return SM_charge_asymmetry_bin_mtt_1500_3000 + (NP_charge_asymmetry_num_bin_mtt_1500_3000
                    - SM_charge_asymmetry_bin_mtt_1500_3000*NP_charge_asymmetry_deno_bin_mtt_1500_3000
                    )/SM_charge_asymmetry_deno_bin_mtt_1500_3000;
            
            
        }
        else{
            
            
            double sigma_pos_bin_mtt_1500_3000_NP = 204.01*C_tG + -2.08*C_Qd1 + 21.2*C_Qd8 + 23.76*C_Qq11 + 8.49*C_Qq31 + 172.2*C_Qq18 + 
            59.21*C_Qq38 + -3.52*C_Qu1 + 39.69*C_Qu8 + 7.84*C_td1 + 56.64*C_td8 + -5.58*C_tq1 + 61.85*C_tq8 + 16.14*C_tu1 + 115.53*C_tu8;
            double sigma_neg_bin_mtt_1500_3000_NP = 200.5*C_tG + -11.33*C_Qd1 + 47.73*C_Qd8 + -0.59*C_Qq11 + -1.24*C_Qq31 + 63.93*C_Qq18 + 
            18.88*C_Qq38 + -24.09*C_Qu1 + 98.88*C_Qu8 + 0.43*C_td1 + 22.17*C_td8 + -36.32*C_tq1 + 147.64*C_tq8 + -0.96*C_tu1 + 40.92*C_tu8;
            
            
            double sigma_pos_bin_mtt_1500_3000_NP_Corrected = SM_sigma_pos_bin_mtt_1500_3000*sigma_pos_bin_mtt_1500_3000_NP/sigma_pos_bin_mtt_1500_3000_madgraph_NLO;
            double sigma_neg_bin_mtt_1500_3000_NP_Corrected = SM_sigma_neg_bin_mtt_1500_3000*sigma_neg_bin_mtt_1500_3000_NP/sigma_neg_bin_mtt_1500_3000_madgraph_NLO;
            
            double NP_charge_asymmetry_num_bin_mtt_1500_3000 = sigma_pos_bin_mtt_1500_3000_NP_Corrected - sigma_neg_bin_mtt_1500_3000_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_1500_3000 = sigma_pos_bin_mtt_1500_3000_NP_Corrected + sigma_neg_bin_mtt_1500_3000_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_1500_3000*(1+(NP_charge_asymmetry_num_bin_mtt_1500_3000-NP_charge_asymmetry_deno_bin_mtt_1500_3000)/SM_charge_asymmetry_deno_bin_mtt_1500_3000);            

            return SM_charge_asymmetry_bin_mtt_1500_3000 + (NP_charge_asymmetry_num_bin_mtt_1500_3000
                    - SM_charge_asymmetry_bin_mtt_1500_3000*NP_charge_asymmetry_deno_bin_mtt_1500_3000
                    )/SM_charge_asymmetry_deno_bin_mtt_1500_3000;
            
        }
        
    }     else if(b_min == 1500 && b_max == 2000){
        

        
        double SM_charge_asymmetry_bin_mtt_1500_2000 = SM.getOptionalParameter("SM_charge_asymmetry_bin_mtt_1500_2000");
        double SM_charge_asymmetry_deno_bin_mtt_1500_2000 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_1500_2000");
        double SM_charge_asymmetry_num_bin_mtt_1500_2000 = SM_charge_asymmetry_bin_mtt_1500_2000*SM_charge_asymmetry_deno_bin_mtt_1500_2000;
        
        double SM_sigma_pos_bin_mtt_1500_2000 =0.5*(SM_charge_asymmetry_num_bin_mtt_1500_2000+SM_charge_asymmetry_deno_bin_mtt_1500_2000);
        double SM_sigma_neg_bin_mtt_1500_2000 =0.5*(SM_charge_asymmetry_deno_bin_mtt_1500_2000-SM_charge_asymmetry_num_bin_mtt_1500_2000);
        
        
        //double Delta_Madgraph = -1.400000000000091;
        //double Total_Madgraph = 1519.0;
        //double pos_Madgraph = 758.8;
        //double neg_Madgraph = 760.2;
        
        double sigma_pos_bin_mtt_1500_2000_madgraph = 623.06;
        double sigma_neg_bin_mtt_1500_2000_madgraph = 613.37;
        
        
        
        
        if(flag_Quadratic){
            
        
            //double Delta = Delta_Madgraph  +0.6*C_tG  +127.3*C_td1*C_td1  +26.33*C_td8  +28.28*C_td8*C_td8  +322.4*C_tu1*C_tu1  +65.61*C_tu8  +71.49*C_tu8*C_tu8  +449.8*C_Qq11*C_Qq11  +449.6*C_Qq31*C_Qq31  +96.74*C_Qq18  +99.97*C_Qq18*C_Qq18  +46.49*C_Qq38  +99.92*C_Qq38*C_Qq38  +2.194*C_tG*C_td8  +6.783*C_tG*C_tu8  +9.33*C_tG*C_Qq18  +3.806*C_tG*C_Qq38  +0.67*C_Qd1*C_td1  +0.163*C_Qu8*C_tu8  +389.5*C_Qq11*C_Qq31  +86.86*C_Qq18*C_Qq38 ; 
            //double Total = Total_Madgraph  +365.6*C_tG  +193.64*C_tG*C_tG  +81.69*C_Qd8  +86.77*C_Qd8*C_Qd8  +390.8*C_Qd1*C_Qd1  +729.9*C_Qu1*C_Qu1  +161.18*C_Qu8  +162.26*C_Qu8*C_Qu8  +390.7*C_td1*C_td1  +81.07*C_td8  +86.88*C_td8*C_td8  +1118.8*C_tq1*C_tq1  +237.84*C_tq8  +248.78*C_tq8*C_tq8  +729.6*C_tu1*C_tu1  +150.79*C_tu8  +162.11*C_tu8*C_tu8  +1119.2*C_Qq11*C_Qq11  +1119.2*C_Qq31*C_Qq31  +230.46*C_Qq18  +248.63*C_Qq18*C_Qq18  +72.81*C_Qq38  +248.68*C_Qq38*C_Qq38  +12.235*C_tG*C_Qd8  +21.962*C_tG*C_Qu8  +11.538*C_tG*C_td8  +34.19*C_tG*C_tq8  +22.517*C_tG*C_tu8  +34.89*C_tG*C_Qq18  +10.308*C_tG*C_Qq38  +5.379*C_Qd8*C_td8  +23.45*C_Qd1*C_td1  +45.23*C_Qu1*C_tu1  +9.957*C_Qu8*C_tu8  +68.33*C_tq1*C_Qq11  +21.29*C_tq1*C_Qq31  +15.024*C_tq8*C_Qq18  +4.593*C_tq8*C_Qq38  +681.1*C_Qq11*C_Qq31  +151.74*C_Qq18*C_Qq38 ; 
            
            
            double sigma_pos_bin_mtt_1500_2000_NP_quad = 92.57*C_tG*C_tG + 118.16*C_Qd1*C_Qd1 + 24.32*C_Qd8*C_Qd8 + 734.31*C_Qq11*C_Qq11 + 
            732.63*C_Qq31*C_Qq31 + 148.8*C_Qq18*C_Qq18 + 152.52*C_Qq38*C_Qq38 + 201.03*C_Qu1*C_Qu1 + 41.1*C_Qu8*C_Qu8 
            + 245.73*C_td1*C_td1 + 50.24*C_td8*C_td8 + 317.54*C_tq1*C_tq1 + 65.51*C_tq8*C_tq8 + 489.18*C_tu1*C_tu1 + 104.44*C_tu8*C_tu8;
            double sigma_neg_bin_mtt_1500_2000_NP_quad = 91.95*C_tG*C_tG + 245.59*C_Qd1*C_Qd1 + 48.46*C_Qd8*C_Qd8 + 318.95*C_Qq11*C_Qq11 +
            320.29*C_Qq31*C_Qq31 + 63.03*C_Qq18*C_Qq18 + 65.21*C_Qq38*C_Qq38 + 490.42*C_Qu1*C_Qu1 + 96.83*C_Qu8*C_Qu8 +
            118.57*C_td1*C_td1 + 24.0*C_td8*C_td8 + 733.93*C_tq1*C_tq1 + 145.07*C_tq8*C_tq8 + 202.0*C_tu1*C_tu1 + 40.84*C_tu8*C_tu8;
            
            
            double sigma_pos_bin_mtt_1500_2000_NP_lin = 170.81*C_tG + -1.77*C_Qd1 + 15.41*C_Qd8 + 12.18*C_Qq11 + 2.58*C_Qq31 + 114.92*C_Qq18 +
            37.21*C_Qq38 + -1.56*C_Qu1 + 25.46*C_Qu8 + 5.03*C_td1 + 39.19*C_td8 + -3.63*C_tq1 + 42.28*C_tq8 + 10.0*C_tu1 + 75.73*C_tu8;
            double sigma_neg_bin_mtt_1500_2000_NP_lin = 167.48*C_tG + -6.39*C_Qd1 + 32.97*C_Qd8 + -0.51*C_Qq11 + 0.27*C_Qq31 + 42.69*C_Qq18 + 
            9.01*C_Qq38 + -15.76*C_Qu1 + 65.88*C_Qu8 + 0.08*C_td1 + 16.22*C_td8 + -23.19*C_tq1 + 100.18*C_tq8 + -2.51*C_tu1 + 25.44*C_tu8;

            
            double sigma_pos_bin_mtt_1500_2000_NP = sigma_pos_bin_mtt_1500_2000_NP_lin + sigma_pos_bin_mtt_1500_2000_NP_quad;
            double sigma_neg_bin_mtt_1500_2000_NP = sigma_neg_bin_mtt_1500_2000_NP_lin + sigma_neg_bin_mtt_1500_2000_NP_quad;
            
            
            double sigma_pos_bin_mtt_1500_2000_NP_Corrected = SM_sigma_pos_bin_mtt_1500_2000*sigma_pos_bin_mtt_1500_2000_NP/sigma_pos_bin_mtt_1500_2000_madgraph;
            double sigma_neg_bin_mtt_1500_2000_NP_Corrected = SM_sigma_neg_bin_mtt_1500_2000*sigma_neg_bin_mtt_1500_2000_NP/sigma_neg_bin_mtt_1500_2000_madgraph;
            
            double NP_charge_asymmetry_num_bin_mtt_1500_2000 = sigma_pos_bin_mtt_1500_2000_NP_Corrected - sigma_neg_bin_mtt_1500_2000_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_1500_2000 = sigma_pos_bin_mtt_1500_2000_NP_Corrected + sigma_neg_bin_mtt_1500_2000_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_1500_2000*(1+(NP_charge_asymmetry_num_bin_mtt_1500_2000-NP_charge_asymmetry_deno_bin_mtt_1500_2000)/SM_charge_asymmetry_deno_bin_mtt_1500_2000);            

            return SM_charge_asymmetry_bin_mtt_1500_2000 + (NP_charge_asymmetry_num_bin_mtt_1500_2000
                    - SM_charge_asymmetry_bin_mtt_1500_2000*NP_charge_asymmetry_deno_bin_mtt_1500_2000
                    )/SM_charge_asymmetry_deno_bin_mtt_1500_2000;
            
        }
        else{
            
            
            double sigma_pos_bin_mtt_1500_2000_NP = 170.81*C_tG + -1.77*C_Qd1 + 15.41*C_Qd8 + 12.18*C_Qq11 + 2.58*C_Qq31 + 114.92*C_Qq18 +
            37.21*C_Qq38 + -1.56*C_Qu1 + 25.46*C_Qu8 + 5.03*C_td1 + 39.19*C_td8 + -3.63*C_tq1 + 42.28*C_tq8 + 10.0*C_tu1 + 75.73*C_tu8;
            double sigma_neg_bin_mtt_1500_2000_NP = 167.48*C_tG + -6.39*C_Qd1 + 32.97*C_Qd8 + -0.51*C_Qq11 + 0.27*C_Qq31 + 42.69*C_Qq18 + 
            9.01*C_Qq38 + -15.76*C_Qu1 + 65.88*C_Qu8 + 0.08*C_td1 + 16.22*C_td8 + -23.19*C_tq1 + 100.18*C_tq8 + -2.51*C_tu1 + 25.44*C_tu8;

            
            double sigma_pos_bin_mtt_1500_2000_NP_Corrected = SM_sigma_pos_bin_mtt_1500_2000*sigma_pos_bin_mtt_1500_2000_NP/sigma_pos_bin_mtt_1500_2000_madgraph;
            double sigma_neg_bin_mtt_1500_2000_NP_Corrected = SM_sigma_neg_bin_mtt_1500_2000*sigma_neg_bin_mtt_1500_2000_NP/sigma_neg_bin_mtt_1500_2000_madgraph;
            
            double NP_charge_asymmetry_num_bin_mtt_1500_2000 = sigma_pos_bin_mtt_1500_2000_NP_Corrected - sigma_neg_bin_mtt_1500_2000_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_1500_2000 = sigma_pos_bin_mtt_1500_2000_NP_Corrected + sigma_neg_bin_mtt_1500_2000_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_1500_2000*(1+(NP_charge_asymmetry_num_bin_mtt_1500_2000-NP_charge_asymmetry_deno_bin_mtt_1500_2000)/SM_charge_asymmetry_deno_bin_mtt_1500_2000);            
            
            return SM_charge_asymmetry_bin_mtt_1500_2000 + (NP_charge_asymmetry_num_bin_mtt_1500_2000
                    - SM_charge_asymmetry_bin_mtt_1500_2000*NP_charge_asymmetry_deno_bin_mtt_1500_2000
                    )/SM_charge_asymmetry_deno_bin_mtt_1500_2000;
            
        }
        
    }     else if(b_min == 2000 && b_max == 2500){
        

        double SM_charge_asymmetry_bin_mtt_2000_2500 = SM.getOptionalParameter("SM_charge_asymmetry_bin_mtt_2000_2500");
        double SM_charge_asymmetry_deno_bin_mtt_2000_2500 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_2000_2500");
        double SM_charge_asymmetry_num_bin_mtt_2000_2500 = SM_charge_asymmetry_bin_mtt_2000_2500*SM_charge_asymmetry_deno_bin_mtt_2000_2500;
        
        double SM_sigma_pos_bin_mtt_2000_2500 =0.5*(SM_charge_asymmetry_num_bin_mtt_2000_2500+SM_charge_asymmetry_deno_bin_mtt_2000_2500);
        double SM_sigma_neg_bin_mtt_2000_2500 =0.5*(SM_charge_asymmetry_deno_bin_mtt_2000_2500-SM_charge_asymmetry_num_bin_mtt_2000_2500);
        
        
        //double Delta_Madgraph = 9.899999999999977;
        //double Total_Madgraph = 280.5;
        //double pos_Madgraph = 145.2;
        //double neg_Madgraph = 135.3;
        
        double sigma_pos_bin_mtt_2000_2500_madgraph = 91.28;
        double sigma_neg_bin_mtt_2000_2500_madgraph = 89.45;
        
        
        if(flag_Quadratic){

            
            double sigma_pos_bin_mtt_2000_2500_NP_quad = 19.16*C_tG*C_tG + 52.13*C_Qd1*C_Qd1 + 11.31*C_Qd8*C_Qd8 + 380.79*C_Qq11*C_Qq11 + 
            380.23*C_Qq31*C_Qq31 + 80.68*C_Qq18*C_Qq18 + 83.66*C_Qq38*C_Qq38 + 111.34*C_Qu1*C_Qu1 + 24.39*C_Qu8*C_Qu8 +
            121.22*C_td1*C_td1 + 26.15*C_td8*C_td8 + 162.15*C_tq1*C_tq1 + 35.82*C_tq8*C_tq8 + 259.37*C_tu1*C_tu1 + 55.42*C_tu8*C_tu8;
            double sigma_neg_bin_mtt_2000_2500_NP_quad = 18.96*C_tG*C_tG + 122.29*C_Qd1*C_Qd1 + 25.74*C_Qd8*C_Qd8 + 161.29*C_Qq11*C_Qq11 + 
            160.46*C_Qq31*C_Qq31 + 33.11*C_Qq18*C_Qq18 + 33.92*C_Qq38*C_Qq38 + 259.13*C_Qu1*C_Qu1 + 54.48*C_Qu8*C_Qu8 +
            52.33*C_td1*C_td1 + 11.05*C_td8*C_td8 + 379.61*C_tq1*C_tq1 + 79.62*C_tq8*C_tq8 + 110.6*C_tu1*C_tu1 + 22.87*C_tu8*C_tu8;
            
            
            double sigma_pos_bin_mtt_2000_2500_NP_lin = 24.76*C_tG + -0.52*C_Qd1 + 3.67*C_Qd8 + 3.98*C_Qq11 + -0.13*C_Qq31 + 36.66*C_Qq18 + 
            12.61*C_Qq38 + -0.59*C_Qu1 + 8.95*C_Qu8 + 1.25*C_td1 + 11.32*C_td8 + -0.29*C_tq1 + 13.01*C_tq8 + 1.89*C_tu1 + 24.59*C_tu8;
            double sigma_neg_bin_mtt_2000_2500_NP_lin = 24.28*C_tG + -2.52*C_Qd1 + 10.16*C_Qd8 + 0.47*C_Qq11 + 0.78*C_Qq31 + 13.11*C_Qq18 + 
            4.73*C_Qq38 + -4.64*C_Qu1 + 21.13*C_Qu8 + 0.9*C_td1 + 4.33*C_td8 + -7.25*C_tq1 + 31.39*C_tq8 + 1.68*C_tu1 + 8.81*C_tu8;
            
            double sigma_pos_bin_mtt_2000_2500_NP = sigma_pos_bin_mtt_2000_2500_NP_lin + sigma_pos_bin_mtt_2000_2500_NP_quad;
            double sigma_neg_bin_mtt_2000_2500_NP = sigma_neg_bin_mtt_2000_2500_NP_lin + sigma_neg_bin_mtt_2000_2500_NP_quad;
            
            
            double sigma_pos_bin_mtt_2000_2500_NP_Corrected = SM_sigma_pos_bin_mtt_2000_2500*sigma_pos_bin_mtt_2000_2500_NP/sigma_pos_bin_mtt_2000_2500_madgraph;
            double sigma_neg_bin_mtt_2000_2500_NP_Corrected = SM_sigma_neg_bin_mtt_2000_2500*sigma_neg_bin_mtt_2000_2500_NP/sigma_neg_bin_mtt_2000_2500_madgraph;
            
            double NP_charge_asymmetry_num_bin_mtt_2000_2500 = sigma_pos_bin_mtt_2000_2500_NP_Corrected - sigma_neg_bin_mtt_2000_2500_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_2000_2500 = sigma_pos_bin_mtt_2000_2500_NP_Corrected + sigma_neg_bin_mtt_2000_2500_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_2000_2500*(1+(NP_charge_asymmetry_num_bin_mtt_2000_2500-NP_charge_asymmetry_deno_bin_mtt_2000_2500)/SM_charge_asymmetry_deno_bin_mtt_2000_2500);            
            
            return SM_charge_asymmetry_bin_mtt_2000_2500 + (NP_charge_asymmetry_num_bin_mtt_2000_2500
                    - SM_charge_asymmetry_bin_mtt_2000_2500*NP_charge_asymmetry_deno_bin_mtt_2000_2500
                    )/SM_charge_asymmetry_deno_bin_mtt_2000_2500;
            
        }
        else{
            

            double sigma_pos_bin_mtt_2000_2500_NP = 24.76*C_tG + -0.52*C_Qd1 + 3.67*C_Qd8 + 3.98*C_Qq11 + -0.13*C_Qq31 + 36.66*C_Qq18 +
            12.61*C_Qq38 + -0.59*C_Qu1 + 8.95*C_Qu8 + 1.25*C_td1 + 11.32*C_td8 + -0.29*C_tq1 + 13.01*C_tq8 + 1.89*C_tu1 + 24.59*C_tu8;
            double sigma_neg_bin_mtt_2000_2500_NP = 24.28*C_tG + -2.52*C_Qd1 + 10.16*C_Qd8 + 0.47*C_Qq11 + 0.78*C_Qq31 + 13.11*C_Qq18 + 
            4.73*C_Qq38 + -4.64*C_Qu1 + 21.13*C_Qu8 + 0.9*C_td1 + 4.33*C_td8 + -7.25*C_tq1 + 31.39*C_tq8 + 1.68*C_tu1 + 8.81*C_tu8;
            
            
            
            double sigma_pos_bin_mtt_2000_2500_NP_Corrected = SM_sigma_pos_bin_mtt_2000_2500*sigma_pos_bin_mtt_2000_2500_NP/sigma_pos_bin_mtt_2000_2500_madgraph;
            double sigma_neg_bin_mtt_2000_2500_NP_Corrected = SM_sigma_neg_bin_mtt_2000_2500*sigma_neg_bin_mtt_2000_2500_NP/sigma_neg_bin_mtt_2000_2500_madgraph;
            
            double NP_charge_asymmetry_num_bin_mtt_2000_2500 = sigma_pos_bin_mtt_2000_2500_NP_Corrected - sigma_neg_bin_mtt_2000_2500_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_2000_2500 = sigma_pos_bin_mtt_2000_2500_NP_Corrected + sigma_neg_bin_mtt_2000_2500_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_2000_2500*(1+(NP_charge_asymmetry_num_bin_mtt_2000_2500-NP_charge_asymmetry_deno_bin_mtt_2000_2500)/SM_charge_asymmetry_deno_bin_mtt_2000_2500);            
            
            return SM_charge_asymmetry_bin_mtt_2000_2500 + (NP_charge_asymmetry_num_bin_mtt_2000_2500
                    - SM_charge_asymmetry_bin_mtt_2000_2500*NP_charge_asymmetry_deno_bin_mtt_2000_2500
                    )/SM_charge_asymmetry_deno_bin_mtt_2000_2500;
            
            
            
        }
        
    }     else if(b_min == 2500 && b_max == 3000){
        
        
        double SM_charge_asymmetry_bin_mtt_2500_3000 = SM.getOptionalParameter("SM_charge_asymmetry_bin_mtt_2500_3000");
        double SM_charge_asymmetry_deno_bin_mtt_2500_3000 = SM.getOptionalParameter("SM_charge_asymmetry_deno_bin_mtt_2500_3000");
        double SM_charge_asymmetry_num_bin_mtt_2500_3000 = SM_charge_asymmetry_bin_mtt_2500_3000*SM_charge_asymmetry_deno_bin_mtt_2500_3000;
        
        double SM_sigma_pos_bin_mtt_2500_3000 =0.5*(SM_charge_asymmetry_num_bin_mtt_2500_3000+SM_charge_asymmetry_deno_bin_mtt_2500_3000);
        double SM_sigma_neg_bin_mtt_2500_3000 =0.5*(SM_charge_asymmetry_deno_bin_mtt_2500_3000-SM_charge_asymmetry_num_bin_mtt_2500_3000);
        
        
        
        
        //double Delta_Madgraph = -2.8999999999999986;
        //double Total_Madgraph = 61.12;
        //double pos_Madgraph = 29.11;
        //double neg_Madgraph = 32.01;
        
        
        double sigma_pos_bin_mtt_2500_3000_madgraph = 17.54;
        double sigma_neg_bin_mtt_2500_3000_madgraph = 17.24;
        
        
        if(flag_Quadratic){
        
            //double Delta = Delta_Madgraph  +51.68*C_td1*C_td1  +3.262*C_td8  +11.478*C_td8*C_td8  +120.07*C_tu1*C_tu1  +9.458*C_tu8  +26.67*C_tu8*C_tu8  +171.8*C_Qq11*C_Qq11  +172.2*C_Qq31*C_Qq31  +13.305*C_Qq18  +38.22*C_Qq18*C_Qq18  +5.559*C_Qq38  +38.21*C_Qq38*C_Qq38  +0.564*C_tG*C_td8  +0.901*C_tG*C_tu8  +1.423*C_tG*C_Qq18  +0.386*C_tG*C_Qq38  +0.006*C_Qd8*C_td8  +0.624*C_Qu1*C_tu1  +0.232*C_Qu8*C_tu8  +0.243*C_tq1*C_Qq11  +0.001*C_tq8*C_Qq38  +135.8*C_Qq11*C_Qq31  +30.27*C_Qq18*C_Qq38 ; 
            //double Total = Total_Madgraph  +14.606*C_tG  +15.074*C_tG*C_tG  +11.4*C_Qd8  +29.416*C_Qd8*C_Qd8  +132.52*C_Qd1*C_Qd1  +302.15*C_Qu1*C_Qu1  +21.356*C_Qu8  +67.13*C_Qu8*C_Qu8  +132.38*C_td1*C_td1  +9.424*C_td8  +29.482*C_td8*C_td8  +434.2*C_tq1*C_tq1  +38.76*C_tq8  +96.46*C_tq8*C_tq8  +302.13*C_tu1*C_tu1  +27.042*C_tu8  +67.15*C_tu8*C_tu8  +434.2*C_Qq11*C_Qq11  +434.4*C_Qq31*C_Qq31  +32.115*C_Qq18  +96.56*C_Qq18*C_Qq18  +12.195*C_Qq38  +96.51*C_Qq38*C_Qq38  +1.869*C_tG*C_Qd8  +3.491*C_tG*C_Qu8  +1.346*C_tG*C_td8  +5.49*C_tG*C_tq8  +3.599*C_tG*C_tu8  +5.263*C_tG*C_Qq18  +2.061*C_tG*C_Qq38  +0.653*C_Qd8*C_td8  +3.018*C_Qd1*C_td1  +7.298*C_Qu1*C_tu1  +1.591*C_Qu8*C_tu8  +10.329*C_tq1*C_Qq11  +3.613*C_tq1*C_Qq31  +2.332*C_tq8*C_Qq18  +0.958*C_tq8*C_Qq38  +339.2*C_Qq11*C_Qq31  +75.31*C_Qq18*C_Qq38 ; 


            
            double sigma_pos_bin_mtt_2500_3000_NP_quad = 4.98*C_tG*C_tG + 21.28*C_Qd1*C_Qd1 + 4.74*C_Qd8*C_Qd8 + 192.55*C_Qq11*C_Qq11 + 
            194.11*C_Qq31*C_Qq31 + 41.97*C_Qq18*C_Qq18 + 43.67*C_Qq38*C_Qq38 + 61.34*C_Qu1*C_Qu1 + 13.52*C_Qu8*C_Qu8 + 
            58.14*C_td1*C_td1 + 12.9*C_td8*C_td8 + 82.9*C_tq1*C_tq1 + 18.17*C_tq8*C_tq8 + 134.43*C_tu1*C_tu1 + 29.76*C_tu8*C_tu8;
            double sigma_neg_bin_mtt_2500_3000_NP_quad = 4.85*C_tG*C_tG + 58.14*C_Qd1*C_Qd1 + 12.44*C_Qd8*C_Qd8 + 83.33*C_Qq11*C_Qq11 +
            82.69*C_Qq31*C_Qq31 + 17.09*C_Qq18*C_Qq18 + 17.76*C_Qq38*C_Qq38 + 135.1*C_Qu1*C_Qu1 + 28.74*C_Qu8*C_Qu8 + 
            21.41*C_td1*C_td1 + 4.69*C_td8*C_td8 + 194.18*C_tq1*C_tq1 + 40.96*C_tq8*C_tq8 + 62.11*C_tu1*C_tu1 + 12.69*C_tu8*C_tu8;
            
            
            double sigma_pos_bin_mtt_2500_3000_NP_lin = 4.73*C_tG + -0.44*C_Qd1 + 1.06*C_Qd8 + 0.43*C_Qq11 + 0.7*C_Qq31 + 12.57*C_Qq18 +
            4.94*C_Qq38 + -0.45*C_Qu1 + 3.36*C_Qu8 + 0.91*C_td1 + 3.77*C_td8 + 0.43*C_tq1 + 4.39*C_tq8 + 1.14*C_tu1 + 8.57*C_tu8;
            double sigma_neg_bin_mtt_2500_3000_NP_lin = 4.61*C_tG + -0.79*C_Qd1 + 3.23*C_Qd8 + 0.53*C_Qq11 + -0.51*C_Qq31 + 4.9*C_Qq18 + 
            2.42*C_Qq38 + -1.35*C_Qu1 + 7.23*C_Qu8 + -0.29*C_td1 + 1.18*C_td8 + -2.82*C_tq1 + 10.79*C_tq8 + 0.59*C_tu1 + 3.24*C_tu8;
            
            
            double sigma_pos_bin_mtt_2500_3000_NP = sigma_pos_bin_mtt_2500_3000_NP_lin + sigma_pos_bin_mtt_2500_3000_NP_quad;
            double sigma_neg_bin_mtt_2500_3000_NP = sigma_neg_bin_mtt_2500_3000_NP_lin + sigma_neg_bin_mtt_2500_3000_NP_quad;
            
            
            double sigma_pos_bin_mtt_2500_3000_NP_Corrected = SM_sigma_pos_bin_mtt_2500_3000*sigma_pos_bin_mtt_2500_3000_NP/sigma_pos_bin_mtt_2500_3000_madgraph;
            double sigma_neg_bin_mtt_2500_3000_NP_Corrected = SM_sigma_neg_bin_mtt_2500_3000*sigma_neg_bin_mtt_2500_3000_NP/sigma_neg_bin_mtt_2500_3000_madgraph;
            
            double NP_charge_asymmetry_num_bin_mtt_2500_3000 = sigma_pos_bin_mtt_2500_3000_NP_Corrected - sigma_neg_bin_mtt_2500_3000_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_2500_3000 = sigma_pos_bin_mtt_2500_3000_NP_Corrected + sigma_neg_bin_mtt_2500_3000_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_2500_3000*(1+(NP_charge_asymmetry_num_bin_mtt_2500_3000-NP_charge_asymmetry_deno_bin_mtt_2500_3000)/SM_charge_asymmetry_deno_bin_mtt_2500_3000);            
            
            return SM_charge_asymmetry_bin_mtt_2500_3000 + (NP_charge_asymmetry_num_bin_mtt_2500_3000
                    - SM_charge_asymmetry_bin_mtt_2500_3000*NP_charge_asymmetry_deno_bin_mtt_2500_3000
                    )/SM_charge_asymmetry_deno_bin_mtt_2500_3000;
            
            
        }
        else{
            
            
            double sigma_pos_bin_mtt_2500_3000_NP = 4.73*C_tG + -0.44*C_Qd1 + 1.06*C_Qd8 + 0.43*C_Qq11 + 0.7*C_Qq31 + 12.57*C_Qq18 +
            4.94*C_Qq38 + -0.45*C_Qu1 + 3.36*C_Qu8 + 0.91*C_td1 + 3.77*C_td8 + 0.43*C_tq1 + 4.39*C_tq8 + 1.14*C_tu1 + 8.57*C_tu8;
            double sigma_neg_bin_mtt_2500_3000_NP = 4.61*C_tG + -0.79*C_Qd1 + 3.23*C_Qd8 + 0.53*C_Qq11 + -0.51*C_Qq31 + 4.9*C_Qq18 + 
            2.42*C_Qq38 + -1.35*C_Qu1 + 7.23*C_Qu8 + -0.29*C_td1 + 1.18*C_td8 + -2.82*C_tq1 + 10.79*C_tq8 + 0.59*C_tu1 + 3.24*C_tu8;
            
            
            
            double sigma_pos_bin_mtt_2500_3000_NP_Corrected = SM_sigma_pos_bin_mtt_2500_3000*sigma_pos_bin_mtt_2500_3000_NP/sigma_pos_bin_mtt_2500_3000_madgraph;
            double sigma_neg_bin_mtt_2500_3000_NP_Corrected = SM_sigma_neg_bin_mtt_2500_3000*sigma_neg_bin_mtt_2500_3000_NP/sigma_neg_bin_mtt_2500_3000_madgraph;
            
            double NP_charge_asymmetry_num_bin_mtt_2500_3000 = sigma_pos_bin_mtt_2500_3000_NP_Corrected - sigma_neg_bin_mtt_2500_3000_NP_Corrected;
            double NP_charge_asymmetry_deno_bin_mtt_2500_3000 = sigma_pos_bin_mtt_2500_3000_NP_Corrected + sigma_neg_bin_mtt_2500_3000_NP_Corrected;
            
            //return  SM_charge_asymmetry_bin_mtt_2500_3000*(1+(NP_charge_asymmetry_num_bin_mtt_2500_3000-NP_charge_asymmetry_deno_bin_mtt_2500_3000)/SM_charge_asymmetry_deno_bin_mtt_2500_3000);            
            
            
            return SM_charge_asymmetry_bin_mtt_2500_3000 + (NP_charge_asymmetry_num_bin_mtt_2500_3000
                    - SM_charge_asymmetry_bin_mtt_2500_3000*NP_charge_asymmetry_deno_bin_mtt_2500_3000
                    )/SM_charge_asymmetry_deno_bin_mtt_2500_3000;
            
            
        }
    }

    
    else {
        throw std::runtime_error("\nERROR: Please specify a correct binning range for charge_asymmetry_tt_diff_mtt_NLO.\n");
    }

}

/*
Charge_Asymmetry_bin_tt_0_500::Charge_Asymmetry_bin_tt_0_500(const StandardModel& SM_i)
 : ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i)) 
{}

double    Charge_Asymmetry_bin_tt_0_500::computeThValue() 
{ 
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double SM_Charge_Asymmetry_bin_tt_0_500 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_0_500();
    bool   flag_Quadratic = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic(); 
    double Delta_Madgraph = -165.81770000001416;
    double Total_Madgraph = 297014.3397;


if(flag_Quadratic){ 

    double Delta = Delta_Madgraph  +277.16*C_tG  +39.956*C_tG*C_tG  +39.489*C_td1*C_td1  +119.887*C_td8  +8.478*C_td8*C_td8  +90.721*C_tu1*C_tu1  +282.549*C_tu8  +19.913*C_tu8*C_tu8  +127.866*C_Qq11*C_Qq11  +128.447*C_Qq31*C_Qq31  +405.634*C_Qq18  +28.703*C_Qq18*C_Qq18  +154.036*C_Qq38  +28.948*C_Qq38*C_Qq38  +0.526*C_Qd8*C_Qu8  +0.122*C_Qd8*C_td8  +0.249*C_Qd8*C_tq8  +0.525*C_Qd8*C_tu8  +1.443*C_Qd1*C_Qu1  +1.655*C_Qd1*C_td1  +3.896*C_Qd1*C_tq1  +4.615*C_Qd1*C_Qq11  +0.035*C_Qd1*C_Qq31  +3.634*C_Qu1*C_td1  +2.953*C_Qu1*C_tq1  +1.274*C_Qu1*C_tu1  +6.568*C_Qu1*C_Qq31  +0.326*C_Qu8*C_td8  +0.626*C_Qu8*C_tq8  +0.773*C_Qu8*C_tu8  +0.089*C_Qu8*C_Qq18  +0.159*C_Qu8*C_Qq38  +1.409*C_td1*C_Qq11  +1.573*C_td1*C_Qq31  +0.572*C_td8*C_tu8  +0.234*C_td8*C_Qq38  +4.073*C_tq1*C_Qq11  +7.01*C_tq1*C_Qq31  +1.222*C_tq8*C_tu8  +3.926*C_tu1*C_Qq11  +1.623*C_tu1*C_Qq31  +0.264*C_tu8*C_Qq18  +0.928*C_tu8*C_Qq38  +105.538*C_Qq11*C_Qq31  +21.858*C_Qq18*C_Qq38 ; 

    double Total = Total_Madgraph  +82478.882*C_tG  +10564.131*C_tG*C_tG  +1101.588*C_Qd8  +46.089*C_Qd8*C_Qd8  +204.659*C_Qd1*C_Qd1  +292.967*C_Qu1*C_Qu1  +1547.089*C_Qu8  +65.096*C_Qu8*C_Qu8  +206.313*C_td1*C_td1  +1097.863*C_td8  +45.901*C_td8*C_td8  +492.8*C_tq1*C_tq1  +2615.608*C_tq8  +109.994*C_tq8*C_tq8  +292.982*C_tu1*C_tu1  +1548.282*C_tu8  +65.134*C_tu8*C_tu8  +493.673*C_Qq11*C_Qq11  +494.118*C_Qq31*C_Qq31  +2613.206*C_Qq18  +109.656*C_Qq18*C_Qq18  +484.73*C_Qq38  +110.159*C_Qq38*C_Qq38  +101.448*C_tG*C_Qd8  +132.612*C_tG*C_Qu8  +113.258*C_tG*C_td8  +296.285*C_tG*C_tq8  +156.073*C_tG*C_tu8  +275.004*C_tG*C_Qq18  +41.635*C_tG*C_Qq38  +53.622*C_Qd8*C_td8  +0.252*C_Qd8*C_tu8  +0.425*C_Qd8*C_Qq18  +0.823*C_Qd1*C_Qu1  +241.252*C_Qd1*C_td1  +2.989*C_Qd1*C_tq1  +3.231*C_Qd1*C_tu1  +0.765*C_Qd1*C_Qq11  +0.16*C_Qu1*C_td1  +2.883*C_Qu1*C_tq1  +339.918*C_Qu1*C_tu1  +0.251*C_Qu8*C_tq8  +74.004*C_Qu8*C_tu8  +0.457*C_Qu8*C_Qq38  +4.932*C_td1*C_tq1  +2.635*C_td1*C_Qq11  +0.371*C_td8*C_Qq18  +0.661*C_tq1*C_tu1  +572.182*C_tq1*C_Qq11  +107.968*C_tq1*C_Qq31  +127.46*C_tq8*C_Qq18  +23.662*C_tq8*C_Qq38  +0.936*C_tu1*C_Qq11  +1.495*C_tu1*C_Qq31  +0.216*C_tu8*C_Qq18  +184.842*C_Qq11*C_Qq31  +41.032*C_Qq18*C_Qq38 ; 

     return SM_Charge_Asymmetry_bin_tt_0_500 + Delta/Total;
 } 
else{ 

    double Delta = Delta_Madgraph+277.16*C_tG+119.887*C_td8+282.549*C_tu8+405.634*C_Qq18+154.036*C_Qq38 ; 

    double Total = Total_Madgraph+154.036*C_tG+154.036*C_Qd8+154.036*C_Qu8+154.036*C_td8+154.036*C_tq8+154.036*C_tu8+154.036*C_Qq18+154.036*C_Qq38 ; 

    return SM_Charge_Asymmetry_bin_tt_0_500 + Delta/Total;
 } 

} 

Charge_Asymmetry_bin_tt_500_750::Charge_Asymmetry_bin_tt_500_750(const StandardModel& SM_i)
 : ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i)) 
{}

double    Charge_Asymmetry_bin_tt_500_750::computeThValue() 
{ 
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double SM_Charge_Asymmetry_bin_tt_500_750 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_500_750();
    bool   flag_Quadratic = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic(); 
    double Delta_Madgraph = 230.0;
    double Total_Madgraph = 176710.0;


if(flag_Quadratic){ 

    double Delta = Delta_Madgraph  +60.0*C_tG  +1.0*C_tG*C_tG  +120.2*C_td1*C_td1  +190.0*C_td8  +26.73*C_td8*C_td8  +285.0*C_tu1*C_tu1  +462.2*C_tu8  +63.45*C_tu8*C_tu8  +405.3*C_Qq11*C_Qq11  +405.0*C_Qq31*C_Qq31  +659.8*C_Qq18  +90.15*C_Qq18*C_Qq18  +264.5*C_Qq38  +90.02*C_Qq38*C_Qq38  +18.2*C_tG*C_td8  +44.18*C_tG*C_tu8  +63.6*C_tG*C_Qq18  +25.41*C_tG*C_Qq38  +0.3*C_tq1*C_Qq11  +0.16*C_tq8*C_Qq18  +330.45*C_Qq11*C_Qq31  +73.02*C_Qq18*C_Qq38 ; 

    double Total = Total_Madgraph  +43260.0*C_tG  +6835.0*C_tG*C_tG  +943.0*C_Qd8  +101.99*C_Qd8*C_Qd8  +458.8*C_Qd1*C_Qd1  +693.2*C_Qu1*C_Qu1  +1415.3*C_Qu8  +154.02*C_Qu8*C_Qu8  +459.0*C_td1*C_td1  +935.2*C_td8  +102.01*C_td8*C_td8  +1143.0*C_tq1*C_tq1  +2324.4*C_tq8  +254.06*C_tq8*C_tq8  +693.2*C_tu1*C_tu1  +1418.4*C_tu8  +154.15*C_tu8*C_tu8  +1143.1*C_Qq11*C_Qq11  +1143.2*C_Qq31*C_Qq31  +2354.2*C_Qq18  +254.05*C_Qq18*C_Qq18  +468.5*C_Qq38  +253.98*C_Qq38*C_Qq38  +119.27*C_tG*C_Qd8  +179.2*C_tG*C_Qu8  +118.94*C_tG*C_td8  +295.1*C_tG*C_tq8  +178.62*C_tG*C_tu8  +295.8*C_tG*C_Qq18  +62.33*C_tG*C_Qq38  +52.56*C_Qd8*C_td8  +237.4*C_Qd1*C_td1  +355.4*C_Qu1*C_tu1  +79.06*C_Qu8*C_tu8  +586.9*C_tq1*C_Qq11  +123.75*C_tq1*C_Qq31  +130.86*C_tq8*C_Qq18  +27.54*C_tq8*C_Qq38  +486.95*C_Qq11*C_Qq31  +108.06*C_Qq18*C_Qq38 ; 

     return SM_Charge_Asymmetry_bin_tt_500_750 + Delta/Total;
 } 
else{ 

    double Delta = Delta_Madgraph+60.0*C_tG+190.0*C_td8+462.2*C_tu8+659.8*C_Qq18+264.5*C_Qq38 ; 

    double Total = Total_Madgraph+264.5*C_tG+264.5*C_Qd8+264.5*C_Qu8+264.5*C_td8+264.5*C_tq8+264.5*C_tu8+264.5*C_Qq18+264.5*C_Qq38 ; 

    return SM_Charge_Asymmetry_bin_tt_500_750 + Delta/Total;
 } 

} 

Charge_Asymmetry_bin_tt_750_1000::Charge_Asymmetry_bin_tt_750_1000(const StandardModel& SM_i)
 : ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i)) 
{}

double    Charge_Asymmetry_bin_tt_750_1000::computeThValue() 
{ 
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double SM_Charge_Asymmetry_bin_tt_750_1000 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_750_1000();
    bool   flag_Quadratic = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic(); 
    double Delta_Madgraph = 30.0;
    double Total_Madgraph = 36990.0;


if(flag_Quadratic){ 

    double Delta = Delta_Madgraph  +7.0*C_tG  +120.4*C_td1*C_td1  +93.6*C_td8  +26.69*C_td8*C_td8  +296.8*C_tu1*C_tu1  +235.3*C_tu8  +65.94*C_tu8*C_tu8  +417.1*C_Qq11*C_Qq11  +417.2*C_Qq31*C_Qq31  +342.4*C_Qq18  +92.76*C_Qq18*C_Qq18  +146.81*C_Qq38  +92.82*C_Qq38*C_Qq38  +9.81*C_tG*C_td8  +24.33*C_tG*C_tu8  +33.13*C_tG*C_Qq18  +14.078*C_tG*C_Qq38  +0.09*C_Qd8*C_td8  +0.28*C_Qu1*C_tu1  +0.76*C_tq1*C_Qq31  +0.06*C_tq8*C_Qq18  +353.41*C_Qq11*C_Qq31  +78.23*C_Qq18*C_Qq38 ; 

    double Total = Total_Madgraph  +8901.0*C_tG  +1954.6*C_tG*C_tG  +378.8*C_Qd8  +92.52*C_Qd8*C_Qd8  +416.2*C_Qd1*C_Qd1  +668.8*C_Qu1*C_Qu1  +616.3*C_Qu8  +148.61*C_Qu8*C_Qu8  +416.2*C_td1*C_td1  +376.6*C_td8  +92.51*C_td8*C_td8  +1079.6*C_tq1*C_tq1  +983.9*C_tq8  +240.02*C_tq8*C_tq8  +668.6*C_tu1*C_tu1  +599.5*C_tu8  +148.66*C_tu8*C_tu8  +1079.7*C_Qq11*C_Qq11  +1079.4*C_Qq31*C_Qq31  +985.0*C_Qq18  +240.04*C_Qq18*C_Qq18  +232.19*C_Qq38  +239.98*C_Qq38*C_Qq38  +51.85*C_tG*C_Qd8  +83.37*C_tG*C_Qu8  +51.77*C_tG*C_td8  +134.95*C_tG*C_tq8  +83.25*C_tG*C_tu8  +134.17*C_tG*C_Qq18  +31.762*C_tG*C_Qq38  +22.99*C_Qd8*C_td8  +103.38*C_Qd1*C_td1  +166.54*C_Qu1*C_tu1  +36.88*C_Qu8*C_tu8  +267.4*C_tq1*C_Qq11  +63.32*C_tq1*C_Qq31  +59.28*C_tq8*C_Qq18  +14.128*C_tq8*C_Qq38  +516.59*C_Qq11*C_Qq31  +114.39*C_Qq18*C_Qq38 ; 

     return SM_Charge_Asymmetry_bin_tt_750_1000 + Delta/Total;
 } 
else{ 

    double Delta = Delta_Madgraph+7.0*C_tG+93.6*C_td8+235.3*C_tu8+342.4*C_Qq18+146.81*C_Qq38 ; 

    double Total = Total_Madgraph+146.81*C_tG+146.81*C_Qd8+146.81*C_Qu8+146.81*C_td8+146.81*C_tq8+146.81*C_tu8+146.81*C_Qq18+146.81*C_Qq38 ; 

    return SM_Charge_Asymmetry_bin_tt_750_1000 + Delta/Total;
 } 

} 

Charge_Asymmetry_bin_tt_1000_1500::Charge_Asymmetry_bin_tt_1000_1500(const StandardModel& SM_i)
 : ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i)) 
{}

double    Charge_Asymmetry_bin_tt_1000_1500::computeThValue() 
{ 
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double SM_Charge_Asymmetry_bin_tt_1000_1500 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_1000_1500();
    bool   flag_Quadratic = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic(); 
    double Delta_Madgraph = 44.0;
    double Total_Madgraph = 12862.0;


if(flag_Quadratic){ 

    double Delta = Delta_Madgraph  +6.0*C_tG  +0.8*C_tG*C_tG  +189.8*C_td1*C_td1  +78.02*C_td8  +42.2*C_td8*C_td8  +484.0*C_tu1*C_tu1  +197.0*C_tu8  +107.49*C_tu8*C_tu8  +673.5*C_Qq11*C_Qq11  +673.5*C_Qq31*C_Qq31  +271.0*C_Qq18  +149.6*C_Qq18*C_Qq18  +133.21*C_Qq38  +149.7*C_Qq38*C_Qq38  +7.21*C_tG*C_td8  +19.26*C_tG*C_tu8  +27.0*C_tG*C_Qq18  +11.544*C_tG*C_Qq38  +0.3*C_Qu1*C_tu1  +0.08*C_Qu8*C_tu8  +0.01*C_tq8*C_Qq18  +0.038*C_tq8*C_Qq38  +589.2*C_Qq11*C_Qq31  +130.91*C_Qq18*C_Qq38 ; 

    double Total = Total_Madgraph  +3080.0*C_tG  +980.4*C_tG*C_tG  +278.52*C_Qd8  +138.24*C_Qd8*C_Qd8  +621.7*C_Qd1*C_Qd1  +1070.2*C_Qu1*C_Qu1  +465.0*C_Qu8  +237.85*C_Qu8*C_Qu8  +621.8*C_td1*C_td1  +272.58*C_td8  +138.16*C_td8*C_td8  +1686.6*C_tq1*C_tq1  +742.6*C_tq8  +374.9*C_tq8*C_tq8  +1070.2*C_tu1*C_tu1  +460.4*C_tu8  +237.71*C_tu8*C_tu8  +1686.5*C_Qq11*C_Qq11  +1686.5*C_Qq31*C_Qq31  +732.4*C_Qq18  +374.8*C_Qq18*C_Qq18  +223.59*C_Qq38  +374.9*C_Qq38*C_Qq38  +38.97*C_tG*C_Qd8  +66.7*C_tG*C_Qu8  +39.13*C_tG*C_td8  +105.75*C_tG*C_tq8  +67.12*C_tG*C_tu8  +105.66*C_tG*C_Qq18  +27.736*C_tG*C_Qq38  +17.184*C_Qd8*C_td8  +77.58*C_Qd1*C_td1  +133.1*C_Qu1*C_tu1  +29.46*C_Qu8*C_tu8  +209.8*C_tq1*C_Qq11  +56.09*C_tq1*C_Qq31  +46.43*C_tq8*C_Qq18  +12.266*C_tq8*C_Qq38  +908.4*C_Qq11*C_Qq31  +201.29*C_Qq18*C_Qq38 ; 

     return SM_Charge_Asymmetry_bin_tt_1000_1500 + Delta/Total;
 } 
else{ 

    double Delta = Delta_Madgraph+6.0*C_tG+78.02*C_td8+197.0*C_tu8+271.0*C_Qq18+133.21*C_Qq38 ; 

    double Total = Total_Madgraph+133.21*C_tG+133.21*C_Qd8+133.21*C_Qu8+133.21*C_td8+133.21*C_tq8+133.21*C_tu8+133.21*C_Qq18+133.21*C_Qq38 ; 

    return SM_Charge_Asymmetry_bin_tt_1000_1500 + Delta/Total;
 } 

} 

Charge_Asymmetry_bin_tt_1500_2000::Charge_Asymmetry_bin_tt_1500_2000(const StandardModel& SM_i)
 : ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i)) 
{}

double    Charge_Asymmetry_bin_tt_1500_2000::computeThValue() 
{ 
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double SM_Charge_Asymmetry_bin_tt_1500_2000 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_1500_2000();
    bool   flag_Quadratic = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic(); 
    double Delta_Madgraph = -1.400000000000091;
    double Total_Madgraph = 1519.0;


if(flag_Quadratic){ 

    double Delta = Delta_Madgraph  +0.6*C_tG  +127.3*C_td1*C_td1  +26.33*C_td8  +28.28*C_td8*C_td8  +322.4*C_tu1*C_tu1  +65.61*C_tu8  +71.49*C_tu8*C_tu8  +449.8*C_Qq11*C_Qq11  +449.6*C_Qq31*C_Qq31  +96.74*C_Qq18  +99.97*C_Qq18*C_Qq18  +46.49*C_Qq38  +99.92*C_Qq38*C_Qq38  +2.194*C_tG*C_td8  +6.783*C_tG*C_tu8  +9.33*C_tG*C_Qq18  +3.806*C_tG*C_Qq38  +0.67*C_Qd1*C_td1  +0.163*C_Qu8*C_tu8  +389.5*C_Qq11*C_Qq31  +86.86*C_Qq18*C_Qq38 ; 

    double Total = Total_Madgraph  +365.6*C_tG  +193.64*C_tG*C_tG  +81.69*C_Qd8  +86.77*C_Qd8*C_Qd8  +390.8*C_Qd1*C_Qd1  +729.9*C_Qu1*C_Qu1  +161.18*C_Qu8  +162.26*C_Qu8*C_Qu8  +390.7*C_td1*C_td1  +81.07*C_td8  +86.88*C_td8*C_td8  +1118.8*C_tq1*C_tq1  +237.84*C_tq8  +248.78*C_tq8*C_tq8  +729.6*C_tu1*C_tu1  +150.79*C_tu8  +162.11*C_tu8*C_tu8  +1119.2*C_Qq11*C_Qq11  +1119.2*C_Qq31*C_Qq31  +230.46*C_Qq18  +248.63*C_Qq18*C_Qq18  +72.81*C_Qq38  +248.68*C_Qq38*C_Qq38  +12.235*C_tG*C_Qd8  +21.962*C_tG*C_Qu8  +11.538*C_tG*C_td8  +34.19*C_tG*C_tq8  +22.517*C_tG*C_tu8  +34.89*C_tG*C_Qq18  +10.308*C_tG*C_Qq38  +5.379*C_Qd8*C_td8  +23.45*C_Qd1*C_td1  +45.23*C_Qu1*C_tu1  +9.957*C_Qu8*C_tu8  +68.33*C_tq1*C_Qq11  +21.29*C_tq1*C_Qq31  +15.024*C_tq8*C_Qq18  +4.593*C_tq8*C_Qq38  +681.1*C_Qq11*C_Qq31  +151.74*C_Qq18*C_Qq38 ; 

     return SM_Charge_Asymmetry_bin_tt_1500_2000 + Delta/Total;
 } 
else{ 

    double Delta = Delta_Madgraph+0.6*C_tG+26.33*C_td8+65.61*C_tu8+96.74*C_Qq18+46.49*C_Qq38 ; 

    double Total = Total_Madgraph+46.49*C_tG+46.49*C_Qd8+46.49*C_Qu8+46.49*C_td8+46.49*C_tq8+46.49*C_tu8+46.49*C_Qq18+46.49*C_Qq38 ; 

    return SM_Charge_Asymmetry_bin_tt_1500_2000 + Delta/Total;
 } 

} 

Charge_Asymmetry_bin_tt_2000_2500::Charge_Asymmetry_bin_tt_2000_2500(const StandardModel& SM_i)
 : ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i)) 
{}

double    Charge_Asymmetry_bin_tt_2000_2500::computeThValue() 
{ 
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double SM_Charge_Asymmetry_bin_tt_2000_2500 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_2000_2500();
    bool   flag_Quadratic = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic(); 
    double Delta_Madgraph = 9.899999999999977;
    double Total_Madgraph = 280.5;


if(flag_Quadratic){ 

    double Delta = Delta_Madgraph  +2.58*C_tG  +83.31*C_td1*C_td1  +8.81*C_td8  +18.56*C_td8*C_td8  +200.2*C_tu1*C_tu1  +24.69*C_tu8  +44.51*C_tu8*C_tu8  +283.6*C_Qq11*C_Qq11  +283.3*C_Qq31*C_Qq31  +32.86*C_Qq18  +62.97*C_Qq18*C_Qq18  +9.395*C_Qq38  +62.95*C_Qq38*C_Qq38  +0.799*C_tG*C_td8  +2.474*C_tG*C_tu8  +3.302*C_tG*C_Qq18  +1.43*C_tG*C_Qq38  +0.284*C_Qu1*C_tu1  +1.339*C_tq1*C_Qq31  +0.009*C_tq8*C_Qq38  +234.5*C_Qq11*C_Qq31  +52.2*C_Qq18*C_Qq38 ; 

    double Total = Total_Madgraph  +66.9*C_tG  +50.09*C_tG*C_tG  +24.618*C_Qd8  +51.86*C_Qd8*C_Qd8  +233.32*C_Qd1*C_Qd1  +476.2*C_Qu1*C_Qu1  +61.12*C_Qu8  +105.75*C_Qu8*C_Qu8  +233.29*C_td1*C_td1  +28.23*C_td8  +51.88*C_td8*C_td8  +708.8*C_tq1*C_tq1  +89.13*C_tq8  +157.56*C_tq8*C_tq8  +476.2*C_tu1*C_tu1  +53.47*C_tu8  +105.91*C_tu8*C_tu8  +708.8*C_Qq11*C_Qq11  +708.9*C_Qq31*C_Qq31  +91.58*C_Qq18  +157.43*C_Qq18*C_Qq18  +28.885*C_Qq38  +157.45*C_Qq38*C_Qq38  +4.216*C_tG*C_Qd8  +8.61*C_tG*C_Qu8  +4.123*C_tG*C_td8  +12.681*C_tG*C_tq8  +8.282*C_tG*C_tu8  +12.864*C_tG*C_Qq18  +4.542*C_tG*C_Qq38  +1.845*C_Qd8*C_td8  +8.754*C_Qd1*C_td1  +16.088*C_Qu1*C_tu1  +3.876*C_Qu8*C_tu8  +25.7*C_tq1*C_Qq11  +8.963*C_tq1*C_Qq31  +5.87*C_tq8*C_Qq18  +1.99*C_tq8*C_Qq38  +486.5*C_Qq11*C_Qq31  +108.14*C_Qq18*C_Qq38 ; 

     return SM_Charge_Asymmetry_bin_tt_2000_2500 + Delta/Total;
 } 
else{ 

    double Delta = Delta_Madgraph+2.58*C_tG+8.81*C_td8+24.69*C_tu8+32.86*C_Qq18+9.395*C_Qq38 ; 

    double Total = Total_Madgraph+9.395*C_tG+9.395*C_Qd8+9.395*C_Qu8+9.395*C_td8+9.395*C_tq8+9.395*C_tu8+9.395*C_Qq18+9.395*C_Qq38 ; 

    return SM_Charge_Asymmetry_bin_tt_2000_2500 + Delta/Total;
 } 

} 

Charge_Asymmetry_bin_tt_2500_3000::Charge_Asymmetry_bin_tt_2500_3000(const StandardModel& SM_i)
 : ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i)) 
{}

double    Charge_Asymmetry_bin_tt_2500_3000::computeThValue() 
{ 
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double SM_Charge_Asymmetry_bin_tt_2500_3000 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_Charge_Asymmetry_bin_tt_2500_3000();
    bool   flag_Quadratic = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic(); 
    double Delta_Madgraph = -2.8999999999999986;
    double Total_Madgraph = 61.12;


if(flag_Quadratic){ 

    double Delta = Delta_Madgraph  +51.68*C_td1*C_td1  +3.262*C_td8  +11.478*C_td8*C_td8  +120.07*C_tu1*C_tu1  +9.458*C_tu8  +26.67*C_tu8*C_tu8  +171.8*C_Qq11*C_Qq11  +172.2*C_Qq31*C_Qq31  +13.305*C_Qq18  +38.22*C_Qq18*C_Qq18  +5.559*C_Qq38  +38.21*C_Qq38*C_Qq38  +0.564*C_tG*C_td8  +0.901*C_tG*C_tu8  +1.423*C_tG*C_Qq18  +0.386*C_tG*C_Qq38  +0.006*C_Qd8*C_td8  +0.624*C_Qu1*C_tu1  +0.232*C_Qu8*C_tu8  +0.243*C_tq1*C_Qq11  +0.001*C_tq8*C_Qq38  +135.8*C_Qq11*C_Qq31  +30.27*C_Qq18*C_Qq38 ; 

    double Total = Total_Madgraph  +14.606*C_tG  +15.074*C_tG*C_tG  +11.4*C_Qd8  +29.416*C_Qd8*C_Qd8  +132.52*C_Qd1*C_Qd1  +302.15*C_Qu1*C_Qu1  +21.356*C_Qu8  +67.13*C_Qu8*C_Qu8  +132.38*C_td1*C_td1  +9.424*C_td8  +29.482*C_td8*C_td8  +434.2*C_tq1*C_tq1  +38.76*C_tq8  +96.46*C_tq8*C_tq8  +302.13*C_tu1*C_tu1  +27.042*C_tu8  +67.15*C_tu8*C_tu8  +434.2*C_Qq11*C_Qq11  +434.4*C_Qq31*C_Qq31  +32.115*C_Qq18  +96.56*C_Qq18*C_Qq18  +12.195*C_Qq38  +96.51*C_Qq38*C_Qq38  +1.869*C_tG*C_Qd8  +3.491*C_tG*C_Qu8  +1.346*C_tG*C_td8  +5.49*C_tG*C_tq8  +3.599*C_tG*C_tu8  +5.263*C_tG*C_Qq18  +2.061*C_tG*C_Qq38  +0.653*C_Qd8*C_td8  +3.018*C_Qd1*C_td1  +7.298*C_Qu1*C_tu1  +1.591*C_Qu8*C_tu8  +10.329*C_tq1*C_Qq11  +3.613*C_tq1*C_Qq31  +2.332*C_tq8*C_Qq18  +0.958*C_tq8*C_Qq38  +339.2*C_Qq11*C_Qq31  +75.31*C_Qq18*C_Qq38 ; 

     return SM_Charge_Asymmetry_bin_tt_2500_3000 + Delta/Total;
 } 
else{ 

    double Delta = Delta_Madgraph+3.262*C_td8+9.458*C_tu8+13.305*C_Qq18+5.559*C_Qq38 ; 

    double Total = Total_Madgraph+5.559*C_tG+5.559*C_Qd8+5.559*C_Qu8+5.559*C_td8+5.559*C_tq8+5.559*C_tu8+5.559*C_Qq18+5.559*C_Qq38 ; 

    return SM_Charge_Asymmetry_bin_tt_2500_3000 + Delta/Total;
 } 

} 
*/



/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// ttll differential cross section for different bins //////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

sigma_ttll_diff_LO::sigma_ttll_diff_LO(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_ttll_bin_100_120" << "SM_sigma_ttll_bin_120_140" 
            << "SM_sigma_ttll_bin_140_180" << "SM_sigma_ttll_bin_180_500" );
                
}

double sigma_ttll_diff_LO::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_lQ1 = (C_lqP + C_lqM)/2.0;
    double C_lQ3 = (C_lqP - C_lqM)/2.0;
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    
    
    //std::cout<<"\033[1;31m C_lu = \033[0m "<< C_lu << std::endl;
    //std::cout<<"\033[1;31m C_eq = \033[0m "<< C_eq << std::endl;
    
    
    if(b_min == 100 && b_max == 120){
        
        double SM_ttll_bin_100_120 = SM.getOptionalParameter("SM_sigma_ttll_bin_100_120");
        
        //double ttll_bin_100_120_Madgraph = 1.;//fb
        
        double total;

        if(flag_Quadratic){
            
            total =  SM_ttll_bin_100_120 * 1/100 * ( 100 - 4*C_lQ1 + 0.11*C_lQ1*C_lQ1 + 4*C_lQ3 + 0.11*C_lQ3*C_lQ3 - 2.6*C_eu + 0.11 *C_eu*C_eu + 1.3*C_lu + 0.11*C_lu*C_lu + 2*C_eq + 0.11*C_eq*C_eq);  //Xsec mg5 included in coefficients
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
                     
        }
        else{
            
            total = SM_ttll_bin_100_120 * 1/100 * ( 100 - 4*C_lQ1 + 4*C_lQ3 - 2.6*C_eu + 1.3*C_lu + 2*C_eq );  //Xsec mg5 included in coefficients
        
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
                    
        }
 
    } else if(b_min == 120 && b_max == 140){
        
        double SM_ttll_bin_120_140 = SM.getOptionalParameter("SM_sigma_ttll_bin_120_140");
        
        //double ttll_bin_120_140_Madgraph = 0.34;//fb
        
        double total;
        
        if(flag_Quadratic){
            
            total = SM_ttll_bin_120_140 * 1/100 * ( 100 - 10*C_lQ1 + 0.5*C_lQ1*C_lQ1 + 10*C_lQ3 + 0.5*C_lQ3*C_lQ3 - 6*C_eu + 0.5 *C_eu*C_eu + 0.6*C_lu + 0.5*C_lu*C_lu + 2.3*C_eq + 0.5*C_eq*C_eq);  //Xsec mg5 included in coefficients
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        else{
            
            total = SM_ttll_bin_120_140 * 1/100 * ( 100 - 10*C_lQ1 + 10*C_lQ3 - 6*C_eu + 0.6*C_lu + 2.3*C_eq );  //Xsec mg5 included in coefficients
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
        }

    } else if(b_min == 140 && b_max == 180){
        
        double SM_ttll_bin_140_180 = SM.getOptionalParameter("SM_sigma_ttll_bin_140_180");
        
        //double ttll_bin_140_180_Madgraph = 0.23;//fb

        double total;
        
        if(flag_Quadratic){
            
            total = SM_ttll_bin_140_180 * 1/100 * ( 100 - 17*C_lQ1 + 1.6*C_lQ1*C_lQ1 + 17*C_lQ3 + 1.6*C_lQ3*C_lQ3 - 11*C_eu + 1.6 *C_eu*C_eu - 2*C_lu + 1.6*C_lu*C_lu + 1.2*C_eq + 1.6*C_eq*C_eq);  //Xsec mg5 included in coefficients
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        else{
            
            total = SM_ttll_bin_140_180 * 1/100 * ( 100 - 17*C_lQ1 + 17*C_lQ3 - 11*C_eu - 2*C_lu + 1.2*C_eq );  //Xsec mg5 included in coefficients
        
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }

    } else if(b_min == 180 && b_max == 500){
        
        double SM_ttll_bin_180_500 = SM.getOptionalParameter("SM_sigma_ttll_bin_180_500");
        
        //double ttll_bin_180_500_Madgraph = 0.26;//fb

        double total;
        
        if(flag_Quadratic){
            
            total = SM_ttll_bin_180_500 * 1/100 * ( 100 - 60*C_lQ1 + 70*C_lQ1*C_lQ1 + 60*C_lQ3 + 70*C_lQ3*C_lQ3 - 40*C_eu + 60 *C_eu*C_eu - 18*C_lu + 50*C_lu*C_lu - 7*C_eq + 60*C_eq*C_eq);  //Xsec mg5 included in coefficients
        
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }
        else{
            
            total = SM_ttll_bin_180_500 * 1/100 * ( 100 - 60*C_lQ1 + 60*C_lQ3 - 40*C_eu - 18*C_lu - 7*C_eq );  //Xsec mg5 included in coefficients        
            
            if (total < 0) return std::numeric_limits<double>::quiet_NaN();
            
            return total;
            
        }

    } else {
        throw std::runtime_error("\nERROR: Please specify a correct binning range for sigma_ttll_diff_LO.\n");
    }
    
}


entang_D_threshold::entang_D_threshold(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_entang_D_threshold");
}

double entang_D_threshold::computeThValue()
{
    
    double entang_D_threshold_SM = SM.getOptionalParameter("SM_entang_D_threshold");
    double entang_D_threshold_madgraph =  -0.467;//This is actually the best SM prediction
    
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();

    
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    //
    if(flag_Quadratic){
        return entang_D_threshold_SM + (-0.12093*C_tG + 0.03404*C_tG*C_tG + -0.00346*C_Qd1 + 
                0.00164*C_Qd1*C_Qd1 + 0.00294*C_Qd8 + 0.00114*C_Qd8*C_Qd8 + 0.00121*C_Qq11 + 
                0.00156*C_Qq11*C_Qq11 + -0.00095*C_Qq31 + -0.00044*C_Qq31*C_Qq31 + 0.00534*C_Qq18 + 
                -0.00074*C_Qq18*C_Qq18 + -0.00153*C_Qq38 + 9e-05*C_Qq38*C_Qq38 + -0.00014*C_Qu1 + 
                0.00041*C_Qu1*C_Qu1 + -0.00156*C_Qu8 + 0.00031*C_Qu8*C_Qu8 + -0.00217*C_td1 + 
                -0.00186*C_td1*C_td1 + 0.00669*C_td8 + 0.00143*C_td8*C_td8 + -7e-05*C_tq1 + 
                0.00019*C_tq1*C_tq1 + 0.00393*C_tq8 + 0.00042*C_tq8*C_tq8 + -0.00413*C_tu1 + 
                -0.00234*C_tu1*C_tu1 + 0.0036*C_tu8 + -0.00157*C_tu8*C_tu8)
                *entang_D_threshold_SM/entang_D_threshold_madgraph;

    }
    else{
        
        return  entang_D_threshold_SM + (-0.12093*C_tG + -0.00346*C_Qd1 + 0.00294*C_Qd8 + 
                0.00121*C_Qq11 + -0.00095*C_Qq31 + 0.00534*C_Qq18 + -0.00153*C_Qq38 + 
                -0.00014*C_Qu1 + -0.00156*C_Qu8 + -0.00217*C_td1 + 0.00669*C_td8 + 
                -7e-05*C_tq1 + 0.00393*C_tq8 + -0.00413*C_tu1 + 0.0036*C_tu8)
                *entang_D_threshold_SM/entang_D_threshold_madgraph;
        
    }
}



entang_Dn_boosted::entang_Dn_boosted(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_entang_Dn_boosted");
}

double entang_Dn_boosted::computeThValue()
{
    
    double entang_Dn_boosted_SM = SM.getOptionalParameter("SM_entang_Dn_boosted");
    double entang_Dn_boosted_madgraph =  -0.5873;//This is actually the best SM prediction
    
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double C_Qd1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd1();
    double C_Qu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu1();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_td1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td1();
    double C_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tq1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq1();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_tu1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu1();
    double C_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double C_Qq11 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq11();
    double C_Qq31 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq31();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();

    
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    //
    if(flag_Quadratic){
        return entang_Dn_boosted_SM + (-0.00324*C_tG + 0.11502*C_tG*C_tG + -0.0019*C_Qd1 + 
                0.01701*C_Qd1*C_Qd1 + -0.02474*C_Qd8 + 0.00494*C_Qd8*C_Qd8 + 0.01159*C_Qq11 + 
                0.03676*C_Qq11*C_Qq11 + -0.00447*C_Qq31 + 0.04024*C_Qq31*C_Qq31 + -0.02589*C_Qq18 + 
                0.01057*C_Qq18*C_Qq18 + -0.01194*C_Qq38 + 0.00976*C_Qq38*C_Qq38 + 0.00041*C_Qu1 + 
                0.03365*C_Qu1*C_Qu1 + -0.01233*C_Qu8 + 0.00434*C_Qu8*C_Qu8 + -0.0068*C_td1 + 
                0.01915*C_td1*C_td1 + -0.02005*C_td8 + 0.01168*C_td8*C_td8 + 0.00458*C_tq1 + 
                0.0354*C_tq1*C_tq1 + -0.02202*C_tq8 + 0.01814*C_tq8*C_tq8 + 0.00614*C_tu1 + 
                0.0362*C_tu1*C_tu1 + -0.00861*C_tu8 + 0.01378*C_tu8*C_tu8)
                *entang_Dn_boosted_SM/entang_Dn_boosted_madgraph;

    }
    else{
        
        return  entang_Dn_boosted_SM + (-0.00324*C_tG + -0.0019*C_Qd1 + -0.02474*C_Qd8 + 
                0.01159*C_Qq11 + -0.00447*C_Qq31 + -0.02589*C_Qq18 + -0.01194*C_Qq38 + 
                0.00041*C_Qu1 + -0.01233*C_Qu8 + -0.0068*C_td1 + -0.02005*C_td8 + 
                0.00458*C_tq1 + -0.02202*C_tq8 + 0.00614*C_tu1 + -0.00861*C_tu8)
                *entang_Dn_boosted_SM/entang_Dn_boosted_madgraph;
        
    }
}



/*
ttll_bin_100_120::ttll_bin_100_120(const StandardModel& SM_i)
 : ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i)) 
{}

double ttll_bin_100_120::computeThValue()
{
    
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_lQ1 = (C_lqP + C_lqM)/2.0;
    double C_lQ3 = (C_lqP - C_lqM)/2.0;
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();

    double SM_ttll_bin_100_120= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttll_bin_100_120(); 
 
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
//    double ttll_bin_100_120_madgraph_NLO = XXXXXXXXXXXX;//fb
//    double ttll_bin_100_120_madgraph_LO = XXXXXXXXXXX;//fb

    
    
    if(flag_Quadratic){
        return  SM_ttll_bin_100_120 * 1/100 * ( 100 - 4*C_lQ1 + 0.11*C_lQ1*C_lQ1 + 4*C_lQ3 + 0.11*C_lQ3*C_lQ3 - 2.6*C_eu + 0.11 *C_eu*C_eu + 1.3*C_lu + 0.11*C_lu*C_lu + 2*C_eq + 0.11*C_eq*C_eq);  //Xsec mg5 included in coefficients
    }
    else{
        return  SM_ttll_bin_100_120 * 1/100 * ( 100 - 4*C_lQ1 + 4*C_lQ3 - 2.6*C_eu + 1.3*C_lu + 2*C_eq );  //Xsec mg5 included in coefficients
    }
}


ttll_bin_120_140::ttll_bin_120_140(const StandardModel& SM_i)
 : ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i)) 
{}
double ttll_bin_120_140::computeThValue()
{
    
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_lQ1 = (C_lqP + C_lqM)/2.0;
    double C_lQ3 = (C_lqP - C_lqM)/2.0;
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();

    double SM_ttll_bin_120_140= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttll_bin_120_140(); 
 
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
//    double ttll_bin_120_140_madgraph_NLO = XXXXXXXXXXXX;//fb
//    double ttll_bin_120_140_madgraph_LO = XXXXXXXXXXX;//fb

    
    
    if(flag_Quadratic){
        return  SM_ttll_bin_120_140 * 1/100 * ( 100 - 10*C_lQ1 + 0.5*C_lQ1*C_lQ1 + 10*C_lQ3 + 0.5*C_lQ3*C_lQ3 - 6*C_eu + 0.5 *C_eu*C_eu + 0.6*C_lu + 0.5*C_lu*C_lu + 2.3*C_eq + 0.5*C_eq*C_eq);  //Xsec mg5 included in coefficients
    }
    else{
        return  SM_ttll_bin_120_140 * 1/100 * ( 100 - 10*C_lQ1 + 10*C_lQ3 - 6*C_eu + 0.6*C_lu + 2.3*C_eq );  //Xsec mg5 included in coefficients
    }
}


ttll_bin_140_180::ttll_bin_140_180(const StandardModel& SM_i)
 : ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i)) 
{}
double ttll_bin_140_180::computeThValue()
{
    
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_lQ1 = (C_lqP + C_lqM)/2.0;
    double C_lQ3 = (C_lqP - C_lqM)/2.0;
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();

    double SM_ttll_bin_140_180= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttll_bin_140_180(); 
 
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
//    double ttll_bin_140_180_madgraph_NLO = XXXXXXXXXXXX;//fb
//    double ttll_bin_140_180_madgraph_LO = XXXXXXXXXXX;//fb

    
    
    if(flag_Quadratic){
        return  SM_ttll_bin_140_180 * 1/100 * ( 100 - 17*C_lQ1 + 1.6*C_lQ1*C_lQ1 + 17*C_lQ3 + 1.6*C_lQ3*C_lQ3 - 11*C_eu + 1.6 *C_eu*C_eu - 2*C_lu + 1.6*C_lu*C_lu + 1.2*C_eq + 1.6*C_eq*C_eq);  //Xsec mg5 included in coefficients
    }
    else{
        return  SM_ttll_bin_140_180 * 1/100 * ( 100 - 17*C_lQ1 + 17*C_lQ3 - 11*C_eu - 2*C_lu + 1.2*C_eq );  //Xsec mg5 included in coefficients
    }
}



ttll_bin_180_500::ttll_bin_180_500(const StandardModel& SM_i)
 : ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i)) 
{}
double ttll_bin_180_500::computeThValue()
{
    
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_lQ1 = (C_lqP + C_lqM)/2.0;
    double C_lQ3 = (C_lqP - C_lqM)/2.0;
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();

    double SM_ttll_bin_180_500= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttll_bin_180_500(); 
 
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
//    double ttll_bin_180_500_madgraph_NLO = XXXXXXXXXXXX;//fb
//    double ttll_bin_180_500_madgraph_LO = XXXXXXXXXXX;//fb

       
    if(flag_Quadratic){
        return  SM_ttll_bin_180_500 * 1/100 * ( 100 - 60*C_lQ1 + 70*C_lQ1*C_lQ1 + 60*C_lQ3 + 70*C_lQ3*C_lQ3 - 40*C_eu + 60 *C_eu*C_eu - 18*C_lu + 50*C_lu*C_lu - 7*C_eq + 60*C_eq*C_eq);  //Xsec mg5 included in coefficients
    }
    else{
        return  SM_ttll_bin_180_500 * 1/100 * ( 100 - 60*C_lQ1 + 60*C_lQ3 - 40*C_eu - 18*C_lu - 7*C_eq );  //Xsec mg5 included in coefficients        
    }
}
*/










/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// OBSERVABLES FOR PROSPECTS OF FUTURE COLLIDERS ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/////// Prospects of muon Colliders at 3 TeV //////////////////////////////////////////////////////////



sigma_mumu_VBF_3TeV_tt::sigma_mumu_VBF_3TeV_tt(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_mumu_VBF_3TeV_tt");
}

double sigma_mumu_VBF_3TeV_tt::computeThValue()
{
    //WRITTEN IN PICO BARNS!!!
    double sigma_mumu_VBF_3TeV_tt_SM = SM.getOptionalParameter("SM_sigma_mumu_VBF_3TeV_tt");
    double sigma_mumu_VBF_3TeV_tt_madgraph = 0.013282;//pb
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    //
    if(flag_Quadratic){
        return sigma_mumu_VBF_3TeV_tt_SM + (-0.000360*C_phiQm+0.000234*C_phiQm*C_phiQm+-0.000392*C_phit+0.000235*C_phit*C_phit
                +0.018636*C_tW+0.057008*C_tW*C_tW+0.000071*C_tphi+0.000017*C_tphi*C_tphi+-0.015499*C_tZ+0.038741*C_tZ*C_tZ)
                *sigma_mumu_VBF_3TeV_tt_SM/sigma_mumu_VBF_3TeV_tt_madgraph;

    }
    else{
        
        return  sigma_mumu_VBF_3TeV_tt_SM + (-0.000360*C_phiQm+-0.000392*C_phit+0.018636*C_tW+0.000071*C_tphi+-0.015499*C_tZ)
                *sigma_mumu_VBF_3TeV_tt_SM/sigma_mumu_VBF_3TeV_tt_madgraph;
        
    }
}


sigma_mumu_3TeV_ttH::sigma_mumu_3TeV_ttH(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_mumu_3TeV_ttH");
}

double sigma_mumu_3TeV_ttH::computeThValue()
{

    //WRITTEN IN FEMTO BARNS!!!
    double sigma_mumu_3TeV_ttH_SM = SM.getOptionalParameter("SM_sigma_mumu_3TeV_ttH");
    double sigma_mumu_VBF_3TeV_ttH_madgraph = 0.031;//fb
    double sigma_mumu_sch_3TeV_ttH_madgraph = 0.412504;//fb
    double sigma_mumu_3TeV_ttH_madgraph = sigma_mumu_VBF_3TeV_ttH_madgraph + sigma_mumu_sch_3TeV_ttH_madgraph;//fb
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    //
    if(flag_Quadratic){
        double sigma_mumu_VBF_3TeV_ttH_NP = (+0.000914*C_phiQm+0.004713*C_phiQm*C_phiQm+-0.000548*C_phit+
        0.004561*C_phit*C_phit+0.137647*C_tW+3.129220*C_tW*C_tW+0.013882*C_phiQ3+0.026469*C_phiQ3*C_phiQ3+
        -0.003543*C_tphi+0.000744*C_tphi*C_tphi+-0.119718*C_tZ+2.331860*C_tZ*C_tZ);
        
        double sigma_mumu_sch_3TeV_ttH_NP = 0.412504+-0.052181*C_phiQm+0.023083*C_phiQm*C_phiQm+0.010818*C_phit+
        0.023120*C_phit*C_phit+-15.499400*C_lu+1103.620000*C_lu*C_lu+-23.593800*C_eu+1102.480000*C_eu*C_eu+
        12.385600*C_tW+262.948000*C_tW*C_tW+-0.043577*C_tphi+0.001296*C_tphi*C_tphi+-9.077840*C_tZ+214.569000*C_tZ*C_tZ;
        
        return sigma_mumu_3TeV_ttH_SM + (sigma_mumu_VBF_3TeV_ttH_NP + sigma_mumu_sch_3TeV_ttH_NP)*(sigma_mumu_3TeV_ttH_SM/sigma_mumu_3TeV_ttH_madgraph);

    }
    else{
        
        double sigma_mumu_VBF_3TeV_ttH_NP = +0.000914*C_phiQm+-0.000548*C_phit+
        0.137647*C_tW+0.013882*C_phiQ3+-0.003543*C_tphi+-0.119718*C_tZ;
        
        double sigma_mumu_sch_3TeV_ttH_NP = +-0.052181*C_phiQm+0.010818*C_phit+
        -15.499400*C_lu+-23.593800*C_eu+12.385600*C_tW+-0.043577*C_tphi+-9.077840*C_tZ;
        
        return sigma_mumu_3TeV_ttH_SM + (sigma_mumu_VBF_3TeV_ttH_NP + sigma_mumu_sch_3TeV_ttH_NP)*(sigma_mumu_3TeV_ttH_SM/sigma_mumu_3TeV_ttH_madgraph);

        
    }
    
    
}




sigma_mumu_3TeV_bb::sigma_mumu_3TeV_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_mumu_3TeV_bb");
}

double sigma_mumu_3TeV_bb::computeThValue()
{

    //WRITTEN IN FEMTO BARNS!!!
    double sigma_mumu_3TeV_bb_SM = SM.getOptionalParameter("SM_sigma_mumu_3TeV_bb");
    double sigma_mumu_VBF_3TeV_bb_madgraph =  0.753955;//fb
    double sigma_mumu_sch_3TeV_bb_madgraph = 10.155200;//fb
    double sigma_mumu_3TeV_bb_madgraph = sigma_mumu_VBF_3TeV_bb_madgraph + sigma_mumu_sch_3TeV_bb_madgraph;//fb
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    
    double C_bphi = 0;

    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    //
    if(flag_Quadratic){
        double sigma_mumu_VBF_3TeV_bb_NP = 0.753955+0.077635*C_bB+4.121930*C_bB*C_bB+-0.041003*C_bW+2.844950*C_bW*C_bW+
        0.009871*C_phiQ3+0.086978*C_phiQ3*C_phiQ3+//0.000480*C_bphi+0.003749*C_bphi*C_bphi
        +-0.004980*C_phib+0.105823*C_phib*C_phib+0.202236*(C_phiQm + C_phiQ3)+0.106541*(C_phiQm + C_phiQ3)*(C_phiQm + C_phiQ3);
        
        double sigma_mumu_sch_3TeV_bb_NP = +634.822000*C_ed+69725.400000*C_ed*C_ed+1497.570000*C_lqP
        +69725.500000*C_lqP*C_lqP+313.447000*C_ld+69722.600000*C_ld*C_ld+0.897097*C_phiQ3+
        0.030020*C_phiQ3*C_phiQ3+-0.097988*C_phib+0.030127*C_phib*C_phib+0.897097*(C_phiQm + C_phiQ3)
        +0.030020*(C_phiQm + C_phiQ3)*(C_phiQm + C_phiQ3);
        
        return sigma_mumu_3TeV_bb_SM + (sigma_mumu_VBF_3TeV_bb_NP + sigma_mumu_sch_3TeV_bb_NP)*(sigma_mumu_3TeV_bb_SM/sigma_mumu_3TeV_bb_madgraph);

    }
    else{
        
        double sigma_mumu_VBF_3TeV_bb_NP = 0.077635*C_bB+
        -0.041003*C_bW+0.009871*C_phiQ3+0.000480*C_bphi
        +-0.004980*C_phib+0.202236*(C_phiQm + C_phiQ3);
        
        double sigma_mumu_sch_3TeV_bb_NP = 634.822000*C_ed+
        1497.570000*C_lqP+313.447000*C_ld+0.897097*C_phiQ3+
        -0.097988*C_phib+0.897097*(C_phiQm + C_phiQ3);
        
        return sigma_mumu_3TeV_bb_SM + (sigma_mumu_VBF_3TeV_bb_NP + sigma_mumu_sch_3TeV_bb_NP)*(sigma_mumu_3TeV_bb_SM/sigma_mumu_3TeV_bb_madgraph);

        
    }
    
    
}




a_3TeV_mumu_bb::a_3TeV_mumu_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_3TeV_mumu_bb");
}

double a_3TeV_mumu_bb::computeThValue()
{
    
    //double a_250_bb_eP_P30_eM_M80 = 69.8;//in percentage//OLD
    //double a_250_bb_eP_P30_eM_M80 = 0.7038
    
    double asym_SM = SM.getOptionalParameter("SM_a_3TeV_mumu_bb");
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 8.3563;
    double NP_pos = 543.0139*ceb + -198.5954*ceQ + -0.2080*cHb + 0.6883*cHQ1 + 0.6883*cHQ3 + -106.3333*clb + 1296.7195*clq1 + 1296.7195*clq3;
    
    double SM_neg = 1.7989;
    double NP_neg = 91.8085*ceb + -120.4171*ceQ + 0.1100*cHb + 0.2088*cHQ1 + 0.2088*cHQ3 + 419.7807*clb + 200.8504*clq1 + 200.8504*clq3;    
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
    
}



sigma_mumu_VBF_10TeV_tt::sigma_mumu_VBF_10TeV_tt(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_mumu_VBF_10TeV_tt");
}

double sigma_mumu_VBF_10TeV_tt::computeThValue()
{
    //WRITTEN IN PICO BARNS!!!
    double sigma_mumu_VBF_10TeV_tt_SM = SM.getOptionalParameter("SM_sigma_mumu_VBF_10TeV_tt");
    double sigma_mumu_VBF_10TeV_tt_madgraph = 0.035865;//pb
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    //
    if(flag_Quadratic){
        return sigma_mumu_VBF_10TeV_tt_SM + (-0.000905*C_phiQm+0.002833*C_phiQm*C_phiQm+-0.002129*C_phit+0.002846*C_phit*C_phit
                +0.047246*C_tW+0.437635*C_tW*C_tW+0.000198*C_tphi+0.000084*C_tphi*C_tphi+-0.037538*C_tZ+0.240757*C_tZ*C_tZ)
                *sigma_mumu_VBF_10TeV_tt_SM/sigma_mumu_VBF_10TeV_tt_madgraph;

    }
    else{
        
        return  sigma_mumu_VBF_10TeV_tt_SM + (-0.000905*C_phiQm+-0.002129*C_phit+0.047246*C_tW+0.000198*C_tphi+-0.037538*C_tZ)
                *sigma_mumu_VBF_10TeV_tt_SM/sigma_mumu_VBF_10TeV_tt_madgraph;
        
    }
}





sigma_mumu_10TeV_ttH::sigma_mumu_10TeV_ttH(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_mumu_10TeV_ttH");
}

double sigma_mumu_10TeV_ttH::computeThValue()
{

    //WRITTEN IN FEMTO BARNS!!!
    
    double sigma_mumu_10TeV_ttH_SM = SM.getOptionalParameter("SM_sigma_mumu_10TeV_ttH");
    double sigma_mumu_VBF_10TeV_ttH_madgraph = 0.179932;//fb
    double sigma_mumu_sch_10TeV_ttH_madgraph = 0.053774;//fb
    double sigma_mumu_10TeV_ttH_madgraph = sigma_mumu_VBF_10TeV_ttH_madgraph + sigma_mumu_sch_10TeV_ttH_madgraph;//fb
    

    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    //
    if(flag_Quadratic){
        
        
        double sigma_mumu_VBF_10TeV_ttH_NP = 0.020219*C_phiQm+0.162560*C_phiQm*C_phiQm+
        -0.027879*C_phit+0.155512*C_phit*C_phit+0.551162*C_tW+70.819200*C_tW*C_tW+
        0.068450*C_phiQ3+0.991230*C_phiQ3*C_phiQ3+-0.012363*C_tphi+0.008227*C_tphi*C_tphi+
        -0.401485*C_tZ+51.466900*C_tZ*C_tZ;
        
        double sigma_mumu_sch_10TeV_ttH_NP = -0.008053*C_phiQm+0.021429*C_phiQm*C_phiQm+
        0.002200*C_phit+0.021416*C_phit*C_phit+-45.807400*C_lu+17609.100000*C_lu*C_lu+
        -69.072700*C_eu+17607.600000*C_eu*C_eu+17.775600*C_tW+2176.630000*C_tW*C_tW+
        -0.005346*C_tphi+0.000172*C_tphi*C_tphi+-3.075590*C_tZ+1773.620000*C_tZ*C_tZ;
        
        return sigma_mumu_10TeV_ttH_SM + (sigma_mumu_VBF_10TeV_ttH_NP + sigma_mumu_sch_10TeV_ttH_NP)*(sigma_mumu_10TeV_ttH_SM/sigma_mumu_10TeV_ttH_madgraph);

    }
    else{
        
        
        double sigma_mumu_VBF_10TeV_ttH_NP = 0.020219*C_phiQm+-0.027879*C_phit+
        0.551162*C_tW+0.068450*C_phiQ3+-0.012363*C_tphi+-0.401485*C_tZ;
        
        double sigma_mumu_sch_10TeV_ttH_NP = -0.008053*C_phiQm+0.002200*C_phit+
        -45.807400*C_lu+-69.072700*C_eu+17.775600*C_tW+-0.005346*C_tphi+-3.075590*C_tZ;
        
        return sigma_mumu_10TeV_ttH_SM + (sigma_mumu_VBF_10TeV_ttH_NP + sigma_mumu_sch_10TeV_ttH_NP)*(sigma_mumu_10TeV_ttH_SM/sigma_mumu_10TeV_ttH_madgraph);

        
        
    }
    
    
}




sigma_mumu_10TeV_bb::sigma_mumu_10TeV_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_mumu_10TeV_bb");
}

double sigma_mumu_10TeV_bb::computeThValue()
{

    //WRITTEN IN FEMTO BARNS!!!
    double sigma_mumu_10TeV_bb_SM = SM.getOptionalParameter("SM_sigma_mumu_10TeV_bb");
    double sigma_mumu_VBF_10TeV_bb_madgraph = 3.383220;//fb
    double sigma_mumu_sch_10TeV_bb_madgraph = 0.913030;//fb
    double sigma_mumu_10TeV_bb_madgraph = sigma_mumu_VBF_10TeV_bb_madgraph + sigma_mumu_sch_10TeV_bb_madgraph;//fb
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    
    double C_bphi = 0;
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    //
    if(flag_Quadratic){
        double sigma_mumu_VBF_10TeV_bb_NP = +-1.021790*C_bB+85.573800*C_bB*C_bB+
        -0.738307*C_bW+83.137200*C_bW*C_bW+0.213586*C_phiQ3+2.371940*C_phiQ3*C_phiQ3+
        0.002960*C_bphi+0.052864*C_bphi*C_bphi+
        -0.033040*C_phib+2.558190*C_phib*C_phib+1.644800*(C_phiQm + C_phiQ3)+
        2.558890*(C_phiQm + C_phiQ3)*(C_phiQm + C_phiQ3);
        
        double sigma_mumu_sch_10TeV_bb_NP = +568.500000*C_ed+774764.000000*C_ed*C_ed+
        1428.520000*C_lqP+774764.000000*C_lqP*C_lqP+267.750000*C_ld+774698.000000*C_ld*C_ld+
        0.080614*C_phiQ3+0.002692*C_phiQ3*C_phiQ3+-0.008836*C_phib+0.002719*C_phib*C_phib+
        0.080614*(C_phiQm + C_phiQ3)+0.002692*(C_phiQm + C_phiQ3)*(C_phiQm + C_phiQ3);
        
        return sigma_mumu_10TeV_bb_SM + (sigma_mumu_VBF_10TeV_bb_NP + sigma_mumu_sch_10TeV_bb_NP)*(sigma_mumu_10TeV_bb_SM/sigma_mumu_10TeV_bb_madgraph);

    }
    else{
        
        double sigma_mumu_VBF_10TeV_bb_NP = +-1.021790*C_bB+
        -0.738307*C_bW+0.213586*C_phiQ3+0.002960*C_bphi+
        -0.033040*C_phib+1.644800*(C_phiQm + C_phiQ3);
        
        double sigma_mumu_sch_10TeV_bb_NP = +568.500000*C_ed+1428.520000*C_lqP+
        267.750000*C_ld+0.080614*C_phiQ3+-0.008836*C_phib+0.080614*(C_phiQm + C_phiQ3);
        
        return sigma_mumu_10TeV_bb_SM + (sigma_mumu_VBF_10TeV_bb_NP + sigma_mumu_sch_10TeV_bb_NP)*(sigma_mumu_10TeV_bb_SM/sigma_mumu_10TeV_bb_madgraph);

        
    }
    
    
}





a_10TeV_mumu_bb::a_10TeV_mumu_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_10TeV_mumu_bb");
}

double a_10TeV_mumu_bb::computeThValue()
{
    
    //double a_250_bb_eP_P30_eM_M80 = 69.8;//in percentage//OLD
    //double a_250_bb_eP_P30_eM_M80 = 0.7038
    
    double asym_SM = SM.getOptionalParameter("SM_a_10TeV_mumu_bb");
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 0.7513;
    double NP_pos = 473.5391*ceb + -198.5954*ceQ + -0.0188*cHb + 0.0618*cHQ1 + 
    0.0618*cHQ3 + -1354.2705*clb + 1296.7195*clq1 + 1226.6425*clq3;    
    
    double SM_neg = 0.1617;
    double NP_neg = 94.9606*ceb + -120.4171*ceQ + 0.0099*cHb + 0.0188*cHQ1 + 
    0.0188*cHQ3 + 1622.0219*clb + 200.8504*clq1 + 201.8824*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
    
}





sigma_mumu_VBF_30TeV_tt::sigma_mumu_VBF_30TeV_tt(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_mumu_VBF_30TeV_tt");
}

double sigma_mumu_VBF_30TeV_tt::computeThValue()
{
    //WRITTEN IN PICO BARNS!!!
    double sigma_mumu_VBF_30TeV_tt_SM = SM.getOptionalParameter("SM_sigma_mumu_VBF_30TeV_tt");
    double sigma_mumu_VBF_30TeV_tt_madgraph = 0.060984;//pb
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    //
    if(flag_Quadratic){
        return sigma_mumu_VBF_30TeV_tt_SM + (-0.001369*C_phiQm+0.025047*C_phiQm*C_phiQm+-0.005211*C_phit+0.025065*C_phit*C_phit
                +0.080061*C_tW+3.675940*C_tW*C_tW+0.000316*C_tphi+0.000197*C_tphi*C_tphi+-0.058553*C_tZ+1.657510*C_tZ*C_tZ)
                *sigma_mumu_VBF_30TeV_tt_SM/sigma_mumu_VBF_30TeV_tt_madgraph;

    }
    else{
        
        return  sigma_mumu_VBF_30TeV_tt_SM + (-0.001369*C_phiQm+-0.005211*C_phit+0.080061*C_tW+0.000316*C_tphi+-0.058553*C_tZ)
                *sigma_mumu_VBF_30TeV_tt_SM/sigma_mumu_VBF_30TeV_tt_madgraph;
        
    }
}





sigma_mumu_30TeV_ttH::sigma_mumu_30TeV_ttH(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_mumu_30TeV_ttH");
}

double sigma_mumu_30TeV_ttH::computeThValue()
{

    //WRITTEN IN FEMTO BARNS!!!
    
    
    double sigma_mumu_30TeV_ttH_SM = SM.getOptionalParameter("SM_sigma_mumu_30TeV_ttH");
    double sigma_mumu_VBF_30TeV_ttH_madgraph = 0.493227;//fb
    double sigma_mumu_sch_30TeV_ttH_madgraph = 0.007052;//fb
    double sigma_mumu_30TeV_ttH_madgraph = sigma_mumu_VBF_30TeV_ttH_madgraph + sigma_mumu_sch_30TeV_ttH_madgraph;//fb
    
    
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    //
    if(flag_Quadratic){
        
        double sigma_mumu_VBF_30TeV_ttH_NP = 0.104163*C_phiQm+2.662680*C_phiQm*C_phiQm+-0.137693*C_phit+
        2.507270*C_phit*C_phit+0.519615*C_tW+1038.360000*C_tW*C_tW+0.086263*C_phiQ3+17.533500*C_phiQ3*C_phiQ3+
        -0.021297*C_tphi+0.055945*C_tphi*C_tphi+-1.452590*C_tZ+741.725000*C_tZ*C_tZ;
        
        double sigma_mumu_sch_30TeV_ttH_NP = -0.001124*C_phiQm+0.021091*C_phiQm*C_phiQm+0.000453*C_phit+
        0.021097*C_phit*C_phit+-179.774000*C_lu+203368.000000*C_lu*C_lu+-569.842000*C_eu+203518.000000*C_eu*C_eu+
        77.347500*C_tW+18710.900000*C_tW*C_tW+-0.000683*C_tphi+0.000026*C_tphi*C_tphi+5.975920*C_tZ+15250.800000*C_tZ*C_tZ;
        
        return sigma_mumu_30TeV_ttH_SM + (sigma_mumu_VBF_30TeV_ttH_NP + sigma_mumu_sch_30TeV_ttH_NP)*(sigma_mumu_30TeV_ttH_SM/sigma_mumu_30TeV_ttH_madgraph);


    }
    else{
        
        double sigma_mumu_VBF_30TeV_ttH_NP = 0.104163*C_phiQm+-0.137693*C_phit+
        0.519615*C_tW+0.086263*C_phiQ3+-0.021297*C_tphi+-1.452590*C_tZ;
        
        double sigma_mumu_sch_30TeV_ttH_NP = -0.001124*C_phiQm+0.000453*C_phit+
        -179.774000*C_lu+-569.842000*C_eu+77.347500*C_tW+-0.000683*C_tphi+5.975920*C_tZ;
        
        return sigma_mumu_30TeV_ttH_SM + (sigma_mumu_VBF_30TeV_ttH_NP + sigma_mumu_sch_30TeV_ttH_NP)*(sigma_mumu_30TeV_ttH_SM/sigma_mumu_30TeV_ttH_madgraph);

        
    }
    
    
}




sigma_mumu_30TeV_bb::sigma_mumu_30TeV_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_mumu_30TeV_bb");
}

double sigma_mumu_30TeV_bb::computeThValue()
{

    //WRITTEN IN FEMTO BARNS!!!
    double sigma_mumu_30TeV_bb_SM = SM.getOptionalParameter("SM_sigma_mumu_30TeV_bb");
    double sigma_mumu_VBF_30TeV_bb_madgraph = 6.613620;//fb
    double sigma_mumu_sch_30TeV_bb_madgraph = 0.101433;//fb
    double sigma_mumu_30TeV_bb_madgraph = sigma_mumu_VBF_30TeV_bb_madgraph + sigma_mumu_sch_30TeV_bb_madgraph;//fb
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    
    double C_bphi = 0;
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    //
    if(flag_Quadratic){
        double sigma_mumu_VBF_30TeV_bb_NP = +-222.935000*C_bB+832.885000*C_bB*C_bB+
        -16.151000*C_bW+1146.330000*C_bW*C_bW+7.830080*C_phiQ3+24.713400*C_phiQ3*C_phiQ3+
        -0.018559*C_bphi+0.131470*C_bphi*C_bphi+
        -0.084030*C_phib+25.041400*C_phib*C_phib+4.409190*(C_phiQm + C_phiQ3)+
        25.090300*(C_phiQm + C_phiQ3)*(C_phiQm + C_phiQ3);
        
        double sigma_mumu_sch_30TeV_bb_NP = +624.000000*C_ed+6972700.000000*C_ed*C_ed+
        1497.000000*C_lqP+6972710.000000*C_lqP*C_lqP+805.500000*C_ld+6971970.000000*C_ld*C_ld+
        0.008951*C_phiQ3+0.000303*C_phiQ3*C_phiQ3+-0.000978*C_phib+0.000301*C_phib*C_phib+
        0.008951*(C_phiQm + C_phiQ3)+0.000303*(C_phiQm + C_phiQ3)*(C_phiQm + C_phiQ3);
        
        return sigma_mumu_30TeV_bb_SM + (sigma_mumu_VBF_30TeV_bb_NP + sigma_mumu_sch_30TeV_bb_NP)*(sigma_mumu_30TeV_bb_SM/sigma_mumu_30TeV_bb_madgraph);

    }
    else{
        
        double sigma_mumu_VBF_30TeV_bb_NP = -222.935000*C_bB+-16.151000*C_bW+7.830080*C_phiQ3+
        -0.018559*C_bphi+
        -0.084030*C_phib+4.409190*(C_phiQm + C_phiQ3);
        
        double sigma_mumu_sch_30TeV_bb_NP = 624.000000*C_ed+1497.000000*C_lqP+
        805.500000*C_ld+0.008951*C_phiQ3+-0.000978*C_phib+0.008951*(C_phiQm + C_phiQ3);
        
        return sigma_mumu_30TeV_bb_SM + (sigma_mumu_VBF_30TeV_bb_NP + sigma_mumu_sch_30TeV_bb_NP)*(sigma_mumu_30TeV_bb_SM/sigma_mumu_30TeV_bb_madgraph);

        
    }
    
    
}



a_30TeV_mumu_bb::a_30TeV_mumu_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_30TeV_mumu_bb");
}

double a_30TeV_mumu_bb::computeThValue()
{
    
    //double a_250_bb_eP_P30_eM_M80 = 69.8;//in percentage//OLD
    //double a_250_bb_eP_P30_eM_M80 = 0.7038
    
    double asym_SM = SM.getOptionalParameter("SM_a_30TeV_mumu_bb");
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 0.0834;
    double NP_pos = -1522.0036*ceb + -13407.0537*ceQ + -0.0021*cHb + 0.0068*cHQ1 + 0.0068*cHQ3 + -15468.3957*clb + -380.9551*clq1 + -380.9551*clq3;    
    
    double SM_neg = 0.0180;
    double NP_neg = 2145.9985*ceb + 13632.0573*ceQ + 0.0011*cHb + 0.0021*cHQ1 + 0.0021*cHQ3 + 16273.8973*clb + 1877.9604*clq1 + 1877.9604*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
    
}



/////// Prospects of e+e- Colliders at 250 GeV //////////////////////////////////////////////////////////




sigma_250_bb_eP_M30_eM_M80::sigma_250_bb_eP_M30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_250_bb_eP_M30_eM_M80");
}

double sigma_250_bb_eP_M30_eM_M80::computeThValue()
{
    

    double sigma_250_bb_eP_M30_eM_M80_Madgraph = 1902.69;//fb
    
    double sigma_250_bb_eP_M30_eM_M80 = SM.getOptionalParameter("SM_sigma_250_bb_eP_M30_eM_M80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_250_bb_eP_M30_eM_M80 +(169.682*ceb + 125.2929*ceb*ceb + -113.88*ceQ + 
                125.2565*ceQ*ceQ + 19.548*cHb + 4.5671*cHb*cHb + 184.36*cHQ1 + 4.4895*cHQ1*cHQ1 + 
                184.36*cHQ3 + 4.4895*cHQ3*cHQ3 + 365.492*clb + 608.953*clb*clb + 2074.5568*clq1 + 
                609.1449*clq1*clq1 + 2074.5568*clq3 + 609.1449*clq3*clq3
                )*sigma_250_bb_eP_M30_eM_M80/sigma_250_bb_eP_M30_eM_M80_Madgraph;
    }
    else{
        return  sigma_250_bb_eP_M30_eM_M80 +(169.682*ceb + -113.88*ceQ + 19.548*cHb + 
                184.36*cHQ1 + 184.36*cHQ3 + 365.492*clb + 2074.5568*clq1 + 2074.5568*clq3
                )*sigma_250_bb_eP_M30_eM_M80/sigma_250_bb_eP_M30_eM_M80_Madgraph;
    }
}








sigma_250_bb_eP_M30_eM_P80::sigma_250_bb_eP_M30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_250_bb_eP_M30_eM_P80");
}

double sigma_250_bb_eP_M30_eM_P80::computeThValue()
{
    
    //double sigma_250_bb_eP_M30_eM_P80_Madgraph = 1018;//fb
    double sigma_250_bb_eP_M30_eM_P80_Madgraph = 949.58;//fb
    
    double sigma_250_bb_eP_M30_eM_P80 = SM.getOptionalParameter("SM_sigma_250_bb_eP_M30_eM_P80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_250_bb_eP_M30_eM_P80 +(1526.061*ceb + 1132.0629*ceb*ceb + -1026.08*ceQ + 
                1131.8329*ceQ*ceQ + -100.866*cHb + 5.7706*cHb*cHb + 89.832*cHQ1 + 5.7788*cHQ1*cHQ1 +
                89.832*cHQ3 + 5.7788*cHQ3*cHQ3 + 40.66*clb + 67.68*clb*clb + 230.568*clq1 + 
                67.7153*clq1*clq1 + 230.568*clq3 + 67.7153*clq3*clq3
                )*sigma_250_bb_eP_M30_eM_P80/sigma_250_bb_eP_M30_eM_P80_Madgraph;
    }
    else{
        return  sigma_250_bb_eP_M30_eM_P80 +(1526.061*ceb + -1026.08*ceQ + -100.866*cHb + 
                89.832*cHQ1 + 89.832*cHQ3 + 40.66*clb + 230.568*clq1 + 230.568*clq3
                )*sigma_250_bb_eP_M30_eM_P80/sigma_250_bb_eP_M30_eM_P80_Madgraph;
    }
}



sigma_250_bb_eP_P30_eM_M80::sigma_250_bb_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_250_bb_eP_P30_eM_M80");
}

double sigma_250_bb_eP_P30_eM_M80::computeThValue()
{
 
    
    //double sigma_250_bb_eP_P30_eM_M80_Madgraph = 3291;//fb
    double sigma_250_bb_eP_P30_eM_M80_Madgraph = 3424.2;//fb
    
    double sigma_250_bb_eP_P30_eM_M80 = SM.getOptionalParameter("SM_sigma_250_bb_eP_P30_eM_M80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    if(flag_Quadratic){
        return  sigma_250_bb_eP_P30_eM_M80 +(91.34*ceb + 66.8377*ceb*ceb + -61.31*ceQ + 66.7895*ceQ*ceQ + 
                51.316*cHb + 7.8376*cHb*cHb + 331.972*cHQ1 + 7.56*cHQ1*cHQ1 + 331.972*cHQ3 + 7.56*cHQ3*cHQ3 +
                678.772*clb + 1130.7153*clb*clb + 3852.8*clq1 + 1131.0753*clq1*clq1 + 3852.8*clq3 + 
                1131.0753*clq3*clq3
                )*sigma_250_bb_eP_P30_eM_M80/sigma_250_bb_eP_P30_eM_M80_Madgraph;
    }
    else{
        return  sigma_250_bb_eP_P30_eM_M80 +(91.34*ceb + -61.31*ceQ + 51.316*cHb + 
                331.972*cHQ1 + 331.972*cHQ3 + 678.772*clb + 3852.8*clq1 + 3852.8*clq3
                )*sigma_250_bb_eP_P30_eM_M80/sigma_250_bb_eP_P30_eM_M80_Madgraph;
    }
}





sigma_250_bb_eP_P30_eM_P80::sigma_250_bb_eP_P30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_250_bb_eP_P30_eM_P80");
}

double sigma_250_bb_eP_P30_eM_P80::computeThValue()
{
    

    double sigma_250_bb_eP_P30_eM_P80_Madgraph = 777.95;//fb
    
    double sigma_250_bb_eP_P30_eM_P80 = SM.getOptionalParameter("SM_sigma_250_bb_eP_P30_eM_P80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_250_bb_eP_P30_eM_P80 +(821.734*ceb + 609.4859*ceb*ceb + -552.428*ceQ + 
                609.2824*ceQ*ceQ + -49.746*cHb + 3.6365*cHb*cHb + 74.238*cHQ1 + 3.6459*cHQ1*cHQ1 +
                74.238*cHQ3 + 3.6459*cHQ3*cHQ3 + 75.452*clb + 125.6447*clb*clb + 428.1096*clq1 + 
                125.6898*clq1*clq1 + 428.1096*clq3 + 125.6898*clq3*clq3
                )*sigma_250_bb_eP_P30_eM_P80/sigma_250_bb_eP_P30_eM_P80_Madgraph;
    }
    else{
        return  sigma_250_bb_eP_P30_eM_P80 +(821.734*ceb + -552.428*ceQ + -49.746*cHb + 
                74.238*cHQ1 + 74.238*cHQ3 + 75.452*clb + 428.1096*clq1 + 428.1096*clq3
                )*sigma_250_bb_eP_P30_eM_P80/sigma_250_bb_eP_P30_eM_P80_Madgraph;
    }
}








a_250_bb_eP_M30_eM_M80::a_250_bb_eP_M30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_250_bb_eP_M30_eM_M80");
    
}

double a_250_bb_eP_M30_eM_M80::computeThValue()
{
    
    // double a_250_bb_eP_M30_eM_M80 = 0.69292
    double asym_SM = SM.getOptionalParameter("SM_a_250_bb_eP_M30_eM_M80");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 1610.5510;
    double NP_pos = 147.5905*ceb + -14.4829*ceQ + -5.7728*cHb + 153.9797*cHQ1 +
    153.9797*cHQ3 + 44.9577*clb + 1815.7625*clq1 + 1815.7625*clq3;
    
    double SM_neg = 292.1390;
    double NP_neg = 22.0915*ceb + -99.3971*ceQ + 25.3208*cHb + 30.3803*cHQ1 +
    30.3803*cHQ3 + 320.5343*clb + 258.7943*clq1 + 258.7943*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}









a_250_bb_eP_M30_eM_P80::a_250_bb_eP_M30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_250_bb_eP_M30_eM_P80");
    
}

double a_250_bb_eP_M30_eM_P80::computeThValue()
{
    
    // double a_250_bb_eP_M30_eM_P80 = 36.5;// In percentage //OLD
    // double a_250_bb_eP_M30_eM_P80 = 0.36826
    double asym_SM = SM.getOptionalParameter("SM_a_250_bb_eP_M30_eM_P80");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos =  649.6362;
    double NP_pos = 1336.7206*ceb + -123.8939*ceQ + -89.5699*cHb + 25.1274*cHQ1 + 25.1274*cHQ3 + 7.1719*clb + 201.8444*clq1 + 201.8444*clq3;
    
    double SM_neg = 299.9438;
    double NP_neg = 189.3404*ceb + -902.1861*ceQ + -11.2961*cHb + 64.7046*cHQ1 + 64.7046*cHQ3 + 33.4881*clb + 28.7236*clq1 + 28.7236*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}



a_250_bb_eP_P30_eM_M80::a_250_bb_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_250_bb_eP_P30_eM_M80");
}

double a_250_bb_eP_P30_eM_M80::computeThValue()
{
    
    //double a_250_bb_eP_P30_eM_M80 = 69.8;//in percentage//OLD
    //double a_250_bb_eP_P30_eM_M80 = 0.7038
    
    double asym_SM = SM.getOptionalParameter("SM_a_250_bb_eP_P30_eM_M80");
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 2917.1787;
    double NP_pos = 81.1444*ceb + -6.1462*ceQ + -0.9202*cHb + 289.6031*cHQ1 + 289.6031*cHQ3 + 89.2599*clb + 3373.1925*clq1 + 3373.1925*clq3;
    
    double SM_neg = 507.0213;
    double NP_neg = 10.1956*ceb + -55.1638*ceQ + 52.2362*cHb + 42.3689*cHQ1 + 42.3689*cHQ3 + 589.5121*clb + 479.6075*clq1 + 479.6075*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
    
}




a_250_bb_eP_P30_eM_P80::a_250_bb_eP_P30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_250_bb_eP_P30_eM_P80");
    
}

double a_250_bb_eP_P30_eM_P80::computeThValue()
{
    
    // double a_250_bb_eP_P30_eM_P80 = 0.48564 
    double asym_SM = SM.getOptionalParameter("SM_a_250_bb_eP_P30_eM_P80");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 577.8768;
    double NP_pos = 720.5272*ceb + -70.4976*ceQ + -47.8308*cHb + 
    37.3405*cHQ1 + 37.3405*cHQ3 + 9.5122*clb + 376.0983*clq1 + 376.0983*clq3;
    
    double SM_neg = 200.0732;
    double NP_neg = 101.2068*ceb + -481.9304*ceQ + -1.9152*cHb + 
    36.8975*cHQ1 + 36.8975*cHQ3 + 65.9398*clb + 52.0113*clq1 + 52.0113*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}






/////// Prospects of e+e- Colliders at 500 GeV //////////////////////////////////////////////////////////


sigma_500_bb_eP_M30_eM_M80::sigma_500_bb_eP_M30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_500_bb_eP_M30_eM_M80");
}

double sigma_500_bb_eP_M30_eM_M80::computeThValue()
{
    

    double sigma_500_bb_eP_M30_eM_M80_Madgraph = 413.01;//fb
    
    double sigma_500_bb_eP_M30_eM_M80 = SM.getOptionalParameter("SM_sigma_500_bb_eP_M30_eM_M80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    if(flag_Quadratic){
        return  sigma_500_bb_eP_M30_eM_M80 +(165.3448*ceb + 503.4372*ceb*ceb + -89.1808*ceQ + 
                503.1713*ceQ*ceQ + 4.9118*cHb + 0.9093*cHb*cHb + 38.088*cHQ1 + 0.9103*cHQ1*cHQ1 + 
                38.088*cHQ3 + 0.9103*cHQ3*cHQ3 + 390.696*clb + 2439.4*clb*clb + 1925.1396*clq1 + 
                2439.9085*clq1*clq1 + 1925.1396*clq3 + 2439.9085*clq3*clq3
                )*sigma_500_bb_eP_M30_eM_M80/sigma_500_bb_eP_M30_eM_M80_Madgraph;
    }
    else{
        return  sigma_500_bb_eP_M30_eM_M80 +(165.3448*ceb + -89.1808*ceQ + 4.9118*cHb + 
                38.088*cHQ1 + 38.088*cHQ3 + 390.696*clb + 1925.1396*clq1 + 1925.1396*clq3
                )*sigma_500_bb_eP_M30_eM_M80/sigma_500_bb_eP_M30_eM_M80_Madgraph;
    }
}





sigma_500_bb_eP_M30_eM_P80::sigma_500_bb_eP_M30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_500_bb_eP_M30_eM_P80");
}

double sigma_500_bb_eP_M30_eM_P80::computeThValue()
{
    
    //double sigma_500_bb_eP_M30_eM_P80 = 216.1;//fb
    double sigma_500_bb_eP_M30_eM_P80_Madgraph = 201.825;//fb
    
    double sigma_500_bb_eP_M30_eM_P80 = SM.getOptionalParameter("SM_sigma_500_bb_eP_M30_eM_P80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_500_bb_eP_M30_eM_P80 +(1488.424*ceb + 4531.4317*ceb*ceb + -803.144*ceQ + 
                4530.6906*ceQ*ceQ + -21.9908*cHb + 1.1384*cHb*cHb + 16.403*cHQ1 + 1.1397*cHQ1*cHQ1 + 
                16.403*cHQ3 + 1.1397*cHQ3*cHQ3 + 43.4182*clb + 270.9813*clb*clb + 213.9166*clq1 + 
                271.0719*clq1*clq1 + 213.9166*clq3 + 271.0719*clq3*clq3
                )*sigma_500_bb_eP_M30_eM_P80/sigma_500_bb_eP_M30_eM_P80_Madgraph;
    }
    else{
        return  sigma_500_bb_eP_M30_eM_P80 +(1488.424*ceb + -803.144*ceQ + -21.9908*cHb + 
                16.403*cHQ1 + 16.403*cHQ3 + 43.4182*clb + 213.9166*clq1 + 213.9166*clq3
                )*sigma_500_bb_eP_M30_eM_P80/sigma_500_bb_eP_M30_eM_P80_Madgraph;
    }
}







sigma_500_bb_eP_P30_eM_M80::sigma_500_bb_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_500_bb_eP_P30_eM_M80");
}

double sigma_500_bb_eP_P30_eM_M80::computeThValue()
{
    
    //double sigma_500_bb_eP_P30_eM_M80 = 718.5;//fb
    double sigma_500_bb_eP_P30_eM_M80_Madgraph = 743.86;//fb
    
    double sigma_500_bb_eP_P30_eM_M80 = SM.getOptionalParameter("SM_sigma_500_bb_eP_P30_eM_M80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    if(flag_Quadratic){
        return  sigma_500_bb_eP_P30_eM_M80 +(89.054*ceb + 270.9682*ceb*ceb + -47.964*ceQ + 
                270.8141*ceQ*ceQ + 12.484*cHb + 1.56*cHb*cHb + 68.968*cHQ1 + 1.5247*cHQ1*cHQ1 +
                68.968*cHQ3 + 1.5247*cHQ3*cHQ3 + 725.4*clb + 4530.5811*clb*clb + 3575.3236*clq1 +
                4531.2214*clq1*clq1 + 3575.3236*clq3 + 4531.2214*clq3*clq3
                )*sigma_500_bb_eP_P30_eM_M80/sigma_500_bb_eP_P30_eM_M80_Madgraph;
    }
    else{
        return  sigma_500_bb_eP_P30_eM_M80 +(89.054*ceb + -47.964*ceQ + 12.484*cHb + 
                68.968*cHQ1 + 68.968*cHQ3 + 725.4*clb + 3575.3236*clq1 + 3575.3236*clq3
                )*sigma_500_bb_eP_P30_eM_M80/sigma_500_bb_eP_P30_eM_M80_Madgraph;
    }
}







sigma_500_bb_eP_P30_eM_P80::sigma_500_bb_eP_P30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_500_bb_eP_P30_eM_P80");
}

double sigma_500_bb_eP_P30_eM_P80::computeThValue()
{
    

    double sigma_500_bb_eP_P30_eM_P80_Madgraph = 166.608;//fb
    
    double sigma_500_bb_eP_P30_eM_P80 = SM.getOptionalParameter("SM_sigma_500_bb_eP_P30_eM_P80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    if(flag_Quadratic){
        return  sigma_500_bb_eP_P30_eM_P80 +(801.364*ceb + 2440.1294*ceb*ceb + 
                -432.48*ceQ + 2439.5459*ceQ*ceQ + -10.7542*cHb + 0.7463*cHb*cHb + 
                14.2124*cHQ1 + 0.7428*cHQ1*cHQ1 + 14.2124*cHQ3 + 0.7428*cHQ3*cHQ3 +
                80.6244*clb + 503.3597*clb*clb + 397.259*clq1 + 503.481*clq1*clq1 +
                397.259*clq3 + 503.481*clq3*clq3
                )*sigma_500_bb_eP_P30_eM_P80/sigma_500_bb_eP_P30_eM_P80_Madgraph;
    }
    else{
        return  sigma_500_bb_eP_P30_eM_P80 +(801.364*ceb + -432.48*ceQ + -10.7542*cHb + 
                14.2124*cHQ1 + 14.2124*cHQ3 + 80.6244*clb + 397.259*clq1 + 397.259*clq3
                )*sigma_500_bb_eP_P30_eM_P80/sigma_500_bb_eP_P30_eM_P80_Madgraph;
    }
}





a_500_bb_eP_M30_eM_M80::a_500_bb_eP_M30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_500_bb_eP_M30_eM_M80");
}

double a_500_bb_eP_M30_eM_M80::computeThValue()
{
    
    
    //double a_500_bb_eP_M30_eM_M80 = 0.68064
    
    double asym_SM = SM.getOptionalParameter("SM_a_500_bb_eP_M30_eM_M80");
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 347.0606;
    double NP_pos = 144.2529*ceb + -9.4157*ceQ + -1.5618*cHb + 31.9791*cHQ1 +
    31.9791*cHQ3 + 46.5714*clb + 1686.8523*clq1 + 1686.8523*clq3;
    
    double SM_neg = 65.9494;
    double NP_neg = 21.0919*ceb + -79.7651*ceQ + 6.4736*cHb + 6.1089*cHQ1 + 
    6.1089*cHQ3 + 344.1246*clb + 238.2873*clq1 + 238.2873*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}




a_500_bb_eP_M30_eM_P80::a_500_bb_eP_M30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_500_bb_eP_M30_eM_P80");
}

double a_500_bb_eP_M30_eM_P80::computeThValue()
{
    
    //double a_500_bb_eP_M30_eM_P80 = 46.9; // in percentage 
    // double a_500_bb_eP_M30_eM_P80 = 0.47686
    double asym_SM = SM.getOptionalParameter("SM_a_500_bb_eP_M30_eM_P80");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos =  149.0336;
    double NP_pos = 1300.4371*ceb + -102.5159*ceQ + -19.9600*cHb + 
    5.0558*cHQ1 + 5.0558*cHQ3 + 5.3046*clb + 187.8792*clq1 + 187.8792*clq3;
    
    double SM_neg = 52.7914;
    double NP_neg = 187.9869*ceb + -700.6281*ceQ + -2.0308*cHb + 
    11.3472*cHQ1 + 11.3472*cHQ3 + 38.1136*clb + 26.0374*clq1 + 26.0374*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}


a_500_bb_eP_P30_eM_M80::a_500_bb_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_500_bb_eP_P30_eM_M80");
}

double a_500_bb_eP_P30_eM_M80::computeThValue()
{
    
    //double a_500_bb_eP_P30_eM_M80 = 68.4; // in percentage 
    //double a_500_bb_eP_P30_eM_M80 = 0.69352
    
    double asym_SM = SM.getOptionalParameter("SM_a_500_bb_eP_P30_eM_M80");
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 629.8709;
    double NP_pos = 78.2693*ceb + -5.7245*ceQ + 0.6528*cHb + 59.1902*cHQ1 + 
    59.1902*cHQ3 + 85.8688*clb + 3133.4382*clq1 + 3133.4382*clq3;
    
    double SM_neg = 113.9891;
    double NP_neg = 10.7847*ceb + -42.2395*ceQ + 11.8312*cHb + 9.7778*cHQ1 + 
    9.7778*cHQ3 + 639.5312*clb + 441.8854*clq1 + 441.8854*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}








a_500_bb_eP_P30_eM_P80::a_500_bb_eP_P30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_500_bb_eP_P30_eM_P80");
}

double a_500_bb_eP_P30_eM_P80::computeThValue()
{
    
    
    //double a_500_bb_eP_P30_eM_P80 = 0.55026
    
    double asym_SM = SM.getOptionalParameter("SM_a_500_bb_eP_P30_eM_P80");
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 129.1429;
    double NP_pos = 700.9852*ceb + -53.4024*ceQ + -10.7908*cHb + 
    7.2791*cHQ1 + 7.2791*cHQ3 + 9.2642*clb + 346.8097*clq1 + 346.8097*clq3;
    
    double SM_neg = 37.4651;
    double NP_neg = 100.3788*ceb + -379.0776*ceQ + 0.0366*cHb + 
    6.9333*cHQ1 + 6.9333*cHQ3 + 71.3602*clb + 50.4493*clq1 + 50.4493*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}



sigma_500_ttH_eP_P30_eM_M80::sigma_500_ttH_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_500_ttH_eP_P30_eM_M80");
}

double sigma_500_ttH_eP_P30_eM_M80::computeThValue()
{
    
    //double sigma_500_ttH_eP_P30_eM_M80 = 0.5176;//fb
    double sigma_500_ttH_eP_P30_eM_M80 = SM.getOptionalParameter("SM_sigma_500_ttH_eP_P30_eM_M80");
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    //double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phil1=0;  //For the moment we neglect these operators
    double C_phi3l1=0; //For the moment we neglect these operators
    double C_phi3l2=0; //For the moment we neglect these operators
    double C_phie=0;   //For the moment we neglect these operators

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+11.710899*C_phil1+1.523439*C_phil1*C_phil1-6.453904*C_phi3l1+1.282291*C_phi3l1*C_phi3l1-18.14602*C_phi3l2+0.768241*C_phi3l2*C_phi3l2+0.328835*C_phie+0.086283*C_phie*C_phie
                        -12.165806*C_tphi+0.43543*C_tphi*C_tphi-10.135556*C_phiQm+0.229592*C_phiQm*C_phiQm-9.743914*C_phit+0.224825*C_phit*C_phit-224.503686*C_tZ+140.841494*C_tZ*C_tZ+391.155313*C_tW+
                        389.310478*C_tW*C_tW-304.471127*C_lqM+239.750591*C_lqM*C_lqM-11.187638*C_eq+14.124559*C_eq*C_eq-295.376291*C_lu+240.656233*C_lu*C_lu-11.642502*C_eu+14.034138*C_eu*C_eu)*sigma_500_ttH_eP_P30_eM_M80/100;

            }
            else{
                return  (100+11.710899*C_phil1-6.453904*C_phi3l1-18.14602*C_phi3l2+0.328835*C_phie
                        -12.165806*C_tphi-10.135556*C_phiQm-9.743914*C_phit-224.503686*C_tZ+
                        391.155313*C_tW-304.471127*C_lqM-11.187638*C_eq-295.376291*C_lu-
                        11.642502*C_eu)*sigma_500_ttH_eP_P30_eM_M80/100;

            
            }
    }
//    ///////Needs to be rewritten in the old basis. NOT finished!!!
    else{
            if(flag_Quadratic){
                return  (100+11.710899*C_phil1-6.453904*C_phi3l1-18.14602*C_phi3l2+0.328835*C_phie
                        -12.165806*C_tphi-10.135556*C_phiQm-9.743914*C_phit-224.503686*C_tZ+
                        391.155313*C_tW-304.471127*C_lqM-11.187638*C_eq-295.376291*C_lu-
                        11.642502*C_eu)*sigma_500_ttH_eP_P30_eM_M80/100;
            }
            else{
                return  (100+11.710899*C_phil1-6.453904*C_phi3l1-18.14602*C_phi3l2+0.328835*C_phie
                        -12.165806*C_tphi-10.135556*C_phiQm-9.743914*C_phit-224.503686*C_tZ+
                        391.155313*C_tW-304.471127*C_lqM-11.187638*C_eq-295.376291*C_lu-
                        11.642502*C_eu)*sigma_500_ttH_eP_P30_eM_M80/100;
            }
    }
}


sigma_500_ttH_eP_M30_eM_P80::sigma_500_ttH_eP_M30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_500_ttH_eP_M30_eM_P80");
}

double sigma_500_ttH_eP_M30_eM_P80::computeThValue()
{
    
    //double sigma_500_ttH_eP_M30_eM_P80 = 0.2442;//fb
    double sigma_500_ttH_eP_M30_eM_P80 = SM.getOptionalParameter("SM_sigma_500_ttH_eP_M30_eM_P80");
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    //double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phil1=0;  //For the moment we neglect these operators
    double C_phi3l1=0; //For the moment we neglect these operators
    double C_phi3l2=0; //For the moment we neglect these operators
    double C_phie=0;   //For the moment we neglect these operators

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+1.450862*C_phil1+0.205376*C_phil1*C_phil1-16.664074*C_phi3l1+0.910369*C_phi3l1*C_phi3l1-18.102474*C_phi3l2+0.837449*C_phi3l2*C_phi3l2+11.966374*C_phie+3.38758*C_phie*C_phie
                        -12.365729*C_tphi+0.445364*C_tphi*C_tphi+10.783755*C_phiQm+0.547252*C_phiQm*C_phiQm+11.436094*C_phit+
                        0.546187*C_phit*C_phit-627.346445*C_tZ+1061.805178*C_tZ*C_tZ+566.038952*C_tW+825.846302*C_tW*C_tW-
                        38.501166*C_lqM+30.35687*C_lqM*C_lqM-403.957822*C_eq+509.751595*C_eq*C_eq-37.269149*C_lu+
                        30.315824*C_lu*C_lu-421.296882*C_eu+508.252377*C_eu*C_eu)*sigma_500_ttH_eP_M30_eM_P80/100;

            }
            else{
                return  (100+1.450862*C_phil1-16.664074*C_phi3l1-18.102474*C_phi3l2+11.966374*C_phie
                        -12.365729*C_tphi+10.783755*C_phiQm+11.436094*C_phit-627.346445*C_tZ+
                        566.038952*C_tW-38.501166*C_lqM-403.957822*C_eq-37.269149*C_lu-
                        421.296882*C_eu)*sigma_500_ttH_eP_M30_eM_P80/100;            
            }
    }
    
    ///////Needs to be rewritten in the old basis. NOT finished!!!
    else{
            if(flag_Quadratic){
                return  (100+1.450862*C_phil1-16.664074*C_phi3l1-18.102474*C_phi3l2+11.966374*C_phie
                        -12.365729*C_tphi+10.783755*C_phiQm+11.436094*C_phit-627.346445*C_tZ+
                        566.038952*C_tW-38.501166*C_lqM-403.957822*C_eq-37.269149*C_lu-
                        421.296882*C_eu)*sigma_500_ttH_eP_M30_eM_P80/100;            
            }
            else{
                return  (100+1.450862*C_phil1-16.664074*C_phi3l1-18.102474*C_phi3l2+11.966374*C_phie
                        -12.365729*C_tphi+10.783755*C_phiQm+11.436094*C_phit-627.346445*C_tZ+
                        566.038952*C_tW-38.501166*C_lqM-403.957822*C_eq-37.269149*C_lu-
                        421.296882*C_eu)*sigma_500_ttH_eP_M30_eM_P80/100;            
            }
    }
}



/////// Prospects of e+e- Colliders at 1000 GeV //////////////////////////////////////////////////////////



sigma_1000_bb_eP_P30_eM_M80::sigma_1000_bb_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1000_bb_eP_P30_eM_M80");
}

double sigma_1000_bb_eP_P30_eM_M80::computeThValue()
{
    
    //double sigma_1000_bb_eP_P30_eM_M80 = 173.9;//fb
    double sigma_1000_bb_eP_P30_eM_M80_Madgraph = 180.148;//fb
    
    double sigma_1000_bb_eP_P30_eM_M80 = SM.getOptionalParameter("SM_sigma_1000_bb_eP_P30_eM_M80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    if(flag_Quadratic){
        return  sigma_1000_bb_eP_P30_eM_M80 +(88.518*ceb + 1084.6219*ceb*ceb + -45.1278*ceQ +
                1084.3773*ceQ*ceQ + 3.0804*cHb + 0.3708*cHb*cHb + 16.508*cHQ1 + 0.3588*cHQ1*cHQ1 +
                16.508*cHQ3 + 0.3588*cHQ3*cHQ3 + 735.946*clb + 18126.6812*clb*clb + 3514.67*clq1 +
                18127.8318*clq1*clq1 + 3514.67*clq3 + 18127.8318*clq3*clq3
                )*sigma_1000_bb_eP_P30_eM_M80/sigma_1000_bb_eP_P30_eM_M80_Madgraph;
    }
    else{
        return  sigma_1000_bb_eP_P30_eM_M80 +(88.518*ceb + -45.1278*ceQ + 3.0804*cHb + 
                16.508*cHQ1 + 16.508*cHQ3 + 735.946*clb + 3514.67*clq1 + 3514.67*clq3
                )*sigma_1000_bb_eP_P30_eM_M80/sigma_1000_bb_eP_P30_eM_M80_Madgraph;
    }
}



sigma_1000_bb_eP_M30_eM_P80::sigma_1000_bb_eP_M30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1000_bb_eP_M30_eM_P80");
}

double sigma_1000_bb_eP_M30_eM_P80::computeThValue()
{
    
    //double sigma_1000_bb_eP_M30_eM_P80 = 52.2;//fb
    double sigma_1000_bb_eP_M30_eM_P80_Madgraph = 48.719;//fb
    
    double sigma_1000_bb_eP_M30_eM_P80 = SM.getOptionalParameter("SM_sigma_1000_bb_eP_M30_eM_P80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_1000_bb_eP_M30_eM_P80 +(1480.272*ceb + 18127.8882*ceb*ceb + -753.812*ceQ +
                18126.6577*ceQ*ceQ + -5.3244*cHb + 0.2682*cHb*cHb + 3.8*cHQ1 + 0.2674*cHQ1*cHQ1 +
                3.8*cHQ3 + 0.2674*cHQ3*cHQ3 + 44.0332*clb + 1084.456*clb*clb + 213.9166*clq1 + 
                451.1966*clq1*clq1 + 210.21*clq3 + 1084.6901*clq3*clq3
                )*sigma_1000_bb_eP_M30_eM_P80/sigma_1000_bb_eP_M30_eM_P80_Madgraph;
    }
    else{
        return  sigma_1000_bb_eP_M30_eM_P80 +(1480.272*ceb + -753.812*ceQ + -5.3244*cHb +
                3.8*cHQ1 + 3.8*cHQ3 + 44.0332*clb + 213.9166*clq1 + 210.21*clq3
                )*sigma_1000_bb_eP_M30_eM_P80/sigma_1000_bb_eP_M30_eM_P80_Madgraph;
    }
}






a_1000_bb_eP_M20_eM_M80::a_1000_bb_eP_M20_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_1000_bb_eP_M20_eM_M80");
}

double a_1000_bb_eP_M20_eM_M80::computeThValue()
{
    

    //double a_1000_bb_eP_M20_eM_M80 = 0.68274
    
    double asym_SM = SM.getOptionalParameter("SM_a_1000_bb_eP_M20_eM_M80");
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 95.38443;
    double NP_pos = 133.41633*ceb + -8.99238*ceQ + -0.37146*cHb + 8.83380*cHQ1 +
    8.83380*cHQ3 + 53.93455*clb + 1898.03403*clq1 + 1898.03403*clq3;
    
    double SM_neg = 17.98357;
    double NP_neg = 18.32967*ceb + -68.35362*ceQ + 1.90766*cHb + 1.50920*cHQ1 +
    1.50920*cHQ3 + 398.93545*clb + 264.55597*clq1 + 264.55597*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}



a_1000_bb_eP_M20_eM_P80::a_1000_bb_eP_M20_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_1000_bb_eP_M20_eM_P80");
}

double a_1000_bb_eP_M20_eM_P80::computeThValue()
{
    
    // double a_1000_bb_eP_M20_eM_P80 = 0.50804
    double asym_SM = SM.getOptionalParameter("SM_a_1000_bb_eP_M20_eM_P80");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 35.67269;
    double NP_pos = 1187.77272*ceb + -68.28617*ceQ + -4.38749*cHb + 
    1.32403*cHQ1 + 1.32403*cHQ3 + 6.99540*clb + 210.14683*clq1 + 210.14683*clq3;
    
    double SM_neg = 11.63731;
    double NP_neg = 178.66128*ceb + -627.56583*ceQ + -0.48285*cHb + 
    2.39583*cHQ1 + 2.39583*cHQ3 + 43.33340*clb + 30.07837*clq1 + 30.07837*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}




a_1000_bb_eP_P20_eM_M80::a_1000_bb_eP_P20_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_1000_bb_eP_P20_eM_M80");
}

double a_1000_bb_eP_P20_eM_M80::computeThValue()
{
    

    //double a_1000_bb_eP_P20_eM_M80 = 0.6857
    
    double asym_SM = SM.getOptionalParameter("SM_a_1000_bb_eP_P20_eM_M80");
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 140.57052;
    double NP_pos = 88.00200*ceb + -6.84170*ceQ + 0.11383*cHb + 13.20309*cHQ1 +
    13.20309*cHQ3 + 85.26107*clb + 2833.42810*clq1 + 2833.42810*clq3;
    
    double SM_neg = 26.20948;
    double NP_neg = 13.18960*ceb + -44.73970*ceQ + 2.66737*cHb + 2.07251*cHQ1 +
    2.07251*cHQ3 + 594.09093*clb + 410.82990*clq1 + 410.82990*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}









a_1000_bb_eP_P20_eM_P80::a_1000_bb_eP_P20_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_1000_bb_eP_P20_eM_P80");
}

double a_1000_bb_eP_P20_eM_P80::computeThValue()
{
    

    //double a_1000_bb_eP_P20_eM_P80 = 0.55008
    
    double asym_SM = SM.getOptionalParameter("SM_a_1000_bb_eP_P20_eM_P80");
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 32.29824;
    double NP_pos = 792.24276*ceb + -54.64439*ceQ + -3.01977*cHb + 
    1.63100*cHQ1 + 1.63100*cHQ3 + 8.83698*clb + 315.89886*clq1 + 315.89886*clq3;
    
    double SM_neg = 9.37476;
    double NP_neg = 118.63724*ceb + -409.34961*ceQ + -0.03449*cHb + 
    1.78022*cHQ1 + 1.78022*cHQ3 + 66.67342*clb + 44.49374*clq1 + 44.49374*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}





sigma_1000_bb_eP_M20_eM_M80::sigma_1000_bb_eP_M20_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1000_bb_eP_M20_eM_M80");
}

double sigma_1000_bb_eP_M20_eM_M80::computeThValue()
{
    
    
    double sigma_1000_bb_eP_M20_eM_M80_Madgraph = 113.368;//fb
    
    double sigma_1000_bb_eP_M20_eM_M80 = SM.getOptionalParameter("SM_sigma_1000_bb_eP_M20_eM_M80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_1000_bb_eP_M20_eM_M80 +(151.746*ceb + 1859.3494*ceb*ceb + 
                -77.346*ceQ + 1859.0294*ceQ*ceQ + 1.5362*cHb + 0.2349*cHb*cHb + 
                10.343*cHQ1 + 0.2354*cHQ1*cHQ1 + 10.343*cHQ3 + 0.2354*cHQ3*cHQ3 + 
                452.87*clb + 11154.9024*clb*clb + 2162.59*clq1 + 11156.0247*clq1*clq1 + 
                2162.59*clq3 + 11156.0247*clq3*clq3
                )*sigma_1000_bb_eP_M20_eM_M80/sigma_1000_bb_eP_M20_eM_M80_Madgraph;
    }
    else{
        return  sigma_1000_bb_eP_M20_eM_M80 +(151.746*ceb + -77.346*ceQ + 1.5362*cHb + 
                10.343*cHQ1 + 10.343*cHQ3 + 452.87*clb + 2162.59*clq1 + 2162.59*clq3
                )*sigma_1000_bb_eP_M20_eM_M80/sigma_1000_bb_eP_M20_eM_M80_Madgraph;
    }
}


sigma_1000_bb_eP_M20_eM_P80::sigma_1000_bb_eP_M20_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1000_bb_eP_M20_eM_P80");
}

double sigma_1000_bb_eP_M20_eM_P80::computeThValue()
{
    
    
    double sigma_1000_bb_eP_M20_eM_P80_Madgraph = 47.31;//fb
    
    double sigma_1000_bb_eP_M20_eM_P80 = SM.getOptionalParameter("SM_sigma_1000_bb_eP_M20_eM_P80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_1000_bb_eP_M20_eM_P80 +(1366.434*ceb + 16733.3329*ceb*ceb + 
                -695.852*ceQ + 16732.3694*ceQ*ceQ + -4.8703*cHb + 0.2519*cHb*cHb +
                3.7199*cHQ1 + 0.2543*cHQ1*cHQ1 + 3.7199*cHQ3 + 0.2543*cHQ3*cHQ3 + 
                50.3288*clb + 1239.392*clb*clb + 240.2252*clq1 + 1239.6402*clq1*clq1 + 
                240.2252*clq3 + 1239.6402*clq3*clq3
                )*sigma_1000_bb_eP_M20_eM_P80/sigma_1000_bb_eP_M20_eM_P80_Madgraph;
    }
    else{
        return  sigma_1000_bb_eP_M20_eM_P80 +(1366.434*ceb + -695.852*ceQ + -4.8703*cHb + 
                3.7199*cHQ1 + 3.7199*cHQ3 + 50.3288*clb + 240.2252*clq1 + 240.2252*clq3
                )*sigma_1000_bb_eP_M20_eM_P80/sigma_1000_bb_eP_M20_eM_P80_Madgraph;
    }
}



sigma_1000_bb_eP_P20_eM_M80::sigma_1000_bb_eP_P20_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1000_bb_eP_P20_eM_M80");
}

double sigma_1000_bb_eP_P20_eM_M80::computeThValue()
{
    
    
    double sigma_1000_bb_eP_P20_eM_M80_Madgraph = 166.78;//fb
    
    double sigma_1000_bb_eP_P20_eM_M80 = SM.getOptionalParameter("SM_sigma_1000_bb_eP_P20_eM_M80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    if(flag_Quadratic){
        return  sigma_1000_bb_eP_P20_eM_M80 +(101.1916*ceb + 1239.5508*ceb*ceb + 
                -51.5814*ceQ + 1239.3309*ceQ*ceQ + 2.7812*cHb + 0.3471*cHb*cHb + 
                15.2756*cHQ1 + 0.3485*cHQ1*cHQ1 + 15.2756*cHQ3 + 0.3485*cHQ3*cHQ3 + 
                679.352*clb + 16732.4235*clb*clb + 3244.258*clq1 + 16733.3306*clq1*clq1 + 
                3244.258*clq3 + 16733.3306*clq3*clq3
                )*sigma_1000_bb_eP_P20_eM_M80/sigma_1000_bb_eP_P20_eM_M80_Madgraph;
    }
    else{
        return  sigma_1000_bb_eP_P20_eM_M80 +(101.1916*ceb + -51.5814*ceQ + 2.7812*cHb +
                15.2756*cHQ1 + 15.2756*cHQ3 + 679.352*clb + 3244.258*clq1 + 3244.258*clq3
                )*sigma_1000_bb_eP_P20_eM_M80/sigma_1000_bb_eP_P20_eM_M80_Madgraph;
    }
}





sigma_1000_bb_eP_P20_eM_P80::sigma_1000_bb_eP_P20_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1000_bb_eP_P20_eM_P80");
}

double sigma_1000_bb_eP_P20_eM_P80::computeThValue()
{
    
    
    double sigma_1000_bb_eP_P20_eM_P80_Madgraph = 41.673;//fb
    
    double sigma_1000_bb_eP_P20_eM_P80 = SM.getOptionalParameter("SM_sigma_1000_bb_eP_P20_eM_P80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_1000_bb_eP_P20_eM_P80 +(910.88*ceb + 11155.5824*ceb*ceb + 
                -463.994*ceQ + 11154.8965*ceQ*ceQ + -3.0543*cHb + 0.1928*cHb*cHb + 
                3.4112*cHQ1 + 0.1919*cHQ1*cHQ1 + 3.4112*cHQ3 + 0.1919*cHQ3*cHQ3 + 
                75.5104*clb + 1859.1527*clb*clb + 360.3926*clq1 + 1859.3867*clq1*clq1 +
                360.3926*clq3 + 1859.3867*clq3*clq3
                )*sigma_1000_bb_eP_P20_eM_P80/sigma_1000_bb_eP_P20_eM_P80_Madgraph;
    }
    else{
        return  sigma_1000_bb_eP_P20_eM_P80 +( 910.88*ceb + -463.994*ceQ + -3.0543*cHb + 
                3.4112*cHQ1 + 3.4112*cHQ3 + 75.5104*clb + 360.3926*clq1 + 360.3926*clq3
                )*sigma_1000_bb_eP_P20_eM_P80/sigma_1000_bb_eP_P20_eM_P80_Madgraph;
    }
}




a_1000_bb_eP_P30_eM_M80::a_1000_bb_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_1000_bb_eP_P30_eM_M80");
}

double a_1000_bb_eP_P30_eM_M80::computeThValue()
{
    
    //double a_1000_bb_eP_P30_eM_M80 = 67.5 ;//in percentage OLD
    //double a_1000_bb_eP_P30_eM_M80 = 0.6873
    
    double asym_SM = SM.getOptionalParameter("SM_a_1000_bb_eP_P30_eM_M80");
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 151.98186;
    double NP_pos = 77.08833*ceb + -4.81075*ceQ + 0.24021*cHb + 14.33980*cHQ1 + 
    14.33980*cHQ3 + 86.15264*clb + 3077.62694*clq1 + 3077.62694*clq3;
    
    double SM_neg = 28.16614;
    double NP_neg = 11.42967*ceb + -40.31705*ceQ + 2.84019*cHb + 2.16820*cHQ1 + 
    2.16820*cHQ3 + 649.79336*clb + 437.04306*clq1 + 437.04306*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}


a_1000_bb_eP_M30_eM_P80::a_1000_bb_eP_M30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_1000_bb_eP_M30_eM_P80");
}

double a_1000_bb_eP_M30_eM_P80::computeThValue()
{
    
    // double a_1000_bb_eP_M30_eM_P80 = 49.5;//in percentage OLD
    // double a_1000_bb_eP_M30_eM_P80 = 0.49428
    double asym_SM = SM.getOptionalParameter("SM_a_1000_bb_eP_M30_eM_P80");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 36.39991;
    double NP_pos = 1287.44260*ceb + -77.12024*ceQ + -4.75187*cHb + 
    1.19735*cHQ1 + 1.19735*cHQ3 + 5.85840*clb + 187.87917*clq1 + 183.07434*clq3;
    
    double SM_neg = 12.31909;
    double NP_neg = 192.82940*ceb + -676.69176*ceQ + -0.57253*cHb + 
    2.60263*cHQ1 + 2.60263*cHQ3 + 38.17480*clb + 26.03743*clq1 + 27.13566*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}


sigma_1000_ttH_eP_P30_eM_M80::sigma_1000_ttH_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1000_ttH_eP_P30_eM_M80");
}

double sigma_1000_ttH_eP_P30_eM_M80::computeThValue()
{
    
    //double sigma_1000_ttH_eP_P30_eM_M80 = 3.354;//fb
    double sigma_1000_ttH_eP_P30_eM_M80 = SM.getOptionalParameter("SM_sigma_1000_ttH_eP_P30_eM_M80");
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    //double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phil1=0;  //For the moment we neglect these operators
    double C_phi3l1=0; //For the moment we neglect these operators
    double C_phi3l2=0; //For the moment we neglect these operators
    double C_phie=0;   //For the moment we neglect these operators

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+126.696859*C_phil1+582.261806*C_phil1*C_phil1+108.527029*C_phi3l1+571.576784*C_phi3l1*C_phi3l1-18.230597*C_phi3l2+0.768095*C_phi3l2*C_phi3l2-4.754208*C_phie+33.591323*C_phie*C_phie
                        -11.708703*C_tphi+0.314572*C_tphi*C_tphi-16.331793*C_phiQm+0.858348*C_phiQm*C_phiQm-10.393621*C_phit+0.859152*C_phit*C_phit-474.907842*C_tZ+804.498626*C_tZ*C_tZ+810.277599*C_tW+2162.847936*C_tW*C_tW-1275.215459*C_lqM+4577.46177*C_lqM*C_lqM-30.307657*C_eq+264.842701*C_eq*C_eq-875.417447*C_lu+4596.835149*C_lu*C_lu-52.208967*C_eu+271.062735*C_eu*C_eu)*sigma_1000_ttH_eP_P30_eM_M80/100;

            }
            else{
                return  (100+126.696859*C_phil1+108.527029*C_phi3l1-18.230597*C_phi3l2-4.754208*C_phie
                        -11.708703*C_tphi-16.331793*C_phiQm-10.393621*C_phit-474.907842*C_tZ+810.277599*C_tW-1275.215459*C_lqM-30.307657*C_eq-875.417447*C_lu-52.208967*C_eu)*sigma_1000_ttH_eP_P30_eM_M80/100;
            }
    }
    ///////Needs to be rewritten in the old basis. NOT finished!!!
    else{
            if(flag_Quadratic){
                return  (100+126.696859*C_phil1+108.527029*C_phi3l1-18.230597*C_phi3l2-4.754208*C_phie
                        -11.708703*C_tphi-16.331793*C_phiQm-10.393621*C_phit-474.907842*C_tZ+810.277599*C_tW-1275.215459*C_lqM-30.307657*C_eq-875.417447*C_lu-52.208967*C_eu)*sigma_1000_ttH_eP_P30_eM_M80/100;
            }
            else{
                return  (100+126.696859*C_phil1+108.527029*C_phi3l1-18.230597*C_phi3l2-4.754208*C_phie
                        -11.708703*C_tphi-16.331793*C_phiQm-10.393621*C_phit-474.907842*C_tZ+810.277599*C_tW-1275.215459*C_lqM-30.307657*C_eq-875.417447*C_lu-52.208967*C_eu)*sigma_1000_ttH_eP_P30_eM_M80/100;
            }
    }
}


sigma_1000_ttH_eP_M30_eM_P80::sigma_1000_ttH_eP_M30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1000_ttH_eP_M30_eM_P80");
}

double sigma_1000_ttH_eP_M30_eM_P80::computeThValue()
{
    //double sigma_1000_ttH_eP_M30_eM_P80 = 1.740;//fb
    double sigma_1000_ttH_eP_M30_eM_P80 = SM.getOptionalParameter("SM_sigma_1000_ttH_eP_M30_eM_P80");
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    //double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phil1=0;  //For the moment we neglect these operators
    double C_phi3l1=0; //For the moment we neglect these operators
    double C_phi3l2=0; //For the moment we neglect these operators
    double C_phie=0;   //For the moment we neglect these operators

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+14.761071*C_phil1+61.693549*C_phil1*C_phil1-3.51419*C_phi3l1+61.225925*C_phi3l1*C_phi3l1-18.228451*C_phi3l2+0.935709*C_phi3l2*C_phi3l2-154.221151*C_phie+1135.13874*C_phie*C_phie
                        -11.595772*C_tphi+0.40181*C_tphi*C_tphi+9.168558*C_phiQm+1.471898*C_phiQm*C_phiQm+18.396741*C_phit+1.48195*C_phit*C_phit-1189.261469*C_tZ+5295.292352*C_tZ*C_tZ+1081.736982*C_tW+4193.88076*C_tW*C_tW-148.867258*C_lqM+535.647553*C_lqM*C_lqM-1026.667338*C_eq+8803.1345*C_eq*C_eq-99.185412*C_lu+508.815379*C_lu*C_lu-1705.860701*C_eu+8822.898534*C_eu*C_eu)*sigma_1000_ttH_eP_M30_eM_P80/100;

            }
            else{
                return  (100+14.761071*C_phil1-3.51419*C_phi3l1-18.228451*C_phi3l2-154.221151*C_phie
                        -11.595772*C_tphi+9.168558*C_phiQm+18.396741*C_phit-1189.261469*C_tZ+1081.736982*C_tW-148.867258*C_lqM-1026.667338*C_eq-99.185412*C_lu-1705.860701*C_eu)*sigma_1000_ttH_eP_M30_eM_P80/100;
            }
    }
    
    ///////Needs to be rewritten in the old basis. NOT finished!!!
    else{  
            if(flag_Quadratic){
                return  (100+14.761071*C_phil1-3.51419*C_phi3l1-18.228451*C_phi3l2-154.221151*C_phie
                        -11.595772*C_tphi+9.168558*C_phiQm+18.396741*C_phit-1189.261469*C_tZ+1081.736982*C_tW-148.867258*C_lqM-1026.667338*C_eq-99.185412*C_lu-1705.860701*C_eu)*sigma_1000_ttH_eP_M30_eM_P80/100;
            }
            else{
                return  (100+14.761071*C_phil1-3.51419*C_phi3l1-18.228451*C_phi3l2-154.221151*C_phie
                        -11.595772*C_tphi+9.168558*C_phiQm+18.396741*C_phit-1189.261469*C_tZ+1081.736982*C_tW-148.867258*C_lqM-1026.667338*C_eq-99.185412*C_lu-1705.860701*C_eu)*sigma_1000_ttH_eP_M30_eM_P80/100;
            }
    }
}



/////// Prospects of Linear Colders at 380 GeV //////////////////////////////////////////////////////////



sigma_380_bb_eP_0_eM_P80::sigma_380_bb_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_380_bb_eP_0_eM_P80");
}

double sigma_380_bb_eP_0_eM_P80::computeThValue()
{
    
    //double sigma_380_bb_eP_0_eM_P80 = 348;//fb
    double sigma_380_bb_eP_0_eM_P80_Madgraph = 330.45;//fb
    
    double sigma_380_bb_eP_0_eM_P80 = SM.getOptionalParameter("SM_sigma_380_bb_eP_0_eM_P80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_380_bb_eP_0_eM_P80 +(1151.4114*ceb + 2013.0309*ceb*ceb + 
                -656.354*ceQ + 2012.7518*ceQ*ceQ + -29.273*cHb + 1.7434*cHb*cHb +
                28.4148*cHQ1 + 1.7574*cHQ1*cHQ1 + 28.4148*cHQ3 + 1.7574*cHQ3*cHQ3 +
                61.1278*clb + 223.5679*clb*clb + 310.908*clq1 + 223.6499*clq1*clq1 +
                310.908*clq3 + 223.6499*clq3*clq3
                )*sigma_380_bb_eP_0_eM_P80/sigma_380_bb_eP_0_eM_P80_Madgraph;
    }
    else{
        return  sigma_380_bb_eP_0_eM_P80 +(1151.4114*ceb + -656.354*ceQ + -29.273*cHb +
                28.4148*cHQ1 + 28.4148*cHQ3 + 61.1278*clb + 310.908*clq1 + 310.908*clq3
                )*sigma_380_bb_eP_0_eM_P80/sigma_380_bb_eP_0_eM_P80_Madgraph;
    }
}




sigma_380_bb_eP_0_eM_M80::sigma_380_bb_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_380_bb_eP_0_eM_M80");
}

double sigma_380_bb_eP_0_eM_M80::computeThValue()
{
    
    //double sigma_380_bb_eP_0_eM_M80 = 998;//fb OLD
    double sigma_380_bb_eP_0_eM_M80_Madgraph = 1034.28;//fb
    
    double sigma_380_bb_eP_0_eM_M80 = SM.getOptionalParameter("SM_sigma_380_bb_eP_0_eM_M80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_380_bb_eP_0_eM_M80 +(127.954*ceb + 223.4694*ceb*ceb + -72.84*ceQ + 
                223.3082*ceQ*ceQ + 15.09*cHb + 2.3165*cHb*cHb + 96.824*cHQ1 + 2.2047*cHQ1*cHQ1 + 
                96.824*cHQ3 + 2.2047*cHQ3*cHQ3 + 550.04*clb + 2012.1106*clb*clb + 2798.244*clq1 +
                2012.8322*clq1*clq1 + 2798.244*clq3 + 2012.8322*clq3*clq3
                )*sigma_380_bb_eP_0_eM_M80/sigma_380_bb_eP_0_eM_M80_Madgraph;
    }
    else{
        return  sigma_380_bb_eP_0_eM_M80 +(127.954*ceb + -72.84*ceQ + 15.09*cHb +
                96.824*cHQ1 + 96.824*cHQ3 + 550.04*clb + 2798.244*clq1 + 2798.244*clq3
                )*sigma_380_bb_eP_0_eM_M80/sigma_380_bb_eP_0_eM_M80_Madgraph;
    }
}



a_380_bb_eP_0_eM_P80::a_380_bb_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_380_bb_eP_0_eM_P80");
}

double a_380_bb_eP_0_eM_P80::computeThValue()
{
    
    //double a_380_bb_eP_0_eM_P80 = 47.9;//in percentage OLD
    // double a_380_bb_eP_0_eM_P80 = 0.48654
    double asym_SM = SM.getOptionalParameter("SM_a_380_bb_eP_0_eM_P80");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 245.6135715;
    double NP_pos = 1006.1668132*ceb + -82.0819923*ceQ + -27.0470725*cHb + 
    10.9831380*cHQ1 + 10.9831380*cHQ3 + 7.7659007*clb + 272.5095266*clq1 + 272.5095266*clq3;
    
    double SM_neg = 84.8364285;
    double NP_neg = 145.2445868*ceb + -574.2720077*ceQ + -2.2259275*cHb + 
    17.4316620*cHQ1 + 17.4316620*cHQ3 + 53.3618992*clb + 38.3984734*clq1 + 38.3984734*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}




a_380_bb_eP_0_eM_M80::a_380_bb_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_380_bb_eP_0_eM_M80");
}

double a_380_bb_eP_0_eM_M80::computeThValue()
{
    
    //double a_380_bb_eP_0_eM_M80 = 67.93;//in percentage OLD
    // double a_380_bb_eP_0_eM_M80 = 0.6929;
    double asym_SM = SM.getOptionalParameter("SM_a_380_bb_eP_0_eM_M80");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 875.4663060;
    double NP_pos = 113.9170800*ceb + -8.0975402*ceQ + -0.3958228*cHb + 83.1980184*cHQ1 + 
    83.1980184*cHQ3 + 69.3365063*clb + 2449.3218223*clq1 + 2449.3218223*clq3;
    
    double SM_neg = 158.8136940;
    double NP_neg = 14.0369200*ceb + -64.7424600*ceQ + 15.4858227*cHb + 13.6259818*cHQ1 + 
    13.6259818*cHQ3 + 480.7034933*clb + 348.9221777*clq1 + 348.9221777*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;

}







/////// Prospects of Linear Colders at 1500 GeV //////////////////////////////////////////////////////////


sigma_1500_bb_eP_0_eM_P80::sigma_1500_bb_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1500_bb_eP_0_eM_P80");
}

double sigma_1500_bb_eP_0_eM_P80::computeThValue()
{
    
    //double sigma_1400_bb_eP_0_eM_P80 = 23.77;//fb OLD
    double sigma_1500_bb_eP_0_eM_P80_Madgraph = 19.6502;//fb
    
    double sigma_1500_bb_eP_0_eM_P80 = SM.getOptionalParameter("SM_sigma_1500_bb_eP_0_eM_P80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_1500_bb_eP_0_eM_P80 +(1137.42*ceb + 31376.2821*ceb*ceb + -572.12*ceQ + 
                31374.4586*ceQ*ceQ + -1.75*cHb + 0.0995*cHb*cHb + 1.5643*cHQ1 + 0.1006*cHQ1*cHQ1 + 
                1.5643*cHQ3 + 0.1006*cHQ3*cHQ3 + 63.178*clb + 3486.0762*clb*clb + 299.348*clq1 + 
                3486.3786*clq1*clq1 + 299.348*clq3 + 3486.3786*clq3*clq3
                )*sigma_1500_bb_eP_0_eM_P80/sigma_1500_bb_eP_0_eM_P80_Madgraph;
    }
    else{
        return  sigma_1500_bb_eP_0_eM_P80 +(1137.42*ceb + -572.12*ceQ + -1.75*cHb + 
                1.5643*cHQ1 + 1.5643*cHQ3 + 63.178*clb + 299.348*clq1 + 299.348*clq3
                )*sigma_1500_bb_eP_0_eM_P80/sigma_1500_bb_eP_0_eM_P80_Madgraph;
    }
}




sigma_1500_bb_eP_0_eM_M80::sigma_1500_bb_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1500_bb_eP_0_eM_M80");
}

double sigma_1500_bb_eP_0_eM_M80::computeThValue()
{
    
    //double sigma_1400_bb_eP_0_eM_M80 = 68.85;//fb
    double sigma_1500_bb_eP_0_eM_M80_Madgraph = 61.893;//fb
    
    double sigma_1500_bb_eP_0_eM_M80 = SM.getOptionalParameter("SM_sigma_1500_bb_eP_0_eM_M80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_1500_bb_eP_0_eM_M80 +(126.338*ceb + 3486.2941*ceb*ceb + -63.672*ceQ +
                3486.0176*ceQ*ceQ + 0.9598*cHb + 0.1316*cHb*cHb + 5.6478*cHQ1 + 0.1342*cHQ1*cHQ1 +
                5.6478*cHQ3 + 0.1342*cHQ3*cHQ3 + 568.5*clb + 31374.0435*clb*clb + 2695.22*clq1 + 
                31376.2553*clq1*clq1 + 2695.22*clq3 + 31376.2553*clq3*clq3
                )*sigma_1500_bb_eP_0_eM_M80/sigma_1500_bb_eP_0_eM_M80_Madgraph;
    }
    else{
        return  sigma_1500_bb_eP_0_eM_M80 +(126.338*ceb + -63.672*ceQ + 0.9598*cHb + 
                5.6478*cHQ1 + 5.6478*cHQ3 + 568.5*clb + 2695.22*clq1 + 2695.22*clq3
                )*sigma_1500_bb_eP_0_eM_M80/sigma_1500_bb_eP_0_eM_M80_Madgraph;
    }
}







a_1500_bb_eP_0_eM_P80::a_1500_bb_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_1500_bb_eP_0_eM_P80");
}

double a_1500_bb_eP_0_eM_P80::computeThValue()
{
    
    //double a_1400_bb_eP_0_eM_P80 = 51.8; //in percentage
    //double a_1500_bb_eP_0_eM_P80 = 0.5283
    double asym_SM = SM.getOptionalParameter("SM_a_1500_bb_eP_0_eM_P80");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 15.0157003;
    double NP_pos = 996.0566886*ceb + -70.8619634*ceQ + -1.6207960*cHb + 
    0.6387019*cHQ1 + 0.6387019*cHQ3 + 7.8897679*clb + 261.3672417*clq1 + 261.3672417*clq3;
    
    double SM_neg = 4.6344997;
    double NP_neg = 141.3633100*ceb + -501.2580366*ceQ + -0.1292040*cHb + 
    0.9256181*cHQ1 + 0.9256181*cHQ3 + 55.2882321*clb + 37.9807583*clq1 + 37.9807583*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}




a_1500_bb_eP_0_eM_M80::a_1500_bb_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_1500_bb_eP_0_eM_M80");
}

double a_1500_bb_eP_0_eM_M80::computeThValue()
{
    
    // double a_1400_bb_eP_0_eM_M80 = 67.4;// in percentage OLD
    // double a_1500_bb_eP_0_eM_M80 = 0.68746;
    double asym_SM = SM.getOptionalParameter("SM_a_1500_bb_eP_0_eM_M80");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 52.2209809;
    double NP_pos = 110.9946828*ceb + -8.3707272*ceQ + -0.0152505*cHb + 
    4.8653457*cHQ1 + 4.8653457*cHQ3 + 94.5528562*clb + 2355.0420720*clq1 + 2355.0420720*clq3;
    
    double SM_neg = 9.6720191;
    double NP_neg = 15.3433172*ceb + -55.3012728*ceQ + 0.9750505*cHb + 
    0.7824543*cHQ1 + 0.7824543*cHQ3 + 473.9471406*clb + 340.1779276*clq1 + 340.1779276*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}




sigma_1500_ttH_eP_0_eM_M80::sigma_1500_ttH_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1500_ttH_eP_0_eM_M80");
}

double sigma_1500_ttH_eP_0_eM_M80::computeThValue()
{
    
    //double sigma_1500_ttH_eP_0_eM_M80 = 1.549;//fb
    double sigma_1500_ttH_eP_0_eM_M80 = SM.getOptionalParameter("SM_sigma_1500_ttH_eP_0_eM_M80");
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    //double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phil1=0;  //For the moment we neglect these operators
    double C_phi3l1=0; //For the moment we neglect these operators
    double C_phi3l2=0; //For the moment we neglect these operators
    double C_phie=0;   //For the moment we neglect these operators

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+441.382932*C_phil1+5209.436564*C_phil1*C_phil1+423.232085*C_phi3l1+5170.240479*C_phi3l1*C_phi3l1-18.191681*C_phi3l2+0.878199*C_phi3l2*C_phi3l2-36.270805*C_phie+560.705729*C_phie*C_phie
                        -11.347306*C_tphi+0.36017*C_tphi*C_tphi-19.05008*C_phiQm+1.763122*C_phiQm*C_phiQm-9.612348*C_phit+1.759418*C_phit*C_phit-745.590684*C_tZ+2748.781721*C_tZ*C_tZ+1227.882916*C_tW+6530.341638*C_tW*C_tW-2810.538297*C_lqM+23675.21407*C_lqM*C_lqM-106.213862*C_eq+2568.916762*C_eq*C_eq-1618.490497*C_lu+23329.80876*C_lu*C_lu-218.873489*C_eu+2606.818362*C_eu*C_eu)*sigma_1500_ttH_eP_0_eM_M80/100;

            }
            else{
                return  (100+441.382932*C_phil1+423.232085*C_phi3l1-18.191681*C_phi3l2-36.270805*C_phie
                        -11.347306*C_tphi-19.05008*C_phiQm-9.612348*C_phit-745.590684*C_tZ+1227.882916*C_tW-2810.538297*C_lqM-106.213862*C_eq-1618.490497*C_lu-218.873489*C_eu)*sigma_1500_ttH_eP_0_eM_M80/100;
            }
    }
    
    ///////Needs to be rewritten in the old basis. NOT finished!!!
    else{
            if(flag_Quadratic){
                return  (100+441.382932*C_phil1+423.232085*C_phi3l1-18.191681*C_phi3l2-36.270805*C_phie
                        -11.347306*C_tphi-19.05008*C_phiQm-9.612348*C_phit-745.590684*C_tZ+1227.882916*C_tW-2810.538297*C_lqM-106.213862*C_eq-1618.490497*C_lu-218.873489*C_eu)*sigma_1500_ttH_eP_0_eM_M80/100;
            }
            else{
                return  (100+441.382932*C_phil1+423.232085*C_phi3l1-18.191681*C_phi3l2-36.270805*C_phie
                        -11.347306*C_tphi-19.05008*C_phiQm-9.612348*C_phit-745.590684*C_tZ+1227.882916*C_tW-2810.538297*C_lqM-106.213862*C_eq-1618.490497*C_lu-218.873489*C_eu)*sigma_1500_ttH_eP_0_eM_M80/100;
            }
    }
}



sigma_1500_ttH_eP_0_eM_P80::sigma_1500_ttH_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1500_ttH_eP_0_eM_P80");
}

double sigma_1500_ttH_eP_0_eM_P80::computeThValue()
{
    
    //double sigma_1500_ttH_eP_0_eM_P80 = 0.8954;//fb
    double sigma_1500_ttH_eP_0_eM_P80 = SM.getOptionalParameter("SM_sigma_1500_ttH_eP_0_eM_P80");
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    //double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phil1=0;  //For the moment we neglect these operators
    double C_phi3l1=0; //For the moment we neglect these operators
    double C_phi3l2=0; //For the moment we neglect these operators
    double C_phie=0;   //For the moment we neglect these operators

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+82.353422*C_phil1+959.488904*C_phil1*C_phil1+64.100178*C_phi3l1+952.358921*C_phi3l1*C_phi3l1-18.1692*C_phi3l2+0.784647*C_phi3l2*C_phi3l2-593.957784*C_phie+8875.780908*C_phie*C_phie
                        -11.047302*C_tphi+0.321624*C_tphi*C_tphi+5.914107*C_phiQm+2.41172*C_phiQm*C_phiQm+19.146706*C_phit+2.438465*C_phit*C_phit-1634.017986*C_tZ+13631.75009*C_tZ*C_tZ+1536.8943*C_tW+11228.0096*C_tW*C_tW-537.673973*C_lqM+4525.102152*C_lqM*C_lqM-1666.140519*C_eq+40471.80881*C_eq*C_eq-307.143464*C_lu+4389.278435*C_lu*C_lu-3450.82006*C_eu+40863.55189*C_eu*C_eu)*sigma_1500_ttH_eP_0_eM_P80/100;

            }
            else{
                return  (100+82.353422*C_phil1+64.100178*C_phi3l1-18.1692*C_phi3l2-593.957784*C_phie
                        -11.047302*C_tphi+5.914107*C_phiQm+19.146706*C_phit-1634.017986*C_tZ+1536.8943*C_tW-537.673973*C_lqM-1666.140519*C_eq-307.143464*C_lu-3450.82006*C_eu)*sigma_1500_ttH_eP_0_eM_P80/100;
    
            }
    }
    ///////Needs to be rewritten in the old basis. NOT finished!!!
    else{
            if(flag_Quadratic){
                return  (100+82.353422*C_phil1+64.100178*C_phi3l1-18.1692*C_phi3l2-593.957784*C_phie
                        -11.047302*C_tphi+5.914107*C_phiQm+19.146706*C_phit-1634.017986*C_tZ+1536.8943*C_tW-537.673973*C_lqM-1666.140519*C_eq-307.143464*C_lu-3450.82006*C_eu)*sigma_1500_ttH_eP_0_eM_P80/100;
            }
            else{
                return  (100+82.353422*C_phil1+64.100178*C_phi3l1-18.1692*C_phi3l2-593.957784*C_phie
                        -11.047302*C_tphi+5.914107*C_phiQm+19.146706*C_phit-1634.017986*C_tZ+1536.8943*C_tW-537.673973*C_lqM-1666.140519*C_eq-307.143464*C_lu-3450.82006*C_eu)*sigma_1500_ttH_eP_0_eM_P80/100;
            }
    }
}





/////// Prospects of Linear Colders at 3000 GeV //////////////////////////////////////////////////////////



sigma_3000_bb_eP_0_eM_P80::sigma_3000_bb_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_3000_bb_eP_0_eM_P80");
}

double sigma_3000_bb_eP_0_eM_P80::computeThValue()
{
    
    //double sigma_3000_bb_eP_0_eM_P80 = 5.158;//fb
    double sigma_3000_bb_eP_0_eM_P80_Madgraph = 4.8944;//fb
    
    double sigma_3000_bb_eP_0_eM_P80 = SM.getOptionalParameter("SM_sigma_3000_bb_eP_0_eM_P80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_3000_bb_eP_0_eM_P80 +(1136.42*ceb + 125508.2066*ceb*ceb + -564.16*ceQ +
                125499.7478*ceQ*ceQ + -0.436*cHb + 0.0249*cHb*cHb + 0.388*cHQ1 + 0.0252*cHQ1*cHQ1 +
                0.388*cHQ3 + 0.0252*cHQ3*cHQ3 + 63.736*clb + 13944.4584*clb*clb + 298.824*clq1 +
                13945.3384*clq1*clq1 + 298.824*clq3 + 13945.3384*clq3*clq3
                )*sigma_3000_bb_eP_0_eM_P80/sigma_3000_bb_eP_0_eM_P80_Madgraph;
    }
    else{
        return  sigma_3000_bb_eP_0_eM_P80 +(1136.42*ceb + -564.16*ceQ + -0.436*cHb +
                0.388*cHQ1 + 0.388*cHQ3 + 63.736*clb + 298.824*clq1 + 298.824*clq3
                )*sigma_3000_bb_eP_0_eM_P80/sigma_3000_bb_eP_0_eM_P80_Madgraph;
    }
}




sigma_3000_bb_eP_0_eM_M80::sigma_3000_bb_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_3000_bb_eP_0_eM_M80");
}

double sigma_3000_bb_eP_0_eM_M80::computeThValue()
{
    
    //double sigma_3000_bb_eP_0_eM_M80 = 14.92;//fb
    double sigma_3000_bb_eP_0_eM_M80_Madgraph = 15.4199;//fb
    
    double sigma_3000_bb_eP_0_eM_M80 = SM.getOptionalParameter("SM_sigma_3000_bb_eP_0_eM_M80");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_3000_bb_eP_0_eM_M80 +(126.272*ceb + 13945.2636*ceb*ceb + -62.706*ceQ + 
                13944.3236*ceQ*ceQ + 0.2395*cHb + 0.0333*cHb*cHb + 1.4052*cHQ1 + 0.0336*cHQ1*cHQ1 +
                1.4052*cHQ3 + 0.0336*cHQ3*cHQ3 + 572.94*clb + 125500.3648*clb*clb + 2690.64*clq1 + 
                125507.2472*clq1*clq1 + 2690.64*clq3 + 125507.2472*clq3*clq3
                )*sigma_3000_bb_eP_0_eM_M80/sigma_3000_bb_eP_0_eM_M80_Madgraph;
    }
    else{
        return  sigma_3000_bb_eP_0_eM_M80 +(126.272*ceb + -62.706*ceQ + 0.2395*cHb + 
                1.4052*cHQ1 + 1.4052*cHQ3 + 572.94*clb + 2690.64*clq1 + 2690.64*clq3
                )*sigma_3000_bb_eP_0_eM_M80/sigma_3000_bb_eP_0_eM_M80_Madgraph;
    }
}




a_3000_bb_eP_0_eM_P80::a_3000_bb_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_3000_bb_eP_0_eM_P80");
}

double a_3000_bb_eP_0_eM_P80::computeThValue()
{
    
    //double a_3000_bb_eP_0_eM_P80 = 52.8; //in percentage OLD
    //double a_3000_bb_eP_0_eM_P80 = 0.53118;
    double asym_SM = SM.getOptionalParameter("SM_a_3000_bb_eP_0_eM_P80");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 3.7471037;
    double NP_pos = 931.9530409*ceb + 7.9854048*ceQ + -0.4033776*cHb + 0.1647738*cHQ1 +
    0.1647738*cHQ3 + 7.9010564*clb + 262.4199070*clq1 + 262.4199070*clq3;
    
    double SM_neg = 1.1472963;
    double NP_neg = 204.4669803*ceb + -572.1454049*ceQ + -0.0326344*cHb + 0.2232302*cHQ1 +
    0.2232302*cHQ3 + 55.8349439*clb + 36.4040942*clq1 + 36.4040942*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}




a_3000_bb_eP_0_eM_M80::a_3000_bb_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_3000_bb_eP_0_eM_M80");
}

double a_3000_bb_eP_0_eM_M80::computeThValue()
{
    
    // double a_3000_bb_eP_0_eM_M80 = 67.6;// in percentage
    // double a_3000_bb_eP_0_eM_M80 = 0.68582;
    double asym_SM = SM.getOptionalParameter("SM_a_3000_bb_eP_0_eM_M80");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 12.9975879;
    double NP_pos = 111.0489210*ceb + -14.4162828*ceQ + -0.0100958*cHb + 1.2080351*cHQ1 + 
    1.2080351*cHQ3 + 167.1403272*clb + 2294.6886074*clq1 + 2294.6886074*clq3;
    
    double SM_neg = 2.4223121;
    double NP_neg = 15.2230773*ceb + -48.2897121*ceQ + 0.2495958*cHb + 0.1971449*cHQ1 + 
    0.1971449*cHQ3 + 405.7996696*clb + 395.9513926*clq1 + 395.9513926*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}





sigma_3000_ttH_eP_0_eM_M80::sigma_3000_ttH_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_3000_ttH_eP_0_eM_M80");
}

double sigma_3000_ttH_eP_0_eM_M80::computeThValue()
{
    
    //double sigma_3000_ttH_eP_0_eM_M80 = 0.5169;//fb
    double sigma_3000_ttH_eP_0_eM_M80 = SM.getOptionalParameter("SM_sigma_3000_ttH_eP_0_eM_M80");
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phil1=0;  //For the moment we neglect these operators
    double C_phi3l1=0; //For the moment we neglect these operators
    double C_phi3l2=0; //For the moment we neglect these operators
    double C_phie=0;   //For the moment we neglect these operators

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  sigma_3000_ttH_eP_0_eM_M80+(2803.82053*C_phil1+150027.1441*C_phil1*C_phil1+2785.533785*C_phi3l1+149773.1056*C_phi3l1*C_phi3l1-18.185847*C_phi3l2+0.853409*C_phi3l2*C_phi3l2-262.301253*C_phie+15636.68674*C_phie*C_phie
                        -10.76737*C_tphi+0.311377*C_tphi*C_tphi+0.311377*C_tphi*C_tphi-23.259425*C_phiQm+4.94913*C_phiQm*C_phiQm-9.576023*C_phit+4.787416*C_phit*C_phit-1607.494378*C_tZ+20765.31896*C_tZ*C_tZ+2637.92154*C_tW+49845.60589*C_tW*C_tW-10973.34078*C_lqM+385229.2693*C_lqM*C_lqM-333.929239*C_eq+41266.51772*C_eq*C_eq-5292.725653*C_lu+381085.3286*C_lu*C_lu-872.258364*C_eu+42018.82314*C_eu*C_eu)*sigma_3000_ttH_eP_0_eM_M80/100;

            }
            else{
                return  sigma_3000_ttH_eP_0_eM_M80+(2803.82053*C_phil1+2785.533785*C_phi3l1-18.185847*C_phi3l2-262.301253*C_phie
                        -10.76737*C_tphi-23.259425*C_phiQm-9.576023*C_phit-1607.494378*C_tZ+2637.92154*C_tW-10973.34078*C_lqM-333.929239*C_eq-5292.725653*C_lu-872.258364*C_eu)*sigma_3000_ttH_eP_0_eM_M80/100;
    
            }
    }
    
    ///////Needs to be rewritten in the old basis. NOT finished!!!
    else{
            if(flag_Quadratic){
                return  sigma_3000_ttH_eP_0_eM_M80+(2803.82053*C_phil1+2785.533785*C_phi3l1-18.185847*C_phi3l2-262.301253*C_phie
                        -10.76737*C_tphi-23.259425*C_phiQm-9.576023*C_phit-1607.494378*C_tZ+2637.92154*C_tW-10973.34078*C_lqM-333.929239*C_eq-5292.725653*C_lu-872.258364*C_eu)*sigma_3000_ttH_eP_0_eM_M80/100;
            }
            else{
                return  sigma_3000_ttH_eP_0_eM_M80+(2803.82053*C_phil1+2785.533785*C_phi3l1-18.185847*C_phi3l2-262.301253*C_phie
                        -10.76737*C_tphi-23.259425*C_phiQm-9.576023*C_phit-1607.494378*C_tZ+2637.92154*C_tW-10973.34078*C_lqM-333.929239*C_eq-5292.725653*C_lu-872.258364*C_eu)*sigma_3000_ttH_eP_0_eM_M80/100;
        }
    }
}



sigma_3000_ttH_eP_0_eM_P80::sigma_3000_ttH_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_3000_ttH_eP_0_eM_P80");
}

double sigma_3000_ttH_eP_0_eM_P80::computeThValue()
{
    
    //double sigma_3000_ttH_eP_0_eM_P80 = 0.3110;//fb
    double sigma_3000_ttH_eP_0_eM_P80 = SM.getOptionalParameter("SM_sigma_3000_ttH_eP_0_eM_P80");
    //double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    //double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    //double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phil1=0;  //For the moment we neglect these operators
    double C_phi3l1=0; //For the moment we neglect these operators
    double C_phi3l2=0; //For the moment we neglect these operators
    double C_phie=0;   //For the moment we neglect these operators

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  sigma_3000_ttH_eP_0_eM_P80+(492.978081*C_phil1+25686.93298*C_phil1*C_phil1+473.716107*C_phi3l1+25643.5765*C_phi3l1*C_phi3l1-18.189704*C_phi3l2+0.826774*C_phi3l2*C_phi3l2-3971.57722*C_phie+247583.0468*C_phie*C_phie
                        -10.289448*C_tphi+0.301091*C_tphi*C_tphi+4.521097*C_phiQm+6.468113*C_phiQm*C_phiQm+22.881545*C_phit+6.541524*C_phit*C_phit-3355.068539*C_tZ+100237.9032*C_tZ*C_tZ+3162.843439*C_tW+83057.23853*C_tW*C_tW-2042.514331*C_lqM+71708.84251*C_lqM*C_lqM-4938.435131*C_eq+635487.3558*C_eq*C_eq-976.17578*C_lu+69628.19228*C_lu*C_lu-13174.05422*C_eu+638535.0073*C_eu*C_eu)*sigma_3000_ttH_eP_0_eM_P80/100;

            }
            else{
                return  sigma_3000_ttH_eP_0_eM_P80+(492.978081*C_phil1+473.716107*C_phi3l1-18.189704*C_phi3l2-3971.57722*C_phie
                        -10.289448*C_tphi+4.521097*C_phiQm+22.881545*C_phit-3355.068539*C_tZ+3162.843439*C_tW-2042.514331*C_lqM-4938.435131*C_eq-976.17578*C_lu-13174.05422*C_eu)*sigma_3000_ttH_eP_0_eM_P80/100;
            }
    }
    
    ///////Needs to be rewritten in the old basis. NOT finished!!!
    else{
            if(flag_Quadratic){
                return  sigma_3000_ttH_eP_0_eM_P80+(492.978081*C_phil1+473.716107*C_phi3l1-18.189704*C_phi3l2-3971.57722*C_phie
                        -10.289448*C_tphi+4.521097*C_phiQm+22.881545*C_phit-3355.068539*C_tZ+3162.843439*C_tW-2042.514331*C_lqM-4938.435131*C_eq-976.17578*C_lu-13174.05422*C_eu)*sigma_3000_ttH_eP_0_eM_P80/100;
            }
            else{
                return  sigma_3000_ttH_eP_0_eM_P80+(492.978081*C_phil1+473.716107*C_phi3l1-18.189704*C_phi3l2-3971.57722*C_phie
                        -10.289448*C_tphi+4.521097*C_phiQm+22.881545*C_phit-3355.068539*C_tZ+3162.843439*C_tW-2042.514331*C_lqM-4938.435131*C_eq-976.17578*C_lu-13174.05422*C_eu)*sigma_3000_ttH_eP_0_eM_P80/100;
            }
    }
}



/////// Prospects of Linear Colders at 240 GeV //////////////////////////////////////////////////////////


sigma_240_bb::sigma_240_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_240_bb" );
}

double sigma_240_bb::computeThValue()
{
    
    //double sigma_240_bb_Madgraph = 1921;//fb OLD
    double sigma_240_bb_Madgraph = 1948.64;//fb
    double sigma_240_bb = SM.getOptionalParameter("SM_sigma_240_bb");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_240_bb +(654.3*ceb + 445.4682*ceb*ceb + -450.448*ceQ + 445.1694*ceQ*ceQ + 
                -22.284*cHb + 6.1388*cHb*cHb + 189.356*cHQ1 + 6.0824*cHQ1*cHQ1 + 189.356*cHQ3 + 
                6.0824*cHQ3*cHQ3 + 287.594*clb + 445.3588*clb*clb + 1661.668*clq1 + 
                445.4776*clq1*clq1 + 1661.668*clq3 + 445.4776*clq3*clq3
                )*sigma_240_bb/sigma_240_bb_Madgraph;
    }
    else{
        return  sigma_240_bb +(654.3*ceb + -450.448*ceQ + -22.284*cHb + 189.356*cHQ1 + 
                189.356*cHQ3 + 287.594*clb + 1661.668*clq1 + 1661.668*clq3
                )*sigma_240_bb/sigma_240_bb_Madgraph;
    }
}





a_240_bb::a_240_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_240_bb" );
}

double a_240_bb::computeThValue()
{
    
    //double a_240_bb = 58.4;//in percentage OLD
    // double a_240_bb = 0.62636;
    double asym_SM = SM.getOptionalParameter("SM_a_240_bb");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 1584.5950752;
    double NP_pos = 573.7625397*ceb + -55.6448403*ceQ + -39.0855201*cHb + 141.2653044*cHQ1 +
    141.2653044*cHQ3 + 34.2530098*clb + 1452.6897692*clq1 + 1452.6897692*clq3;
    
    double SM_neg = 364.0449248;
    double NP_neg = 80.5374604*ceb + -394.8031597*ceQ + 16.8015205*cHb + 48.0906956*cHQ1 + 
    48.0906956*cHQ3 + 253.3409900*clb + 208.9782308*clq1 + 208.9782308*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}




/////// Prospects of Linear Colders at 360 GeV //////////////////////////////////////////////////////////


sigma_360_bb::sigma_360_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_360_bb" );
}

double sigma_360_bb::computeThValue()
{
    
    //double sigma_360_bb_Madgraph = 757.8;//fb
    double sigma_360_bb_Madgraph = 767.33;//fb
    double sigma_360_bb = SM.getOptionalParameter("SM_sigma_360_bb");
    
    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    
    if(flag_Quadratic){
        return  sigma_360_bb +(640.698*ceb + 1003.5376*ceb*ceb + -370.464*ceQ + 
                1003.1506*ceQ*ceQ + -8.02*cHb + 2.2212*cHb*cHb + 70.676*cHQ1 + 
                2.2282*cHQ1*cHQ1 + 70.676*cHQ3 + 2.2282*cHQ3*cHQ3 + 304.324*clb + 
                1003.1153*clb*clb + 1561.9484*clq1 + 1003.4786*clq1*clq1 + 
                1561.9484*clq3 + 1003.4786*clq3*clq3
                )*sigma_360_bb/sigma_360_bb_Madgraph;
    }
    else{
        return  sigma_360_bb +(640.698*ceb + -370.464*ceQ + -8.02*cHb + 70.676*cHQ1 + 
                70.676*cHQ3 + 304.324*clb + 1561.9484*clq1 + 1561.9484*clq3
                )*sigma_360_bb/sigma_360_bb_Madgraph;
    }
}




a_360_bb::a_360_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_360_bb" );
}

double a_360_bb::computeThValue()
{
    
    // double a_360_bb = 62.9;// in percentage
    // double a_360_bb = 0.64086;
    double asym_SM = SM.getOptionalParameter("SM_a_360_bb");

    double cHQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double cHb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double ceb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double ceQ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double clb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double clqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    double clqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double clq1 = (clqP + clqM)/2.0;
    double clq3 = (clqP - clqM)/2.0;
    
    

    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHQ1;
    if(flag_LHC_WG_Basis){
        double cHQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        cHQ1 = cHQ3 + cHQm;
    }
    else{
        cHQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    }
    
    
    //We provide these values at linear order
    double SM_pos = 629.5405519;
    double NP_pos = 559.5420937*ceb + -46.1775177*ceQ + -15.0147661*cHb + 54.0983184*cHQ1 +
    54.0983184*cHQ3 + 36.8141818*clb + 1367.8373494*clq1 + 1367.8373494*clq3;
    
    double SM_neg = 137.7894481;
    double NP_neg = 81.1559063*ceb + -324.2864823*ceQ + 6.9947662*cHb + 16.5776817*cHQ1 + 
    16.5776817*cHQ3 + 267.5098183*clb + 194.1110506*clq1 + 194.1110506*clq3;
    
    double NP_num  = NP_pos-NP_neg;
    double NP_deno = NP_pos+NP_pos;
    
    //double SM_num  = SM_pos-SM_neg;
    double SM_deno = SM_pos+SM_neg;
    
    
    return asym_SM + (NP_num - asym_SM*NP_deno)/SM_deno;
}


/////// Prospects of Linear Colders at the Z pole  91.2 GeV //////////////////////////////////////////////////////////

//////// Cross Section ///////////////


sigma_Z_pole_bb::sigma_Z_pole_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_Z_pole_bb" );
}

double sigma_Z_pole_bb::computeThValue()
{
    
    
    //double sigma_Z_pole_bb_Madgraph = 8659000;//fb
    double sigma_Z_pole_bb = SM.getOptionalParameter("SM_sigma_Z_pole_bb");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  sigma_Z_pole_bb/100*(100+14*((C_phiQm+C_phiQ3)/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));

            }
            else{
                return  sigma_Z_pole_bb/100*(100+14*((C_phiQm+C_phiQ3)/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));
            }
    }
    
    else{
            if(flag_Quadratic){
                return  sigma_Z_pole_bb/100*(100+14*(C_phiQ1/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));
            }
            else{
                return  sigma_Z_pole_bb/100*(100+14*(C_phiQ1/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));
            }
    }
}



sigma_Z_pole_bb_eP_0_eM_P80::sigma_Z_pole_bb_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_Z_pole_bb_eP_0_eM_P80" );
}

double sigma_Z_pole_bb_eP_0_eM_P80::computeThValue()
{
    
    //double sigma_Z_pole_bb_eP_0_eM_P80_Madgraph = 7762000;//fb
    double sigma_Z_pole_bb_eP_0_eM_P80 = SM.getOptionalParameter("SM_sigma_Z_pole_bb_eP_0_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  sigma_Z_pole_bb_eP_0_eM_P80/100*(100+14*((C_phiQm+C_phiQ3)/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));

            }
            else{
                return  sigma_Z_pole_bb_eP_0_eM_P80/100*(100+14*((C_phiQm+C_phiQ3)/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));
            }
    }
    
    else{
            if(flag_Quadratic){
                return  sigma_Z_pole_bb_eP_0_eM_P80/100*(100+14*(C_phiQ1/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));
            }
            else{
                return  sigma_Z_pole_bb_eP_0_eM_P80/100*(100+14*(C_phiQ1/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));
            }
    }
}





sigma_Z_pole_bb_eP_0_eM_M80::sigma_Z_pole_bb_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_Z_pole_bb_eP_0_eM_M80" );
}

double sigma_Z_pole_bb_eP_0_eM_M80::computeThValue()
{
    
    //double sigma_Z_pole_bb_eP_0_eM_M80_Madgraph = 9554000;//fb
    double sigma_Z_pole_bb_eP_0_eM_M80 = SM.getOptionalParameter("SM_sigma_Z_pole_bb_eP_0_eM_M80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  sigma_Z_pole_bb_eP_0_eM_M80/100*(100+14*((C_phiQm+C_phiQ3)/0.998)+14*(C_phiQ3/0.998)-2.3*(C_phib/0.998));

            }
            else{
                return  sigma_Z_pole_bb_eP_0_eM_M80/100*(100+14*((C_phiQm+C_phiQ3)/0.998)+14*(C_phiQ3/0.998)-2.3*(C_phib/0.998));
            }
    }
    
    else{
            if(flag_Quadratic){
                return  sigma_Z_pole_bb_eP_0_eM_M80/100*(100+14*(C_phiQ1/0.998)+14*(C_phiQ3/0.998)-2.3*(C_phib/0.998));
            }
            else{
                return  sigma_Z_pole_bb_eP_0_eM_M80/100*(100+14*(C_phiQ1/0.998)+14*(C_phiQ3/0.998)-2.3*(C_phib/0.998));
            }
    }
}




sigma_Z_pole_bb_eP_M30_eM_P80::sigma_Z_pole_bb_eP_M30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_Z_pole_bb_eP_M30_eM_P80" );
}

double sigma_Z_pole_bb_eP_M30_eM_P80::computeThValue()
{
    
    //double sigma_Z_pole_bb_eP_M30_eM_P80_Madgraph = 9518000;//fb
    double sigma_Z_pole_bb_eP_M30_eM_P80 = SM.getOptionalParameter("SM_sigma_Z_pole_bb_eP_M30_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  sigma_Z_pole_bb_eP_M30_eM_P80/100*(100+14*((C_phiQm+C_phiQ3)/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));

            }
            else{
                return  sigma_Z_pole_bb_eP_M30_eM_P80/100*(100+14*((C_phiQm+C_phiQ3)/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));
            }
    }
    
    else{
            if(flag_Quadratic){
                return  sigma_Z_pole_bb_eP_M30_eM_P80/100*(100+14*(C_phiQ1/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));
            }
            else{
                return  sigma_Z_pole_bb_eP_M30_eM_P80/100*(100+14*(C_phiQ1/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));
            }
    }
}






sigma_Z_pole_bb_eP_P30_eM_M80::sigma_Z_pole_bb_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_Z_pole_bb_eP_P30_eM_M80" );
}

double sigma_Z_pole_bb_eP_P30_eM_M80::computeThValue()
{
    
    //double sigma_Z_pole_bb_eP_P30_eM_M80_Madgraph = 11990000;//fb
    double sigma_Z_pole_bb_eP_P30_eM_M80 = SM.getOptionalParameter("SM_sigma_Z_pole_bb_eP_P30_eM_M80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  sigma_Z_pole_bb_eP_P30_eM_M80/100*(100+14*((C_phiQm+C_phiQ3)/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));

            }
            else{
                return  sigma_Z_pole_bb_eP_P30_eM_M80/100*(100+14*((C_phiQm+C_phiQ3)/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));
            }
    }
    
    else{
            if(flag_Quadratic){
                return  sigma_Z_pole_bb_eP_P30_eM_M80/100*(100+14*(C_phiQ1/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));
            }
            else{
                return  sigma_Z_pole_bb_eP_P30_eM_M80/100*(100+14*(C_phiQ1/0.998)+14*(C_phiQ3/0.998)-2.5*(C_phib/0.998));
            }
    }
}







///////// Asymmetry //////////////////


a_Z_pole_bb::a_Z_pole_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_Z_pole_bb" );
}

double a_Z_pole_bb::computeThValue()
{
    
    double a_Z_pole_bb = SM.getOptionalParameter("SM_a_Z_pole_bb");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_Z_pole_bb + (0.1*((C_phiQm+C_phiQ3)/0.998)+0.1*(C_phiQ3/0.998)+0.5*(C_phib/0.998))/100;

            }
            else{
                return  a_Z_pole_bb + (0.1*((C_phiQm+C_phiQ3)/0.998)+0.1*(C_phiQ3/0.998)+0.5*(C_phib/0.998))/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_Z_pole_bb + (0.1*(C_phiQ1/0.998)+0.1*(C_phiQ3/0.998)+0.5*(C_phib/0.998))/100;
            }
            else{
                return  a_Z_pole_bb + (0.1*(C_phiQ1/0.998)+0.1*(C_phiQ3/0.998)+0.5*(C_phib/0.998))/100;
            }
    }
}




a_Z_pole_bb_eP_0_eM_P80::a_Z_pole_bb_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_Z_pole_bb_eP_0_eM_P80" );
}

double a_Z_pole_bb_eP_0_eM_P80::computeThValue()
{
    
    //double a_Z_pole_bb_eP_0_eM_P80 = -52.4;//in percentage
    double a_Z_pole_bb_eP_0_eM_P80 = SM.getOptionalParameter("SM_a_Z_pole_bb_eP_0_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    /////SOMETHING WEIRD HERE NEEDS TO BE CHECKED
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_Z_pole_bb_eP_0_eM_P80+(-0.5*((C_phiQm+C_phiQ3)/0.998)-0.5*(C_phiQ3/0.998)-3*(C_phib/0.998)-0.5*C_lqM)/100;

            }
            else{
                return  a_Z_pole_bb_eP_0_eM_P80+(-0.5*((C_phiQm+C_phiQ3)/0.998)-0.5*(C_phiQ3/0.998)-3*(C_phib/0.998)-0.5*C_lqM)/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_Z_pole_bb_eP_0_eM_P80+(-0.5*(C_phiQ1/0.998)-0.5*(C_phiQ3/0.998)-3*(C_phib/0.998)-0.5*C_lqM)/100;
            }
            else{
                return  a_Z_pole_bb_eP_0_eM_P80+(-0.5*(C_phiQ1/0.998)-0.5*(C_phiQ3/0.998)-3*(C_phib/0.998)-0.5*C_lqM)/100;
            }
    }
}


a_Z_pole_bb_eP_0_eM_M80::a_Z_pole_bb_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_Z_pole_bb_eP_0_eM_P80" );
}

double a_Z_pole_bb_eP_0_eM_M80::computeThValue()
{
    
    //double a_Z_pole_bb_eP_0_eM_M80 = 60;// in percentage
    double a_Z_pole_bb_eP_0_eM_M80 = SM.getOptionalParameter("SM_a_Z_pole_bb_eP_0_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_Z_pole_bb_eP_0_eM_M80 + (-0.5*((C_phiQm+C_phiQ3)/0.998)-0.5*(C_phiQ3/0.998)+2.8*(C_phib/0.998))/100;

            }
            else{
                return  a_Z_pole_bb_eP_0_eM_M80 + (-0.5*((C_phiQm+C_phiQ3)/0.998)-0.5*(C_phiQ3/0.998)+2.8*(C_phib/0.998))/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_Z_pole_bb_eP_0_eM_M80 + (-0.5*(C_phiQ1/0.998)-0.5*(C_phiQ3/0.998)+2.8*(C_phib/0.998))/100;
            }
            else{
                return  a_Z_pole_bb_eP_0_eM_M80 + (-0.5*(C_phiQ1/0.998)-0.5*(C_phiQ3/0.998)+2.8*(C_phib/0.998))/100;
            }
    }
}




a_Z_pole_bb_eP_M30_eM_P80::a_Z_pole_bb_eP_M30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_Z_pole_bb_eP_M30_eM_P80" );
}

double a_Z_pole_bb_eP_M30_eM_P80::computeThValue()
{
    
    // double a_Z_pole_bb_eP_M30_eM_P80 = -59.7;// in percentage
    double a_Z_pole_bb_eP_M30_eM_P80 = SM.getOptionalParameter("SM_a_Z_pole_bb_eP_M30_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_Z_pole_bb_eP_M30_eM_P80 + (-0.6*((C_phiQm+C_phiQ3)/0.998)-0.6*(C_phiQ3/0.998)-3.1*(C_phib/0.998))/100;

            }
            else{
                return  a_Z_pole_bb_eP_M30_eM_P80 + (-0.6*((C_phiQm+C_phiQ3)/0.998)-0.6*(C_phiQ3/0.998)-3.1*(C_phib/0.998))/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_Z_pole_bb_eP_M30_eM_P80 + (-0.6*(C_phiQ1/0.998)-0.6*(C_phiQ3/0.998)-3.1*(C_phib/0.998))/100;
            }
            else{
                return  a_Z_pole_bb_eP_M30_eM_P80 + (-0.6*(C_phiQ1/0.998)-0.6*(C_phiQ3/0.998)-3.1*(C_phib/0.998))/100;
            }
    }
}






a_Z_pole_bb_eP_P30_eM_M80::a_Z_pole_bb_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_Z_pole_bb_eP_P30_eM_M80" );
}

double a_Z_pole_bb_eP_P30_eM_M80::computeThValue()
{
    
    //double a_Z_pole_bb_eP_P30_eM_M80 = 64.5;// in percentage
    double a_Z_pole_bb_eP_P30_eM_M80 = SM.getOptionalParameter("SM_a_Z_pole_bb_eP_P30_eM_M80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    //double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_Z_pole_bb_eP_P30_eM_M80 + (0.6*((C_phiQm+C_phiQ3)/0.998)+0.6*(C_phiQ3/0.998)+3.2*(C_phib/0.998))/100;

            }
            else{
                return  a_Z_pole_bb_eP_P30_eM_M80 + (0.6*((C_phiQm+C_phiQ3)/0.998)+0.6*(C_phiQ3/0.998)+3.2*(C_phib/0.998))/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_Z_pole_bb_eP_P30_eM_M80 + (0.6*(C_phiQ1/0.998)+0.6*(C_phiQ3/0.998)+3.2*(C_phib/0.998))/100;
            }
            else{
                return  a_Z_pole_bb_eP_P30_eM_M80 + (0.6*(C_phiQ1/0.998)+0.6*(C_phiQ3/0.998)+3.2*(C_phib/0.998))/100;
            }
    }
}




///////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// OPTIMAL OBSERVABLES ////////////
//////////////////////////////////////////////////////////////////////////////////////////////




opt_obs_ilc_500_M30_M80::opt_obs_ilc_500_M30_M80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_ilc_500_M30_M80" << "eff_ilc_500_M30_M80");               
}

double opt_obs_ilc_500_M30_M80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_ilc_500_M30_M80");
    double eff = SM.getOptionalParameter("eff_ilc_500_M30_M80");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    


        
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_ilc_500_M30_M80_1 = 1.6388108355063156*C_phiQm + 0.6087828462366077*C_phit
 + -61.890724618292566*C_tW + 39.75505329432949*C_tZ
 + -4.5375385353076815e-07*CI_tW + 3.358614712277795e-07*CI_tZ
 + 94.4094725375339*C_lqM + 5.68112958507642*C_eq
 + 44.714498905389284*C_lu + 13.081461570106693*C_eu;

        return lumi*eff*opt_obs_ilc_500_M30_M80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_ilc_500_M30_M80_2 = 0.6087828462366078*C_phiQm + 1.2003646209917191*C_phit
 + -49.54909693163696*C_tW + 31.517839093491197*C_tZ
 + 4.5375385353076815e-07*CI_tW + -3.3586147107303444e-07*CI_tZ
 + 40.131540147129115*C_lqM + 7.859649201693296*C_eq
 + 70.14369540821401*C_lu + 5.295931344571681*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_M80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_ilc_500_M30_M80_3 = -61.89072461829256*C_phiQm + -49.54909693163696*C_phit
 + 3149.8689674659568*C_tW + -2023.5923066087373*C_tZ
 + 0.0*CI_tW + 1.215405407765529e-15*CI_tZ
 + -3722.6183854350847*C_lqM + -393.1825271022163*C_eq
 + -3176.6897996411835*C_lu + -538.5779626257699*C_eu;

        return lumi*eff*opt_obs_ilc_500_M30_M80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_ilc_500_M30_M80_4 = 39.755053294329485*C_phiQm + 31.51783909349119*C_phit
 + -2023.5923066087369*C_tW + 1305.854298764936*C_tZ
 + 0.0*CI_tW + -1.0718680202974305e-15*CI_tZ
 + 2394.940694108159*C_lqM + 256.8258832653644*C_eq
 + 2034.7150265358978*C_lu + 358.58651289404975*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_M80_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_ilc_500_M30_M80_5 = -4.5375385353076815e-07*C_phiQm + 4.5375385353076815e-07*C_phit
 + 0.0*C_tW + 0.0*C_tZ
 + 63.98323788950894*CI_tW + -41.20002299823266*CI_tZ
 + 4.571586963749346e-06*C_lqM + -5.776155790041091e-07*C_eq
 + -4.571586993537247e-06*C_lu + 5.776155818381808e-07*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_M80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_ilc_500_M30_M80_6 = 3.358614712277794e-07*C_phiQm + -3.3586147107303444e-07*C_phit
 + 1.215405407765529e-15*C_tW + -1.0718680202974305e-15*C_tZ
 + -41.20002299823267*CI_tW + 26.535300137348194*CI_tZ
 + 0.0001276573399437371*C_lqM + -1.0890583782044935e-06*C_eq
 + -0.00012765733992174254*C_lu + 1.0890583764385416e-06*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_M80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        double opt_obs_ilc_500_M30_M80_7 = 94.40947253753392*C_phiQm + 40.131540147129115*C_phit
 + -3722.6183854350847*C_tW + 2394.940694108159*C_tZ
 + 4.571586963749346e-06*CI_tW + 0.0001276573399437371*CI_tZ
 + 5478.8820315893045*C_lqM + 372.2141194997444*C_eq
 + 2861.7256424690518*C_lu + 765.2784244327142*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_M80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_ilc_500_M30_M80_8 = 5.68112958507642*C_phiQm + 7.859649201693296*C_phit
 + -393.1825271022163*C_tW + 256.8258832653644*C_tZ
 + -5.776155790041091e-07*CI_tW + -1.0890583782044931e-06*CI_tZ
 + 372.2141194997444*C_lqM + 70.84796056092617*C_eq
 + 480.31785699418907*C_lu + 58.56415794507859*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_M80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_ilc_500_M30_M80_9 = 44.71449890538929*C_phiQm + 70.14369540821401*C_phit
 + -3176.6897996411835*C_tW + 2034.7150265358976*C_tZ
 + -4.571586993537247e-06*CI_tW + -0.00012765733992174254*CI_tZ
 + 2861.7256424690518*C_lqM + 480.3178569941891*C_eq
 + 4192.860190165782*C_lu + 415.8196899853501*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_M80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_ilc_500_M30_M80_10 = 13.081461570106692*C_phiQm + 5.295931344571679*C_phit
 + -538.57796262577*C_tW + 358.58651289404986*C_tZ
 + 5.776155818381807e-07*CI_tW + 1.0890583764385414e-06*CI_tZ
 + 765.278424432714*C_lqM + 58.56415794507859*C_eq
 + 415.81968998535007*C_lu + 144.48630803103853*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_M80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_ilc_500_M30_M80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}












opt_obs_ilc_500_M30_P80::opt_obs_ilc_500_M30_P80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_ilc_500_M30_P80" << "eff_ilc_500_M30_P80");               
}

double opt_obs_ilc_500_M30_P80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_ilc_500_M30_P80");
    double eff = SM.getOptionalParameter("eff_ilc_500_M30_P80");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    


        
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_ilc_500_M30_P80_1 = 1.473280350963868*C_phiQm + 0.8558669789576177*C_phit
 + 71.14845273910511*C_tW + -78.32569012092065*C_tZ
 + 1.86496195078167e-05*CI_tW + -1.60363466517847e-05*CI_tZ
 + -5.555199051108081*C_lqM + -97.94963336923034*C_eq
 + -9.22804394132109*C_lu + -63.73711146290252*C_eu;

        return lumi*eff*opt_obs_ilc_500_M30_P80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_ilc_500_M30_P80_2 = 0.8558669789576178*C_phiQm + 2.3707804897443365*C_phit
 + 94.18390707604702*C_tW + -104.00278568083235*C_tZ
 + -1.86496195078167e-05*CI_tW + 1.603634665139295e-05*CI_tZ
 + -12.489635867217741*C_lqM + -67.4535041216653*C_eq
 + -6.177877354440506*C_lu + -154.4707207676446*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_P80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_ilc_500_M30_P80_3 = 71.14845273910508*C_phiQm + 94.18390707604702*C_phit
 + 5046.90680308217*C_tW + -5558.066430368364*C_tZ
 + 7.084551415757191e-15*CI_tW + -6.247877500120491e-15*CI_tZ
 + -556.1363276217832*C_lqM + -5058.148987505395*C_eq
 + -493.40553825194127*C_lu + -6419.144020481242*C_eu;

        return lumi*eff*opt_obs_ilc_500_M30_P80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_ilc_500_M30_P80_4 = -78.32569012092067*C_phiQm + -104.00278568083236*C_phit
 + -5558.066430368364*C_tW + 6124.4215686141415*C_tZ
 + -6.247877500120492e-15*CI_tW + 5.510013403204266e-15*CI_tZ
 + 616.019841811461*C_lqM + 5572.704924711923*C_eq
 + 541.77237616826*C_lu + 7084.853041854161*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_P80_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_ilc_500_M30_P80_5 = 1.86496195078167e-05*C_phiQm + -1.86496195078167e-05*C_phit
 + 7.084551415757191e-15*C_tW + -6.247877500120491e-15*C_tZ
 + 101.9943518363094*CI_tW + -112.47559954705848*CI_tZ
 + 3.9692297755298845e-06*C_lqM + 0.0006027328443855061*C_eq
 + -3.816999477733891e-06*C_lu + -0.0006437145940170607*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_P80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_ilc_500_M30_P80_6 = -1.6036346651784704e-05*C_phiQm + 1.6036346651392955e-05*C_phit
 + -6.247877500120492e-15*C_tW + 5.510013403204266e-15*C_tZ
 + -112.4755995470585*CI_tW + 124.05277666148825*CI_tZ
 + 8.363574642724233e-05*C_lqM + -0.001393494180541751*C_eq
 + -8.537760071380305e-05*C_lu + 0.001429636052681893*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_P80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        double opt_obs_ilc_500_M30_P80_7 = -5.555199051108081*C_phiQm + -12.489635867217741*C_phit
 + -556.1363276217832*C_tW + 616.019841811461*C_tZ
 + 3.9692297755298845e-06*CI_tW + 8.363574642724234e-05*CI_tZ
 + 90.7541172654634*C_lqM + 448.8743339767417*C_eq
 + 46.61154327023644*C_lu + 829.7903864992192*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_P80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_ilc_500_M30_P80_8 = -97.94963336923034*C_phiQm + -67.45350412166529*C_phit
 + -5058.148987505395*C_tW + 5572.704924711923*C_tZ
 + 0.0006027328443855061*CI_tW + -0.001393494180541751*CI_tZ
 + 448.8743339767417*C_lqM + 6602.657549907333*C_eq
 + 626.973975788001*C_lu + 4909.197684000763*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_P80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_ilc_500_M30_P80_9 = -9.228043941321088*C_phiQm + -6.177877354440507*C_phit
 + -493.40553825194115*C_tW + 541.7723761682599*C_tZ
 + -3.816999477733891e-06*CI_tW + -8.537760071380305e-05*CI_tZ
 + 46.61154327023644*C_lqM + 626.9739757880011*C_eq
 + 77.97789163093746*C_lu + 473.03352794471476*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_P80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_ilc_500_M30_P80_10 = -63.73711146290251*C_phiQm + -154.4707207676446*C_phit
 + -6419.144020481242*C_tW + 7084.853041854161*C_tZ
 + -0.0006437145940170606*CI_tW + 0.0014296360526818933*CI_tZ
 + 829.7903864992193*C_lqM + 4909.197684000764*C_eq
 + 473.03352794471476*C_lu + 10144.997076407559*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_M30_P80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_ilc_500_M30_P80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}













opt_obs_ilc_500_P30_M80::opt_obs_ilc_500_P30_M80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_ilc_500_P30_M80" << "eff_ilc_500_P30_M80");               
}

double opt_obs_ilc_500_P30_M80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_ilc_500_P30_M80");
    double eff = SM.getOptionalParameter("eff_ilc_500_P30_M80");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_ilc_500_P30_M80_1 = 3.5186751022116605*C_phiQm + 1.6884988529110831*C_phit
 + -120.19178947404997*C_tW + 71.7175623113869*C_tZ
 + -1.9137900158989966e-05*CI_tW + 1.663558025918004e-05*CI_tZ
 + 195.31327707410614*C_lqM + 3.744578685564284*C_eq
 + 99.20925285715637*C_lu + 8.050235419947409*C_eu;

        return lumi*eff*opt_obs_ilc_500_P30_M80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_ilc_500_P30_M80_2 = 1.6884988529110831*C_phiQm + 2.622092726658831*C_phit
 + -100.70975857166331*C_tW + 59.950992746168325*C_tZ
 + -1.832711599033101e-05*CI_tW + 1.640487886154404e-05*CI_tZ
 + 96.42808477423281*C_lqM + 4.881424311553053*C_eq
 + 146.62569023884342*C_lu + 4.019152435312039*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_M80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_ilc_500_P30_M80_3 = -120.19178947404998*C_phiQm + -100.70975857166331*C_phit
 + 5216.288128089109*C_tW + -3105.404538038601*C_tZ
 + -0.0017452997076787726*CI_tW + 0.0015391826715125362*CI_tZ
 + -6739.022191248597*C_lqM + -204.87560749940198*C_eq
 + -5743.100146346434*C_lu + -281.5133376904457*C_eu;

        return lumi*eff*opt_obs_ilc_500_P30_M80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_ilc_500_P30_M80_4 = 71.71756231138689*C_phiQm + 59.950992746168325*C_phit
 + -3105.404538038601*C_tW + 1853.800393838675*C_tZ
 + 0.0015391826715104118*CI_tW + -0.0013574077196350964*CI_tZ
 + 4022.3459880439113*C_lqM + 123.62510656026484*C_eq
 + 3423.085539113187*C_lu + 172.5015816680891*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_M80_4;

        
    } else if(b_min == 5 && b_max == 5){
        

        double opt_obs_ilc_500_P30_M80_5 = -1.9137900158989966e-05*C_phiQm + -1.832711599033101e-05*C_phit
 + -0.0017452997076787726*C_tW + 0.0015391826715104118*C_tZ
 + 106.16313454567828*CI_tW + -63.34234176112494*CI_tZ
 + 0.0004128493513643155*C_lqM + 2.6561312859269752e-06*C_eq
 + -0.0004128493513715997*C_lu + -2.6655129889878537e-06*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_M80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_ilc_500_P30_M80_6 = 1.6635580259180036e-05*C_phiQm + 1.6404878861544037e-05*C_phit
 + 0.0015391826715125362*C_tW + -0.0013574077196350964*C_tZ
 + -63.34234176112493*CI_tW + 37.78737131087368*CI_tZ
 + -0.00030281986190683336*C_lqM + -1.864096884308424e-05*C_eq
 + 0.00030281986189902667*C_lu + 1.869372174121568e-05*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_M80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        double opt_obs_ilc_500_P30_M80_7 = 195.31327707410614*C_phiQm + 96.42808477423281*C_phit
 + -6739.022191248596*C_tW + 4022.3459880439113*C_tZ
 + 0.00041284935136431554*CI_tW + -0.0003028198619068334*CI_tZ
 + 10847.388031348375*C_lqM + 214.78702490809235*C_eq
 + 5655.585549512605*C_lu + 448.0790418520734*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_M80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_ilc_500_P30_M80_8 = 3.744578685564284*C_phiQm + 4.881424311553054*C_phit
 + -204.87560749940198*C_tW + 123.6251065602649*C_tZ
 + 2.6561312859269757e-06*CI_tW + -1.864096884308424e-05*CI_tZ
 + 214.78702490809235*C_lqM + 11.881968590943949*C_eq
 + 275.13336795660103*C_lu + 9.96477827438408*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_M80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_ilc_500_P30_M80_9 = 99.20925285715637*C_phiQm + 146.62569023884342*C_phit
 + -5743.100146346434*C_tW + 3423.085539113187*C_tZ
 + -0.0004128493513715997*CI_tW + 0.0003028198618990267*CI_tZ
 + 5655.585549512606*C_lqM + 275.13336795660103*C_eq
 + 8213.590978054763*C_lu + 241.40239924727*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_M80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_ilc_500_P30_M80_10 = 8.050235419947409*C_phiQm + 4.019152435312038*C_phit
 + -281.5133376904457*C_tW + 172.5015816680891*C_tZ
 + -2.6655129889878537e-06*CI_tW + 1.869372174121568e-05*CI_tZ
 + 448.0790418520734*C_lqM + 9.964778274384082*C_eq
 + 241.40239924727004*C_lu + 25.14956896249474*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_M80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_ilc_500_P30_M80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}






opt_obs_ilc_500_P30_P80::opt_obs_ilc_500_P30_P80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_ilc_500_P30_P80" << "eff_ilc_500_P30_P80");               
}

double opt_obs_ilc_500_P30_P80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_ilc_500_P30_P80");
    double eff = SM.getOptionalParameter("eff_ilc_500_P30_P80");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
        
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_ilc_500_P30_P80_1 = 0.4733668460101026*C_phiQm + 0.005756343626120958*C_phit
 + 22.762006896207158*C_tW + -23.17958375877933*C_tZ
 + -3.925645983740149e-06*CI_tW + 4.946022417778193e-06*CI_tZ
 + -2.3638086699094676*C_lqM + -32.13085739637971*C_eq
 + -9.774645391519162*C_lu + -11.495703966221601*C_eu;

        return lumi*eff*opt_obs_ilc_500_P30_P80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_ilc_500_P30_P80_2 = 0.0057563436261209855*C_phiQm + 0.8601413576896537*C_phit
 + 36.89957118648607*C_tW + -37.956826882338994*C_tZ
 + 4.0130915566293896e-06*CI_tW + -1.7311533488048075e-06*CI_tZ
 + -14.962958794240743*C_lqM + -17.408009480232202*C_eq
 + -4.337439418110696*C_lu + -58.43192804314783*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_P80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_ilc_500_P30_P80_3 = 22.762006896207158*C_phiQm + 36.89957118648607*C_phit
 + 3048.89891247043*C_tW + -3085.8057697355357*C_tZ
 + 0.0*CI_tW + -1.2882803255075677e-15*CI_tZ
 + -961.9281659113891*C_lqM + -2511.601189704324*C_eq
 + -842.9794543770362*C_lu + -3255.219637505407*C_eu;

        return lumi*eff*opt_obs_ilc_500_P30_P80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_ilc_500_P30_P80_4 = -23.17958375877933*C_phiQm + -37.956826882339*C_phit
 + -3085.8057697355357*C_tW + 3127.708096860304*C_tZ
 + -4.397918814979122e-05*CI_tW + 3.7457640266459524e-05*CI_tZ
 + 968.9409490551667*C_lqM + 2545.550545973207*C_eq
 + 842.1533257191585*C_lu + 3320.0348961155846*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_P80_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_ilc_500_P30_P80_5 = -3.92564598374015e-06*C_phiQm + 4.0130915566293896e-06*C_phit
 + 0.0*C_tW + -4.3979188149791225e-05*C_tZ
 + 61.425393657016826*CI_tW + -62.27247793767658*CI_tZ
 + 5.2319165417257535e-06*C_lqM + -5.82493810743151e-06*C_eq
 + -5.242656511477227e-06*C_lu + 5.82493810743151e-06*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_P80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_ilc_500_P30_P80_6 = 4.946022417778193e-06*C_phiQm + -1.7311533488048075e-06*C_phit
 + -1.2882803255075677e-15*C_tW + 3.745764026645953e-05*C_tZ
 + -62.27247793767658*CI_tW + 63.17102901190249*CI_tZ
 + -2.2809564712517208e-05*C_lqM + 5.2646698065230034e-05*C_eq
 + 2.281903630943624e-05*C_lu + -5.2646698065230034e-05*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_P80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        double opt_obs_ilc_500_P30_P80_7 = -2.363808669909469*C_phiQm + -14.962958794240743*C_phit
 + -961.928165911389*C_tW + 968.9409490551666*C_tZ
 + 5.2319165417257535e-06*CI_tW + -2.280956471251721e-05*CI_tZ
 + 421.3332324199011*C_lqM + 627.1073242906411*C_eq
 + 218.57294795931676*C_lu + 1179.5572262441685*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_P80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_ilc_500_P30_P80_8 = -32.13085739637971*C_phiQm + -17.408009480232195*C_phit
 + -2511.601189704324*C_tW + 2545.5505459732067*C_tZ
 + -5.82493810743151e-06*CI_tW + 5.2646698065230034e-05*CI_tZ
 + 627.107324290641*C_lqM + 2712.675104725252*C_eq
 + 856.9108996593474*C_lu + 2059.0064893056597*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_P80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_ilc_500_P30_P80_9 = -9.774645391519162*C_phiQm + -4.337439418110695*C_phit
 + -842.9794543770362*C_tW + 842.1533257191585*C_tZ
 + -5.242656511477227e-06*CI_tW + 2.281903630943624e-05*CI_tZ
 + 218.5729479593168*C_lqM + 856.9108996593474*C_eq
 + 350.4964024241973*C_lu + 669.1274417171611*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_P80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_ilc_500_P30_P80_10 = -11.495703966221598*C_phiQm + -58.43192804314782*C_phit
 + -3255.219637505407*C_tW + 3320.034896115584*C_tZ
 + 5.82493810743151e-06*CI_tW + -5.2646698065230034e-05*CI_tZ
 + 1179.5572262441688*C_lqM + 2059.006489305659*C_eq
 + 669.1274417171611*C_lu + 4396.112489563018*C_eu;
        
        return lumi*eff*opt_obs_ilc_500_P30_P80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_ilc_500_P30_P80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}






opt_obs_ilc_1000_M20_M80::opt_obs_ilc_1000_M20_M80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_ilc_1000_M20_M80" << "eff_ilc_1000_M20_M80");               
}

double opt_obs_ilc_1000_M20_M80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_ilc_1000_M20_P80");
    double eff = SM.getOptionalParameter("eff_ilc_1000_M20_P80");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
        
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_ilc_1000_M20_M80_1 = 0.6860255612551367*C_phiQm + 0.06288087439962653*C_phit
 + -23.41498140175173*C_tW + 14.942617472141432*C_tZ
 + -1.5850460242713856e-07*CI_tW + 6.819953588639327e-08*CI_tZ
 + 156.16560963066826*C_lqM + 2.902415391808151*C_eq
 + 29.946100561240897*C_lu + 18.079647292882644*C_eu;

        return lumi*eff*opt_obs_ilc_1000_M20_M80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_ilc_1000_M20_M80_2 = 0.06288087439962661*C_phiQm + 0.4592130249379536*C_phit
 + -16.747559833106205*C_tW + 10.546031168151007*C_tZ
 + 1.4171817091957683e-07*CI_tW + -1.0506971653627906e-07*CI_tZ
 + 21.413271565666292*C_lqM + 8.3541457325761*C_eq
 + 104.95434479416713*C_lu + 2.420719478427417*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_M80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_ilc_1000_M20_M80_3 = -23.414981401751735*C_phiQm + -16.747559833106205*C_phit
 + 1463.9138521156908*C_tW + -931.7561776504119*C_tZ
 + 0.0*CI_tW + -3.379419914274885e-15*CI_tZ
 + -5583.078766211686*C_lqM + -388.27174840632733*C_eq
 + -4357.824116447796*C_lu + -691.5044987503966*C_eu;

        return lumi*eff*opt_obs_ilc_1000_M20_M80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_ilc_1000_M20_M80_4 = 14.942617472141432*C_phiQm + 10.546031168151007*C_phit
 + -931.7561776504119*C_tW + 594.8517664330902*C_tZ
 + 4.663994735490485e-15*CI_tW + -1.482830675412827e-15*CI_tZ
 + 3564.5522743690162*C_lqM + 249.86136660430208*C_eq
 + 2764.720983522906*C_lu + 458.659745221885*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_M80_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_ilc_1000_M20_M80_5 = -1.5850460242713856e-07*C_phiQm + 1.4171817091957683e-07*C_phit
 + 0.0*C_tW + 4.663994735490485e-15*C_tZ
 + 182.80996244705983*CI_tW + -116.77066944091251*CI_tZ
 + -0.002957143824068356*C_lqM + 2.2117614580499612e-05*C_eq
 + 0.002957143824114569*C_lu + 2.3275677716787556e-05*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_M80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_ilc_1000_M20_M80_6 = 6.819953588639328e-08*C_phiQm + -1.0506971653627905e-07*C_phit
 + -3.379419914274885e-15*C_tW + -1.4828306754128266e-15*C_tZ
 + -116.77066944091251*CI_tW + 74.62198395131131*CI_tZ
 + -0.0001915373055074597*C_lqM + -1.2459929722146745e-05*C_eq
 + 0.00019162955762175153*C_lu + -1.4330068098814526e-05*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_M80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        double opt_obs_ilc_1000_M20_M80_7 = 156.16560963066826*C_phiQm + 21.413271565666307*C_phit
 + -5583.078766211686*C_tW + 3564.552274369016*C_tZ
 + -0.002957143824068356*CI_tW + -0.00019153730550745948*CI_tZ
 + 35714.55119597334*C_lqM + 849.3374780265794*C_eq
 + 8426.06923900336*C_lu + 4137.744662499108*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_M80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_ilc_1000_M20_M80_8 = 2.902415391808149*C_phiQm + 8.3541457325761*C_phit
 + -388.2717484063273*C_tW + 249.86136660430202*C_tZ
 + 2.211761458049961e-05*CI_tW + -1.2459929722146745e-05*CI_tZ
 + 849.3374780265793*C_lqM + 227.18844157301316*C_eq
 + 1960.0742786114297*C_lu + 101.7852406606287*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_M80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_ilc_1000_M20_M80_9 = 29.94610056124089*C_phiQm + 104.95434479416713*C_phit
 + -4357.824116447795*C_tW + 2764.720983522905*C_tZ
 + 0.0029571438241145693*CI_tW + 0.00019162955762175174*CI_tZ
 + 8426.069239003356*C_lqM + 1960.0742786114301*C_eq
 + 24522.349727420136*C_lu + 1163.0189013968397*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_M80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_ilc_1000_M20_M80_10 = 18.079647292882644*C_phiQm + 2.4207194784274177*C_phit
 + -691.5044987503966*C_tW + 458.659745221885*C_tZ
 + 2.3275677716787556e-05*CI_tW + -1.4330068098814524e-05*CI_tZ
 + 4137.744662499108*C_lqM + 101.78524066062869*C_eq
 + 1163.0189013968397*C_lu + 707.7057375048124*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_M80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_ilc_1000_M20_M80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}






opt_obs_ilc_1000_M20_P80::opt_obs_ilc_1000_M20_P80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_ilc_1000_M20_P80" << "eff_ilc_1000_M20_P80");               
}

double opt_obs_ilc_1000_M20_P80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_ilc_1000_M20_P80");
    double eff = SM.getOptionalParameter("eff_ilc_1000_M20_P80");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
        
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_ilc_1000_M20_P80_1 = 0.38815750859537757*C_phiQm + 0.017725299710632858*C_phit
 + 16.917847736363274*C_tW + -18.334092520201246*C_tZ
 + -5.236205865034183e-06*CI_tW + 4.4702132972494e-06*CI_tZ
 + -1.8604178070130581*C_lqM + -101.16393351009786*C_eq
 + -12.47548388369432*C_lu + -18.739068237217694*C_eu;

        return lumi*eff*opt_obs_ilc_1000_M20_P80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_ilc_1000_M20_P80_2 = 0.017725299710632803*C_phiQm + 0.7629996081368549*C_phit
 + 27.644347768486423*C_tW + -30.142537021454572*C_tZ
 + -5.508611194859263e-06*CI_tW + 4.748558754616871e-06*CI_tZ
 + -18.87698560535277*C_lqM + -26.03270566173262*C_eq
 + -3.079437486950912*C_lu + -198.19718964853132*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_P80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_ilc_1000_M20_P80_3 = 16.917847736363274*C_phiQm + 27.644347768486423*C_phit
 + 1955.2076088573008*C_tW + -2120.609705874259*C_tZ
 + -2.7751319746930538e-14*CI_tW + 2.447393434958441e-14*CI_tZ
 + -819.3657072312832*C_lqM + -5250.318815369343*C_eq
 + -695.5831807727862*C_lu + -7846.641091336725*C_eu;

        return lumi*eff*opt_obs_ilc_1000_M20_P80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_ilc_1000_M20_P80_4 = -18.33409252020125*C_phiQm + -30.142537021454576*C_phit
 + -2120.609705874259*C_tW + 2300.7201369929126*C_tZ
 + 1.810046696739932e-14*CI_tW + 2.0485689567875043e-07*CI_tZ
 + 896.1078376764397*C_lqM + 5699.131446787922*C_eq
 + 748.8847000776334*C_lu + 8544.776506684268*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_P80_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_ilc_1000_M20_P80_5 = -5.236205865034183e-06*C_phiQm + -5.508611194859263e-06*C_phit
 + -2.7751319746930538e-14*C_tW + 1.810046696739932e-14*C_tZ
 + 241.39470069314348*CI_tW + -262.4673637876004*CI_tZ
 + -1.9190774479269234e-06*C_lqM + 0.0022836962741094964*C_eq
 + 7.61776900034376e-05*C_lu + 0.0021438044514300618*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_P80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_ilc_1000_M20_P80_6 = 4.4702132972494e-06*C_phiQm + 4.748558754616871e-06*C_phit
 + 2.447393434958441e-14*C_tW + 2.0485689567875046e-07*C_tZ
 + -262.46736378760033*CI_tW + 285.4211030820254*CI_tZ
 + 0.00019157677869885328*C_lqM + -0.0018807003184547518*C_eq
 + 0.0003749395117035614*C_lu + -0.0014934101986907253*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_P80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        double opt_obs_ilc_1000_M20_P80_7 = -1.860417807013059*C_phiQm + -18.876985605352772*C_phit
 + -819.3657072312831*C_tW + 896.1078376764397*C_tZ
 + -1.9190774479269217e-06*CI_tW + 0.0001915767786988533*CI_tZ
 + 700.5228232572056*C_lqM + 1272.9465621157246*C_eq
 + 154.7348659814437*C_lu + 4993.080386550328*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_P80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_ilc_1000_M20_P80_8 = -101.16393351009786*C_phiQm + -26.03270566173262*C_phit
 + -5250.318815369343*C_tW + 5699.131446787922*C_tZ
 + 0.0022836962741094964*CI_tW + -0.0018807003184547516*CI_tZ
 + 1272.9465621157246*C_lqM + 27263.454945510202*C_eq
 + 3359.6641265001895*C_lu + 10471.605923028983*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_P80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_ilc_1000_M20_P80_9 = -12.47548388369432*C_phiQm + -3.079437486950913*C_phit
 + -695.5831807727862*C_tW + 748.8847000776334*C_tZ
 + 7.617769000343759e-05*CI_tW + 0.0003749395117035614*CI_tZ
 + 154.73486598144376*C_lqM + 3359.6641265001895*C_eq
 + 603.2690567268976*C_lu + 1473.2800933050007*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_P80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_ilc_1000_M20_P80_10 = -18.73906823721768*C_phiQm + -198.19718964853132*C_phit
 + -7846.641091336726*C_tW + 8544.776506684268*C_tZ
 + 0.0021438044514300618*CI_tW + -0.0014934101986907253*CI_tZ
 + 4993.080386550328*C_lqM + 10471.605923028983*C_eq
 + 1473.2800933050007*C_lu + 52250.58196060555*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M20_P80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_ilc_1000_M20_P80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}





opt_obs_ilc_1000_P20_M80::opt_obs_ilc_1000_P20_M80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_ilc_1000_P20_M80" << "eff_ilc_1000_P20_M80");               
}

double opt_obs_ilc_1000_P20_M80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_ilc_1000_P20_M80");
    double eff = SM.getOptionalParameter("eff_ilc_1000_P20_M80");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
        
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_ilc_1000_P20_M80_1 = 1.09540877542705*C_phiQm + 0.18862148058577427*C_phit
 + -35.860481188499*C_tW + 21.837694350089254*C_tZ
 + -3.8556927569159297e-07*CI_tW + -4.208355006267873e-05*CI_tZ
 + 247.35940816466837*C_lqM + 2.361489532060558*C_eq
 + 53.608109714112345*C_lu + 12.958176063415348*C_eu;

        return lumi*eff*opt_obs_ilc_1000_P20_M80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_ilc_1000_P20_M80_2 = 0.18862148058577427*C_phiQm + 0.7311115129696788*C_phit
 + -26.932680566372618*C_tW + 16.32311969002878*C_tZ
 + 2.4581259926028614e-07*CI_tW + -4.248116619643736e-05*CI_tZ
 + 47.431814339684934*C_lqM + 5.919016723447084*C_eq
 + 166.141524013229*C_lu + 2.7626531851219776*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_M80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_ilc_1000_P20_M80_3 = -35.860481188499*C_phiQm + -26.932680566372618*C_phit
 + 2044.5059429291803*C_tW + -1241.3701949747465*C_tZ
 + 0.0*CI_tW + 0.0*CI_tZ
 + -8252.314633901886*C_lqM + -253.2486053570436*C_eq
 + -6429.266922940911*C_lu + -453.90775475959526*C_eu;

        return lumi*eff*opt_obs_ilc_1000_P20_M80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_ilc_1000_P20_M80_4 = 21.837694350089254*C_phiQm + 16.32311969002878*C_phit
 + -1241.3701949747465*C_tW + 755.3703536746532*C_tZ
 + 0.0*CI_tW + 1.1273319498160023e-15*CI_tZ
 + 5026.268211757163*C_lqM + 155.28617443744682*C_eq
 + 3905.8535280221327*C_lu + 285.6331744755713*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_M80_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_ilc_1000_P20_M80_5 = -3.8556927569159297e-07*C_phiQm + 2.4581259926028614e-07*C_phit
 + 0.0*C_tW + 0.0*C_tZ
 + 255.6325928253371*CI_tW + -155.71018707273277*CI_tZ
 + -1.0234999869608925e-05*C_lqM + -0.0007765274107141076*C_eq
 + 1.0234999958369636e-05*C_lu + -0.0007757483322086861*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_M80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_ilc_1000_P20_M80_6 = -4.208355006267874e-05*C_phiQm + -4.2481166196437365e-05*C_phit
 + 0.0*C_tW + 1.1273319498160023e-15*C_tZ
 + -155.71018707273274*CI_tW + 94.85691448041831*CI_tZ
 + 8.383098582487106e-06*C_lqM + 0.0014309239557751135*C_eq
 + -7.905455998246563e-06*C_lu + 0.001429619036931999*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_M80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        double opt_obs_ilc_1000_P20_M80_7 = 247.3594081646683*C_phiQm + 47.43181433968492*C_phit
 + -8252.314633901886*C_tW + 5026.268211757162*C_tZ
 + -1.0234999869608927e-05*CI_tW + 8.383098582487104e-06*CI_tZ
 + 55908.33771830471*C_lqM + 590.3421847455065*C_eq
 + 13192.129236062276*C_lu + 2929.598534938029*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_M80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_ilc_1000_P20_M80_8 = 2.3614895320605584*C_phiQm + 5.9190167234470845*C_phit
 + -253.24860535704366*C_tW + 155.28617443744682*C_tZ
 + -0.0007765274107141076*CI_tW + 0.0014309239557751135*CI_tZ
 + 590.3421847455065*C_lqM + 70.13333227786778*C_eq
 + 1353.3567213303177*C_lu + 31.820356432446374*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_M80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_ilc_1000_P20_M80_9 = 53.608109714112324*C_phiQm + 166.141524013229*C_phit
 + -6429.266922940911*C_tW + 3905.853528022133*C_tZ
 + 1.0234999958369636e-05*CI_tW + -7.905455998246563e-06*CI_tZ
 + 13192.12923606227*C_lqM + 1353.3567213303177*C_eq
 + 37920.70357525274*C_lu + 818.0937122149071*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_M80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_ilc_1000_P20_M80_10 = 12.958176063415348*C_phiQm + 2.7626531851219776*C_phit
 + -453.9077547595953*C_tW + 285.63317447557125*C_tZ
 + -0.0007757483322086861*CI_tW + 0.0014296190369319987*CI_tZ
 + 2929.598534938029*C_lqM + 31.820356432446374*C_eq
 + 818.0937122149071*C_lu + 227.3111085740191*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_M80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_ilc_1000_P20_M80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}







opt_obs_ilc_1000_P20_P80::opt_obs_ilc_1000_P20_P80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_ilc_1000_P20_P80" << "eff_ilc_1000_P20_P80");               
}

double opt_obs_ilc_1000_P20_P80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_ilc_1000_P20_P80");
    double eff = SM.getOptionalParameter("eff_ilc_1000_P20_P80");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
        
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_ilc_1000_P20_P80_1 = 0.2133873547953592*C_phiQm + -0.08204510089055535*C_phit
 + 7.609915207933011*C_tW + -7.823409130970538*C_tZ
 + 4.73627890037032e-05*CI_tW + -4.1803568806235226e-05*CI_tZ
 + 2.3473006572698276*C_lqM + -51.76068979603156*C_eq
 + -13.874493035865026*C_lu + 5.112067117457322*C_eu;

        return lumi*eff*opt_obs_ilc_1000_P20_P80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_ilc_1000_P20_P80_2 = -0.08204510089055536*C_phiQm + 0.4362606739526286*C_phit
 + 15.456186935555392*C_tW + -16.12746803190067*C_tZ
 + 3.9251341347755306e-05*CI_tW + -3.4610178296704273e-05*CI_tZ
 + -22.833811820556413*C_lqM + -5.071716399604718*C_eq
 + -0.3002445465021264*C_lu + -111.67142836604901*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_P80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_ilc_1000_P20_P80_3 = 7.609915207933011*C_phiQm + 15.456186935555394*C_phit
 + 1400.4489833258497*C_tW + -1439.379070377874*C_tZ
 + 1.7858100018417536e-14*CI_tW + -1.7115801162713176e-14*CI_tZ
 + -1172.8478580194842*C_lqM + -3278.2485728256843*C_eq
 + -976.7877661608331*C_lu + -5056.925263873175*C_eu;

        return lumi*eff*opt_obs_ilc_1000_P20_P80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_ilc_1000_P20_P80_4 = -7.823409130970539*C_phiQm + -16.12746803190067*C_phit
 + -1439.379070377874*C_tW + 1480.4727992282471*C_tZ
 + -1.8170017133778093e-14*CI_tW + 1.4306312182279883e-06*CI_tZ
 + 1206.3380574246942*C_lqM + 3370.8570331613587*C_eq
 + 989.4288463610189*C_lu + 5242.634125639532*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_P80_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_ilc_1000_P20_P80_5 = 4.7362789003703204e-05*C_phiQm + 3.925134134775531e-05*C_phit
 + 1.7858100018417536e-14*C_tW + -1.8170017133778096e-14*C_tZ
 + 172.85250585026296*CI_tW + -178.17406927719543*CI_tZ
 + -1.378805527384992e-05*C_lqM + -0.0007729415946775807*C_eq
 + -3.858758917118177e-06*C_lu + 0.00028074072062016453*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_P80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_ilc_1000_P20_P80_6 = -4.180356880623523e-05*C_phiQm + -3.4610178296704273e-05*C_phit
 + -1.7115801162713176e-14*C_tW + 1.4306312182279883e-06*C_tZ
 + -178.17406927719543*CI_tW + 183.75237985496273*CI_tZ
 + 3.359556558028491e-06*C_lqM + 0.0006354777989346247*C_eq
 + -9.060090807748995e-05*C_lu + -0.0002014050417294052*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_P80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        double opt_obs_ilc_1000_P20_P80_7 = 2.347300657269825*C_phiQm + -22.833811820556413*C_phit
 + -1172.8478580194842*C_tW + 1206.3380574246942*C_tZ
 + -1.3788055273849919e-05*CI_tW + 3.359556558028491e-06*CI_tZ
 + 1905.1342195184993*C_lqM + 1571.9428762098114*C_eq
 + 429.0079598942002*C_lu + 6315.299860874898*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_P80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_ilc_1000_P20_P80_8 = -51.76068979603156*C_phiQm + -5.071716399604718*C_phit
 + -3278.2485728256847*C_tW + 3370.8570331613587*C_tZ
 + -0.0007729415946775806*CI_tW + 0.0006354777989346247*CI_tZ
 + 1571.9428762098114*C_lqM + 15001.898179534715*C_eq
 + 4029.175309281545*C_lu + 5888.744683321876*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_P80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_ilc_1000_P20_P80_9 = -13.874493035865024*C_phiQm + -0.3002445465021264*C_phit
 + -976.7877661608333*C_tW + 989.428846361019*C_tZ
 + -3.858758917118177e-06*CI_tW + -9.060090807748996e-05*CI_tZ
 + 429.00795989420044*C_lqM + 4029.175309281545*C_eq
 + 1566.5073741320496*C_lu + 1861.67840204391*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_P80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_ilc_1000_P20_P80_10 = 5.112067117457315*C_phiQm + -111.67142836604903*C_phit
 + -5056.925263873175*C_tW + 5242.634125639533*C_tZ
 + 0.0002807407206201645*CI_tW + -0.00020140504172940518*CI_tZ
 + 6315.299860874898*C_lqM + 5888.744683321876*C_eq
 + 1861.67840204391*C_lu + 30612.688463381062*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_P20_P80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_ilc_1000_P20_P80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}






opt_obs_ilc_1000_M30_P80::opt_obs_ilc_1000_M30_P80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_ilc_1000_M30_P80" << "eff_ilc_1000_M30_P80");               
}

double opt_obs_ilc_1000_M30_P80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_ilc_1000_M30_P80");
    double eff = SM.getOptionalParameter("eff_ilc_1000_M30_P80");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_ilc_1000_M30_P80_1 = 0.43927175900014614*C_phiQm + 0.04933441339020023*C_phit
 + 19.34055613625289*C_tW + -21.161997441957084*C_tZ
 + 1.954224664942196e-06*CI_tW + -4.3904978309551674e-05*CI_tZ
 + -2.2187989430850106*C_lqM + -114.61434800174523*C_eq
 + -11.496246670061486*C_lu + -25.68949672116183*C_eu;

        return lumi*eff*opt_obs_ilc_1000_M30_P80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_ilc_1000_M30_P80_2 = 0.04933441339020023*C_phiQm + 0.8514654351901008*C_phit
 + 30.771855385912957*C_tW + -33.82854174700178*C_tZ
 + -2.6017446768034072e-06*CI_tW + -3.98543108449264e-05*CI_tZ
 + -17.14307762660929*C_lqM + -32.12245971859595*C_eq
 + -3.205997389379503*C_lu + -220.91488873896276*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M30_P80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_ilc_1000_M30_P80_3 = 19.34055613625289*C_phiQm + 30.77185538591296*C_phit
 + 2095.2849101593974*C_tW + -2293.7630902938713*C_tZ
 + 4.45194887783085e-14*CI_tW + -5.3905963051372404e-14*CI_tZ
 + -723.4534617927372*C_lqM + -5759.250645394546*C_eq
 + -616.6371418717267*C_lu + -8554.72525423054*C_eu;

        return lumi*eff*opt_obs_ilc_1000_M30_P80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_ilc_1000_M30_P80_4 = -21.161997441957084*C_phiQm + -33.82854174700178*C_phit
 + -2293.763090293871*C_tW + 2511.8208071903314*C_tZ
 + -3.926181008951456e-14*CI_tW + 4.8080422270969955e-14*CI_tZ
 + 799.9934078600339*C_lqM + 6311.267870544609*C_eq
 + 671.2206722256877*C_lu + 9396.899882666941*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M30_P80_4;

        
    } else if(b_min == 5 && b_max == 5){
        
        double opt_obs_ilc_1000_M30_P80_5 = 1.954224664942196e-06*C_phiQm + -2.6017446768034072e-06*C_phit
 + 4.45194887783085e-14*C_tW + -3.926181008951457e-14*C_tZ
 + 258.7284225630976*CI_tW + -283.9081797552015*CI_tZ
 + 1.984926683376991e-05*C_lqM + -0.0008377744554338955*C_eq
 + 2.5102474309473477e-05*C_lu + 0.0005860120479975471*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M30_P80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_ilc_1000_M30_P80_6 = -4.390497830955167e-05*C_phiQm + -3.98543108449264e-05*C_phit
 + -5.3905963051372404e-14*C_tW + 4.808042227096996e-14*C_tZ
 + -283.9081797552015*CI_tW + 311.5703244407138*CI_tZ
 + 5.6137770578445584e-05*C_lqM + -0.0009105711352818717*C_eq
 + -0.00011899538895231577*C_lu + 0.0011337309881414232*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M30_P80_6;

        
    } else if(b_min == 7 && b_max == 7){
                
        double opt_obs_ilc_1000_M30_P80_7 = -2.2187989430850124*C_phiQm + -17.14307762660929*C_phit
 + -723.4534617927372*C_tW + 799.993407860034*C_tZ
 + 1.984926683376991e-05*CI_tW + 5.61377705784456e-05*CI_tZ
 + 514.9266357385035*C_lqM + 1152.8930557807962*C_eq
 + 113.28019339128403*C_lu + 4503.356918254478*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M30_P80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_ilc_1000_M30_P80_8 = -114.61434800174523*C_phiQm + -32.12245971859595*C_phit
 + -5759.250645394545*C_tW + 6311.267870544609*C_tZ
 + -0.0008377744554338955*CI_tW + -0.0009105711352818717*CI_tZ
 + 1152.8930557807957*C_lqM + 30559.050316336587*C_eq
 + 3062.5598538252675*C_lu + 11686.585514343278*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M30_P80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_ilc_1000_M30_P80_9 = -11.496246670061488*C_phiQm + -3.205997389379502*C_phit
 + -616.6371418717266*C_tW + 671.2206722256878*C_tZ
 + 2.5102474309473477e-05*CI_tW + -0.00011899538895231577*CI_tZ
 + 113.28019339128409*C_lqM + 3062.5598538252675*C_eq
 + 448.02426261325695*C_lu + 1328.6522683361013*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M30_P80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_ilc_1000_M30_P80_10 = -25.689496721161845*C_phiQm + -220.9148887389628*C_phit
 + -8554.72525423054*C_tW + 9396.899882666941*C_tZ
 + 0.0005860120479975471*CI_tW + 0.0011337309881414232*CI_tZ
 + 4503.356918254478*C_lqM + 11686.585514343282*C_eq
 + 1328.6522683361018*C_lu + 57883.096059097006*C_eu;
        
        return lumi*eff*opt_obs_ilc_1000_M30_P80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_ilc_1000_M30_P80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}






opt_obs_ilc_1000_P30_M80::opt_obs_ilc_1000_P30_M80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_ilc_1000_P30_M80" << "eff_ilc_1000_P30_M80");                  
}

double opt_obs_ilc_1000_P30_M80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_ilc_1000_P30_M80");
    double eff = SM.getOptionalParameter("eff_ilc_1000_P30_M80");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
                
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_ilc_1000_P30_M80_1 = 1.1987886923974536*C_phiQm + 0.22117351869803004*C_phit
 + -38.96249738223764*C_tW + 23.53768918125916*C_tZ
 + 2.937825507109394e-06*CI_tW + -2.67760771171347e-06*CI_tZ
 + 270.30657017078545*C_lqM + 2.126715524302746*C_eq
 + 59.59753043673251*C_lu + 11.47142726960921*C_eu;

        return lumi*eff*opt_obs_ilc_1000_P30_M80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_ilc_1000_P30_M80_2 = 0.2211735186980302*C_phiQm + 0.8004634375502897*C_phit
 + -29.46841972072286*C_tW + 17.74022521362101*C_tZ
 + 3.60746888280997e-06*CI_tW + -3.0903052729717966e-06*CI_tZ
 + 54.11193638652344*C_lqM + 5.228895174196872*C_eq
 + 181.51150553903307*C_lu + 2.5865412003729755*C_eu;

        
        return lumi*eff*opt_obs_ilc_1000_P30_M80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_ilc_1000_P30_M80_3 = -38.962497382237636*C_phiQm + -29.468419720722856*C_phit
 + 2188.8655666552872*C_tW + -1318.2381311494682*C_tZ
 + 0.0*CI_tW + -8.873447611166219e-15*CI_tZ
 + -8918.469843892013*C_lqM + -220.82454357474063*C_eq
 + -6946.411009837359*C_lu + -396.1261219170531*C_eu;


        return lumi*eff*opt_obs_ilc_1000_P30_M80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_ilc_1000_P30_M80_4 = 23.537689181259154*C_phiQm + 17.740225213621002*C_phit
 + -1318.238131149468*C_tW + 795.4505826087668*C_tZ
 + 1.7786420601446763e-16*CI_tW + 6.997840156803055e-15*CI_tZ
 + 5388.540993832447*C_lqM + 134.29716002121233*C_eq
 + 4189.301836106421*C_lu + 246.9864294852694*C_eu;

        
        return lumi*eff*opt_obs_ilc_1000_P30_M80_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_ilc_1000_P30_M80_5 = 2.9378255071093934e-06*C_phiQm + 3.60746888280997e-06*C_phit
 + 0.0*C_tW + 1.7786420601446766e-16*C_tZ
 + 273.84336969425925*CI_tW + -165.46618543761434*CI_tZ
 + 3.4300843848702727e-06*C_lqM + 0.0007084082711695386*C_eq
 + 2.6420151642797557e-05*C_lu + 0.0007092667233139591*C_eu;

        
        return lumi*eff*opt_obs_ilc_1000_P30_M80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_ilc_1000_P30_M80_6 = -2.6776077117134704e-06*C_phiQm + -3.090305272971797e-06*C_phit
 + -8.873447611166219e-15*C_tW + 6.997840156803056e-15*C_tZ
 + -165.46618543761434*CI_tW + 99.98711858937997*CI_tZ
 + -0.0005245214244543577*C_lqM + -0.0006184121062542568*C_eq
 + 0.0004987865070970718*C_lu + -0.0006315769510948558*C_eu;

        
        return lumi*eff*opt_obs_ilc_1000_P30_M80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        
        double opt_obs_ilc_1000_P30_M80_7 = 270.30657017078545*C_phiQm + 54.11193638652344*C_phit
 + -8918.469843892013*C_tW + 5388.5409938324465*C_tZ
 + 3.4300843848702714e-06*CI_tW + -0.0005245214244543577*CI_tZ
 + 60984.31291778166*C_lqM + 520.006778050391*C_eq
 + 14389.642977270483*C_lu + 2588.651906505266*C_eu;

        
        return lumi*eff*opt_obs_ilc_1000_P30_M80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_ilc_1000_P30_M80_8 = 2.126715524302746*C_phiQm + 5.228895174196872*C_phit
 + -220.82454357474063*C_tW + 134.2971600212124*C_tZ
 + 0.0007084082711695386*CI_tW + -0.0006184121062542568*CI_tZ
 + 520.006778050391*C_lqM + 49.89417802641793*C_eq
 + 1190.7915073724098*C_lu + 22.68946156305393*C_eu;

        
        return lumi*eff*opt_obs_ilc_1000_P30_M80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_ilc_1000_P30_M80_9 = 59.597530436732534*C_phiQm + 181.51150553903307*C_phit
 + -6946.41100983736*C_tW + 4189.3018361064205*C_tZ
 + 2.6420151642797557e-05*CI_tW + 0.0004987865070970719*CI_tZ
 + 14389.642977270483*C_lqM + 1190.7915073724098*C_eq
 + 41278.975813086436*C_lu + 722.0597190881294*C_eu;

        
        return lumi*eff*opt_obs_ilc_1000_P30_M80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_ilc_1000_P30_M80_10 = 11.47142726960921*C_phiQm + 2.5865412003729755*C_phit
 + -396.12612191705307*C_tW + 246.9864294852694*C_tZ
 + 0.000709266723313959*CI_tW + -0.0006315769510948557*CI_tZ
 + 2588.651906505266*C_lqM + 22.68946156305394*C_eq
 + 722.0597190881294*C_lu + 162.8033954813007*C_eu;

        
        return lumi*eff*opt_obs_ilc_1000_P30_M80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_ilc_1000_P30_M80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}









opt_obs_clic_380_0_M80::opt_obs_clic_380_0_M80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_clic_380_0_M80" << "eff_clic_380_0_M80");                             
}

double opt_obs_clic_380_0_M80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_clic_380_0_M80");
    double eff = SM.getOptionalParameter("eff_clic_380_0_M80");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
                
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_clic_380_0_M80_1 = 2.4012344013789924*C_phiQm + 1.905265776860176*C_phit
 + -90.62194879632287*C_tW + 54.933017584564254*C_tZ
 + -1.510179528126631e-07*CI_tW + 1.0410629898214733e-07*CI_tZ
 + 77.74087527712051*C_lqM + 4.506516485060968*C_eq
 + 63.52967970860912*C_lu + 5.67926424309723*C_eu;

        return lumi*eff*opt_obs_clic_380_0_M80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_clic_380_0_M80_2 = 1.9052657768601764*C_phiQm + 2.1375471670361286*C_phit
 + -84.96212312728132*C_tW + 51.415831922461194*C_tZ
 + 1.510179528126631e-07*CI_tW + -1.0410629911245894e-07*CI_tZ
 + 62.72264618582824*C_lqM + 4.759631849474812*C_eq
 + 69.67230618240092*C_lu + 4.545296874962294*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_M80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_clic_380_0_M80_3 = -90.62194879632287*C_phiQm + -84.96212312728133*C_phit
 + 3708.4886305328705*C_tW + -2246.22219875002*C_tZ
 + 1.1737718321964578e-15*CI_tW + 6.167796862743943e-17*CI_tZ
 + -2958.4494235079305*C_lqM + -198.0143710095109*C_eq
 + -2802.628244967366*C_lu + -218.63935047010096*C_eu;


        return lumi*eff*opt_obs_clic_380_0_M80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_clic_380_0_M80_4 = 54.93301758456427*C_phiQm + 51.41583192246122*C_phit
 + -2246.2221987500207*C_tW + 1366.012117675563*C_tZ
 + 5.952705866708424e-16*CI_tW + -1.553058444004611e-16*CI_tZ
 + 1795.3542595860504*C_lqM + 122.31257159197723*C_eq
 + 1699.1594628392736*C_lu + 135.86305981983088*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_M80_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_clic_380_0_M80_5 = -1.510179528126631e-07*C_phiQm + 1.510179528126631e-07*C_phit
 + 1.1737718321964578e-15*C_tW + 5.952705866708425e-16*C_tZ
 + 14.577980749417915*CI_tW + -8.853789038374096*CI_tZ
 + -8.54576825429878e-07*C_lqM + -6.808031644695242e-05*C_eq
 + 8.545768275524168e-07*C_lu + 6.808031644788933e-05*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_M80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_clic_380_0_M80_6 = 1.0410629898214732e-07*C_phiQm + -1.0410629911245892e-07*C_phit
 + 6.167796862743938e-17*C_tW + -1.5530584440046111e-16*C_tZ
 + -8.853789038374096*CI_tW + 5.378408517442349*CI_tZ
 + 6.877120194308497e-05*C_lqM + 4.607666230663461e-05*C_eq
 + -6.877120194495685e-05*C_lu + -4.6076662307460866e-05*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_M80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        
        double opt_obs_clic_380_0_M80_7 = 77.74087527712051*C_phiQm + 62.72264618582823*C_phit
 + -2958.44942350793*C_tW + 1795.3542595860506*C_tZ
 + -8.545768254298782e-07*CI_tW + 6.877120194308497e-05*CI_tZ
 + 2520.38115848052*C_lqM + 149.88262352581376*C_eq;

        
        return lumi*eff*opt_obs_clic_380_0_M80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_clic_380_0_M80_8 = 4.506516485060968*C_phiQm + 4.759631849474812*C_phit
 + -198.0143710095109*C_tW + 122.31257159197723*C_tZ
 + -6.808031644695244e-05*CI_tW + 4.607666230663461e-05*CI_tZ
 + 149.88262352581376*C_lqM + 12.995584678830621*C_eq
 + 157.39902419389858*C_lu + 12.69703053936366*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_M80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_clic_380_0_M80_9 = 63.52967970860913*C_phiQm + 69.67230618240093*C_phit
 + -2802.628244967366*C_tW + 1699.1594628392731*C_tZ
 + 8.545768275524168e-07*CI_tW + -6.877120194495683e-05*CI_tZ
 + 2090.295913596953*C_lqM + 157.39902419389855*C_eq
 + 2276.4585542807317*C_lu + 154.45107386236768*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_M80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_clic_380_0_M80_10 = 5.679264243097231*C_phiQm + 4.545296874962295*C_phit
 + -218.63935047010094*C_tW + 135.86305981983088*C_tZ
 + 6.808031644788933e-05*CI_tW + -4.6076662307460866e-05*CI_tZ
 + 185.6572021876204*C_lqM + 12.69703053936366*C_eq
 + 154.45107386236768*C_lu + 16.843302387395145*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_M80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_clic_380_0_M80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}





opt_obs_clic_380_0_P80::opt_obs_clic_380_0_P80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_clic_380_0_P80" << "eff_clic_380_0_P80");             
}

double opt_obs_clic_380_0_P80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    double lumi = SM.getOptionalParameter("lumi_clic_380_0_P80");
    double eff = SM.getOptionalParameter("eff_clic_380_0_P80");
    
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
            
                
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_clic_380_0_P80_1 = 1.0341761741381041*C_phiQm + 0.8678278844204601*C_phit
 + 56.68893403665082*C_tW + -60.871723947002224*C_tZ
 + -8.783522025136684e-07*CI_tW + 2.4612823960247657e-06*CI_tZ
 + -6.020033435067013*C_lqM + -43.071797425225654*C_eq
 + -6.984378895746486*C_lu + -38.34578687274316*C_eu;

        return lumi*eff*opt_obs_clic_380_0_P80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_clic_380_0_P80_2 = 0.8678278844204601*C_phiQm + 1.2800040237569856*C_phit
 + 63.24723568308892*C_tW + -68.04839784248016*C_tZ
 + 8.783522026203955e-07*CI_tW + -2.461282396160697e-06*CI_tZ
 + -7.974170072066974*C_lqM + -39.473715433271835*C_eq
 + -6.2889840875507055*C_lu + -51.985960903184775*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_P80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_clic_380_0_P80_3 = 56.68893403665082*C_phiQm + 63.24723568308892*C_phit
 + 3637.768424262441*C_tW + -3898.3399469046262*C_tZ
 + 7.157912155269471e-15*CI_tW + -2.5181380615073766e-15*CI_tZ
 + -445.3252512598177*C_lqM + -2492.4269919420503*C_eq
 + -426.14005905002784*C_lu + -2700.2034311245907*C_eu;


        return lumi*eff*opt_obs_clic_380_0_P80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_clic_380_0_P80_4 = -60.87172394700221*C_phiQm + -68.04839784248014*C_phit
 + -3898.339946904626*C_tW + 4181.759795045329*C_tZ
 + -6.312574456482685e-15*CI_tW + 3.909617197126763e-05*CI_tZ
 + 475.9784437679133*C_lqM + 2673.819071838232*C_eq
 + 454.24304876930967*C_lu + 2900.3381287370885*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_P80_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_clic_380_0_P80_5 = -8.783522025136682e-07*C_phiQm + 8.783522026203954e-07*C_phit
 + 7.157912155269471e-15*C_tW + -6.312574456482686e-15*C_tZ
 + 14.234915995947016*CI_tW + -15.268652477864318*CI_tZ
 + -1.8531641033033386e-06*C_lqM + 0.019574895100669603*C_eq
 + 1.8531641049344373e-06*C_lu + -0.019574895100669603*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_P80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_clic_380_0_P80_6 = 2.461282396024766e-06*C_phiQm + -2.461282396160697e-06*C_phit
 + -2.518138061507376e-15*C_tW + 3.909617197126763e-05*C_tZ
 + -15.26865247786432*CI_tW + 16.381275342798084*CI_tZ
 + -6.937747711679283e-06*C_lqM + -0.01726360807033705*C_eq
 + 6.937747710479603e-06*C_lu + 0.01726360807033528*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_P80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        
        double opt_obs_clic_380_0_P80_7 = -6.020033435067013*C_phiQm + -7.974170072066974*C_phit
 + -445.32525125981766*C_tW + 475.9784437679133*C_tZ
 + -1.8531641033033386e-06*CI_tW + -6.937747711679283e-06*CI_tZ
 + 65.74424987845913*C_lqM + 285.70874598570134*C_eq
 + 54.46067778721985*C_lu + 341.2755822181571*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_P80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_clic_380_0_P80_8 = -43.071797425225654*C_phiQm + -39.473715433271835*C_phit
 + -2492.4269919420503*C_tW + 2673.819071838232*C_tZ
 + 0.019574895100669603*CI_tW + -0.01726360807033705*CI_tZ
 + 285.70874598570134*C_lqM + 1833.7299868067848*C_eq
 + 306.6170831239439*C_lu + 1731.5740587852574*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_P80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_clic_380_0_P80_9 = -6.984378895746486*C_phiQm + -6.2889840875507055*C_phit
 + -426.1400590500278*C_tW + 454.2430487693097*C_tZ
 + 1.8531641049344373e-06*CI_tW + 6.937747710479601e-06*CI_tZ
 + 54.46067778721985*C_lqM + 306.617083123944*C_eq
 + 61.28659874963077*C_lu + 290.0472349059146*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_P80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_clic_380_0_P80_10 = -38.34578687274317*C_phiQm + -51.98596090318478*C_phit
 + -2700.20343112459*C_tW + 2900.338128737089*C_tZ
 + -0.019574895100669603*CI_tW + 0.01726360807033528*CI_tZ
 + 341.2755822181572*C_lqM + 1731.5740587852576*C_eq
 + 290.0472349059146*C_lu + 2150.8037263028136*C_eu;

        
        return lumi*eff*opt_obs_clic_380_0_P80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_clic_380_0_P80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}











opt_obs_clic_1500_0_M80::opt_obs_clic_1500_0_M80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_clic_1500_0_M80" << "eff_clic_1500_0_M80");             
}

double opt_obs_clic_1500_0_M80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    double lumi = SM.getOptionalParameter("lumi_clic_1500_0_M80");
    double eff = SM.getOptionalParameter("eff_clic_1500_0_M80");
    
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    
     
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_clic_1500_0_M80_1 = 0.41494184774765785*C_phiQm + 0.03656147200861573*C_phit
 + -13.65373715985302*C_tW + 8.498266992604943*C_tZ
 + 9.819943953004157e-07*CI_tW + -1.4466765561111347e-06*CI_tZ
 + 211.882605565019*C_lqM + 2.1043144313805966*C_eq
 + 33.06916555854299*C_lu + 16.59456121668658*C_eu;

        return lumi*eff*opt_obs_clic_1500_0_M80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_clic_1500_0_M80_2 = 0.036561472008615745*C_phiQm + 0.27522284271636205*C_phit
 + -9.90826667162875*C_tW + 6.113613180121138*C_tZ
 + -1.0544602060381405e-06*CI_tW + 1.5233798660043683e-06*CI_tZ
 + 24.92859321957831*C_lqM + 7.318547415731294*C_eq
 + 141.24661742701258*C_lu + 2.202934119038063*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_M80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_clic_1500_0_M80_3 = -13.65373715985302*C_phiQm + -9.908266671628748*C_phit
 + 983.4038953477036*C_tW + -610.3086040133214*C_tZ
 + 0.0*CI_tW + -6.126433762720552e-16*CI_tZ
 + -7190.616243375866*C_lqM + -318.4023654269937*C_eq
 + -5543.014431477027*C_lu + -601.1466121685295*C_eu;


        return lumi*eff*opt_obs_clic_1500_0_M80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_clic_1500_0_M80_4 = 8.498266992604943*C_phiQm + 6.113613180121138*C_phit
 + -610.3086040133215*C_tW + 379.5373865291553*C_tZ
 + -0.0006419359240135993*CI_tW + 0.0007007006238664515*CI_tZ
 + 4476.625314511583*C_lqM + 199.44349700731678*C_eq
 + 3435.179815184874*C_lu + 387.99197436867405*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_M80_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_clic_1500_0_M80_5 = 9.819943953004157e-07*C_phiQm + -1.0544602060381403e-06*C_phit
 + 0.0*C_tW + -0.0006419359240135993*C_tZ
 + 253.44366233030155*CI_tW + -157.80606655450362*CI_tZ
 + -3.920848322499579e-06*C_lqM + -9.55163508443112e-07*C_eq
 + 3.920847736100021e-06*C_lu + 2.9054969738236944e-07*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_M80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_clic_1500_0_M80_6 = -1.4466765561111345e-06*C_phiQm + 1.523379866004368e-06*C_phit
 + -6.126433762720552e-16*C_tW + 0.0007007006238664515*C_tZ
 + -157.8060665545036*CI_tW + 98.28503464206874*CI_tZ
 + 5.709251467507465e-06*C_lqM + 4.889102804076378e-05*C_eq
 + -3.153715364212929e-06*C_lu + 4.7432612791321025e-05*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_M80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        
        double opt_obs_clic_1500_0_M80_7 = 211.882605565019*C_phiQm + 24.928593219578303*C_phit
 + -7190.616243375867*C_tW + 4476.625314511584*C_tZ
 + -3.9208483224995796e-06*CI_tW + 5.709251467507464e-06*CI_tZ
 + 108409.56096180479*C_lqM + 1322.4110518226498*C_eq
 + 20068.56442424916*C_lu + 8488.520946719624*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_M80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_clic_1500_0_M80_8 = 2.104314431380597*C_phiQm + 7.318547415731295*C_phit
 + -318.4023654269937*C_tW + 199.4434970073168*C_tZ
 + -9.55163508443112e-07*CI_tW + 4.889102804076377e-05*CI_tZ
 + 1322.4110518226503*C_lqM + 293.12835700404423*C_eq
 + 3793.0698383271533*C_lu + 100.99041914328296*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_M80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_clic_1500_0_M80_9 = 33.06916555854296*C_phiQm + 141.2466174270126*C_phit
 + -5543.014431477027*C_tW + 3435.1798151848743*C_tZ
 + 3.920847736100021e-06*CI_tW + -3.1537153642129295e-06*CI_tZ
 + 20068.56442424916*C_lqM + 3793.069838327153*C_eq
 + 73242.65698580653*C_lu + 1991.2434471603299*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_M80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_clic_1500_0_M80_10 = 16.59456121668658*C_phiQm + 2.202934119038064*C_phit
 + -601.1466121685295*C_tW + 387.9919743686741*C_tZ
 + 2.9054969738236944e-07*CI_tW + 4.7432612791321025e-05*CI_tZ
 + 8488.520946719624*C_lqM + 100.9904191432829*C_eq
 + 1991.2434471603297*C_lu + 998.6055794619215*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_M80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_clic_1500_0_M80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}





opt_obs_clic_1500_0_P80::opt_obs_clic_1500_0_P80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_clic_1500_0_P80" << "eff_clic_1500_0_P80");              
}

double opt_obs_clic_1500_0_P80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_clic_1500_0_P80");
    double eff = SM.getOptionalParameter("eff_clic_1500_0_P80");
    
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_clic_1500_0_P80_1 = 0.13765536346506738*C_phiQm + -0.03150972988950435*C_phit
 + 5.415254729903891*C_tW + -5.727837805228482*C_tZ
 + 5.008320032421266e-06*CI_tW + -4.41594799747584e-06*CI_tZ
 + 0.8494657886761665*C_lqM + -78.43138802720634*C_eq
 + -14.457640452249116*C_lu + 1.6993611998350175*C_eu;

        return lumi*eff*opt_obs_clic_1500_0_P80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_clic_1500_0_P80_2 = -0.03150972988950436*C_phiQm + 0.2809542282743418*C_phit
 + 9.904548569070963*C_tW + -10.586930607190805*C_tZ
 + 2.7422791882167675e-06*CI_tW + -2.4103893586063074e-06*CI_tZ
 + -22.713729787618487*C_lqM + -7.7078658864503495*C_eq
 + -0.99787333491129*C_lu + -163.18698886414813*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_P80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_clic_1500_0_P80_3 = 5.415254729903891*C_phiQm + 9.904548569070965*C_phit
 + 937.5081186264679*C_tW + -993.1759797993705*C_tZ
 + 0.0*CI_tW + 0.0*CI_tZ
 + -1039.7580915262515*C_lqM + -4308.6983721686265*C_eq
 + -875.2355489629199*C_lu + -6710.154605340167*C_eu;


        return lumi*eff*opt_obs_clic_1500_0_P80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_clic_1500_0_P80_4 = -5.727837805228482*C_phiQm + -10.586930607190805*C_phit
 + -993.1759797993707*C_tW + 1052.5262689769156*C_tZ
 + 0.00010553220750169323*CI_tW + -9.306902668688913e-05*CI_tZ
 + 1107.1336951691812*C_lqM + 4565.800796489775*C_eq
 + 916.241453056005*C_lu + 7150.431571269555*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_P80_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_clic_1500_0_P80_5 = 5.008320032421266e-06*C_phiQm + 2.742279188216767e-06*C_phit
 + 0.0*C_tW + 0.00010553220750169322*C_tZ
 + 238.85827237183304*CI_tW + -253.6819319472234*CI_tZ
 + -0.0025324700874964135*C_lqM + -0.024448486436486467*C_eq
 + -0.002750089357457568*C_lu + 0.020420744007673384*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_P80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_clic_1500_0_P80_6 = -4.41594799747584e-06*C_phiQm + -2.4103893586063074e-06*C_phit
 + 0.0*C_tW + -9.306902668688913e-05*C_tZ
 + -253.6819319472234*CI_tW + 269.5043011809208*CI_tZ
 + 0.002265787682125453*C_lqM + 0.02697293325709436*C_eq
 + 0.002469140376901939*C_lu + -0.022900450054608775*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_P80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        
        double opt_obs_clic_1500_0_P80_7 = 0.8494657886761665*C_phiQm + -22.713729787618487*C_phit
 + -1039.7580915262515*C_tW + 1107.1336951691812*C_tZ
 + -0.0025324700874964135*CI_tW + 0.002265787682125453*CI_tZ
 + 2838.2963573053908*C_lqM + 2744.380306293122*C_eq
 + 483.2540222910784*C_lu + 13651.851199371788*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_P80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_clic_1500_0_P80_8 = -78.43138802720632*C_phiQm + -7.7078658864503495*C_phit
 + -4308.6983721686265*C_tW + 4565.800796489775*C_tZ
 + -0.024448486436486467*CI_tW + 0.02697293325709436*CI_tZ
 + 2744.380306293121*C_lqM + 48364.33282617488*C_eq
 + 8889.612536017306*C_lu + 14577.239175981536*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_P80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_clic_1500_0_P80_9 = -14.457640452249116*C_phiQm + -0.9978733349112892*C_phit
 + -875.23554896292*C_tW + 916.2414530560051*C_tZ
 + -0.002750089357457568*CI_tW + 0.002469140376901939*CI_tZ
 + 483.2540222910784*C_lqM + 8889.612536017306*C_eq
 + 2424.7196180296273*C_lu + 3338.5023956991186*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_P80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_clic_1500_0_P80_10 = 1.6993611998350104*C_phiQm + -163.18698886414813*C_phit
 + -6710.154605340167*C_tW + 7150.431571269555*C_tZ
 + 0.020420744007673387*CI_tW + -0.02290045005460878*CI_tZ
 + 13651.851199371788*C_lqM + 14577.239175981544*C_eq
 + 3338.502395699118*C_lu + 97928.88435870936*C_eu;

        
        return lumi*eff*opt_obs_clic_1500_0_P80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_clic_1500_0_P80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}







opt_obs_clic_3000_0_M80::opt_obs_clic_3000_0_M80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_clic_3000_0_M80" << "eff_clic_3000_0_M80");               
}

double opt_obs_clic_3000_0_M80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    double lumi = SM.getOptionalParameter("lumi_clic_3000_0_M80");
    double eff = SM.getOptionalParameter("eff_clic_3000_0_M80");
    
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    
                
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_clic_3000_0_M80_1 = 0.10665675620537152*C_phiQm + 0.006120265198131997*C_phit
 + -3.481924152921489*C_tW + 2.169906794631394*C_tZ
 + -4.789674072427836e-07*CI_tW + 1.1585576988050365e-06*CI_tZ
 + 218.04039029903777*C_lqM + 1.6984007351286863*C_eq
 + 27.448088960492342*C_lu + 17.117247022155606*C_eu;

        return lumi*eff*opt_obs_clic_3000_0_M80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_clic_3000_0_M80_2 = 0.006120265198131997*C_phiQm + 0.07063538850516685*C_phit
 + -2.506685190387274*C_tW + 1.5480484890040143*C_tZ
 + -4.680659731560242e-07*CI_tW + 1.150641800849326e-06*CI_tZ
 + 18.907743964383243*C_lqM + 7.385104812715184*C_eq
 + 144.97604010771937*C_lu + 1.7808526707912509*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_M80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_clic_3000_0_M80_3 = -3.481924152921489*C_phiQm + -2.506685190387275*C_phit
 + 465.8432929446417*C_tW + -289.7891333524571*C_tZ
 + 0.0*CI_tW + 3.033572847024531e-15*CI_tZ
 + -7347.8538201048395*C_lqM + -317.1353164125313*C_eq
 + -5632.792600979696*C_lu + -619.274494399182*C_eu;


        return lumi*eff*opt_obs_clic_3000_0_M80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_clic_3000_0_M80_4 = 2.169906794631394*C_phiQm + 1.5480484890040143*C_phit
 + -289.7891333524571*C_tW + 180.46352022095024*C_tZ
 + -2.088732261424171e-05*CI_tW + 3.434654155607854e-05*CI_tZ
 + 4580.121750561242*C_lqM + 198.76084341562716*C_eq
 + 3494.3285731977794*C_lu + 400.3736201208188*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_M80_4;

        
    } else if(b_min == 5 && b_max == 5){
   
        double opt_obs_clic_3000_0_M80_5 = -4.789674072427836e-07*C_phiQm + -4.6806597315602414e-07*C_phit
 + 0.0*C_tW + -2.088732261424171e-05*C_tZ
 + 275.8250360851673*CI_tW + -171.93730810726458*CI_tZ
 + -9.747407510863218e-07*C_lqM + -0.012506306684315698*C_eq
 + 9.747401722121282e-07*C_lu + -0.0033981157188001463*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_M80_5;

        
    } else if(b_min == 6 && b_max == 6){
        
        double opt_obs_clic_3000_0_M80_6 = 1.1585576988050365e-06*C_phiQm + 1.150641800849326e-06*C_phit
 + 3.0335728470245314e-15*C_tW + 3.434654155607854e-05*C_tZ
 + -171.93730810726458*CI_tW + 107.20570866252996*CI_tZ
 + 0.0016272086607826257*C_lqM + 0.010846896608532733*C_eq
 + -0.001626890644555206*C_lu + 0.003182417010832638*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_M80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        
        double opt_obs_clic_3000_0_M80_7 = 218.04039029903777*C_phiQm + 18.907743964383258*C_phit
 + -7347.8538201048395*C_tW + 4580.121750561242*C_tZ
 + -9.747407510863218e-07*CI_tW + 0.0016272086607826257*CI_tZ
 + 446634.16778215754*C_lqM + 4482.935777858627*C_eq
 + 69136.5320758395*C_lu + 35038.35047738055*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_M80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_clic_3000_0_M80_8 = 1.6984007351286863*C_phiQm + 7.385104812715186*C_phit
 + -317.1353164125313*C_tW + 198.76084341562716*C_tZ
 + -0.012506306684315698*CI_tW + 0.010846896608532733*CI_tZ
 + 4482.935777858625*C_lqM + 1178.9043320654328*C_eq
 + 15279.184826511586*C_lu + 323.94601814530847*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_M80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_clic_3000_0_M80_9 = 27.44808896049235*C_phiQm + 144.97604010771937*C_phit
 + -5632.792600979696*C_tW + 3494.3285731977803*C_tZ
 + 9.747401722121282e-07*CI_tW + -0.001626890644555206*CI_tZ
 + 69136.5320758395*C_lqM + 15279.184826511586*C_eq
 + 300735.0999930231*C_lu + 7274.877296483532*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_M80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_clic_3000_0_M80_10 = 17.117247022155606*C_phiQm + 1.7808526707912513*C_phit
 + -619.274494399182*C_tW + 400.3736201208188*C_tZ
 + -0.003398115718800146*CI_tW + 0.003182417010832637*CI_tZ
 + 35038.35047738055*C_lqM + 323.94601814530824*C_eq
 + 7274.877296483532*C_lu + 4169.67504212753*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_M80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_clic_3000_0_M80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}





opt_obs_clic_3000_0_P80::opt_obs_clic_3000_0_P80(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_clic_3000_0_P80" << "eff_clic_3000_0_P80");               
}

double opt_obs_clic_3000_0_P80::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_clic_3000_0_P80");
    double eff = SM.getOptionalParameter("eff_clic_3000_0_P80");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_clic_3000_0_P80_1 = 0.03557555683714805*C_phiQm + -0.01019018670122251*C_phit
 + 1.3548183058442778*C_tW + -1.4319477759528147*C_tZ
 + 7.586879084393368e-06*CI_tW + -7.551997580429217e-06*CI_tZ
 + 1.6429280154454258*C_lqM + -80.4314129733981*C_eq
 + -14.933801174827105*C_lu + 6.559065280624168*C_eu;

        return lumi*eff*opt_obs_clic_3000_0_P80_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_clic_3000_0_P80_2 = -0.010190186701222508*C_phiQm + 0.07285212051112055*C_phit
 + 2.5268062513821437*C_tW + -2.700117566348261*C_tZ
 + -8.044958901495904e-06*CI_tW + 6.217092338917716e-06*CI_tZ
 + -23.418581880871088*C_lqM + -3.1112272836851673*C_eq
 + -0.20474634142293802*C_lu + -168.78343584875955*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_P80_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_clic_3000_0_P80_3 = 1.3548183058442778*C_phiQm + 2.5268062513821437*C_phit
 + 441.6306475259417*C_tW + -468.094575317429*C_tZ
 + -2.2416945973823423e-14*CI_tW + 0.00018589028351812944*CI_tZ
 + -1060.972953308112*C_lqM + -4343.545384483194*C_eq
 + -895.1839379354278*C_lu + -6866.053235950871*C_eu;


        return lumi*eff*opt_obs_clic_3000_0_P80_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_clic_3000_0_P80_4 = -1.4319477759528145*C_phiQm + -2.700117566348261*C_phit
 + -468.094575317429*C_tW + 496.2993870967572*C_tZ
 + 2.627146624412578e-14*CI_tW + -0.00017595413617225793*CI_tZ
 + 1129.5794504410328*C_lqM + 4600.191855822186*C_eq
 + 936.1159110158246*C_lu + 7313.668222089809*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_P80_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_clic_3000_0_P80_5 = 7.586879084393368e-06*C_phiQm + -8.044958901495906e-06*C_phit
 + -2.2416945973823423e-14*C_tW + 2.6271466244125783e-14*C_tZ
 + 259.6509178980169*CI_tW + -275.6324331373924*CI_tZ
 + -0.002502097678094232*C_lqM + 0.0005387824547244251*C_eq
 + 0.0024130670808982227*C_lu + 0.0005287454515653952*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_P80_5;

        
    } else if(b_min == 6 && b_max == 6){
        

        double opt_obs_clic_3000_0_P80_6 = -7.551997580429217e-06*C_phiQm + 6.217092338917716e-06*C_phit
 + 0.00018589028351812944*C_tW + -0.00017595413617225793*C_tZ
 + -275.6324331373924*CI_tW + 292.6766259172592*CI_tZ
 + 0.0023794638657654243*C_lqM + -0.0014188107569577735*C_eq
 + -0.0022971815962022214*C_lu + 0.00048047355137390913*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_P80_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        
        double opt_obs_clic_3000_0_P80_7 = 1.6429280154454258*C_phiQm + -23.418581880871088*C_phit
 + -1060.972953308112*C_tW + 1129.5794504410326*C_tZ
 + -0.002502097678094232*CI_tW + 0.0023794638657654247*CI_tZ
 + 11710.576087728981*C_lqM + 9542.496636605832*C_eq
 + 1625.3092019318083*C_lu + 56033.1732001331*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_P80_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_clic_3000_0_P80_8 = -80.43141297339808*C_phiQm + -3.1112272836851673*C_phit
 + -4343.5453844831945*C_tW + 4600.191855822186*C_tZ
 + 0.0005387824547244251*CI_tW + -0.0014188107569577735*CI_tZ
 + 9542.496636605832*C_lqM + 196966.39327873895*C_eq
 + 36403.12574542785*C_lu + 48675.23130193319*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_P80_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_clic_3000_0_P80_9 = -14.933801174827105*C_phiQm + -0.2047463414229398*C_phit
 + -895.1839379354279*C_tW + 936.1159110158246*C_tZ
 + 0.0024130670808982227*CI_tW + -0.0022971815962022214*CI_tZ
 + 1625.3092019318083*C_lqM + 36403.12574542785*C_eq
 + 10114.537052966794*C_lu + 11999.145597546125*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_P80_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_clic_3000_0_P80_10 = 6.559065280624168*C_phiQm + -168.78343584875955*C_phit
 + -6866.053235950871*C_tW + 7313.668222089808*C_tZ
 + 0.0005287454515653952*CI_tW + 0.0004804735513739091*CI_tZ
 + 56033.1732001331*C_lqM + 48675.23130193321*C_eq
 + 11999.145597546125*C_lu + 404167.6611825897*C_eu;

        
        return lumi*eff*opt_obs_clic_3000_0_P80_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_clic_3000_0_P80, it goes from 1 to 10 in consecutive pairs.\n");
    }

}






opt_obs_fcc_350::opt_obs_fcc_350(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_fcc_350" << "eff_fcc_350");               
}

double opt_obs_fcc_350::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_fcc_350");
    double eff = SM.getOptionalParameter("eff_fcc_350");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    
    
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_fcc_350_1 = 0.14011277771651373*C_phiQm + 0.11431227773816992*C_phit
 + -12.185763553299475*C_tW + 9.288530644404148*C_tZ
 + -6.301835467966674e-07*CI_tW + 8.784285905278263e-07*CI_tZ
 + 6.184787905699479*C_lqM + 2.931086905517934*C_eq
 + 5.75905215840804*C_lu + 3.2039699273299567*C_eu;

        return lumi*eff*opt_obs_fcc_350_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_fcc_350_2 = 0.1143122777381699*C_phiQm + 0.1367395724202781*C_phit
 + -11.705408832184528*C_tW + 8.87988728268522*C_tZ
 + -3.9672492297987094e-07*CI_tW + 6.84935576119742e-07*CI_tZ
 + 5.595438527411547*C_lqM + 3.017496409752115*C_eq
 + 5.962696662165202*C_lu + 2.7770546798487117*C_eu;

        
        return lumi*eff*opt_obs_fcc_350_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_fcc_350_3 = -12.185763553299477*C_phiQm + -11.705408832184528*C_phit
 + 1633.3032632995635*C_tW + -1307.3223957864402*C_tZ
 + 1.8864190160300214e-15*CI_tW + -1.4944179914089339e-15*CI_tZ
 + -736.0393300979259*C_lqM + -480.70029686695653*C_eq
 + -729.967393550436*C_lu + -487.8938208019186*C_eu;


        return lumi*eff*opt_obs_fcc_350_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_fcc_350_4 = 9.288530644404148*C_phiQm + 8.879887282685218*C_phit
 + -1307.3223957864402*C_tW + 1053.9960263492844*C_tZ
 + 0.0001524403949571921*CI_tW + -0.00013437451691864874*CI_tZ
 + 583.6675446449821*C_lqM + 392.19020948005766*C_eq
 + 578.6225727646719*C_lu + 398.4466256008428*C_eu;

        
        return lumi*eff*opt_obs_fcc_350_4;

        
    } else if(b_min == 5 && b_max == 5){
                
        double opt_obs_fcc_350_5 = -6.301835467966675e-07*C_phiQm + -3.96724922979871e-07*C_phit
 + 1.8864190160300214e-15*C_tW + 0.0001524403949571921*C_tZ
 + 0.8726123248001882*CI_tW + -0.698618345434601*CI_tZ
 + 4.1951726349847224e-05*C_lqM + 1.0624318346612374e-07*C_eq
 + -4.195172635131853e-05*C_lu + -1.0624318334552496e-07*C_eu;

        
        return lumi*eff*opt_obs_fcc_350_5;

        
    } else if(b_min == 6 && b_max == 6){
        

        double opt_obs_fcc_350_6 = 8.784285905278264e-07*C_phiQm + 6.849355761197421e-07*C_phit
 + -1.4944179914089339e-15*C_tW + -0.00013437451691864874*C_tZ
 + -0.6986183454346009*CI_tW + 0.5601153022131046*CI_tZ
 + 0.0007083607546546914*C_lqM + 6.902259140289349e-06*C_eq
 + 0.0007907258901739203*C_lu + 7.893951274792737e-06*C_eu;

        
        return lumi*eff*opt_obs_fcc_350_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        
        double opt_obs_fcc_350_7 = 6.184787905699479*C_phiQm + 5.595438527411546*C_phit
 + -736.0393300979259*C_tW + 583.6675446449822*C_tZ
 + 4.1951726349847224e-05*CI_tW + 0.0007083607546546913*CI_tZ
 + 342.1073510838338*C_lqM + 208.12251497925953*C_eq
 + 332.53911409852356*C_lu + 214.53504825609784*C_eu;

        
        return lumi*eff*opt_obs_fcc_350_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_fcc_350_8 = 2.9310869055179345*C_phiQm + 3.017496409752115*C_phit
 + -480.7002968669566*C_tW + 392.19020948005766*C_tZ
 + 1.0624318346612374e-07*CI_tW + 6.902259140289349e-06*CI_tZ
 + 208.12251497925953*C_lqM + 151.04132911264546*C_eq
 + 209.727403326688*C_lu + 150.33102927041702*C_eu;

        
        return lumi*eff*opt_obs_fcc_350_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_fcc_350_9 = 5.75905215840804*C_phiQm + 5.962696662165201*C_phit
 + -729.967393550436*C_tW + 578.6225727646719*C_tZ
 + -4.195172635131853e-05*CI_tW + 0.0007907258901739204*CI_tZ
 + 332.53911409852356*C_lqM + 209.727403326688*C_eq
 + 337.03832768284326*C_lu + 208.87182725637823*C_eu;

        
        return lumi*eff*opt_obs_fcc_350_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_fcc_350_10 = 3.203969927329957*C_phiQm + 2.7770546798487117*C_phit
 + -487.89382080191854*C_tW + 398.4466256008428*C_tZ
 + -1.0624318334552494e-07*CI_tW + 7.893951274792737e-06*CI_tZ
 + 214.53504825609784*C_lqM + 150.33102927041702*C_eq
 + 208.87182725637823*C_lu + 156.42028181004318*C_eu;

        
        return lumi*eff*opt_obs_fcc_350_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_fcc_350, it goes from 1 to 10 in consecutive pairs.\n");
    }

}







opt_obs_fcc_365::opt_obs_fcc_365(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_fcc_365" << "eff_fcc_365");               
}

double opt_obs_fcc_365::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    double lumi = SM.getOptionalParameter("lumi_fcc_365");
    double eff = SM.getOptionalParameter("eff_fcc_365");
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    
    
            
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_fcc_365_1 = 0.3096101166222116*C_phiQm + 0.1374996869775552*C_phit
 + -22.582867702025293*C_tW + 17.312761125744444*C_tZ
 + -3.1100818808860036e-07*CI_tW + 7.94747946323648e-07*CI_tZ
 + 13.215301662149164*C_lqM + 5.108151080225039*C_eq
 + 10.109979583497916*C_lu + 7.100868603670604*C_eu;

        return lumi*eff*opt_obs_fcc_365_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_fcc_365_2 = 0.1374996869775552*C_phiQm + 0.28831813648246973*C_phit
 + -19.388357633528106*C_tW + 14.594320806408248*C_tZ
 + -2.992831646787645e-07*CI_tW + -5.632834217509583e-07*CI_tZ
 + 8.94409029528181*C_lqM + 5.772549777343706*C_eq
 + 11.64370528931688*C_lu + 4.001592478968842*C_eu;

        
        return lumi*eff*opt_obs_fcc_365_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_fcc_365_3 = -22.582867702025293*C_phiQm + -19.388357633528106*C_phit
 + 2968.9335825898656*C_tW + -2377.824339806803*C_tZ
 + 0.0*CI_tW + 5.434739628219844e-16*CI_tZ
 + -1431.744492620986*C_lqM + -905.6658689043426*C_eq
 + -1387.793890979745*C_lu + -958.2560508354894*C_eu;


        return lumi*eff*opt_obs_fcc_365_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_fcc_365_4 = 17.312761125744444*C_phiQm + 14.59432080640825*C_phit
 + -2377.824339806803*C_tW + 1917.8300139679118*C_tZ
 + -1.6058961970037645e-05*CI_tW + 8.85549201789431e-05*CI_tZ
 + 1136.4088048979982*C_lqM + 738.5098104772077*C_eq
 + 1099.8616487812044*C_lu + 784.2347754227291*C_eu;

        
        return lumi*eff*opt_obs_fcc_365_4;

        
    } else if(b_min == 5 && b_max == 5){

        double opt_obs_fcc_365_5 = -3.110081880886003e-07*C_phiQm + -2.992831646787645e-07*C_phit
 + 0.0*C_tW + -1.6058961970037645e-05*C_tZ
 + 6.494043805075946*CI_tW + -5.203461449601805*CI_tZ
 + 2.7405901287701535e-06*C_lqM + -1.348991197969875e-07*C_eq
 + -2.7405901287701535e-06*C_lu + 1.348991197969875e-07*C_eu;

        
        return lumi*eff*opt_obs_fcc_365_5;

        
    } else if(b_min == 6 && b_max == 6){
        
    
        double opt_obs_fcc_365_6 = 7.947479463236479e-07*C_phiQm + -5.632834217509582e-07*C_phit
 + 5.434739628219844e-16*C_tW + 8.855492017894309e-05*C_tZ
 + -5.203461449601805*CI_tW + 4.175591519765497*CI_tZ
 + 0.002724351057717014*C_lqM + 2.9369050480910674e-06*C_eq
 + 0.002482700959665445*C_lu + -2.9432552978997855e-06*C_eu;

        
        return lumi*eff*opt_obs_fcc_365_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        
        double opt_obs_fcc_365_7 = 13.215301662149166*C_phiQm + 8.94409029528181*C_phit
 + -1431.7444926209862*C_tW + 1136.4088048979982*C_tZ
 + 2.7405901287701535e-06*CI_tW + 0.0027243510577170138*CI_tZ
 + 732.6709969552045*C_lqM + 410.14022323814953*C_eq
 + 656.7345420730154*C_lu + 460.8862102342187*C_eu;

        
        return lumi*eff*opt_obs_fcc_365_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_fcc_365_8 = 5.1081510802250385*C_phiQm + 5.772549777343706*C_phit
 + -905.6658689043426*C_tW + 738.5098104772076*C_tZ
 + -1.348991197969875e-07*CI_tW + 2.9369050480910678e-06*CI_tZ
 + 410.14022323814953*C_lqM + 303.1792753366912*C_eq
 + 423.41835775708824*C_lu + 296.9563045900016*C_eu;

        
        return lumi*eff*opt_obs_fcc_365_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_fcc_365_9 = 10.109979583497916*C_phiQm + 11.643705289316882*C_phit
 + -1387.793890979745*C_tW + 1099.8616487812042*C_tZ
 + -2.7405901287701535e-06*CI_tW + 0.002482700959665445*CI_tZ
 + 656.7345420730154*C_lqM + 423.41835775708824*C_eq
 + 693.4202984063455*C_lu + 415.9272934918138*C_eu;

        
        return lumi*eff*opt_obs_fcc_365_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_fcc_365_10 = 7.100868603670604*C_phiQm + 4.001592478968843*C_phit
 + -958.2560508354894*C_tW + 784.2347754227291*C_tZ
 + 1.348991197969875e-07*CI_tW + -2.9432552978997855e-06*CI_tZ
 + 460.88621023421877*C_lqM + 296.9563045900016*C_eq
 + 415.92729349181377*C_lu + 345.34132925320677*C_eu;

        
        return lumi*eff*opt_obs_fcc_365_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_fcc_365, it goes from 1 to 10 in consecutive pairs.\n");
    }

}










opt_obs_cepc_350::opt_obs_cepc_350(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_cepc_350" << "eff_cepc_350" );               
}

double opt_obs_cepc_350::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    double lumi = SM.getOptionalParameter("lumi_cepc_350");
    double eff = SM.getOptionalParameter("eff_cepc_350");
    
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    
    
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_cepc_350_1 = 0.14011277771651373*C_phiQm + 0.11431227773816992*C_phit + 
        -12.185763553299475*C_tW + 9.288530644404148*C_tZ + -6.301835467966674e-07*CI_tW + 
        8.784285905278263e-07*CI_tZ + 6.184787905699479*C_lqM + 2.931086905517934*C_eq + 
        5.75905215840804*C_lu + 3.2039699273299567*C_eu;

        return lumi*eff*opt_obs_cepc_350_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_cepc_350_2 = 0.1143122777381699*C_phiQm + 0.1367395724202781*C_phit + 
        -11.705408832184528*C_tW + 8.87988728268522*C_tZ + -3.9672492297987094e-07*CI_tW + 
        6.84935576119742e-07*CI_tZ + 5.595438527411547*C_lqM + 3.017496409752115*C_eq + 
        5.962696662165202*C_lu + 2.7770546798487117*C_eu;

        
        return lumi*eff*opt_obs_cepc_350_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_cepc_350_3 = -12.185763553299477*C_phiQm + -11.705408832184528*C_phit + 
        1633.3032632995635*C_tW + -1307.3223957864402*C_tZ + 1.8864190160300214e-15*CI_tW + 
        -1.4944179914089339e-15*CI_tZ + -736.0393300979259*C_lqM + -480.70029686695653*C_eq + 
        -729.967393550436*C_lu + -487.8938208019186*C_eu;


        return lumi*eff*opt_obs_cepc_350_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_cepc_350_4 = 9.288530644404148*C_phiQm + 8.879887282685218*C_phit + 
        -1307.3223957864402*C_tW + 1053.9960263492844*C_tZ + 0.0001524403949571921*CI_tW +
        -0.00013437451691864874*CI_tZ + 583.6675446449821*C_lqM + 392.19020948005766*C_eq +
        578.6225727646719*C_lu + 398.4466256008428*C_eu;

        
        return lumi*eff*opt_obs_cepc_350_4;

        
    } else if(b_min == 5 && b_max == 5){
                
        double opt_obs_cepc_350_5 = -6.301835467966675e-07*C_phiQm + -3.96724922979871e-07*C_phit + 
        1.8864190160300214e-15*C_tW + 0.0001524403949571921*C_tZ + 0.8726123248001882*CI_tW + 
        -0.698618345434601*CI_tZ + 4.1951726349847224e-05*C_lqM + 1.0624318346612374e-07*C_eq + 
        -4.195172635131853e-05*C_lu + -1.0624318334552496e-07*C_eu;

        
        return lumi*eff*opt_obs_cepc_350_5;

        
    } else if(b_min == 6 && b_max == 6){

        double opt_obs_cepc_350_6 = 8.784285905278264e-07*C_phiQm + 6.849355761197421e-07*C_phit +
        -1.4944179914089339e-15*C_tW + -0.00013437451691864874*C_tZ + -0.6986183454346009*CI_tW + 
        0.5601153022131046*CI_tZ + 0.0007083607546546914*C_lqM + 6.902259140289349e-06*C_eq + 
        0.0007907258901739203*C_lu + 7.893951274792737e-06*C_eu;

        
        return lumi*eff*opt_obs_cepc_350_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        
        double opt_obs_cepc_350_7 = 6.184787905699479*C_phiQm + 5.595438527411546*C_phit + 
        -736.0393300979259*C_tW + 583.6675446449822*C_tZ + 4.1951726349847224e-05*CI_tW + 
        0.0007083607546546913*CI_tZ + 342.1073510838338*C_lqM + 208.12251497925953*C_eq + 
        332.53911409852356*C_lu + 214.53504825609784*C_eu;

        
        return lumi*eff*opt_obs_cepc_350_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_cepc_350_8 = 2.9310869055179345*C_phiQm + 3.017496409752115*C_phit + 
        -480.7002968669566*C_tW + 392.19020948005766*C_tZ + 1.0624318346612374e-07*CI_tW + 
        6.902259140289349e-06*CI_tZ + 208.12251497925953*C_lqM + 151.04132911264546*C_eq + 
        209.727403326688*C_lu + 150.33102927041702*C_eu;

        
        return lumi*eff*opt_obs_cepc_350_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_cepc_350_9 = 5.75905215840804*C_phiQm + 5.962696662165201*C_phit + 
        -729.967393550436*C_tW + 578.6225727646719*C_tZ + -4.195172635131853e-05*CI_tW + 
        0.0007907258901739204*CI_tZ + 332.53911409852356*C_lqM + 209.727403326688*C_eq + 
        337.03832768284326*C_lu + 208.87182725637823*C_eu;

        
        return lumi*eff*opt_obs_cepc_350_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_cepc_350_10 = 3.203969927329957*C_phiQm + 2.7770546798487117*C_phit + 
        -487.89382080191854*C_tW + 398.4466256008428*C_tZ + -1.0624318334552494e-07*CI_tW + 
        7.893951274792737e-06*CI_tZ + 214.53504825609784*C_lqM + 150.33102927041702*C_eq + 
        208.87182725637823*C_lu + 156.42028181004318*C_eu;

        
        return lumi*eff*opt_obs_cepc_350_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_cepc_350, it goes from 1 to 10 in consecutive pairs.\n");
    }

}








opt_obs_cepc_360::opt_obs_cepc_360(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_cepc_360" << "eff_cepc_360" );              
}

double opt_obs_cepc_360::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    double lumi = SM.getOptionalParameter("lumi_cepc_360");
    double eff = SM.getOptionalParameter("eff_cepc_360");
    
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_cepc_360_1 = 0.2613241621954116*C_phiQm + 0.1426733201447283*C_phit +
        -20.064010398417633*C_tW + 15.354456154006607*C_tZ + 2.2233893311168813e-06*CI_tW +
        -2.17343351945387e-06*CI_tZ + 11.220847604453926*C_lqM + 4.6311827882863135*C_eq + 
        9.141898861570013*C_lu + 5.964710688938101*C_eu;

        return lumi*eff*opt_obs_cepc_360_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_cepc_360_2 = 0.14267332014472828*C_phiQm + 0.2463791974661928*C_phit + 
        -17.85939561371129*C_tW + 13.478652708822949*C_tZ + 2.2321190489044256e-06*CI_tW + 
        -2.0736169422038726e-06*CI_tZ + 8.355437781487353*C_lqM + 5.068584956663598*C_eq + 
        10.158229783498602*C_lu + 3.886657692490702*C_eu;

        
        return lumi*eff*opt_obs_cepc_360_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_cepc_360_3 = -20.06401039841763*C_phiQm + -17.859395613711285*C_phit +
        2653.2306615355387*C_tW + -2124.576155828475*C_tZ + 3.9195595110845996e-15*CI_tW + 
        -3.481369316155128e-15*CI_tZ + -1251.3127887303644*C_lqM + -799.7349591134495*C_eq + 
        -1221.812197269644*C_lu + -834.9182338644936*C_eu;


        return lumi*eff*opt_obs_cepc_360_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_cepc_360_4 = 15.35445615400661*C_phiQm + 13.478652708822949*C_phit + 
        -2124.576155828475*C_tW + 1713.3646136694374*C_tZ + -4.2083644579048585e-08*CI_tW + 
        1.2931486456754067e-05*CI_tZ + 992.9035403701065*C_lqM + 652.2447902353307*C_eq + 
        968.3786440723969*C_lu + 682.8383138416785*C_eu;

        
        return lumi*eff*opt_obs_cepc_360_4;

        
    } else if(b_min == 5 && b_max == 5){
                

        double opt_obs_cepc_360_5 = 2.2233893311168813e-06*C_phiQm + 2.232119048904426e-06*C_phit +
        3.9195595110845996e-15*C_tW + -4.208364457904858e-08*C_tZ + 4.320076095510152*CI_tW + 
        -3.4606905205694356*CI_tZ + -1.7745027352550902e-05*C_lqM + 4.609720390592397e-05*C_eq + 
        1.7745027346969064e-05*C_lu + -4.609720390592397e-05*C_eu;

        
        return lumi*eff*opt_obs_cepc_360_5;

        
    } else if(b_min == 6 && b_max == 6){
        

        double opt_obs_cepc_360_6 = -2.17343351945387e-06*C_phiQm + -2.0736169422038726e-06*C_phit +
        -3.4813693161551284e-15*C_tW + 1.2931486456754065e-05*C_tZ + -3.460690520569435*CI_tW + 
        2.7763522906532816*CI_tZ + -2.2137317576642757e-05*C_lqM + -3.473687864018645e-05*C_eq +
        1.5433020527597982e-05*C_lu + 4.7339506921695865e-05*C_eu;

        
        return lumi*eff*opt_obs_cepc_360_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        
        double opt_obs_cepc_360_7 = 11.220847604453926*C_phiQm + 8.35543778148735*C_phit + 
        -1251.3127887303644*C_tW + 992.9035403701064*C_tZ + -1.7745027352550902e-05*CI_tW + 
        -2.2137317576642753e-05*CI_tZ + 620.4800958800838*C_lqM + 356.8696842997229*C_eq +
        571.0327724522122*C_lu + 389.9461374063771*C_eu;

        
        return lumi*eff*opt_obs_cepc_360_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_cepc_360_8 = 4.631182788286314*C_phiQm + 5.0685849566635985*C_phit + 
        -799.7349591134495*C_tW + 652.2447902353307*C_tZ + 4.609720390592397e-05*CI_tW + 
        -3.473687864018645e-05*CI_tZ + 356.86968429972285*C_lqM + 262.13125655995003*C_eq + 
        365.4021218466032*C_lu + 258.202217005985*C_eu;

        
        return lumi*eff*opt_obs_cepc_360_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_cepc_360_9 = 9.141898861570015*C_phiQm + 10.158229783498603*C_phit + 
        -1221.812197269644*C_tW + 968.3786440723967*C_tZ + 1.7745027346969064e-05*CI_tW + 
        1.5433020527597982e-05*CI_tZ + 571.032772452212*C_lqM + 365.4021218466032*C_eq + 
        594.7145421768979*C_lu + 360.67311623081036*C_eu;

        
        return lumi*eff*opt_obs_cepc_360_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_cepc_360_10 = 5.964710688938101*C_phiQm + 3.8866576924907026*C_phit + 
        -834.9182338644936*C_tW + 682.8383138416785*C_tZ + -4.609720390592397e-05*CI_tW +
        4.7339506921695865e-05*CI_tZ + 389.94613740637715*C_lqM + 258.20221700598506*C_eq +
        360.67311623081036*C_lu + 289.6974251448306*C_eu;

        
        return lumi*eff*opt_obs_cepc_360_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_cepc_360, it goes from 1 to 10 in consecutive pairs.\n");
    }

}



















opt_obs_muon_3TeV::opt_obs_muon_3TeV(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_muon_3TeV" << "eff_muon_3TeV");             
}


double opt_obs_muon_3TeV::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    double lumi = SM.getOptionalParameter("lumi_muon_3TeV");
    double eff = SM.getOptionalParameter("eff_muon_3TeV");  
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_muon_3TeV_1 = 0.04385217227970124*C_phiQm + -0.028148548666231166*C_phit
 + -1.4092608584531616*C_tW + 1.124019139686304*C_tZ
 + -7.977038258363135e-06*CI_tW + 7.036793152737949e-06*CI_tZ
 + 81.4818764767476*C_lqM + -8.604788220305771*C_eq
 + -11.840306688250017*C_lu + 51.63086027090416*C_eu;

        return lumi*eff*opt_obs_muon_3TeV_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_muon_3TeV_2 = -0.028148548666231166*C_phiQm + 0.04190396734486602*C_phit
 + -0.2780150129561469*C_tW + 0.15115927218233643*C_tZ
 + -5.633978599381841e-06*CI_tW + 4.966064908062675e-06*CI_tZ
 + -35.77560383864257*C_lqM + 24.35667155369132*C_eq
 + 54.59369968751945*C_lu + -34.73898919422445*C_eu;

        
        return lumi*eff*opt_obs_muon_3TeV_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_muon_3TeV_3 = -1.4092608584531616*C_phiQm + -0.27801501295614695*C_phit
 + 441.18242192091446*C_tW + -355.4222818137299*C_tZ
 + 0.0001726210905694075*CI_tW + -0.00015223482257482505*CI_tZ
 + -4479.79183572281*C_lqM + -1844.1840900063494*C_eq
 + -3517.0510674854777*C_lu + -3364.5790285723892*C_eu;


        return lumi*eff*opt_obs_muon_3TeV_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_muon_3TeV_4 = 1.1240191396863042*C_phiQm + 0.15115927218233632*C_phit
 + -355.42228181372985*C_tW + 287.2187291027542*C_tZ
 + -0.00015047760467409963*CI_tW + 0.0001362792618805949*CI_tZ
 + 3590.793894397487*C_lqM + 1491.1074227140994*C_eq
 + 2771.6673406191285*C_lu + 2808.6943190871875*C_eu;

        
        return lumi*eff*opt_obs_muon_3TeV_4;

        
    } else if(b_min == 5 && b_max == 5){

        
        double opt_obs_muon_3TeV_5 = -7.977038258363135e-06*C_phiQm + -5.633978599381841e-06*C_phit
 + 0.0001726210905694075*C_tW + -0.00015047760467409963*C_tZ
 + 260.17152020222505*CI_tW + -210.00727926700463*CI_tZ
 + 0.0015099379091762296*C_lqM + 0.0034543598159384697*C_eq
 + -0.0001589410331599803*C_lu + 0.0031486271865535068*C_eu;

        
        return lumi*eff*opt_obs_muon_3TeV_5;

        
    } else if(b_min == 6 && b_max == 6){
        

        double opt_obs_muon_3TeV_6 = 7.0367931527379485e-06*C_phiQm + 4.9660649080626745e-06*C_phit
 + -0.00015223482257482505*C_tW + 0.0001362792618805949*C_tZ
 + -210.00727926700463*CI_tW + 169.84737755084333*CI_tZ
 + -0.00045823347557478393*C_lqM + -0.005046201139186361*C_eq
 + -0.00024864009377309113*C_lu + -0.0008886617540116975*C_eu;

        
        return lumi*eff*opt_obs_muon_3TeV_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        
        double opt_obs_muon_3TeV_7 = 81.4818764767476*C_phiQm + -35.775603838642574*C_phit
 + -4479.79183572281*C_tW + 3590.793894397487*C_tZ
 + 0.0015099379091762296*CI_tW + -0.00045823347557478393*CI_tZ
 + 180021.00658073605*C_lqM + 16617.773389917675*C_eq
 + 27234.44940953821*C_lu + 113804.03974324113*C_eu;

        
        return lumi*eff*opt_obs_muon_3TeV_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_muon_3TeV_8 = -8.604788220305773*C_phiQm + 24.35667155369132*C_phit
 + -1844.1840900063492*C_tW + 1491.1074227140994*C_tZ
 + 0.0034543598159384697*CI_tW + -0.0050462011391863605*CI_tZ
 + 16617.773389917675*C_lqM + 38846.35764251843*C_eq
 + 58429.53098862293*C_lu + 10227.76492675033*C_eu;

        
        return lumi*eff*opt_obs_muon_3TeV_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_muon_3TeV_9 = -11.84030668825002*C_phiQm + 54.59369968751945*C_phit
 + -3517.0510674854777*C_tW + 2771.667340619128*C_tZ
 + -0.00015894103315998033*CI_tW + -0.000248640093773091*CI_tZ
 + 27234.44940953821*C_lqM + 58429.53098862294*C_eq
 + 132087.34957660464*C_lu + 24206.416428780652*C_eu;

        
        return lumi*eff*opt_obs_muon_3TeV_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_muon_3TeV_10 = 51.63086027090416*C_phiQm + -34.73898919422445*C_phit
 + -3364.5790285723892*C_tW + 2808.694319087188*C_tZ
 + 0.0031486271865535063*CI_tW + -0.0008886617540116973*CI_tZ
 + 113804.03974324113*C_lqM + 10227.76492675033*C_eq
 + 24206.416428780656*C_lu + 107954.38625978182*C_eu;

        
        return lumi*eff*opt_obs_muon_3TeV_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_muon_3TeV, it goes from 1 to 10 in consecutive pairs.\n");
    }

}









opt_obs_muon_10TeV::opt_obs_muon_10TeV(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_muon_10TeV" << "eff_muon_10TeV");             
}

double opt_obs_muon_10TeV::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    double lumi = SM.getOptionalParameter("lumi_muon_10TeV");
    double eff = SM.getOptionalParameter("eff_muon_10TeV");  
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
    
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_muon_10TeV_1 = 0.004009591465404328*C_phiQm + -0.0025955162780930347*C_phit
 + -0.1279413945534101*C_tW + 0.10208027766790802*C_tZ
 + -3.126517690344216e-08*CI_tW + 2.8042751565421526e-08*CI_tZ
 + 82.56159769990425*C_lqM + -9.078661839817244*C_eq
 + -12.645295235141134*C_lu + 52.36745183601147*C_eu;

        return lumi*eff*opt_obs_muon_10TeV_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_muon_10TeV_2 = -0.0025955162780930347*C_phiQm + 0.003842251716906442*C_phit
 + -0.024786515652280197*C_tW + 0.013334592054865113*C_tZ
 + 1.4031866438365768e-08*CI_tW + -5.130136388713993e-08*CI_tZ
 + -36.87550634000039*C_lqM + 24.7558327927607*C_eq
 + 55.49642160324448*C_lu + -35.61466057807351*C_eu;

        
        return lumi*eff*opt_obs_muon_10TeV_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_muon_10TeV_3 = -0.1279413945534101*C_phiQm + -0.024786515652280204*C_phit
 + 283.25824765621496*C_tW + -228.5965992035851*C_tZ
 + 3.253024792087326e-14*CI_tW + 2.9637807221894017e-06*CI_tZ
 + -4506.791869657863*C_lqM + -1843.3403169753803*C_eq
 + -3535.242771315796*C_lu + -3391.2047870117485*C_eu;


        return lumi*eff*opt_obs_muon_10TeV_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_muon_10TeV_4 =  0.10208027766790802*C_phiQm + 0.013334592054865119*C_phit
 + -228.5965992035851*C_tW + 184.8744569891588*C_tZ
 + -2.8688478935464505e-14*CI_tW + -1.890427431205348e-06*CI_tZ
 + 3612.9242583598993*C_lqM + 1490.2041522727304*C_eq
 + 2785.6031058809654*C_lu + 2831.5724126148025*C_eu;

        
        return lumi*eff*opt_obs_muon_10TeV_4;

        
    } else if(b_min == 5 && b_max == 5){
                

        double opt_obs_muon_10TeV_5 = -3.126517690344216e-08*C_phiQm + 1.403186643836577e-08*C_phit
 + 3.253024792087326e-14*C_tW + -2.8688478935464505e-14*C_tZ
 + 266.58240013802356*CI_tW + -215.1848633518403*CI_tZ
 + 0.001861721329438426*C_lqM + 0.001049389567764944*C_eq
 + -0.000773498529909885*C_lu + 0.0018733695680428456*C_eu;

        
        return lumi*eff*opt_obs_muon_10TeV_5;

        
    } else if(b_min == 6 && b_max == 6){
        

        double opt_obs_muon_10TeV_6 = 2.8042751565421523e-08*C_phiQm + -5.130136388713993e-08*C_phit
 + 2.963780722189402e-06*C_tW + -1.8904274312053481e-06*C_tZ
 + -215.18486335184028*CI_tW + 174.03687752436392*CI_tZ
 + -0.0022641706423453044*C_lqM + -0.001851416742736093*C_eq
 + -0.0006531182581423842*C_lu + -0.0018750003691621938*C_eu;

        
        return lumi*eff*opt_obs_muon_10TeV_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        
        double opt_obs_muon_10TeV_7 = 82.56159769990425*C_phiQm + -36.8755063400004*C_phit
 + -4506.791869657863*C_tW + 3612.9242583598993*C_tZ
 + 0.0018617213294384262*CI_tW + -0.0022641706423453044*CI_tZ
 + 2017659.163435116*C_lqM + 174935.5369042937*C_eq
 + 286344.173467685*C_lu + 1275060.3700656057*C_eu;

        
        return lumi*eff*opt_obs_muon_10TeV_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_muon_10TeV_8 = -9.078661839817244*C_phiQm + 24.7558327927607*C_phit
 + -1843.3403169753803*C_tW + 1490.2041522727307*C_tZ
 + 0.001049389567764944*CI_tW + -0.001851416742736093*CI_tZ
 + 174935.5369042937*C_lqM + 432924.2401327633*C_eq
 + 651637.3556322666*C_lu + 105591.6078334171*C_eu;

        
        return lumi*eff*opt_obs_muon_10TeV_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_muon_10TeV_9 = -12.64529523514113*C_phiQm + 55.49642160324447*C_phit
 + -3535.242771315796*C_tW + 2785.603105880965*C_tZ
 + -0.0007734985299098851*CI_tW + -0.0006531182581423844*CI_tZ
 + 286344.173467685*C_lqM + 651637.3556322666*C_eq
 + 1481331.0071603623*C_lu + 260167.47923858245*C_eu;

        
        return lumi*eff*opt_obs_muon_10TeV_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_muon_10TeV_10 = 52.367451836011476*C_phiQm + -35.61466057807352*C_phit
 + -3391.2047870117485*C_tW + 2831.572412614803*C_tZ
 + 0.0018733695680428456*CI_tW + -0.0018750003691621938*CI_tZ
 + 1275060.3700656057*C_lqM + 105591.60783341713*C_eq
 + 260167.47923858243*C_lu + 1212777.6020041464*C_eu;

        
        return lumi*eff*opt_obs_muon_10TeV_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_muon_10TeV, it goes from 1 to 10 in consecutive pairs.\n");
    }

}









opt_obs_muon_30TeV::opt_obs_muon_30TeV(const StandardModel& SM_i)
: ThObservable(SM_i), myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "lumi_muon_30TeV" << "eff_muon_30TeV");             
}

double opt_obs_muon_30TeV::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    double lumi = SM.getOptionalParameter("lumi_muon_30TeV");
    double eff = SM.getOptionalParameter("eff_muon_30TeV");  
    
    //bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    //bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //['ClqA', 'CeqA', 'CpqA', 'ClqV', 'CeqV', 'CpqV', 'CuZR', 'CuAR', 'CuZI', 'CuAI']
    //['cpQM','cpt','ctW','ctZ','ctWI','ctZI','cQlM1','cQe1','ctl1','cte1']
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double CI_tW = 0.;
    double CI_tZ = 0.;
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    
       
    if(b_min == 1 && b_max == 1){
        
        double opt_obs_muon_30TeV_1 = 0.0004461256530174988*C_phiQm + -0.00028900246330781226*C_phit
 + -0.014228308284791904*C_tW + 0.011352784758373114*C_tZ
 + -1.1816937810350695e-09*CI_tW + 3.854912971639567e-08*CI_tZ
 + 82.65934156076267*C_lqM + -9.121602198856497*C_eq
 + -12.718396391937766*C_lu + 52.433273612505545*C_eu;

        return lumi*eff*opt_obs_muon_30TeV_1;
        
        
    } else if(b_min == 2 && b_max == 2){
        
        double opt_obs_muon_30TeV_2 = -0.00028900246330781226*C_phiQm + 0.0004276109473491372*C_phit
 + -0.002750126494113585*C_tW + 0.0014775347515775697*C_tZ
 + 5.040646735668308e-10*CI_tW + -3.859728099488311e-08*CI_tZ
 + -36.97554779688711*C_lqM + 24.79165483170672*C_eq
 + 55.57774124725651*C_lu + -35.69229666080014*C_eu;

        
        return lumi*eff*opt_obs_muon_30TeV_2;

        
    } else if(b_min == 3 && b_max == 3){
        
        double opt_obs_muon_30TeV_3 = -0.014228308284791906*C_phiQm + -0.002750126494113587*C_phit
 + 269.19388278260504*C_tW + -217.3029378226469*C_tZ
 + 3.185952115961814e-15*CI_tW + 2.1050341735418113e-06*CI_tZ
 + -4509.21767007115*C_lqM + -1843.434786215226*C_eq
 + -3536.991405756283*C_lu + -3393.278552271313*C_eu;


        return lumi*eff*opt_obs_muon_30TeV_3;

    } else if(b_min == 4 && b_max == 4){
        
        double opt_obs_muon_30TeV_4 = 0.011352784758373115*C_phiQm + 0.0014775347515775684*C_phit
 + -217.3029378226469*C_tW + 175.76231573209878*C_tZ
 + -2.80969639058673e-15*CI_tW + -1.8118717493423767e-06*CI_tZ
 + 3614.8998799890082*C_lqM + 1490.257010091355*C_eq
 + 2786.9374266503196*C_lu + 2833.3470680536734*C_eu;

        
        return lumi*eff*opt_obs_muon_30TeV_4;

        
    } else if(b_min == 5 && b_max == 5){
                

        double opt_obs_muon_30TeV_5 = -1.1816937810350695e-09*C_phiQm + 5.040646735668308e-10*C_phit
 + 3.185952115961814e-15*C_tW + -2.8096963905867296e-15*C_tZ
 + 267.1652755741726*CI_tW + -215.65650170132454*CI_tZ
 + 0.0015534013544804435*C_lqM + -6.681327479592764e-05*C_eq
 + 0.0006738972790641265*C_lu + 6.635383658870812e-05*C_eu;

        
        return lumi*eff*opt_obs_muon_30TeV_5;

        
    } else if(b_min == 6 && b_max == 6){
        

        double opt_obs_muon_30TeV_6 = 3.854912971639568e-08*C_phiQm + -3.859728099488311e-08*C_phit
 + 2.1050341735418113e-06*C_tW + -1.8118717493423765e-06*C_tZ
 + -215.65650170132457*CI_tW + 174.41927946458918*CI_tZ
 + 0.011139318264304546*C_lqM + -0.0005712903651506302*C_eq
 + 0.014551832618290027*C_lu + -0.0005135301116716116*C_eu;

        
        return lumi*eff*opt_obs_muon_30TeV_6;

        
    } else if(b_min == 7 && b_max == 7){
        
        double opt_obs_muon_30TeV_7 = 82.65934156076267*C_phiQm + -36.97554779688711*C_phit
 + -4509.21767007115*C_tW + 3614.8998799890082*C_tZ
 + 0.0015534013544804435*CI_tW + 0.011139318264304546*CI_tZ
 + 18172430.059509046*C_lqM + 1566603.2045368692*C_eq
 + 2564049.7374233436*C_lu + 11483908.709193476*C_eu;

        
        return lumi*eff*opt_obs_muon_30TeV_7;

        
    } else if(b_min == 8 && b_max == 8){
        
        double opt_obs_muon_30TeV_8 = -9.121602198856497*C_phiQm + 24.79165483170672*C_phit
 + -1843.4347862152263*C_tW + 1490.2570100913547*C_tZ
 + -6.681327479592764e-05*CI_tW + -0.0005712903651506302*CI_tZ
 + 1566603.2045368692*C_lqM + 3897392.3117358787*C_eq
 + 5866674.762341308*C_lu + 943905.3767505556*C_eu;

        
        return lumi*eff*opt_obs_muon_30TeV_8;

        
    } else if(b_min == 9 && b_max == 9){
        
        double opt_obs_muon_30TeV_9 = -12.718396391937766*C_phiQm + 55.577741247256505*C_phit
 + -3536.9914057562823*C_tW + 2786.937426650319*C_tZ
 + 0.0006738972790641266*CI_tW + 0.014551832618290027*CI_tZ
 + 2564049.7374233445*C_lqM + 5866674.762341308*C_eq
 + 13342993.830230948*C_lu + 2334423.9987671273*C_eu;

        
        return lumi*eff*opt_obs_muon_30TeV_9;

        
    } else if(b_min == 10 && b_max == 10){
        
        
        double opt_obs_muon_30TeV_10 = 52.433273612505545*C_phiQm + -35.69229666080014*C_phit
 + -3393.2785522713134*C_tW + 2833.347068053673*C_tZ
 + 6.63538365887081e-05*CI_tW + -0.0005135301116716115*CI_tZ
 + 11483908.709193476*C_lqM + 943905.3767505558*C_eq
 + 2334423.9987671273*C_lu + 10925596.17041314*C_eu;

        
        return lumi*eff*opt_obs_muon_30TeV_10;

        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct opt_obs_muon_30TeV, it goes from 1 to 10 in consecutive pairs.\n");
    }

}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////// Observables for CP-violation at tth and thj (prospects and proposed asymmetries) ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



b4_ttH_LO::b4_ttH_LO(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_b4_ttH_LO_bin_m1_m0p5" << "SM_b4_ttH_LO_bin_m0p5_0"
            << "SM_b4_ttH_LO_bin_0_0p5" << "SM_b4_ttH_LO_bin_0p5_1");
    
}

double b4_ttH_LO::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiG();
    double cHGtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiGtil();
    //double cHW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiW();
    //double cHWtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiWtil();
    double ctHRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double ctHIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphiIm();
    double ctGRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double ctGIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tGIm();
    double ctWRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double ctWIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tWIm();
    double ctBRe = 0.;
    double ctBIm = 0.;
    
        
    
    
    //Expressions are written in the basis of SMEFTsim
    //let's change to the basis of dim6top (same as
    //SMEFT@NLO except for the gs on CtG). These
    //models define the covariance derivative 
    //with a different convention
    
    if(flag_LHC_WG_Basis){
        double ctZRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
        double ctZIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZIm();
        double sw2 = 0.22305;
        double cw2 = 1 - sw2;
        double tw2 = sw2/cw2;
        double sw  = sqrt(sw2);
        double tw  = sqrt(tw2);
        ctBRe = ctZRe/sw - ctWRe/tw;
        ctBIm = ctZIm/sw - ctWIm/tw;
        ctWRe = -ctWRe;
        ctWIm = -ctWIm;
        ctGRe = -ctGRe;
        ctGIm = -ctGIm;
    }
    else{
        ctBRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        ctBIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tBIm();
    }
    
    
   
    if(b_min == -1 && b_max == -0.5){
        
        double SM_b4_tth_bin_m1_m0p5 = SM.getOptionalParameter("SM_b4_ttH_LO_bin_m1_m0p5");
        double b4_tth_bin_m1_m0p5_madgraph_LO = 1.0632053747005017;//pb
        double total;
        

        
        if(flag_Quadratic){
                total =  SM_b4_tth_bin_m1_m0p5 + (
                        0.910009531724715 * cHG + -0.0010318009149136 * cHGtil
                    + -0.0067981802286535 * ctGIm + -1.118548376749951 * ctGRe
                    + 0.0001145055069108 * ctHIm + -0.1309188173876159 * ctHRe
                    + -0.0020499059754391 * ctWIm + -0.0013815505152491 * ctWRe
                    + 0.8719616200347676 * cHG * cHG + 0.8386370679735822 * cHGtil * cHGtil
                    + 0.9774807816545604 * ctGIm * ctGIm + 0.99600530493426 * ctGRe * ctGRe
                    + 0.0055871762794901 * ctHIm * ctHIm + 0.0046644155584449 * ctHRe * ctHRe
                    + 0.1513855248714063 * ctWIm * ctWIm + 0.0223113988052159 * ctWRe * ctWRe
                    + -0.004353109781129 * cHG * cHGtil + -0.0059424753457291 * cHG * ctGIm
                    + -1.1589183159553875 * cHG * ctGRe + -0.0001031669131571 * cHG * ctHIm
                    + -0.0808290830781629 * cHG * ctHRe + 0.0023298372266031 * cHG * ctWIm
                    + 0.0007355367034189 * cHG * ctWRe + -0.8622884150085206 * cHGtil * ctGIm
                    + -0.0085612044843103 * cHGtil * ctGRe + 0.0893168984995095 * cHGtil * ctHIm
                    + -0.0007889628027357 * cHGtil * ctHRe + -0.0049693652224635 * cHGtil * ctWIm
                    + -0.003888785938789 * cHGtil * ctWRe + 0.0044373833815463 * ctGIm * ctGRe
                    + -0.0398407363498014 * ctGIm * ctHIm + 0.0001267839505427 * ctGIm * ctHRe
                    + 0.0068555328240843 * ctGIm * ctWIm + -0.0045838333288629 * ctGIm * ctWRe
                    + -0.0003431226109376 * ctGRe * ctHIm + 0.0759842998525177 * ctGRe * ctHRe
                    + 0.0412398880264256 * ctGRe * ctWIm + 0.0037903964632825 * ctGRe * ctWRe
                    + -0.0001532894777914 * ctHIm * ctHRe + -0.0002395101986699 * ctHIm * ctWIm
                    + -0.0003104171792067 * ctHIm * ctWRe + 0.0001367863317469 * ctHRe * ctWIm
                    + -0.0008153939464238 * ctHRe * ctWRe + -0.0022824930758392 * ctWIm * ctWRe
                        -0.0007691736987786996 * ctBRe + 0.002162289681735499 * ctBRe * ctBRe+ -0.00010516033086452325 * ctBIm + 0.0005771520876739905 * ctBIm * ctBIm
                        )*(SM_b4_tth_bin_m1_m0p5/b4_tth_bin_m1_m0p5_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        else{
                total =  SM_b4_tth_bin_m1_m0p5 + (0.9099876415629798*cHG + 
                        -1.1208249512813937*ctGRe + -0.1310184227297144*ctHRe
                        -0.0007691736987786996 * ctBRe -0.0013815505152491303 * ctWRe
                        )*(SM_b4_tth_bin_m1_m0p5/b4_tth_bin_m1_m0p5_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
        }
       
    } else if(b_min == -0.5 && b_max == 0.){
        
        double SM_b4_tth_bin_m0p5_0 = SM.getOptionalParameter("SM_b4_ttH_LO_bin_m0p5_0");
        double b4_tth_bin_m0p5_0_madgraph_LO = 1.932424901417113;
        double total;
        

        if(flag_Quadratic){
                total =  SM_b4_tth_bin_m0p5_0 + (2.564246291660702 * cHG + -0.005823342315884 * cHGtil
                    + -0.0016743198961839 * ctGIm + -2.778424493033977 * ctGRe
                    + 0.0007195133202335 * ctHIm + -0.2344245527833774 * ctHRe
                    + -0.0064220762451969 * ctWIm + -0.001215693713534 * ctWRe
                    + 2.034099761021467 * cHG * cHG + 1.9612673310063 * cHGtil * cHGtil
                    + 2.5341046159273257 * ctGIm * ctGIm + 2.544939547558986 * ctGRe * ctGRe
                    + 0.0039923939634707 * ctHIm * ctHIm + 0.0077910494317363 * ctHRe * ctHRe
                    + 0.4588069842581513 * ctWIm * ctWIm + 0.065854083346678 * ctWRe * ctWRe
                    + -0.0071569239960523 * cHG * cHGtil + 0.0145336259449214 * cHG * ctGIm
                    + -2.277554920224764 * cHG * ctGRe + -4.978015231518844e-05 * cHG * ctHIm
                    + -0.2068188648732378 * cHG * ctHRe + -0.0040071191552686 * cHG * ctWIm
                    + 0.0002867576099352 * cHG * ctWRe + -1.4549114499408065 * cHGtil * ctGIm
                    + -0.026659343007864 * cHGtil * ctGRe + 0.1064186672197977 * cHGtil * ctHIm
                    + -0.002336704675229 * cHGtil * ctHRe + 0.0024569029783631 * cHGtil * ctWIm
                    + -0.0026429889231636 * cHGtil * ctWRe + -0.0026918504123321 * ctGIm * ctGRe
                    + 0.0026133954963231 * ctGIm * ctHIm + -0.0008204230232737 * ctGIm * ctHRe
                    + -0.0313376168904284 * ctGIm * ctWIm + 0.0058758543148654 * ctGIm * ctWRe
                    + 4.366654352272459e-06 * ctGRe * ctHIm + 0.1771923049354946 * ctGRe * ctHRe
                    + 0.0063163081295789 * ctGRe * ctWIm + -0.0038756154471681 * ctGRe * ctWRe
                    + -1.0308033893416904e-05 * ctHIm * ctHRe + -0.0012053193562682 * ctHIm * ctWIm
                    + -8.001925082197133e-05 * ctHIm * ctWRe + 0.0002104946911641 * ctHRe * ctWIm
                    + -0.0013179836923686 * ctHRe * ctWRe + 0.0093561728027508 * ctWIm * ctWRe
                        -0.00331238196421145 * ctBRe + 0.008086112326838674 * ctBRe * ctBRe+ 0.00015893289734032756 * ctBIm + 0.003681537592577791 * ctBIm * ctBIm
                        )*(SM_b4_tth_bin_m0p5_0/b4_tth_bin_m0p5_0_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        else{
                total =  SM_b4_tth_bin_m0p5_0 + (2.564414022112631*cHG + 
                        -2.7816281236042366*ctGRe -0.2344736662557363*ctHRe
                        -0.00331238196421145 * ctBRe -0.0012156937135340828 * ctWRe
                        )*(SM_b4_tth_bin_m0p5_0/b4_tth_bin_m0p5_0_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
        }

    } else if(b_min == 0. && b_max == 0.5){
        
        double SM_b4_tth_bin_0_0p5 = SM.getOptionalParameter("SM_b4_ttH_LO_bin_0_0p5");
        double b4_tth_bin_0_0p5_madgraph_LO = 2.48186478266304;
        double total;
        

            if(flag_Quadratic){
                total = SM_b4_tth_bin_0_0p5 + (3.2093136673424736 * cHG + 0.004029159850137 * cHGtil
                    + -0.0213738731216522 * ctGIm + -3.504682622127001 * ctGRe
                    + -0.0002875376687773 * ctHIm + -0.3032103265837639 * ctHRe
                    + 0.0095376518003493 * ctWIm + -0.0058800129619509 * ctWRe
                    + 2.508732816346594 * cHG * cHG + 2.437345459874909 * cHGtil * cHGtil
                    + 3.1881280812208166 * ctGIm * ctGIm + 3.165825489942664 * ctGRe * ctGRe
                    + 0.0036148729182692 * ctHIm * ctHIm + 0.0099998321780285 * ctHRe * ctHRe
                    + 0.458074038263705 * ctWIm * ctWIm + 0.0758786285553478 * ctWRe * ctWRe
                    + 0.0053996492403103 * cHG * cHGtil + -0.0236162521426582 * cHG * ctGIm
                    + -2.7111214421094987 * cHG * ctGRe + -0.001177880066321 * cHG * ctHIm
                    + -0.2508174145511634 * cHG * ctHRe + 0.0027312222369129 * cHG * ctWIm
                    + -0.0005506402375662 * cHG * ctWRe + -1.712920295286303 * cHGtil * ctGIm
                    + 0.0148464169032713 * cHGtil * ctGRe + 0.1255660004129038 * cHGtil * ctHIm
                    + 0.0014699068562402 * cHGtil * ctHRe + 0.0055682588248907 * cHGtil * ctWIm
                    + -0.0006447342027894 * cHGtil * ctWRe + -0.0113505890042336 * ctGIm * ctGRe
                    + 0.0159020808765088 * ctGIm * ctHIm + -0.0004755818091291 * ctGIm * ctHRe
                    + 0.025203510325142 * ctGIm * ctWIm + -0.0009750172709239 * ctGIm * ctWRe
                    + -0.000673228149598 * ctGRe * ctHIm + 0.2243402305478286 * ctGRe * ctHRe
                    + -0.0176639654745808 * ctGRe * ctWIm + -0.0002779594299604 * ctGRe * ctWRe
                    + -3.642230847214099e-05 * ctHIm * ctHRe + -0.0017585157882703 * ctHIm * ctWIm
                    + 0.0004703647702423 * ctHIm * ctWRe + -0.000285286634482 * ctHRe * ctWIm
                    + -0.0026938458799784 * ctHRe * ctWRe + -0.0017642667261674 * ctWIm * ctWRe
                        -0.005253102367525533 * ctBRe + 0.012034842573137972 * ctBRe * ctBRe+ -0.00011904416973186616 * ctBIm + 0.005659950685029516 * ctBIm * ctBIm
                        )*(SM_b4_tth_bin_0_0p5/b4_tth_bin_0_0p5_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

            }
            else{
                total = SM_b4_tth_bin_0_0p5 + (3.2093661829551725*cHG +
                        -3.501689241803734*ctGRe + -0.3031442522037286*ctHRe
                        -0.005253102367525533 * ctBRe - 0.005880012961950926 * ctWRe
                        )*(SM_b4_tth_bin_0_0p5/b4_tth_bin_0_0p5_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
            }

    } else if(b_min == 0.5 && b_max == 1.){
        
        double SM_b4_tth_bin_0p5_1 = SM.getOptionalParameter("SM_b4_ttH_LO_bin_0p5_1");
        double b4_tth_bin_0p5_1_madgraph_LO = 4.587536187686661;//pb
        double total;
        

            if(flag_Quadratic){
                total = SM_b4_tth_bin_0p5_1 + (3.8624206533859455 * cHG + 0.0025652805438462 * cHGtil
                    + 0.0380812173024499 * ctGIm + -4.967609706819283 * ctGRe
                    + -0.0006965496592073 * ctHIm + -0.5544826010765769 * ctHRe
                    + 0.0006073357310019 * ctWIm + -0.0256463227849223 * ctWRe
                    + 2.347939739733353 * cHG * cHG + 2.3092071630867377 * cHGtil * cHGtil
                    + 3.199130314657721 * ctGIm * ctGIm + 3.427832512642536 * ctGRe * ctGRe
                    + 0.004453926977181 * ctHIm * ctHIm + 0.0167718441092212 * ctHRe * ctHRe
                    + 0.159767945592875 * ctWIm * ctWIm + 0.059120496919097 * ctWRe * ctWRe
                    + 0.0072252855175872 * cHG * cHGtil + 0.03975689764912 * cHG * ctGIm
                    + -3.1711629756704776 * cHG * ctGRe + 0.0014716333628268 * cHG * ctHIm
                    + -0.2420802478163346 * cHG * ctHRe + 0.0018271931377158 * cHG * ctWIm
                    + -0.0002326603969194 * cHG * ctWRe + -1.898959594407268 * cHGtil * ctGIm
                    + 0.014375573266062 * cHGtil * ctGRe + 0.188418633982538 * cHGtil * ctHIm
                    + 0.0007511425222509 * cHGtil * ctHRe + -0.0014741276289671 * cHGtil * ctWIm
                    + 0.007150320736038 * cHGtil * ctWRe + 0.0037375122760155 * ctGIm * ctGRe
                    + 0.0081175256170974 * ctGIm * ctHIm + 0.0007985463884474 * ctGIm * ctHRe
                    + -0.0001131618468446 * ctGIm * ctWIm + 0.0008483977113818 * ctGIm * ctWRe
                    + 9.220570277163054e-05 * ctGRe * ctHIm + 0.3707695293252685 * ctGRe * ctHRe
                    + -0.0295125287850878 * ctGRe * ctWIm + 0.0053821483274255 * ctGRe * ctWRe
                    + 0.0001626805165413 * ctHIm * ctHRe + -0.0008105961019141 * ctHIm * ctWIm
                    + 0.0001035970931805 * ctHIm * ctWRe + 6.611959525185118e-06 * ctHRe * ctWIm
                    + 0.0003098374575341 * ctHRe * ctWRe + -0.0025092774493894 * ctWIm * ctWRe
                        -0.009187820270389016 * ctBRe + 0.016843223885503147 * ctBRe * ctBRe+ 5.500801554077228e-05 * ctBIm + 0.005015742054571379 * ctBIm * ctBIm
                        )*(SM_b4_tth_bin_0p5_1/b4_tth_bin_0p5_1_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

            }
            else{
                total = SM_b4_tth_bin_0p5_1 + (3.862512212641887*cHG +
                        -4.967994348065772*ctGRe + -0.554551144413014*ctHRe
                        -0.009187820270389016 * ctBRe -0.02564632278492235 * ctWRe
                        )*(SM_b4_tth_bin_0p5_1/b4_tth_bin_0p5_1_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
            }

    } else{
        throw std::runtime_error("\nERROR: Please specify a correct binning for b4_tth_LO. \n");
    }    
    
}





Asymmetry_Dazi_ord_ttH::Asymmetry_Dazi_ord_ttH(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

//    setParametersForObservable(make_vector<std::string>() << "SM_Asymmetry_Dazi_ord_ttH");
    
}

double Asymmetry_Dazi_ord_ttH::computeThValue()
{
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiG();
    double cHGtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiGtil();
    //double cHW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiW();
    //double cHWtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiWtil();
    double ctHRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double ctHIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphiIm();
    double ctGRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double ctGIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tGIm();
    double ctWRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double ctWIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tWIm();
    double ctBRe = 0.;
    double ctBIm = 0.;
    
    //Expressions are written in the basis of SMEFTsim
    //let's change to the basis of dim6top (same as
    //SMEFT@NLO except for the gs on CtG). These
    //models define the covariance derivative 
    //with a different convention
    
    if(flag_LHC_WG_Basis){
        double ctZRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
        double ctZIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZIm();
        double sw2 = 0.22305;
        double cw2 = 1 - sw2;
        double tw2 = sw2/cw2;
        double sw  = sqrt(sw2);
        double tw  = sqrt(tw2);
        ctBRe = ctZRe/sw - ctWRe/tw;
        ctBIm = ctZIm/sw - ctWIm/tw;
        ctWRe = -ctWRe;
        ctWIm = -ctWIm;
        ctGRe = -ctGRe;
        ctGIm = -ctGIm;
    }
    else{
        ctBRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        ctBIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tBIm();
    }
    
    
    
    
    
    //double SM_Asymmetry_Dazi_ord = SM.getOptionalParameter("SM_Asymmetry_Dazi_ord");
    
    double Dazi_ord_neg_Mad = 0.801238412607049;
    double Dazi_ord_pos_Mad = 0.8006171719284476;
    


    if(flag_Quadratic){

            double Dazi_ord_neg_NP = 0.8409726917581889 * cHG + 0.0239521993589093 * cHGtil
                + -0.0727650754320362 * ctGIm + -0.9799091583330292 * ctGRe
                + -0.0001997761653249 * ctHIm + -0.0974439820691808 * ctHRe
                + 0.0040841173511341 * ctWIm + -0.0026346877927101 * ctWRe
                + 0.6180724071879639 * cHG * cHG + 0.6006261707767361 * cHGtil * cHGtil
                + 0.7856366780703115 * ctGIm * ctGIm + 0.8059251252505129 * ctGRe * ctGRe
                + 0.0013971371538088 * ctHIm * ctHIm + 0.0031357537523052 * ctHRe * ctHRe
                + 0.0977883829166236 * ctWIm * ctWIm + 0.0177777400561392 * ctWRe * ctWRe
                + 0.0157290488652358 * cHG * cHGtil + -0.0033975293377011 * cHG * ctGIm
                + -0.7471772164696022 * cHG * ctGRe + 0.0006005292045924 * cHG * ctHIm
                + -0.0620994359176081 * cHG * ctHRe + 0.0014581468930783 * cHG * ctWIm
                + -5.417222834249148e-05 * cHG * ctWRe + -0.4707346806682253 * cHGtil * ctGIm
                + -0.0258435683811506 * cHGtil * ctGRe + 0.0402500484934529 * cHGtil * ctHIm
                + -0.0011440632677678 * cHGtil * ctHRe + 0.0003948705556123 * cHGtil * ctWIm
                + -3.223394386769663e-05 * cHGtil * ctWRe + 0.019267926306919 * ctGIm * ctGRe
                + -0.0011549478510061 * ctGIm * ctHIm + 0.0043931390740712 * ctGIm * ctHRe
                + -0.0026955097559458 * ctGIm * ctWIm + -2.6134787065639213e-05 * ctGIm * ctWRe
                + -0.0039137945238467 * ctGRe * ctHIm + 0.0672010069425599 * ctGRe * ctHRe
                + -0.0042875376558723 * ctGRe * ctWIm + 0.0009131837278532 * ctGRe * ctWRe
                + -1.0463403882919558e-06 * ctHIm * ctHRe + -0.0003603290153253 * ctHIm * ctWIm
                + 0.0002103912520551 * ctHIm * ctWRe + -3.3601141254107565e-05 * ctHRe * ctWIm
                + -0.0003184009878738 * ctHRe * ctWRe + 9.336705314867688e-05 * ctWIm * ctWRe
            -0.0013981083710685738 * ctBRe + 0.0031128276299091096 * ctBRe * ctBRe+ -6.751260631875089e-05 * ctBIm + 0.0011880155398607123 * ctBIm * ctBIm
            ;
            double Dazi_ord_pos_NP = 0.8374737681273691 * cHG + -0.0239936915366412 * cHGtil
                + 0.0740756917886777 * ctGIm + -0.9887205393601546 * ctGRe
                + 0.0001758920216041 * ctHIm + -0.0972082903115259 * ctHRe
                + -0.0038178577929426 * ctWIm + -0.0027960955347098 * ctWRe
                + 0.6174050704291627 * cHG * cHG + 0.6004297673773595 * cHGtil * cHGtil
                + 0.7898132427425815 * ctGIm * ctGIm + 0.8070470149243939 * ctGRe * ctGRe
                + 0.0014116881912348 * ctHIm * ctHIm + 0.0031074396853647 * ctHRe * ctHRe
                + 0.0976538671279632 * ctWIm * ctWIm + 0.0177390091018495 * ctWRe * ctWRe
                + -0.0155516480557833 * cHG * cHGtil + 0.00733235897668 * cHG * ctGIm
                + -0.7355705059116835 * cHG * ctGRe + -0.000578124797484 * cHG * ctHIm
                + -0.0620909954882131 * cHG * ctHRe + -0.0009997175215992 * cHG * ctWIm
                + 9.219615703674744e-05 * cHG * ctWRe + -0.4727237329541328 * cHGtil * ctGIm
                + 0.0248890286688689 * cHGtil * ctGRe + 0.0408622177561639 * cHGtil * ctHIm
                + 0.0010001226021191 * cHGtil * ctHRe + -0.0001432347113879 * cHGtil * ctWIm
                + 2.80668204017398e-05 * cHGtil * ctWRe + -0.0202015021970715 * ctGIm * ctGRe
                + -0.0009467910567448 * ctGIm * ctHIm + -0.0044521278787526 * ctGIm * ctHRe
                + 0.0027922955345794 * ctGIm * ctWIm + 0.0002115650012277 * ctGIm * ctWRe
                + 0.0037674177682461 * ctGRe * ctHIm + 0.0677935966148256 * ctGRe * ctHRe
                + 0.0043479597737826 * ctGRe * ctWIm + -0.0001145085088655 * ctGRe * ctWRe
                + -4.895422227624291e-06 * ctHIm * ctHRe + -0.0002782677303533 * ctHIm * ctWIm
                + -0.0001811936976982 * ctHIm * ctWRe + 4.451757689860614e-05 * ctHRe * ctWIm
                + -0.0004004144226487 * ctHRe * ctWRe + 0.0003520788300325 * ctWIm * ctWRe
            -0.001549419314305453 * ctBRe + 0.0031134638609745385 * ctBRe * ctBRe+ 6.587933630411436e-05 * ctBIm + 0.0011885295925029155 * ctBIm * ctBIm
            ;
    
    
            //double num = (Dazi_ord_pos_Mad + Dazi_ord_pos_NP) - (Dazi_ord_neg_Mad + Dazi_ord_neg_NP);
            //double deno = (Dazi_ord_pos_Mad + Dazi_ord_pos_NP) + (Dazi_ord_neg_Mad + Dazi_ord_neg_NP);
    
            
            //We set the SM difference to zero
            double num = (Dazi_ord_pos_NP) - (Dazi_ord_neg_NP);
            double deno = (Dazi_ord_pos_Mad + Dazi_ord_pos_NP) + (Dazi_ord_neg_Mad + Dazi_ord_neg_NP);
    
            double total = num/deno;
            return total;
            
    }
    else{

            
            
            double Dazi_ord_neg_NP = 0.8409727694620143*cHG + 0.0239522015720341*cHGtil + -0.0727650821867264*ctGIm +
            -0.9799092492967638*ctGRe + -0.0001997762114054*ctHIm + -0.0974440045456804*ctHRe
            -0.0013981083710685738 * ctBRe+ -6.751260631875089e-05 * ctBIm
            -0.0026346877927101164 * ctWRe+ 0.004084117351134148 * ctWIm
            ;
            double Dazi_ord_pos_NP = 0.8374738455079025*cHG + -0.0239936937535997*cHGtil + 0.0740756986650311*ctGIm +
            -0.9887206311418384*ctGRe + 0.0001758920621755*ctHIm + -0.0972083127336606*ctHRe
            -0.001549419314305453 * ctBRe+ 6.587933630411436e-05 * ctBIm
            -0.0027960955347098027 * ctWRe+ -0.003817857792942636 * ctWIm
            ;
    
    
            //double num = (Dazi_ord_pos_Mad + Dazi_ord_pos_NP) - (Dazi_ord_neg_Mad + Dazi_ord_neg_NP);
            //double deno = (Dazi_ord_pos_Mad + Dazi_ord_pos_NP) + (Dazi_ord_neg_Mad + Dazi_ord_neg_NP);
    
            
            //We set the SM difference to zero
            double num = (Dazi_ord_pos_NP) - (Dazi_ord_neg_NP);
            double deno = (Dazi_ord_pos_Mad + Dazi_ord_pos_NP) + (Dazi_ord_neg_Mad + Dazi_ord_neg_NP);
    
            
            
            double total = num/deno;
            return total;
            
    }

}


Asymmetry_Dazi_ord_ttH_ee::Asymmetry_Dazi_ord_ttH_ee(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

//    setParametersForObservable(make_vector<std::string>() << "SM_Asymmetry_Dazi_ord_ttH_ee");
    
}

double Asymmetry_Dazi_ord_ttH_ee::computeThValue()
{
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiG();
    double cHGtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiGtil();
    //double cHW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiW();
    //double cHWtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiWtil();
    double ctHRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double ctHIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphiIm();
    double ctGRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double ctGIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tGIm();
    double ctWRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double ctWIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tWIm();
    double ctBRe = 0.;
    double ctBIm = 0.;
    
    
    //Expressions are written in the basis of SMEFTsim
    //let's change to the basis of dim6top (same as
    //SMEFT@NLO except for the gs on CtG). These
    //models define the covariance derivative 
    //with a different convention
    
    if(flag_LHC_WG_Basis){
        double ctZRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
        double ctZIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZIm();
        double sw2 = 0.22305;
        double cw2 = 1 - sw2;
        double tw2 = sw2/cw2;
        double sw  = sqrt(sw2);
        double tw  = sqrt(tw2);
        ctBRe = ctZRe/sw - ctWRe/tw;
        ctBIm = ctZIm/sw - ctWIm/tw;
        ctWRe = -ctWRe;
        ctWIm = -ctWIm;
        ctGRe = -ctGRe;
        ctGIm = -ctGIm;
    }
    else{
        ctBRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        ctBIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tBIm();
    }
    
    
    
    
    
    
    //double SM_Asymmetry_Dazi_ord = SM.getOptionalParameter("SM_Asymmetry_Dazi_ord");
    
    double Dazi_ord_ee_neg_Mad = 0.8006766500797665;
    double Dazi_ord_ee_pos_Mad = 0.8011790447500986;
    

    if(flag_Quadratic){

        
            double Dazi_ord_ee_neg_NP = 0.8404287198682079 * cHG + 0.0208422614129312 * cHGtil
                + -0.0550089310565959 * ctGIm + -0.983778963439606 * ctGRe
                + 0.0016550942929671 * ctHIm + -0.0975507985435286 * ctHRe
                + 0.0239873280104936 * ctWIm + -0.0026831297740093 * ctWRe
                + 0.6176130449506895 * cHG * cHG + 0.6001347124313424 * cHGtil * cHGtil
                + 0.7865719489740035 * ctGIm * ctGIm + 0.80611827533555 * ctGRe * ctGRe
                + 0.0014229577967017 * ctHIm * ctHIm + 0.0031487266504941 * ctHRe * ctHRe
                + 0.0977314966404424 * ctWIm * ctWIm + 0.0177935600401541 * ctWRe * ctWRe
                + 0.0040914045299683 * cHG * cHGtil + -0.0219899064314808 * cHG * ctGIm
                + -0.7388365912386904 * cHG * ctGRe + -0.0003736312654786 * cHG * ctHIm
                + -0.061822929264468 * cHG * ctHRe + 4.755266873937569e-05 * cHG * ctWIm
                + 0.0004908479318843 * cHG * ctWRe + -0.4730595292823339 * cHGtil * ctGIm
                + 0.0067909859506482 * cHGtil * ctGRe + 0.0404155989650598 * cHGtil * ctHIm
                + -0.0006071631925375 * cHGtil * ctHRe + 0.000410970632705 * cHGtil * ctWIm
                + -0.0010873287770317 * cHGtil * ctWRe + 0.0067907908337478 * ctGIm * ctGRe
                + -0.0012276151256363 * ctGIm * ctHIm + 0.0042570184763379 * ctGIm * ctHRe
                + -0.0077141578106353 * ctGIm * ctWIm + -0.0004329878738224 * ctGIm * ctWRe
                + -0.0026512985824996 * ctGRe * ctHIm + 0.0673226914661984 * ctGRe * ctHRe
                + 0.0058230088343248 * ctGRe * ctWIm + -5.741735732788912e-05 * ctGRe * ctWRe
                + -0.0001332105072905 * ctHIm * ctHRe + -0.0003511364464264 * ctHIm * ctWIm
                + 0.000125880399987 * ctHIm * ctWRe + -0.0004405570611716 * ctHRe * ctWIm
                + -0.0004815913190063 * ctHRe * ctWRe + -0.000362755095259 * ctWIm * ctWRe
            -0.0016684521282414322 * ctBRe + 0.0031109123237574587 * ctBRe * ctBRe+ -4.204517990341605e-05 * ctBIm + 0.0011882036111250938 * ctBIm * ctBIm
            ;
            double Dazi_ord_ee_pos_NP = 0.8380177400173502 * cHG + -0.0208837535906631 * cHGtil
                + 0.0563195474132376 * ctGIm + -0.9848507342535776 * ctGRe
                + -0.0016789784366879 * ctHIm + -0.0971014738371782 * ctHRe
                + -0.0237210684368537 * ctWIm + -0.0027476538685015 * ctWRe
                + 0.6178644326664374 * cHG * cHG + 0.6009212257227535 * cHGtil * cHGtil
                + 0.7888779718388894 * ctGIm * ctGIm + 0.8068538648393566 * ctGRe * ctGRe
                + 0.0013858675483419 * ctHIm * ctHIm + 0.0030944667871758 * ctHRe * ctHRe
                + 0.097710764743593 * ctWIm * ctWIm + 0.0177231911784963 * ctWRe * ctWRe
                + -0.0039140034663799 * cHG * cHGtil + 0.0259247403583757 * cHG * ctGIm
                + -0.7439138030158363 * cHG * ctGRe + 0.0003960357139343 * cHG * ctHIm
                + -0.0623674908935729 * cHG * ctHRe + 0.0004108770005809 * cHG * ctWIm
                + -0.0004528240008638 * cHG * ctWRe + -0.4703992585666175 * cHGtil * ctGIm
                + -0.0077455265899209 * cHGtil * ctGRe + 0.0406966990783182 * cHGtil * ctHIm
                + 0.0004632224275368 * cHGtil * ctHRe + -0.0001593347843943 * cHGtil * ctWIm
                + 0.0010831616538765 * cHGtil * ctWRe + -0.0077243666834214 * ctGIm * ctGRe
                + -0.0008741250260713 * ctGIm * ctHIm + -0.0043160072941473 * ctGIm * ctHRe
                + 0.0078109435808307 * ctGIm * ctWIm + 0.0006184181227847 * ctGIm * ctWRe
                + 0.0025049217993331 * ctGRe * ctHIm + 0.0676719742864739 * ctGRe * ctHRe
                + -0.0057625867117164 * ctGRe * ctWIm + 0.0008560926014874 * ctGRe * ctWRe
                + 0.0001272687320934 * ctHIm * ctHRe + -0.0002874604068827 * ctHIm * ctWIm
                + -9.668285932849868e-05 * ctHIm * ctWRe + 0.0004514735017376 * ctHRe * ctWIm
                + -0.0002372243269866 * ctHRe * ctWRe + 0.0008082013407264 * ctWIm * ctWRe
            -0.0012790754463039983 * ctBRe + 0.003115378933016111 * ctBRe * ctBRe+ 4.0411909950716085e-05 * ctBIm + 0.0011883414318798607 * ctBIm * ctBIm
            ;
    
    
            //double num = (Dazi_ord_ee_pos_Mad + Dazi_ord_ee_pos_NP) - (Dazi_ord_ee_neg_Mad + Dazi_ord_ee_neg_NP);
            //double deno = (Dazi_ord_ee_pos_Mad + Dazi_ord_ee_pos_NP) + (Dazi_ord_ee_neg_Mad + Dazi_ord_ee_neg_NP);
            
            
            //We set the SM difference to zero
            double num = (Dazi_ord_ee_pos_NP) - (Dazi_ord_ee_neg_NP);
            double deno = (Dazi_ord_ee_pos_Mad + Dazi_ord_ee_pos_NP) + (Dazi_ord_ee_neg_Mad + Dazi_ord_ee_neg_NP);
    
            double total = num/deno;
            return total;
            
        }
        else{

            
    
            double Dazi_ord_ee_neg_NP = 0.8404288553887629*cHG + 0.0208422647737813*cHGtil + -0.0550100485100096*ctGIm +
            -0.9837800350478916*ctGRe + 0.0016550964830614*ctHIm + -0.0975509276270819*ctHRe
            -0.0016684521282414322 * ctBRe+ -4.204517990341605e-05 * ctBIm
            -0.002683129774009377 * ctWRe+ 0.02398732801049369 * ctWIm
            ;
            double Dazi_ord_ee_pos_NP = 0.8380178751491306*cHG + -0.0208837569582036*cHGtil + 0.0563206662942763*ctGIm +
            -0.9848518070293188*ctGRe + -0.0016789806583868*ctHIm + -0.0971016023261651*ctHRe
            -0.0012790754463039983 * ctBRe+ 4.0411909950716085e-05 * ctBIm
            -0.002747653868501579 * ctWRe+ -0.023721068436853746 * ctWIm
            ;
    
    
            //double num = (Dazi_ord_ee_pos_Mad + Dazi_ord_ee_pos_NP) - (Dazi_ord_ee_neg_Mad + Dazi_ord_ee_neg_NP);
            //double deno = (Dazi_ord_ee_pos_Mad + Dazi_ord_ee_pos_NP) + (Dazi_ord_ee_neg_Mad + Dazi_ord_ee_neg_NP);
            
            
            //We set the SM difference to zero
            double num = (Dazi_ord_ee_pos_NP) - (Dazi_ord_ee_neg_NP);
            double deno = (Dazi_ord_ee_pos_Mad + Dazi_ord_ee_pos_NP) + (Dazi_ord_ee_neg_Mad + Dazi_ord_ee_neg_NP);
            
            
            
            double total = num/deno;
            return total;

    }
}


sigma_ttH_diff_LO_mtth::sigma_ttH_diff_LO_mtth(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_ttH_diff_LO_mtth_bin_450_655" << "SM_sigma_ttH_diff_LO_mtth_bin_655_860"
            << "SM_sigma_ttH_diff_LO_mtth_bin_860_1270" << "SM_sigma_ttH_diff_LO_mtth_bin_1270_2500");

    
}

double sigma_ttH_diff_LO_mtth::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiG();
    double cHGtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiGtil();
    //double cHW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiW();
    //double cHWtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiWtil();
    double ctHRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double ctHIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphiIm();
    double ctGRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double ctGIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tGIm();
    double ctWRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double ctWIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tWIm();
    double ctBRe = 0.;
    double ctBIm = 0.;
    
    
    
    //Expressions are written in the basis of SMEFTsim
    //let's change to the basis of dim6top (same as
    //SMEFT@NLO except for the gs on CtG). These
    //models define the covariance derivative 
    //with a different convention
    
    if(flag_LHC_WG_Basis){
        double ctZRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
        double ctZIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZIm();
        double sw2 = 0.22305;
        double cw2 = 1 - sw2;
        double tw2 = sw2/cw2;
        double sw  = sqrt(sw2);
        double tw  = sqrt(tw2);
        ctBRe = ctZRe/sw - ctWRe/tw;
        ctBIm = ctZIm/sw - ctWIm/tw;
        ctWRe = -ctWRe;
        ctWIm = -ctWIm;
        ctGRe = -ctGRe;
        ctGIm = -ctGIm;
    }
    else{
        ctBRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        ctBIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tBIm();
    }
    
    
   
    if(b_min == 450 && b_max == 655){
        
        double SM_ttH_bin_450_655 = SM.getOptionalParameter("SM_sigma_ttH_diff_LO_mtth_bin_450_655");
        double ttH_bin_450_655_madgraph_LO = 0.0085968515036142;//pb/GeV 
        double total;
        
        
        if(flag_Quadratic){
                total =  SM_ttH_bin_450_655 + (0.0045181955745771 * cHG + 2.092434682771094e-06 * cHGtil
                    + 3.35721352610896e-06 * ctGIm + -0.0073319817382764 * ctGRe
                    + 2.178950715664429e-06 * ctHIm + -0.0010488957886131 * ctHRe
                    + -6.397754803461053e-06 * ctWIm + -3.0007469619147686e-05 * ctWRe
                    + 0.001129479117399 * cHG * cHG + 0.0012542998928326 * cHGtil * cHGtil
                    + 0.0027597465911086 * ctGIm * ctGIm + 0.0031944543434764 * ctGRe * ctGRe
                    + 5.658099728467203e-06 * ctHIm * ctHIm + 3.208560834544569e-05 * ctHRe * ctHRe
                    + 3.1330082486985544e-05 * ctWIm * ctWIm + 3.6107605324769976e-05 * ctWRe * ctWRe
                    + 5.753962177779126e-05 * cHG * cHGtil + 2.8298667156714807e-05 * cHG * ctGIm
                    + -0.0029193903114386 * cHG * ctGRe + 3.5475454384648053e-06 * cHG * ctHIm
                    + -0.0002322100406697 * cHG * ctHRe + -4.0326775995974525e-08 * cHG * ctWIm
                    + 9.549539522979432e-06 * cHG * ctWRe + -0.0013688338065862 * cHGtil * ctGIm
                    + -7.438424425927569e-05 * cHGtil * ctGRe + 0.0003213688692185 * cHGtil * ctHIm
                    + -4.67918246424609e-07 * cHGtil * ctHRe + -3.209779568677406e-06 * cHGtil * ctWIm
                    + -1.756648507190347e-05 * cHGtil * ctWRe + 9.308637288288346e-05 * ctGIm * ctGRe
                    + -1.0905995485980285e-05 * ctGIm * ctHIm + 9.160715861988856e-07 * ctGIm * ctHRe
                    + 4.055744317929544e-06 * ctGIm * ctWIm + 4.462759282716777e-06 * ctGIm * ctWRe
                    + 4.508771085515173e-07 * ctGRe * ctHIm + 0.0006295005037564 * ctGRe * ctHRe
                    + 1.8551257867788214e-05 * ctGRe * ctWIm + 9.718015508280995e-07 * ctGRe * ctWRe
                    + -2.1921421677216978e-07 * ctHIm * ctHRe + 4.255271382203052e-07 * ctHIm * ctWIm
                    + -5.998001470667131e-07 * ctHIm * ctWRe + 1.6381698403344913e-07 * ctHRe * ctWIm
                    + 2.2094272996431297e-06 * ctHRe * ctWRe + -4.901345972218524e-10 * ctWIm * ctWRe
                    -8.307825152903403e-06 * ctBRe + 1.1261681936616807e-05 * ctBRe * ctBRe + 9.963448805855313e-07 * ctBIm * ctBIm
                    )*(SM_ttH_bin_450_655/ttH_bin_450_655_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        else{
                total =  SM_ttH_bin_450_655 + (0.0045181955745771*cHG +
                        -0.0073319817382764*ctGRe + -0.0010488957886131*ctHRe
                        -8.307825152903403e-06 * ctBRe-3.0007469619147686e-05 * ctWRe
                        )*(SM_ttH_bin_450_655/ttH_bin_450_655_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
        
        }
        
    } else if(b_min == 655 && b_max == 860){
        
        double SM_ttH_bin_655_860 = SM.getOptionalParameter("SM_sigma_ttH_diff_LO_mtth_bin_655_860");
        double ttH_bin_655_860_madgraph_LO = 0.0089470371846531;//pb/GeV
        double total;
 

        if(flag_Quadratic){
                total =  SM_ttH_bin_655_860 + (0.0068342312470782 * cHG + 6.516776100529143e-06 * cHGtil
                    + 1.7655998727092448e-05 * ctGIm + -0.0096851624939012 * ctGRe
                    + -8.016268841642971e-07 * ctHIm + -0.0010829961514144 * ctHRe
                    + 4.5660108292516766e-06 * ctWIm + -3.266069476552669e-05 * ctWRe
                    + 0.0032625736708071 * cHG * cHG + 0.0034043208850053 * cHGtil * cHGtil
                    + 0.0056434013751323 * ctGIm * ctGIm + 0.0061425691862662 * ctGRe * ctGRe
                    + 1.34073766868692e-05 * ctHIm * ctHIm + 3.306703774944704e-05 * ctHRe * ctHRe
                    + 0.0001699478599239 * ctWIm * ctWIm + 7.677430471397795e-05 * ctWRe * ctWRe
                    + -2.364760048584498e-05 * cHG * cHGtil + -1.4463105183893285e-05 * cHG * ctGIm
                    + -0.005825556907969 * cHG * ctGRe + -1.0233814119407383e-05 * cHG * ctHIm
                    + -0.0003562114438536 * cHG * ctHRe + 4.041138353352847e-06 * cHG * ctWIm
                    + 2.0052515994350983e-06 * cHG * ctWRe + -0.003378026398263 * cHGtil * ctGIm
                    + 8.076642084585417e-05 * cHGtil * ctGRe + 0.0004541628998359 * cHGtil * ctHIm
                    + 5.62908023971076e-07 * cHGtil * ctHRe + 6.153746747686317e-06 * cHGtil * ctWIm
                    + 4.800465880375193e-06 * cHGtil * ctWRe + -7.617018209962787e-05 * ctGIm * ctGRe
                    + -3.316085847073677e-05 * ctGIm * ctHIm + -2.770873236021253e-06 * ctGIm * ctHRe
                    + -7.011545265733163e-06 * ctGIm * ctWIm + -1.1505958574542414e-05 * ctGIm * ctWRe
                    + 1.7917779605858275e-06 * ctGRe * ctHIm + 0.000739929641391 * ctGRe * ctHRe
                    + -3.541643525756911e-05 * ctGRe * ctWIm + -1.0605493261818856e-05 * ctGRe * ctWRe
                    + -1.2537585183551948e-07 * ctHIm * ctHRe + -1.5696577463876328e-06 * ctHIm * ctWIm
                    + 3.765760764606885e-07 * ctHIm * ctWRe + 3.604191362999165e-07 * ctHRe * ctWIm
                    + -2.1530467315235542e-06 * ctHRe * ctWRe + 4.389160141851645e-06 * ctWIm * ctWRe
                    -1.5467073895367167e-05 * ctBRe + 2.201484616365379e-05 * ctBRe * ctBRe+ 4.364785857635667e-06 * ctBIm * ctBIm
                    )*(SM_ttH_bin_655_860/ttH_bin_655_860_madgraph_LO);
                        //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        else{
                total =  SM_ttH_bin_655_860 + (0.0068342312470782*cHG +
                        -0.0096851624939012*ctGRe + -0.0010829961514144*ctHRe
                        -1.5467073895367167e-05 * ctBRe -3.266069476552669e-05 * ctWRe
                        )*(SM_ttH_bin_655_860/ttH_bin_655_860_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
        
        }
    
    } else if(b_min == 860 && b_max == 1270){
        
        double SM_ttH_bin_860_1270 = SM.getOptionalParameter("SM_sigma_ttH_diff_LO_mtth_bin_860_1270");
        double ttH_bin_860_1270_madgraph_LO = 0.0026949395920587;//pb/GeV
        double total;
 

        if(flag_Quadratic){
                total =  SM_ttH_bin_860_1270 + (0.0035482569627141 * cHG + -8.220476758191708e-07 * cHGtil
                    + 1.2052119047978316e-05 * ctGIm + -0.0041580771412648 * ctGRe
                    + -1.45859601326237e-06 * ctHIm + -0.0003266623947329 * ctHRe
                    + 7.515840323157885e-07 * ctWIm + -1.1085493448650952e-05 * ctWRe
                    + 0.0031788183581889 * cHG * cHG + 0.0030970512164488 * cHGtil * cHGtil
                    + 0.0038709717867162 * ctGIm * ctGIm + 0.0040041851524892 * ctGRe * ctGRe
                    + 7.55327927619287e-06 * ctHIm * ctHIm + 1.044536285490516e-05 * ctHRe * ctHRe
                    + 0.0003022025414976 * ctWIm * ctWIm + 7.053103081616774e-05 * ctWRe * ctWRe
                    + -6.845532348889237e-06 * cHG * cHGtil + -3.0658735341057295e-05 * cHG * ctGIm
                    + -0.0038039141269859 * cHG * ctGRe + 4.037792277630792e-06 * cHG * ctHIm
                    + -0.000253969240318 * cHG * ctHRe + 5.585076583089488e-06 * cHG * ctWIm
                    + -8.115494553287108e-07 * cHG * ctWRe + -0.0025610325162679 * cHGtil * ctGIm
                    + -4.903421501601487e-06 * cHGtil * ctGRe + 0.0001789011931272 * cHGtil * ctHIm
                    + 1.4482656857522924e-06 * cHGtil * ctHRe + 2.1733719725858247e-06 * cHGtil * ctWIm
                    + 8.391042189262454e-06 * cHGtil * ctWRe + -1.5487075151731095e-06 * ctGIm * ctGRe
                    + -1.6017969617732204e-05 * ctGIm * ctHIm + 1.1205379927686776e-06 * ctGIm * ctHRe
                    + 1.9900791938254336e-06 * ctGIm * ctWIm + -5.365572614642063e-06 * ctGIm * ctWRe
                    + -1.8151564085000028e-06 * ctGRe * ctHIm + 0.0002528327953381 * ctGRe * ctHRe
                    + -4.974839844277338e-06 * ctGRe * ctWIm + 9.916092663458942e-06 * ctGRe * ctWRe
                    + 1.123734370417351e-07 * ctHIm * ctHRe + -1.931109460764816e-06 * ctHIm * ctWIm
                    + 1.5128957953081112e-07 * ctHIm * ctWRe + -1.3678888535800165e-07 * ctHRe * ctWIm
                    + -2.746731171583885e-06 * ctHRe * ctWRe + 6.285203592998334e-07 * ctWIm * ctWRe
                    -7.332564342765746e-06 * ctBRe + 1.488317517997495e-05 * ctBRe * ctBRe + 5.322544053956713e-06 * ctBIm * ctBIm
                    )*(SM_ttH_bin_860_1270/ttH_bin_860_1270_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

            }
            else{
                total =  SM_ttH_bin_860_1270 + (0.0035482569627141*cHG +
                        -0.0041580771412648*ctGRe + -0.0003266623947329*ctHRe
                        -7.332564342765746e-06 * ctBRe-1.1085493448650952e-05 * ctWRe
                        )*(SM_ttH_bin_860_1270/ttH_bin_860_1270_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
        }
    
    } else if(b_min == 1270 && b_max == 2500){
        
        double SM_ttH_bin_1270_2500 = SM.getOptionalParameter("SM_sigma_ttH_diff_LO_mtth_bin_1270_2500");
        double ttH_bin_1270_2500_madgraph_LO = 0.0002602108651872;//pb/GeV
        double total;
 

        if(flag_Quadratic){
                total =  SM_ttH_bin_1270_2500 + (0.0010427269178651 * cHG + -9.856161004434538e-07 * cHGtil
                    + 4.622973593784718e-06 * ctGIm + -0.0007284069850997 * ctGRe
                    + 2.406669648430393e-07 * ctHIm + -3.129120993544568e-05 * ctHRe
                    + 1.558576659940747e-06 * ctWIm + -5.950496569095965e-08 * ctWRe
                    + 0.0011892767595462 * cHG * cHG + 0.0011046291764241 * cHGtil * cHGtil
                    + 0.0011424484585363 * ctGIm * ctGIm + 0.0010765623784725 * ctGRe * ctGRe
                    + 1.175058315219801e-06 * ctHIm * ctHIm + 1.2422431746160378e-06 * ctHRe * ctHRe
                    + 0.0002574437383637 * ctWIm * ctWIm + 3.634533253690786e-05 * ctWRe * ctWRe
                    + -3.889698786486935e-06 * cHG * cHGtil + 1.663871428932151e-05 * cHG * ctGIm
                    + -0.0009516653801491 * cHG * ctGRe + -8.764646223761842e-07 * cHG * ctHIm
                    + -0.0001115383166857 * cHG * ctHRe + -2.849926045572032e-06 * cHG * ctWIm
                    + -3.1134703074364467e-06 * cHG * ctWRe + -0.0006835996252276 * cHGtil * ctGIm
                    + 2.6865465597297744e-06 * cHGtil * ctGRe + 1.839047460509835e-05 * cHGtil * ctHIm
                    + -9.444369397582044e-07 * cHGtil * ctHRe + -4.2031651278482206e-06 * cHGtil * ctWIm
                    + -1.902334553772218e-06 * cHGtil * ctWRe + -3.7867570888689487e-06 * ctGIm * ctGRe
                    + 4.67918943547474e-06 * ctGIm * ctHIm + -3.3779355185015254e-08 * ctGIm * ctHRe
                    + -9.90107948205884e-06 * ctGIm * ctWIm + 3.2439146454010803e-06 * ctGIm * ctWRe
                    + -1.4688339893080404e-08 * ctGRe * ctHIm + 2.912940987744028e-05 * ctGRe * ctHRe
                    + -2.044738314961812e-06 * ctGRe * ctWIm + 1.1890602603092937e-06 * ctGRe * ctWRe
                    + 8.774074439221735e-09 * ctHIm * ctHRe + -6.938904354794587e-07 * ctHIm * ctWIm
                    + 9.554777502529532e-08 * ctHIm * ctWRe + -3.3366971401629874e-08 * ctHRe * ctWIm
                    + -9.04292431679249e-07 * ctHRe * ctWRe + -1.1428837119814367e-06 * ctWIm * ctWRe
                    -1.1656806445381473e-06 * ctBRe + 4.4686841467467264e-06 * ctBRe * ctBRe+ 2.6432943491666583e-06 * ctBIm * ctBIm
                        )*(SM_ttH_bin_1270_2500/ttH_bin_1270_2500_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        else{
                total =  SM_ttH_bin_1270_2500 + (0.0010427269178651*cHG +
                        -0.0007284069850997*ctGRe + -3.129120993544568e-05*ctHRe
                        -1.1656806445381473e-06 * ctBRe -5.9504965690959644e-08 * ctWRe
                        )*(SM_ttH_bin_1270_2500/ttH_bin_1270_2500_madgraph_LO);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
        
        }
    
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct binning for sigma_ttH_diff_LO_mtth. \n");
    }

}
    
    
    
    
    


Asymmetry_trip_prod_pt_pe_pp_ttH::Asymmetry_trip_prod_pt_pe_pp_ttH(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

//    setParametersForObservable(make_vector<std::string>() << "SM_Asymmetry_trip_prod_pt_pe_pp_ttH");
    
}

double Asymmetry_trip_prod_pt_pe_pp_ttH::computeThValue()
{
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiG();
    double cHGtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiGtil();
    //double cHW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiW();
    //double cHWtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiWtil();
    double ctHRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double ctHIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphiIm();
    double ctGRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double ctGIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tGIm();
    double ctWRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double ctWIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tWIm();
    double ctBRe = 0.;
    double ctBIm = 0.;
    
    
        //Expressions are written in the basis of SMEFTsim
    //let's change to the basis of dim6top (same as
    //SMEFT@NLO except for the gs on CtG). These
    //models define the covariance derivative 
    //with a different convention
    
    if(flag_LHC_WG_Basis){
        double ctZRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
        double ctZIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZIm();
        double sw2 = 0.22305;
        double cw2 = 1 - sw2;
        double tw2 = sw2/cw2;
        double sw  = sqrt(sw2);
        double tw  = sqrt(tw2);
        ctBRe = ctZRe/sw - ctWRe/tw;
        ctBIm = ctZIm/sw - ctWIm/tw;
        ctWRe = -ctWRe;
        ctWIm = -ctWIm;
        ctGRe = -ctGRe;
        ctGIm = -ctGIm;
    }
    else{
        ctBRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        ctBIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tBIm();
    }
    
    
    
    //double SM_Asymmetry_Dazi_ord = SM.getOptionalParameter("SM_Asymmetry_Dazi_ord");
    
    double trip_prod_neg_Mad = 2.5127944405346114;
    double trip_prod_pos_Mad = 2.520483758548904;
    


    if(flag_Quadratic){


    
            double trip_prod_neg_NP = 2.640183516559485 * cHG + -0.291135613883015 * cHGtil
                + 0.3678572986205774 * ctGIm + -3.087690836545025 * ctGRe
                + 0.0094009278222139 * ctHIm + -0.3056879017920962 * ctHRe
                + 0.0028779039125397 * ctWIm + -0.0090384597387344 * ctWRe
                + 1.9428373183076488 * cHG * cHG + 1.8884134166653317 * cHGtil * cHGtil
                + 2.470902330930285 * ctGIm * ctGIm + 2.5390721751356464 * ctGRe * ctGRe
                + 0.004393313833058 * ctHIm * ctHIm + 0.0098526634558913 * ctHRe * ctHRe
                + 0.3066205731868689 * ctWIm * ctWIm + 0.0560668334542675 * ctWRe * ctWRe
                + -0.0611273430332702 * cHG * cHGtil + -0.1562901590147846 * cHG * ctGIm
                + -2.3280072653525665 * cHG * ctGRe + -0.0332829035749282 * cHG * ctHIm
                + -0.1922967483887002 * cHG * ctHRe + 0.001261653307765 * cHG * ctWIm
                + -0.0028384867252228 * cHG * ctWRe + -1.4936120158370252 * cHGtil * ctGIm
                + 0.2195689708810037 * cHGtil * ctGRe + 0.1284110619808202 * cHGtil * ctHIm
                + 0.0435618768488701 * cHGtil * ctHRe + -0.0005958003085194 * cHGtil * ctWIm
                + 0.0010590191160643 * cHGtil * ctWRe + 0.0180277458513596 * ctGIm * ctGRe
                + -0.0037356858757553 * ctGIm * ctHIm + -0.0251806671367735 * ctGIm * ctHRe
                + 0.0145409140622336 * ctGIm * ctWIm + 0.0006672516568214 * ctGIm * ctWRe
                + -0.0129283651774673 * ctGRe * ctHIm + 0.2122410869542897 * ctGRe * ctHRe
                + 0.0148165894863743 * ctGRe * ctWIm + -0.0048128893054671 * ctGRe * ctWRe
                + -0.0003906827817827 * ctHIm * ctHRe + -0.0012818830236245 * ctHIm * ctWIm
                + 0.0001845529183801 * ctHIm * ctWRe + -0.0001007080896056 * ctHRe * ctWIm
                + -0.0016692696078997 * ctHRe * ctWRe + 0.0110862621259655 * ctWIm * ctWRe
            -0.004354637597045269 * ctBRe + 0.009775150516485295 * ctBRe * ctBRe+ -0.00028492026513413093 * ctBIm + 0.0037313612606456603 * ctBIm * ctBIm
            ;
            double trip_prod_pos_NP = 2.6328115532520737 * cHG + 0.2910052623622711 * cHGtil
                + -0.363739875902877 * ctGIm + -3.096941761896972 * ctGRe
                + -0.0094759620726647 * ctHIm + -0.305830247123692 * ctHRe
                + -0.0020413671276346 * ctWIm + -0.0080240264123283 * ctWRe
                + 1.938529650870392 * cHG * cHG + 1.884815095098858 * cHGtil * cHGtil
                + 2.4785195653397807 * ctGIm * ctGIm + 2.5282292515897424 * ctGRe * ctGRe
                + 0.0044308712361479 * ctHIm * ctHIm + 0.0097609071828327 * ctHRe * ctHRe
                + 0.3074217267298444 * ctWIm * ctWIm + 0.0555200231929575 * ctWRe * ctWRe
                + 0.0616850233108645 * cHG * cHGtil + 0.1686557442326578 * cHG * ctGIm
                + -2.3329631336855035 * cHG * ctGRe + 0.0333533207345481 * cHG * ctHIm
                + -0.1980141552978979 * cHG * ctHRe + 0.0001794692513655 * cHG * ctWIm
                + 0.0029579718310278 * cHG * ctWRe + -1.4707955683599927 * cHGtil * ctGIm
                + -0.2225681716920635 * cHGtil * ctGRe + 0.1264439339654254 * cHGtil * ctHIm
                + -0.0440141389916522 * cHGtil * ctHRe + 0.0013869829103071 * cHGtil * ctWIm
                + -0.0010721150335592 * cHGtil * ctWRe + -0.020961743091158 * ctGIm * ctGRe
                + -0.0028677842945724 * ctGIm * ctHIm + 0.0249953324058435 * ctGIm * ctHRe
                + -0.0142367301713332 * ctGIm * ctWIm + -8.440301745165969e-05 * ctGIm * ctWRe
                + 0.0124684222616443 * ctGRe * ctHIm + 0.2119156676645677 * ctGRe * ctHRe
                + -0.0146267046439642 * ctGRe * ctWIm + 0.0073226429336149 * ctGRe * ctWRe
                + 0.0003720140582797 * ctHIm * ctHRe + -0.0007259142884621 * ctHIm * ctWIm
                + -9.280611641033528e-05 * ctHIm * ctWRe + 0.0001350325322032 * ctHRe * ctWIm
                + -0.0005892796643464 * ctHRe * ctWRe + -0.0096859099971627 * ctWIm * ctWRe
            -0.004906596203596264 * ctBRe + 0.009788072415787982 * ctBRe * ctBRe+ 0.0002797884743307305 * ctBIm + 0.0037358256356170294 * ctBIm * ctBIm
            ;
    
    
            //double num = (trip_prod_pos_Mad + trip_prod_pos_NP) - (trip_prod_neg_Mad + trip_prod_neg_NP);
            //double deno = (trip_prod_pos_Mad + trip_prod_pos_NP) + (trip_prod_neg_Mad + trip_prod_neg_NP);
            
            
            //We set the SM difference to zero
            double num = (trip_prod_pos_NP) - (trip_prod_neg_NP);
            double deno = (trip_prod_pos_Mad + trip_prod_pos_NP) + (trip_prod_neg_Mad + trip_prod_neg_NP);
    
            double total = num/deno;
            return total;
            
    }
    else{

            
    
            double trip_prod_neg_NP = 2.640095238800265*cHG + -0.291215708923054*cHGtil + 0.3672351161865572*ctGIm +
            -3.086551544820942*ctGRe + 0.0093896874318005*ctHIm + -0.3056973190252862*ctHRe
            -0.004354637597045269 * ctBRe+ -0.00028492026513413093 * ctBIm
            -0.009038459738734433 * ctWRe+ 0.002877903912539781 * ctWIm
            ;
            double trip_prod_pos_NP = 2.6338438330285854*cHG + 0.2910853340659681*cHGtil + -0.3631167628767298*ctGIm +
            -3.0994782998687884*ctGRe + -0.0094647329097138*ctHIm + -0.3059123320650544*ctHRe
            -0.004906596203596264 * ctBRe+ 0.0002797884743307305 * ctBIm
            -0.008024026412328383 * ctWRe+ -0.0020413671276346412 * ctWIm
            ;
    
    
            //double num = (trip_prod_pos_Mad + trip_prod_pos_NP) - (trip_prod_neg_Mad + trip_prod_neg_NP);
            //double deno = (trip_prod_pos_Mad + trip_prod_pos_NP) + (trip_prod_neg_Mad + trip_prod_neg_NP);
            
            
            //We set the SM difference to zero
            double num = (trip_prod_pos_NP) - (trip_prod_neg_NP);
            double deno = (trip_prod_pos_Mad + trip_prod_pos_NP) + (trip_prod_neg_Mad + trip_prod_neg_NP);
            

            
            double total = num/deno;
            return total;
            
    }
    
}






Asymmetry_sign_trip_prod_pe_pp_ttH::Asymmetry_sign_trip_prod_pe_pp_ttH(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

//    setParametersForObservable(make_vector<std::string>() << "SM_Asymmetry_sign_trip_prod_pe_pp_ttH");
    
}

double Asymmetry_sign_trip_prod_pe_pp_ttH::computeThValue()
{
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    double cHG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiG();
    double cHGtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiGtil();
    //double cHW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiW();
    //double cHWtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiWtil();
    double ctHRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double ctHIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphiIm();
    double ctGRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double ctGIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tGIm();
    double ctWRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double ctWIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tWIm();
    double ctBRe = 0.;
    double ctBIm = 0.;
    
    //Expressions are written in the basis of SMEFTsim
    //let's change to the basis of dim6top (same as
    //SMEFT@NLO except for the gs on CtG). These
    //models define the covariance derivative 
    //with a different convention
    
    if(flag_LHC_WG_Basis){
        double ctZRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
        double ctZIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZIm();
        double sw2 = 0.22305;
        double cw2 = 1 - sw2;
        double tw2 = sw2/cw2;
        double sw  = sqrt(sw2);
        double tw  = sqrt(tw2);
        ctBRe = ctZRe/sw - ctWRe/tw;
        ctBIm = ctZIm/sw - ctWIm/tw;
        ctWRe = -ctWRe;
        ctWIm = -ctWIm;
        ctGRe = -ctGRe;
        ctGIm = -ctGIm;
    }
    else{
        ctBRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        ctBIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tBIm();
    }
    
    
    
    //double SM_Asymmetry_Dazi_ord = SM.getOptionalParameter("SM_Asymmetry_Dazi_ord");
    
    double sign_trip_prod_pe_pp_neg_Mad = 2.5199987657363723;
    double sign_trip_prod_pe_pp_pos_Mad = 2.512378708453129;
    
    

    if(flag_Quadratic){
            
            
            
    
            double sign_trip_prod_pe_pp_neg_NP = 2.642738820737528 * cHG + 0.1857689148975252 * cHGtil
                + -0.2201130538105204 * ctGIm + -3.1022541249960875 * ctGRe
                + -0.0059790521485915 * ctHIm + -0.3068221716991218 * ctHRe
                + -0.0015026106634554 * ctWIm + -0.0083629306863192 * ctWRe
                + 1.943560394107679 * cHG * cHG + 1.8868727613590808 * cHGtil * cHGtil
                + 2.4782475378888056 * ctGIm * ctGIm + 2.537205329210237 * ctGRe * ctGRe
                + 0.0044286014949246 * ctHIm * ctHIm + 0.0099118147633419 * ctHRe * ctHRe
                + 0.3067008919750776 * ctWIm * ctWIm + 0.0556639803507113 * ctWRe * ctWRe
                + 0.0508984896879067 * cHG * cHGtil + 0.1476151527744092 * cHG * ctGIm
                + -2.340190729813459 * cHG * ctGRe + 0.0299484591443408 * cHG * ctHIm
                + -0.1959685579715413 * cHG * ctHRe + -0.0009799542926578 * cHG * ctWIm
                + 0.0013594420154272 * cHG * ctWRe + -1.4761384971318905 * cHGtil * ctGIm
                + -0.1713637905247774 * cHGtil * ctGRe + 0.1279944459298787 * cHGtil * ctHIm
                + -0.0292542423455988 * cHGtil * ctHRe + 0.0008465992940455 * cHGtil * ctWIm
                + 0.0003448470109748 * cHGtil * ctWRe + -0.0261746428921166 * ctGIm * ctGRe
                + -0.0032607955481276 * ctGIm * ctHIm + 0.0140985825882586 * ctGIm * ctHRe
                + -0.0085553818589635 * ctGIm * ctWIm + -0.0012864686136748 * ctGIm * ctWRe
                + 0.0074007867049598 * ctGRe * ctHIm + 0.2117097917766713 * ctGRe * ctHRe
                + -0.0207629155173062 * ctGRe * ctWIm + 0.0067564123668787 * ctGRe * ctWRe
                + 6.818562664689942e-05 * ctHIm * ctHRe + -0.0007987813897161 * ctHIm * ctWIm
                + 4.501543383594875e-05 * ctHIm * ctWRe + 0.0002535284543679 * ctHRe * ctWIm
                + -0.0008213347686312 * ctHRe * ctWRe + -0.0081842768917958 * ctWIm * ctWRe
            -0.004607323137037447 * ctBRe + 0.009777333000546243 * ctBRe * ctBRe+ 0.00016620851227672206 * ctBIm + 0.003731326965602036 * ctBIm * ctBIm
            ;
            double sign_trip_prod_pe_pp_pos_NP = 2.630256251912717 * cHG + -0.185899266418269 * cHGtil
                + 0.2242304765282196 * ctGIm + -3.082378464959626 * ctGRe
                + 0.0059040178981408 * ctHIm + -0.3046959772166663 * ctHRe
                + 0.002339089715825 * ctWIm + -0.0086983779187188 * ctWRe
                + 1.937806573270836 * cHG * cHG + 1.886355750337498 * cHGtil * cHGtil
                + 2.471174359351984 * ctGIm * ctGIm + 2.5300960959070182 * ctGRe * ctGRe
                + 0.0043955835742812 * ctHIm * ctHIm + 0.0097017558753822 * ctHRe * ctHRe
                + 0.307299030584504 * ctWIm * ctWIm + 0.0559151752699298 * ctWRe * ctWRe
                + -0.0503411686164696 * cHG * cHGtil + -0.1352535186964854 * cHG * ctGIm
                + -2.317996993087569 * cHG * ctGRe + -0.0298780736340808 * cHG * ctHIm
                + -0.1941871296364016 * cHG * ctHRe + 0.0024201523602079 * cHG * ctWIm
                + -0.0012399863287341 * cHG * ctWRe + -1.4878249507037815 * cHGtil * ctGIm
                + 0.1683650160418539 * cHGtil * ctGRe + 0.1268273709363991 * cHGtil * ctHIm
                + 0.0288020394197951 * cHGtil * ctHRe + -5.606197840776873e-05 * cHGtil * ctWIm
                + -0.0003579384133383 * cHGtil * ctWRe + 0.0232417256197431 * ctGIm * ctGRe
                + -0.0033420150171774 * ctGIm * ctHIm + -0.0142839014035242 * ctGIm * ctHRe
                + 0.0088594432873936 * ctGIm * ctWIm + 0.00186901482171 * ctGIm * ctWRe
                + -0.0078606428282776 * ctGRe * ctHIm + 0.2123883426407318 * ctGRe * ctHRe
                + 0.0209527371632488 * ctGRe * ctWIm + -0.004247300263279 * ctGRe * ctWRe
                + -8.685221915353248e-05 * ctHIm * ctHRe + -0.001207429352545 * ctHIm * ctWIm
                + 4.671135156023021e-05 * ctHIm * ctWRe + -0.000219233464493 * ctHRe * ctWIm
                + -0.0014368897467772 * ctHRe * ctWRe + 0.0095836863532764 * ctWIm * ctWRe
            -0.0046526081687880555 * ctBRe + 0.0097831385708423 * ctBRe * ctBRe+ -0.00017133958134900318 * ctBIm + 0.0037348097496134464 * ctBIm * ctBIm
            ;
            
    
            //double num = (sign_trip_prod_pe_pp_pos_Mad + sign_trip_prod_pe_pp_pos_NP) - (sign_trip_prod_pe_pp_neg_Mad + sign_trip_prod_pe_pp_neg_NP);
            //double deno = (sign_trip_prod_pe_pp_pos_Mad + sign_trip_prod_pe_pp_pos_NP) + (sign_trip_prod_pe_pp_neg_Mad + sign_trip_prod_pe_pp_neg_NP);
            
            //We set the SM difference to zero
            double num = (sign_trip_prod_pe_pp_pos_NP) - (sign_trip_prod_pe_pp_neg_NP);
            double deno = (sign_trip_prod_pe_pp_pos_Mad + sign_trip_prod_pe_pp_pos_NP) + (sign_trip_prod_pe_pp_neg_Mad + sign_trip_prod_pe_pp_neg_NP);

            double total = num/deno;
            return total;
            
        }
        else{

            
    
            double sign_trip_prod_pe_pp_neg_NP = 2.642738924499983*cHG + 0.1857689223794195*cHGtil + -0.2201130894378797*ctGIm +
            -3.102254615722175*ctGRe + -0.0059790531590874*ctHIm + -0.3068222235539151*ctHRe
            -0.004607323137037447 * ctBRe+ 0.00016620851227672206 * ctBIm
            -0.008362930686319298 * ctWRe+ -0.0015026106634554587 * ctWIm
            ;
            double sign_trip_prod_pe_pp_pos_NP = 2.6302563644690693*cHG + -0.1858992739054133*cHGtil + 0.2242305128220224*ctGIm + 
            -3.0823789676020823*ctGRe + 0.0059040188959551*ctHIm + -0.30469602871212*ctHRe
            -0.0046526081687880555 * ctBRe+ -0.00017133958134900318 * ctBIm
            -0.008698377918718886 * ctWRe+ 0.002339089715825096 * ctWIm
            ;
    
    
            //double num = (sign_trip_prod_pe_pp_pos_Mad + sign_trip_prod_pe_pp_pos_NP) - (sign_trip_prod_pe_pp_neg_Mad + sign_trip_prod_pe_pp_neg_NP);
            //double deno = (sign_trip_prod_pe_pp_pos_Mad + sign_trip_prod_pe_pp_pos_NP) + (sign_trip_prod_pe_pp_neg_Mad + sign_trip_prod_pe_pp_neg_NP);
            
            //We set the SM difference to zero
            double num = (sign_trip_prod_pe_pp_pos_NP) - (sign_trip_prod_pe_pp_neg_NP);
            double deno = (sign_trip_prod_pe_pp_pos_Mad + sign_trip_prod_pe_pp_pos_NP) + (sign_trip_prod_pe_pp_neg_Mad + sign_trip_prod_pe_pp_neg_NP);

            
            double total = num/deno;
            return total;
            
    }
    
}










Asymmetry_cos_je_tHj::Asymmetry_cos_je_tHj(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

//    setParametersForObservable(make_vector<std::string>() << "SM_Asymmetry_cos_je_tHj");

}

double Asymmetry_cos_je_tHj::computeThValue()
{
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    //double cHG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiG();
    //double cHGtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiGtil();
    double cHW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiW();
    double cHWtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiWtil();
    double ctHRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double ctHIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphiIm();
    //double ctGRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    //double ctGIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tGIm();
    double ctWRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double ctWIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tWIm();
    
    
    
    //Expressions are written in the basis of SMEFTsim
    //let's change to the basis of dim6top (same as
    //SMEFT@NLO except for the gs on CtG). These
    //models define the covariance derivative 
    //with a different convention
    
    if(flag_LHC_WG_Basis){
//        double ctZRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
//        double ctZIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZIm();
//        sw2 = 0.22305;
//        cw2 = 1 - sw2;
//        tw2 = sw2/cw2;
//        sw  = sqrt(sw2);
//        tw  = sqrt(tw2);
//        ctBRe = ctZRe/sw - ctWRe/tw;
//        ctBIm = ctZIm/sw - ctWIm/tw;
        ctWRe = -ctWRe;
        ctWIm = -ctWIm;
//        ctGRe = -ctGRe;
//        ctGIm = -ctGIm;
    }
//    else{
//        double ctBRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
//        double ctBIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tBIm();
//    }
    
    
    //double SM_Asymmetry_cos_je_tHj = SM.getOptionalParameter("SM_Asymmetry_cos_je_tHj");
    
    double cos_je_neg_Mad = 1.802755193194093;
    double cos_je_pos_Mad = 4.185564082240573;
    

    if(flag_Quadratic){


    
            double cos_je_neg_NP = 0.3648966515528267*cHW + 0.001727621412455*cHWtil + 0.0002137853404453*ctHIm + 
            -0.0682956856995048*ctHRe + -0.0030877434371623*ctWIm + -1.329003891718189*ctWRe + 
            0.1866736021132992*cHW*cHW + 0.1452346608094477*cHWtil*cHWtil + 0.0116273887361596*ctHIm*ctHIm + 
            0.0231503175651335*ctHRe*ctHRe + 3.83451885101712*ctWIm*ctWIm + 3.8332542149381097*ctWRe*ctWRe + 
            -0.0007900978367795*cHW*cHWtil + 1.0544889795918722e-05*cHW*ctHIm + -0.0327623256577681*cHW*ctHRe +
            0.0007437915907154*cHW*ctWIm + -0.6726490763976125*cHW*ctWRe + -0.0049881868923777*cHWtil*ctHIm + 
            -0.0005470517961899*cHWtil*ctHRe + -0.4769755806909019*cHWtil*ctWIm + -0.0006823264112964*cHWtil*ctWRe + 
            0.0002621344978575*ctHIm*ctHRe + 0.0748309927213606*ctHIm*ctWIm + -0.0018785485471096*ctHIm*ctWRe +
            -0.0018820801261675*ctHRe*ctWIm + 0.1058363336301706*ctHRe*ctWRe + -0.0047330354036098*ctWIm*ctWRe;
            double cos_je_pos_NP = 0.5847019499353336*cHW + -0.0007184422927661*cHWtil + -0.0003465092977812*ctHIm +
            -0.0887773510001904*ctHRe + 0.0062306067229668*ctWIm + -0.9345195813265864*ctWRe + 
            0.427307590720827*cHW*cHW + 0.3073774610006161*cHWtil*cHWtil + 0.0090770914570092*ctHIm*ctHIm + 
            0.0430639890806705*ctHRe*ctHRe + 1.5822244524211948*ctWIm*ctWIm + 1.583096674934071*ctWRe*ctWRe + 
            0.0005480173007809*cHW*cHWtil + 4.277811723747693e-05*cHW*ctHIm + -0.0874650266649605*cHW*ctHRe + 
            -0.0024641068137472*cHW*ctWIm + -0.6056247049583481*cHW*ctWRe + -0.0108252574635077*cHWtil*ctHIm + 
            0.0006241202527166*cHWtil*ctHRe + -0.4081977848482742*cHWtil*ctWIm + 0.0020255764722883*cHWtil*ctWRe + 
            -0.0001016003457396*ctHIm*ctHRe + 0.0744221969307464*ctHIm*ctWIm + 0.0015065993293668*ctHIm*ctWRe + 
            -0.0001967249821085*ctHRe*ctWIm + 0.1432442378938073*ctHRe*ctWRe + -0.0018638514243001*ctWIm*ctWRe;
            
            
            
            double num = (cos_je_pos_Mad + cos_je_pos_NP) - (cos_je_neg_Mad + cos_je_neg_NP);
            double deno = (cos_je_pos_Mad + cos_je_pos_NP) + (cos_je_neg_Mad + cos_je_neg_NP);
    
            
            double total = num/deno;
            return total;
            
    }
    else{
    
            double cos_je_neg_NP = 0.3648966515528267*cHW + 0.001727621412455*cHWtil + 0.0002137853404453*ctHIm + 
            -0.0682956856995048*ctHRe + -0.0030877434371623*ctWIm + -1.329003891718189*ctWRe;
            double cos_je_pos_NP = 0.5847019499353336*cHW + -0.0007184422927661*cHWtil + -0.0003465092977812*ctHIm +
            -0.0887773510001904*ctHRe + 0.0062306067229668*ctWIm + -0.9345195813265864*ctWRe;
    
    
            double num = (cos_je_pos_Mad + cos_je_pos_NP) - (cos_je_neg_Mad + cos_je_neg_NP);
            double deno = (cos_je_pos_Mad + cos_je_pos_NP) + (cos_je_neg_Mad + cos_je_neg_NP);
    
            double total = num/deno;
            return total;
            
    
    }  
}





Asymmetry_cos_se_tHj::Asymmetry_cos_se_tHj(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

//    setParametersForObservable(make_vector<std::string>() << "SM_Asymmetry_cos_se_tHj");

}

double Asymmetry_cos_se_tHj::computeThValue()
{
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    //double cHG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiG();
    //double cHGtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiGtil();
    double cHW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiW();
    double cHWtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiWtil();
    double ctHRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double ctHIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphiIm();
    //double ctGRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    //double ctGIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tGIm();
    double ctWRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double ctWIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tWIm();
    
    
        //Expressions are written in the basis of SMEFTsim
    //let's change to the basis of dim6top (same as
    //SMEFT@NLO except for the gs on CtG). These
    //models define the covariance derivative 
    //with a different convention
    
    if(flag_LHC_WG_Basis){
//        double ctZRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
//        double ctZIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZIm();
//        sw2 = 0.22305;
//        cw2 = 1 - sw2;
//        tw2 = sw2/cw2;
//        sw  = sqrt(sw2);
//        tw  = sqrt(tw2);
//        ctBRe = ctZRe/sw - ctWRe/tw;
//        ctBIm = ctZIm/sw - ctWIm/tw;
        ctWRe = -ctWRe;
        ctWIm = -ctWIm;
//        ctGRe = -ctGRe;
//        ctGIm = -ctGIm;
    }
//    else{
//        double ctBRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
//        double ctBIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tBIm();
//    }
    
    
    //double SM_Asymmetry_cos_se_tHj = SM.getOptionalParameter("SM_Asymmetry_cos_se_tHj");
    
    double cos_se_neg_Mad = 2.990807050970236;
    double cos_se_pos_Mad = 2.9975122244644314;
    

    if(flag_Quadratic){


    
            double cos_se_neg_NP = 0.4728914373547805*cHW + 0.0118407416464167*cHWtil + -0.0033338074845438*ctHIm +
            -0.0804656940482105*ctHRe + -0.689834314040347*ctWIm + -1.132189187920139*ctWRe + 
            0.3068652120228926*cHW*cHW + 0.225623165775055*cHWtil*cHWtil + 0.0104250772428093*ctHIm*ctHIm + 
            0.0338402455559234*ctHRe*ctHRe + 2.71088765118971*ctWIm*ctWIm + 2.7161130000666813*ctWRe*ctWRe + 
            -0.0055312276935633*cHW*cHWtil + 0.0096577534668971*cHW*ctHIm + -0.0594580614259007*cHW*ctHRe + 
            -0.1979153488691058*cHW*ctWIm + -0.6372052845724958*cHW*ctWRe + -0.0078440619072392*cHWtil*ctHIm +
            -0.0031289075947225*cHWtil*ctHRe + -0.4433548116003082*cHWtil*ctWIm + 0.2461428388687654*cHWtil*ctWRe + 
            -0.0067462516967434*ctHIm*ctHRe + 0.0758507838776238*ctHIm*ctWIm + -0.008892021283742*ctHIm*ctWRe + 
            -0.0085734101641831*ctHRe*ctWIm + 0.1236903103660209*ctHRe*ctWRe + -0.0123076830730909*ctWIm*ctWRe;
            double cos_se_pos_NP = 0.4767071641321322*cHW + -0.0108315625270247*cHWtil + 0.0032010835270088*ctHIm + 
            -0.0766073426520603*ctHRe + 0.6929771773260096*ctWIm + -1.131334285096125*ctWRe + 
            0.3071159808107537*cHW*cHW + 0.2269889560350069*cHWtil*cHWtil + 0.0102794029503121*ctHIm*ctHIm + 
            0.0323740610901135*ctHRe*ctHRe + 2.705855652237456*ctWIm*ctWIm + 2.7002378898037698*ctWRe*ctWRe +
            0.0052891474393248*cHW*cHWtil + -0.00960443047912*cHW*ctHIm + -0.0607692401807437*cHW*ctHRe + 
            0.1961950351073473*cHW*ctWIm + -0.6410678813860773*cHW*ctWRe + -0.0079693317308205*cHWtil*ctHIm + 
            0.0032059760035094*cHWtil*ctHRe + -0.4418189733437092*cHWtil*ctWIm + -0.2447995889997185*cHWtil*ctWRe + 
            0.006906785245703*ctHIm*ctHRe + 0.0734022421909097*ctHIm*ctWIm + 0.0085200721419787*ctHIm*ctWRe + 
            0.0064946056775615*ctHRe*ctWIm + 0.1253901576970864*ctHRe*ctWRe + 0.0057108033786834*ctWIm*ctWRe;
    
    
            //double num = (cos_se_pos_Mad + cos_se_pos_NP) - (cos_se_neg_Mad + cos_se_neg_NP);
            //double deno = (cos_se_pos_Mad + cos_se_pos_NP) + (cos_se_neg_Mad + cos_se_neg_NP);
    
            
            //We set the SM difference to zero
            double num = (cos_se_pos_NP) - (cos_se_neg_NP);
            double deno = (cos_se_pos_Mad + cos_se_pos_NP) + (cos_se_neg_Mad + cos_se_neg_NP);
    
            double total = num/deno;
            return total;
            
    }
    else{

            
    
            double cos_se_neg_NP = 0.4728914373547805*cHW + 0.0118407416464167*cHWtil + -0.0033338074845438*ctHIm +
            -0.0804656940482105*ctHRe + -0.689834314040347*ctWIm + -1.132189187920139*ctWRe;
            double cos_se_pos_NP = 0.4767071641321322*cHW + -0.0108315625270247*cHWtil + 0.0032010835270088*ctHIm + 
            -0.0766073426520603*ctHRe + 0.6929771773260096*ctWIm + -1.131334285096125*ctWRe;
    
    
            //double num = (cos_se_pos_Mad + cos_se_pos_NP) - (cos_se_neg_Mad + cos_se_neg_NP);
            //double deno = (cos_se_pos_Mad + cos_se_pos_NP) + (cos_se_neg_Mad + cos_se_neg_NP);
    
            
            //We set the SM difference to zero
            double num = (cos_se_pos_NP) - (cos_se_neg_NP);
            double deno = (cos_se_pos_Mad + cos_se_pos_NP) + (cos_se_neg_Mad + cos_se_neg_NP);
    
            double total = num/deno;
            return total;
            
    }  
}





Asymmetry_cos_ye_tHj::Asymmetry_cos_ye_tHj(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

//    setParametersForObservable(make_vector<std::string>() << "SM_Asymmetry_cos_ye_tHj");

}

double Asymmetry_cos_ye_tHj::computeThValue()
{
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    //double cHG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiG();
    //double cHGtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiGtil();
    double cHW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiW();
    double cHWtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiWtil();
    double ctHRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double ctHIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphiIm();
    //double ctGRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    //double ctGIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tGIm();
    double ctWRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double ctWIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tWIm();
    
    
    //Expressions are written in the basis of SMEFTsim
    //let's change to the basis of dim6top (same as
    //SMEFT@NLO except for the gs on CtG). These
    //models define the covariance derivative 
    //with a different convention
    
    if(flag_LHC_WG_Basis){
//        double ctZRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
//        double ctZIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZIm();
//        sw2 = 0.22305;
//        cw2 = 1 - sw2;
//        tw2 = sw2/cw2;
//        sw  = sqrt(sw2);
//        tw  = sqrt(tw2);
//        ctBRe = ctZRe/sw - ctWRe/tw;
//        ctBIm = ctZIm/sw - ctWIm/tw;
        ctWRe = -ctWRe;
        ctWIm = -ctWIm;
//        ctGRe = -ctGRe;
//        ctGIm = -ctGIm;
    }
//    else{
//        double ctBRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
//        double ctBIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tBIm();
//    }
    
    //double SM_Asymmetry_cos_ye_tHj = SM.getOptionalParameter("SM_Asymmetry_cos_ye_tHj");
    
    double cos_ye_neg_Mad = 2.9854504039888377;
    double cos_ye_pos_Mad = 3.002868871445828;
    

    if(flag_Quadratic){


    
            double cos_ye_neg_NP = 0.4728811316749511*cHW + -0.0743571487036044*cHWtil + 0.0134093106580527*ctHIm +
            -0.0788765537664126*ctHRe + -0.5965260044476208*ctWIm + -1.133879975967292*ctWRe + 
            0.3066300756612584*cHW*cHW + 0.2265944646346798*cHWtil*cHWtil + 0.0105197938347638*ctHIm*ctHIm +
            0.0348166658751635*ctHRe*ctHRe + 2.7103779842736944*ctWIm*ctWIm + 2.7128731828613617*ctWRe*ctWRe + 
            0.0067989490892239*cHW*cHWtil + 0.0062514079134364*cHW*ctHIm + -0.0596622342861595*cHW*ctHRe + 
            -0.6504445487426931*cHW*ctWIm + -0.6408303732204232*cHW*ctWRe + -0.0077879672123913*cHWtil*ctHIm + 
            -0.0045140089094199*cHWtil*ctHRe + -0.4445489789604543*cHWtil*ctWIm + 0.4626096180508141*cHWtil*ctWRe +
            -0.005594974377675*ctHIm*ctHRe + 0.0768049243145821*ctHIm*ctWIm + -0.0251294260739384*ctHIm*ctWRe +
            0.0703429679105624*ctHRe*ctWIm + 0.1273317894395585*ctHRe*ctWRe + -0.0157222879603941*ctWIm*ctWRe;
            double cos_ye_pos_NP = 0.476717469811962*cHW + 0.0753663278229962*cHWtil + -0.0135420346155876*ctHIm + 
            -0.078196482933858*ctHRe + 0.599668867733271*ctWIm + -1.1296434970927198*ctWRe + 
            0.307351117172388*cHW*cHW + 0.2260176571753823*cHWtil*cHWtil + 0.0101846863583576*ctHIm*ctHIm +
            0.0313976407708736*ctHRe*ctHRe + 2.706365319141492*ctWIm*ctWIm + 2.703477707011085*ctWRe*ctWRe + 
            -0.0070410294141127*cHW*cHWtil + -0.0061980849317276*cHW*ctHIm + -0.0605650927340252*cHW*ctHRe +
            0.6487242343490647*cHW*ctWIm + -0.6374431506853131*cHW*ctWRe + -0.0080254293734157*cHWtil*ctHIm + 
            0.004591077351006*cHWtil*ctHRe + -0.4406240019648604*cHWtil*ctWIm + -0.4612663681943691*cHWtil*ctWRe + 
            0.0057555079260233*ctHIm*ctHRe + 0.0724481406951879*ctHIm*ctWIm + 0.0247574767566629*ctHIm*ctWRe + 
            -0.0724217725511453*ctHRe*ctWIm + 0.1217486012314901*ctHRe*ctWRe + 0.009125406250685*ctWIm*ctWRe;
    
    
            //double num = (cos_ye_pos_Mad + cos_ye_pos_NP) - (cos_ye_neg_Mad + cos_ye_neg_NP);
            //double deno = (cos_ye_pos_Mad + cos_ye_pos_NP) + (cos_ye_neg_Mad + cos_ye_neg_NP);
    
            
            //We set the SM difference to zero
            double num = (cos_ye_pos_NP) - (cos_ye_neg_NP);
            double deno = (cos_ye_pos_Mad + cos_ye_pos_NP) + (cos_ye_neg_Mad + cos_ye_neg_NP);
    
            double total = num/deno;
            return total;
            
    }
    else{

            
    
            double cos_ye_neg_NP = 0.4728811316749511*cHW + -0.0743571487036044*cHWtil + 0.0134093106580527*ctHIm +
            -0.0788765537664126*ctHRe + -0.5965260044476208*ctWIm + -1.133879975967292*ctWRe;
            double cos_ye_pos_NP = 0.476717469811962*cHW + 0.0753663278229962*cHWtil + -0.0135420346155876*ctHIm + 
            -0.078196482933858*ctHRe + 0.599668867733271*ctWIm + -1.1296434970927198*ctWRe;
    
    
            //double num = (cos_ye_pos_Mad + cos_ye_pos_NP) - (cos_ye_neg_Mad + cos_ye_neg_NP);
            //double deno = (cos_ye_pos_Mad + cos_ye_pos_NP) + (cos_ye_neg_Mad + cos_ye_neg_NP);
    
            
            //We set the SM difference to zero
            double num = (cos_ye_pos_NP) - (cos_ye_neg_NP);
            double deno = (cos_ye_pos_Mad + cos_ye_pos_NP) + (cos_ye_neg_Mad + cos_ye_neg_NP);
    
            double total = num/deno;
            return total;
            
    }  
}





Asymmetry_trip_prod_ph_pt_pj_tHj::Asymmetry_trip_prod_ph_pt_pj_tHj(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

//    setParametersForObservable(make_vector<std::string>() << "SM_Asymmetry_trip_prod_ph_pt_pj_tHj");

}

double Asymmetry_trip_prod_ph_pt_pj_tHj::computeThValue()
{
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    //double cHG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiG();
    //double cHGtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiGtil();
    double cHW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiW();
    double cHWtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiWtil();
    double ctHRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double ctHIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphiIm();
    //double ctGRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    //double ctGIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tGIm();
    double ctWRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double ctWIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tWIm();
    
    
    //Expressions are written in the basis of SMEFTsim
    //let's change to the basis of dim6top (same as
    //SMEFT@NLO except for the gs on CtG). These
    //models define the covariance derivative 
    //with a different convention
    
    if(flag_LHC_WG_Basis){
//        double ctZRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
//        double ctZIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZIm();
//        sw2 = 0.22305;
//        cw2 = 1 - sw2;
//        tw2 = sw2/cw2;
//        sw  = sqrt(sw2);
//        tw  = sqrt(tw2);
//        ctBRe = ctZRe/sw - ctWRe/tw;
//        ctBIm = ctZIm/sw - ctWIm/tw;
        ctWRe = -ctWRe;
        ctWIm = -ctWIm;
//        ctGRe = -ctGRe;
//        ctGIm = -ctGIm;
    }
//    else{
//        double ctBRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
//        double ctBIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tBIm();
//    }
    
    
    
    //double SM_Asymmetry_cos_ye_tHj = SM.getOptionalParameter("SM_Asymmetry_cos_ye_tHj");
    
    double trip_prod_neg_Mad = 2.997424819968232;
    double trip_prod_pos_Mad = 2.990894455466435;
    

    if(flag_Quadratic){

            
            double trip_prod_neg_NP = 0.4731751789769144*cHW + -0.2062173155393779*cHWtil + -0.0210420414673162*ctHIm + 
            -0.0784772415015998*ctHRe + 0.0035049087021471*ctWIm + -1.1126408682479452*ctWRe + 
            0.3062680149098*cHW*cHW + 0.2256041264624489*cHWtil*cHWtil + 0.0103033595638413*ctHIm*ctHIm + 
            0.0327718453878561*ctHRe*ctHRe + 2.706594862259218*ctWIm*ctWIm + 2.705927327451248*ctWRe*ctWRe +
            0.0580361245603444*cHW*cHWtil + -0.0028341721569468*cHW*ctHIm + -0.0603326748771091*cHW*ctHRe + 
            -0.0341451929677045*cHW*ctWIm + -0.6366438156186173*cHW*ctWRe + -0.0075429332474376*cHWtil*ctHIm + 
            -0.0255348605499368*cHWtil*ctHRe + -0.4407813603683847*cHWtil*ctWIm + 0.0116076094967981*cHWtil*ctWRe + 
            -0.0001820004036614*ctHIm*ctHRe + 0.0756006805298208*ctHIm*ctWIm + 0.0513525574146963*ctHIm*ctWRe + 
            -0.0407550906968653*ctHRe*ctWIm + 0.1280372222916361*ctHRe*ctWRe + 0.013498171229231*ctWIm*ctWRe;
            double trip_prod_pos_NP = 0.4764234225099984*cHW + 0.20722649465877*cHWtil + 0.0209093175097812*ctHIm + 
            -0.078595795198671*ctHRe + -0.000362045416562*ctWIm + -1.150882604826262*ctWRe +
            0.3077131779238463*cHW*cHW + 0.2270079953476132*cHWtil*cHWtil + 0.0104011206292801*ctHIm*ctHIm + 
            0.0334424612581808*ctHRe*ctHRe + 2.710148441129685*ctWIm*ctWIm + 2.710423562420865*ctWRe*ctWRe + 
            -0.0582782715491244*cHW*cHWtil + 0.0028875063300902*cHW*ctHIm + -0.0599189477103146*cHW*ctHRe + 
            0.0324241655935762*cHW*ctWIm + -0.6418354389210809*cHW*ctWRe + -0.0082755885390069*cHWtil*ctHIm +
            0.0256119540036406*cHWtil*ctHRe + -0.4445427933974889*cHWtil*ctWIm + -0.0102638339471822*cHWtil*ctWRe + 
            0.0003425431837288*ctHIm*ctHRe + 0.0736905682201566*ctHIm*ctWIm + -0.0517246077639719*ctHIm*ctWRe + 
            0.038675770013493*ctHRe*ctWIm + 0.1210837943540063*ctHRe*ctWRe + -0.0200958660055103*ctWIm*ctWRe;
    
    
            //double num = (trip_prod_pos_Mad + trip_prod_pos_NP) - (trip_prod_neg_Mad + trip_prod_neg_NP);
            //double deno = (trip_prod_pos_Mad + trip_prod_pos_NP) + (trip_prod_neg_Mad + trip_prod_neg_NP);
    
            
            //We set the SM difference to zero
            double num = (trip_prod_pos_NP) - (trip_prod_neg_NP);
            double deno = (trip_prod_pos_Mad + trip_prod_pos_NP) + (trip_prod_neg_Mad + trip_prod_neg_NP);
    
            double total = num/deno;
            return total;
            
        }
        else{

            double trip_prod_neg_NP = 0.4731751789769144*cHW + -0.2062173155393779*cHWtil + -0.0210420414673162*ctHIm + 
            -0.0784772415015998*ctHRe + 0.0035049087021471*ctWIm + -1.1126408682479452*ctWRe;
            double trip_prod_pos_NP = 0.4764234225099984*cHW + 0.20722649465877*cHWtil + 0.0209093175097812*ctHIm + 
            -0.078595795198671*ctHRe + -0.000362045416562*ctWIm + -1.150882604826262*ctWRe;
    
    
            //double num = (trip_prod_pos_Mad + trip_prod_pos_NP) - (trip_prod_neg_Mad + trip_prod_neg_NP);
            //double deno = (trip_prod_pos_Mad + trip_prod_pos_NP) + (trip_prod_neg_Mad + trip_prod_neg_NP);
    
            
            //We set the SM difference to zero
            double num = (trip_prod_pos_NP) - (trip_prod_neg_NP);
            double deno = (trip_prod_pos_Mad + trip_prod_pos_NP) + (trip_prod_neg_Mad + trip_prod_neg_NP);
    
            double total = num/deno;
            return total;
            
    }  
}








sigma_tHj_diff_LO_Del_R_th::sigma_tHj_diff_LO_Del_R_th(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tHj_diff_LO_Del_R_th_bin_0_pi" 
            << "SM_sigma_tHj_diff_LO_Del_R_th_bin_pi_8");
    
}

double sigma_tHj_diff_LO_Del_R_th::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    //double cHG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiG();
    //double cHGtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiGtil();
    double cHW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiW();
    double cHWtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiWtil();
    double ctHRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double ctHIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphiIm();
    //double ctGRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    //double ctGIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tGIm();
    double ctWRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double ctWIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tWIm();
    
    
        
    //Expressions are written in the basis of SMEFTsim
    //let's change to the basis of dim6top (same as
    //SMEFT@NLO except for the gs on CtG). These
    //models define the covariance derivative 
    //with a different convention
    
    if(flag_LHC_WG_Basis){
//        double ctZRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
//        double ctZIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZIm();
//        sw2 = 0.22305;
//        cw2 = 1 - sw2;
//        tw2 = sw2/cw2;
//        sw  = sqrt(sw2);
//        tw  = sqrt(tw2);
//        ctBRe = ctZRe/sw - ctWRe/tw;
//        ctBIm = ctZIm/sw - ctWIm/tw;
        ctWRe = -ctWRe;
        ctWIm = -ctWIm;
//        ctGRe = -ctGRe;
//        ctGIm = -ctGIm;
    }
//    else{
//        double ctBRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
//        double ctBIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tBIm();
//    }
    
    
   
    if(b_min == 0 && b_max == 3.14){
        
        double SM_Del_R_th_tHj_bin_0_pi = SM.getOptionalParameter("SM_sigma_tHj_diff_LO_Del_R_th_bin_0_pi");
        double Del_R_th_tHj_bin_0_pi_madgraph = 0.5198924690483601;
        double total;
        
        
        if(flag_Quadratic){
                total = SM_Del_R_th_tHj_bin_0_pi + (0.220246481678596*cHW + 0.0001485728632859*cHWtil + 0.0001008891760229*ctHIm +
                        -0.0600527834531858*ctHRe + 0.0004390818255602*ctWIm + -0.526439687822913*ctWRe + 
                        0.0878777287723505*cHW*cHW + 0.0537182291189895*cHWtil*cHWtil + 0.00257295232299*ctHIm*ctHIm + 
                        0.0100192193344844*ctHRe*ctHRe + 0.8714545891772361*ctWIm*ctWIm + 0.8683513697290252*ctWRe*ctWRe + 
                        1.1954223612303805e-05*cHW*cHWtil + -0.0001657857250467*cHW*ctHIm + -0.024684972339762*cHW*ctHRe + 
                        1.755473359284221e-05*cHW*ctWIm + -0.2590069019348239*cHW*ctWRe + -0.0030608640639377*cHWtil*ctHIm + 
                        0.0001448608174933*cHWtil*ctHRe + -0.1377822969273901*cHWtil*ctWIm + -0.0002095049296119*cHWtil*ctWRe + 
                        0.0001531599176172*ctHIm*ctHRe + 0.0233551969992832*ctHIm*ctWIm + 0.0001221412647315*ctHIm*ctWRe + 
                        -0.0004981550659563*ctHRe*ctWIm + 0.0464732517727849*ctHRe*ctWRe + -0.0017079266931867*ctWIm*ctWRe
                        )*(SM_Del_R_th_tHj_bin_0_pi/Del_R_th_tHj_bin_0_pi_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        else{
                total = SM_Del_R_th_tHj_bin_0_pi + (0.220246481678596*cHW + 
                        -0.0600527834531858*ctHRe + -0.526439687822913*ctWRe
                        )*(SM_Del_R_th_tHj_bin_0_pi/Del_R_th_tHj_bin_0_pi_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        
    } else if(b_min == 3.14 && b_max == 8.){
        
        double SM_Del_R_th_tHj_bin_pi_8 = SM.getOptionalParameter("SM_sigma_tHj_diff_LO_Del_R_th_bin_pi_8");
        double Del_R_th_tHj_bin_pi_8_madgraph = 0.9286855031587032;
        double total;
        
        
        if(flag_Quadratic){
                
                total = SM_Del_R_th_tHj_bin_pi_8 + (0.0897896143595042*cHW + 2.0137369438635863e-05*cHWtil + -0.0001184136133698*ctHIm +
                        0.0022046467428809*ctHRe + 0.000143163526555*ctWIm + -0.2092845561240396*ctWRe + 
                        0.0844443789341863*cHW*cHW + 0.0623462174258319*cHWtil*cHWtil + 0.0026189688979234*ctHIm*ctHIm + 
                        0.0072496081836699*ctHRe*ctHRe + 0.6634699666908053*ctWIm*ctWIm + 0.6652835432621799*ctWRe*ctWRe + 
                        8.08782734279534e-06*cHW*cHWtil + 9.311417809410028e-05*cHW*ctHIm + -0.0106075370870455*cHW*ctHRe + 
                        -0.0003069951264093*cHW*ctWIm + -0.1574825552102092*cHW*ctWRe + -0.0015264174728027*cHWtil*ctHIm + 
                        -8.004726796563944e-05*cHWtil*ctHRe + -0.0987965665599762*cHWtil*ctWIm + 0.0001784457437636*cHWtil*ctWRe + 
                        -4.8553041115131374e-05*ctHIm*ctHRe + 0.0183802012849704*ctHIm*ctWIm + -0.0001651313443605*ctHIm*ctWRe + 
                        9.2546382412656e-05*ctHRe*ctWIm + 0.0280427583650578*ctHRe*ctWRe + 0.0003215575263969*ctWIm*ctWRe
                        )*(SM_Del_R_th_tHj_bin_pi_8/Del_R_th_tHj_bin_pi_8_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        else{
                
                total = SM_Del_R_th_tHj_bin_pi_8 + (0.0897896143595042*cHW +
                        0.0022046467428809*ctHRe + -0.2092845561240396*ctWRe
                        )*(SM_Del_R_th_tHj_bin_pi_8/Del_R_th_tHj_bin_pi_8_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
        }
        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct binning for sigma_tHj_diff_LO_Del_R_th. \n");
    }    
    
}









sigma_tHj_diff_LO_mth::sigma_tHj_diff_LO_mth(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tHj_diff_LO_mth_bin_200_340" 
            << "SM_sigma_tHj_diff_LO_mth_bin_340_424" << "SM_sigma_tHj_diff_LO_mth_bin_424_620"
            << "SM_sigma_tHj_diff_LO_mth_bin_620_1600" );
    
}

double sigma_tHj_diff_LO_mth::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    //double cHG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiG();
    //double cHGtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiGtil();
    double cHW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiW();
    double cHWtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiWtil();
    double ctHRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double ctHIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphiIm();
    //double ctGRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    //double ctGIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tGIm();
    double ctWRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double ctWIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tWIm();
    
    
    
    //Expressions are written in the basis of SMEFTsim
    //let's change to the basis of dim6top (same as
    //SMEFT@NLO except for the gs on CtG). These
    //models define the covariance derivative 
    //with a different convention
    
    if(flag_LHC_WG_Basis){
//        double ctZRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
//        double ctZIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZIm();
//        sw2 = 0.22305;
//        cw2 = 1 - sw2;
//        tw2 = sw2/cw2;
//        sw  = sqrt(sw2);
//        tw  = sqrt(tw2);
//        ctBRe = ctZRe/sw - ctWRe/tw;
//        ctBIm = ctZIm/sw - ctWIm/tw;
        ctWRe = -ctWRe;
        ctWIm = -ctWIm;
//        ctGRe = -ctGRe;
//        ctGIm = -ctGIm;
    }
//    else{
//        double ctBRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
//        double ctBIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tBIm();
//    }
    
    
   
    if(b_min == 200 && b_max == 340){
        
        double SM_sigma_tHj_diff_LO_mth_bin_200_340 = SM.getOptionalParameter("SM_sigma_tHj_diff_LO_mth_bin_200_340");
        double sigma_tHj_diff_LO_mth_bin_200_340_madgraph = 4.8312389333236;
        double total;
        
        
        if(flag_Quadratic){
                
                total = SM_sigma_tHj_diff_LO_mth_bin_200_340 + (2.044625422571*cHW + -0.0019697321068240514*cHWtil + -0.0002770594241113811*ctHIm +
                        -0.9319723301799*ctHRe + -0.00476895114764897*ctWIm + -4.0956850184244*ctWRe + 
                        0.2928532476119*cHW*cHW + 0.09546992012006428*cHWtil*cHWtil + 0.01454612189076554*ctHIm*ctHIm + 
                        0.1209362969154*ctHRe*ctHRe + 2.7343139034392*ctWIm*ctWIm + 2.6995685969779997*ctWRe*ctWRe + 
                        -0.005934498577654642*cHW*cHWtil + -0.0018803268831102064*cHW*ctHIm + -0.2194966254055*cHW*ctHRe + 
                        0.003760700951024876*cHW*ctWIm + -1.2677935653004*cHW*ctWRe + -0.02264745989419263*cHWtil*ctHIm + 
                        -0.000901401797751721*cHWtil*ctHRe + -0.5172919603443*cHWtil*ctWIm + 0.0005135710164916318*cHWtil*ctWRe + 
                        0.001622458841596397*ctHIm*ctHRe + 0.1287365371493*ctHIm*ctWIm + -0.005794226903136551*ctHIm*ctWRe + 
                        -0.0036847586393026423*ctHRe*ctWIm + 0.37264722012039997*ctHRe*ctWRe + 0.00219928660555313*ctWIm*ctWRe
                        )*(SM_sigma_tHj_diff_LO_mth_bin_200_340/sigma_tHj_diff_LO_mth_bin_200_340_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        else{
                
                total = SM_sigma_tHj_diff_LO_mth_bin_200_340 + (2.044625422571*cHW +
                        -0.9319723301799*ctHRe + -4.0956850184244*ctWRe
                        )*(SM_sigma_tHj_diff_LO_mth_bin_200_340/sigma_tHj_diff_LO_mth_bin_200_340_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
        
        }
        
    } else if(b_min == 340 && b_max == 424){
        
        double SM_sigma_tHj_diff_LO_mth_bin_340_424 = SM.getOptionalParameter("SM_sigma_tHj_diff_LO_mth_bin_340_424");
        double sigma_tHj_diff_LO_mth_bin_340_424_madgraph = 15.7487145559461;
        double total;
        
        
        if(flag_Quadratic){
                
                total = SM_sigma_tHj_diff_LO_mth_bin_340_424 + (4.9981313432903995*cHW + -0.00433369039049003*cHWtil + 0.00039224905580559366*ctHIm +
                        -1.0434688823926*ctHRe + 0.003958581262702987*ctWIm + -11.327199630097*ctWRe + 
                        1.3719672914796999*cHW*cHW + 0.6464375218769001*cHWtil*cHWtil + 0.06497929284932091*ctHIm*ctHIm + 
                        0.2766726496071*ctHRe*ctHRe + 10.297289637753002*ctWIm*ctWIm + 10.2481785950402*ctWRe*ctWRe + 
                        -0.015556904787081122*cHW*cHWtil + 0.00012779846093879946*cHW*ctHIm + -0.5974518855705999*cHW*ctHRe + 
                        0.01289784306163222*cHW*ctWIm + -4.6693468708844*cHW*ctWRe + -0.07474505546859601*cHWtil*ctHIm + 
                        -0.0017800266970389987*cHWtil*ctHRe + -2.2400298269615*cHWtil*ctWIm + -0.018183520573615563*cHWtil*ctWRe + 
                        -0.0003375385250080826*ctHIm*ctHRe + 0.5066845255559*ctHIm*ctWIm + 0.02001141484961333*ctHIm*ctWRe +
                        0.02423769252751784*ctHRe*ctWIm + 1.0459736936309*ctHRe*ctWRe + 0.0640694541759329*ctWIm*ctWRe
                        )*(SM_sigma_tHj_diff_LO_mth_bin_340_424/sigma_tHj_diff_LO_mth_bin_340_424_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        else{
                
                total = SM_sigma_tHj_diff_LO_mth_bin_340_424 + (4.9981313432903995*cHW +
                        -1.0434688823926*ctHRe + -11.327199630097*ctWRe
                        )*(SM_sigma_tHj_diff_LO_mth_bin_340_424/sigma_tHj_diff_LO_mth_bin_340_424_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
                
        
        }
        
        
    }  else if(b_min == 424 && b_max == 620){
        
        double SM_sigma_tHj_diff_LO_mth_bin_424_620 = SM.getOptionalParameter("SM_sigma_tHj_diff_LO_mth_bin_424_620");
        double sigma_tHj_diff_LO_mth_bin_424_620_madgraph = 10.9357549691572;
        double total;
        
        
        if(flag_Quadratic){
                
                total = SM_sigma_tHj_diff_LO_mth_bin_424_620 + (1.5666660066113*cHW + 0.009873964296414073*cHWtil + -0.003017092917417808*ctHIm +
                        0.08256721872489103*ctHRe + -0.008516060692950744*ctWIm + -4.136280835257099*ctWRe + 
                        1.0890377895256*cHW*cHW + 0.7007065924062*cHWtil*cHWtil + 0.0390975212028466*ctHIm*ctHIm + 
                        0.09289750818528009*ctHRe*ctHRe + 7.2496001463796*ctWIm*ctWIm + 7.2836843498583*ctWRe*ctWRe +
                        0.009983943811839311*cHW*cHWtil + -0.0012669630490372552*cHW*ctHIm + -0.1910815830423*cHW*ctHRe +
                        -0.0110446887577412*cHW*ctWIm + -2.7163356332474002*cHW*ctWRe + -0.02914841005248382*cHWtil*ctHIm + 
                        -0.0011939154296200202*cHWtil*ctHRe + -1.6033076535583*cHWtil*ctWIm + -0.003923502799203852*cHWtil*ctWRe +
                        0.0007781238255122336*ctHIm*ctHRe + 0.2845742436244*ctHIm*ctWIm + -0.006589865014594622*ctHIm*ctWRe + 
                        -0.014562482561483663*ctHRe*ctWIm + 0.4522172735994*ctHRe*ctWRe + -0.03856803764934983*ctWIm*ctWRe
                        )*(SM_sigma_tHj_diff_LO_mth_bin_424_620/sigma_tHj_diff_LO_mth_bin_424_620_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        else{
                
                total = SM_sigma_tHj_diff_LO_mth_bin_424_620 + (1.5666660066113*cHW +
                        0.08256721872489103*ctHRe -4.136280835257099*ctWRe
                        )*(SM_sigma_tHj_diff_LO_mth_bin_424_620/sigma_tHj_diff_LO_mth_bin_424_620_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
                
        }
        
        
    }  else if(b_min == 620 && b_max == 1600){
        
        double SM_sigma_tHj_diff_LO_mth_bin_620_1600 = SM.getOptionalParameter("SM_sigma_tHj_diff_LO_mth_bin_620_1600");
        double sigma_tHj_diff_LO_mth_bin_620_1600_madgraph = 2.0062389761137;
        double total;
        
        
        if(flag_Quadratic){
                
                total = SM_sigma_tHj_diff_LO_mth_bin_620_1600 + (0.1111424096413*cHW + -0.0014860638211933177*cHWtil + 0.0003135970091938667*ctHIm +
                        0.02440742106196714*ctHRe + 0.0029584200534060745*ctWIm + -0.3390724228498*ctWRe + 
                        0.29202717282470003*cHW*cHW + 0.2430289940472*cHWtil*cHWtil + 0.005608877718665806*ctHIm*ctHIm + 
                        0.008411593347051882*ctHRe*ctHRe + 2.6193177426385*ctWIm*ctWIm + 2.6275190748015*ctWRe*ctWRe + 
                        0.0003166402726267892*cHW*cHWtil + 0.0003708859349638161*cHW*ctHIm + -0.010811319724541921*cHW*ctHRe + 
                        -0.004126259538358132*cHW*ctWIm + -0.4598013457796*cHW*ctWRe + -0.0019298379837962567*cHWtil*ctHIm +
                        0.0005878696857919419*cHWtil*ctHRe + -0.3198592839751*cHWtil*ctWIm + 0.0026820414812528397*cHWtil*ctWRe +
                        -8.366197674275488e-05*ctHIm*ctHRe + 0.045894456459680266*ctHIm*ctWIm + -0.0014537499925816463*ctHIm*ctWRe + 
                        -0.0007806602952328085*ctHRe*ctWIm + 0.05445691902987437*ctHRe*ctWRe + -0.0006867750649745359*ctWIm*ctWRe
                        )*(SM_sigma_tHj_diff_LO_mth_bin_620_1600/sigma_tHj_diff_LO_mth_bin_620_1600_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        else{
                
                total = SM_sigma_tHj_diff_LO_mth_bin_620_1600 + (0.1111424096413*cHW +
                        0.02440742106196714*ctHRe + -0.3390724228498*ctWRe
                        )*(SM_sigma_tHj_diff_LO_mth_bin_620_1600/sigma_tHj_diff_LO_mth_bin_620_1600_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
                
        }
        
        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct binning for sigma_tHj_diff_LO_mth. \n");
    }    
    
}








sigma_tHj_diff_LO_trip_prod_z_pt_pj::sigma_tHj_diff_LO_trip_prod_z_pt_pj(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{

    setParametersForObservable(make_vector<std::string>() << "SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m1_m0p1" 
            << "SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m0p1_0" << "SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0_0p1"
            << "SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0p1_1" );
    
}

double sigma_tHj_diff_LO_trip_prod_z_pt_pj::computeThValue()
{
    
    b_min = getBinMin();
    b_max = getBinMax();
    
    
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    //double cHG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiG();
    //double cHGtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiGtil();
    double cHW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiW();
    double cHWtil = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiWtil();
    double ctHRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double ctHIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphiIm();
    //double ctGRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    //double ctGIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tGIm();
    double ctWRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double ctWIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tWIm();
    
    
    //Expressions are written in the basis of SMEFTsim
    //let's change to the basis of dim6top (same as
    //SMEFT@NLO except for the gs on CtG). These
    //models define the covariance derivative 
    //with a different convention
    
    if(flag_LHC_WG_Basis){
//        double ctZRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
//        double ctZIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZIm();
//        sw2 = 0.22305;
//        cw2 = 1 - sw2;
//        tw2 = sw2/cw2;
//        sw  = sqrt(sw2);
//        tw  = sqrt(tw2);
//        ctBRe = ctZRe/sw - ctWRe/tw;
//        ctBIm = ctZIm/sw - ctWIm/tw;
        ctWRe = -ctWRe;
        ctWIm = -ctWIm;
//        ctGRe = -ctGRe;
//        ctGIm = -ctGIm;
    }
//    else{
//        double ctBRe = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
//        double ctBIm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tBIm();
//    }
    
   
    if(b_min == -1 && b_max == -0.1){
        
        double SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m1_m0p1 = SM.getOptionalParameter("SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m1_m0p1");
        double sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m1_m0p1_madgraph = 0.8264858170834004;
        double total;
        
        if(flag_Quadratic){
                
                total = SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m1_m0p1 + (0.3646362275869991*cHW + -0.0166895627470372*cHWtil + 
                        0.0256209852360194*ctHIm + -0.128543600482164*ctHRe + -0.0939877611285754*ctWIm + -0.7936049243097774*ctWRe +
                        0.1278696411053465*cHW*cHW + 0.0731620368677026*cHWtil*cHWtil + 0.0041888307920575*ctHIm*ctHIm + 
                        0.0193907320385548*ctHRe*ctHRe + 1.0675602775250457*ctWIm*ctWIm + 1.0687772237286093*ctWRe*ctWRe + 
                        -0.0782150127334121*cHW*cHWtil + 0.0025943196337468*cHW*ctHIm + -0.0496091674851035*cHW*ctHRe + 
                        0.0254989881300741*cHW*ctWIm + -0.334101544416724*cHW*ctWRe + -0.0039346191573913*cHWtil*ctHIm + 
                        0.0282265393923323*cHWtil*ctHRe + -0.2033655508825295*cHWtil*ctWIm + -0.0165254794118828*cHWtil*ctWRe + 
                        -2.621088637231366e-05*ctHIm*ctHRe + 0.0231618106640849*ctHIm*ctWIm + -0.0427921919641906*ctHIm*ctWRe + 
                        0.0307631048922972*ctHRe*ctWIm + 0.0632664748372589*ctHRe*ctWRe + -0.0083383017737166*ctWIm*ctWRe
                        )*(SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m1_m0p1/sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m1_m0p1_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        else{
                
                total = SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m1_m0p1 + (0.3646362275869991*cHW + -0.0166895627470372*cHWtil + 
                        0.0256209852360194*ctHIm + -0.128543600482164*ctHRe + -0.0939877611285754*ctWIm + -0.7936049243097774*ctWRe
                        )*(SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m1_m0p1/sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m1_m0p1_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
        
        }
        
    } else if(b_min == -0.1 && b_max == 0){
        
        double SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m0p1_0 = SM.getOptionalParameter("SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m0p1_0");
        double sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m0p1_0_madgraph = 22.45084375748576;
        double total;
        
        if(flag_Quadratic){
                
                total = SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m0p1_0 + (1.442341421722291*cHW + 4.019977713678601*cHWtil + 
                        0.153196279604299*ctHIm + 0.3625647084148817*ctHRe + 0.7530568161311931*ctWIm + -4.19557216302995*ctWRe +
                        1.908492498127476*cHW*cHW + 1.614134281181288*cHWtil*cHWtil + 0.0669586130040683*ctHIm*ctHIm + 
                        0.1603252933722704*ctHRe*ctHRe + 17.50734966952049*ctWIm*ctWIm + 17.514736467909568*ctWRe*ctWRe + 
                        -0.3525825156328727*cHW*cHWtil + 0.0197738212419989*cHW*ctHIm + -0.1592811472687367*cHW*ctHRe + 
                        0.2387483055297098*cHW*ctWIm + -3.380438409950413*cHW*ctWRe + -0.0451921162715584*cHWtil*ctHIm +
                        0.2065480016189984*cHWtil*ctHRe + -2.6391380244958267*cHWtil*ctWIm + -0.0237387770588975*cHWtil*ctWRe +
                        0.0011320829420689*ctHIm*ctHRe + 0.5343061646884646*ctHIm*ctWIm + -0.4926921646822787*ctHIm*ctWRe +
                        0.4285201539717223*ctHRe*ctWIm + 0.6738494834458951*ctHRe*ctWRe + -0.1121784940758641*ctWIm*ctWRe
                        )*(SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m0p1_0/sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m0p1_0_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

        }
        else{
                
                total = SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m0p1_0 + (1.442341421722291*cHW + 4.019977713678601*cHWtil + 
                        0.153196279604299*ctHIm + 0.3625647084148817*ctHRe + 0.7530568161311931*ctWIm + -4.19557216302995*ctWRe
                        )*(SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m0p1_0/sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_m0p1_0_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
                
        }
        
        
    }  else if(b_min == 0 && b_max == 0.1){
        
        double SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0_0p1 = SM.getOptionalParameter("SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0_0p1");
        double sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0_0p1_madgraph = 22.534252619573525;
        double total;
        
        if(flag_Quadratic){
                
                total = SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0_0p1 + (1.473504604282017*cHW + -3.9759593578800665*cHWtil +
                        -0.1574375509821199*ctHIm + 0.3831143236353689*ctHRe + -0.6988322279880105*ctWIm + -4.175337801382581*ctWRe +
                        1.920195946765777*cHW*cHW + 1.604637203954196*cHWtil*cHWtil + 0.0641534545419859*ctHIm*ctHIm + 
                        0.1565135244341224*ctHRe*ctHRe + 17.428569162642475*ctWIm*ctWIm + 17.46609058059385*ctWRe*ctWRe +
                        0.3551360188717259*cHW*cHWtil + -0.015651082292107*cHW*ctHIm + -0.146335080030712*cHW*ctHRe + 
                        -0.2747404294361077*cHW*ctWIm + -3.394518055348159*cHW*ctWRe + -0.0430494668268945*cHWtil*ctHIm + 
                        -0.2027440689736365*cHWtil*ctHRe + -2.603236269869225*cHWtil*ctWIm + 0.0682574309145335*cHWtil*ctWRe +
                        0.0009754633482168*ctHIm*ctHRe + 0.5455637449209495*ctHIm*ctWIm + 0.4770873023247127*ctHIm*ctWRe + 
                        -0.441382381279638*ctHRe*ctWIm + 0.7163879472097825*ctHRe*ctWRe + 0.0441737265426905*ctWIm*ctWRe
                        )*(SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0_0p1/sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0_0p1_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

            }
            else{
                
                total = SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0_0p1 + (1.473504604282017*cHW + -3.9759593578800665*cHWtil +
                        -0.1574375509821199*ctHIm + 0.3831143236353689*ctHRe + -0.6988322279880105*ctWIm + -4.175337801382581*ctWRe
                        )*(SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0_0p1/sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0_0p1_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
        }
        
    }  else if(b_min == 0.1 && b_max == 1){
        
        double SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0p1_1 = SM.getOptionalParameter("SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0p1_1");
        double sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0p1_1_madgraph = 0.8288582248370959;
        double total;
        
        if(flag_Quadratic){
                
                total = SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0p1_1 + (0.3664904378118876*cHW + 0.0129199443509313*cHWtil + 
                        -0.0252972039241013*ctHIm + -0.1288352216337421*ctHRe + 0.091454877214733*ctWIm + -0.7913200499790009*ctWRe +
                        0.1289218570571187*cHW*cHW + 0.0720990445685774*cHWtil*cHWtil + 0.0042481396953987*ctHIm*ctHIm + 
                        0.018976406699871*ctHRe*ctHRe + 1.069274633793857*ctWIm*ctWIm + 1.062631870741071*ctWRe*ctWRe +
                        0.0776623121682901*cHW*cHWtil + -0.0029931540061552*cHW*ctHIm + -0.0500193309817597*cHW*ctHRe + 
                        -0.0234113227498942*cHW*ctWIm + -0.3334286909771868*cHW*ctWRe + -0.0038311889777161*cHWtil*ctHIm +
                        -0.0285635670320629*cHWtil*ctHRe + -0.1976737308512463*cHWtil*ctWIm + 0.013071461654528*cHWtil*ctWRe + 
                        -2.959037628247512e-05*ctHIm*ctHRe + 0.022689293687742*ctHIm*ctWIm + 0.0441127888606079*ctHIm*ctWRe +
                        -0.0316437501216862*ctHRe*ctWIm + 0.0590186720983172*ctHRe*ctWRe + 0.0085645224919978*ctWIm*ctWRe
                        )*(SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0p1_1/sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0p1_1_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;

            }
            else{
                
                total = SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0p1_1 + (0.3664904378118876*cHW + 0.0129199443509313*cHWtil + 
                        -0.0252972039241013*ctHIm + -0.1288352216337421*ctHRe + 0.091454877214733*ctWIm + -0.7913200499790009*ctWRe
                        )*(SM_sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0p1_1/sigma_tHj_diff_LO_trip_prod_z_pt_pj_bin_0p1_1_madgraph);
                //if (total < 0) return std::numeric_limits<double>::quiet_NaN();
                return total;
        }
        
        
    } else{
        throw std::runtime_error("\nERROR: Please specify a correct binning for sigma_tHj_diff_LO_mthj. \n");
    }    
    
}





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////
////////////////
//////////////// FROM NOW ON => OLD, NOT USED ANYMORE. CHECK AND DELETE!!!
////////////////
////////////////



/////// Prospects of Linear Colders at 250 GeV //////////////////////////////////////////////////////////


sigma_250_bb_eLpR::sigma_250_bb_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double sigma_250_bb_eLpR::computeThValue()
{
    // double sigma_250_bb_eLpR_madgraph = 3290;//fb
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  1000*((0.31*(C_phiQ3/0.998)+0.31*((C_phiQm+C_phiQ3)/0.998)+0.05*(C_phib/(0.998))+0.27*(C_bW)*(-1)*(C_bW)*(-1)+0.077*((0.881533*C_bW-C_bZ)/0.47213)*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)
                        -0.25*(C_bW)*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+0.091*C_ed-0.064*C_eq+0.71*C_ld+3.77*C_lqP));
            }
            else{
                return  1000*(0.31*(C_phiQ3/0.998)+0.31*((C_phiQm+C_phiQ3)/0.998)+0.05*(C_phib/(0.998))
                        +0.091*C_ed-0.064*C_eq+0.71*C_ld+3.77*C_lqP);
            }
    }
    
    else{
            if(flag_Quadratic){
                return  1000*((0.31*(C_phiQ3/0.998)+0.31*(((C_phiQ1-C_phiQ3)+C_phiQ3)/0.998)+0.05*(C_phib/(0.998))+0.27*(C_bW)*(-1)*(C_bW)*(-1)+0.077*C_bB*(-1)*C_bB*(-1)
                        -0.25*(C_bW)*(-1)*C_bB*(-1)+0.091*C_ed-0.064*C_eq+0.71*C_ld+3.77*C_lqP));
            }
            else{
                return  1000*(0.31*(C_phiQ3/0.998)+0.31*(((C_phiQ1-C_phiQ3)+C_phiQ3)/0.998)+0.05*(C_phib/(0.998))
                        +0.091*C_ed-0.064*C_eq+0.71*C_ld+3.77*C_lqP);
            }
    }
}



a_250_bb_eLpR::a_250_bb_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double a_250_bb_eLpR::computeThValue()
{
    //double a_250_bb_eLpR_madgraph = 69.6; in percentage
    //double a_250_bb_eLpR_madgraph = 0.696; 
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_LHC_WG_Basis){
        if(flag_Quadratic){
        return ((0.4*(C_phiQ3/0.998)+0.3*(((C_phiQm+C_phiQ3))/0.998)-2.2*(C_phib/(0.998))-5.1*C_bW*(-1)*C_bW*(-1)-1.29*((0.881533*C_bW-C_bZ)/0.47213)*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)
                +4.46*C_bW*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+3.8*C_eq-29.5*C_ld+8.57*C_lqP))/100;
    }
    else{
        return (0.4*(C_phiQ3/0.998)+0.3*(((C_phiQm+C_phiQ3))/0.998)-2.2*(C_phib/(0.998))+3.8*C_eq-29.5*C_ld+8.57*C_lqP)/100;
    }
    }
    else{
        if(flag_Quadratic){
        return ((0.4*(C_phiQ3/0.998)+0.3*((C_phiQ1)/0.998)-2.2*(C_phib/(0.998))-5.1*C_bW*(-1)*C_bW*(-1)-1.29*C_bB*(-1)*C_bB*(-1)
                +4.46*C_bW*(-1)*C_bB*(-1)+3.8*C_eq-29.5*C_ld+8.57*C_lqP))/100;
    }
    else{
        return (0.4*(C_phiQ3/0.998)+0.3*((C_phiQ1)/0.998)-2.2*(C_phib/(0.998))+3.8*C_eq-29.5*C_ld+8.57*C_lqP)/100;
    }
    }
}


sigma_250_bb_eRpL::sigma_250_bb_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double sigma_250_bb_eRpL::computeThValue()
{
    //double sigma_250_bb_eRpL_madgraph = 1020;//fb
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_LHC_WG_Basis){
        if(flag_Quadratic){
        return 1000*(0.094*(C_phiQ3/0.998)+0.094*(((C_phiQm+C_phiQ3))/0.998)-0.11*(C_phib/(0.998))+0.018*C_bW*(-1)*C_bW*(-1)
                +0.31*((0.881533*C_bW-C_bZ)/0.47213)*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+0.023*C_bW*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+1.61*C_ed-1.08*C_eq+0.04*C_ld+0.23*C_lqP);
        }
        else{
            return 1000*(0.094*(C_phiQ3/0.998)+0.094*(((C_phiQm+C_phiQ3))/0.998)-0.11*(C_phib/(0.998))+1.61*C_ed-1.08*C_eq+0.04*C_ld+0.23*C_lqP);
        }
    }
    else{
        if(flag_Quadratic){
        return 1000*(0.094*(C_phiQ3/0.998)+0.094*((C_phiQ1)/0.998)-0.11*(C_phib/(0.998))+0.018*C_bW*(-1)*C_bW*(-1)
                +0.31*C_bB*(-1)*C_bB*(-1)+0.023*C_bW*(-1)*C_bB*(-1)+1.61*C_ed-1.08*C_eq+0.04*C_ld+0.23*C_lqP);
        }
        else{
            return 1000*(0.094*(C_phiQ3/0.998)+0.094*((C_phiQ1)/0.998)-0.11*(C_phib/(0.998))+1.61*C_ed-1.08*C_eq+0.04*C_ld+0.23*C_lqP);
        }
    }
}



a_250_bb_eRpL::a_250_bb_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double a_250_bb_eRpL::computeThValue()
{
    //double a_250_bb_eRpL_madgraph = 35.9; In percentage
    //double a_250_bb_eRpL_madgraph = 0.359; over one
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_LHC_WG_Basis){
        if(flag_Quadratic){
            return ((-7.8*(C_phiQ3/0.998)-7.7*(((C_phiQm+C_phiQ3))/0.998)-4.5*(C_phib/(0.998))-0.23*C_bW*(-1)*C_bW*(-1)-7.75*((0.881533*C_bW-C_bZ)/0.47213)*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)
                -0.61*C_bW*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+62*C_ed+119*C_eq-4.6*C_ld+7.9*C_lqP))/100;
        }
        else{
            return ((-7.8*(C_phiQ3/0.998)-7.7*(((C_phiQm+C_phiQ3))/0.998)-4.5*(C_phib/(0.998))+62*C_ed+119*C_eq-4.6*C_ld+7.9*C_lqP))/100;
        }
    }
    else{
        if(flag_Quadratic){
            return ((-7.8*(C_phiQ3/0.998)-7.7*((C_phiQ1)/0.998)-4.5*(C_phib/(0.998))-0.23*C_bW*(-1)*C_bW*(-1)-7.75*C_bB*(-1)*C_bB*(-1)
                -0.61*C_bW*(-1)*C_bB*(-1)+62*C_ed+119*C_eq-4.6*C_ld+7.9*C_lqP))/100;
        }
        else{
            return ((-7.8*(C_phiQ3/0.998)-7.7*((C_phiQ1)/0.998)-4.5*(C_phib/(0.998))+62*C_ed+119*C_eq-4.6*C_ld+7.9*C_lqP))/100;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////



/////// Prospects of Linear Colders at 500 GeV //////////////////////////////////////////////////////////

/////// 500 GeV bb observables //////////////////////////////////////////////////////////////////////////



sigma_500_bb_eLpR::sigma_500_bb_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double sigma_500_bb_eLpR::computeThValue()
{
   // double sigma_500_bb_eLpR_madgraph = 720;//fb
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_LHC_WG_Basis){
        if(flag_Quadratic){
        return 1000*(0.064*(C_phiQ3/0.998)+0.064*(((C_phiQm+C_phiQ3))/0.998)+0.012*(C_phib/(0.998))+0.24*C_bW*(-1)*C_bW*(-1)+0.085*((0.881533*C_bW-C_bZ)/0.47213)*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)
                -0.255*C_bW*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+0.095*C_ed-0.05*C_eq+0.76*C_ld+3.5*C_lqP);
        }
        else{
            return 1000*(0.064*(C_phiQ3/0.998)+0.064*(((C_phiQm+C_phiQ3))/0.998)+0.012*(C_phib/(0.998))+0.095*C_ed-0.05*C_eq+0.76*C_ld+3.5*C_lqP);
        }
    }
    else{
        if(flag_Quadratic){
        return 1000*(0.064*(C_phiQ3/0.998)+0.064*((C_phiQ1)/0.998)+0.012*(C_phib/(0.998))+0.24*C_bW*(-1)*C_bW*(-1)+0.085*C_bB*(-1)*C_bB*(-1)
                -0.255*C_bW*(-1)*C_bB*(-1)+0.095*C_ed-0.05*C_eq+0.76*C_ld+3.5*C_lqP);
        }
        else{
            return 1000*(0.064*(C_phiQ3/0.998)+0.064*((C_phiQ1)/0.998)+0.012*(C_phib/(0.998))+0.095*C_ed-0.05*C_eq+0.76*C_ld+3.5*C_lqP);
        }
    }
}


a_500_bb_eLpR::a_500_bb_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double a_500_bb_eLpR::computeThValue()
{
    //double a_500_bb_eLpR_madgraph = 67.7; In percentage
    //double a_500_bb_eLpR_madgraph = 0.677; over one
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_LHC_WG_Basis){
        if(flag_Quadratic){
            return (0.2*(C_phiQ3/0.998)+0.2*(((C_phiQm+C_phiQ3))/0.998)-2.3*(C_phib/(0.998))-13.7*C_bW*(-1)*C_bW*(-1)
                -3.8*((0.881533*C_bW-C_bZ)/0.47213)*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+13.5*C_bW*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+8*C_eq-139*C_ld+40.3*C_lqP)/100;
        }
        else{
            return (0.2*(C_phiQ3/0.998)+0.2*(((C_phiQm+C_phiQ3))/0.998)-2.3*(C_phib/(0.998))+8*C_eq-139*C_ld+40.3*C_lqP)/100;
        }
    }
    else{
        if(flag_Quadratic){
            return (0.2*(C_phiQ3/0.998)+0.2*((C_phiQ1)/0.998)-2.3*(C_phib/(0.998))-13.7*C_bW*(-1)*C_bW*(-1)
                -3.8*C_bB*(-1)*C_bB*(-1)+13.5*C_bW*(-1)*C_bB*(-1)+8*C_eq-139*C_ld+40.3*C_lqP)/100;
        }
        else{
            return (0.2*(C_phiQ3/0.998)+0.2*((C_phiQ1)/0.998)-2.3*(C_phib/(0.998))+8*C_eq-139*C_ld+40.3*C_lqP)/100;
        }
    }
}


sigma_500_bb_eRpL::sigma_500_bb_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double sigma_500_bb_eRpL::computeThValue()
{
   // double sigma_500_bb_eRpL_madgraph = 220;//fb
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_LHC_WG_Basis){
        if(flag_Quadratic){
            return 1000*(0.02*(C_phiQ3/0.998)+0.02*(((C_phiQm+C_phiQ3))/0.998)-0.024*(C_phib/(0.998))+0.014*C_bW*(-1)*C_bW*(-1)+0.29*((0.881533*C_bW-C_bZ)/0.47213)*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)
                -0.007*C_bW*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+1.56*C_ed-0.84*C_eq+0.046*C_ld+0.2*C_lqP);
        }
        else{
            return 1000*(0.02*(C_phiQ3/0.998)+0.02*(((C_phiQm+C_phiQ3))/0.998)-0.024*(C_phib/(0.998))+1.56*C_ed-0.84*C_eq+0.046*C_ld+0.2*C_lqP);
        }
    }
    else{
        if(flag_Quadratic){
            return 1000*(0.02*(C_phiQ3/0.998)+0.02*((C_phiQ1)/0.998)-0.024*(C_phib/(0.998))+0.014*C_bW*(-1)*C_bW*(-1)+0.29*C_bB*(-1)*C_bB*(-1)
                -0.007*C_bW*(-1)*C_bB*(-1)+1.56*C_ed-0.84*C_eq+0.046*C_ld+0.2*C_lqP);
        }
        else{
            return 1000*(0.02*(C_phiQ3/0.998)+0.02*((C_phiQ1)/0.998)-0.024*(C_phib/(0.998))+1.56*C_ed-0.84*C_eq+0.046*C_ld+0.2*C_lqP);
        }
    }
}



a_500_bb_eRpL::a_500_bb_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double a_500_bb_eRpL::computeThValue()
{
    //double a_500_bb_eRpL_madgraph = 46.7; in percentage
    //double a_500_bb_eRpL_madgraph = 0.467; over one
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_LHC_WG_Basis){
        if(flag_Quadratic){
            return (-6.9*(C_phiQ3/0.998)-7.1*(((C_phiQm+C_phiQ3))/0.998)-3.5*(C_phib/(0.998))-1.57*C_bW*(-1)*C_bW*(-1)-23*((0.881533*C_bW-C_bZ)/0.47213)*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)
                +0.39*C_bW*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+219*C_ed+380*C_eq-26.6*C_ld+28.4*C_lqP)/100;
        }
        else{
            return (-6.9*(C_phiQ3/0.998)-7.1*(((C_phiQm+C_phiQ3))/0.998)-3.5*(C_phib/(0.998))+219*C_ed+380*C_eq-26.6*C_ld+28.4*C_lqP)/100;
        }
    }
    else{
        if(flag_Quadratic){
            return (-6.9*(C_phiQ3/0.998)-7.1*((C_phiQ1)/0.998)-3.5*(C_phib/(0.998))-1.57*C_bW*(-1)*C_bW*(-1)-23*C_bB*(-1)*C_bB*(-1)
                +0.39*C_bW*(-1)*C_bB*(-1)+219*C_ed+380*C_eq-26.6*C_ld+28.4*C_lqP)/100;
        }
        else{
            return (-6.9*(C_phiQ3/0.998)-7.1*((C_phiQ1)/0.998)-3.5*(C_phib/(0.998))+219*C_ed+380*C_eq-26.6*C_ld+28.4*C_lqP)/100;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////// 500 GeV tt observables //////////////////////////////////////////////////////////////////////////
//These observables should be used when the optimal observables are used (they are redundant)


sigma_500_tt_eLpR::sigma_500_tt_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double sigma_500_tt_eLpR::computeThValue()
{
    //double sigma_500_tt_eLpR_madgraph = 930754;//fb
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    if(flag_LHC_WG_Basis){
        return 1000*(-35.373*(C_phit/0.998)+54.9562*(C_phiQ3/0.998)-54.9562*(((C_phiQm+C_phiQ3))/0.998)+639.01*(C_tW/(0.999*0.6532))+204.081*(((0.881533*C_tW-C_tZ)/0.47213))/(0.998*0.3492));
    }
    else{
        return 1000*(-35.373*(C_phit/0.998)+54.9562*(C_phiQ3/0.998)-54.9562*((C_phiQ1)/0.998)+639.01*(C_tW/(0.999*0.6532))+204.081*(C_tB)/(0.998*0.3492));
    }
}



a_500_tt_eLpR::a_500_tt_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double a_500_tt_eLpR::computeThValue()
{
    //double a_500_tt_eLpR_madgraph = -0.380287;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    if(flag_LHC_WG_Basis){
        return (-0.0127129*(C_phit/0.998)-0.0326318*(C_phiQ3/0.998)+0.0326318*(((C_phiQm+C_phiQ3))/0.998)-0.118262*(C_tW/(0.999*0.6532))-0.0384058*((((0.881533*C_tW-C_tZ)/0.47213) )/(0.998*0.3492)));
    }
    else{
        return (-0.0127129*(C_phit/0.998)-0.0326318*(C_phiQ3/0.998)+0.0326318*((C_phiQ1)/0.998)-0.118262*(C_tW/(0.999*0.6532))-0.0384058*((C_tB )/(0.998*0.3492)));
    }
}


sigma_500_tt_eRpL::sigma_500_tt_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double sigma_500_tt_eRpL::computeThValue()
{
   // double sigma_500_tt_eRpL_madgraph = 481105;//fb
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    if(flag_LHC_WG_Basis){
        return 1000*(32.3079*(C_phit/0.998)-16.7611*(C_phiQ3/0.998)+16.7611*(((C_phiQm+C_phiQ3))/0.998)+31.5853*(C_tW/(0.999*0.6532))+267.885*((((0.881533*C_tW-C_tZ)/0.47213) )/(0.998*0.3492)));
    }
    else{
        return 1000*(32.3079*(C_phit/0.998)-16.7611*(C_phiQ3/0.998)+16.7611*((C_phiQ1)/0.998)+31.5853*(C_tW/(0.999*0.6532))+267.885*((C_tB )/(0.998*0.3492)));
    }
}




a_500_tt_eRpL::a_500_tt_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double a_500_tt_eRpL::computeThValue()
{
   // double a_500_tt_eRpL_madgraph = -0.457824;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    if(flag_LHC_WG_Basis){
        return (-0.0413699*(C_phit/0.998)-0.0125082*(C_phiQ3/0.998)+0.0125082*(((C_phiQm+C_phiQ3))/0.998)-0.0109327*(C_tW/(0.999*0.6532))-0.125174*((((0.881533*C_tW-C_tZ)/0.47213) )/(0.998*0.3492)));
    }
    else{
        return (-0.0413699*(C_phit/0.998)-0.0125082*(C_phiQ3/0.998)+0.0125082*((C_phiQ1)/0.998)-0.0109327*(C_tW/(0.999*0.6532))-0.125174*((C_tB )/(0.998*0.3492)));
    }
}





pt_500_tt_eLpR::pt_500_tt_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double pt_500_tt_eLpR::computeThValue()
{
    //double pt_500_tt_eLpR_madgraph = (0.570477+0.573802+0.576393+0.570551+0.584227)/5;
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    if(flag_LHC_WG_Basis){
        return (-0.0196093*(C_phit/0.998)+0.215508*(C_tW/(0.999*0.6532))+0.0336945*((((0.881533*C_tW-C_tZ)/0.47213) )/(0.998*0.3492)));
    }
    else{
        return (-0.0196093*(C_phit/0.998)+0.215508*(C_tW/(0.999*0.6532))+0.0336945*((C_tB )/(0.998*0.3492)));
    }
}


pt_500_tt_eRpL::pt_500_tt_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double pt_500_tt_eRpL::computeThValue()
{
   // double pt_500_tt_eRpL_madgraph = (-0.432304-0.428132-0.423239-0.430406-0.431086)/5;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    if(flag_LHC_WG_Basis){
        return (-0.00550366*(C_phit/0.998)+0.0176743*(C_phiQ3/0.998)-0.0302046*(((C_phiQm+C_phiQ3))/0.998)+0.104522*(C_tW/(0.999*0.6532))-0.204084*((((0.881533*C_tW-C_tZ)/0.47213) )/(0.998*0.3492)));
    }
    else{
        return (-0.00550366*(C_phit/0.998)+0.0176743*(C_phiQ3/0.998)-0.0302046*((C_phiQ1)/0.998)+0.104522*(C_tW/(0.999*0.6532))-0.204084*((C_tB )/(0.998*0.3492)));
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////// 500 GeV Optimal observables //////////////////////////////////////////////////////////////////////////
//[arxiv:1807.02121]

op1::op1(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op1::computeThValue()
{
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
        return C_phit;
}


op2::op2(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op2::computeThValue()
{
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    if(flag_LHC_WG_Basis)
    {
        double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        return (C_phiQm);
    }
    else
    {
        double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
        double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
        
        return (C_phiQ1-C_phiQ3);
    }
}

op3::op3(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op3::computeThValue()
{
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
        return C_tW;
}

op4::op4(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op4::computeThValue()
{
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    
    if(flag_LHC_WG_Basis)
    {
        double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
        double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
        return (0.881533*C_tW-C_tZ)/0.47213;
    }
    else
    {
        double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        
        return C_tB;
    }
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////



/////// Prospects of Linear Colders at 1000 GeV /////////////////////////////////////////////////////////


/////// 1000 GeV bb observables /////////////////////////////////////////////////////////////////////////

sigma_1000_bb_eLpR::sigma_1000_bb_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double sigma_1000_bb_eLpR::computeThValue()
{
   // double sigma_1000_bb_eLpR_madgraph = 720;//fb
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_LHC_WG_Basis){
        if(flag_Quadratic){
            return 1000*(0.015*(C_phiQ3/0.998)+0.015*(((C_phiQm+C_phiQ3))/0.998)+0.003*(C_phib/(0.998))+0.23*C_bW*(-1)*C_bW*(-1)+0.087*((0.881533*C_bW-C_bZ)/0.47213)*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)
                -0.255*C_bW*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+0.093*C_ed-0.048*C_eq+0.77*C_ld+3.44*C_lqP);
        }
        else{
            return 1000*(0.015*(C_phiQ3/0.998)+0.015*(((C_phiQm+C_phiQ3))/0.998)+0.003*(C_phib/(0.998))+0.093*C_ed-0.048*C_eq+0.77*C_ld+3.44*C_lqP);
        }
    }
    else{
        if(flag_Quadratic){
            return 1000*(0.015*(C_phiQ3/0.998)+0.015*((C_phiQ1)/0.998)+0.003*(C_phib/(0.998))+0.23*C_bW*(-1)*C_bW*(-1)+0.087*C_bB*(-1)*C_bB*(-1)
                -0.255*C_bW*(-1)*C_bB*(-1)+0.093*C_ed-0.048*C_eq+0.77*C_ld+3.44*C_lqP);
        }
        else{
            return 1000*(0.015*(C_phiQ3/0.998)+0.015*((C_phiQ1)/0.998)+0.003*(C_phib/(0.998))+0.093*C_ed-0.048*C_eq+0.77*C_ld+3.44*C_lqP);
        }
    }
}


a_1000_bb_eLpR::a_1000_bb_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double a_1000_bb_eLpR::computeThValue()
{
    //double a_1000_bb_eLpR_madgraph = 67.7; In percentage
    //double a_1000_bb_eLpR_madgraph = 0.677; In percentage
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_LHC_WG_Basis){
        if(flag_Quadratic){
            return (0.48*(C_phiQ3/0.998)+0.48*(((C_phiQm+C_phiQ3))/0.998)-2.7*(C_phib/(0.998))-23.6*C_bW*(-1)*C_bW*(-1)-5.76*((0.881533*C_bW-C_bZ)/0.47213)*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)
                +21.2*C_bW*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+37.4*C_eq-593*C_ld+145*C_lqP)/100;
        }
        else{
            return (0.48*(C_phiQ3/0.998)+0.48*(((C_phiQm+C_phiQ3))/0.998)-2.7*(C_phib/(0.998))+37.4*C_eq-593*C_ld+145*C_lqP)/100;
        }
    }
    else{
        if(flag_Quadratic){
            return (0.48*(C_phiQ3/0.998)+0.48*((C_phiQ1)/0.998)-2.7*(C_phib/(0.998))-23.6*C_bW*(-1)*C_bW*(-1)-5.76*C_bB*(-1)*C_bB*(-1)
                +21.2*C_bW*(-1)*C_bB*(-1)+37.4*C_eq-593*C_ld+145*C_lqP)/100;
        }
        else{
            return (0.48*(C_phiQ3/0.998)+0.48*((C_phiQ1)/0.998)-2.7*(C_phib/(0.998))+37.4*C_eq-593*C_ld+145*C_lqP)/100;
        }
    }
}


sigma_1000_bb_eRpL::sigma_1000_bb_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double sigma_1000_bb_eRpL::computeThValue()
{
   // double sigma_1000_bb_eRpL_madgraph = 220;//fb
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_LHC_WG_Basis){
        if(flag_Quadratic){
            return 1000*(0.004*(C_phiQ3/0.998)+0.004*(((C_phiQm+C_phiQ3))/0.998)-0.006*(C_phib/(0.998))+0.014*C_bW*(-1)*C_bW*(-1)+0.29*((0.881533*C_bW-C_bZ)/0.47213)*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)
                -0.013*C_bW*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+1.55*C_ed-0.79*C_eq+0.046*C_ld+0.21*C_lqP);
        }
        else{
            return 1000*(0.004*(C_phiQ3/0.998)+0.004*(((C_phiQm+C_phiQ3))/0.998)-0.006*(C_phib/(0.998))+1.55*C_ed-0.79*C_eq+0.046*C_ld+0.21*C_lqP);
        }
    }
    else{
        if(flag_Quadratic){
            return 1000*(0.004*(C_phiQ3/0.998)+0.004*((C_phiQ1)/0.998)-0.006*(C_phib/(0.998))+0.014*C_bW*(-1)*C_bW*(-1)+0.29*C_bB*(-1)*C_bB*(-1)
                -0.013*C_bW*(-1)*C_bB*(-1)+1.55*C_ed-0.79*C_eq+0.046*C_ld+0.21*C_lqP);
        }
        else{
            return 1000*(0.004*(C_phiQ3/0.998)+0.004*((C_phiQ1)/0.998)-0.006*(C_phib/(0.998))+1.55*C_ed-0.79*C_eq+0.046*C_ld+0.21*C_lqP);
        }
    }
}



a_1000_bb_eRpL::a_1000_bb_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double a_1000_bb_eRpL::computeThValue()
{
    //double a_1000_bb_eRpL_madgraph = 46.7; In percentage
    //double a_1000_bb_eRpL_madgraph = 0.467;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_LHC_WG_Basis){
        if(flag_Quadratic){
            return (-6.8*(C_phiQ3/0.998)-6.8*(((C_phiQm+C_phiQ3))/0.998)-3.6*(C_phib/(0.998))-3.15*C_bW*(-1)*C_bW*(-1)-27.5*((0.881533*C_bW-C_bZ)/0.47213)*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)
                +0.44*C_bW*(-1)*((0.881533*C_bW-C_bZ)/0.47213)*(-1)+826*C_ed+1611*C_eq-120*C_ld+6.5*C_lqP)/100;
        }
        else{
            return (-6.8*(C_phiQ3/0.998)-6.8*(((C_phiQm+C_phiQ3))/0.998)-3.6*(C_phib/(0.998))+826*C_ed+1611*C_eq-120*C_ld+6.5*C_lqP)/100;
        }
    }
    else{
        if(flag_Quadratic){
            return (-6.8*(C_phiQ3/0.998)-6.8*((C_phiQ1)/0.998)-3.6*(C_phib/(0.998))-3.15*C_bW*(-1)*C_bW*(-1)-27.5*C_bB*(-1)*C_bB*(-1)
                +0.44*C_bW*(-1)*C_bB*(-1)+826*C_ed+1611*C_eq-120*C_ld+6.5*C_lqP)/100;
        }
        else{
            return (-6.8*(C_phiQ3/0.998)-6.8*((C_phiQ1)/0.998)-3.6*(C_phib/(0.998))+826*C_ed+1611*C_eq-120*C_ld+6.5*C_lqP)/100;
        }
    }
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////// 1000 GeV Optimal observables /////////////////////////////////////////////////////////////////////////
//[arxiv:1807.02121]




op_1000_1::op_1000_1(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op_1000_1::computeThValue()
{
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    if(flag_LHC_WG_Basis)
    {
        double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
        return (C_phiQm);
    }
    else
    {
        double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
        double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
        
        return (C_phiQ1-C_phiQ3);
    }
}



op_1000_2::op_1000_2(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op_1000_2::computeThValue()
{
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
        return C_phit;
}







op_1000_3::op_1000_3(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op_1000_3::computeThValue()
{
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
        return C_tW;
}


op_1000_4::op_1000_4(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op_1000_4::computeThValue()
{
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    
    if(flag_LHC_WG_Basis)
    {
        double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
        
        return C_tZ;
    }
    else
    {
        double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
        
        return -0.47213*C_tB+0.881533*C_tW;
    }
}



op_1000_5::op_1000_5(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op_1000_5::computeThValue()
{
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
        return C_lqM;
}



op_1000_6::op_1000_6(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op_1000_6::computeThValue()
{
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
        return C_eq;
}




op_1000_7::op_1000_7(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op_1000_7::computeThValue()
{
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
        return C_lu;
}



op_1000_8::op_1000_8(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op_1000_8::computeThValue()
{
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
        return C_eu;
}



op_eigen_ttll_1::op_eigen_ttll_1(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op_eigen_ttll_1::computeThValue()
{
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    
    return 0.815497*C_lqM -0.540171*C_eu -0.201143*C_lu -0.0521625*C_eq;
}



op_eigen_ttll_2::op_eigen_ttll_2(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op_eigen_ttll_2::computeThValue()
{
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    
    return 0.168067*C_lqM - 0.0696008*C_eu + 0.68548*C_lu + 0.705001*C_eq;
}




op_eigen_ttll_3::op_eigen_ttll_3(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op_eigen_ttll_3::computeThValue()
{
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    
    return 0.548904*C_lqM + 0.83473*C_eu - 0.00492088*C_lu - 0.0436616*C_eq;
}




op_eigen_ttll_4::op_eigen_ttll_4(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double op_eigen_ttll_4::computeThValue()
{
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    
    return 0.0736368*C_lqM - 0.0812221*C_eu + 0.699739*C_lu - 0.705936*C_eq;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// OBSERVABLES FOR PROSPECTS OF FUTURE COLLIDERS ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Translation to other basis //////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


gLt::gLt(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double gLt::computeThValue()
{
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double vev = 246.22;
   
    return -C_phiQm*(vev*vev/2.0);
}

gLb::gLb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double gLb::computeThValue()
{
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double vev = 246.22;
   
    return -(C_phiQm + 2.0*C_phiQ3)*(vev*vev/2.0);
}

gRt::gRt(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double gRt::computeThValue()
{
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double vev = 246.22;
   
    return -C_phit*(vev*vev/2.0);
}

gRb::gRb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{}

double gRb::computeThValue()
{
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double vev = 246.22;
   
    return -C_phib*(vev*vev/2.0);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


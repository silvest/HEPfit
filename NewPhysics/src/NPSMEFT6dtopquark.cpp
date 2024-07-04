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
    
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
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
        }
        
    } else throw std::runtime_error("wrong bin choice in sigma_ttH_diff_LO_ATLAS_231204450");
   
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
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
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
            
            return  SM_charge_asymmetry_bin_mtt_0_500*(1+(NP_charge_asymmetry_num_bin_mtt_0_500-NP_charge_asymmetry_deno_bin_mtt_0_500)/SM_charge_asymmetry_deno_bin_mtt_0_500);            
        
            
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
            
            return  SM_charge_asymmetry_bin_mtt_0_500*(1+(NP_charge_asymmetry_num_bin_mtt_0_500-NP_charge_asymmetry_deno_bin_mtt_0_500)/SM_charge_asymmetry_deno_bin_mtt_0_500);            
        
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
            
            return  SM_charge_asymmetry_bin_mtt_500_750*(1+(NP_charge_asymmetry_num_bin_mtt_500_750-NP_charge_asymmetry_deno_bin_mtt_500_750)/SM_charge_asymmetry_deno_bin_mtt_500_750);            
        
            
            
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
            
            return  SM_charge_asymmetry_bin_mtt_500_750*(1+(NP_charge_asymmetry_num_bin_mtt_500_750-NP_charge_asymmetry_deno_bin_mtt_500_750)/SM_charge_asymmetry_deno_bin_mtt_500_750);            
        
            
            
            
            
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
            
            return  SM_charge_asymmetry_bin_mtt_750_1000*(1+(NP_charge_asymmetry_num_bin_mtt_750_1000-NP_charge_asymmetry_deno_bin_mtt_750_1000)/SM_charge_asymmetry_deno_bin_mtt_750_1000);            
        
            
            
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
            
            return  SM_charge_asymmetry_bin_mtt_750_1000*(1+(NP_charge_asymmetry_num_bin_mtt_750_1000-NP_charge_asymmetry_deno_bin_mtt_750_1000)/SM_charge_asymmetry_deno_bin_mtt_750_1000);            
        
            
            
            
            
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
            
            return  SM_charge_asymmetry_bin_mtt_1000_1500*(1+(NP_charge_asymmetry_num_bin_mtt_1000_1500-NP_charge_asymmetry_deno_bin_mtt_1000_1500)/SM_charge_asymmetry_deno_bin_mtt_1000_1500);            
        
            
            
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
            
            return  SM_charge_asymmetry_bin_mtt_1000_1500*(1+(NP_charge_asymmetry_num_bin_mtt_1000_1500-NP_charge_asymmetry_deno_bin_mtt_1000_1500)/SM_charge_asymmetry_deno_bin_mtt_1000_1500);            
        
            
            
            
            
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
            
            return  SM_charge_asymmetry_bin_mtt_1500_3000*(1+(NP_charge_asymmetry_num_bin_mtt_1500_3000-NP_charge_asymmetry_deno_bin_mtt_1500_3000)/SM_charge_asymmetry_deno_bin_mtt_1500_3000);            

            
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
            
            return  SM_charge_asymmetry_bin_mtt_1500_3000*(1+(NP_charge_asymmetry_num_bin_mtt_1500_3000-NP_charge_asymmetry_deno_bin_mtt_1500_3000)/SM_charge_asymmetry_deno_bin_mtt_1500_3000);            

        }
        
    }     else if(b_min == 1500 && b_max == 2000){
        
        ////////////////////////////////
        //NEEDS TO BE UPDATED!!!!!!!!!//
        ////////////////////////////////
        
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
            
            return  SM_charge_asymmetry_bin_mtt_1500_2000*(1+(NP_charge_asymmetry_num_bin_mtt_1500_2000-NP_charge_asymmetry_deno_bin_mtt_1500_2000)/SM_charge_asymmetry_deno_bin_mtt_1500_2000);            

            
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
            
            return  SM_charge_asymmetry_bin_mtt_1500_2000*(1+(NP_charge_asymmetry_num_bin_mtt_1500_2000-NP_charge_asymmetry_deno_bin_mtt_1500_2000)/SM_charge_asymmetry_deno_bin_mtt_1500_2000);            

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
            
            return  SM_charge_asymmetry_bin_mtt_2000_2500*(1+(NP_charge_asymmetry_num_bin_mtt_2000_2500-NP_charge_asymmetry_deno_bin_mtt_2000_2500)/SM_charge_asymmetry_deno_bin_mtt_2000_2500);            

            
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
            
            return  SM_charge_asymmetry_bin_mtt_2000_2500*(1+(NP_charge_asymmetry_num_bin_mtt_2000_2500-NP_charge_asymmetry_deno_bin_mtt_2000_2500)/SM_charge_asymmetry_deno_bin_mtt_2000_2500);            

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
            
            return  SM_charge_asymmetry_bin_mtt_2500_3000*(1+(NP_charge_asymmetry_num_bin_mtt_2500_3000-NP_charge_asymmetry_deno_bin_mtt_2500_3000)/SM_charge_asymmetry_deno_bin_mtt_2500_3000);            

            
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
            
            return  SM_charge_asymmetry_bin_mtt_2500_3000*(1+(NP_charge_asymmetry_num_bin_mtt_2500_3000-NP_charge_asymmetry_deno_bin_mtt_2500_3000)/SM_charge_asymmetry_deno_bin_mtt_2500_3000);            

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
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
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





/////// Prospects of e+e- Colliders at 250 GeV //////////////////////////////////////////////////////////



sigma_250_bb_eP_P30_eM_M80::sigma_250_bb_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_250_bb_eP_P30_eM_M80");
}

double sigma_250_bb_eP_P30_eM_M80::computeThValue()
{
 
    
    //double sigma_250_bb_eP_P30_eM_M80_Madgraph = 3291;//fb
    double sigma_250_bb_eP_P30_eM_M80 = SM.getOptionalParameter("SM_sigma_250_bb_eP_P30_eM_M80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    //
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+9.3*((C_phiQm+C_phiQ3)/0.998)+9.3*(C_phiQ3/0.998)+1.5*(C_phib/0.998)+114.5*C_lqP+20.6*C_ld+2.9*C_ed)*sigma_250_bb_eP_P30_eM_M80/100;

            }
            else{
                return  (100+9.3*((C_phiQm+C_phiQ3)/0.998)+9.3*(C_phiQ3/0.998)+1.5*(C_phib/0.998)+114.5*C_lqP+20.6*C_ld+2.9*C_ed)*sigma_250_bb_eP_P30_eM_M80/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  (100+9.3*(C_phiQ1/0.998)+9.3*(C_phiQ3/0.998)+1.5*(C_phib/0.998)+114.5*C_lqP+20.6*C_ld+2.9*C_ed)*sigma_250_bb_eP_P30_eM_M80/100;
            }
            else{
                return  (100+9.3*(C_phiQ1/0.998)+9.3*(C_phiQ3/0.998)+1.5*(C_phib/0.998)+114.5*C_lqP+20.6*C_ld+2.9*C_ed)*sigma_250_bb_eP_P30_eM_M80/100;
            }
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
    double sigma_250_bb_eP_M30_eM_P80 = SM.getOptionalParameter("SM_sigma_250_bb_eP_M30_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+9.2*((C_phiQm+C_phiQ3)/0.998)+9.2*(C_phiQ3/0.998)-10.7*(C_phib/0.998)+22*C_lqP+4.1*C_ld+158*C_ed)*sigma_250_bb_eP_M30_eM_P80/100;

            }
            else{
                return  (100+9.2*((C_phiQm+C_phiQ3)/0.998)+9.2*(C_phiQ3/0.998)-10.7*(C_phib/0.998)+22*C_lqP+4.1*C_ld+158*C_ed)*sigma_250_bb_eP_M30_eM_P80/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  (100+9.2*(C_phiQ1/0.998)+9.2*(C_phiQ3/0.998)-10.7*(C_phib/0.998)+22*C_lqP+4.1*C_ld+158*C_ed)*sigma_250_bb_eP_M30_eM_P80/100;
            }
            else{
                return  (100+9.2*(C_phiQ1/0.998)+9.2*(C_phiQ3/0.998)-10.7*(C_phib/0.998)+22*C_lqP+4.1*C_ld+158*C_ed)*sigma_250_bb_eP_M30_eM_P80/100;
            }
    }
}



a_250_bb_eP_P30_eM_M80::a_250_bb_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_250_bb_eP_P30_eM_M80");
}

double a_250_bb_eP_P30_eM_M80::computeThValue()
{
    
    //double a_250_bb_eP_P30_eM_M80 = 69.8;//in percentage
    double a_250_bb_eP_P30_eM_M80 = SM.getOptionalParameter("SM_a_250_bb_eP_P30_eM_M80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    //We divide by 100 because the values introduced are in percentage
    
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_250_bb_eP_P30_eM_M80+(0.3*((C_phiQm+C_phiQ3)/0.998)+0.3*(C_phiQ3/0.998)-2.4*(C_phib/0.998)+12.9*C_lqP-17.4*C_ld+0.1*C_ed)/100;

            }
            else{
                return  a_250_bb_eP_P30_eM_M80+(0.3*((C_phiQm+C_phiQ3)/0.998)+0.3*(C_phiQ3/0.998)-2.4*(C_phib/0.998)+12.9*C_lqP-17.4*C_ld+0.1*C_ed)/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_250_bb_eP_P30_eM_M80+(0.3*(C_phiQ1/0.998)+0.3*(C_phiQ3/0.998)-2.4*(C_phib/0.998)+12.9*C_lqP-17.4*C_ld+0.1*C_ed)/100;
            }
            else{
                return  a_250_bb_eP_P30_eM_M80+(0.3*(C_phiQ1/0.998)+0.3*(C_phiQ3/0.998)-2.4*(C_phib/0.998)+12.9*C_lqP-17.4*C_ld+0.1*C_ed)/100;
            }
    }
}



a_250_bb_eP_M30_eM_P80::a_250_bb_eP_M30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_250_bb_eP_M30_eM_P80");
    
}

double a_250_bb_eP_M30_eM_P80::computeThValue()
{
    
    // double a_250_bb_eP_M30_eM_P80 = 36.5;// In percentage
    double a_250_bb_eP_M30_eM_P80 = SM.getOptionalParameter("SM_a_250_bb_eP_M30_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_250_bb_eP_M30_eM_P80+(-7.6*((C_phiQm+C_phiQ3)/0.998)-7.6*(C_phiQ3/0.998)-4.7*(C_phib/0.998)+7.8*C_lqP-4.1*C_ld+30*C_ed)/100;

            }
            else{
                return  a_250_bb_eP_M30_eM_P80+(-7.6*((C_phiQm+C_phiQ3)/0.998)-7.6*(C_phiQ3/0.998)-4.7*(C_phib/0.998)+7.8*C_lqP-4.1*C_ld+30*C_ed)/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_250_bb_eP_M30_eM_P80+(-7.6*(C_phiQ1/0.998)-7.6*(C_phiQ3/0.998)-4.7*(C_phib/0.998)+7.8*C_lqP-4.1*C_ld+30*C_ed)/100;
            }
            else{
                return  a_250_bb_eP_M30_eM_P80+(-7.6*(C_phiQ1/0.998)-7.6*(C_phiQ3/0.998)-4.7*(C_phib/0.998)+7.8*C_lqP-4.1*C_ld+30*C_ed)/100;
            }
    }
}



/////// Prospects of e+e- Colliders at 500 GeV //////////////////////////////////////////////////////////


sigma_500_bb_eP_P30_eM_M80::sigma_500_bb_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_500_bb_eP_P30_eM_M80");
}

double sigma_500_bb_eP_P30_eM_M80::computeThValue()
{
    
    //double sigma_500_bb_eP_P30_eM_M80 = 718.5;//fb
    double sigma_500_bb_eP_P30_eM_M80 = SM.getOptionalParameter("SM_sigma_500_bb_eP_P30_eM_M80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+8.9*((C_phiQm+C_phiQ3)/0.998)+8.9*(C_phiQ3/0.998)+1.6*(C_phib/0.998)+487*C_lqP+102*C_ld+13*C_ed)*sigma_500_bb_eP_P30_eM_M80/100;

            }
            else{
                return  (100+8.9*((C_phiQm+C_phiQ3)/0.998)+8.9*(C_phiQ3/0.998)+1.6*(C_phib/0.998)+487*C_lqP+102*C_ld+13*C_ed)*sigma_500_bb_eP_P30_eM_M80/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  (100+8.9*(C_phiQ1/0.998)+8.9*(C_phiQ3/0.998)+1.6*(C_phib/0.998)+487*C_lqP+102*C_ld+13*C_ed)*sigma_500_bb_eP_P30_eM_M80/100;
            }
            else{
                return  (100+8.9*(C_phiQ1/0.998)+8.9*(C_phiQ3/0.998)+1.6*(C_phib/0.998)+487*C_lqP+102*C_ld+13*C_ed)*sigma_500_bb_eP_P30_eM_M80/100;
            }
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
    double sigma_500_bb_eP_M30_eM_P80 = SM.getOptionalParameter("SM_sigma_500_bb_eP_M30_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+7.9*((C_phiQm+C_phiQ3)/0.998)+7.9*(C_phiQ3/0.998)-11*(C_phib/0.998)+97*C_lqP+20*C_ld+721*C_ed)*sigma_500_bb_eP_M30_eM_P80/100;

            }
            else{
                return  (100+7.9*((C_phiQm+C_phiQ3)/0.998)+7.9*(C_phiQ3/0.998)-11*(C_phib/0.998)+97*C_lqP+20*C_ld+721*C_ed)*sigma_500_bb_eP_M30_eM_P80/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  (100+7.9*(C_phiQ1/0.998)+7.9*(C_phiQ3/0.998)-11*(C_phib/0.998)+97*C_lqP+20*C_ld+721*C_ed)*sigma_500_bb_eP_M30_eM_P80/100;
            }
            else{
                return  (100+7.9*(C_phiQ1/0.998)+7.9*(C_phiQ3/0.998)-11*(C_phib/0.998)+97*C_lqP+20*C_ld+721*C_ed)*sigma_500_bb_eP_M30_eM_P80/100;
            }
    }
}



a_500_bb_eP_P30_eM_M80::a_500_bb_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_500_bb_eP_P30_eM_M80");
}

double a_500_bb_eP_P30_eM_M80::computeThValue()
{
    
    //double a_500_bb_eP_P30_eM_M80 = 68.4; // in percentage 
    double a_500_bb_eP_P30_eM_M80 = SM.getOptionalParameter("SM_a_500_bb_eP_P30_eM_M80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_500_bb_eP_P30_eM_M80+(0.45*((C_phiQm+C_phiQ3)/0.998)+0.45*(C_phiQ3/0.998)-2.6*(C_phib/0.998)+1.1*C_lqP-2.9*C_ld+0.46*C_ed)/100;

            }
            else{
                return  a_500_bb_eP_P30_eM_M80+(0.45*((C_phiQm+C_phiQ3)/0.998)+0.45*(C_phiQ3/0.998)-2.6*(C_phib/0.998)+1.1*C_lqP-2.9*C_ld+0.46*C_ed)/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_500_bb_eP_P30_eM_M80+(0.45*(C_phiQ1/0.998)+0.45*(C_phiQ3/0.998)-2.6*(C_phib/0.998)+1.1*C_lqP-2.9*C_ld+0.46*C_ed)/100;
            }
            else{
                return  a_500_bb_eP_P30_eM_M80+(0.45*(C_phiQ1/0.998)+0.45*(C_phiQ3/0.998)-2.6*(C_phib/0.998)+1.1*C_lqP-2.9*C_ld+0.46*C_ed)/100;
            }
    }
}



a_500_bb_eP_M30_eM_P80::a_500_bb_eP_M30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_500_bb_eP_M30_eM_P80");
}

double a_500_bb_eP_M30_eM_P80::computeThValue()
{
    
    //double a_500_bb_eP_M30_eM_P80 = 46.9; // in percentage 
    double a_500_bb_eP_M30_eM_P80 = SM.getOptionalParameter("SM_a_500_bb_eP_M30_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    

    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_500_bb_eP_M30_eM_P80+(-6.9*((C_phiQm+C_phiQ3)/0.998)-6.9*(C_phiQ3/0.998)-3.7*(C_phib/0.998)+6.7*C_lqP-5*C_ld+0.48*C_ed)/100;

            }
            else{
                return  a_500_bb_eP_M30_eM_P80+(-6.9*((C_phiQm+C_phiQ3)/0.998)-6.9*(C_phiQ3/0.998)-3.7*(C_phib/0.998)+6.7*C_lqP-5*C_ld+0.48*C_ed)/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_500_bb_eP_M30_eM_P80+(-6.9*(C_phiQ1/0.998)-6.9*(C_phiQ3/0.998)-3.7*(C_phib/0.998)+6.7*C_lqP-5*C_ld+0.48*C_ed)/100;
            }
            else{
                return  a_500_bb_eP_M30_eM_P80+(-6.9*(C_phiQ1/0.998)-6.9*(C_phiQ3/0.998)-3.7*(C_phib/0.998)+6.7*C_lqP-5*C_ld+0.48*C_ed)/100;
            }
    }
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
    double sigma_1000_bb_eP_P30_eM_M80 = SM.getOptionalParameter("SM_sigma_1000_bb_eP_P30_eM_M80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+8.8*((C_phiQm+C_phiQ3)/0.998)+8.8*(C_phiQ3/0.998)+1.7*(C_phib/0.998)+
                        1974*C_lqP+428*C_ld+54*C_ed)*sigma_1000_bb_eP_P30_eM_M80/100;

            }
            else{
                return  (100+8.8*((C_phiQm+C_phiQ3)/0.998)+8.8*(C_phiQ3/0.998)+1.7*(C_phib/0.998)+
                        1974*C_lqP+428*C_ld+54*C_ed)*sigma_1000_bb_eP_P30_eM_M80/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  (100+8.8*(C_phiQ1/0.998)+8.8*(C_phiQ3/0.998)+1.7*(C_phib/0.998)+
                        1974*C_lqP+428*C_ld+54*C_ed)*sigma_1000_bb_eP_P30_eM_M80/100;
            }
            else{
                return  (100+8.8*(C_phiQ1/0.998)+8.8*(C_phiQ3/0.998)+1.7*(C_phib/0.998)+
                        1974*C_lqP+428*C_ld+54*C_ed)*sigma_1000_bb_eP_P30_eM_M80/100;
            }
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
    double sigma_1000_bb_eP_M30_eM_P80 = SM.getOptionalParameter("SM_sigma_1000_bb_eP_M30_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+7.5*((C_phiQm+C_phiQ3)/0.998)+7.5*(C_phiQ3/0.998)-11*(C_phib/0.998)+398*C_lqP+89*C_ld+2975*C_ed)*sigma_1000_bb_eP_M30_eM_P80/100;

            }
            else{
                return  (100+7.5*((C_phiQm+C_phiQ3)/0.998)+7.5*(C_phiQ3/0.998)-11*(C_phib/0.998)+398*C_lqP+89*C_ld+2975*C_ed)*sigma_1000_bb_eP_M30_eM_P80/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  (100+7.5*(C_phiQ1/0.998)+7.5*(C_phiQ3/0.998)-11*(C_phib/0.998)+398*C_lqP+89*C_ld+2975*C_ed)*sigma_1000_bb_eP_M30_eM_P80/100;
            }
            else{
                return  (100+7.5*(C_phiQ1/0.998)+7.5*(C_phiQ3/0.998)-11*(C_phib/0.998)+398*C_lqP+89*C_ld+2975*C_ed)*sigma_1000_bb_eP_M30_eM_P80/100;
            }
    }
}



a_1000_bb_eP_P30_eM_M80::a_1000_bb_eP_P30_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_1000_bb_eP_P30_eM_M80");
}

double a_1000_bb_eP_P30_eM_M80::computeThValue()
{
    
    //double a_1000_bb_eP_P30_eM_M80 = 67.5 ;//in percentage
    double a_1000_bb_eP_P30_eM_M80 = SM.getOptionalParameter("SM_a_1000_bb_eP_P30_eM_M80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    //double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_1000_bb_eP_P30_eM_M80+(0.48*((C_phiQm+C_phiQ3)/0.998)+0.48*(C_phiQ3/0.998)-2.6*(C_phib/0.998)+0.06*C_ld+0.08*C_ed)/100;

            }
            else{
                return  a_1000_bb_eP_P30_eM_M80+(0.48*((C_phiQm+C_phiQ3)/0.998)+0.48*(C_phiQ3/0.998)-2.6*(C_phib/0.998)+0.06*C_ld+0.08*C_ed)/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_1000_bb_eP_P30_eM_M80+(0.48*(C_phiQ1/0.998)+0.48*(C_phiQ3/0.998)-2.6*(C_phib/0.998)+0.06*C_ld+0.08*C_ed)/100;
            }
            else{
                return  a_1000_bb_eP_P30_eM_M80+(0.48*(C_phiQ1/0.998)+0.48*(C_phiQ3/0.998)-2.6*(C_phib/0.998)+0.06*C_ld+0.08*C_ed)/100;
            }
    }
}


a_1000_bb_eP_M30_eM_P80::a_1000_bb_eP_M30_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_1000_bb_eP_M30_eM_P80");
}

double a_1000_bb_eP_M30_eM_P80::computeThValue()
{
    
    // double a_1000_bb_eP_M30_eM_P80 = 49.5;//in percentage
    double a_1000_bb_eP_M30_eM_P80 = SM.getOptionalParameter("SM_a_1000_bb_eP_M30_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    //double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_1000_bb_eP_M30_eM_P80+(-6.6*((C_phiQm+C_phiQ3)/0.998)-6.6*(C_phiQ3/0.998)-3.5*(C_phib/0.998)+0.22*C_lqP-0.23*C_ld)/100;

            }
            else{
                return  a_1000_bb_eP_M30_eM_P80+(-6.6*((C_phiQm+C_phiQ3)/0.998)-6.6*(C_phiQ3/0.998)-3.5*(C_phib/0.998)+0.22*C_lqP-0.23*C_ld)/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_1000_bb_eP_M30_eM_P80+(-6.6*(C_phiQ1/0.998)-6.6*(C_phiQ3/0.998)-3.5*(C_phib/0.998)+0.22*C_lqP-0.23*C_ld)/100;
            }
            else{
                return  a_1000_bb_eP_M30_eM_P80+(-6.6*(C_phiQ1/0.998)-6.6*(C_phiQ3/0.998)-3.5*(C_phib/0.998)+0.22*C_lqP-0.23*C_ld)/100;
            }
    }
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
    double sigma_380_bb_eP_0_eM_P80 = SM.getOptionalParameter("SM_sigma_380_bb_eP_0_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+8.3*((C_phiQm+C_phiQ3)/0.998)+8.3*(C_phiQ3/0.998)-9.1*(C_phib/0.998)+87.1*C_lqP+18.7*C_ld+347.1*C_ed)*sigma_380_bb_eP_0_eM_P80/100;

            }
            else{
                return  (100+8.3*((C_phiQm+C_phiQ3)/0.998)+8.3*(C_phiQ3/0.998)-9.1*(C_phib/0.998)+87.1*C_lqP+18.7*C_ld+347.1*C_ed)*sigma_380_bb_eP_0_eM_P80/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  (100+8.3*(C_phiQ1/0.998)+8.3*(C_phiQ3/0.998)-9.1*(C_phib/0.998)+87.1*C_lqP+18.7*C_ld+347.1*C_ed)*sigma_380_bb_eP_0_eM_P80/100;
            }
            else{
                return  (100+8.3*(C_phiQ1/0.998)+8.3*(C_phiQ3/0.998)-9.1*(C_phib/0.998)+87.1*C_lqP+18.7*C_ld+347.1*C_ed)*sigma_380_bb_eP_0_eM_P80/100;
            }
    }
}




sigma_380_bb_eP_0_eM_M80::sigma_380_bb_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_380_bb_eP_0_eM_M80");
}

double sigma_380_bb_eP_0_eM_M80::computeThValue()
{
    
    //double sigma_380_bb_eP_0_eM_M80 = 998;//fb
    double sigma_380_bb_eP_0_eM_M80 = SM.getOptionalParameter("SM_sigma_380_bb_eP_0_eM_M80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+9*((C_phiQm+C_phiQ3)/0.998)+9*(C_phiQ3/0.998)+1.5*(C_phib/0.998)+274*C_lqP+57*C_ld+13*C_ed)*sigma_380_bb_eP_0_eM_M80/100;

            }
            else{
                return  (100+9*((C_phiQm+C_phiQ3)/0.998)+9*(C_phiQ3/0.998)+1.5*(C_phib/0.998)+274*C_lqP+57*C_ld+13*C_ed)*sigma_380_bb_eP_0_eM_M80/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  (100+9*(C_phiQ1/0.998)+9*(C_phiQ3/0.998)+1.5*(C_phib/0.998)+274*C_lqP+57*C_ld+13*C_ed)*sigma_380_bb_eP_0_eM_M80/100;
            }
            else{
                return  (100+9*(C_phiQ1/0.998)+9*(C_phiQ3/0.998)+1.5*(C_phib/0.998)+274*C_lqP+57*C_ld+13*C_ed)*sigma_380_bb_eP_0_eM_M80/100;
            }
    }
}



a_380_bb_eP_0_eM_P80::a_380_bb_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_380_bb_eP_0_eM_P80");
}

double a_380_bb_eP_0_eM_P80::computeThValue()
{
    
    //double a_380_bb_eP_0_eM_P80 = 47.9;//in percentage
    double a_380_bb_eP_0_eM_P80 = SM.getOptionalParameter("SM_a_380_bb_eP_0_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_380_bb_eP_0_eM_P80+(-6.1*((C_phiQm+C_phiQ3)/0.998)-6.1*(C_phiQ3/0.998)-3.4*(C_phib/0.998)+12.2*C_lqP-8.4*C_ld+2.8*C_ed)/100;

            }
            else{
                return  a_380_bb_eP_0_eM_P80+(-6.1*((C_phiQm+C_phiQ3)/0.998)-6.1*(C_phiQ3/0.998)-3.4*(C_phib/0.998)+12.2*C_lqP-8.4*C_ld+2.8*C_ed)/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_380_bb_eP_0_eM_P80+(-6.1*(C_phiQ1/0.998)-6.1*(C_phiQ3/0.998)-3.4*(C_phib/0.998)+12.2*C_lqP-8.4*C_ld+2.8*C_ed)/100;
            }
            else{
                return  a_380_bb_eP_0_eM_P80+(-6.1*(C_phiQ1/0.998)-6.1*(C_phiQ3/0.998)-3.4*(C_phib/0.998)+12.2*C_lqP-8.4*C_ld+2.8*C_ed)/100;
            }
    }
}




a_380_bb_eP_0_eM_M80::a_380_bb_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_380_bb_eP_0_eM_M80");
}

double a_380_bb_eP_0_eM_M80::computeThValue()
{
    
    //double a_380_bb_eP_0_eM_M80 = 67.93;//in percentage
    double a_380_bb_eP_0_eM_M80 = SM.getOptionalParameter("SM_a_380_bb_eP_0_eM_M80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    


    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_380_bb_eP_0_eM_M80+(0.34*((C_phiQm+C_phiQ3)/0.998)+0.34*(C_phiQ3/0.998)-2.6*(C_phib/0.998)+12*C_lqP-9.4*C_ld+0.5*C_ed)/100;

            }
            else{
                return  a_380_bb_eP_0_eM_M80+(0.34*((C_phiQm+C_phiQ3)/0.998)+0.34*(C_phiQ3/0.998)-2.6*(C_phib/0.998)+12*C_lqP-9.4*C_ld+0.5*C_ed)/100;            
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_380_bb_eP_0_eM_M80+(0.34*(C_phiQ1/0.998)+0.34*(C_phiQ3/0.998)-2.6*(C_phib/0.998)+12*C_lqP-9.4*C_ld+0.5*C_ed)/100;            
            }
            else{
                return  a_380_bb_eP_0_eM_M80+(0.34*(C_phiQ1/0.998)+0.34*(C_phiQ3/0.998)-2.6*(C_phib/0.998)+12*C_lqP-9.4*C_ld+0.5*C_ed)/100;            
            }
    }

}







/////// Prospects of Linear Colders at 1400 GeV //////////////////////////////////////////////////////////


sigma_1400_bb_eP_0_eM_P80::sigma_1400_bb_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1400_bb_eP_0_eM_P80");
}

double sigma_1400_bb_eP_0_eM_P80::computeThValue()
{
    
    //double sigma_1400_bb_eP_0_eM_P80 = 23.77;//fb
    double sigma_1400_bb_eP_0_eM_P80 = SM.getOptionalParameter("SM_sigma_1400_bb_eP_0_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+7.7*((C_phiQm+C_phiQ3)/0.998)+7.7*(C_phiQ3/0.998)-9.2*(C_phib/0.998)+1235*C_lqP+277*C_ld+4990*C_ed)*sigma_1400_bb_eP_0_eM_P80/100;

            }
            else{
                return  (100+7.7*((C_phiQm+C_phiQ3)/0.998)+7.7*(C_phiQ3/0.998)-9.2*(C_phib/0.998)+1235*C_lqP+277*C_ld+4990*C_ed)*sigma_1400_bb_eP_0_eM_P80/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  (100+7.7*(C_phiQ1/0.998)+7.7*(C_phiQ3/0.998)-9.2*(C_phib/0.998)+1235*C_lqP+277*C_ld+4990*C_ed)*sigma_1400_bb_eP_0_eM_P80/100;
            }
            else{
                return  (100+7.7*(C_phiQ1/0.998)+7.7*(C_phiQ3/0.998)-9.2*(C_phib/0.998)+1235*C_lqP+277*C_ld+4990*C_ed)*sigma_1400_bb_eP_0_eM_P80/100;
            }
    }
}




sigma_1400_bb_eP_0_eM_M80::sigma_1400_bb_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_sigma_1400_bb_eP_0_eM_M80");
}

double sigma_1400_bb_eP_0_eM_M80::computeThValue()
{
    
    //double sigma_1400_bb_eP_0_eM_M80 = 68.85;//fb
    double sigma_1400_bb_eP_0_eM_M80 = SM.getOptionalParameter("SM_sigma_1400_bb_eP_0_eM_M80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
        
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+8.7*((C_phiQm+C_phiQ3)/0.998)+8.7*(C_phiQ3/0.998)+1.5*(C_phib/0.998)+3843*C_lqP+846*C_ld+182.5*C_ed)*sigma_1400_bb_eP_0_eM_M80/100;

            }
            else{
                return  (100+8.7*((C_phiQm+C_phiQ3)/0.998)+8.7*(C_phiQ3/0.998)+1.5*(C_phib/0.998)+3843*C_lqP+846*C_ld+182.5*C_ed)*sigma_1400_bb_eP_0_eM_M80/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  (100+8.7*(C_phiQ1/0.998)+8.7*(C_phiQ3/0.998)+1.5*(C_phib/0.998)+3843*C_lqP+846*C_ld+182.5*C_ed)*sigma_1400_bb_eP_0_eM_M80/100;
            }
            else{
                return  (100+8.7*(C_phiQ1/0.998)+8.7*(C_phiQ3/0.998)+1.5*(C_phib/0.998)+3843*C_lqP+846*C_ld+182.5*C_ed)*sigma_1400_bb_eP_0_eM_M80/100;
            }
    }
}







a_1400_bb_eP_0_eM_P80::a_1400_bb_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_1400_bb_eP_0_eM_P80");
}

double a_1400_bb_eP_0_eM_P80::computeThValue()
{
    
    //double a_1400_bb_eP_0_eM_P80 = 51.8; //in percentage
    double a_1400_bb_eP_0_eM_P80 = SM.getOptionalParameter("SM_a_1400_bb_eP_0_eM_P80");
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
                return  a_1400_bb_eP_0_eM_P80+(-5.6*((C_phiQm+C_phiQ3)/0.998)-5.6*(C_phiQ3/0.998)-3.1*(C_phib/0.998))/100;

            }
            else{
                return  a_1400_bb_eP_0_eM_P80+(-5.6*((C_phiQm+C_phiQ3)/0.998)-5.6*(C_phiQ3/0.998)-3.1*(C_phib/0.998))/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_1400_bb_eP_0_eM_P80+(-5.6*(C_phiQ1/0.998)-5.6*(C_phiQ3/0.998)-3.1*(C_phib/0.998))/100;
            }
            else{
                return  a_1400_bb_eP_0_eM_P80+(-5.6*(C_phiQ1/0.998)-5.6*(C_phiQ3/0.998)-3.1*(C_phib/0.998))/100;
            }
    }
}




a_1400_bb_eP_0_eM_M80::a_1400_bb_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_1400_bb_eP_0_eM_M80");
}

double a_1400_bb_eP_0_eM_M80::computeThValue()
{
    
    // double a_1400_bb_eP_0_eM_M80 = 67.4;// in percentage
    double a_1400_bb_eP_0_eM_M80 = SM.getOptionalParameter("SM_a_1400_bb_eP_0_eM_M80");
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
                return  a_1400_bb_eP_0_eM_M80+(0.4*((C_phiQm+C_phiQ3)/0.998)+0.4*(C_phiQ3/0.998)-2.7*(C_phib/0.998))/100;

            }
            else{
                return  a_1400_bb_eP_0_eM_M80+(0.4*((C_phiQm+C_phiQ3)/0.998)+0.4*(C_phiQ3/0.998)-2.7*(C_phib/0.998))/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_1400_bb_eP_0_eM_M80+(0.4*(C_phiQ1/0.998)+0.4*(C_phiQ3/0.998)-2.7*(C_phib/0.998))/100;
            }
            else{
                return  a_1400_bb_eP_0_eM_M80+(0.4*(C_phiQ1/0.998)+0.4*(C_phiQ3/0.998)-2.7*(C_phib/0.998))/100;
            }
    }
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
    double sigma_3000_bb_eP_0_eM_P80 = SM.getOptionalParameter("SM_sigma_3000_bb_eP_0_eM_P80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    //double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    //double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    //double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    //double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    


    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+7.6*((C_phiQm+C_phiQ3)/0.998)+7.6*(C_phiQ3/0.998)-9.3*(C_phib/0.998)+5684*C_lqP+1273*C_ld+25000*C_ed)*sigma_3000_bb_eP_0_eM_P80/100;

            }
            else{
                return  (100+7.6*((C_phiQm+C_phiQ3)/0.998)+7.6*(C_phiQ3/0.998)-9.3*(C_phib/0.998)+5684*C_lqP+1273*C_ld+25000*C_ed)*sigma_3000_bb_eP_0_eM_P80/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  (100+7.6*(C_phiQ1/0.998)+7.6*(C_phiQ3/0.998)-9.3*(C_phib/0.998)+5684*C_lqP+1273*C_ld+25000*C_ed)*sigma_3000_bb_eP_0_eM_P80/100;
            }
            else{
                return  (100+7.6*(C_phiQ1/0.998)+7.6*(C_phiQ3/0.998)-9.3*(C_phib/0.998)+5684*C_lqP+1273*C_ld+25000*C_ed)*sigma_3000_bb_eP_0_eM_P80/100;
            }
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
    double sigma_3000_bb_eP_0_eM_M80 = SM.getOptionalParameter("SM_sigma_3000_bb_eP_0_eM_M80");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    //double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    //double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    //double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    //double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (100+8.8*((C_phiQm+C_phiQ3)/0.998)+8.8*(C_phiQ3/0.998)+1.6*(C_phib/0.998)+17193*C_lqP+4002*C_ld+858*C_ed)*sigma_3000_bb_eP_0_eM_M80/100;

            }
            else{
                return  (100+8.8*((C_phiQm+C_phiQ3)/0.998)+8.8*(C_phiQ3/0.998)+1.6*(C_phib/0.998)+17193*C_lqP+4002*C_ld+858*C_ed)*sigma_3000_bb_eP_0_eM_M80/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  (100+8.8*(C_phiQ1/0.998)+8.8*(C_phiQ3/0.998)+1.6*(C_phib/0.998)+17193*C_lqP+4002*C_ld+858*C_ed)*sigma_3000_bb_eP_0_eM_M80/100;
            }
            else{
                return  (100+8.8*(C_phiQ1/0.998)+8.8*(C_phiQ3/0.998)+1.6*(C_phib/0.998)+17193*C_lqP+4002*C_ld+858*C_ed)*sigma_3000_bb_eP_0_eM_M80/100;
            }
    }
}




a_3000_bb_eP_0_eM_P80::a_3000_bb_eP_0_eM_P80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_3000_bb_eP_0_eM_P80");
}

double a_3000_bb_eP_0_eM_P80::computeThValue()
{
    
    //double a_3000_bb_eP_0_eM_P80 = 52.8; //in percentage
    double a_3000_bb_eP_0_eM_P80 = SM.getOptionalParameter("SM_a_3000_bb_eP_0_eM_P80");
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
                return  a_3000_bb_eP_0_eM_P80+(-5.5*((C_phiQm+C_phiQ3)/0.998)-5.5*(C_phiQ3/0.998)-3.1*(C_phib/0.998))/100;

            }
            else{
                return  a_3000_bb_eP_0_eM_P80+(-5.5*((C_phiQm+C_phiQ3)/0.998)-5.5*(C_phiQ3/0.998)-3.1*(C_phib/0.998))/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_3000_bb_eP_0_eM_P80+(-5.5*(C_phiQ1/0.998)-5.5*(C_phiQ3/0.998)-3.1*(C_phib/0.998))/100;
            }
            else{
                return  a_3000_bb_eP_0_eM_P80+(-5.5*(C_phiQ1/0.998)-5.5*(C_phiQ3/0.998)-3.1*(C_phib/0.998))/100;
            }
    }
}




a_3000_bb_eP_0_eM_M80::a_3000_bb_eP_0_eM_M80(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_3000_bb_eP_0_eM_M80");
}

double a_3000_bb_eP_0_eM_M80::computeThValue()
{
    
    // double a_3000_bb_eP_0_eM_M80 = 67.6;// in percentage
    double a_3000_bb_eP_0_eM_M80 = SM.getOptionalParameter("SM_a_3000_bb_eP_0_eM_M80");
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
                return  a_3000_bb_eP_0_eM_M80+(0.42*((C_phiQm+C_phiQ3)/0.998)+0.42*(C_phiQ3/0.998)-2.8*(C_phib/0.998))/100;

            }
            else{
                return  a_3000_bb_eP_0_eM_M80+(0.42*((C_phiQm+C_phiQ3)/0.998)+0.42*(C_phiQ3/0.998)-2.8*(C_phib/0.998))/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_3000_bb_eP_0_eM_M80+(0.42*(C_phiQ1/0.998)+0.42*(C_phiQ3/0.998)-2.8*(C_phib/0.998))/100;
            }
            else{
                return  a_3000_bb_eP_0_eM_M80+(0.42*(C_phiQ1/0.998)+0.42*(C_phiQ3/0.998)-2.8*(C_phib/0.998))/100;
            }
    }
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
    
    //double sigma_240_bb_Madgraph = 1921;//fb
    double sigma_240_bb = SM.getOptionalParameter("SM_sigma_240_bb");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  sigma_240_bb/100*(100+9.3*((C_phiQm+C_phiQ3)/0.998)+9.3*(C_phiQ3/0.998)-1.4*(C_phib/0.998)+85*C_lqP+15*C_ld+36*C_ed);

            }
            else{
                return  sigma_240_bb/100*(100+9.3*((C_phiQm+C_phiQ3)/0.998)+9.3*(C_phiQ3/0.998)-1.4*(C_phib/0.998)+85*C_lqP+15*C_ld+36*C_ed);
            }
    }
    
    else{
            if(flag_Quadratic){
                return  sigma_240_bb/100*(100+9.3*(C_phiQ1/0.998)+9.3*(C_phiQ3/0.998)-1.4*(C_phib/0.998)+85*C_lqP+15*C_ld+36*C_ed);
            }
            else{
                return  sigma_240_bb/100*(100+9.3*(C_phiQ1/0.998)+9.3*(C_phiQ3/0.998)-1.4*(C_phib/0.998)+85*C_lqP+15*C_ld+36*C_ed);
            }
    }
}





a_240_bb::a_240_bb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{
    setParametersForObservable(make_vector<std::string>() << "SM_a_240_bb" );
}

double a_240_bb::computeThValue()
{
    
    //double a_240_bb = 58.4;//in percentage
    double a_240_bb = SM.getOptionalParameter("SM_a_240_bb");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_240_bb + (-1.5*((C_phiQm+C_phiQ3)/0.998)-1.5*(C_phiQ3/0.998)-2.3*(C_phib/0.998)+14*C_lqP-14*C_ld+3*C_ed)/100;

            }
            else{
                return  a_240_bb + (-1.5*((C_phiQm+C_phiQ3)/0.998)-1.5*(C_phiQ3/0.998)-2.3*(C_phib/0.998)+14*C_lqP-14*C_ld+3*C_ed)/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_240_bb + (-1.5*(C_phiQ1/0.998)-1.5*(C_phiQ3/0.998)-2.3*(C_phib/0.998)+14*C_lqP-14*C_ld+3*C_ed)/100;
            }
            else{
                return  a_240_bb + (-1.5*(C_phiQ1/0.998)-1.5*(C_phiQ3/0.998)-2.3*(C_phib/0.998)+14*C_lqP-14*C_ld+3*C_ed)/100;
            }
    }
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
    double sigma_360_bb = SM.getOptionalParameter("SM_sigma_360_bb");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  sigma_360_bb/100*(100+8.9*((C_phiQm+C_phiQ3)/0.998)+8.9*(C_phiQ3/0.998)-1.3*(C_phib/0.998)+202*C_lqP+41*C_ld+87*C_ed);

            }
            else{
                return  sigma_360_bb/100*(100+8.9*((C_phiQm+C_phiQ3)/0.998)+8.9*(C_phiQ3/0.998)-1.3*(C_phib/0.998)+202*C_lqP+41*C_ld+87*C_ed);
            }
    }
    
    else{
            if(flag_Quadratic){
                return  sigma_360_bb/100*(100+8.9*(C_phiQ1/0.998)+8.9*(C_phiQ3/0.998)-1.3*(C_phib/0.998)+202*C_lqP+41*C_ld+87*C_ed);
            }
            else{
                return  sigma_360_bb/100*(100+8.9*(C_phiQ1/0.998)+8.9*(C_phiQ3/0.998)-1.3*(C_phib/0.998)+202*C_lqP+41*C_ld+87*C_ed);
            }
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
    double a_360_bb = SM.getOptionalParameter("SM_a_360_bb");
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();

    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  a_360_bb + (-1.3*((C_phiQm+C_phiQ3)/0.998)-1.3*(C_phiQ3/0.998)-2.4*(C_phib/0.998)+19*C_lqP-11*C_ld+2.3*C_ed)/100;

            }
            else{
                return  a_360_bb + (-1.3*((C_phiQm+C_phiQ3)/0.998)-1.3*(C_phiQ3/0.998)-2.4*(C_phib/0.998)+19*C_lqP-11*C_ld+2.3*C_ed)/100;
            }
    }
    
    else{
            if(flag_Quadratic){
                return  a_360_bb + (-1.3*(C_phiQ1/0.998)-1.3*(C_phiQ3/0.998)-2.4*(C_phib/0.998)+19*C_lqP-11*C_ld+2.3*C_ed)/100;
            }
            else{
                return  a_360_bb + (-1.3*(C_phiQ1/0.998)-1.3*(C_phiQ3/0.998)-2.4*(C_phib/0.998)+19*C_lqP-11*C_ld+2.3*C_ed)/100;
            }
    }
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


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
#include <boost/bind.hpp>


const std::string NPSMEFT6dtopquark::NPSMEFT6dtopquarkVars[NNPSMEFT6dtopquarkVars]
        = {"C_phit","C_phiQ3","C_phiQ1","C_tW","C_tB","C_tphi","C_phib","C_bW","C_bB"};


NPSMEFT6dtopquark::NPSMEFT6dtopquark()
: NPbase()
{
}



void NPSMEFT6dtopquark::setParameter(const std::string name, const double& value)
{
    if (name.compare("C_phit") == 0)
        C_phit = value;
    else if (name.compare("C_phiQ3") == 0)
        C_phiQ3 = value;
    else if (name.compare("C_phiQ1") == 0)
        C_phiQ1 = value;
    else if (name.compare("C_tW") == 0)
        C_tW = value;
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
    else
        res = StandardModel::setFlag(name,value);

    return(res);
}

//Observables from LEP1



Rb_NPSMEFT6dtopquark::Rb_NPSMEFT6dtopquark(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};


double Rb_NPSMEFT6dtopquark::computeThValue()
{
    double smlep_bb = 0.21579;
    double lep_bb_madgraph = 0.22;
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    if(flag_Quadratic){
        return smlep_bb + (0.023*C_phiQ3+0.023*C_phiQ1-0.005*C_phib-0.0033*C_bW+0.0018*C_bW*C_bW+0.0015*C_bW*C_bB-0.0019*C_bB)*(smlep_bb/lep_bb_madgraph);
    }
    else{
        return smlep_bb + (0.023*C_phiQ3+0.023*C_phiQ1-0.005*C_phib-0.0033*C_bW-0.0019*C_bB)*(smlep_bb/lep_bb_madgraph);
    }
}


AFBLR::AFBLR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double AFBLR::computeThValue()
{
    double smlep_alr = 0.9347*0.75;
    double lep_alr_madgraph = 0.66;
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    if(flag_Quadratic){
    return smlep_alr + (0.008*C_phiQ3+0.008*C_phiQ1+0.034*C_phib+0.0056*C_bW*C_bW-0.002*C_phib*C_bB+0.002*C_bB*C_bB+0.0065*C_bW*C_bB)*(smlep_alr/lep_alr_madgraph);
    }
    else{
        return smlep_alr + (0.008*C_phiQ3+0.008*C_phiQ1+0.034*C_phib)*(smlep_alr/lep_alr_madgraph);
    }
}

//Observables from LHC Run 2 (also prospects at High lumniosity)


sigmattZ::sigmattZ(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmattZ::computeThValue()
{
    double smxttz = 0.88;
    double xttz_madgraph = 0.5887;
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    if(flag_Quadratic){
        return ( smxttz + (0.041*C_phit+0.0024*C_phit*C_phit+0.0018*C_phit*C_phiQ3+0.0018*C_phit*C_phiQ1-0.066*C_phiQ3
            +0.0035*C_phiQ3*C_phiQ3+0.066*C_phiQ1+0.005*C_phiQ1*C_phiQ1-0.0001*C_phiQ3*C_tB+0.0001*C_phiQ1*C_tB+0.00068*C_tW
            +0.018*C_tW*C_tW+0.0001*C_tW*C_phiQ3-0.0002*C_tW*C_phiQ1+0.00024*C_tB+0.0016*C_tB*C_tB-0.01*C_tW*C_tB)*(smxttz/xttz_madgraph));
    }
    else{
        return ( smxttz + (0.041*C_phit-0.066*C_phiQ3+0.066*C_phiQ1+0.00068*C_tW+0.00024*C_tB)*(smxttz/xttz_madgraph));
    }
}

sigmattA_1::sigmattA_1(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmattA_1::computeThValue()
{   
    double xtta_madgraph = 2.18;
    double smxtta = 0.495;
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    if(flag_Quadratic){
        return (smxtta + (0.015*C_tW+0.007*C_tW*C_tW+0.015*C_tB
            +0.007*C_tB*C_tB+0.014*C_tW*C_tB)*(smxtta/xtta_madgraph));
    }
    else{
        return (smxtta + (0.015*C_tW+0.015*C_tB)*(smxtta/xtta_madgraph)); 
    }
}

sigmattA_2::sigmattA_2(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmattA_2::computeThValue()
{
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double xtta_madgraph_180_pt_300 = 0.068;
    double xtta_madgraph = 2.18;
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    if(flag_Quadratic){
        return (xtta_madgraph_180_pt_300 +(0.0015*C_tW+0.0017*C_tW*C_tW+0.0014*C_tB+
                0.0016*C_tB*C_tB+0.0032*C_tW*C_tB))/(xtta_madgraph + (0.015*C_tW+0.007*C_tW*C_tW+0.015*C_tB
            +0.007*C_tB*C_tB+0.014*C_tW*C_tB));
    }
    else{
        return (xtta_madgraph_180_pt_300 +(0.0015*C_tW+0.0014*C_tB))/(xtta_madgraph + (0.015*C_tW+0.015*C_tB));
    }
}


sigmattH::sigmattH(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmattH::computeThValue()
{
    double smxtth = 0.507;
    double xtth_madgraph = 0.4;
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    if(flag_Quadratic){
        return  smxtth + (0.0007*C_tW+0.00067*C_tW*C_tW-0.049*C_tphi+0.0015*C_tphi*C_tphi)*(smxtth/xtth_madgraph);
    }
    else{
        return  smxtth + (0.0007*C_tW-0.049*C_tphi)*(smxtth/xtth_madgraph);
    }
}

sigmattW::sigmattW(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmattW::computeThValue()
{
    double smxttw = 0.60;
    double xttw_madgraph = 0.35;
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    if(flag_Quadratic){
        return smxttw + (0.00036*C_phiQ3-0.0003*C_phiQ1+0.00012*C_phiQ3*C_phiQ1+0.0027*C_tW
            +0.0032*C_tW*C_tW-0.00011*C_tW*C_phiQ1)*(smxttw/xttw_madgraph);
    }
    else{
        return smxttw + (0.00036*C_phiQ3-0.0003*C_phiQ1+0.0027*C_tW)*(smxttw/xttw_madgraph);
    }
}

sigmatq::sigmatq(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmatq::computeThValue()
{
    double smxst = 216.99;
    double xst_madgraph = 44.14;
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    if(flag_Quadratic){
        return smxst + (5.26*C_phiQ3+0.16*C_phiQ3*C_phiQ3+1.52*C_tW+0.31*C_tW*C_tW+0.1*C_phiQ3*C_tW)*(smxst/xst_madgraph);
    }
    else{
        return smxst + (5.26*C_phiQ3+1.52*C_tW)*(smxst/xst_madgraph);
    }
}

sigmatW::sigmatW(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmatW::computeThValue()
{
    double smxgbtw = 71.7;
    double xgbtw_madgraph = 13.5;
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    if(flag_Quadratic){
        return smxgbtw + (1.61*C_phiQ3+0.05*C_phiQ3*C_phiQ3-0.74*C_tW+0.135*C_tW*C_tW-0.046*C_phiQ3*C_tW)*(smxgbtw/xgbtw_madgraph);
    }
    else{
        return smxgbtw + (1.61*C_phiQ3-0.74*C_tW)*(smxgbtw/xgbtw_madgraph);
    }
}

sigmatqZ::sigmatqZ(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmatqZ::computeThValue()
{
    double smxztq = 0.8;
    double xztq_madgraph = 0.48;
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    if(flag_Quadratic){
        return smxztq + (0.0029*C_phit+0.0005*C_phit*C_phiQ3-0.0005*C_phit*C_phiQ1
            +0.092*C_phiQ3+0.014*C_phiQ3*C_phiQ3+0.0003*C_phiQ3*C_bW+0.001*C_phiQ3*C_tW+0.01*C_phiQ1
            +0.0003*C_phiQ1*C_tW+0.0017*C_phiQ3*C_phiQ1+0.007*C_tW+0.016*C_tW*C_tW-0.0003*C_phib
            -0.0002*C_phib*C_bW+0.026*C_bW*C_bW)*(smxztq/xztq_madgraph);
    }
    else{
        return smxztq + (0.0029*C_phit+0.092*C_phiQ3+0.01*C_phiQ1+0.007*C_tW-0.0003*C_phib)*(smxztq/xztq_madgraph);
    }
}




F0::F0(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double F0::computeThValue()
{
    double smF0 = 0.6978;
    double F0_madgraph = 0.699;
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    if(flag_Quadratic){
        return smF0 + (-0.04*C_tW+0.0025*C_tW*C_tW)*(smF0/F0_madgraph);
    }
    else{
        return smF0 + (-0.04*C_tW)*(smF0/F0_madgraph);
    }
}


FL::FL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double FL::computeThValue()
{
    double smFL = 0.3022;
    double FL_madgraph = 0.301;
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    if(flag_Quadratic){
        return smFL + (0.04*C_tW-0.0025*C_tW*C_tW)*(smFL/FL_madgraph);
    }
    else{
        return smFL + (0.04*C_tW)*(smFL/FL_madgraph);
    }
}


//Prospects of Linear Colders at 250 GeV
//250 bb observables

sigma_250_bb_eLpR::sigma_250_bb_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigma_250_bb_eLpR::computeThValue()
{
   // double sigma_250_bb_eLpR_madgraph = 3.29;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
        return (0.31*C_phiQ3+0.31*C_phiQ1+0.05*C_phib+0.11*C_bW+0.052*C_bB);
}



a_250_bb_eLpR::a_250_bb_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double a_250_bb_eLpR::computeThValue()
{
    //double a_250_bb_eLpR_madgraph = 69.6;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
        return ((0.4*C_phiQ3+0.3*C_phiQ1-2.2*C_phib-2.3*C_bW-0.6504*C_bB));
}


sigma_250_bb_eRpL::sigma_250_bb_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigma_250_bb_eRpL::computeThValue()
{
    //double sigma_250_bb_eRpL_madgraph = 1.02;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
        return (0.094*C_phiQ3+0.094*C_phiQ1-0.11*C_phib+0.007*C_bW+0.017*C_bB);
}




a_250_bb_eRpL::a_250_bb_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double a_250_bb_eRpL::computeThValue()
{
    //double a_250_bb_eRpL_madgraph = 35.9;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
        return (-7.8*C_phiQ3-7.7*C_phiQ1-4.5*C_phib-0.5*C_bW+2.4304*C_bB);
}


//Prospects of Linear Colders at 500 GeV
//500 bb observables


sigma_500_bb_eLpR::sigma_500_bb_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigma_500_bb_eLpR::computeThValue()
{
   // double sigma_500_bb_eLpR_madgraph = 0.72;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
        return (0.064*C_phiQ3+0.064*C_phiQ1+0.012*C_phib+0.1*C_bW+0.0132856*C_bB);
}



a_500_bb_eLpR::a_500_bb_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double a_500_bb_eLpR::computeThValue()
{
    //double a_500_bb_eLpR_madgraph = 67.7;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
        return (0.2*C_phiQ3+0.2*C_phiQ1-2.3*C_phib-8*C_bW-0.5204*C_bB);
}


sigma_500_bb_eRpL::sigma_500_bb_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigma_500_bb_eRpL::computeThValue()
{
   // double sigma_500_bb_eRpL_madgraph = 0.22;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
        return (0.02*C_phiQ3+0.02*C_phiQ1-0.024*C_phib+0.006*C_bW+0.00530063*C_bB);
}




a_500_bb_eRpL::a_500_bb_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double a_500_bb_eRpL::computeThValue()
{
    //double a_500_bb_eRpL_madgraph = 46.7;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
        return (-6.9*C_phiQ3-7.1*C_phiQ1-3.5*C_phib-1.5*C_bW+1.256*C_bB);
}



//500 tt observables


sigma_500_tt_eLpR::sigma_500_tt_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigma_500_tt_eLpR::computeThValue()
{
    //double sigma_500_tt_eLpR_madgraph = 930.754;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        return (-35.373*C_phit+54.9562*C_phiQ3-54.9562*C_phiQ1+639.01*C_tW+204.081*C_tB);
}



a_500_tt_eLpR::a_500_tt_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double a_500_tt_eLpR::computeThValue()
{
    //double a_500_tt_eLpR_madgraph = -0.380287;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        return (-0.0127129*C_phit-0.0326318*C_phiQ3+0.0326318*C_phiQ1-0.118262*C_tW-0.0384058*C_tB);
}


sigma_500_tt_eRpL::sigma_500_tt_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigma_500_tt_eRpL::computeThValue()
{
   // double sigma_500_tt_eRpL_madgraph = 481.105;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        return (32.3079*C_phit-16.7611*C_phiQ3+16.7611*C_phiQ1+31.5853*C_tW+267.885*C_tB);
}




a_500_tt_eRpL::a_500_tt_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double a_500_tt_eRpL::computeThValue()
{
   // double a_500_tt_eRpL_madgraph = -0.457824;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        return (-0.0413699*C_phit-0.0125082*C_phiQ3+0.0125082*C_phiQ1-0.0109327*C_tW-0.125174*C_tB);
}





pt_500_tt_eLpR::pt_500_tt_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double pt_500_tt_eLpR::computeThValue()
{
    //double pt_500_tt_eLpR_madgraph = (0.570477+0.573802+0.576393+0.570551+0.584227)/5;
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        return (-0.0196093*C_phit+0.215508*C_tW+0.0336945*C_tB);
}


pt_500_tt_eRpL::pt_500_tt_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double pt_500_tt_eRpL::computeThValue()
{
   // double pt_500_tt_eRpL_madgraph = (-0.432304-0.428132-0.423239-0.430406-0.431086)/5;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        return (-0.00550366*C_phit+0.0176743*C_phiQ3-0.0302046*C_phiQ1+0.104522*C_tW-0.204084*C_tB);
}

//double NPSMEFT6dtopquark::sigmattbarZ(const double sqrt_s) const
//{
//}


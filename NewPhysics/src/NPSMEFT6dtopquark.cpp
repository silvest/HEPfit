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
        = {"C_phit","C_phiQ3","C_phiQ1","C_tW","C_tB","C_tphi","C_phib","C_bW","C_bB","C_phitb","C_ed","C_eq","C_ld","C_lqP","C_eu","C_lu","C_lqM"};


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
    else if (name.compare("C_phitb") == 0)
        C_phitb = value;
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




C_phit::C_phit(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_phit::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
}

C_phiQ3::C_phiQ3(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_phiQ3::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
}


C_phiQ1::C_phiQ1(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_phiQ1::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
}

C_tW::C_tW(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_tW::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
}

C_tB::C_tB(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_tB::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
}


C_tphi::C_tphi(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_tphi::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
}

C_phib::C_phib(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_phib::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
}


C_bW::C_bW(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_bW::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
}


C_bB::C_bB(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_bB::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
}


C_phitb::C_phitb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_phitb::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
}


C_ed::C_ed(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_ed::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
}


C_eq::C_eq(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_eq::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
}

C_ld::C_ld(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_ld::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
}

C_lqP::C_lqP(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_lqP::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
}

C_eu::C_eu(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_eu::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
}

C_lu::C_lu(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_lu::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
}

C_lqM::C_lqM(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_lqM::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
}


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
        return smlep_bb + (0.023*C_phiQ3+0.023*C_phiQ1-0.005*C_phib+0.0018*C_bW*(-1*0.99*0.653228500107)*C_bW*(-1*0.99*0.653228500107)+0.002*C_bW*(-1*0.99*0.653228500107)*C_bB*(-1*0.99*0.349192873528))*(smlep_bb/lep_bb_madgraph);
    }
    else{
        return smlep_bb + (0.023*C_phiQ3+0.023*C_phiQ1-0.005*C_phib)*(smlep_bb/lep_bb_madgraph);
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
        return smlep_alr + (0.008*C_phiQ3+0.008*C_phiQ1+0.034*C_phib+0.0056*C_bW*(-1*0.99*0.653228500107)*C_bW*(-1*0.99*0.653228500107)-0.002*C_phib*C_bB*(-1*0.99*0.349192873528)+0.0015*C_bB*(-1*0.99*0.349192873528)*C_bB*(-1*0.99*0.349192873528)+0.0076*C_bW*(-1*0.99*0.653228500107)*C_bB*(-1*0.99*0.349192873528))*(smlep_alr/lep_alr_madgraph);
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
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_phitb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    if(flag_Quadratic){
        return ( smxttz + (0.041*C_phit+0.0024*C_phit*C_phit+0.0018*C_phit*C_phiQ3+0.0018*C_phit*C_phiQ1-0.066*C_phiQ3+0.0035*C_phiQ3*C_phiQ3
                +0.066*C_phiQ1+0.005*C_phiQ1*C_phiQ1-0.0001*C_phiQ3*C_tB+0.0001*C_phiQ1*C_tB+0.00068*C_tW+0.018*C_tW*C_tW+0.0001*C_tW*C_phiQ3
                -0.0002*C_tW*C_phiQ1+0.00024*C_tB+0.0016*C_tB*C_tB-0.01*C_tW*C_tB+0.002*C_bW*C_bW+0.0003*C_phitb+0.0003*C_phitb*C_phitb)*(smxttz/xttz_madgraph));
    }
    else{
        return ( smxttz + (0.041*C_phit-0.066*C_phiQ3+0.066*C_phiQ1+0.00068*C_tW+0.00024*C_tB+0.0003*C_phitb)*(smxttz/xttz_madgraph));
    }
}

sigmattA_1::sigmattA_1(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmattA_1::computeThValue()
{   
    double xtta_madgraph = 2.18;
    double smxtta = 0.063;
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    if(flag_Quadratic){
        return (smxtta + (0.0034*C_phiQ3+0.00035*C_phiQ3*C_phiQ3+0.015*C_tW+0.007*C_tW*C_tW
                +0.015*C_tB+0.007*C_tB*C_tB+0.014*C_tW*C_tB+0.001*C_bW*C_bW)*(smxtta/xtta_madgraph));
    }
    else{
        return (smxtta + (0.0034*C_phiQ3+0.015*C_tW+0.015*C_tB)*(smxtta/xtta_madgraph)); 
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
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    if(flag_Quadratic){
        return  smxtth + (0.0002*C_phiQ3+0.0007*C_tW+0.00014*C_tB+0.00067*C_tW*C_tW
                -0.049*C_tphi+0.0015*C_tphi*C_tphi+0.0001*C_bW*C_bW)*(smxtth/xtth_madgraph);
    }
    else{
        return  smxtth + (0.0002*C_phiQ3+0.0007*C_tW+0.00014*C_tB-0.049*C_tphi)*(smxtth/xtth_madgraph);
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
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    if(flag_Quadratic){
        return smxttw + (-0.00013*C_phit+0.00036*C_phiQ3-0.0003*C_phiQ1+0.00012*C_phiQ3*C_phiQ1
                +0.0027*C_tW+0.0032*C_tW*C_tW-0.00011*C_tW*C_phiQ1)*(smxttw/xttw_madgraph);
    }
    else{
        return smxttw + (-0.00013*C_phit+0.00036*C_phiQ3-0.0003*C_phiQ1+0.0027*C_tW)*(smxttw/xttw_madgraph);
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
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_phitb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    if(flag_Quadratic){
        return smxst + (5.26*C_phiQ3+0.16*C_phiQ3*C_phiQ3+1.52*C_tW+0.31*C_tW*C_tW+0.1*C_phiQ3*C_tW
                +0.19*C_bW*C_bW-0.012*C_bW*C_phitb+0.019*C_phitb*C_phitb)*(smxst/xst_madgraph);
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
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_phitb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    if(flag_Quadratic){
        return smxgbtw + (1.61*C_phiQ3+0.05*C_phiQ3*C_phiQ3-0.74*C_tW+0.135*C_tW*C_tW-0.046*C_phiQ3*C_tW
                +0.14*C_bW*C_bW-0.022*C_bW*C_phitb+0.017*C_phitb*C_phitb)*(smxgbtw/xgbtw_madgraph);
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
    double smxztq = 0.0942;
    double xztq_madgraph = 0.48;
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB =  myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_phitb =  myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    if(flag_Quadratic){
        return smxztq + (0.0029*C_phit+0.0005*C_phit*C_phiQ3-0.0005*C_phit*C_phiQ1+0.092*C_phiQ3
                +0.014*C_phiQ3*C_phiQ3+0.0003*C_phiQ3*C_bW*(-1*(0.99)*(0.6532))+0.001*C_phiQ3*C_tW+0.01*C_phiQ1
                +0.0003*C_phiQ1*C_tW+0.0017*C_phiQ3*C_phiQ1+0.007*C_tW+0.016*C_tW*C_tW-0.0003*C_phib
                +0.028*C_bW*(-1*(0.99)*(0.6532))*C_bW*(-1*(0.99)*(0.6532))+0.0024*C_bB*(-1*(0.99)*(0.3492))*C_bB*(-1*(0.99)*(0.3492))
                +0.002*C_phitb*C_phitb)*(smxztq/xztq_madgraph);
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
        return smF0 +(-0.04*C_tW)*(smF0/F0_madgraph);
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
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_Quadratic){
        return ((0.31*C_phiQ3+0.31*C_phiQ1+0.05*C_phib+0.27*C_bW*(-1*(0.99)*(0.653228500107))*C_bW*(-1*(0.99)*(0.653228500107))+0.077*C_bB*(-1*(0.99)*(0.349192873528))*C_bB*(-1*(0.99)*(0.3492))
                -0.25*C_bW*(-1*(0.99)*(0.653228500107))*C_bB*(-1*(0.99)*(0.349192873528))+0.091*C_ed-0.064*C_eq+0.71*C_ld+3.77*C_lqP));
    }
    else
    {
        return (0.31*C_phiQ3+0.31*C_phiQ1+0.05*C_phib
                +0.091*C_ed-0.064*C_eq+0.71*C_ld+3.77*C_lqP);
    }
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
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_Quadratic){
        return ((0.4*C_phiQ3+0.3*C_phiQ1-2.2*C_phib-5.1*C_bW*(-1*(0.99)*(0.653228500107))*C_bW*(-1*(0.99)*(0.653228500107))-1.29*C_bB*(-1*(0.99)*(0.349192873528))*C_bB*(-1*(0.99)*(0.349192873528))
                +4.46*C_bW*(-1*(0.99)*(0.653228500107))*C_bB*(-1*(0.99)*(0.349192873528))+3.8*C_eq-29.5*C_ld+8.57*C_lqP));
    }
    else{
        return (0.4*C_phiQ3+0.3*C_phiQ1-2.2*C_phib+3.8*C_eq-29.5*C_ld+8.57*C_lqP);
    }
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
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_Quadratic){
        return (0.094*C_phiQ3+0.094*C_phiQ1-0.11*C_phib+0.018*C_bW*(-1*(0.99)*(0.653228500107))*C_bW*(-1*(0.99)*(0.653228500107))
                +0.31*C_bB*(-1*(0.99)*(0.349192873528))*C_bB*(-1*(0.99)*(0.349192873528))+0.023*C_bW*(-1*(0.99)*(0.653228500107))*C_bB*(-1*(0.99)*(0.349192873528))+1.61*C_ed-1.08*C_eq+0.04*C_ld+0.23*C_lqP);
    }
    else{
        return (0.094*C_phiQ3+0.094*C_phiQ1-0.11*C_phib+1.61*C_ed-1.08*C_eq+0.04*C_ld+0.23*C_lqP);
    }
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
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_Quadratic){
        return ((-7.8*C_phiQ3-7.7*C_phiQ1-4.5*C_phib-0.23*C_bW*(-1*(0.99)*(0.653228500107))*C_bW*(-1*(0.99)*(0.653228500107))-7.75*C_bB*(-1*(0.99)*(0.349192873528))*C_bB*(-1*(0.99)*(0.349192873528))
                -0.61*C_bW*(-1*(0.99)*(0.653228500107))*C_bB*(-1*(0.99)*(0.349192873528))+62*C_ed+119*C_eq-4.6*C_ld+7.9*C_lqP));
    }
    else{
        return ((-7.8*C_phiQ3-7.7*C_phiQ1-4.5*C_phib+62*C_ed+119*C_eq-4.6*C_ld+7.9*C_lqP));
    }
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
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_Quadratic){
        return (0.064*C_phiQ3+0.064*C_phiQ1+0.012*C_phib+0.24*C_bW*(-1*(0.99)*(0.653228500107))*C_bW*(-1*(0.99)*(0.653228500107))+0.085*C_bB*(-1*(0.99)*(0.349192873528))*C_bB*(-1*(0.99)*(0.349192873528))
                -0.255*C_bW*(-1*(0.99)*(0.653228500107))*C_bB*(-1*(0.99)*(0.349192873528))+0.095*C_ed-0.05*C_eq+0.76*C_ld+3.5*C_lqP);
    }
    else{
        return (0.064*C_phiQ3+0.064*C_phiQ1+0.012*C_phib+0.095*C_ed-0.05*C_eq+0.76*C_ld+3.5*C_lqP);
    }
    
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
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_Quadratic){
        return (0.2*C_phiQ3+0.2*C_phiQ1-2.3*C_phib-13.7*C_bW*(-1*(0.99)*(0.653228500107))*C_bW*(-1*(0.99)*(0.653228500107))
                -3.8*C_bB*(-1*(0.99)*(0.349192873528))*C_bB*(-1*(0.99)*(0.349192873528))+13.5*C_bW*(-1*(0.99)*(0.653228500107))*C_bB*(-1*(0.99)*(0.349192873528))+8*C_eq-139*C_ld+40.3*C_lqP);
    }
    else{
        return (0.2*C_phiQ3+0.2*C_phiQ1-2.3*C_phib+8*C_eq-139*C_ld+40.3*C_lqP);
    }
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
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_Quadratic){
        return (0.02*C_phiQ3+0.02*C_phiQ1-0.024*C_phib+0.014*C_bW*(-1*(0.99)*(0.653228500107))*C_bW*(-1*(0.99)*(0.653228500107))+0.29*C_bB*(-1*(0.99)*(0.349192873528))*C_bB*(-1*(0.99)*(0.349192873528))
                -0.007*C_bW*(-1*(0.99)*(0.653228500107))*C_bB*(-1*(0.99)*(0.349192873528))+1.56*C_ed-0.84*C_eq+0.046*C_ld+0.2*C_lqP);
    }
    else{
        return (0.02*C_phiQ3+0.02*C_phiQ1-0.024*C_phib+1.56*C_ed-0.84*C_eq+0.046*C_ld+0.2*C_lqP);
    }
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
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_Quadratic){
        return (-6.9*C_phiQ3-7.1*C_phiQ1-3.5*C_phib-1.57*C_bW*(-1*(0.99)*(0.653228500107))*C_bW*(-1*(0.99)*(0.653228500107))-23*C_bB*(-1*(0.99)*(0.349192873528))*C_bB*(-1*(0.99)*(0.349192873528))
                +0.39*C_bW*(-1*(0.99)*(0.653228500107))*C_bB*(-1*(0.99)*(0.349192873528))+219*C_ed+380*C_eq-26.6*C_ld+28.4*C_lqP);
    }
    else{
        return (-6.9*C_phiQ3-7.1*C_phiQ1-3.5*C_phib+219*C_ed+380*C_eq-26.6*C_ld+28.4*C_lqP);
    }
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








//Prospects of Linear Colders at 1000 GeV
//1000 bb observables


sigma_1000_bb_eLpR::sigma_1000_bb_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigma_1000_bb_eLpR::computeThValue()
{
   // double sigma_1000_bb_eLpR_madgraph = 0.72;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_Quadratic){
        return (0.015*C_phiQ3+0.015*C_phiQ1+0.003*C_phib+0.23*C_bW*(-1*(0.99)*(0.6532))*C_bW*(-1*(0.99)*(0.6532))+0.087*C_bB*(-1*(0.99)*(0.3492))*C_bB*(-1*(0.99)*(0.3492))
                -0.255*C_bW*(-1*(0.99)*(0.6532))*C_bB*(-1*(0.99)*(0.3492))+0.093*C_ed-0.048*C_eq+0.77*C_ld+3.44*C_lqP);
    }
    else{
        return (0.015*C_phiQ3+0.015*C_phiQ1+0.003*C_phib+0.093*C_ed-0.048*C_eq+0.77*C_ld+3.44*C_lqP);
    }
    
}


a_1000_bb_eLpR::a_1000_bb_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double a_1000_bb_eLpR::computeThValue()
{
    //double a_1000_bb_eLpR_madgraph = 67.7;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_Quadratic){
        return (0.48*C_phiQ3+0.48*C_phiQ1-2.7*C_phib-23.6*C_bW*(-1*(0.99)*(0.6532))*C_bW*(-1*(0.99)*(0.6532))-5.76*C_bB*(-1*(0.99)*(0.3492))*C_bB*(-1*(0.99)*(0.3492))
                +21.2*C_bW*(-1*(0.99)*(0.6532))*C_bB*(-1*(0.99)*(0.3492))+37.4*C_eq-593*C_ld+145*C_lqP);
    }
    else{
        return (0.48*C_phiQ3+0.48*C_phiQ1-2.7*C_phib+37.4*C_eq-593*C_ld+145*C_lqP);
    }
}


sigma_1000_bb_eRpL::sigma_1000_bb_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigma_1000_bb_eRpL::computeThValue()
{
   // double sigma_1000_bb_eRpL_madgraph = 0.22;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_Quadratic){
        return (0.004*C_phiQ3+0.004*C_phiQ1-0.006*C_phib+0.014*C_bW*(-1*(0.99)*(0.6532))*C_bW*(-1*(0.99)*(0.6532))+0.29*C_bB*(-1*(0.99)*(0.3492))*C_bB*(-1*(0.99)*(0.3492))
                -0.013*C_bW*(-1*(0.99)*(0.6532))*C_bB*(-1*(0.99)*(0.3492))+1.55*C_ed-0.79*C_eq+0.046*C_ld+0.21*C_lqP);
    }
    else{
        return (0.004*C_phiQ3+0.004*C_phiQ1-0.006*C_phib+1.55*C_ed-0.79*C_eq+0.046*C_ld+0.21*C_lqP);
    }
}



a_1000_bb_eRpL::a_1000_bb_eRpL(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double a_1000_bb_eRpL::computeThValue()
{
    //double a_1000_bb_eRpL_madgraph = 46.7;
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_bB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bB();
    double C_ed = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ed();
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
    double C_ld = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_ld();
    double C_lqP = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqP();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    if(flag_Quadratic){
        return (-6.8*C_phiQ3-6.8*C_phiQ1-3.6*C_phib-3.15*C_bW*(-1*(0.99)*(0.6532))*C_bW*(-1*(0.99)*(0.6532))-27.5*C_bB*(-1*(0.99)*(0.3492))*C_bB*(-1*(0.99)*(0.3492))
                +0.44*C_bW*(-1*(0.99)*(0.6532))*C_bB*(-1*(0.99)*(0.3492))+826*C_ed+1611*C_eq-120*C_ld+6.5*C_lqP);
    }
    else{
        return (-6.8*C_phiQ3-6.8*C_phiQ1-3.6*C_phib+826*C_ed+1611*C_eq-120*C_ld+6.5*C_lqP);
    }
}



//OPTIMAL OBSERVABLES


op1::op1(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double op1::computeThValue()
{
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
        return C_phit;
}


op2::op2(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double op2::computeThValue()
{
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
        return (C_phiQ1-C_phiQ3);
}

op3::op3(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double op3::computeThValue()
{
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
        return C_tW;
}

op4::op4(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double op4::computeThValue()
{
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        return C_tB;
}


op_1000_1::op_1000_1(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double op_1000_1::computeThValue()
{
    double C_eu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eu();
        return C_eu;
}


op_1000_2::op_1000_2(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double op_1000_2::computeThValue()
{
    double C_eq = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_eq();
        return C_eq;
}


op_1000_3::op_1000_3(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double op_1000_3::computeThValue()
{
    double C_lu = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lu();
        return C_lu;
}


op_1000_4::op_1000_4(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double op_1000_4::computeThValue()
{
    double C_lqM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_lqM();
        return C_lqM;
}


op_1000_5::op_1000_5(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double op_1000_5::computeThValue()
{
    double C_phiQ1 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ1();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();

        return C_phiQ1-C_phiQ3;
}

op_1000_6::op_1000_6(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double op_1000_6::computeThValue()
{
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
        return C_phit;
}


op_1000_7::op_1000_7(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double op_1000_7::computeThValue()
{
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
        return C_tW;
}


op_1000_8::op_1000_8(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double op_1000_8::computeThValue()
{
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
        return C_tB;
}




//double NPSMEFT6dtopquark::sigmattbarZ(const double sqrt_s) const
//{
//}


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


const std::string NPSMEFT6dtopquark::NPSMEFT6dtopquarkVars[NNPSMEFT6dtopquarkVars]
        = {"C_phit","C_phiQ3","C_phiQ1","C_phiQm","C_tW","C_tZ","C_tB","C_tphi","C_phib","C_bW","C_bB","C_bZ" ,"C_phitb","C_tG","C_ed","C_eq","C_ld","C_lqP","C_eu","C_lu","C_lqM","C_tu8","C_td8", "C_Qq18","C_tq8","C_Qq38","C_Qu8","C_Qd8",
        "SM_ttZ_bin_0_40","SM_ttZ_bin_40_70","SM_ttZ_bin_70_110","SM_ttZ_bin_110_160","SM_ttZ_bin_160_220","SM_ttZ_bin_220_290","SM_ttZ_bin_290_400", "SM_ttA_bin_20_25", "SM_ttA_bin_25_30",
        "SM_ttA_bin_30_35", "SM_ttA_bin_35_40", "SM_ttA_bin_40_47", "SM_ttA_bin_47_55", "SM_ttA_bin_55_70", "SM_ttA_bin_70_85", "SM_ttA_bin_85_132", "SM_ttA_bin_132_180", "SM_ttA_bin_180_300",
        "SM_tAq_inc", "SM_tZQ_inc", "SM_ttA_inc", "SM_ttZ_inc", "SM_ttW_inc", "SM_tW_inc", "SM_tW_inc_8TeV", "SM_ttH_inc", "SM_sigmatchannel13", "SM_sigmatchannel8", "SM_sigmaschannel8", "SM_sigmaschannelTev",
        "SM_tH_tchan_value", "SM_ttbar_LHC13", "SM_ttbar_LHC8", "SM_ttbar_Tev", "F0_SM", "FL_SM","Rb_SM", "AFBLR_SM", "ttWqEM_SM"};


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
    else if (name.compare("SM_tAq_inc") == 0)
        SM_tAq_inc = value;
    else if (name.compare("SM_ttZ_bin_0_40") == 0)
        SM_ttZ_bin_0_40 = value;
    else if (name.compare("SM_ttZ_bin_40_70") == 0)
        SM_ttZ_bin_40_70 = value;
    else if (name.compare("SM_ttZ_bin_70_110") == 0)
        SM_ttZ_bin_70_110 = value;
    else if (name.compare("SM_ttZ_bin_110_160") == 0)
        SM_ttZ_bin_110_160 = value;
    else if (name.compare("SM_ttZ_bin_160_220") == 0)
        SM_ttZ_bin_160_220 = value;
    else if (name.compare("SM_ttZ_bin_220_290") == 0)
        SM_ttZ_bin_220_290 = value;
    else if (name.compare("SM_ttZ_bin_290_400") == 0)
        SM_ttZ_bin_290_400 = value;
    
    
    
    
    else if (name.compare("SM_ttA_bin_20_25") == 0)
        SM_ttA_bin_20_25 = value;
    else if (name.compare("SM_ttA_bin_25_30") == 0)
        SM_ttA_bin_25_30 = value;
    else if (name.compare("SM_ttA_bin_30_35") == 0)
        SM_ttA_bin_30_35 = value;
    else if (name.compare("SM_ttA_bin_35_40") == 0)
        SM_ttA_bin_35_40 = value;
    else if (name.compare("SM_ttA_bin_40_47") == 0)
        SM_ttA_bin_40_47 = value;
    else if (name.compare("SM_ttA_bin_47_55") == 0)
        SM_ttA_bin_47_55 = value;
    else if (name.compare("SM_ttA_bin_55_70") == 0)
        SM_ttA_bin_55_70 = value;
    else if (name.compare("SM_ttA_bin_70_85") == 0)
        SM_ttA_bin_70_85 = value;
    else if (name.compare("SM_ttA_bin_85_132") == 0)
        SM_ttA_bin_85_132 = value;
    else if (name.compare("SM_ttA_bin_132_180") == 0)
        SM_ttA_bin_132_180 = value;
    else if (name.compare("SM_ttA_bin_180_300") == 0)
        SM_ttA_bin_180_300 = value;
    
    
    
    
    else if (name.compare("SM_tZQ_inc") == 0)
        SM_tZQ_inc = value;
    else if (name.compare("SM_ttA_inc") == 0)
        SM_ttA_inc = value;
    else if (name.compare("SM_ttZ_inc") == 0)
        SM_ttZ_inc = value;
    else if (name.compare("SM_ttH_inc") == 0)
        SM_ttH_inc = value;
    else if (name.compare("SM_ttW_inc") == 0)
        SM_ttW_inc = value;
    else if (name.compare("SM_tW_inc") == 0)
        SM_tW_inc = value;
    else if (name.compare("SM_tW_inc_8TeV") == 0)
        SM_tW_inc_8TeV = value;
    else if (name.compare("SM_sigmatchannel13") == 0)
        SM_sigmatchannel13=value;
    else if (name.compare("SM_sigmatchannel8") == 0)
        SM_sigmatchannel8=value;
    else if (name.compare("SM_sigmaschannel8") == 0)
        SM_sigmaschannel8=value;
    
    else if (name.compare("SM_sigmaschannelTev") == 0)
        SM_sigmaschannelTev=value;
    
    else if(name.compare("SM_tH_tchan_value")==0)
        SM_tH_tchan_value=value;
    
    else if(name.compare("SM_ttbar_LHC13")==0)
        SM_ttbar_LHC13=value;
    else if(name.compare("SM_ttbar_LHC8")==0)
        SM_ttbar_LHC8=value;
    else if(name.compare("SM_ttbar_Tev")==0)
        SM_ttbar_Tev=value;
    else if(name.compare("F0_SM")==0)
        F0_SM=value;
    else if(name.compare("FL_SM")==0)
        FL_SM=value;
    else if(name.compare("Rb_SM")==0)
        Rb_SM=value;
    else if(name.compare("AFBLR_SM")==0)
        AFBLR_SM=value;
     else if(name.compare("ttWqEM_SM")==0)
        ttWqEM_SM=value;
    else
        NPbase::setParameter(name, value);
}




// Flag to swith on/off the quadratic terms. flag = true means quadratic terms
//are switched on
// Flag to swith between LO and NLO terms. flag = true means NLO terms are switched used

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


C_phiQm::C_phiQm(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_phiQm::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
}

C_tW::C_tW(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_tW::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
}

C_tZ::C_tZ(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_tZ::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
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

C_tG::C_tG(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_tG::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
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

C_bZ::C_bZ(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_bZ::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bZ();
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


C_tu8::C_tu8(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};


double C_tu8::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
}



C_td8::C_td8(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_td8::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
}

C_Qq18::C_Qq18(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_Qq18::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
}

C_tq8::C_tq8(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_tq8::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
}

C_Qq38::C_Qq38(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_Qq38::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
}

C_Qu8::C_Qu8(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_Qu8::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
}

C_Qd8::C_Qd8(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double C_Qd8::computeThValue()
{
    return myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Observables from LEP ////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

Rb_NPSMEFT6dtopquark::Rb_NPSMEFT6dtopquark(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double Rb_NPSMEFT6dtopquark::computeThValue()
{
    double Rb_SM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_Rb_SM();
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
{};

double AFBLR::computeThValue()
{
    double AFBLR_SM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_AFBLR_SM();
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
///////Observables from Tevatron ////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


sigmattbarTev::sigmattbarTev(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmattbarTev::computeThValue()
{
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double   xttbarTev_madgraph_LO = 4948.8;//fb
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double mC_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double mC_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    double SM_ttbar_tev = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttbar_Tev_value();

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




sigmaschannelTev::sigmaschannelTev(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmaschannelTev::computeThValue()
    {
    double SM_sigmaschannelTev = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_sigmaschannelTev();
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



/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Observables from LHC ////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////



/////// Helicities  /////////////////////////////////////////////////////////////////////////////////////


F0::F0(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double F0::computeThValue()
{
    double F0_madgraph = 0.70381;
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double F0_SM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_F0_SM();
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
{};

double FL::computeThValue()
{
    double FL_madgraph = 0.29619;
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double FL_SM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_FL_SM();
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
{};

double sigmattbarLHC13::computeThValue()
{
 
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
    double SM_ttbar_LHC13 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttbar_LHC13_value();


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
{};

double sigmattbarLHC8::computeThValue()
{
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

    double SM_ttbar_LHC8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttbar_LHC8_value();

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
{};

double sigmattZ::computeThValue()
{
    double xttzNLO_madgraph = 740;//fb
    double xttzLO_madgraph = 510.23;//fb
    double mC_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double mC_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
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
    double SM_ttZ_inc = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttZ_inc();
    
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic) {
                return  SM_ttZ_inc+ (-3.25*C_tZ+ 288.50*C_tG + 48.59*C_phit -76.50*C_phiQm 
                        + 62.33*C_tZ*C_tZ+ 291.80*C_tG*C_tG +3.20*C_phit*C_phit
                        +3.75*C_phiQm*C_phiQm-0.0016*C_phit*C_phiQm)*(SM_ttZ_inc/(xttzNLO_madgraph)) 
                        +(0.0096*mC_tu8+0.0060*mC_tu8*mC_tu8+0.0081*mC_td8+0.0039*mC_td8*mC_td8)*(1000*SM_ttZ_inc/xttzLO_madgraph);
            }
            else{
                return  SM_ttZ_inc + (-3.25*C_tZ+288.50*C_tG+48.59*C_phit
                        -76.50*C_phiQm)*(SM_ttZ_inc/(xttzNLO_madgraph))+
                        (0.0096*mC_tu8+0.0081*mC_td8)*(1000*SM_ttZ_inc/xttzLO_madgraph);
            }
    }
    else{
        if(flag_Quadratic){
                return  SM_ttZ_inc+ (-3.25*(-0.472123*C_tB + 0.881533*C_tW)+ 288.50*C_tG + 48.59*C_phit -76.50*(C_phiQ1-C_phiQ3) 
                        + 62.33*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW)+ 291.80*C_tG*C_tG +3.20*C_phit*C_phit
                        +3.75*(C_phiQ1-C_phiQ3)*(C_phiQ1-C_phiQ3)-0.0016*C_phit*(C_phiQ1-C_phiQ3))*(SM_ttZ_inc/(xttzNLO_madgraph)) 
                        +(0.0096*mC_tu8+0.0060*mC_tu8*mC_tu8+0.0081*mC_td8+0.0039*mC_td8*mC_td8)*(1000*SM_ttZ_inc/xttzLO_madgraph);
            }
            else{
                return  SM_ttZ_inc + (-3.25*(-0.472123*C_tB + 0.881533*C_tW)+288.50*C_tG+48.59*C_phit
                        -76.50*(C_phiQ1-C_phiQ3))*(SM_ttZ_inc/(xttzNLO_madgraph))+
                        (0.0096*mC_tu8+0.0081*mC_td8)*(1000*SM_ttZ_inc/xttzLO_madgraph);
            }
    }
}


sigmattA::sigmattA(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmattA::computeThValue()
{   
    double xtta_madgraph_NLO = 1986;//fb
    double xtta_madgraph_LO = 1317.6;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double mC_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double mC_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_tZ = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tZ();
    double C_tB = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tB();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double SM_ttA_inc = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttA_inc();
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttA_inc + (-88.4173*C_tZ + 110.514*C_tW + 578*C_tG  + 291.96*C_tZ*C_tZ + 380.115*C_tW*C_tW +152.0*C_tG*C_tG - 660.969*C_tZ*C_tW)*(SM_ttA_inc/xtta_madgraph_NLO)+
                        (0.0360*mC_tu8+0.0215*mC_td8+0.0147*mC_tu8*mC_tu8+0.0102*mC_td8*mC_td8)*(1000*SM_ttA_inc/xtta_madgraph_LO);
            }
            else{
                return  SM_ttA_inc + (-88.4173*C_tZ + 110.514*C_tW + 578*C_tG  )*(SM_ttA_inc/xtta_madgraph_NLO)+
                        (0.0360*mC_tu8+0.0215*mC_td8+0.0147*mC_tu8*mC_tu8)*(1000*SM_ttA_inc/xtta_madgraph_LO); 
            }
    }    
    
    else{
        if(flag_Quadratic){
            return  SM_ttA_inc + (-88.4173*(-0.472123*C_tB + 0.881533*C_tW) + 110.514*C_tW + 578*C_tG  + 291.96*(-0.472123*C_tB + 0.881533*C_tW)*(-0.472123*C_tB + 0.881533*C_tW) + 380.115*C_tW*C_tW +152.0*C_tG*C_tG - 660.969*(-0.472123*C_tB + 0.881533*C_tW)*C_tW)*(SM_ttA_inc/xtta_madgraph_NLO)+
                        (0.0360*mC_tu8+0.0215*mC_td8+0.0147*mC_tu8*mC_tu8+0.0102*mC_td8*mC_td8)*(1000*SM_ttA_inc/xtta_madgraph_LO);
        }
        else{
            return  SM_ttA_inc + (-88.4173*(-0.472123*C_tB + 0.881533*C_tW) + 110.514*C_tW + 578*C_tG  )*(SM_ttA_inc/xtta_madgraph_NLO)+
                        (0.0360*mC_tu8+0.0215*mC_td8)*(1000*SM_ttA_inc/xtta_madgraph_LO); 
        }
    }
        
}


sigmattW::sigmattW(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmattW::computeThValue()
{
    double SM_ttW_inc = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttW_inc();
    double NLOxttw_madgraph = 543.3;//fb
    double LOxttw_madgraph = 361.16;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double mC_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double mC_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();


    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  (SM_ttW_inc + (158.3*C_tG  + 23.4*C_tG*C_tG)*(SM_ttW_inc/NLOxttw_madgraph))+
                        (-0.00001*mC_tu8-0.00021*mC_td8-0.00013*mC_tu8*mC_tu8-0.00007*mC_td8*mC_td8)*(1000*SM_ttW_inc/LOxttw_madgraph);
            }
            else{
                return  (SM_ttW_inc + (158.3*C_tG)*(SM_ttW_inc/NLOxttw_madgraph))+
                        (-0.00001*mC_tu8-0.00021*mC_td8)*(1000*SM_ttW_inc/LOxttw_madgraph) ;
            }      
    }
    
    else{
            if(flag_Quadratic){
                return  (SM_ttW_inc + (158.3*C_tG  + 23.4*C_tG*C_tG)*(SM_ttW_inc/NLOxttw_madgraph))+
                        (-0.00001*mC_tu8-0.00021*mC_td8-0.00013*mC_tu8*mC_tu8-0.00007*mC_td8*mC_td8)*(1000*SM_ttW_inc/LOxttw_madgraph);
            }
            else{
                return  (SM_ttW_inc + (158.3*C_tG)*(SM_ttW_inc/NLOxttw_madgraph))+
                        (-0.00001*mC_tu8-0.00021*mC_td8)*(1000*SM_ttW_inc/LOxttw_madgraph) ;
            }      
    }
}


sigmatchannel13::sigmatchannel13(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double sigmatchannel13::computeThValue()
{
    double SM_sigmatchannel13 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_sigmatchannel13();
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
{};


double sigmatchannel8::computeThValue()
    {
    double SM_sigmatchannel8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_sigmatchannel8();
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
{};


double sigmaschannel8::computeThValue()
    {
    double SM_sigmaschannel8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_sigmaschannel8();
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
{};

double sigmatW::computeThValue()
{
    double NLOxgbtw_madgraph = 56188.5;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double SM_tW_inc = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_tW_inc();
    
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
{};

double sigmatW_8TeV::computeThValue()
{
    double NLOxgbtw_madgraph_8TeV = 17.3969;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_bW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_bW();
    double C_phitb = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phitb();
    double SM_tW_inc_8TeV = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_tW_inc_8TeV();

    
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
{};

double sigmatqZ::computeThValue()
{
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
    double SM_tZQ_inc = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_tZQ_inc();

    
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
{};

double sigmatAq::computeThValue()
{
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
    double SM_tAq_inc= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_tAq_inc();
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
{};

double tH_tchan::computeThValue()
{
    double C_tW = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tW();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double SM_tH_tchan_value= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_tH_tchan_value();
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
{};

double sigmattH::computeThValue()
{
    double xtth_madgraph_NLO = 459.7;//fb
    double xtth_madgraph_LO = 354.57;//fb
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double C_tG = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tG();
    double C_tphi = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tphi();
    double mC_tu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tu8();
    double mC_td8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_td8();
    double SM_ttH_inc = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttH_inc();
    
    if(flag_LHC_WG_Basis){
            if(flag_Quadratic){
                return  SM_ttH_inc + (-59.3*C_tphi+ 458.9*C_tG +1.9*C_tphi*C_tphi+ 598.6*C_tG*C_tG)*(SM_ttH_inc/xtth_madgraph_NLO)+
                        (0.0203*mC_tu8+0.0126*mC_td8+0.0088*mC_tu8*mC_tu8+0.0046*mC_td8*mC_td8)*(1000*SM_ttH_inc/xtth_madgraph_LO) ;
            }
            else{
                return  SM_ttH_inc + (-59.3*C_tphi+ 458.9*C_tG )*(SM_ttH_inc/xtth_madgraph_NLO)+
                        (0.0203*mC_tu8+0.0126*mC_td8)*(1000*SM_ttH_inc/xtth_madgraph_LO);
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
{};

double ttHSUM::computeThValue()
{
    
    
    double xtth_madgraph_NLO = 459.7;//fb
    double xtth_madgraph_LO = 354.57;//fb
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
    double SM_ttH_inc = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttH_inc();
    double C_Qq18 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq18();
    double C_tq8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_tq8();
    double C_Qq38 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qq38();
    double C_Qu8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qu8();
    double C_Qd8 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_Qd8();
    


    double SM_tH_tchan_value= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_tH_tchan_value();
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





ttWqEM::ttWqEM(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double ttWqEM::computeThValue()
{
    double ttWqEM_SM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_ttWqEM_SM();
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
{};

double ttWqSUM::computeThValue()
{
    double ttWqEM_SM = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_ttWqEM_SM();

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
    
    double SM_ttW_inc = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttW_inc();
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


ttZ_bin_0_40::ttZ_bin_0_40(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};



double ttZ_bin_0_40::computeThValue()
{
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
    double SM_ttZ_bin_0_40= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttZ_bin_0_40();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
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
}


ttZ_bin_40_70::ttZ_bin_40_70(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};



double ttZ_bin_40_70::computeThValue()
{
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
    double SM_ttZ_bin_40_70= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttZ_bin_40_70();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
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
}

ttZ_bin_70_110::ttZ_bin_70_110(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};


double ttZ_bin_70_110::computeThValue()
{
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
    double SM_ttZ_bin_70_110= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttZ_bin_70_110();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
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
}

ttZ_bin_110_160::ttZ_bin_110_160(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};


double ttZ_bin_110_160::computeThValue()
{
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
    double SM_ttZ_bin_110_160= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttZ_bin_110_160();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();  
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
}


ttZ_bin_160_220::ttZ_bin_160_220(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double ttZ_bin_160_220::computeThValue()
{
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
    double SM_ttZ_bin_160_220= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttZ_bin_160_220();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
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
}

ttZ_bin_220_290::ttZ_bin_220_290(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};


double ttZ_bin_220_290::computeThValue()
{

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
    double SM_ttZ_bin_220_290= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttZ_bin_220_290();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic(); 
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
}

ttZ_bin_290_400::ttZ_bin_290_400(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};


double ttZ_bin_290_400::computeThValue()
{
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
    double SM_ttZ_bin_290_400= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttZ_bin_290_400();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////



/////// ttA differential cross section for different bins ///////////////////////////////////////////////


ttA_bin_20_25::ttA_bin_20_25(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double ttA_bin_20_25::computeThValue()
{
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
    double SM_ttA_bin_20_25= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttA_bin_20_25();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
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
}


ttA_bin_25_30::ttA_bin_25_30(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double ttA_bin_25_30::computeThValue()
{
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
    double SM_ttA_bin_25_30= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttA_bin_25_30();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
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
}


ttA_bin_30_35::ttA_bin_30_35(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double ttA_bin_30_35::computeThValue()
{
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
    double SM_ttA_bin_30_35= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttA_bin_30_35();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double ttA_bin_30_35_madgraph_NLO=0.898765;//fb
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
}







ttA_bin_35_40::ttA_bin_35_40(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double ttA_bin_35_40::computeThValue()
{
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
    double SM_ttA_bin_35_40= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttA_bin_35_40();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
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
}







ttA_bin_40_47::ttA_bin_40_47(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double ttA_bin_40_47::computeThValue()
{
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
    double SM_ttA_bin_40_47= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttA_bin_40_47();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
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
}







ttA_bin_47_55::ttA_bin_47_55(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double ttA_bin_47_55::computeThValue()
{
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
    double SM_ttA_bin_47_55= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttA_bin_47_55();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
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
}







ttA_bin_55_70::ttA_bin_55_70(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double ttA_bin_55_70::computeThValue()
{
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
    double SM_ttA_bin_55_70= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttA_bin_55_70();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
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
}







ttA_bin_70_85::ttA_bin_70_85(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double ttA_bin_70_85::computeThValue()
{
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
    double SM_ttA_bin_70_85= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttA_bin_70_85();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double ttA_bin_70_85_madgraph_NLO=0.22535;//fb
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
}




ttA_bin_85_132::ttA_bin_85_132(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double ttA_bin_85_132::computeThValue()
{
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
    double SM_ttA_bin_85_132= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttA_bin_85_132();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
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
}





ttA_bin_132_180::ttA_bin_132_180(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double ttA_bin_132_180::computeThValue()
{
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
    double SM_ttA_bin_132_180= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttA_bin_132_180();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double ttA_bin_132_180_madgraph_NLO=0.0582305;//fb
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
}






ttA_bin_180_300::ttA_bin_180_300(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double ttA_bin_180_300::computeThValue()
{
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
    double SM_ttA_bin_180_300= myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_SM_ttA_bin_180_300();
    bool   flag_LHC_WG_Basis=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_LHC_WG_Basis();
    bool   flag_Quadratic=myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_flag_Quadratic();
    double ttA_bin_180_300_madgraph_NLO=0.0256399;//fb
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
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// OBSERVABLES FOR PROSPECTS OF FUTURE COLLIDERS ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////// Prospects of Linear Colders at 250 GeV //////////////////////////////////////////////////////////


sigma_250_bb_eLpR::sigma_250_bb_eLpR(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

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
{};

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
{};

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
{};

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
{};

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
{};

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
{};

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
{};

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
{};

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
{};

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
{};

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
{};

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
{};

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
{};

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
{};

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
{};

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
{};

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
{};

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
{};

double gLt::computeThValue()
{
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double vev = 246.22;
   
    return -C_phiQm*(vev*vev/2.0);
}

gLb::gLb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double gLb::computeThValue()
{
    double C_phiQm = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQm();
    double C_phiQ3 = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phiQ3();
    double vev = 246.22;
   
    return -(C_phiQm + 2.0*C_phiQ3)*(vev*vev/2.0);
}

gRt::gRt(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double gRt::computeThValue()
{
    double C_phit = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phit();
    double vev = 246.22;
   
    return -C_phit*(vev*vev/2.0);
}

gRb::gRb(const StandardModel& SM_i)
: ThObservable(SM_i),myNPSMEFT6dtopquark(static_cast<const NPSMEFT6dtopquark&> (SM_i))
{};

double gRb::computeThValue()
{
    double C_phib = myNPSMEFT6dtopquark.getNPSMEFT6dtopquark_C_phib();
    double vev = 246.22;
   
    return -C_phib*(vev*vev/2.0);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


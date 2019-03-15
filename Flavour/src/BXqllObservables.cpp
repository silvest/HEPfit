/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BXqllObservables.h"
#include "StandardModel.h"
#include <TMath.h>

/*******************************************************************************
 * Observables                                                                 *
 * ****************************************************************************/


R_BXqll::R_BXqll(const StandardModel& SM_i, QCD::quark quark_i, QCD::lepton lep_i) 
: ThObservable(SM_i), myBXqll(SM_i, quark_i, lep_i) 
{  
//    if (SM.getModelName().compare("StandardModel") != 0) std::cout << "\nWARNING: B to Xq l+ l-: R_BXqll not implemented in: " + SM.getModelName() + " model, returning Standard Model value.\n" << std::endl;
//    lep = lep_i;
//    quark = quark_i;
//    
    setParametersForObservable(myBXqll.initializeBXqllParameters());
}

double R_BXqll::computeThValue()
{ 
//    double q_min = getBinMin();
//    double q_max = getBinMax();
//    return (myBXqll.integrate_Rquark(q_min, q_max, LOWQ2));
    return myBXqll.getR_LOWQ2(0.15);
}


Rlow_BXqll::Rlow_BXqll(const StandardModel& SM_i, QCD::quark quark_i, QCD::lepton lep_i) 
: ThObservable(SM_i), myBXqll(SM_i, quark_i, lep_i) 
{  
//    if (SM.getModelName().compare("StandardModel") != 0) std::cout << "\nWARNING: B to Xq l+ l-: Rlow_BXqll not implemented in: " + SM.getModelName() + " model, returning Standard Model value.\n" << std::endl;
//    lep = lep_i;
//    quark = quark_i;
//    
    setParametersForObservable(myBXqll.initializeBXqllParameters());
}

double Rlow_BXqll::computeThValue()
{
    double q_min = getBinMin();
    double q_max = getBinMax();
    return (myBXqll.integrate_Rquark(q_min, q_max, LOWQ2));
}


Rhigh_BXqll::Rhigh_BXqll(const StandardModel& SM_i, QCD::quark quark_i, QCD::lepton lep_i) 
: ThObservable(SM_i)
//, myBXqll(SM_i, quark_i, lep_i) 
{  
//    if (SM.getModelName().compare("StandardModel") != 0) std::cout << "\nWARNING: B to Xq l+ l-: Rhigh_BXqll not implemented in: " + SM.getModelName() + " model, returning Standard Model value.\n" << std::endl;
//    lep = lep_i;
//    quark = quark_i;
//    
//    setParametersForObservable(myBXqll.initializeBXqllParameters());
}

double Rhigh_BXqll::computeThValue()
{   
//    double sh_min = 14.4/4.9/4.9; // Isidori range
//    double sh_max = 1.;
//    return (myBXqll.integrate_Rquark(sh_min, sh_max, HIGHQ2));
    return 0;
}

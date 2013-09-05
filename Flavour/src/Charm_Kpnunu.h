/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CHARM_KPNUNU_H
#define	CHARM_KPNUNU_H

#include "StandardModel.h"
#include <sstream>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_clausen.h>

using namespace gslpp;

class Charm_Kpnunu {

public:
    
    Charm_Kpnunu(const StandardModel& model_i);
    
    ~Charm_Kpnunu();
    
    vector<double> Cp(orders order);
    
    matrix<double> RGevolP(int nf, orders order);

    vector<double> MatchingCp(orders order);    
    
    vector<double> C_p(orders order);
    
    double C_P(orders order);
    
    vector<double> Cb(orders order);
    
    matrix<double> RGevolB(int nf, orders order);
    
    vector<double> MatchingCb(orders order);
    
    vector<double> C_b(orders order);
    
    double C_Be(orders order);
    
    double C_Bt(orders order);
    
    double P_C(orders order);
    
    double C_TOT(orders order, orders_ew order);
    
private:
    const StandardModel& model;
    const StandardModelMatching& modelmatching;
    vector<double> cp, dcp, c_p, cpmuW0, cpmuW1, cpmuW2, cb, dcb, c_b, cbmuW0,
                   cbmuW1, cbmuW2;
    matrix<double> U4p, U5p, J5p1, J4p1, J5p2, J4p2, dc_p, 
                   U4b, U5b, J5b1, J4b1, J5b2, J4b2, dc_b;
    double etab, etacb, etac, mc, kc ,xi1c, xi2c, xc, CP, CBe, CBt;
};

#endif	/* CHARM_KPNUNU_H */


/* 
 * File:   EWSM.cpp
 * Author: mishima
 */

#include "EWSM.h"
#include <cmath>


EWSM::EWSM(const StandardModel& SM_i) : SM(SM_i) {
    EWSMC = new EWSMcommon(SM);
    myOneLoopEW = new OneLoopEW(*EWSMC);
    myTwoLoopQCD = new TwoLoopQCD(*EWSMC);
    myThreeLoopQCD = new ThreeLoopQCD(*EWSMC);
    myTwoLoopEW = new TwoLoopEW(*EWSMC);
    myThreeLoopEW2QCD = new ThreeLoopEW2QCD(*EWSMC);
    myThreeLoopEW = new ThreeLoopEW(*EWSMC);
}

EWSM::EWSM(const EWSM& orig) : SM(orig.SM) {
    EWSMC = new EWSMcommon(SM);
    myOneLoopEW = new OneLoopEW(*EWSMC);
    myTwoLoopQCD = new TwoLoopQCD(*EWSMC);
    myThreeLoopQCD = new ThreeLoopQCD(*EWSMC);
    myTwoLoopEW = new TwoLoopEW(*EWSMC);
    myThreeLoopEW2QCD = new ThreeLoopEW2QCD(*EWSMC);
    myThreeLoopEW = new ThreeLoopEW(*EWSMC);
    
    for (int j=0; j<orders_EW_size; j++) {
        EWSM::orders_EW j_order = (EWSM::orders_EW) j;
        DeltaAlpha_l[j] = orig.getDeltaAlpha_l(j_order);
        DeltaAlpha_t[j] = orig.getDeltaAlpha_t(j_order);
        DeltaRho[j] = orig.getDeltaRho(j_order);
        DeltaR_rem[j] = orig.getDeltaR_rem(j_order);
        for (int i=0; i<6; i++) {
            StandardModel::lepton i_l = (StandardModel::lepton) i;
            StandardModel::quark i_q = (StandardModel::quark) i;
            deltaRho_rem_l[i][j] = orig.getDeltaRho_rem_l(i_l, j_order);
            deltaRho_rem_q[i][j] = orig.getDeltaRho_rem_q(i_q, j_order);
            deltaKappa_rem_l[i][j] = orig.getDeltaKappa_rem_l(i_l, j_order);
            deltaKappa_rem_q[i][j] = orig.getDeltaKappa_rem_q(i_q, j_order);
        }
    }
    DeltaRbar_rem = orig.getDeltaRbar_rem();
    rho_GammaW_leptons = orig.getRho_GammaW_leptons();
    rho_GammaW_quarks = orig.getRho_GammaW_quarks();
}

EWSM::~EWSM() {
    if (myOneLoopEW != NULL) delete myOneLoopEW;
    if (myTwoLoopQCD != NULL) delete myTwoLoopQCD;
    if (myThreeLoopQCD != NULL) delete myThreeLoopQCD;
    if (myTwoLoopEW != NULL) delete myTwoLoopEW;
    if (myThreeLoopEW2QCD != NULL) delete myThreeLoopEW2QCD;
    if (myThreeLoopEW != NULL) delete myThreeLoopEW;
    if (EWSMC != NULL) delete EWSMC;
}


////////////////////////////////////////////////////////////////////////

void EWSM::ComputeDeltaAlpha(const bool flag_order[orders_EW_size]) {
    
    if (flag_order[EW1]) 
        DeltaAlpha_l[EW1] = myOneLoopEW->DeltaAlpha_l();
    if (flag_order[EW1QCD1]) 
        DeltaAlpha_l[EW1QCD1] = myTwoLoopQCD->DeltaAlpha_l();
    if (flag_order[EW1QCD2]) 
        DeltaAlpha_l[EW1QCD2] = myThreeLoopQCD->DeltaAlpha_l();
    if (flag_order[EW2]) 
        DeltaAlpha_l[EW2] = myTwoLoopEW->DeltaAlpha_l();
    if (flag_order[EW2QCD1]) 
        DeltaAlpha_l[EW2QCD1] = myThreeLoopEW2QCD->DeltaAlpha_l();
    if (flag_order[EW3]) 
        DeltaAlpha_l[EW3] = myThreeLoopEW->DeltaAlpha_l();
    
    if (flag_order[EW1]) 
        DeltaAlpha_t[EW1] = myOneLoopEW->DeltaAlpha_t();
    if (flag_order[EW1QCD1]) 
        DeltaAlpha_t[EW1QCD1] = myTwoLoopQCD->DeltaAlpha_t();
    if (flag_order[EW1QCD2]) 
        DeltaAlpha_t[EW1QCD2] = myThreeLoopQCD->DeltaAlpha_t();
    if (flag_order[EW2]) 
        DeltaAlpha_t[EW2] = myTwoLoopEW->DeltaAlpha_t();
    if (flag_order[EW2QCD1]) 
        DeltaAlpha_t[EW2QCD1] = myThreeLoopEW2QCD->DeltaAlpha_t();
    if (flag_order[EW3]) 
        DeltaAlpha_t[EW3] = myThreeLoopEW->DeltaAlpha_t();
}
 
void EWSM::ComputeCC(double Mw_i, const bool flag_order[orders_EW_size]) {    
    EWSMC->ComputeForCC(Mw_i);
        
    if (flag_order[EW1]) 
        DeltaRho[EW1] = myOneLoopEW->DeltaRho();
    if (flag_order[EW1QCD1]) 
        DeltaRho[EW1QCD1] = myTwoLoopQCD->DeltaRho();
    if (flag_order[EW1QCD2]) 
        DeltaRho[EW1QCD2] = myThreeLoopQCD->DeltaRho();
    if (flag_order[EW2]) 
        DeltaRho[EW2] = myTwoLoopEW->DeltaRho();
    if (flag_order[EW2QCD1]) 
        DeltaRho[EW2QCD1] = myThreeLoopEW2QCD->DeltaRho();
    if (flag_order[EW3]) 
        DeltaRho[EW3] = myThreeLoopEW->DeltaRho();
    
    if (flag_order[EW1]) 
        DeltaR_rem[EW1] = myOneLoopEW->DeltaR_rem();
    if (flag_order[EW1QCD1]) 
        DeltaR_rem[EW1QCD1] = myTwoLoopQCD->DeltaR_rem();
    if (flag_order[EW1QCD2]) 
        DeltaR_rem[EW1QCD2] = myThreeLoopQCD->DeltaR_rem();
    if (flag_order[EW2]) 
        DeltaR_rem[EW2] = myTwoLoopEW->DeltaR_rem();
    if (flag_order[EW2QCD1]) 
        DeltaR_rem[EW2QCD1] = myThreeLoopEW2QCD->DeltaR_rem();
    if (flag_order[EW3]) 
        DeltaR_rem[EW3] = myThreeLoopEW->DeltaR_rem();    
}    

void EWSM::ComputeNC(double Mw_i, const bool flag_order[orders_EW_size]) {
    EWSMC->ComputeForNC(Mw_i);
    
    DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem();    

    for (int i=0; i<6; i++) {
        StandardModel::lepton i_l = (StandardModel::lepton) i;
        StandardModel::quark i_q = (StandardModel::quark) i;
        
        if (flag_order[EW1]) 
            deltaRho_rem_l[i][EW1] = myOneLoopEW->deltaRho_rem_l(i_l);
        if (flag_order[EW1QCD1]) 
            deltaRho_rem_l[i][EW1QCD1] = myTwoLoopQCD->deltaRho_rem_l(i_l);
        if (flag_order[EW1QCD2]) 
            deltaRho_rem_l[i][EW1QCD2] = myThreeLoopQCD->deltaRho_rem_l(i_l);
        if (flag_order[EW2]) 
            deltaRho_rem_l[i][EW2] = myTwoLoopEW->deltaRho_rem_l(i_l);
        if (flag_order[EW2QCD1]) 
            deltaRho_rem_l[i][EW2QCD1] = myThreeLoopEW2QCD->deltaRho_rem_l(i_l);
        if (flag_order[EW3]) 
            deltaRho_rem_l[i][EW3] = myThreeLoopEW->deltaRho_rem_l(i_l);
        
        if (flag_order[EW1]) 
            deltaRho_rem_q[i][EW1] = myOneLoopEW->deltaRho_rem_q(i_q);
        if (flag_order[EW1QCD1]) 
            deltaRho_rem_q[i][EW1QCD1] = myTwoLoopQCD->deltaRho_rem_q(i_q);
        if (flag_order[EW1QCD2]) 
            deltaRho_rem_q[i][EW1QCD2] = myThreeLoopQCD->deltaRho_rem_q(i_q);
        if (flag_order[EW2]) 
            deltaRho_rem_q[i][EW2] = myTwoLoopEW->deltaRho_rem_q(i_q);
        if (flag_order[EW2QCD1]) 
            deltaRho_rem_q[i][EW2QCD1] = myThreeLoopEW2QCD->deltaRho_rem_q(i_q);
        if (flag_order[EW3]) 
            deltaRho_rem_q[i][EW3] = myThreeLoopEW->deltaRho_rem_q(i_q);
        
        if (flag_order[EW1]) 
            deltaKappa_rem_l[i][EW1] = myOneLoopEW->deltaKappa_rem_l(i_l);
        if (flag_order[EW1QCD1]) 
            deltaKappa_rem_l[i][EW1QCD1] = myTwoLoopQCD->deltaKappa_rem_l(i_l);
        if (flag_order[EW1QCD2]) 
            deltaKappa_rem_l[i][EW1QCD2] = myThreeLoopQCD->deltaKappa_rem_l(i_l);
        if (flag_order[EW2]) 
            deltaKappa_rem_l[i][EW2] = myTwoLoopEW->deltaKappa_rem_l(i_l);
        if (flag_order[EW2QCD1]) 
            deltaKappa_rem_l[i][EW2QCD1] = myThreeLoopEW2QCD->deltaKappa_rem_l(i_l);
        if (flag_order[EW3]) 
            deltaKappa_rem_l[i][EW3] = myThreeLoopEW->deltaKappa_rem_l(i_l);
        
        if (flag_order[EW1]) 
            deltaKappa_rem_q[i][EW1] = myOneLoopEW->deltaKappa_rem_q(i_q);
        if (flag_order[EW1QCD1]) 
            deltaKappa_rem_q[i][EW1QCD1] = myTwoLoopQCD->deltaKappa_rem_q(i_q);
        if (flag_order[EW1QCD2]) 
            deltaKappa_rem_q[i][EW1QCD2] = myThreeLoopQCD->deltaKappa_rem_q(i_q);
        if (flag_order[EW2]) 
            deltaKappa_rem_q[i][EW2] = myTwoLoopEW->deltaKappa_rem_q(i_q);
        if (flag_order[EW2QCD1]) 
            deltaKappa_rem_q[i][EW2QCD1] = myThreeLoopEW2QCD->deltaKappa_rem_q(i_q);
        if (flag_order[EW3]) 
            deltaKappa_rem_q[i][EW3] = myThreeLoopEW->deltaKappa_rem_q(i_q);
    }
}

void EWSM::ComputeRhoWij(double Mw_i, const bool flag_order[orders_EW_size]) {
    EWSMC->ComputeForRhoWij(Mw_i);

    if (flag_order[EW1]) {
        /* independent of fermion masses */
        rho_GammaW_leptons = myOneLoopEW->rho_GammaW_l(StandardModel::NEUTRINO_1,
                                                       StandardModel::ELECTRON);
        rho_GammaW_quarks = myOneLoopEW->rho_GammaW_q(StandardModel::UP, 
                                                      StandardModel::DOWN);
    }
}
    








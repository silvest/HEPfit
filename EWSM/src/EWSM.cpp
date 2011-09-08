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

    /* Initializations */
    for (int j=0; j<=orders_EW_size; j++) {
        DeltaAlpha_l[j] = 0.0;
        DeltaAlpha_t[j] = 0.0;
        DeltaRho[j] = 0.0;
        DeltaR_rem[j] = 0.0;
        for (int i=0; i<6; i++) {
            deltaRho_rem_l[i][j] = 0.0;
            deltaRho_rem_q[i][j] = 0.0;
            deltaKappa_rem_l[i][j] = 0.0;
            deltaKappa_rem_q[i][j] = 0.0;
        }
    }
    DeltaRbar_rem = 0.0;
}

EWSM::EWSM(const EWSM& orig) : SM(orig.SM) {
    EWSMC = new EWSMcommon(SM);
    myOneLoopEW = new OneLoopEW(*EWSMC);
    myTwoLoopQCD = new TwoLoopQCD(*EWSMC);
    myThreeLoopQCD = new ThreeLoopQCD(*EWSMC);
    myTwoLoopEW = new TwoLoopEW(*EWSMC);
    myThreeLoopEW2QCD = new ThreeLoopEW2QCD(*EWSMC);
    myThreeLoopEW = new ThreeLoopEW(*EWSMC);
    
    for (int j=0; j<=orders_EW_size; j++) {
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
    if (flag_order[EW1]==true) 
        DeltaAlpha_l[EW1] = myOneLoopEW->DeltaAlpha_l();
    if (flag_order[EW1QCD1]==true) 
        DeltaAlpha_l[EW1QCD1] = myTwoLoopQCD->DeltaAlpha_l();
    if (flag_order[EW1QCD2]==true) 
        DeltaAlpha_l[EW1QCD2] = myThreeLoopQCD->DeltaAlpha_l();
    if (flag_order[EW2]==true) 
        DeltaAlpha_l[EW2] = myTwoLoopEW->DeltaAlpha_l();
    if (flag_order[EW2QCD1]==true) 
        DeltaAlpha_l[EW2QCD1] = myThreeLoopEW2QCD->DeltaAlpha_l();
    if (flag_order[EW3]==true) 
        DeltaAlpha_l[EW3] = myThreeLoopEW->DeltaAlpha_l();
    
    if (flag_order[EW1]==true) 
        DeltaAlpha_t[EW1] = myOneLoopEW->DeltaAlpha_t();
    if (flag_order[EW1QCD1]==true) 
        DeltaAlpha_t[EW1QCD1] = myTwoLoopQCD->DeltaAlpha_t();
    if (flag_order[EW1QCD2]==true) 
        DeltaAlpha_t[EW1QCD2] = myThreeLoopQCD->DeltaAlpha_t();
    if (flag_order[EW2]==true) 
        DeltaAlpha_t[EW2] = myTwoLoopEW->DeltaAlpha_t();
    if (flag_order[EW2QCD1]==true) 
        DeltaAlpha_t[EW2QCD1] = myThreeLoopEW2QCD->DeltaAlpha_t();
    if (flag_order[EW3]==true) 
        DeltaAlpha_t[EW3] = myThreeLoopEW->DeltaAlpha_t();
    
    /* The sum of the corrections */
    for (int j=0; j<orders_EW_size; j++) {
        DeltaAlpha_l[orders_EW_size] += DeltaAlpha_l[j];
        DeltaAlpha_t[orders_EW_size] += DeltaAlpha_t[j];
    }
}
 
void EWSM::ComputeCC(double Mw_i, const bool flag_order[orders_EW_size]) {    
    EWSMC->Compute(Mw_i);
        
    if (flag_order[EW1]==true) 
        DeltaRho[EW1] = myOneLoopEW->DeltaRho();
    if (flag_order[EW1QCD1]==true) 
        DeltaRho[EW1QCD1] = myTwoLoopQCD->DeltaRho();
    if (flag_order[EW1QCD2]==true) 
        DeltaRho[EW1QCD2] = myThreeLoopQCD->DeltaRho();
    if (flag_order[EW2]==true) 
        DeltaRho[EW2] = myTwoLoopEW->DeltaRho();
    if (flag_order[EW2QCD1]==true) 
        DeltaRho[EW2QCD1] = myThreeLoopEW2QCD->DeltaRho();
    if (flag_order[EW3]==true) 
        DeltaRho[EW3] = myThreeLoopEW->DeltaRho();
    
    if (flag_order[EW1]==true) 
        DeltaR_rem[EW1] = myOneLoopEW->DeltaR_rem();
    if (flag_order[EW1QCD1]==true) 
        DeltaR_rem[EW1QCD1] = myTwoLoopQCD->DeltaR_rem();
    if (flag_order[EW1QCD2]==true) 
        DeltaR_rem[EW1QCD2] = myThreeLoopQCD->DeltaR_rem();
    if (flag_order[EW2]==true) 
        DeltaR_rem[EW2] = myTwoLoopEW->DeltaR_rem();
    if (flag_order[EW2QCD1]==true) 
        DeltaR_rem[EW2QCD1] = myThreeLoopEW2QCD->DeltaR_rem();
    if (flag_order[EW3]==true) 
        DeltaR_rem[EW3] = myThreeLoopEW->DeltaR_rem();    
    
    /* The sum of the corrections */
    for (int j=0; j<orders_EW_size; j++) {
        DeltaRho[orders_EW_size] += DeltaRho[j];
        DeltaR_rem[orders_EW_size] += DeltaR_rem[j];
    }
    
}    

void EWSM::ComputeNC(double Mw_i, const bool flag_order[orders_EW_size]) {
    EWSMC->Compute(Mw_i);
    
    DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem();    

    for (int i=0; i<6; i++) {
        StandardModel::lepton i_l = (StandardModel::lepton) i;
        StandardModel::quark i_q = (StandardModel::quark) i;
        
        if (flag_order[EW1]==true) 
            deltaRho_rem_l[i][EW1] = myOneLoopEW->deltaRho_rem_l(i_l);
        if (flag_order[EW1QCD1]==true) 
            deltaRho_rem_l[i][EW1QCD1] = myTwoLoopQCD->deltaRho_rem_l(i_l);
        if (flag_order[EW1QCD2]==true) 
            deltaRho_rem_l[i][EW1QCD2] = myThreeLoopQCD->deltaRho_rem_l(i_l);
        if (flag_order[EW2]==true) 
            deltaRho_rem_l[i][EW2] = myTwoLoopEW->deltaRho_rem_l(i_l);
        if (flag_order[EW2QCD1]==true) 
            deltaRho_rem_l[i][EW2QCD1] = myThreeLoopEW2QCD->deltaRho_rem_l(i_l);
        if (flag_order[EW3]==true) 
            deltaRho_rem_l[i][EW3] = myThreeLoopEW->deltaRho_rem_l(i_l);
        
        if (flag_order[EW1]==true) 
            deltaRho_rem_q[i][EW1] = myOneLoopEW->deltaRho_rem_q(i_q);
        if (flag_order[EW1QCD1]==true) 
            deltaRho_rem_q[i][EW1QCD1] = myTwoLoopQCD->deltaRho_rem_q(i_q);
        if (flag_order[EW1QCD2]==true) 
            deltaRho_rem_q[i][EW1QCD2] = myThreeLoopQCD->deltaRho_rem_q(i_q);
        if (flag_order[EW2]==true) 
            deltaRho_rem_q[i][EW2] = myTwoLoopEW->deltaRho_rem_q(i_q);
        if (flag_order[EW2QCD1]==true) 
            deltaRho_rem_q[i][EW2QCD1] = myThreeLoopEW2QCD->deltaRho_rem_q(i_q);
        if (flag_order[EW3]==true) 
            deltaRho_rem_q[i][EW3] = myThreeLoopEW->deltaRho_rem_q(i_q);
        
        if (flag_order[EW1]==true) 
            deltaKappa_rem_l[i][EW1] = myOneLoopEW->deltaKappa_rem_l(i_l);
        if (flag_order[EW1QCD1]==true) 
            deltaKappa_rem_l[i][EW1QCD1] = myTwoLoopQCD->deltaKappa_rem_l(i_l);
        if (flag_order[EW1QCD2]==true) 
            deltaKappa_rem_l[i][EW1QCD2] = myThreeLoopQCD->deltaKappa_rem_l(i_l);
        if (flag_order[EW2]==true) 
            deltaKappa_rem_l[i][EW2] = myTwoLoopEW->deltaKappa_rem_l(i_l);
        if (flag_order[EW2QCD1]==true) 
            deltaKappa_rem_l[i][EW2QCD1] = myThreeLoopEW2QCD->deltaKappa_rem_l(i_l);
        if (flag_order[EW3]==true) 
            deltaKappa_rem_l[i][EW3] = myThreeLoopEW->deltaKappa_rem_l(i_l);
        
        if (flag_order[EW1]==true) 
            deltaKappa_rem_q[i][EW1] = myOneLoopEW->deltaKappa_rem_q(i_q);
        if (flag_order[EW1QCD1]==true) 
            deltaKappa_rem_q[i][EW1QCD1] = myTwoLoopQCD->deltaKappa_rem_q(i_q);
        if (flag_order[EW1QCD2]==true) 
            deltaKappa_rem_q[i][EW1QCD2] = myThreeLoopQCD->deltaKappa_rem_q(i_q);
        if (flag_order[EW2]==true) 
            deltaKappa_rem_q[i][EW2] = myTwoLoopEW->deltaKappa_rem_q(i_q);
        if (flag_order[EW2QCD1]==true) 
            deltaKappa_rem_q[i][EW2QCD1] = myThreeLoopEW2QCD->deltaKappa_rem_q(i_q);
        if (flag_order[EW3]==true) 
            deltaKappa_rem_q[i][EW3] = myThreeLoopEW->deltaKappa_rem_q(i_q);
    }
    
    /* The sum of the corrections */
    for (int j=0; j<orders_EW_size; j++) {
        for (int i=0; i<6; i++) {
            deltaRho_rem_l[i][orders_EW_size] += deltaRho_rem_l[i][j];
            deltaRho_rem_q[i][orders_EW_size] += deltaRho_rem_q[i][j];
            deltaKappa_rem_l[i][orders_EW_size] += deltaKappa_rem_l[i][j];
            deltaKappa_rem_q[i][orders_EW_size] += deltaKappa_rem_q[i][j];
        }
    }
}
    
    
////////////////////////////////////////////////////////////////////////


// double EWSM::mcMz() const {
//     double mc_at_mb = Mrun(quarks[BOTTOM].getMass(), quarks[CHARM].getMass(), 4.0);
//     double mc_at_Mz = Mrun(Mz, mc_at_mb, 5.0);

//     /* TEST */
//     //std::cout << "mc(mc)= " << quarks[CHARM].getMass() << std::endl;
//     //std::cout << "mc(mb)_LO+NLO= " << mc_at_mb << std::endl;
//     //std::cout << "mc(Mz)_LO+NLO= " << mc_at_Mz << std::endl;

//     return ( mc_at_Mz);
//     //return ( 0.563817 );// <--- used in ZFITTER with the effective mass mc=1.5
// }

// double EWSM::mbMz() const {
//     /* TEST */
//     //std::cout << "mb(mb)= " << quarks[BOTTOM].getMass() << std::endl;
//     //std::cout << "mb(Mz)_LO= " << mrun(Mz, quarks[BOTTOM].getMass(), 5.0, 0) << std::endl;
//     //std::cout << "mb(Mz)_LO+NLO= " << mrun(Mz, quarks[BOTTOM].getMass(), 5.0, 1) << std::endl;

//     return ( Mrun(Mz, quarks[BOTTOM].getMass(), 5.0));
//     //return (2.819440);// <--- used in ZFITTER with the effective mass mb=4.7
// }



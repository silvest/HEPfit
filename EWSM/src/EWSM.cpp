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
    
    for (int j=0; j<=orders_EW_size; j++) {
        EWSM::orders_EW j_order = (EWSM::orders_EW) j;
        DeltaAlpha_l[j] = orig.getDeltaAlpha_l(j_order);
        DeltaAlpha_t[j] = orig.getDeltaAlpha_t(j_order);
        DeltaRho[j] = orig.getDeltaRho(j_order);
        DeltaR_rem[j] = orig.getDeltaR_rem(j_order);
        for (int i=0; i<6; i++) {
            StandardModel::lepton i_l = (StandardModel::lepton) i;
            StandardModel::quark i_q = (StandardModel::quark) i;
            deltaRho_rem_l[j][i] = orig.getDeltaRho_rem_l(j_order,i_l);
            deltaRho_rem_q[j][i] = orig.getDeltaRho_rem_q(j_order,i_q);
            deltaKappa_rem_l[j][i] = orig.getDeltaKappa_rem_l(j_order,i_l);
            deltaKappa_rem_q[j][i] = orig.getDeltaKappa_rem_q(j_order,i_q);
        }
    }
    DeltaRbar_rem = orig.getDeltaRbar_rem();

    alphaMZ = orig.getAlphaMZ();
    Mw = orig.getMw();
    DeltaR = orig.getDeltaR();
    for (int i=0; i<6; i++) {
        StandardModel::lepton i_l = (StandardModel::lepton) i;
        StandardModel::quark i_q = (StandardModel::quark) i;
        rhoZ_l[i] = orig.getRhoZ_l(i_l);
        rhoZ_q[i] = orig.getRhoZ_q(i_q);
        kappaZ_l[i] = orig.getKappaZ_l(i_l);
        kappaZ_q[i] = orig.getKappaZ_q(i_q);
    }
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

void EWSM::TEST() const {
    std::cout << "Mw_tree = " << EWSMC->GetMW_tree() << std::endl;    
}

void EWSM::Compute() {
    
    Mw = EWSMC->GetMW_tree(); /* !! TEST !! */
    EWSMC->Compute(Mw);
    
    DeltaAlpha_l[TREE] = 0.0;
    DeltaAlpha_l[EW1] = myOneLoopEW->DeltaAlpha_l();
    DeltaAlpha_l[EW1QCD] = myTwoLoopQCD->DeltaAlpha_l();
    DeltaAlpha_l[EW1QCD2] = myThreeLoopQCD->DeltaAlpha_l();
    DeltaAlpha_l[EW2] = myTwoLoopEW->DeltaAlpha_l();
    DeltaAlpha_l[EW2QCD1] = myThreeLoopEW2QCD->DeltaAlpha_l();
    DeltaAlpha_l[EW3] = myThreeLoopEW->DeltaAlpha_l();
    
    DeltaAlpha_t[TREE] = 0.0;
    DeltaAlpha_t[EW1] = myOneLoopEW->DeltaAlpha_t();
    DeltaAlpha_t[EW1QCD] = myTwoLoopQCD->DeltaAlpha_t();
    DeltaAlpha_t[EW1QCD2] = myThreeLoopQCD->DeltaAlpha_t();
    DeltaAlpha_t[EW2] = myTwoLoopEW->DeltaAlpha_t();
    DeltaAlpha_t[EW2QCD1] = myThreeLoopEW2QCD->DeltaAlpha_t();
    DeltaAlpha_t[EW3] = myThreeLoopEW->DeltaAlpha_t();
    
    DeltaRho[TREE] = 0.0;
    DeltaRho[EW1] = myOneLoopEW->DeltaRho();
    DeltaRho[EW1QCD] = myTwoLoopQCD->DeltaRho();
    DeltaRho[EW1QCD2] = myThreeLoopQCD->DeltaRho();
    DeltaRho[EW2] = myTwoLoopEW->DeltaRho();
    DeltaRho[EW2QCD1] = myThreeLoopEW2QCD->DeltaRho();
    DeltaRho[EW3] = myThreeLoopEW->DeltaRho();
    
    DeltaR_rem[TREE] = 0.0;    
    DeltaR_rem[EW1] = myOneLoopEW->DeltaR_rem();
    DeltaR_rem[EW1QCD] = myTwoLoopQCD->DeltaR_rem();
    DeltaR_rem[EW1QCD2] = myThreeLoopQCD->DeltaR_rem();
    DeltaR_rem[EW2] = myTwoLoopEW->DeltaR_rem();
    DeltaR_rem[EW2QCD1] = myThreeLoopEW2QCD->DeltaR_rem();
    DeltaR_rem[EW3] = myThreeLoopEW->DeltaR_rem();    
    
    DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem();

    for (int i=0; i<6; i++) {
        StandardModel::lepton i_l = (StandardModel::lepton) i;
        StandardModel::quark i_q = (StandardModel::quark) i;

        deltaRho_rem_l[TREE][i_l] = 0.0;
        deltaRho_rem_l[EW1][i_l] = myOneLoopEW->deltaRho_rem_l(i_l);
        deltaRho_rem_l[EW1QCD][i_l] = myTwoLoopQCD->deltaRho_rem_l(i_l);
        deltaRho_rem_l[EW1QCD2][i_l] = myThreeLoopQCD->deltaRho_rem_l(i_l);
        deltaRho_rem_l[EW2][i_l] = myTwoLoopEW->deltaRho_rem_l(i_l);
        deltaRho_rem_l[EW2QCD1][i_l] = myThreeLoopEW2QCD->deltaRho_rem_l(i_l);
        deltaRho_rem_l[EW3][i_l] = myThreeLoopEW->deltaRho_rem_l(i_l);
        
        deltaRho_rem_q[TREE][i_q] = 0.0;
        deltaRho_rem_q[EW1][i_q] = myOneLoopEW->deltaRho_rem_q(i_q);
        deltaRho_rem_q[EW1QCD][i_q] = myTwoLoopQCD->deltaRho_rem_q(i_q);
        deltaRho_rem_q[EW1QCD2][i_q] = myThreeLoopQCD->deltaRho_rem_q(i_q);
        deltaRho_rem_q[EW2][i_q] = myTwoLoopEW->deltaRho_rem_q(i_q);
        deltaRho_rem_q[EW2QCD1][i_q] = myThreeLoopEW2QCD->deltaRho_rem_q(i_q);
        deltaRho_rem_q[EW3][i_q] = myThreeLoopEW->deltaRho_rem_q(i_q);
        
        deltaKappa_rem_l[TREE][i_l] = 0.0;
        deltaKappa_rem_l[EW1][i_l] = myOneLoopEW->deltaKappa_rem_l(i_l);
        deltaKappa_rem_l[EW1QCD][i_l] = myTwoLoopQCD->deltaKappa_rem_l(i_l);
        deltaKappa_rem_l[EW1QCD2][i_l] = myThreeLoopQCD->deltaKappa_rem_l(i_l);
        deltaKappa_rem_l[EW2][i_l] = myTwoLoopEW->deltaKappa_rem_l(i_l);
        deltaKappa_rem_l[EW2QCD1][i_l] = myThreeLoopEW2QCD->deltaKappa_rem_l(i_l);
        deltaKappa_rem_l[EW3][i_l] = myThreeLoopEW->deltaKappa_rem_l(i_l);
        
        deltaKappa_rem_q[TREE][i_q] = 0.0;
        deltaKappa_rem_q[EW1][i_q] = myOneLoopEW->deltaKappa_rem_q(i_q);
        deltaKappa_rem_q[EW1QCD][i_q] = myTwoLoopQCD->deltaKappa_rem_q(i_q);
        deltaKappa_rem_q[EW1QCD2][i_q] = myThreeLoopQCD->deltaKappa_rem_q(i_q);
        deltaKappa_rem_q[EW2][i_q] = myTwoLoopEW->deltaKappa_rem_q(i_q);
        deltaKappa_rem_q[EW2QCD1][i_q] = myThreeLoopEW2QCD->deltaKappa_rem_q(i_q);
        deltaKappa_rem_q[EW3][i_q] = myThreeLoopEW->deltaKappa_rem_q(i_q);
    }

    /* Total */
    DeltaAlpha_l[orders_EW_size] = 0.0;
    DeltaAlpha_t[orders_EW_size] = 0.0;
    DeltaRho[orders_EW_size] = 0.0;
    DeltaR_rem[orders_EW_size] = 0.0;
    for (int i=0; i<6; i++) {
        deltaRho_rem_l[orders_EW_size][i] = 0.0;
        deltaRho_rem_q[orders_EW_size][i] = 0.0;
        deltaKappa_rem_l[orders_EW_size][i] = 0.0;
        deltaKappa_rem_q[orders_EW_size][i] = 0.0;
    }
    for (int j=0; j<orders_EW_size; j++) {
        DeltaAlpha_l[orders_EW_size] += DeltaAlpha_l[j];
        DeltaAlpha_t[orders_EW_size] += DeltaAlpha_t[j];
        DeltaRho[orders_EW_size] += DeltaRho[j];
        DeltaR_rem[orders_EW_size] += DeltaR_rem[j];
        for (int i=0; i<6; i++) {
            deltaRho_rem_l[orders_EW_size][i] += deltaRho_rem_l[j][i];
            deltaRho_rem_q[orders_EW_size][i] += deltaRho_rem_q[j][i];
            deltaKappa_rem_l[orders_EW_size][i] += deltaKappa_rem_l[j][i];
            deltaKappa_rem_q[orders_EW_size][i] += deltaKappa_rem_q[j][i];
        }
    }

    double DeltaAlpha = DeltaAlpha_l[orders_EW_size] 
                        + DeltaAlpha_t[orders_EW_size] + SM.getDAle5Mz();
    alphaMZ = SM.getAle()/(1.0 - DeltaAlpha);

    //Mw
    //DeltaR    
    //rhoZ_l[6]
    //rhoZ_q[6]
    //kappaZ_l[6]
    //kappaZ_q[6]

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



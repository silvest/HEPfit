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
            deltaRho_rem_l[j][i] = 0.0;
            deltaRho_rem_q[j][i] = 0.0;
            deltaKappa_rem_l[j][i] = 0.0;
            deltaKappa_rem_q[j][i] = 0.0;
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

    schemeMw = orig.getSchemeMw();
    schemeRhoZ = orig.getSchemeRhoZ();
    schemeKappaZ = orig.getSchemeKappaZ();
    for (int j=0; j<orders_EW_size; j++) {
        EWSM::orders_EW j_order = (EWSM::orders_EW) j;
        flag_order[j] = orig.getFlag_order(j_order);
    }    
    
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

    alphaMz = orig.getAlphaMz();
    DeltaAlpha = orig.getDeltaAlpha();
    DeltaAlpha_l5q = orig.getDeltaAlpha_l5q();
    Mw_tree = orig.getMw_tree();
    Mw = orig.getMw();
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

void EWSM::setFlags(const schemes_EW schemeMw_i, 
                    const schemes_EW schemeRhoZ_i, const schemes_EW schemeKappaZ_i, 
                    const bool flag_order_i[orders_EW_size]) {
    
    schemeMw = schemeMw_i;
    schemeRhoZ = schemeRhoZ_i;
    schemeKappaZ = schemeKappaZ_i;
    for (int j=0; j<orders_EW_size; j++) {
        flag_order[j] = flag_order_i[j];
    }
}

void EWSM::Compute() {
    ComputeAlphaMz();
    ComputeMw();
    ComputeRhoZ();
    ComputeKappaZ();
}


////////////////////////////////////////////////////////////////////////

void EWSM::ComputeAlphaMz() {
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

    DeltaAlpha = DeltaAlpha_l[orders_EW_size] 
                 + DeltaAlpha_t[orders_EW_size] + SM.getDAle5Mz();
    DeltaAlpha_l5q = DeltaAlpha_l[orders_EW_size] + SM.getDAle5Mz(); 
    
    
#define TEST_DEBUG
#ifdef TEST_DEBUG
    DeltaAlpha = 0.10;
#endif    
    
    alphaMz = SM.getAle()/(1.0 - DeltaAlpha);
}

void EWSM::ComputeMw() {    
    double tmp = 4.0*M_PI*SM.getAle()/sqrt(2.0)/SM.getGF()/SM.getMz()/SM.getMz();
    Mw_tree = SM.getMz()/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp));
    
    if (schemeMw==APPROXIMATEFORMULA) {
        myApproximateFormulae = new ApproximateFormulae(SM, DeltaAlpha);
        Mw = myApproximateFormulae->Mw();        
        delete myApproximateFormulae;
    } else {
        EWSMC->Compute(Mw_tree);
        
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

        resumMw();
        
        // iterations 

        
    }

    EWSMC->Compute(Mw);
}    

void EWSM::ComputeRhoZ() {
    DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem();    
    
    for (int i=0; i<6; i++) {
        StandardModel::lepton i_l = (StandardModel::lepton) i;
        StandardModel::quark i_q = (StandardModel::quark) i;
        
        if (flag_order[EW1]==true) 
            deltaRho_rem_l[EW1][i_l] = myOneLoopEW->deltaRho_rem_l(i_l);
        if (flag_order[EW1QCD1]==true) 
            deltaRho_rem_l[EW1QCD1][i_l] = myTwoLoopQCD->deltaRho_rem_l(i_l);
        if (flag_order[EW1QCD2]==true) 
            deltaRho_rem_l[EW1QCD2][i_l] = myThreeLoopQCD->deltaRho_rem_l(i_l);
        if (flag_order[EW2]==true) 
            deltaRho_rem_l[EW2][i_l] = myTwoLoopEW->deltaRho_rem_l(i_l);
        if (flag_order[EW2QCD1]==true) 
            deltaRho_rem_l[EW2QCD1][i_l] = myThreeLoopEW2QCD->deltaRho_rem_l(i_l);
        if (flag_order[EW3]==true) 
            deltaRho_rem_l[EW3][i_l] = myThreeLoopEW->deltaRho_rem_l(i_l);
        
        if (flag_order[EW1]==true) 
            deltaRho_rem_q[EW1][i_q] = myOneLoopEW->deltaRho_rem_q(i_q);
        if (flag_order[EW1QCD1]==true) 
            deltaRho_rem_q[EW1QCD1][i_q] = myTwoLoopQCD->deltaRho_rem_q(i_q);
        if (flag_order[EW1QCD2]==true) 
            deltaRho_rem_q[EW1QCD2][i_q] = myThreeLoopQCD->deltaRho_rem_q(i_q);
        if (flag_order[EW2]==true) 
            deltaRho_rem_q[EW2][i_q] = myTwoLoopEW->deltaRho_rem_q(i_q);
        if (flag_order[EW2QCD1]==true) 
            deltaRho_rem_q[EW2QCD1][i_q] = myThreeLoopEW2QCD->deltaRho_rem_q(i_q);
        if (flag_order[EW3]==true) 
            deltaRho_rem_q[EW3][i_q] = myThreeLoopEW->deltaRho_rem_q(i_q);
    }

    /* The sum of the corrections */
    for (int j=0; j<orders_EW_size; j++) {
        for (int i=0; i<6; i++) {
            deltaRho_rem_l[orders_EW_size][i] += deltaRho_rem_l[j][i];
            deltaRho_rem_q[orders_EW_size][i] += deltaRho_rem_q[j][i];
        }
    }

    resumRhoZ();
}
    
void EWSM::ComputeKappaZ() {     
    DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem();    
    
    for (int i=0; i<6; i++) {
        StandardModel::lepton i_l = (StandardModel::lepton) i;
        StandardModel::quark i_q = (StandardModel::quark) i;
        
        if (flag_order[EW1]==true) 
            deltaKappa_rem_l[EW1][i_l] = myOneLoopEW->deltaKappa_rem_l(i_l);
        if (flag_order[EW1QCD1]==true) 
            deltaKappa_rem_l[EW1QCD1][i_l] = myTwoLoopQCD->deltaKappa_rem_l(i_l);
        if (flag_order[EW1QCD2]==true) 
            deltaKappa_rem_l[EW1QCD2][i_l] = myThreeLoopQCD->deltaKappa_rem_l(i_l);
        if (flag_order[EW2]==true) 
            deltaKappa_rem_l[EW2][i_l] = myTwoLoopEW->deltaKappa_rem_l(i_l);
        if (flag_order[EW2QCD1]==true) 
            deltaKappa_rem_l[EW2QCD1][i_l] = myThreeLoopEW2QCD->deltaKappa_rem_l(i_l);
        if (flag_order[EW3]==true) 
            deltaKappa_rem_l[EW3][i_l] = myThreeLoopEW->deltaKappa_rem_l(i_l);
        
        if (flag_order[EW1]==true) 
            deltaKappa_rem_q[EW1][i_q] = myOneLoopEW->deltaKappa_rem_q(i_q);
        if (flag_order[EW1QCD1]==true) 
            deltaKappa_rem_q[EW1QCD1][i_q] = myTwoLoopQCD->deltaKappa_rem_q(i_q);
        if (flag_order[EW1QCD2]==true) 
            deltaKappa_rem_q[EW1QCD2][i_q] = myThreeLoopQCD->deltaKappa_rem_q(i_q);
        if (flag_order[EW2]==true) 
            deltaKappa_rem_q[EW2][i_q] = myTwoLoopEW->deltaKappa_rem_q(i_q);
        if (flag_order[EW2QCD1]==true) 
            deltaKappa_rem_q[EW2QCD1][i_q] = myThreeLoopEW2QCD->deltaKappa_rem_q(i_q);
        if (flag_order[EW3]==true) 
            deltaKappa_rem_q[EW3][i_q] = myThreeLoopEW->deltaKappa_rem_q(i_q);
    }
    
    /* The sum of the corrections */
    for (int j=0; j<orders_EW_size; j++) {
        for (int i=0; i<6; i++) {
            deltaKappa_rem_l[orders_EW_size][i] += deltaKappa_rem_l[j][i];
            deltaKappa_rem_q[orders_EW_size][i] += deltaKappa_rem_q[j][i];
        }
    }
    
    resumKappaZ();

    // Write codes for Im[kappa_Z^f]
    
    
    
    
    
    
    /* Using the approximate formula for the real parts */
    if (schemeKappaZ==APPROXIMATEFORMULA) {
        myApproximateFormulae = new ApproximateFormulae(SM, DeltaAlpha);
        double sin2thetaEff_l[6], sin2thetaEff_q[6];
        for (int i=0; i<6; i++) {
            StandardModel::lepton i_l = (StandardModel::lepton) i;
            StandardModel::quark i_q = (StandardModel::quark) i;        
            sin2thetaEff_l[i] = myApproximateFormulae->sin2thetaEff(i_l);
            sin2thetaEff_q[i] = myApproximateFormulae->sin2thetaEff(i_q);
            kappaZ_l[i].real() = sin2thetaEff_l[i]/EWSMC->GetSW2(); 
            kappaZ_q[i].real() = sin2thetaEff_q[i]/EWSMC->GetSW2();
        }
        delete myApproximateFormulae;
    }   
}
        
void EWSM::resumMw() {
    double cW2_to_sW2 = EWSMC->GetCW2()/EWSMC->GetSW2();
    double R;

    switch (schemeMw) {
        case NORESUM: 
            R = DeltaAlpha_l5q - cW2_to_sW2*DeltaRho[orders_EW_size] 
                + DeltaR_rem[orders_EW_size];
            break;
        case OMSI:
            throw "Write codes for resummationMw()";
            break;
        case INTERMEDIATE:
            throw "Write codes for resummationMw()";            
            break;        
        case OMSII:
            throw "Write codes for resummationMw()";            
            break;
        default:
            throw "Error in resummationMw()";            
            break;
    }   

    double tmp = 4.0*M_PI*SM.getAle()/sqrt(2.0)/SM.getGF()/SM.getMz()/SM.getMz();
    Mw = SM.getMz()/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp*R));
}

void EWSM::resumRhoZ() {
    /* Real parts */
    switch (schemeRhoZ) {
        case NORESUM: 
            for (int i=0; i<6; i++) {
                rhoZ_l[i].real() = 1.0 + DeltaRho[orders_EW_size] 
                                   + deltaRho_rem_l[orders_EW_size][i].real();
                rhoZ_q[i].real() = 1.0 + DeltaRho[orders_EW_size] 
                                   + deltaRho_rem_q[orders_EW_size][i].real();
            }
            break;
        case OMSI:
            throw "Write codes for resummationRhoZ()";
            break;
        case INTERMEDIATE:
            throw "Write codes for resummationRhoZ()";
            break;        
        case OMSII:
            throw "Write codes for resummationRhoZ()";
            break;
        default:
            throw "Error in resummationRhoZ()";
            break;
    }

    /* Imaginary parts */
    for (int i=0; i<6; i++) {
        rhoZ_l[i].imag() = deltaRho_rem_l[orders_EW_size][i].imag();
        rhoZ_q[i].imag() = deltaRho_rem_q[orders_EW_size][i].imag();    
    }
}

void EWSM::resumKappaZ() {
    double cW2_to_sW2 = EWSMC->GetCW2()/EWSMC->GetSW2(); 
    
    /* Real parts */
    switch (schemeKappaZ) {
        case NORESUM: 
            for (int i=0; i<6; i++) {
                kappaZ_l[i].real() = 1.0 + cW2_to_sW2*DeltaRho[orders_EW_size] 
                                     + deltaKappa_rem_l[orders_EW_size][i].real();
                kappaZ_q[i].real() = 1.0 + cW2_to_sW2*DeltaRho[orders_EW_size] 
                                     + deltaKappa_rem_q[orders_EW_size][i].real();
            }
            break;
        case OMSI:
            throw "Write codes for resummationKappaZ()";
            break;
        case INTERMEDIATE:
            throw "Write codes for resummationKappaZ()";
            break;        
        case OMSII:
            throw "Write codes for resummationKappaZ()";
            break;
        case APPROXIMATEFORMULA:
            /* The real parts are given by the approximate formulae. 
             * See ComputeKappaZ() */
            break;
        default:
            throw "Error in resummationKappaZ()";
            break;
    }

    /* Imaginary parts */
    for (int i=0; i<6; i++) {
        kappaZ_l[i].imag() = deltaKappa_rem_l[orders_EW_size][i].imag();
        kappaZ_q[i].imag() = deltaKappa_rem_q[orders_EW_size][i].imag();    
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



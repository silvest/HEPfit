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
            deltaRho_rem_l[i][j] = orig.getDeltaRho_rem_l(i_l, j_order);
            deltaRho_rem_q[i][j] = orig.getDeltaRho_rem_q(i_q, j_order);
            deltaKappa_rem_l[i][j] = orig.getDeltaKappa_rem_l(i_l, j_order);
            deltaKappa_rem_q[i][j] = orig.getDeltaKappa_rem_q(i_q, j_order);
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

        /* Resummation */
        Mw = resumMw();
        
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
    }

    /* The sum of the corrections */
    for (int j=0; j<orders_EW_size; j++) {
        for (int i=0; i<6; i++) {
            deltaRho_rem_l[i][orders_EW_size] += deltaRho_rem_l[i][j];
            deltaRho_rem_q[i][orders_EW_size] += deltaRho_rem_q[i][j];
        }
    }

    /* Resummations */
    for (int i=0; i<6; i++) {
        rhoZ_l[i] = resumRhoZ(deltaRho_rem_l[i]);
        rhoZ_q[i] = resumRhoZ(deltaRho_rem_q[i]);
    }
}
    
void EWSM::ComputeKappaZ() {     
    DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem();    
    
    for (int i=0; i<6; i++) {
        StandardModel::lepton i_l = (StandardModel::lepton) i;
        StandardModel::quark i_q = (StandardModel::quark) i;
        
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
            deltaKappa_rem_l[i][orders_EW_size] += deltaKappa_rem_l[i][j];
            deltaKappa_rem_q[i][orders_EW_size] += deltaKappa_rem_q[i][j];
        }
    }
    
    /* Resummations */
    for (int i=0; i<6; i++) {
        kappaZ_l[i] = resumKappaZ(deltaKappa_rem_l[i]);
        kappaZ_q[i] = resumKappaZ(deltaKappa_rem_q[i]);
    }

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
        
double EWSM::resumMw() {
    double cW2_to_sW2 = EWSMC->GetCW2()/EWSMC->GetSW2();
    double R;

    if (DeltaR_rem[EW1QCD2]!=0.0) throw "Error in EWSM::resumMw()";
    if (DeltaR_rem[EW2QCD1]!=0.0) throw "Error in EWSM::resumMw()";
    if (DeltaR_rem[EW3]!=0.0) throw "Error in EWSM::resumMw()";
    
    switch (schemeMw) {
        case NORESUM: 
            // R = 1 + DeltaRho
            R = 1.0 + DeltaAlpha_l5q - cW2_to_sW2*DeltaRho[orders_EW_size] 
                + DeltaR_rem[orders_EW_size];
            break;
        case OMSI:
            // R = 1/(1 - DeltaRho)
            R = 1.0/(1.0 + cW2_to_sW2*DeltaRho[orders_EW_size])
                /(1.0 - DeltaAlpha_l5q 
                  - DeltaR_rem[EW1] - DeltaR_rem[EW1QCD1] - DeltaR_rem[EW2]);
            break;
        case INTERMEDIATE:
            // R = 1/(1 - DeltaRho)
            R = 1.0/( (1.0 + cW2_to_sW2*DeltaRho[orders_EW_size])
                      *(1.0 - DeltaAlpha_l5q - DeltaR_rem[EW1]) 
                      - DeltaR_rem[EW1QCD1] - DeltaR_rem[EW2] );
            break;        
        case OMSII:
            // R = 1/(1 - DeltaRho)
            R = 1.0/( (1.0 + cW2_to_sW2*DeltaRho[orders_EW_size])
                      *(1.0 - DeltaAlpha_l5q)
                      - (1.0 + cW2_to_sW2*DeltaRho[EW1])*DeltaR_rem[EW1]
                      - DeltaR_rem[EW1QCD1] - DeltaR_rem[EW2] );            
            break;
        default:
            throw "Error in resummationMw()";            
            break;
    }   

    double tmp = 4.0*M_PI*SM.getAle()/sqrt(2.0)/SM.getGF()/SM.getMz()/SM.getMz();
    return (SM.getMz()/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp*R)));
}

complex EWSM::resumRhoZ(const complex deltaRho_rem[orders_EW_size+1]) {
    complex rhoZ;
    
    double deltaRho_rem_G_real = deltaRho_rem[EW1].real() 
                                 + deltaRho_rem[EW1QCD1].real();
    if (deltaRho_rem[EW1QCD2]!=0.0) throw "Error in EWSM::resumRhoZ()";
    if (deltaRho_rem[EW2QCD1]!=0.0) throw "Error in EWSM::resumRhoZ()";    
    if (deltaRho_rem[EW3]!=0.0) throw "Error in EWSM::resumRhoZ()";    
    
    /* Real parts */
    switch (schemeRhoZ) {
        case NORESUM: 
            rhoZ.real() = 1.0 + DeltaRho[orders_EW_size] 
                          + deltaRho_rem[orders_EW_size].real();
            break;
        case OMSI:
            rhoZ.real() = (1.0 + deltaRho_rem_G_real + deltaRho_rem[EW2].real())
                          /(1.0 - DeltaRho[orders_EW_size]*(1.0 - DeltaRbar_rem));
            break;
        case INTERMEDIATE:
            rhoZ.real() = (1.0 + deltaRho_rem_G_real)
                          /(1.0 - DeltaRho[orders_EW_size]*(1.0 - DeltaRbar_rem))
                          + deltaRho_rem[EW2].real();            
            break;        
        case OMSII:
            rhoZ.real() = 1.0 + DeltaRho[orders_EW_size] 
                          + pow(DeltaRho[EW1], 2.0) - DeltaRho[EW1]*DeltaRbar_rem
                          + deltaRho_rem_G_real*(1.0 + DeltaRho[EW1]) 
                          + deltaRho_rem[EW2].real();  
            break;
        default:
            throw "Error in resummationRhoZ()";
            break;
    }

    /* Imaginary parts */
    rhoZ.imag() = deltaRho_rem[orders_EW_size].imag();

    return rhoZ;
}

complex EWSM::resumKappaZ(const complex deltaKappa_rem[orders_EW_size+1]) {
    double cW2_to_sW2 = EWSMC->GetCW2()/EWSMC->GetSW2(); 
    complex kappaZ;
    
    double deltaKappa_rem_G_real = deltaKappa_rem[EW1].real() 
                                   + deltaKappa_rem[EW1QCD1].real() 
                                   + deltaKappa_rem[EW1QCD2].real() ;
    if (deltaKappa_rem[EW2QCD1]!=0.0) throw "Error in EWSM::resumKappaZ()";    
    if (deltaKappa_rem[EW3]!=0.0) throw "Error in EWSM::resumKappaZ()";     
    
    /* Real parts */
    switch (schemeKappaZ) {
        case NORESUM: 
            kappaZ.real() = 1.0 + cW2_to_sW2*DeltaRho[orders_EW_size] 
                            + deltaKappa_rem[orders_EW_size].real();
            break;
        case OMSI:
            kappaZ.real() = (1.0 + deltaKappa_rem_G_real 
                             + deltaKappa_rem[EW2].real())
                            *(1.0 + cW2_to_sW2*DeltaRho[orders_EW_size]
                                    *(1.0 - DeltaRbar_rem));
            break;
        case INTERMEDIATE:
            kappaZ.real() = (1.0 + deltaKappa_rem_G_real)
                            *(1.0 + cW2_to_sW2*DeltaRho[orders_EW_size]
                                    *(1.0 - DeltaRbar_rem))
                            + deltaKappa_rem[EW2].real();
            break;        
        case OMSII:
            kappaZ.real() = 1.0 + cW2_to_sW2*DeltaRho[orders_EW_size]
                            - cW2_to_sW2*DeltaRho[EW1]*DeltaRbar_rem
                            +  deltaKappa_rem_G_real
                               *(1.0 + cW2_to_sW2*DeltaRho[EW1])
                            + deltaKappa_rem[EW2].real();
            break;
        case APPROXIMATEFORMULA:
            /* The real parts are given by the approximate formulae. 
             * See ComputeKappaZ() */
            kappaZ.real() = 0.0; // dummy
            break;
        default:
            throw "Error in resummationKappaZ()";
            break;
    }

    /* Imaginary parts */
    kappaZ.imag() = deltaKappa_rem[orders_EW_size].imag();

    return kappaZ;
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



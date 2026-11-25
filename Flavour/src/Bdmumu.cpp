/*
 * Copyright (C) 2012 SusyFit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Bdmumu.h"
#include "StandardModel.h"
#include "AmpDB2.h"
#include <limits>
#include "EvolBsmm.h"
#include "HeffDB1.h"
#include "NPSMEFTd6GeneralMatching.h"

Bdmumu::Bdmumu(const StandardModel& SM_i, int obsFlag, QCD::lepton lep_i)
: ThObservable(SM_i)
{
    lep = lep_i;
    if(lep == QCD::MU) leptonindex = 1;
    else if(lep == QCD::ELECTRON) leptonindex = 0;
    else if(lep == QCD::TAU) leptonindex = 2;
    else throw std::runtime_error("Bdmumu::Bdmumu(): lepton type not defined. Can be only MU or ELECTRON");
    if (obsFlag > 0 and obsFlag < 5) obs = obsFlag;
    else throw std::runtime_error("obsFlag in Bdmumu(myFlavour, obsFlag) called from ThFactory::ThFactory() can only be 1 (BR) or 2 (BRbar) or 3 (Amumu) or 4 (Smumu)");
    SM.initializeMeson(QCD::B_D);
};

double Bdmumu::computeThValue()
{
    double FBd = SM.getMesons(QCD::B_D).getDecayconst();

    computeObs(FULLNLO, FULLNLO_QED, obs >= 3);

    double PRF = pow(coupling, 2.) / M_PI /8. / SM.getMesons(QCD::B_D).computeWidth() * pow(FBd, 2.) * pow(mlep, 2.) * mBd * beta;
    // The Bd time-integrated rate uses y_d = 0.
    if (obs == 1 || obs == 2) return PRF * ampSq;
    if (obs == 3) return( Amumu );
    if (obs == 4) return( Smumu );
    
    throw std::runtime_error("Bdmumu::computeThValue(): Observable type not defined. Can be only any of (1,2,3,4)");
    return (EXIT_FAILURE);
}

void Bdmumu::computeObs(orders order, orders_qed order_qed, bool withMixing)
{   
    double mu = SM.getMub();  
    
    mlep = SM.getLeptons(lep).getMass();
    mBd = SM.getMesons(QCD::B_D).getMass();
    mW = SM.Mw();
    const double GF = SM.getGF();
    coupling = GF * GF * mW * mW / M_PI / M_PI;
    mb = SM.Mrun(mu, SM.getQuarks(QCD::BOTTOM).getMass_scale(),
            SM.getQuarks(QCD::BOTTOM).getMass(), QCD::BOTTOM, FULLNNLO);
    md = SM.Mrun(mu, SM.getQuarks(QCD::DOWN).getMass_scale(),
            SM.getQuarks(QCD::DOWN).getMass(), QCD::DOWN, FULLNNLO);
    chiral = pow(mBd, 2.) / 2. / mlep * mb / (mb + md);
    beta = sqrt(1. - pow(2. * mlep / mBd, 2.));
    computeAmpSq(order, order_qed, mu);
    if (!withMixing) return;
    if (ampSq == 0.) {
        Amumu = Smumu = std::numeric_limits<double>::quiet_NaN();
        return;
    }
    const double mixingPhase = AmpDB2::mixingPhase(SM, QCD::B_D);
    Amumu = (absP * absP * cos(2. * argP - mixingPhase) -  absS * absS * cos(2. * argS - mixingPhase)) / (absP * absP + absS * absS);
    /* LHCb convention: S = [|P|^2 sin(2 phi_P - phi^NP) - |S|^2 sin(2 phi_S - phi^NP)]/(|P|^2 + |S|^2) (De Bruyn et al.,
     * arXiv:1204.1737, eqs. (16), (19); NP mixing phase phi^NP as in Buras et al., arXiv:1303.3820), with the CP asymmetry
     * [Gamma(B) - Gamma(Bbar)]/[...] = S sin(Delta m t)/[...] and lambda_SM = -1. M21, and hence mixingPhase, is the
     * conjugate of the standard M12 (mixingPhase - mixingPhase_SM = -phi^NP) and P carries the same CKM phase
     * (2 argP = mixingPhase in the SM), so S = sin(mixingPhase - 2 argP) for the P term. */
    Smumu = (absP * absP * sin(mixingPhase - 2. * argP) -  absS * absS * sin(mixingPhase - 2. * argS)) / (absP * absP + absS * absS);
}

double Bdmumu::computeAmumu(orders order)
{
    computeObs(FULLNLO, FULLNLO_QED);
    return(Amumu);
}

double Bdmumu::computeSmumu(orders order)
{
    computeObs(FULLNLO, FULLNLO_QED);
    return(Smumu);
}

void Bdmumu::computeAmpSq(orders order, orders_qed order_qed, double mu)
{  
    if (SM.getFlavour().getHDB1().getCoeffdmumu().getOrder() < order % 3){
        std::stringstream out;
        out << order;
        throw std::runtime_error("Bdmumu::computeAmpSq(): required cofficient of "
                                 "order " + out.str() + " not computed");
    }
    allcoeff = SM.getFlavour().ComputeCoeffdmumu(mu, NDR);

    EvolBsmm& evolbdmm = SM.getFlavour().getHDB1().getUBdmm();
    double alsmu = evolbdmm.alphatilde_s(mu);
    double alemu = evolbdmm.alphatilde_e(mu);
//    double alemu = SM.ale_OS(mu)/4./M_PI; // to be checked
    gslpp::matrix<gslpp::complex> Vckm = SM.getVCKM();
//    double sw2 =  (M_PI * SM.getAle() ) / ( sqrt(2.) * SM.getGF() * SM.Mw() * SM.Mw()) ;
    
    C_10p = 0.;
    C_S = 0.;
    C_Sp = 0.;
    C_P = 0.;
    C_Pp = 0.;
   
    if((order == FULLNLO) && (order_qed == FULLNLO_QED)){
    
    switch(order_qed) {
        case FULLNLO_QED:
        {
            C_10 = (*(allcoeff[LO]))(7) /alemu  + (*(allcoeff[NLO]))(7) * alsmu/alemu 
                    + (*(allcoeff[NNLO]))(7) * alsmu * alsmu/alemu + (*(allcoeff[LO_QED ]))(7) /alsmu
                    + (*(allcoeff[NLO_QED11]))(7) + (*(allcoeff[NLO_QED02]))(7) * alemu /alsmu /alsmu 
                    + (*(allcoeff[NLO_QED21]))(7) * alsmu 
                    + (*(allcoeff[NLO_QED12]))(7) * alemu /alsmu+ (*(allcoeff[NLO_QED22]))(7) * alemu;
            
            gslpp::complex NPfactor = 2. * coupling;

            if( const NPSMEFTd6GeneralMatching* matching = dynamic_cast<const NPSMEFTd6GeneralMatching*>(&SM.getMatching()))
            {
                /*
                 * Run the NP LEFT coefficients from the EW matching scale to mu.
                 * For QED we keep only the one-loop leading logarithm. In the
                 * axial combinations entering B_d -> l+l-, the QED penguin
                 * terms in the JMS RGEs cancel. With q_e=-1 and q_d=-1/3,
                 * the remaining vector-current mixing is 4 * alpha_e/(4 pi)
                 * * log(muEW/mu).
                 *
                 * The scalar LEFT coefficients have the same QCD anomalous
                 * dimension as the quark mass, so their QCD running is
                 * implemented with the exact HEPfit mass-running machinery.
                 */
                const double muEW = SM.getMuw();
                const double qedLog = evolbdmm.alphatilde_e(muEW) * log(muEW / mu);
                const double vectorQEDMix = 4. * qedLog;

                const double mbEW = SM.Mrun(muEW, SM.getQuarks(QCD::BOTTOM).getMass_scale(),
                        SM.getQuarks(QCD::BOTTOM).getMass(), QCD::BOTTOM, FULLNNLO);
                const double scalarQCD = mb / mbEW;
                const double scalarQED = 1. + 20. / 3. * qedLog;
                const double scalarRunning = scalarQCD * scalarQED;

                const gslpp::complex LdeVLR = matching->getCdeVLR(0,2,leptonindex,leptonindex);
                const gslpp::complex LedVLL = matching->getCedVLL(leptonindex,leptonindex,0,2);
                const gslpp::complex LedVLR = matching->getCedVLR(leptonindex,leptonindex,0,2);
                const gslpp::complex LedVRR = matching->getCedVRR(leptonindex,leptonindex,0,2);

                /* The JMS coefficients with quark indices (0,2) multiply (dbar ... b) operators, i.e. the
                 * Bbar_d decay. The SM coefficients of this class (CMbdmm, lambda_t = computelamt_d() = V_td V_tb^*) and
                 * M21 refer instead to the (bbar ... d) operators, as CddVLL[2][0][2][0] in NPSMEFTd6GeneralMatching::CMdbd2().
                 * The NP shifts, built below in the dbar-b convention (JMS -> C10, C_S, C_P dictionary), therefore
                 * enter the amplitude complex conjugated. The hermitian conjugation maps each shift to its own
                 * conjugate (the sign of (lbar gamma5 l)^dagger is compensated by P_R <-> P_L in the quark density).
                 */
                C_10 = C_10 + ((LdeVLR - LedVLL
                        + vectorQEDMix * (LdeVLR + LedVLL)) / NPfactor).conjugate();
                C_10p = C_10p + ((-LedVLR + LedVRR
                        - vectorQEDMix * (LedVLR + LedVRR)) / NPfactor).conjugate();

                /* JMS also contains T^{RR}_{ed} -> S^{RR}_{ed} QED mixing.
                 * It vanishes for the present tree-level SMEFT matching, but
                 * keeping it here makes the LEFT running complete at O(alpha_e).
                 * The tensor-induced term is kept as a strict QED leading log;
                 * mixed QCD-QED logarithms are beyond the accuracy used here.
                 */
                const double tensorToScalarQED = 32. * qedLog;
                const gslpp::complex LedSRR = scalarRunning
                        * matching->getCedSRR(leptonindex,leptonindex,0,2)
                        + tensorToScalarQED
                        * matching->getCedTRR(leptonindex,leptonindex,0,2);
                const gslpp::complex LedSRL = scalarRunning
                        * matching->getCedSRL(leptonindex,leptonindex,0,2);
                const gslpp::complex LedSRRhc = scalarRunning
                        * matching->getCedSRR(leptonindex,leptonindex,2,0).conjugate()
                        + tensorToScalarQED
                        * matching->getCedTRR(leptonindex,leptonindex,2,0).conjugate();
                const gslpp::complex LedSRLhc = scalarRunning
                        * matching->getCedSRL(leptonindex,leptonindex,2,0).conjugate();

                C_S = C_S + ((LedSRR + LedSRLhc) / NPfactor).conjugate();
                C_Sp = C_Sp + ((LedSRL + LedSRRhc) / NPfactor).conjugate();
                C_P = C_P + ((LedSRR - LedSRLhc) / NPfactor).conjugate();
                C_Pp = C_Pp + ((LedSRL - LedSRRhc) / NPfactor).conjugate();
            }

            gslpp::complex CC_P = (C_10 - C_10p + mBd * mBd / ( 2. * mlep * (mb + md)) * (C_P - C_Pp) );

                
            absP = CC_P.abs(); 
            argP = CC_P.arg();
            
            gslpp::complex CC_S = ( beta * mBd * mBd / ( 2. * mlep * (mb + md)) * (C_S - C_Sp) );

            absS = CC_S.abs();
            argS = CC_S.arg();
    
            
            ampSq = absP * absP + absS * absS;
                  
        }
        break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Bdmumu::computeAmpSq(): order " + out.str() + " not implemented");;
        }
    }
        
}

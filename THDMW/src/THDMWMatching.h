/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef THDMWMATCHING_H
#define THDMWMATCHING_H

#include <Polylogarithms.h>
#include "gslpp.h"
#include "StandardModelMatching.h"

class THDMW;

/**
 * @class THDMWMatching
 * @ingroup THDMW
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class THDMWMatching : public StandardModelMatching {
public:
    THDMWMatching(const THDMW & THDMW_i);


    /**
     * @return C10 Wilson coefficient of  B_q \to l \bar{l}\f$ according to 1504.00839.
     */
    virtual gslpp::complex C10NP(double xt, double xS, gslpp::complex etaU);
    
    /**
     * @return CS Wilson coefficient of  B_q \to l \bar{l}\f$ according to 1504.00839.
     */
    virtual gslpp::complex CSNP(double nu1, double xh, double xt, double xS, gslpp::complex etaU,  gslpp::complex etaD);
    
    
    /**
     * @return CP Wilson coefficient of  B_q \to l \bar{l}\f$ according to 1504.00839.
     */
    virtual gslpp::complex CPNP(double xt, double xS, gslpp::complex etaU,  gslpp::complex etaD);
    
    
    
    
    /**
     * @return THDMW Wilson coefficients for \f$ B_s \to \bar{B_s}\f$ according to 1504.00839 
     */
    virtual  std::vector<WilsonCoefficient>& CMdbs2();
    virtual  std::vector<WilsonCoefficient>& CMdbsp2();

    

    virtual std::vector<WilsonCoefficient>& CMBMll(QCD::lepton lepton);
    
    
    
    
    /**
     * @return C^(NP,ct)_(VLL, Eta_d Eta^*_u) according to 1504.00839
     */
    virtual double CNPVLLctEtadEtasu(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,ct)_(VLL, Eta_u^4) according to 1504.00839
     */
    virtual double CNPVLLctEtau4(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,ct)_(VLL, Eta_u^2) according to 1504.00839
     */
    virtual double CNPVLLctEtau2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,ct)_(VLL) according to 1504.00839
     */
    virtual gslpp::complex CNPVLLct(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD);
    
    /**
     * @return C^(NP,tt)_(VLL, Eta_d Eta^*_u) according to 1504.00839
     */
    virtual double CNPVLLttEtadEtasu(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,tt)_(VLL, Eta_u^4) according to 1504.00839
     */
    virtual double CNPVLLttEtau4(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,tt)_(VLL, Eta_u^2) according to 1504.00839
     */
    virtual double CNPVLLttEtau2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,tt)_(VLL) according to 1504.00839
     */
    virtual gslpp::complex CNPVLLtt(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD);

    /**
     * @return C^(NP,cc)_(VLL, Eta_d Eta^*_u) according to 1504.00839
     */
    virtual double CNPVLLccEtadEtasu(double xc, double xb, double xS);
    
    /**
     * @return C^(NP,cc)_(VLL, Eta_u^2) according to 1504.00839
     */
    virtual double CNPVLLccEtau2(double xc, double xb, double xS);
    
    /**
     * @return C^(NP,cc)_(VLL) according to 1504.00839
     */
    virtual gslpp::complex CNPVLLcc(double xc, double xb, double xS,  gslpp::complex etaU, gslpp::complex etaD);

    /**
     * @return \bar{C^(NP)_(VLL)} according to 1504.00839
     */
    virtual gslpp::complex CNPVLL(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD);
    
    
    
    
    
    
    /**
     * @return C^(NP,ct)_(SRR1, Eta_d Eta^*_u Eta_u^2) according to 1504.00839
     */
    virtual double CNPSRR1ctEtadEtasuEtau2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,ct)_(SRR1, (Eta_d Eta^*_u)^2) according to 1504.00839
     */
    virtual double CNPSRR1ctEtad2Etasu2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,ct)_(SRR1, (Eta_d Eta^*_u)) according to 1504.00839
     */
    virtual double CNPSRR1ctEtadEtasu(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,ct)_(SRR1, (|Eta_u|^4) according to 1504.00839
     */
    virtual double CNPSRR1ctEtau4(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,ct)_(SRR1, (|Eta_u|^2) according to 1504.00839
     */
    virtual double CNPSRR1ctEtau2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,ct)_(SRR1) according to 1504.00839
     */
    virtual gslpp::complex CNPSRR1ct(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD);
 
    
    
    
    
    
    
    /**
     * @return C^(NP,tt)_(SRR1, Eta_d Eta^*_u Eta_u^2) according to 1504.00839
     */
    virtual double CNPSRR1ttEtadEtasuEtau2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,tt)_(SRR1, (Eta_d Eta^*_u)^2) according to 1504.00839
     */
    virtual double CNPSRR1ttEtad2Etasu2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,tt)_(SRR1, (Eta_d Eta^*_u)) according to 1504.00839
     */
    virtual double CNPSRR1ttEtadEtasu(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,tt)_(SRR1, (|Eta_u|^4) according to 1504.00839
     */
    virtual double CNPSRR1ttEtau4(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,tt)_(SRR1, (|Eta_u|^2) according to 1504.00839
     */
    virtual double CNPSRR1ttEtau2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,tt)_(SRR1) according to 1504.00839
     */
    virtual gslpp::complex CNPSRR1tt(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD);
 
    /**
     * @return C^(NP)_(SRR1) according to 1504.00839
     */
    virtual gslpp::complex CNPSRR1(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD);
 
    
 
    
    
    
    
    
    /**
     * @return C^(NP,ct)_(SRR2, Eta_d Eta^*_u Eta_u^2) according to 1504.00839
     */
    virtual double CNPSRR2ctEtadEtasuEtau2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,ct)_(SRR2, (Eta_d Eta^*_u)^2) according to 1504.00839
     */
    virtual double CNPSRR2ctEtad2Etasu2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,ct)_(SRR2, (Eta_d Eta^*_u)) according to 1504.00839
     */
    virtual double CNPSRR2ctEtadEtasu(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,ct)_(SRR2, (|Eta_u|^4) according to 1504.00839
     */
    virtual double CNPSRR2ctEtau4(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,ct)_(SRR2, (|Eta_u|^2) according to 1504.00839
     */
    virtual double CNPSRR2ctEtau2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,ct)_(SRR2) according to 1504.00839
     */
    virtual gslpp::complex CNPSRR2ct(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD);
 
    
    
    
    
    
    /**
     * @return C^(NP,tt)_(SRR2, Eta_d Eta^*_u Eta_u^2) according to 1504.00839
     */
    virtual double CNPSRR2ttEtadEtasuEtau2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,tt)_(SRR2, (Eta_d Eta^*_u)^2) according to 1504.00839
     */
    virtual double CNPSRR2ttEtad2Etasu2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,tt)_(SRR2, (Eta_d Eta^*_u)) according to 1504.00839
     */
    virtual double CNPSRR2ttEtadEtasu(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,tt)_(SRR2, (|Eta_u|^4) according to 1504.00839
     */
    virtual double CNPSRR2ttEtau4(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,tt)_(SRR2, (|Eta_u|^2) according to 1504.00839
     */
    virtual double CNPSRR2ttEtau2(double xc, double xb, double xt, double xS);
    
    /**
     * @return C^(NP,tt)_(SRR2) according to 1504.00839
     */
    virtual gslpp::complex CNPSRR2tt(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD);
 
    /**
     * @return C^(NP)_(SRR1) according to 1504.00839
     */
    virtual gslpp::complex CNPSRR2(double xc, double xb, double xt, double xS,  gslpp::complex etaU, gslpp::complex etaD);
 
    
    
    
    
    
    
    
    
    
    /**
     * @return f1 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f1(double xc, double xt, double xS);
  
    /**
     * @return f2 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f2(double xb, double xt, double xS);
    
    /**
     * @return f3 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f3(double xb, double xt, double xS);
    
    /**
     * @return f4 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f4(double xb, double xt, double xS);
    
    /**
     * @return f5 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f5(double xc, double xt, double xS);
    
    /**
     * @return f6 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f6(double xb, double xt, double xS);
    
    /**
     * @return f7 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f7(double xb, double xt, double xS);
    
    /**
     * @return f8 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f8(double xb, double xt);
    
    /**
     * @return f9 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f9(double xc ,double xb, double xt, double xS);
    
    /**
     * @return f10 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f10( double xt, double xS);
    
    /**
     * @return f11 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f11( double xt, double xS);
    
    /**
     * @return f12 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f12( double xt, double xS);
    
    /**
     * @return f13 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f13( double xt, double xS);
    
    /**
     * @return f14 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f14(double xb, double xt, double xS);
    
    /**
     * @return f15 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f15(double xb, double xt, double xS);
    
    /**
     * @return f16 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f16(double xb, double xt, double xS);
    
    /**
     * @return f17 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f17(double xt, double xS);
    
    /**
     * @return f18 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f18(double xb, double xt, double xS);
    
    /**
     * @return f19 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f19(double xb, double xt, double xS);
    
    /**
     * @return f20 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f20(double xb, double xt, double xS);
    
    /**
     * @return f21 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f21(double xb, double xt);
    
    /**
     * @return f22 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f22(double xb, double xt, double xS);
    
    /**
     * @return f23 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f23(double xb, double xt, double xS);
    
    /**
     * @return f24 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f24(double xt, double xS);
    
    /**
     * @return f25 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f25(double xb, double xt, double xS);
    
    /**
     * @return f26 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f26(double xb, double xt, double xS);
    
    /**
     * @return f27 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f27(double xb, double xt, double xS);
    
    /**
     * @return f28 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f28(double xb, double xt);
    
    /**
     * @return f29 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f29(double xt, double xS);
    
    /**
     * @return f30 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f30(double xb, double xt, double xS);
    
    /**
     * @return f31 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f31(double xb, double xt, double xS);
    
    /**
     * @return f32 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f32(double xc, double xb, double xt, double xS);
    
    /**
     * @return f33 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f33(double xc, double xb, double xt, double xS);
    
    /**
     * @return f34 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f34(double xt, double xS);
    
    /**
     * @return f35 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f35(double xc, double xb, double xt, double xS);
    
    /**
     * @return f36 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f36(double xt, double xS);
    
    /**
     * @return f37 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f37(double xt, double xS);
    
    /**
     * @return f38 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f38(double xb, double xt, double xS);
    
    /**
     * @return f39 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f39(double xb, double xt, double xS);
    
    /**
     * @return f40 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f40(double xb, double xt, double xS);
    
    /**
     * @return f41 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f41(double xt, double xS);
    
    /**
     * @return f42 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f42(double xt, double xS);
    
    /**
     * @return f43 needed to calculate  \f$ B_s \to \bar{B_s}\f$ according to 1504.00839
     */
    virtual  double f43(double xt, double xS);
    
    
    
    
     const Polylogarithms getPolyLog() const
       {
           return PolyLog;
       }

    void updateTHDMWParameters();

private:
    const THDMW & myTHDMW;

    gslpp::matrix<gslpp::complex> myCKM;
    WilsonCoefficient mcBMll,mcbsg , mcdbs2, mcdbsp2;

    //double GF, mMU;
    //gslpp::complex CWbsgArrayLO[8], CWbsgArrayNLO[8], CWbsgArrayNNLO[8];
    //double mtbsg, mhpbsg, mubsg; // caching
    //gslpp::complex su, sd, sl; // caching

    const Polylogarithms PolyLog;
};





#endif /* THDMWMATCHING_H */

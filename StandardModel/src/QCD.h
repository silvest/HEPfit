/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef QCD_H
#define	QCD_H

#include "Model.h"
#include "Meson.h"
#include "OrderScheme.h"
#define MEPS 1.e-10 // mass precision

/**
 * @class QCD
 * @ingroup StandardModel
 * @brief A class for QCD. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class QCD: public Model {
public:
    enum meson {B_D, B_S, B_P, K_0, K_P, D_0, P_0, P_P, MESON_END}; 
    enum quark {UP,DOWN,CHARM,STRANGE,TOP,BOTTOM};

    static const int NQCDvars = 76;//43;//26;

    /**
     * array containing the labels under which all QCD parameters must be
     * stored in a Parameters object
     */
    static const std::string QCDvars[NQCDvars];

    virtual std::string ModelName() const 
    {
        return "QCD";
    }
    
    QCD() 
    : BBs(5), BBd(5), BD(5), BK(5), BKd1(10), BKd3(10)
    {
        Nc=3.;
        CF = Nc/2.-1./(2.*Nc);
        quarks[UP].setCharge(2./3.);
        quarks[UP].setMass_scale(2.);
        quarks[UP].setIsospin(1./2.);
        quarks[CHARM].setCharge(2./3.);
        quarks[CHARM].setIsospin(1./2.);    
        quarks[TOP].setCharge(2./3.);
        quarks[TOP].setIsospin(1./2.);    
        quarks[DOWN].setCharge(-1./3.);
        quarks[DOWN].setMass_scale(2.);
        quarks[DOWN].setIsospin(-1./2.);
        quarks[STRANGE].setCharge(-1./3.);
        quarks[STRANGE].setMass_scale(2.);
        quarks[STRANGE].setIsospin(-1./2.);   
        quarks[BOTTOM].setCharge(-1./3.); 
        quarks[BOTTOM].setIsospin(-1./2.);
        //to be moved to the Als class
        for (int i = 0; i < CacheSize; i++)
        {
            for (int j = 0; j < 8; j++)
                als_cache[j][i] = 0.;
            for (int j = 0; j < 10; j++)
                mrun_cache[j][i] = 0.;
            for (int j = 0; j < 4; j++)
                logLambda5_cache[j][i] = 0.;
            for (int j = 0; j < 4; j++)
                mp2mbar_cache[j][i] = 0.;
        }
    };

    ////////////////////////////////////////////////////////////////////////

    virtual bool SetFlag(const std::string, const bool&);
    
    virtual bool PreUpdate();
     
    virtual bool PostUpdate();      
    
    /**
     * updates the QCD parameters found in the argument
     * @param a map containing the parameters (all as double) to be updated
     */
    virtual bool Update(const std::map<std::string, double>&);
    
    /**
     * Checks that all required parameters are present
     * @param a map containing the parameters (all as double) to be updated
     */
    virtual bool CheckParameters(const std::map<std::string, double>&);
    
    /**
     * Initializes the QCD parameters found in the argument
     * @param a map containing the parameters (all as double) to be updated
     * "AlsMz"
     * "Mz"
     * "mup" @f$\overline{\mathrm{MS}}@f$ mass @f$m_u(2\mathrm{GeV})@f$
     * "mdown" @f$\overline{\mathrm{MS}}@f$ mass @f$m_d(2\mathrm{GeV})@f$
     * "mcharm" @f$\overline{\mathrm{MS}}@f$ mass @f$m_c(m_c)@f$
     * "mstrange" @f$\overline{\mathrm{MS}}@f$ mass @f$m_s(2\mathrm{GeV})@f$
     * "mtop" the top quark pole mass @f$m_t@f$
     * "mbottom" @f$\overline{\mathrm{MS}}@f$ mass @f$m_b(m_b)@f$
     * "mut"
     * "mub"
     * "muc"
     * "MBd"
     * "MBs"
     * "MBp"
     * "MK0"
     * "MKp"
     * "MD"
     * "FBs"
     * "FBsoFBd"
     * "FD"
     * "BBsoBBd"
     * "BBs1"
     * "BBs2"
     * "BBs3"
     * "BBs4"
     * "BBs5"
     * "BBsscale"
     * "BBsscheme"
     * "BD1"
     * "BD2"
     * "BD3"
     * "BD4"
     * "BD5"
     * "BDscale"
     * "BDscheme"
     * "BK1"
     * "BK2"
     * "BK3"
     * "BK4"
     * "BK5"
     * "BKscale"
     * "BKscheme"
     * "BK(1/2)1"
     * "BK(1/2)2"
     * "BK(1/2)3"
     * "BK(1/2)4"
     * "BK(1/2)5"
     * "BK(1/2)6"
     * "BK(1/2)7"
     * "BK(1/2)8"
     * "BK(1/2)9"
     * "BK(1/2)10"
     * "BKd_scale"
     * "BKd_scheme"
     * "BK(3/2)1"
     * "BK(3/2)2"
     * "BK(3/2)3"
     * "BK(3/2)4"
     * "BK(3/2)5"
     * "BK(3/2)6"
     * "BK(3/2)7"
     * "BK(3/2)8"
     * "BK(3/2)9"
     * "BK(3/2)10"
     * "ReA2_Kd"
     * "ReA0_Kd"
     * "DeltaP_cu"
     * "Br_Kp_Ppenu"
     * "IB_Kl"
     * "IB_Kp"
     * "W_Kl"
     * "W_Kp"
     * "Br_Kp_Mupnu"
     * "Br_B_Xcenu"
     */
    virtual bool Init(const std::map<std::string, double>&);

    ////////////////////////////////////////////////////////////////////////

    Meson getMesons(const QCD::meson i) const 
    {
        return mesons[i];
    }

    Particle getQuarks(const QCD::quark i) const 
    {
        return quarks[i];
    }
    
    /**
     * @return @f$\alpha_s(Mz)@f$
     */
    double getAlsMz() const 
    {
        return AlsMz;
    }

    /**
     * set the initial condition @f$\alpha_s(Mz)@f$
     * @param AlsMz the initial condition @f$\alpha_s(Mz)@f$
     */
    void setAlsMz(double AlsMz) 
    {
        this->AlsMz = AlsMz;
    }

    /**
     * @return the scale Mz at which the initial condition for @f$\alpha_s(Mz)@f$ is given
     */
    double getMz() const 
    {
        return Mz;
    }

    /**
     * set the scale M at which the initial condition for @f$\alpha_s(M)@f$ is given
     * @param M the scale M in GeV
     */
    void setMz(double Mz) 
    {
        this->Mz = Mz;
    }

    /**
     * @return the number of colours
     */
    double getNc() const
    {
        return Nc;
    }

    /**
     * set the number of colours
     * @param Nc the number of colours
     */
    void setNc(double Nc)
    {
        this->Nc = Nc;
    }

    /**
     * @return the threshold between six- and five-flavour theory in GeV
     */
    double getMut() const 
    {
        return mut;
    }

    /**
     * set the threshold between six- and five-flavour theory
     * @param mut the threshold between six- and five-flavour theory in GeV
     */
    void setMut(double mut)
    {
        this->mut = mut;
    }

    /**
     * @return the threshold between five- and four-flavour theory in GeV
     */
    double getMub() const 
    {
        return mub;
    }

    /**
     * set the threshold between five- and four-flavour theory
     * @param mub the threshold between five- and four-flavour theory in GeV
     */
    void setMub(double mub) 
    {
        this->mub = mub;
    }

    /**
     * @return the threshold between four- and three-flavour theory in GeV
     */
    double getMuc() const 
    {
        return muc;
    }

    /**
     * set the threshold between four- and three-flavour theory
     * @param muc the threshold between four- and three-flavour theory in GeV
     */
    void setMuc(double muc) 
    {
        this->muc = muc;
    }

    /**
     * @return the pole mass of the top quark
     */
    double getMtpole() const 
    {
        return mtpole;
    }

    double getCF() const 
    {
        return CF;
    }
    
    /**
      * 
      * @return Expirimental value of the real part of the amplitude for the decay
      * K_L in two pions wthout ispspin change
      */
    double getReA0_kd() const{
        return ReA0_kd;
    }
    
    /**
      * 
      * @return Expirimental value of the real part of the amplitude for the decay
      * K_L in two pions with double isospin change
      */
    double getReA2_kd() const{
        return ReA2_kd;
    }
    
    /**
      * 
      * @return value of isospin breacking contribution for the decay
      * K_L in two pions
      */
    double getOmega_eta_etap() const {
        return Omega_eta_etap;
    }
    
    /**
     * 
     * @return the experimental value for the Br of the semileptonic decay of
     * the K+ in pion, electron and neutrino
     */
    double getBr_Kp_ppenu() const{
        return Br_Kp_Ppenu;
    }
    
    /**
     * 
     * @return the experimental value for the Br of the semileptonic decay of
     * the K+ in mu+ and neutrino
     */
    double getBr_Kp_Mupnu() const{
        return Br_Kp_Mupnu;
    }
    
    /**
     * 
     * @return the experimental value for the Br of the semileptonic decay of
     * the B -> Xc electron and neutrino
     */
    double getBr_B_Xcenu() const{
        return Br_B_Xcenu;
    }
    
    /**
     * 
     * @return isospin breaking corrections in relating Br(K_l -> pi0 nu nu) and
     * Br(K+ -> pi0 e nu)
     */
    double getIB_Kl() const{
        return IB_Kl;
    }
    
    /**
     * 
     * @return isospin breaking corrections in relating Br(K+ -> pi+ nu nu) and
     * Br(K+ -> pi0 e nu)
     */
    double getIB_Kp() const{
        return IB_Kp;
    }
    
    /**
     * 
     * @return long distance correction to the charm contribution to 
     * Br(K+ -> P+ nu nu), hep-ph/0503107 and hep-ph/0603079
     */
    double getDeltaP_cu() const{
        return DeltaP_cu;
    }

    BParameter getBBd() const 
    {
        return BBd;
    }

    BParameter getBBs() const 
    {
        return BBs;
    }
    
    BParameter getBD() const 
    {
        return BD;
    }
    
    BParameter getBK() const 
    {
        return BK;
    }
    
    BParameter getBKd1() const {
        return BKd1;
    }
    
    BParameter getBKd3() const {
        return BKd3;
    }
    
    /*BParameter getBD() const 
     {
        return BD;
    }*/
    
    ////////////////////////////////////////////////////////////////////////

    double Thresholds(int i) const;

    double AboveTh(double mu) const;

    double BelowTh(double mu) const;

    /**
     * the number of active flavour at scale @f$\mu@f$
     * @param mu the scale @f$\mu@f$ in GeV
     * @return the number of active flavour at scale @f$\mu@f$
     */
    double Nf(double mu) const;
    
    ////////////////////////////////////////////////////////////////////////

    /**
     * the @f$\beta_0@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\beta_0@f$ coefficient
     */
    double Beta0(double nf) const;

    /**
     * the @f$\beta_1@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\beta_1@f$ coefficient
     */
    double Beta1(double nf) const;

    /**
     * the @f$\beta_2@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\beta_2@f$ coefficient
     */
    double Beta2(double nf) const;

    /**
     * the strong running coupling @f$\alpha_s@f$ in the @f$\overline{\mathrm{MS}}@f$ scheme
     * @param mu the scale @f$\mu@f$ in GeV
     * @param logLambda log(Lambda)
     * @param nf the number of active flavours @f$n_f@f$
     * @param order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return @f$\alpha_s@f$
     */
    double AlsWithLambda(double mu, double logLambda, double nf, orders order) const;

    /**
     * the strong running coupling @f$\alpha_s@f$ in the @f$\overline{\mathrm{MS}}@f$ scheme
     * @param mu the scale @f$\mu@f$ in GeV
     * @param nf the number of active flavours
     * @param alsi the initial condition @f$\alpha_s(m_i)@f$
     * @param mi the scale @f$m_i@f$ in GeV
     * @param order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return @f$\alpha_s@f$
     */
    double Als(double mu, double nf, double alsi, double mi, orders order) const;

    /**
     * the strong running coupling @f$\alpha_s@f$ in the @f$\overline{\mathrm{MS}}@f$ scheme
     * @param mu the scale @f$\mu@f$ in GeV
     * @param nfmu the number of active flavours at the scale @f$\mu@f$
     * @param order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return @f$\alpha_s@f$
     */
    double Als(double mu, double nfmu, orders order) const;

    /**
     * the strong running coupling @f$\alpha_s@f$ in the @f$\overline{\mathrm{MS}}@f$ scheme
     * @param mu the scale @f$\mu@f$ in GeV
     * @param order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return @f$\alpha_s@f$
     */
    double Als(double mu, orders order = FULLNLO) const;

    /**
     * @f$\ln\Lambda_\mathrm{QCD}@f$ with five active flavours in GeV
     * @param order (=LO, FULLNLO, FULLNNLO)
     * @return @f$\ln\Lambda_\mathrm{QCD}@f$ for five active flavours
     */
    double logLambda5(orders order) const;

    double logLambda(double muMatching, double mf, double nfNEW, double nfORG, 
                     double logLambdaORG, orders order) const;

    double logLambda(double mu, orders order) const;  
    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * the @f$\gamma_0@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\gamma_0@f$ coefficient
     */
    double Gamma0(double nf) const;

    /**
     * the @f$\gamma_1@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\gamma_1@f$ coefficient
     */
    double Gamma1(double nf) const;

    /**
     * the @f$\gamma_2@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\gamma_2@f$ coefficient
     */
    double Gamma2(double nf) const;
    
    /**
     * @brief threshold corrections to the running mass
     * @param nf_f
     * @param nf_i
     * @return 
     */
    double threCorrForMass(double nf_f, double nf_i) const;
    
    /**
     * the running quark mass @f$m(\mu)@f$
     * @param mu the scale @f$\mu@f$ in GeV
     * @param m the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$
     * @param order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return the running quark mass @f$m(\mu)@f$
     */
    double Mrun(double mu, double m, orders order = FULLNLO) const;
    
    /**
     * runs the quark mass from @f$\mu_i@f$ to @f$\mu_f@f$
     * @param mu_f the scale @f$\mu_f@f$ in GeV
     * @param mu_i the scale @f$\mu_i@f$ in GeV
     * @param m the @f$\overline{\mathrm{MS}}@f$ mass @f$m(mu_i)@f$
     * @param order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return the running quark mass @f$m(\mu_f)@f$
     */
    double Mrun(double mu_f, double mu_i, double m, orders order = FULLNLO) const;

    /**
     * runs the quark mass from @f$\mu_i@f$ to @f$\mu_f@f$ at fixed nf
     * @param mu_f the scale @f$\mu_f@f$ in GeV
     * @param mu_i the scale @f$\mu_i@f$ in GeV
     * @param m the @f$\overline{\mathrm{MS}}@f$ mass @f$m(mu_i)@f$
     * @param nf the number of active flavours
     * @param order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return the running quark mass @f$m(\mu_f)@f$
     */
    double Mrun(double mu_f, double mu_i, double m, double nf, orders order = FULLNLO) const;

    ////////////////////////////////////////////////////////////////////////    
    
    /**
     * convert the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ to the pole mass
     * @param mbar the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     * @return the pole mass in GeV
     */
    double Mbar2Mp(double mbar) const;

    /**
     * convert the pole mass to the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$
     * @param mp the pole mass in GeV
     * @return the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     */
    double Mp2Mbar(double mp) const;

    double MS2DRqmass(const double& MSscale, const double& MSbar) const;
    
    /**
     * convert @f$\overline{\mathrm{MS}}@f$ to @f$\overline{\mathrm{DR}}@f$ quark masses
     * @param MSbar the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$
     * @return the @f$\overline{\mathrm{DR}}@f$ mass @f$m(m)@f$
     */
    double MS2DRqmass(const double& MSbar) const;
    
    ////////////////////////////////////////////////////////////////////////

protected:
    double Nc, CF, AlsMz, Mz, mut, mub, muc, mtpole;
    double ReA0_kd, ReA2_kd, Omega_eta_etap;
    double Br_Kp_Ppenu, IB_Kl, IB_Kp, DeltaP_cu, Br_Kp_Mupnu, Br_B_Xcenu;
    Particle quarks[6];
    Meson mesons[MESON_END];
    BParameter BBs, BBd, BD, BK, BKd1, BKd3;
    virtual void SetParameter(const std::string, const double&);
    bool computeYu, computeYd;

private:
    static const int CacheSize = 5;
    mutable double als_cache[8][CacheSize], logLambda5_cache[4][CacheSize], 
                   mp2mbar_cache[4][CacheSize], mrun_cache[10][CacheSize];
    bool computeFBd, computeBd, computemt;
    double BBsoBBd, FBsoFBd;
    double ZeroNf5(double *x, double *) const;
    double Mp2Mbara(double * mu, double * mp) const;
    void CacheShift(double cache[][5], int n) const;
};

#endif	/* QCD_H */
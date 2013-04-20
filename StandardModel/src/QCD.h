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

    static const int NQCDvars = 76;

    /**
     * array containing the labels under which all QCD parameters must be
     * stored in a Parameters object
     */
    static const std::string QCDvars[NQCDvars];

    QCD();

    virtual std::string ModelName() const 
    {
        return "QCD";
    }
    
    std::string orderToString(const orders order) const;
    
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
      * @return Expirimental value of the real part of the amplitude for the decay
      * K_L in two pions wthout ispspin change
      */
    double getReA0_kd() const
    {
        return ReA0_kd;
    }
    
    /**
      * @return Expirimental value of the real part of the amplitude for the decay
      * K_L in two pions with double isospin change
      */
    double getReA2_kd() const
    {
        return ReA2_kd;
    }
    
    /**
      * @return value of isospin breacking contribution for the decay
      * K_L in two pions
      */
    double getOmega_eta_etap() const 
    {
        return Omega_eta_etap;
    }
    
    /**
     * @return the experimental value for the Br of the semileptonic decay of
     * the K+ in pion, electron and neutrino
     */
    double getBr_Kp_ppenu() const
    {
        return Br_Kp_Ppenu;
    }
    
    /**
     * @return the experimental value for the Br of the semileptonic decay of
     * the K+ in mu+ and neutrino
     */
    double getBr_Kp_Mupnu() const
    {
        return Br_Kp_Mupnu;
    }
    
    /**
     * @return the experimental value for the Br of the semileptonic decay of
     * the B -> Xc electron and neutrino
     */
    double getBr_B_Xcenu() const
    {
        return Br_B_Xcenu;
    }
    
    /**
     * @return isospin breaking corrections in relating Br(K_l -> pi0 nu nu) and
     * Br(K+ -> pi0 e nu)
     */
    double getIB_Kl() const
    {
        return IB_Kl;
    }
    
    /**
     * @return isospin breaking corrections in relating Br(K+ -> pi+ nu nu) and
     * Br(K+ -> pi0 e nu)
     */
    double getIB_Kp() const
    {
        return IB_Kp;
    }
    
    /**
     * @return long distance correction to the charm contribution to 
     * Br(K+ -> P+ nu nu), hep-ph/0503107 and hep-ph/0603079
     */
    double getDeltaP_cu() const
    {
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
    
    BParameter getBKd1() const 
    {
        return BKd1;
    }
    
    BParameter getBKd3() const 
    {
        return BKd3;
    }
    
    ////////////////////////////////////////////////////////////////////////

    double Thresholds(const int i) const;

    double AboveTh(const double mu) const;

    double BelowTh(const double mu) const;

    /**
     * the number of active flavour at scale @f$\mu@f$
     * @param mu the scale @f$\mu@f$ in GeV
     * @return the number of active flavour at scale @f$\mu@f$
     */
    double Nf(const double mu) const;
    
    ////////////////////////////////////////////////////////////////////////

    /**
     * the @f$\beta_0@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\beta_0@f$ coefficient
     */
    double Beta0(const double nf) const;

    /**
     * the @f$\beta_1@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\beta_1@f$ coefficient
     */
    double Beta1(const double nf) const;

    /**
     * the @f$\beta_2@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\beta_2@f$ coefficient
     */
    double Beta2(const double nf) const;

    /**
     * the strong running coupling @f$\alpha_s@f$ in the @f$\overline{\mathrm{MS}}@f$ scheme
     * @param mu the scale @f$\mu@f$ in GeV
     * @param alsi the initial condition @f$\alpha_s(m_i)@f$
     * @param mu_i the scale @f$m_i@f$ in GeV
     * @param order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return @f$\alpha_s@f$
     */
    double AlsWithInit(const double mu, const double alsi, const double mu_i, 
                       const orders order) const;    
    
    double AlsWithLambda(const double mu, const orders order) const;

    /**
     * the strong running coupling @f$\alpha_s@f$ in the @f$\overline{\mathrm{MS}}@f$ scheme
     * @param mu the scale @f$\mu@f$ in GeV
     * @param order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return @f$\alpha_s@f$
     */
    double Als(const double mu, const orders order = FULLNNLO) const;

    /**
     * @f$\ln\Lambda_\mathrm{QCD}@f$ with nf flavours in GeV
     * @param The number of active flavours. 
     * @param order (=LO, FULLNLO, FULLNNLO)
     * @return @f$\ln\Lambda_\mathrm{QCD}@f$ for five active flavours
     */
    double logLambda(const double nf, orders order) const;  
    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * the @f$\gamma_0@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\gamma_0@f$ coefficient
     */
    double Gamma0(const double nf) const;

    /**
     * the @f$\gamma_1@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\gamma_1@f$ coefficient
     */
    double Gamma1(const double nf) const;

    /**
     * the @f$\gamma_2@f$ coefficient
     * @param nf the number of active flavours
     * @return the @f$\gamma_2@f$ coefficient
     */
    double Gamma2(const double nf) const;
    
    /**
     * the running quark mass @f$m(\mu)@f$
     * @param mu the scale @f$\mu@f$ in GeV
     * @param m the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$
     * @param order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return the running quark mass @f$m(\mu)@f$
     */
    double Mrun(const double mu, const double m, const orders order = FULLNLO) const;
    
    /**
     * runs the quark mass from @f$\mu_i@f$ to @f$\mu_f@f$
     * @param mu_f the scale @f$\mu_f@f$ in GeV
     * @param mu_i the scale @f$\mu_i@f$ in GeV
     * @param m the @f$\overline{\mathrm{MS}}@f$ mass @f$m(mu_i)@f$
     * @param order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return the running quark mass @f$m(\mu_f)@f$
     */
    double Mrun(const double mu_f, const double mu_i, const double m, 
                const orders order = FULLNNLO) const;

    ////////////////////////////////////////////////////////////////////////    
    
    /**
     * convert the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ to the pole mass
     * @param mbar the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     * @param order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return the pole mass in GeV
     */
    double Mbar2Mp(const double mbar, const orders order = FULLNNLO) const;

    /**
     * convert the pole mass to the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$
     * @param mp the pole mass in GeV
     * @param order (=LO, NLO, NNLO, FULLNLO, FULLNNLO)
     * @return the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     */
    double Mp2Mbar(const double mp, const orders order = FULLNNLO) const;

    double MS2DRqmass(const double MSscale, const double MSbar) const;
    
    /**
     * convert @f$\overline{\mathrm{MS}}@f$ to @f$\overline{\mathrm{DR}}@f$ quark masses
     * @param MSbar the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$
     * @return the @f$\overline{\mathrm{DR}}@f$ mass @f$m(m)@f$
     */
    double MS2DRqmass(const double MSbar) const;
    
    ////////////////////////////////////////////////////////////////////////

protected:
    double Nc, CF, mtpole;
    Particle quarks[6];
    Meson mesons[MESON_END];
    bool computeYu, computeYd;

    // model parameters
    double AlsMz, Mz, mut, mub, muc;
    double ReA0_kd, ReA2_kd, Omega_eta_etap;
    double Br_Kp_Ppenu, IB_Kl, IB_Kp, DeltaP_cu, Br_Kp_Mupnu, Br_B_Xcenu;
    double BBsoBBd, FBsoFBd;
    BParameter BBs, BBd, BD, BK, BKd1, BKd3;
    virtual void SetParameter(const std::string, const double&);

private:
    double zeta2, zeta3;
    bool computeFBd, computeBd, computemt;
    double AlsWithLambda(const double mu, const double logLambda, 
                         const orders order) const;
    double ZeroNf5(double *x, double *) const;
    double logLambda5(orders order) const;
    double logLambda(const double muMatching, const double mf, 
                     const double nfNEW, const double nfORG, 
                     const double logLambdaORG, orders order) const;
    double threCorrForMass(const double nf_f, const double nf_i) const;
    double MrunTMP(const double mu_f, const double mu_i, const double m, 
                   const orders order = FULLNNLO) const;
    double Mp2MbarTMP(double *mu, double *params) const;

    // caches
    static const int CacheSize = 5;
    mutable double als_cache[8][CacheSize], logLambda5_cache[4][CacheSize], 
                   mrun_cache[10][CacheSize], mp2mbar_cache[5][CacheSize];
    void CacheShift(double cache[][5], int n) const;
};

#endif	/* QCD_H */

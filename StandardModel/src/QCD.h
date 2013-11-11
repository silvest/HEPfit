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
 * @details This class is .... 
 * 
 * Model parameters: 
 * \li \b AlsMz:&nbsp; the strong coupling constant at the Z-boson mass, @f$\alpha_s(M_Z)@f$,
 * \li \b Mz:&nbsp; the mass of the @f$Z@f$ boson in GeV,
 * \li \b mup:&nbsp; the @f$\overline{\mathrm{MS}}@f$ mass of the up quark at 2 GeV, @f$m_u(2\,\mathrm{GeV})@f$, in GeV,
 * \li \b mdown:&nbsp; the @f$\overline{\mathrm{MS}}@f$ mass of the down quark at 2 GeV, @f$m_d(2\,\mathrm{GeV})@f$, in GeV,
 * \li \b mcharm:&nbsp; the @f$\overline{\mathrm{MS}}@f$ scale-invariant mass of the charm quark, @f$m_c(m_c)@f$, in GeV,
 * \li \b mstrange:&nbsp; the @f$\overline{\mathrm{MS}}@f$ mass of the strange quark at 2 GeV, @f$m_s(2\,\mathrm{GeV})@f$, in GeV,
 * \li \b mtop:&nbsp; the pole mass of the top quark in GeV,
 * \li \b mbottom:&nbsp; the @f$\overline{\mathrm{MS}}@f$ scale-invariant mass of the bottom quark, @f$m_b(m_b)@f$, in GeV,
 * \li \b mut:&nbsp; the threshold between six- and five-flavour theory in GeV,
 * \li \b mub:&nbsp; the threshold between five- and four-flavour theory in GeV, 
 * \li \b muc:&nbsp; the threshold between four- and three-flavour theory in GeV, 
 * \li \b MBd:&nbsp; the mass of the \f$ B_d \f$ meson in \f$ GeV \f$,
 * \li \b tBd:&nbsp; the lifetime of the \f$ B_d \f$ meson in \f$ ps^{-1} \f$,
 * \li \b MBs:&nbsp; the mass of the \f$ B_s \f$ meson in \f$ GeV \f$,
 * \li \b tBs:&nbsp; the lifetime of the \f$ B_s \f$ meson in \f$ ps^{-1} \f$,
 * \li \b MBp:&nbsp; the mass of the \f$ B^\pm \f$ meson in \f$ GeV \f$,
 * \li \b MK0:&nbsp; the mass of the \f$ K^0 \f$ meson in \f$ GeV \f$,
 * \li \b MKp:&nbsp; the mass of the \f$ K^\pm \f$ meson in \f$ GeV \f$,
 * \li \b MD:&nbsp; the mass of the \f$ D^0 \f$ meson in \f$ GeV \f$,
 * \li \b tKl:&nbsp; the lifetime of the \f$ K_L \f$ meson in \f$ ps^{-1} \f$,
 * \li \b tKp:&nbsp; the lifetime of the \f$ K^\pm \f$ meson in \f$ ps^{-1} \f$,
 * \li \b FBs:&nbsp; the decay constant of the \f$ B_s \f$ meson in \f$ GeV \f$,
 * \li \b FBsoFBd:&nbsp; the ratio \f$ F_{B_d}/F_{B_s} \f$ in \f$ GeV \f$ necessary to compute \f$ F_{B_s} \f$,
 * \li \b FD:&nbsp; the decay constant of the \f$ D^0 \f$ meson in \f$ GeV \f$,
 * \li \b BBsoBBd:&nbsp;
 * \li \b BBs1:&nbsp; the bag parameter for \f$ O_1 \f$ in \f$ \Delta b = 2 \f$ processes in \f$ B_s \f$,
 * \li \b BBs2:&nbsp; the bag parameter for \f$ O_2 \f$ in \f$ \Delta b = 2 \f$ processes in \f$ B_s \f$,
 * \li \b BBs3:&nbsp; the bag parameter for \f$ O_3 \f$ in \f$ \Delta b = 2 \f$ processes in \f$ B_s \f$,
 * \li \b BBs4:&nbsp; the bag parameter for \f$ O_4 \f$ in \f$ \Delta b = 2 \f$ processes in \f$ B_s \f$,
 * \li \b BBs5:&nbsp; the bag parameter for \f$ O_5 \f$ in \f$ \Delta b = 2 \f$ processes in \f$ B_s \f$,
 * \li \b BBsscale:&nbsp; the scale at which the bag parameters are specified for the \f$ B_s \f$ system,
 * \li \b BBsscheme:&nbsp; the scheme in which the bag parameters are specified for the \f$ B_s \f$ system,
 * \li \b BD1:&nbsp; the bag parameter for \f$ O_1 \f$ in \f$ \Delta c = 2 \f$ processes in \f$ D^0 \f$,
 * \li \b BD2:&nbsp; the bag parameter for \f$ O_2 \f$ in \f$ \Delta c = 2 \f$ processes in \f$ D^0 \f$,
 * \li \b BD3:&nbsp; the bag parameter for \f$ O_3 \f$ in \f$ \Delta c = 2 \f$ processes in \f$ D^0 \f$,
 * \li \b BD4:&nbsp; the bag parameter for \f$ O_4 \f$ in \f$ \Delta c = 2 \f$ processes in \f$ D^0 \f$,
 * \li \b BD5:&nbsp; the bag parameter for \f$ O_5 \f$ in \f$ \Delta c = 2 \f$ processes in \f$ D^0 \f$,
 * \li \b BDscale:&nbsp; the scale at which the bag parameters are specified for the \f$ D_0 \f$ system,
 * \li \b BDscheme:&nbsp; the scheme in which the bag parameters are specified for the \f$ D_0 \f$ system,
 * \li \b BK1:&nbsp; the bag parameter for \f$ O_1 \f$ in \f$ \Delta c = 2 \f$ processes in \f$ K^0 \f$,
 * \li \b BK2:&nbsp; the bag parameter for \f$ O_2 \f$ in \f$ \Delta c = 2 \f$ processes in \f$ K^0 \f$,
 * \li \b BK3:&nbsp; the bag parameter for \f$ O_3 \f$ in \f$ \Delta c = 2 \f$ processes in \f$ K^0 \f$,
 * \li \b BK4:&nbsp; the bag parameter for \f$ O_4 \f$ in \f$ \Delta c = 2 \f$ processes in \f$ K^0 \f$,
 * \li \b BK5:&nbsp; the bag parameter for \f$ O_5 \f$ in \f$ \Delta c = 2 \f$ processes in \f$ K^0 \f$,
 * \li \b BKscale:&nbsp; the scale at which the bag parameters are specified for the \f$ K^0 \f$ system,
 * \li \b BKscheme:&nbsp; the scheme in which the bag parameters are specified for the \f$ K^0 \f$ system,
 * \li \b BK(1/2)1:&nbsp;
 * \li \b BK(1/2)2:&nbsp;
 * \li \b BK(1/2)3:&nbsp;
 * \li \b BK(1/2)4:&nbsp;
 * \li \b BK(1/2)5:&nbsp;
 * \li \b BK(1/2)6:&nbsp;
 * \li \b BK(1/2)7:&nbsp;
 * \li \b BK(1/2)8:&nbsp;
 * \li \b BK(1/2)9:&nbsp;
 * \li \b BK(1/2)10:&nbsp;
 * \li \b BKd_scale:&nbsp;
 * \li \b BKd_scheme:&nbsp;
 * \li \b BK(3/2)1:&nbsp;
 * \li \b BK(3/2)2:&nbsp;
 * \li \b BK(3/2)3:&nbsp;
 * \li \b BK(3/2)4:&nbsp;
 * \li \b BK(3/2)5:&nbsp;
 * \li \b BK(3/2)6:&nbsp;
 * \li \b BK(3/2)7:&nbsp;
 * \li \b BK(3/2)8:&nbsp;
 * \li \b BK(3/2)9:&nbsp;
 * \li \b BK(3/2)10:&nbsp;
 * \li \b ReA0_Kd:&nbsp; the experimental value of the real part of the amplitude for @f$K^0\to\pi\pi@f$ with @f$\Delta I=0@f$,
 * \li \b ReA2_Kd:&nbsp; the experimental value of the real part of the amplitude for @f$K^0\to\pi\pi@f$ with @f$\Delta I=2@f$, 
 * \li \b Omega_eta_etap:&nbsp; the isospin breaking contribution in @f$K^0\to\pi\pi@f$,
 * \li \b Br_Kp_P0enu:&nbsp; the experimental value for the branching ratio of @f$K^+\to\pi^0e^+\nu@f$,
 * \li \b Br_Kp_munu:&nbsp; the experimental value for the branching ratio of @f$K^+\to\mu^+\nu@f$,
 * \li \b Br_B_Xcenu:&nbsp; the experimental value for the branching ratio of @f$B\to X_c e\nu@f$,
 * \li \b DeltaP_cu:&nbsp; the long-distance correction to the charm contribution of @f$K^+\to\pi^+\nu\bar{\nu}@f$,
 * \li \b IB_Kl:&nbsp; the isospin breaking corrections between @f$K_L\to\pi^0\nu\bar{\nu}@f$ and @f$K^+\to\pi^0 e^+\nu@f$,
 * \li \b IB_Kp:&nbsp; the isospin breaking corrections between @f$K^+\to\pi^+ \nu\bar{\nu}@f$ and @f$K^+\to\pi^0 e^+\nu@f$. 
 * 
 */
class QCD: public Model {
public:

    /**
     * An enum type for mesons. 
     */
    enum meson 
    {
        B_D, /**< @f$B_d@f$ meson */
        B_S, /**< @f$B_s@f$ meson */
        B_P, /**< @f$B^\pm@f$ meson */
        K_0, /**< @f$K^0@f$ meson */
        K_P, /**< @f$K^\pm@f$ meson */
        D_0, /**< @f$D^0@f$ meson */
        P_0, /**< @f$\pi^0@f$ meson */
        P_P, /**< @f$\pi^\pm@f$ meson */
        MESON_END /**< The size of this enum. */
    }; 

    /**
     * An enum type for quarks. 
     */
    enum quark 
    {
        UP, /**< Up quark */
        DOWN, /**< Down quark */
        CHARM, /**< Charm quark */
        STRANGE, /**< Strange quark */
        TOP, /**< Top quark */
        BOTTOM /**< Bottom quark */
    };

    static const int NQCDvars = 78; /**< The number of model parameters in QCD. */

    /**
     * An array containing the labels under which all QCD parameters are stored 
     * in a vector of ModelParameter via InputParser::ReadParameters(). 
     */
    static const std::string QCDvars[NQCDvars];

    /**
     * The default constructor. 
     */
    QCD();

    /**
     * @return The name of the model defined in the current class. 
     */
    virtual std::string ModelName() const 
    {
        return "QCD";
    }
    
    /**
     * Converts an object of the enum type "orders" to the corresponding string. 
     * @param[in] order An object of the enum type "orders". 
     * @return The string of the given "order". 
     */
    std::string orderToString(const orders order) const;
    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Sets flags for QCD. 
     * @param[in] name A name of a flag. 
     * @param[in] value A value of the given flag. 
     * @return A boolean value indicating whether the given flag name is associated 
     * with QCD. 
     */
    virtual bool SetFlag(const std::string name, const bool& value);
    
    /**
     * Pre update. 
     * @return 
     */
    virtual bool PreUpdate();

    /**
     * Post update. 
     * @return 
     */
    virtual bool PostUpdate();      
    
    /**
     * Updates the QCD parameters found in the argument. 
     * @param[in] DPars A map containing the parameters (all as double) to be updated. 
     */
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    /**
     * Checks that all required parameters are present in a given map. 
     * @param[in] DPars A map containing the parameters (all as double) to be used in Monte Carlo. 
     */
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    /**
     * Initializes the QCD parameters found in the argument. 
     * @param[in] DPars A map containing the parameters (all as double) to be used in Monte Carlo. 
     */
    virtual bool Init(const std::map<std::string, double>& DPars);

    ////////////////////////////////////////////////////////////////////////
    // get and set methods for class members

    /**
     * @param[in] m The name of a meson. 
     * @return The object of the meson found in the argument. 
     */
    Meson getMesons(const meson m) const 
    {
        return mesons[m];
    }

    /**
     * @param[in] q The name of a quark. 
     * @return The object of the quark found in the argument. 
     */
    Particle getQuarks(const quark q) const 
    {
        return quarks[q];
    }
    
    /**
     * @return The strong coupling constant at @f$M_Z@f$, @f$\alpha_s(M_Z)@f$. 
     */
    double getAlsMz() const 
    {
        return AlsMz;
    }

    /**
     * Sets the strong coupling constant at @f$M_Z@f$, @f$\alpha_s(M_Z)@f$. 
     * @param[in] AlsMz @f$\alpha_s(M_Z)@f$. 
     */
    void setAlsMz(double AlsMz) 
    {
        this->AlsMz = AlsMz;
    }

    /**
     * @return The @f$Z@f$-boson mass @f$M_Z@f$. 
     */
    double getMz() const 
    {
        return Mz;
    }

    /**
     * Sets the @f$Z@f$-boson mass @f$M_Z@f$. 
     * @param[in] Mz @f$M_Z@f$ in GeV. 
     */
    void setMz(double Mz) 
    {
        this->Mz = Mz;
    }

    /**
     * @return The number of colours. 
     */
    double getNc() const
    {
        return Nc;
    }

    /**
     * Sets the number of colours. 
     * @param[in] Nc The number of colours. 
     */
    void setNc(double Nc)
    {
        this->Nc = Nc;
    }

    /**
     * @return The threshold between six- and five-flavour theory in GeV. 
     */
    double getMut() const 
    {
        return mut;
    }

    /**
     * Sets the threshold between six- and five-flavour theory. 
     * @param[in] mut The threshold between six- and five-flavour theory in GeV. 
     */
    void setMut(double mut)
    {
        this->mut = mut;
    }

    /**
     * @return The threshold between five- and four-flavour theory in GeV. 
     */
    double getMub() const 
    {
        return mub;
    }

    /**
     * Sets the threshold between five- and four-flavour theory. 
     * @param[in] mub The threshold between five- and four-flavour theory in GeV. 
     */
    void setMub(double mub) 
    {
        this->mub = mub;
    }

    /**
     * @return The threshold between four- and three-flavour theory in GeV. 
     */
    double getMuc() const 
    {
        return muc;
    }

    /**
     * Set the threshold between four- and three-flavour theory. 
     * @param[in] muc The threshold between four- and three-flavour theory in GeV. 
     */
    void setMuc(double muc) 
    {
        this->muc = muc;
    }

    /**
     * @return The pole mass of the top quark. 
     */
    double getMtpole() const 
    {
        return mtpole;
    }

    /**
     * @return The Casimir factor of QCD. 
     */
    double getCF() const 
    {
        return CF;
    }
    
    /**
     * @return 
     */
    BParameter getBBd() const 
    {
        return BBd;
    }

    /**
     * @return 
     */
    BParameter getBBs() const 
    {
        return BBs;
    }
    
    /**
     * @return 
     */
    BParameter getBD() const 
    {
        return BD;
    }
    
    /**
     * @return 
     */
    BParameter getBK() const 
    {
        return BK;
    }
    
    /**
     * @return 
     */
    BParameter getBKd1() const 
    {
        return BKd1;
    }
    
    /**
     * @return 
     */
    BParameter getBKd3() const 
    {
        return BKd3;
    }
    
    /**
      * @return The experimental value of the real part of the amplitude for 
      * @f$K^0\to\pi\pi@f$ with @f$\Delta I=0@f$. 
      */
    double getReA0_Kd() const
    {
        return ReA0_Kd;
    }
    
    /**
      * @return The experimental value of the real part of the amplitude for 
      * @f$K^0\to\pi\pi@f$ with @f$\Delta I=2@f$.
      */
    double getReA2_Kd() const
    {
        return ReA2_Kd;
    }
    
    /**
      * @return The isospin breaking contribution in @f$K^0\to\pi\pi@f$. 
      */
    double getOmega_eta_etap() const 
    {
        return Omega_eta_etap;
    }
    
    /**
     * @return The experimental value for the branching ratio of @f$K^+\to\pi^0e^+\nu@f$. 
     */
    double getBr_Kp_P0enu() const
    {
        return Br_Kp_P0enu;
    }
    
    /**
     * @return The experimental value for the branching ratio of @f$K^+\to\mu^+\nu@f$. 
     */
    double getBr_Kp_munu() const
    {
        return Br_Kp_munu;
    }
    
    /**
     * @return The experimental value for the branching ratio of @f$B\to X_c e\nu@f$. 
     */
    double getBr_B_Xcenu() const
    {
        return Br_B_Xcenu;
    }
    
    /**
     * @return The long-distance correction to the charm contribution of @f$K^+\to\pi^+\nu\bar{\nu}@f$. 
     * 
     * References: 
     * [<A HREF="http://inspirehep.net/record/678222?ln=en">Isidori et al.(2005)</A>], 
     * [<A HREF="http://inspirehep.net/record/712083?ln=en">Buras et al.(2006)</A>]
     */
    double getDeltaP_cu() const
    {
        return DeltaP_cu;
    }

    /**
     * @return The isospin breaking corrections between 
     * @f$K_L\to\pi^0\nu\bar{\nu}@f$ and @f$K^+\to\pi^0 e^+\nu@f$. 
     */
    double getIB_Kl() const
    {
        return IB_Kl;
    }
    
    /**
     * @return The isospin breaking corrections between  
     * @f$K^+\to\pi^+ \nu\bar{\nu}@f$ and @f$K^+\to\pi^0 e^+\nu@f$. 
     */
    double getIB_Kp() const
    {
        return IB_Kp;
    }
    
    ////////////////////////////////////////////////////////////////////////

    /**
     * 
     * @param[in] i
     * @return 
     */
    double Thresholds(const int i) const;

    /**
     * 
     * @param[in] mu
     * @return 
     */
    double AboveTh(const double mu) const;

    /**
     * 
     * @param[in] mu
     * @return 
     */
    double BelowTh(const double mu) const;

    /**
     * @param[in] mu A scale @f$\mu@f$ in GeV. 
     * @return The number of active flavour at scale @f$\mu@f$. 
     */
    double Nf(const double mu) const;
    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @param[in] nf The number of active flavours. 
     * @return The @f$\beta_0@f$ coefficient. 
     */
    double Beta0(const double nf) const;

    /**
     * @param[in] nf The number of active flavours. 
     * @return The @f$\beta_1@f$ coefficient. 
     */
    double Beta1(const double nf) const;

    /**
     * @param[in] nf The number of active flavours. 
     * @return The @f$\beta_2@f$ coefficient. 
     */
    double Beta2(const double nf) const;

    /**
     * Computes the running strong coupling @f$\alpha_s(\mu)@f$ from @f$\alpha_s(\mu_i)@f$ 
     * in the @f$\overline{\mathrm{MS}}@f$ scheme, where it is forbidden to across 
     * a flavour threshould in the RG running from @f$\mu_i@f$ to @f$\mu@f$. 
     * @param[in] mu A scale @f$\mu@f$ in GeV. 
     * @param[in] alsi An initial condition for the coupling at the scale given below. 
     * @param[in] mu_i An initial scale @f$\mu_i@f$ in GeV. 
     * @param[in] order LO, NLO or FULLNLO in the @f$\alpha_s@f$ expansion. 
     * @return The strong coupling constant @f$\alpha_s(\mu)@f$ in the 
     * @f$\overline{\mathrm{MS}}@f$ scheme. 
     */
    double AlsWithInit(const double mu, const double alsi, const double mu_i, 
                       const orders order) const;    
    
    /**
     * Computes the running strong coupling @f$\alpha_s(\mu)@f$ in the 
     * @f$\overline{\mathrm{MS}}@f$ scheme with the use of @f$\Lambda_{\rm QCD}@f$. 
     * @param[in] mu A scale @f$\mu@f$ in GeV. 
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion. 
     * @return The strong coupling constant @f$\alpha_s(\mu)@f$ in the 
     * @f$\overline{\mathrm{MS}}@f$ scheme. 
     */
    double AlsWithLambda(const double mu, const orders order) const;

    /**
     * Computes the running strong coupling @f$\alpha_s(\mu)@f$ in the 
     * @f$\overline{\mathrm{MS}}@f$ scheme. In the cases of LO, NLO and FULLNNLO, 
     * the coupling is computed with AlsWithInit(). On the other hand, in the 
     * cases of NNLO and FULLNNLO, the coupling is computed with AlsWithLambda(). 
     * @param[in] mu A scale @f$\mu@f$ in GeV. 
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion. 
     * @return The strong coupling constant @f$\alpha_s(\mu)@f$ in the 
     * @f$\overline{\mathrm{MS}}@f$ scheme. 
     */
    double Als(const double mu, const orders order = FULLNLO) const;

    /**
     * Computes @f$\ln\Lambda_\mathrm{QCD}@f$ with nf flavours in GeV. 
     * @param[in] nf The number of active flavours. 
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion. 
     * @return @f$\ln\Lambda_\mathrm{QCD}@f$ with nf flavours in GeV. 
     */
    double logLambda(const double nf, orders order) const;  
    
    /**
     * temporary function waiting for the implementation of NNLO etact
     * @param mu
     * @return 
     */
    
    double Als4(const double mu) const;
    
    /**
     * temporary function waiting for the implementation of NNLO etact
     * @param mu_f
     * @param mu_i
     * @param m
     * @return 
     */
    double Mrun4(const double mu_f, const double mu_i, const double m) const;
    
    ////////////////////////////////////////////////////////////////////////
    
    /**
     * @param[in] nf The number of active flavours. 
     * @return The @f$\gamma_0@f$ coefficient. 
     */
    double Gamma0(const double nf) const;

    /**
     * @param[in] nf The number of active flavours. 
     * @return The @f$\gamma_1@f$ coefficient. 
     */
    double Gamma1(const double nf) const;

    /**
     * @param[in] nf The number of active flavours. 
     * @return The @f$\gamma_2@f$ coefficient. 
     */
    double Gamma2(const double nf) const;
    
    /**
     * Computes a running quark mass @f$m(\mu)@f$ from @f$m(m)@f$. 
     * @param[in] mu A scale @f$\mu@f$ in GeV
     * @param[in] m The @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV. 
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion. 
     * @return The running quark mass @f$m(\mu)@f$ in GeV. 
     */
    double Mrun(const double mu, const double m, const orders order = FULLNLO) const;
    
    /**
     * Runs a quark mass from @f$\mu_i@f$ to @f$\mu_f@f$. 
     * @param[in] mu_f A scale @f$\mu_f@f$ in GeV. 
     * @param[in] mu_i A scale @f$\mu_i@f$ in GeV. 
     * @param[in] m The @f$\overline{\mathrm{MS}}@f$ mass @f$m(\mu_i)@f$ in GeV. 
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion. 
     * @return The running quark mass @f$m(\mu_f)@f$ in GeV. 
     */
    double Mrun(const double mu_f, const double mu_i, const double m, 
                const orders order = FULLNLO) const;

    ////////////////////////////////////////////////////////////////////////    
    
    /**
     * Converts the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ to the pole mass
     * @param[in] mbar the @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion. 
     * @return The pole mass in GeV
     */
    double Mbar2Mp(const double mbar, const orders order = FULLNLO) const;

    /**
     * Converts a quark pole mass to the corresponding @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$. 
     * @param[in] mp The pole mass of the bottom or top quark in GeV. 
     * @param[in] order LO, NLO, FULLNLO, NNLO or FULLNNLO in the @f$\alpha_s@f$ expansion. 
     * @return The @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV. 
     */
    double Mp2Mbar(const double mp, const orders order = FULLNLO) const;

    /**
     * 
     * @param[in] MSscale
     * @param[in] MSbar
     * @return 
     */
    double MS2DRqmass(const double MSscale, const double MSbar) const;
    
    /**
     * Converts a quark mass from the @f$\overline{\mathrm{MS}}@f$ scheme to 
     * the @f$\overline{\mathrm{DR}}@f$ scheme. 
     * @param[in] MSbar The @f$\overline{\mathrm{MS}}@f$ mass @f$m(m)@f$ in GeV. 
     * @return The @f$\overline{\mathrm{DR}}@f$ mass @f$m(m)@f$ in GeV. 
     */
    double MS2DRqmass(const double MSbar) const;
    
    ////////////////////////////////////////////////////////////////////////

protected:
    double Nc; /**< The number of colours. */
    double CF; /**< The Casimir factor in the SU(Nc) gauge theory. */
    double mtpole;  /**< The pole mass of the top quark. */
    Particle quarks[6];
    Meson mesons[MESON_END];
    bool computeYu, computeYd;
    BParameter BBs, BBd, BD, BK, BKd1, BKd3;

    // model parameters
    double AlsMz, Mz, mut, mub, muc;
    double ReA0_Kd, ReA2_Kd, Omega_eta_etap;
    double Br_Kp_P0enu, IB_Kl, IB_Kp, DeltaP_cu, Br_Kp_munu, Br_B_Xcenu;
    double BBsoBBd, FBsoFBd;
    virtual void SetParameter(const std::string, const double&);

private:
    double zeta2, zeta3;
    bool computeFBd, computeBd, computemt;
    double AlsWithLambda(const double mu, const double logLambda, 
                         const orders order) const;
    double ZeroNf6NLO(double *logLambda6, double *logLambda5_in) const;
    double ZeroNf5(double *logLambda5, double *order) const;
    double ZeroNf4NLO(double *logLambda4, double *logLambda5_in) const;
    double ZeroNf3NLO(double *logLambda3, double *logLambda4_in) const;
    double logLambda5(orders order) const;
    double logLambdaNLO(const double nfNEW, const double nfORG,
                        const double logLambdaORG) const;
    double logLambda(const double muMatching, const double mf, 
                     const double nfNEW, const double nfORG, 
                     const double logLambdaORG, orders order) const;
    double threCorrForMass(const double nf_f, const double nf_i) const;
    double MrunTMP(const double mu_f, const double mu_i, const double m, 
                   const orders order) const;
    double Mp2MbarTMP(double *mu, double *params) const;

    // caches
    static const int CacheSize = 5;
    mutable double als_cache[8][CacheSize], logLambda5_cache[4][CacheSize],
                   logLambdaNLO_cache[9][CacheSize],
                   mrun_cache[10][CacheSize], mp2mbar_cache[5][CacheSize];
    void CacheShift(double cache[][5], int n) const;
};

#endif	/* QCD_H */

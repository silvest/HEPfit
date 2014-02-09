/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMCACHE_H
#define	EWSMCACHE_H

#include <cmath>
#include <PVfunctions.h>
#include <Polylogarithms.h>
#include <ClausenFunctions.h>
#include "StandardModel.h"

using namespace gslpp;

/**
 * @class EWSMcache
 * @ingroup StandardModel
 * @brief A class for cache variables used in computing radiative corrections
 * to the %EW precision observables.
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details This class provides caching methods for a bunch of the functions
 * appearing in EWSMOneLoopEW, EWSMTwoLoopQCD, EWSMTwoLoopEW, EWSMThreeLoopQCD,
 * EWSMThreeLoopEW2QCD and EWSMThreeLoopEW classes. Each cache method calls the
 * private member function CacheCheck() and newCacheForDouble()
 * (or newCacheForComplex()).
 * 
 * Moreover, this class contains methods to access model parameters and functions
 * in QCD and StandardModel class, such that the classes listed above do not call
 * directly the functions in QCD nor StandardModel.
 *
 * The internal flags #FlagDebug and #FlagCacheInEWSMcache, which can be changed
 * with setFlagDebug() and setFlagCacheInEWSMcache(), respectively, are designed 
 * for debugging.
 */
class EWSMcache {

public:

    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    EWSMcache(const StandardModel& SM_i);

    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief A set method to change the internal boolean flag #FlagDebug.
     * @details The flag #FlagDebug=true is used for testing and debugging the
     * codes associated with the current class. The flag #FlagDebug
     * is set to false in the constructor EWSMcache() by default.
     * @param[in] FlagDebug a boolean flag for debugging
     */
    void setFlagDebug(bool FlagDebug)
    {
        this->FlagDebug = FlagDebug;
    }

    /**
     * @brief A set method to change the model flag CacheInEWSMcache in
     * StandardModel.
     * @details Setting CacheInEWSMcache to false, the caching methods
     * defined in the current class are not employed in numerical computations.
     * The flag is set to true in the constructor EWSMcache() by default.
     * 
     * See also @ref StandardModelFlags "the description of the StandardModel flags".
     * @param[in] FlagCacheInEWSMcache a boolean flag for caching
     */
    void setFlagCacheInEWSMcache(bool FlagCacheInEWSMcache)
    {
        this->FlagCacheInEWSMcache = FlagCacheInEWSMcache;
    }

    
    ////////////////////////////////////////////////////////////////////////     
    
    /**
     * @brief A get method to access the member reference to the object of type
     * StandardModel passed to the constructor.
     * @return the reference to the object of type StandardModel passed to the
     * constructor.
     */
    //const StandardModel& getSM() const
    //{
    //    return SM;
    //}
    
    /**
     * @brief A get method to access the member object of type PVfunctions.
     * @return the object of type PVfunctions in the current class
     */
    const PVfunctions getPV() const 
    {
        return PV;
    }

    /**
     * @brief A get method to access the member object of type Polylogarithms. 
     * @return the object of type Polylogarithms in the current class
     */
    const Polylogarithms getPolyLog() const 
    {
        return PolyLog;
    }

    /**
     * @brief A get method to access the member object of type ClausenFunctions. 
     * @return the object of type ClausenFunctions in the current class
     */
    const ClausenFunctions getClausen() const 
    {
        return Clausen;
    }

    
    ////////////////////////////////////////////////////////////////////////         
    // get functions for constants
    
    /**
     * @brief A get method to access the value of the zeta function @f$\zeta(2)@f$.
     * @return @f$\zeta(2)@f$
     */
    double getZeta2() const
    {
        return zeta2;
    }

    /**
     * @brief A get method to access the value of the zeta function @f$\zeta(3)@f$.
     * @return @f$\zeta(3)@f$
     */
    double getZeta3() const
    {
        return zeta3;
    }
    
    /**
     * @brief A get method to access the value of the zeta function @f$\zeta(4)@f$.
     * @return @f$\zeta(4)@f$
     */
    double getZeta4() const
    {
        return zeta4;
    }        
    
    /**
     * @brief A get method to access the value of the zeta function @f$\zeta(5)@f$.
     * @return @f$\zeta(5)@f$
     */
    double getZeta5() const
    {
        return zeta5;
    }

    /**
     * @brief A get method to access the constant @f$S_2@f$.
     * @details The constant @f$S_2@f$ is defined as
     * @f[
     * S_2
     * = \frac{4}{9 \sqrt{3}}\, {\rm Cl}_2 \left( \frac{\pi}{3} \right)
     * \approx 0.260434137632162\,,
     * @f]
     * which appears in three-loop amplitudes. 
     * @return @f$S_2@f$
     */
    double getS2() const
    {
        return S2;
    }
 
    /**
     * @brief A get method to access the constant @f$D_3@f$.
     * @details The constant @f$D_3@f$ is defined as
     * @f[
     * D_3
     * = 6\, \zeta(3) - \frac{15}{4}\, \zeta(4)
     * - 6 \left[ {\rm Cl}_2 \left( \frac{\pi}{3} \right) \right]^2
     * \approx -3.02700949398765\,,
     * @f]
     * which appears in three-loop amplitudes.
     * @return @f$D_3@f$
     */
    double getD3() const
    {
        return D3;
    }   
    
    /**
     * @brief A get method to access the constant @f$B_4@f$.
     * @details The constant @f$B_4@f$ is defined as
     * @f[
     * B_4
     * = 16\, {\rm Li}_4 \left( \frac{1}{2} \right)
     * - 4\, \zeta(2) \ln^2 2 + \frac{2}{3} \ln^4 2 - \frac{13}{2}\, \zeta(4)
     * \approx -1.76280008707377\,,
     * @f]
     * which appears in three-loop amplitudes. 
     * @return @f$B_4@f$
     */
    double getB4() const
    {
        return B4;
    }    

    /**
     * @brief A get method to access the constant @f$\ln 2@f$.
     * @return @f$\ln 2@f$
     */
    double getLog2() const
    {
        return log2;
    }    

    
    //////////////////////////////////////////////////////////////////////// 

    /**
     * @brief The mass of a lepton. 
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$m_l@f$
     */
    double ml(const StandardModel::lepton l) const;

    /**
     * @brief The mass squared of a lepton. 
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$m_l^2@f$
     */
    double ml2(const StandardModel::lepton l) const
    {
        return ( ml(l)*ml(l) );
    }

    /**
     * @brief The mass of a quark. 
     * @param[in] q name of a quark (see QCD::quark)
     * @param[in] mu renormalization scale
     * @param[in] order order in QCD (= LO, FULLNLO, FULLNNLO[defalut])
     * @attention If the flag #FlagDebug is set to true,
     * @f$m_{u,d,s}(2\,\mathrm{GeV})@f$, @f$m_c(m_c)@f$, @f$m_b(m_b)@f$ or
     * @f$m_t^{\mathrm{pole}}@f$ is returned. 
     * @return the MSbar mass of u, d, s, c, b at the scale mu
     * or the pole mass of t
     */
    double mq(const QCD::quark q, const double mu, 
              const orders order=FULLNNLO) const;

    /**
     * @brief The mass squared of a quark. 
     * @param[in] q name of a quark (see QCD::quark)
     * @param[in] mu renormalization scale
     * @param[in] order order in QCD (= LO, FULLNLO, FULLNNLO[defalut])
     * @attention If the flag #FlagDebug is set to true,
     * @f$(m_{u,d,s}(2\,\mathrm{GeV}))^2@f$, @f$(m_c(m_c))^2@f$, @f$(m_b(m_b))^2@f$ or
     * @f$(m_t^{\mathrm{pole}})^2@f$ is returned.
     * @return the MSbar mass squared of u, d, s, c, b at the scale mu
     * or the pole mass squared of t
     */
    double mq2(const QCD::quark q, const double mu,
               const orders order=FULLNNLO) const
    {
        return ( mq(q, mu, order)*mq(q, mu, order) );
    }
    
    /**
     * @brief The pole mass of the top quark.
     * @return @f$m_t^{\mathrm{pole}}@f$
     */
    double Mt() const 
    {
        return SM.getMtpole();
    }
    
    /**
     * @brief The pole mass squared of the top quark.
     * @return @f$(m_t^{\mathrm{pole}})^2@f$
     */
    double Mt2() const 
    {
        return ( SM.getMtpole()*SM.getMtpole() ); 
    }

    /**
     * @brief The strong coupling @f$\alpha_s(M_Z^2)@f$. 
     * @return @f$\alpha_s(M_Z^2)@f$
     */
    double alsMz() const 
    {
        return SM.getAlsMz();
    }
    
    /**
     * @brief The Fermi constant @f$G_\mu@f$.
     * @return @f$G_\mu@f$
     */
    double GF() const
    {
        return SM.getGF();
    }
    
    /**
     * @brief The fine structure constant @f$\alpha@f$.
     * @return @f$\alpha@f$
     */
    double ale() const
    {
        return SM.getAle();
    }

    /**
     * @brief The five-flavour hadronic correction to @f$\alpha@f$ at @f$M_Z^2@f$,
     * @f$\Delta\alpha_{\mathrm{had}}^{(5)}(M_Z^2)@f$
     * @return @f$\Delta\alpha_{\mathrm{had}}^{(5)}(M_Z^2)@f$
     */
    double dAle5Mz() const 
    {
        return SM.getDAle5Mz();
    }

    /**
     * @brief The @f$Z@f$-boson mass @f$M_Z@f$.
     * @return @f$M_Z@f$
     */
    double Mz() const 
    {
        return SM.getMz();
    }

    /**
     * @brief The @f$Z@f$-boson mass squared @f$M_Z^2@f$.
     * @return @f$M_Z^2@f$
     */
    double Mz2() const
    {
        return ( SM.getMz()*SM.getMz() );
    }

    /**
     * @brief The Higgs mass @f$m_h@f$.
     * @return @f$m_h@f$
     */
    double mh() const 
    {
        return SM.getMHl();
    }
    
    /**
     * @brief The Higgs mass squared @f$m_h^2@f$.
     * @return @f$m_h^2@f$
     */
    double mh2() const
    {
        return ( SM.getMHl()*SM.getMHl() );
    }

    /**
     * @brief The @f$W@f$-boson mass @f$M_W@f$.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$M_W@f$
     */
    double Mw(const double Mw_i) const
    {
        return Mw_i;
    }
    
    /**
     * @brief The @f$W@f$-boson mass squared @f$M_W^2@f$.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$M_W^2@f$
     */
    double Mw2(const double Mw_i) const
    {
        return ( Mw_i*Mw_i );
    }

    /**
     * @brief @f$c_W^2@f$.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$c_W^2@f$
     */
    double cW2(const double Mw_i) const 
    {
        return ( Mw(Mw_i)*Mw(Mw_i)/Mz()/Mz() );
    }   
    
    /**
     * @brief @f$s_W^2@f$.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$s_W^2@f$
     */
    double sW2(const double Mw_i) const 
    {
        return ( 1.0 - cW2(Mw_i) );
    } 
    
    /**
     * @brief The charge of a lepton, @f$Q_l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$Q_l@f$
     */
    double Ql(const StandardModel::lepton l) const
    {
        return SM.getLeptons(l).getCharge();
    }    

    /**
     * @brief The charge of a quark, @f$Q_q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$Q_q@f$
     */
    double Qq(const QCD::quark q) const
    {
        return SM.getQuarks(q).getCharge();
    }

    /**
     * @brief The isospin of a lepton, @f$I_3^l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$I_3^l@f$
     */
    double I3l(const StandardModel::lepton l) const
    {
        return SM.getLeptons(l).getIsospin();
    }

    /**
     * @brief The isospin of a quark, @f$I_3^q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$I_3^q@f$
     */
    double I3q(const QCD::quark q) const
    {
        return SM.getQuarks(q).getIsospin();
    }

    /**
     * @brief The tree-level vector coupling for @f$Z\to l\bar{l}@f$, 
     * denoted as @f$v_l@f$. 
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$v_l@f$
     */
    double vl(const StandardModel::lepton l, const double Mw_i) const 
    {
        return ( al(l) - 2.0*Ql(l)*sW2(Mw_i) );
    }

    /**
     * @brief The tree-level vector coupling for @f$Z\to q\bar{q}@f$,
     * denoted as @f$v_q@f$.
     * @param[in] q nname of a quark (see QCD::quark)
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$v_q@f$
     */
    double vq(const QCD::quark q, const double Mw_i) const
    {
        return ( aq(q) - 2.0*Qq(q)*sW2(Mw_i) );
    }
    
    /**
     * @brief The tree-level axial-vector coupling for @f$Z\to l\bar{l}@f$,
     * denoted as @f$a_l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$a_l@f$
     */
    double al(const StandardModel::lepton l) const
    {
        return ( SM.getLeptons(l).getIsospin() );
    }

    /**
     * @brief The tree-level axial-vector coupling for @f$Z\to q\bar{q}@f$,
     * denoted as @f$a_q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$a_q@f$
     */
    double aq(const QCD::quark q) const
    {
        return ( SM.getQuarks(q).getIsospin() );
    }
    
    /**
     * @brief @f$\sigma_l = |v_l+a_l|@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\sigma_l@f$
     */
    double sigmal(const StandardModel::lepton l, const double Mw_i) const
    {
        return ( 1.0 - 2.0*fabs(Ql(l))*sW2(Mw_i) );
    }

    /**
     * @brief @f$\sigma_q = |v_q+a_q|@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\sigma_q@f$
     */
    double sigmaq(const QCD::quark q, const double Mw_i) const 
    {
        return ( 1.0 - 2.0*fabs(Qq(q))*sW2(Mw_i) );
    }  
    
    /**
     * @brief @f$\delta_l = v_l-a_l@f$.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta_l@f$
     */    
    double deltal(const StandardModel::lepton l, const double Mw_i) const 
    {
        return ( - 2.0*Ql(l)*sW2(Mw_i) );   
    } 

    /**
     * @brief @f$\delta_q = v_q-a_q@f$.
     * @param[in] q name of a quark (see QCD::quark)
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\delta_q@f$
     */    
    double deltaq(const QCD::quark q, const double Mw_i) const 
    {
        return ( - 2.0*Qq(q)*sW2(Mw_i) );   
    }  
    
    /**
     * @brief The conversion factor from @f$\alpha@f$ to @f$G_\mu@f$, denoted
     * as @f$f@f$.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$ f = \sqrt{2}G_\mu M_Z^2 s_W^2 c_W^2/(\pi\alpha)@f$.
     */
    double f_AlphaToGF(const double Mw_i) const 
    {
        return ( sqrt(2.0)*GF()*pow(Mz(),2.0)*sW2(Mw_i)*cW2(Mw_i)/M_PI/ale() );
    }
   
    /**
     * @brief The quantity @f$X_t@f$ with the coupling @f$G_\mu@f$.
     * @return @f$X_t^{G_\mu}=G_\mu m_t^2/(8\sqrt{2}\pi^2)@f$
     */
    double Xt_GF() const 
    {
        return ( GF()*Mt()*Mt()/8.0/sqrt(2.0)/M_PI/M_PI );
    }
    
    /**
     * @brief The quantity @f$X_t@f$ with the coupling @f$\alpha@f$.
     * @details See also Xt_GF() and f_AlphaToGF(). 
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$X_t^{\alpha}@f$
     */
    double Xt_alpha(const double Mw_i) const 
    {
        return ( Xt_GF()/f_AlphaToGF(Mw_i) );
    }

    /**
     * @brief The strong coupling @f$\alpha_s(\mu)@f$. 
     * @param[in] mu renormalization scale @f$\mu@f$ in GeV
     * @param[in] order order in QCD (see #orders)
     * @return @f$\alpha_s(\mu)@f$ at the given order
     */
    double Als(const double mu, const orders order) const
    {
        return ( SM.Als(mu, order) );
    }

    /**
     * @brief The strong coupling @f$\alpha_s(m_t^2)@f$ at NNLO.
     * @attention The constant value @f$\alpha_s(m_t^2)=0.1074432788@f$ is
     * returned when the flag #FlagDebug is set to true. 
     * @return @f$\alpha_s(m_t^2)@f$ at NNLO
     */
    double alsMt() const 
    {
        if (FlagDebug)
            return ( 0.1074432788 );// for debug
        else
            return ( SM.Als(Mt(),FULLNNLO) );
    }

    /**
     * @brief The constant @f$\Phi@f$ for two-loop %QCD contribution.
     * @return @f$\Phi = {\rm arcsin}(\sqrt{r})@f$ with @f$r = M_Z^2/(4m_t^2)@f$
     */
    double Phi_QCD2() const 
    {
        double r_QCD2 = Mz()*Mz()/4.0/Mt()/Mt();
        return ( asin(sqrt(r_QCD2)) );
    }

    /**
     * @brief The constant @f$\gamma@f$ for two-loop %QCD contribution. 
     * @return @f$\gamma = \ln (2\sqrt{r})@f$ with @f$r = M_Z^2/(4m_t^2)@f$
     */
    double gamma_QCD2() const 
    {
        double r_QCD2 = Mz()*Mz()/4.0/Mt()/Mt();
        return ( log(2.0*sqrt(r_QCD2)) );
    }

    /**
     * @brief The constant @f$h@f$ for two-loop %QCD contribution.
     * @return @f$h= \ln(2\sqrt{1-r})@f$ with @f$r = M_Z^2/(4m_t^2)@f$
     */
    double h_QCD2() const 
    {
        double r_QCD2 = Mz()*Mz()/4.0/Mt()/Mt();
        return ( log(2.0*sqrt(1.0-r_QCD2)) );
    }
    
    /**
     * @brief A logarithm appearing in the functions @f$V_1'@f$ and @f$A_1'@f$
     * for two-loop %QCD contribution.
     * @return @f$\mathrm{Re}[\ln (1-e^{2i\Phi}) - 2\ln (1-e^{4i\Phi})]@f$ with
     * @f$\Phi=@f$Phi_QCD2()
     */
    double logV1primeAndA1prime() const 
    {
        gsl_complex OneMinusE2Iphi = gsl_complex_rect(1.0-cos(2.0*Phi_QCD2()), 
                                     -sin(2.0*Phi_QCD2()));
        gsl_complex OneMinusE4Iphi = gsl_complex_rect(1.0-cos(4.0*Phi_QCD2()), 
                                     -sin(4.0*Phi_QCD2()));
        return ( GSL_REAL(gsl_complex_log(OneMinusE2Iphi))
                - 2.0*GSL_REAL(gsl_complex_log(OneMinusE4Iphi)) );
    }
    
    /**
     * @brief The constant @f${\rm Cl}_2(2 \Phi)@f$.
     * @details This constant appears in two-loop %QCD contribution, where
     * @f$\Phi=\mathrm{arcsin}\big(M_Z/(2m_t)\big)@f$.
     * @return @f${\rm Cl}_2(2 \Phi)@f$
     */
    double Cl2_2Phi() const 
    {
        double Phi= asin(Mz()/2.0/Mt());
        return ( Clausen.Cl2(2.0*Phi) );
    }

    /**
     * @brief The constant @f${\rm Cl}_2(4 \Phi)@f$.
     * @details This constant appears in two-loop %QCD contribution, where
     * @f$\Phi=\mathrm{arcsin}\big(M_Z/(2m_t)\big)@f$.
     * @return @f${\rm Cl}_2(4 \Phi)@f$
     */
    double Cl2_4Phi() const 
    {
        double Phi= asin(Mz()/2.0/Mt());
        return ( Clausen.Cl2(4.0*Phi) );
    }

    /**
     * @brief The constant @f${\rm Cl}_3(2 \Phi)@f$.
     * @details This constant appears in two-loop %QCD contribution, where
     * @f$\Phi=\mathrm{arcsin}\big(M_Z/(2m_t)\big)@f$.
     * @return @f${\rm Cl}_3(2 \Phi)@f$
     */
    double Cl3_2Phi() const 
    {
        double Phi= asin(Mz()/2.0/Mt());
        return ( Clausen.Cl3(2.0*Phi) );
    }

    /**
     * @brief The constant @f${\rm Cl}_3(4 \Phi)@f$.
     * @details This constant appears in two-loop %QCD contribution, where
     * @f$\Phi=\mathrm{arcsin}\big(\sqrt{M_Z/(2m_t)}\big)@f$.
     * @return @f${\rm Cl}_3(4 \Phi)@f$
     */
    double Cl3_4Phi() const 
    {
        double Phi= asin(Mz()/2.0/Mt());
        return ( Clausen.Cl3(4.0*Phi) );
    }

    
    ////////////////////////////////////////////////////////////////////////  
    // Functions for the caches

    /**
     * @brief A cache method.
     * @return @f$\ln (M_Z/m_e)@f$
     */
    double logMZtoME() const;

    /**
     * @brief A cache method.
     * @return @f$\ln (M_Z/m_\mu)@f$
     */
    double logMZtoMMU() const;

    /**
     * @brief A cache method.
     * @return @f$\ln (M_Z/m_\tau)@f$
     */
    double logMZtoMTAU() const;

    /**
     * @brief A cache method.
     * @return @f$\ln (M_Z/m_t)@f$
     */
    double logMZtoMTOP() const;

    /**
     * @brief A cache method.
     * @return @f$\ln (m_t/m_h)@f$
     */
    double logMTOPtoMH() const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\ln c_W^2@f$
     */
    double log_cW2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\mathrm{Li}_2(M_W^2/m_t^2)@f$
     */
    double Li2_MW2toMTOP2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\mathrm{Li}_3(M_W^2/m_t^2)@f$
     */
    double Li3_MW2toMTOP2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$\mathrm{Li}_3(-M_W^2/m_t^2/(1-M_W^2/m_t^2))@f$
     */
    double Li3_for_F1(const double Mw_i) const;    

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$A_0(M_W^2)@f$ with @f$\mu=M_Z@f$
     */
    double A0_Mz2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @return @f$A_0(m_h^2)@f$ with @f$\mu=M_Z@f$
     */
    double A0_Mz2_mh2() const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$A_0(M_Z^2)@f$ with @f$\mu=M_W@f$
     */
    double A0_Mw2_Mz2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$A_0(m_h^2)@f$ with @f$\mu=M_W@f$
     */
    double A0_Mw2_mh2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @return @f$A_0(M_Z^2)@f$ with @f$\mu=M_Z@f$
     */
    double A0_Mz2_Mz2() const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$A_0(M_W^2)@f$ with @f$\mu=M_W@f$
     */
    double A0_Mw2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_0(M_W^2;m_h^2,M_W^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0_Mz2_Mw2_mh2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_0(0;m_h^2,M_W^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0_Mz2_0_mh2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_0(M_Z^2;m_t^2,m_t^2)@f$ with @f$\mu=M_W@f$
     */
    complex B0_Mw2_Mz2_Mt2_Mt2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_0(M_Z^2;M_W^2,M_W^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0_Mz2_Mz2_Mw2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @return @f$B_0(M_Z^2;m_h^2,M_Z^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0_Mz2_Mz2_mh2_Mz2() const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_0(M_W^2;M_Z^2,M_W^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0_Mz2_Mw2_Mz2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_0(M_W^2;0,M_W^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0_Mz2_Mw2_0_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_0(0;M_Z^2,M_W^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0_Mz2_0_Mz2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_0(0;0,M_W^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0_Mz2_0_0_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_0(M_Z^2;M_W^2,M_W^2)@f$ with @f$\mu=M_W@f$
     */
    complex B0_Mw2_Mz2_Mw2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_0(M_W^2;M_Z^2,M_W^2)@f$ with @f$\mu=M_W@f$
     */
    complex B0_Mw2_Mw2_Mz2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_0(M_W^2;m_h^2,M_W^2)@f$ with @f$\mu=M_W@f$
     */
    complex B0_Mw2_Mw2_mh2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_0(M_W^2;0,M_W^2)@f$ with @f$\mu=M_W@f$
     */
    complex B0_Mw2_Mw2_0_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$B_0(M_Z^2;m_l^2,m_l^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0_Mz2_Mz2_ml2_ml2(const StandardModel::lepton l) const;

    /**
     * @brief A cache method.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$B_0(M_Z^2;m_q^2,m_q^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0_Mz2_Mz2_mq2_mq2(const QCD::quark q) const;    

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_{0p}(0;m_h^2,M_W^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0p_Mz2_0_mh2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @return @f$B_{0p}(M_Z^2;m_h^2,M_Z^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0p_Mz2_Mz2_mh2_Mz2() const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_{0p}(0;M_Z^2,M_W^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0p_Mz2_0_Mz2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_{0p}(M_Z^2;M_W^2,M_W^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0p_Mz2_Mz2_Mw2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_{0p}(M_W^2;M_Z^2,M_W^2)@f$ with @f$\mu=M_W@f$
     */
    complex B0p_Mw2_Mw2_Mz2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_{0p}(M_W^2;m_h^2,M_W^2)@f$ with @f$\mu=M_W@f$
     */
    complex B0p_Mw2_Mw2_mh2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_{0p}(M_W^2;0,M_W^2)@f$ with @f$\mu=M_W@f$
     */
    complex B0p_Mw2_Mw2_0_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$B_{0p}(M_Z^2;m_l^2,m_l^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0p_Mz2_Mz2_ml2_ml2(const StandardModel::lepton l) const;

    /**
     * @brief A cache method.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$B_{0p}(M_Z^2;m_q^2,m_q^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B0p_Mz2_Mz2_mq2_mq2(const QCD::quark q) const;
    
    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a lepton doublet
     * @return @f$B_1(0;m_l^2,m_{l'}^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B1_Mz2_0_ml2_mlprime2(const int gen) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a quark doublet
     * @return @f$B_1(0;m_q^2,m_{q'}^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B1_Mz2_0_mq2_mqprime2(const int gen) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a lepton doublet
     * @return @f$B_1(0;m_{l'}^2,m_l^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B1_Mz2_0_mlprime2_ml2(const int gen) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a quark doublet
     * @return @f$B_1(0;m_{q'}^2,m_q^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B1_Mz2_0_mqprime2_mq2(const int gen) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a lepton doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_1(M_W^2;m_l^2,m_{l'}^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B1_Mz2_Mw2_ml2_mlprime2(const int gen, const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a quark doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_1(M_W^2;m_q^2,m_{q'}^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B1_Mz2_Mw2_mq2_mqprime2(const int gen, const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a lepton doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_1(M_W^2;m_{l'}^2,m_l^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B1_Mz2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a quark doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_1(M_W^2;m_{q'}^2,m_q^2)@f$ with @f$\mu=M_Z@f$
     */
    complex B1_Mz2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a lepton doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_{1p}(M_W^2;m_l^2,m_{l'}^2)@f$ with @f$\mu=M_W@f$
     */
    complex B1p_Mw2_Mw2_ml2_mlprime2(const int gen, const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a quark doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_{1p}(M_W^2;m_q^2,m_{q'}^2)@f$ with @f$\mu=M_W@f$
     */
    complex B1p_Mw2_Mw2_mq2_mqprime2(const int gen, const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a lepton doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_{1p}(M_W^2;m_{l'}^2,m_l^2)@f$ with @f$\mu=M_W@f$
     */
    complex B1p_Mw2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a quark doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_{1p}(M_W^2;m_{q'}^2,m_q^2)@f$ with @f$\mu=M_W@f$
     */
    complex B1p_Mw2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$B_f(M_Z^2;m_l^2,m_l^2)@f$ with @f$\mu=M_Z@f$
     */
    complex Bf_Mz2_Mz2_ml2_ml2(const StandardModel::lepton l) const;

    /**
     * @brief A cache method.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$B_f(M_Z^2;m_q^2,m_q^2)@f$ with @f$\mu=M_Z@f$
     */
    complex Bf_Mz2_Mz2_mq2_mq2(const QCD::quark q) const;

    /**
     * @brief A cache method.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$B_f(0;m_l^2,m_l^2)@f$ with @f$\mu=M_Z@f$
     */
    complex Bf_Mz2_0_ml2_ml2(const StandardModel::lepton l) const;

    /**
     * @brief A cache method.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$B_f(0;m_q^2,m_q^2)@f$ with @f$\mu=M_Z@f$
     */
    complex Bf_Mz2_0_mq2_mq2(const QCD::quark q) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a lepton doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_f(M_W^2;m_{l'}^2,m_l^2)@f$ with @f$\mu=M_Z@f$
     */
    complex Bf_Mz2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a quark doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_f(M_W^2;m_{q'}^2,m_q^2)@f$ with @f$\mu=M_Z@f$
     */
    complex Bf_Mz2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a lepton doublet
     * @return @f$B_f(0;m_{l'}^2,m_l^2)@f$ with @f$\mu=M_Z@f$
     */
    complex Bf_Mz2_0_mlprime2_ml2(const int gen) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a quark doublet
     * @return @f$B_f(0;m_{q'}^2,m_q^2)@f$ with @f$\mu=M_Z@f$
     */
    complex Bf_Mz2_0_mqprime2_mq2(const int gen) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a lepton doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_f(M_W^2;m_{l'}^2,m_l^2)@f$ with @f$\mu=M_W@f$
     */
    complex Bf_Mw2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a quark doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_f(M_W^2;m_{q'}^2,m_q^2)@f$ with @f$\mu=M_W@f$
     */
    complex Bf_Mw2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const;
    
    /**
     * @brief A cache method.
     * @param[in] l name of a lepton (see StandardModel::lepton)
     * @return @f$B_{fp}(M_Z^2;m_l^2,m_l^2)@f$ with @f$\mu=M_Z@f$
     */
    complex Bfp_Mz2_Mz2_ml2_ml2(const StandardModel::lepton l) const;

    /**
     * @brief A cache method.
     * @param[in] q name of a quark (see QCD::quark)
     * @return @f$B_{fp}(M_Z^2;m_q^2,m_q^2)@f$ with @f$\mu=M_Z@f$
     */
    complex Bfp_Mz2_Mz2_mq2_mq2(const QCD::quark q) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a lepton doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_{fp}(M_W^2;m_{l'}^2,m_l^2)@f$ with @f$\mu=M_W@f$
     */
    complex Bfp_Mw2_Mw2_mlprime2_ml2(const int gen, const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] gen the generation index of a quark doublet
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$B_{fp}(M_W^2;m_{q'}^2,m_q^2)@f$ with @f$\mu=M_W@f$
     */
    complex Bfp_Mw2_Mw2_mqprime2_mq2(const int gen, const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$C_0(M_Z^2;M_W^2,m_t^2,M_W^2)@f$
     */
    complex C0_Mz2_Mw2_Mt2_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$C_0(M_Z^2;M_t^2,M_W^2,m_t^2)@f$
     */
    complex C0_Mz2_Mt2_Mw2_Mt2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$C_0(M_Z^2;0,M_W^2,0)@f$
     */
    complex C0_Mz2_0_Mw2_0(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$C_0(M_Z^2;M_W^2,0,M_W^2)@f$
     */
    complex C0_Mz2_Mw2_0_Mw2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$C_0(M_W^2;M_W^2,0,M_Z^2)@f$
     */
    complex C0_Mw2_Mw2_0_Mz2(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @param[in] Mw_i the @f$W@f$-boson mass
     * @return @f$C_0(M_W^2;0,M_Z^2,0)@f$
     */
    complex C0_Mw2_0_Mz2_0(const double Mw_i) const;

    /**
     * @brief A cache method.
     * @return @f$C_0(M_Z^2;0,M_Z^2,0)@f$
     */
    complex C0_Mz2_0_Mz2_0() const;
    
    
    //////////////////////////////////////////////////////////////////////// 

private:
    bool FlagDebug;///< A flag for debugging (false by default).
    bool FlagCacheInEWSMcache;///< A flag for caching (true by default).
    
    const StandardModel& SM;///< A reference to an object of type StandardModel.
    const PVfunctions PV;///< An object of type PVfunctions.
    const ClausenFunctions Clausen;///< An object of type ClausenFunctions.
    const Polylogarithms PolyLog;///< An object of type Polylogarithms.
    
    // Constants 
    double zeta2;///< The constant @f$\zeta(2)@f$.
    double zeta3;///< The constant @f$\zeta(3)@f$.
    double zeta4;///< The constant @f$\zeta(4)@f$.
    double zeta5;///< The constant @f$\zeta(5)@f$.
    double S2;///< The constant @f$S_2=(4/9/\sqrt{3})\mathrm{Cl}_2(\pi/3)@f$.
    double D3;///< The constant @f$D_3=6\zeta(3) - (15/4)\zeta(4) - 6[\mathrm{Cl}_2(\pi/3)]^2@f$.
    double B4;///< The constant @f$B_4=16\mathrm{Li}_4(1/2) - 4\zeta(2)\ln^2(2) + (2/3)\ln^4(2) - (13/2)\zeta(4)@f$.
    double log2;///< The constant @f$\ln 2@f$. 

    
    ////////////////////////////////////////////////////////////////////////     
    // Caches for double variables
    // The size of the caches are equal to the number of relevant parameters
    // plus one. 
    
    /* Logarithms */
    mutable double logMZtoME_cache[3];///< A cache of @f$\ln (M_Z/m_e)@f$.
    mutable double logMZtoMMU_cache[3];///< A cache of @f$\ln (M_Z/m_\mu)@f$.
    mutable double logMZtoMTAU_cache[3];///< A cache of @f$\ln (M_Z/m_\tau)@f$.
    mutable double logMZtoMTOP_cache[3];///< A cache of @f$\ln (M_Z/m_t)@f$.
    mutable double logMTOPtoMH_cache[3];///< A cache of @f$\ln (m_t/m_h)@f$.
    mutable double log_cW2_cache[3];///< A cache of @f$\ln c_W^2@f$.
    
    /* Dilogarithm and Trilogarithm for two-loop QCD corrections */
    mutable double Li2_MW2toMTOP2_cache[3];///< A cache of @f$\mathrm{Li}_2(M_W^2/m_t^2)@f$.
    mutable double Li3_MW2toMTOP2_cache[3];///< A cache of @f$\mathrm{Li}_3(M_W^2/m_t^2)@f$.
    mutable double Li3_for_F1_cache[3];///< A cache of @f$\mathrm{Li}_3(-M_W^2/m_t^2/(1-M_W^2/m_t^2))@f$.
    
    /* One-loop functions */
    mutable double A0_Mz2_Mw2_cache[3];///< A cache of a PV function.
    mutable double A0_Mz2_mh2_cache[3];///< A cache of a PV function.
    mutable double A0_Mw2_Mz2_cache[3];///< A cache of a PV function.
    mutable double A0_Mw2_mh2_cache[3];///< A cache of a PV function.
    mutable double A0_Mz2_Mz2_cache[2];///< A cache of a PV function.
    mutable double A0_Mw2_Mw2_cache[2];///< A cache of a PV function.

    
    ////////////////////////////////////////////////////////////////////////     
    // Caches for complex variables
    // The size of the caches are equal to the number of relevant parameters
    // plus one. 
    
    mutable double B0_Mz2_Mw2_mh2_Mw2_cache[5];///< A cache of a PV function.
    mutable double B0_Mz2_0_mh2_Mw2_cache[5];///< A cache of a PV function.
    mutable double B0_Mw2_Mz2_Mt2_Mt2_cache[5];///< A cache of a PV function.
    mutable double B0_Mz2_Mz2_Mw2_Mw2_cache[4];///< A cache of a PV function.
    mutable double B0_Mz2_Mz2_mh2_Mz2_cache[4];///< A cache of a PV function.
    mutable double B0_Mz2_Mw2_Mz2_Mw2_cache[4];///< A cache of a PV function.
    mutable double B0_Mz2_Mw2_0_Mw2_cache[4];///< A cache of a PV function.
    mutable double B0_Mz2_0_Mz2_Mw2_cache[4];///< A cache of a PV function.
    mutable double B0_Mz2_0_0_Mw2_cache[4];///< A cache of a PV function.
    mutable double B0_Mw2_Mz2_Mw2_Mw2_cache[4];///< A cache of a PV function.
    mutable double B0_Mw2_Mw2_Mz2_Mw2_cache[4];///< A cache of a PV function.
    mutable double B0_Mw2_Mw2_mh2_Mw2_cache[4];///< A cache of a PV function.
    mutable double B0_Mw2_Mw2_0_Mw2_cache[3];///< A cache of a PV function.
    mutable double B0_Mz2_Mz2_ml2_ml2_cache[6][4];///< A cache of a PV function.
    mutable double B0_Mz2_Mz2_mq2_mq2_cache[6][4];///< A cache of a PV function.
    
    mutable double B0p_Mz2_0_mh2_Mw2_cache[5];///< A cache of a PV function.
    mutable double B0p_Mz2_Mz2_mh2_Mz2_cache[4];///< A cache of a PV function.
    mutable double B0p_Mz2_0_Mz2_Mw2_cache[4];///< A cache of a PV function.
    mutable double B0p_Mz2_Mz2_Mw2_Mw2_cache[4];///< A cache of a PV function.
    mutable double B0p_Mw2_Mw2_Mz2_Mw2_cache[4];///< A cache of a PV function.
    mutable double B0p_Mw2_Mw2_mh2_Mw2_cache[4];///< A cache of a PV function.
    mutable double B0p_Mw2_Mw2_0_Mw2_cache[3];///< A cache of a PV function.
    mutable double B0p_Mz2_Mz2_ml2_ml2_cache[6][4];///< A cache of a PV function.
    mutable double B0p_Mz2_Mz2_mq2_mq2_cache[6][4];///< A cache of a PV function.
    
    mutable double B1_Mz2_0_ml2_mlprime2_cache[3][5];///< A cache of a PV function.
    mutable double B1_Mz2_0_mq2_mqprime2_cache[3][5];///< A cache of a PV function.
    mutable double B1_Mz2_0_mlprime2_ml2_cache[3][5];///< A cache of a PV function.
    mutable double B1_Mz2_0_mqprime2_mq2_cache[3][5];///< A cache of a PV function.
    mutable double B1_Mz2_Mw2_ml2_mlprime2_cache[3][6];///< A cache of a PV function.
    mutable double B1_Mz2_Mw2_mq2_mqprime2_cache[3][6];///< A cache of a PV function.
    mutable double B1_Mz2_Mw2_mlprime2_ml2_cache[3][6];///< A cache of a PV function.
    mutable double B1_Mz2_Mw2_mqprime2_mq2_cache[3][6];///< A cache of a PV function.
    
    mutable double B1p_Mw2_Mw2_ml2_mlprime2_cache[3][5];///< A cache of a PV function.
    mutable double B1p_Mw2_Mw2_mq2_mqprime2_cache[3][5];///< A cache of a PV function.
    mutable double B1p_Mw2_Mw2_mlprime2_ml2_cache[3][5];///< A cache of a PV function.
    mutable double B1p_Mw2_Mw2_mqprime2_mq2_cache[3][5];///< A cache of a PV function.
        
    mutable double Bf_Mz2_Mz2_ml2_ml2_cache[6][4];///< A cache of a PV function.
    mutable double Bf_Mz2_Mz2_mq2_mq2_cache[6][4];///< A cache of a PV function.
    mutable double Bf_Mz2_0_ml2_ml2_cache[6][4];///< A cache of a PV function.
    mutable double Bf_Mz2_0_mq2_mq2_cache[6][4];///< A cache of a PV function.
    mutable double Bf_Mz2_Mw2_mlprime2_ml2_cache[3][6];///< A cache of a PV function.
    mutable double Bf_Mz2_Mw2_mqprime2_mq2_cache[3][6];///< A cache of a PV function.
    mutable double Bf_Mz2_0_mlprime2_ml2_cache[3][5];///< A cache of a PV function.
    mutable double Bf_Mz2_0_mqprime2_mq2_cache[3][5];///< A cache of a PV function.
    mutable double Bf_Mw2_Mw2_mlprime2_ml2_cache[3][5];///< A cache of a PV function.
    mutable double Bf_Mw2_Mw2_mqprime2_mq2_cache[3][5];///< A cache of a PV function.

    mutable double Bfp_Mz2_Mz2_ml2_ml2_cache[6][4];///< A cache of a PV function.
    mutable double Bfp_Mz2_Mz2_mq2_mq2_cache[6][4];///< A cache of a PV function.
    mutable double Bfp_Mw2_Mw2_mlprime2_ml2_cache[3][5];///< A cache of a PV function.
    mutable double Bfp_Mw2_Mw2_mqprime2_mq2_cache[3][5];///< A cache of a PV function.
    
    mutable double C0_Mz2_Mw2_Mt2_Mw2_cache[5];///< A cache of a PV function.
    mutable double C0_Mz2_Mt2_Mw2_Mt2_cache[5];///< A cache of a PV function.
    mutable double C0_Mz2_0_Mw2_0_cache[4];///< A cache of a PV function.
    mutable double C0_Mz2_Mw2_0_Mw2_cache[4];///< A cache of a PV function.
    mutable double C0_Mw2_Mw2_0_Mz2_cache[4];///< A cache of a PV function.
    mutable double C0_Mw2_0_Mz2_0_cache[4];///< A cache of a PV function.
    mutable double C0_Mz2_0_Mz2_0_cache[3];///< A cache of a PV function.

    
    ////////////////////////////////////////////////////////////////////////     
    // Methods for caching, which will be compiled as inline functions.

    /**
     * @brief A check method for caching. 
     * @details This function checks if the values of the parameters in the array
     * params[] are all identical to those stored in the array cache[]. When they
     * are identical to each other, the current function returns true. Otherwise,
     * the function returns false.
     * @param[in] cache a cache of the parameters to be checked
     * @param[in] NumPar the number of the parameters to be checked
     * @param[in] params the parameters to be checked
     * @attention If the flag #FlagCacheInEWSMcache is set to false, the current
     * function always returns false.
     * @return true (false) if the parameters in params[] are (not) identical to
     * those in cache[].
     */
    bool CacheCheck(const double cache[], 
                    const int NumPar, const double params[]) const 
    {
        if (!FlagCacheInEWSMcache) return false;
        bool bCache = true;
        for(int i=0; i<NumPar; ++i)
            bCache &= (params[i] == cache[i]);
        return bCache;
    }

    /**
     * @brief A method to update a cache of the parameters and the quantity
     * under consideration.
     * @details This function updates cache[] with params[] and newResult,
     * where newResult for the quantity under consideration depends on the
     * parameters in params[]. Both the values of params[] and that of newResult
     * are stored into cache[], whose last element corresponds to the latter.
     * @param[out] cache a cache of the parameters and the quantity
     * @param[in] NumPar the number of the parameters
     * @param[in] params an array of the parameters
     * @param[in] newResult the new result of the quantity
     * @attention If the flag #FlagCacheInEWSMcache is set to false, the current
     * function does not modify cache[].
     */
    void newCacheForDouble(double cache[], const int NumPar, 
                           const double params[], const double newResult) const 
    {
        if (!FlagCacheInEWSMcache) return;
        for(int i=0; i<NumPar; ++i)
            cache[i] = params[i];
        cache[NumPar] = newResult;
    }

    /**
     * @brief A method to update a cache of the parameters and the quantity
     * under consideration.
     * @details This function updates cache[] with params[] and newResult,
     * where newResult for the quantity under consideration depends on the
     * parameters in params[]. Both the values of params[] and that of newResult
     * are stored into cache[], whose last two elements correspond to the real
     * and imaginary parts of the latter.
     * @param[out] cache a cache of the parameters and the quantity
     * @param[in] NumPar the number of the parameters
     * @param[in] params an array of the parameters
     * @param[in] newResult the new result of the quantity
     * @attention If the flag #FlagCacheInEWSMcache is set to false, the current
     * function does not modify cache[].
     */
    void newCacheForComplex(double cache[], const int NumPar, 
                            const double params[], const complex newResult) const 
    {
        if (!FlagCacheInEWSMcache) return;
        for(int i=0; i<NumPar; ++i)
            cache[i] = params[i];
        cache[NumPar] = newResult.real();
        cache[NumPar+1] = newResult.imag();
    }
    
};

#endif	/* EWSMCACHE_H */


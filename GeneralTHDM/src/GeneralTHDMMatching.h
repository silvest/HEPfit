/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMMATCHING_H
#define	GENERALTHDMMATCHING_H

 #include <Polylogarithms.h>
#include "gslpp.h"
#include "StandardModelMatching.h"

class GeneralTHDM;

/**
 * @class GeneralTHDMMatching
 * @ingroup GeneralTHDM
 * @brief 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class GeneralTHDMMatching : public StandardModelMatching {
public:
    GeneralTHDMMatching(const GeneralTHDM & GeneralTHDM_i);

    
    
    
    
    
    /**
     * @brief Wilson coefficient for \f$ (g-2)_{\mu} \f$.
     * @return
     */
    virtual std::vector<WilsonCoefficient>& CMgminus2mu();

    /** Calculates the muon g-2 at LO**/
    /**
     * @brief Calculates amplitudes for \f$ (g-2)_{\mu} \f$ at one loop from 1502.04199, before \cite Broggio:2014mna was used.
     * @return 
     */    
    virtual double gminus2muLO();
    //FOR THE MOMENT WE LEAVE THIS DEFINITION OF THE LO CONTRIBUTION TO THE g-2, MUST BE CHECKED
    
    
    
    /** One-loop function for the contribution to the g-2 **/
    /**
     * @brief Loop \f$ (g-2)_{\mu} \f$ at one loop from 1502.04199. The complete loop function is \int_0^1 dx \frac{x^2(2-x)}{ratio_sq*x^2-x+1}
     * but for small values of ratio_sq we can expand and take only the first term for which the functions becomes (-7.0/6.0-log(ratio_sq))
     * @return 
     */    
    double F1oneloopgm2(const double ratio_sq);
    
    
    
    /** One-loop function for the contribution to the g-2 **/
    /**
     * @brief Loop \f$ (g-2)_{\mu} \f$ at one loop from 1502.04199. The complete loop function is \int_0^1 dx \frac{-x^3}{ratio_sq*x^2-x+1}
     * but for small values of ratio_sq we can expand and take only the first term for which the functions becomes (11.0/6.0+log(ratio_sq))
     * @return 
     */    
    double F2oneloopgm2(const double ratio_sq);
    
    
    
    /** One-loop function for the contribution to the g-2 **/
    /**
     * @brief Loop \f$ (g-2)_{\mu} \f$ at one loop from 1502.04199. The complete loop function is \int_0^1 dx \frac{x^2(1-x)}{ratio_sq*x*(1-x)-x}
     * but for small values of ratio_sq we can expand and take only the first term for which the functions becomes (-1/12. - ratio_sq/60.)
     * @return 
     */    
    double F3oneloopgm2(const double ratio_sq);
    
    
    
    /** Two-loop (Barr-Zee) function for the contribution to the g-2 **/
    /**
     * @brief Loop \f$ (g-2)_{\mu} \f$ at two loops (Barr-Zee) from 1502.04199. The complete loop function is \frac{ratio_sq}{2}*\int_0^1 dx \frac{2x(1-x)-1}{ratio_sq-x(1-x)}\log(\frac{ratio_sq}{x(1-x)})
     * 
     * @return 
     */    
    double F1twoloopgm2(const double ratio_sq);
    
    
    
    /** Two-loop (Barr-Zee) function for the contribution to the g-2 **/
    /**
     * @brief Loop \f$ (g-2)_{\mu} \f$ at two loops (Barr-Zee) from 1502.04199. The complete loop function is \frac{ratio_sq}{2}*\int_0^1 dx \frac{1}{ratio_sq-x(1-x)}\log(\frac{ratio_sq}{x(1-x)})
     * 
     * @return 
     */    
    double F1tildetwoloopgm2(const double ratio_sq);
    
    
    /** Two-loop (Barr-Zee) function for the contribution to the g-2 **/
    /**
     * @brief Loop \f$ (g-2)_{\mu} \f$ at two loops (Barr-Zee) from 1502.04199. The complete loop function is \frac{1}{2}*\int_0^1 dx \frac{x(x-1)}{ratio_sq-x(1-x)}\log(\frac{ratio_sq}{x(1-x)})
     * 
     * @return 
     */    
    double F2twoloopgm2(const double ratio_sq);
    
    
    
    /** Two-loop (Barr-Zee) function for the contribution to the g-2 **/
    /**
     * @brief Loop \f$ (g-2)_{\mu} \f$ at two loops (Barr-Zee) from 1502.04199. The complete loop function is \frac{1}{2}*\int_0^1 dx \frac{x(3x(4x-1)+10)ratio_sq-x(1-x)}{ratio_sq-x(1-x)}\log(\frac{ratio_sq}{x(1-x)})
     * 
     * @return 
     */    
    double F3twoloopgm2(const double ratio_sq);
    
    
    
    /** Two-loop (Barr-Zee) function for the contribution to the g-2 **/
    /**
     * @brief Loop \f$ (g-2)_{\mu} \f$ at two loops (Barr-Zee) from 1502.04199. In the notation of [1502.04199] this corresponds to (Q_t*x+Q_b*(1-x))x(1+x)G(ratio_sq,0), i.e. we neglect mb/mHp and mb/mW which is a great approximation
     * The complete loop function is \int_0^1 dx (2/3 *x-1/3 *(1-x))x(1+x)*\frac{\log{\frac{ratio_sq*x}{x(1-x)}}}{x(1-x-ratio_sq*x)})
     * 
     * @return 
     */    
    double F4twoloopgm2(const double ratio_sq);
    
    
    
    
    /** Two-loop (Barr-Zee) function for the contribution to the g-2 **/
    /**
     * @brief Loop \f$ (g-2)_{\mu} \f$ at two loops (Barr-Zee) from 1502.04199. In the notation of [1502.04199] this corresponds to (Q_t*x+Q_b*(1-x))x(1-x)G(ratio_sq,0), i.e. we neglect mb/mHp and mb/mW which is a great approximation
     * The complete loop function is \int_0^1 dx (2/3 *x-1/3 *(1-x))x(1-x)*\frac{\log{\frac{ratio_sq*x}{x(1-x)}}}{x(1-x-ratio_sq*x)})
     * 
     * @return 
     */    
    double F5twoloopgm2(const double ratio_sq);
    
    
    
    /** Calculates the NLO contribution to the muon g-2**/
    /**
     * @brief Calculates amplitudes for \f$ (g-2)_{\mu} \f$ at approximate two-loop from \cite Broggio:2014mna.
     * @return 
     */    
    virtual double gminus2muNLO();
    
    
    
    
    
      /** Calculates the bosonic NLO contribution to the muon g-2**/
    /**
     * @brief Calculates amplitudes for \f$ (g-2)_{\mu} \f$ at approximate two-loop from \cite Cherchiglia:2016eui.
     * @return 
     */    
    //virtual double gminus2muNLOF();
    //WE'RE REMOVING THIS DEFINITIONS AND SUBSTITUTING THEM FROM THOSE OF 1502.04199
    
     /** Calculates the bosonic NLO contribution to the muon g-2**/
    /**
     * @brief Calculates amplitudes for \f$ (g-2)_{\mu} \f$ at approximate two-loop from \cite Cherchiglia:2016eui.
     * @return 
     */    
    //virtual double gminus2muNLOB();
    //WE'RE REMOVING THIS DEFINITIONS AND SUBSTITUTING THEM FROM THOSE OF 1502.04199
    
    
    
    

      /** Calculates the square root of a negative number**/
    /**
     * @brief Calculates the square root of a negative number
     * @return  square root of a negative number
     */    
    virtual gslpp::complex negsquareroot(double x);
    
        /** Calculates the power root of a negative number**/
    /**
     * @brief Calculates the power root of a negative number
     * @return  power root of a negative number
     */    
    virtual gslpp::complex negpow(double basis, double exp);
    
           /** Calculates the log of a negative number**/
    /**
     * @brief Calculates the log of a negative number
     * @return  log of a negative number
     */  
    
    virtual gslpp::complex neglog(gslpp::complex argument);

    
    
    
    
    
    
    /** Calculates the function of Eq. (68) of 1607.06292**/
    /**
     * @brief Calculates the function of Eq. (68) of 1607.06292
     * @return  function of Eq. (68) of 1607.06292
     */    
    //virtual gslpp::complex TF(double m1, double m2, double m3);
    //WE'RE REMOVING THIS DEFINITIONS AND SUBSTITUTING THEM FROM THOSE OF 1502.04199
    
    
    
    
    
    
    
    
    /**
     * @return GeneralTHDM Wilson coefficients for \f$ B_s \to \bar{B_s}\f$ according to @cite Geng:1988bq, @cite Deschamps:2009rh
     */
    virtual  std::vector<WilsonCoefficient>& CMdbs2();
//    virtual  std::vector<WilsonCoefficient>& CMdbsp2();

    /**
     * @return GeneralTHDM Wilson coefficient for \f$ B \to \tau \nu \f$ from @cite Hou:1992sy
     */
    virtual  std::vector<WilsonCoefficient>& CMbtaunu(QCD::meson meson_i);
    
    /**
     * 
     * @return Wilson coefficient for \f$ D \rightarrow \lepton \nu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMcleptonnu(QCD::meson meson_i, QCD::lepton lepton_i) ;
    
    
    virtual std::vector<WilsonCoefficient>& CMBMll(QCD::lepton lepton);
    
    
    /**
     * @return C10 Wilson coefficient of  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual double C10Bll(double xt, double xHp, gslpp::complex su);
    
    /**
     * @return Box CS Wilson coefficient of  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual  gslpp::complex CSboxBll(double xt, double xHp, gslpp::complex su, gslpp::complex sd, gslpp::complex sl);
    
       /**
     * @return Box CP Wilson coefficient of  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual  gslpp::complex CPboxBll(double xt, double xHp, gslpp::complex su, gslpp::complex sd, gslpp::complex sl);
    
          /**
     * @return Z-penguin in the unitary gauge CP Wilson coefficient of  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual  gslpp::complex CPZUBll(double xt, double xHp, double sW2, gslpp::complex su, gslpp::complex sd);
    
   
    /**
     * @return f1 needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual  double f1(double xHp, double xt);
    
        
    /**
     * @return f2 needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual  double f2(double xHp, double xt);
    
        
    /**
     * @return f3 needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual  double f3(double xHp, double xt);
    
        
    /**
     * @return f4 needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual  double f4(double xHp, double xt);
    
        
    /**
     * @return f5 needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual  double f5(double xHp, double xt);
    
        
    /**
     * @return f6 needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual  double f6(double xHp, double xt);
    
        
    /**
     * @return f7 needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual  double f7(double xHp, double xt);
    
        
    /**
     * @return f8 needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual  double f8(double xHp, double xt);
   
    
        
    /**
     * @return f9 needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual  double f9(double xHp, double xt);
    
      /**
     * @return f10 needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual  double f10(double xHp, double xt);
    
    
     /**
     * @return g0 needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual gslpp::complex  g0(double xHp, double xt, gslpp::complex su, gslpp::complex sd);
    
      /**
     * @return g1a needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual gslpp::complex  g1a(double xHp, double xt, gslpp::complex su, gslpp::complex sd);
    
      /**
     * @return g2a needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual gslpp::complex  g2a(double xHp, double xt, gslpp::complex su, gslpp::complex sd);
    
      /**
     * @return g3a needed to calculate  \f$ B_q \to l \bar{l}\f$ according to @cite Li:2014fea
     */
    virtual gslpp::complex  g3a(double xHp, double xt, gslpp::complex su, gslpp::complex sd);
    
    
    virtual gslpp::complex lambdaHHphi(double lambda3, double Relambda7,double Imlambda7, double Ri1, double Ri2, double Ri3 );
    
    virtual gslpp::complex CphiU(double xHp, double xt, double vev, double xphi, double mu, double Ri1, double Ri2, double Ri3, double mHi_2, double lambda3, double Relambda7,double Imlambda7, gslpp::complex su, gslpp::complex sd);
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return GeneralTHDM Wilson coefficients, Misiak basis, for \f$ B \rightarrow X_{s} \gamma, l^{+} l^{-} \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMbsg() ;
    
    /**
     * 
     * @param i int, flag for the caching
     * @param sigmau
     * @param sigmad
     * @param mt top mass
     * @param mhp charged Higgs mass
     * @param mu matching scale
     * @param order
     * @return return the value of the Wilson coefficients for \f$ B \rightarrow X_{s} \gamma, l^{+} l^{-} \f$
     */
    gslpp::complex setWCbsg (int i, gslpp::complex sigmau, gslpp::complex sigmad, double mt, double mhp, double mu, orders order);
    
     const Polylogarithms getPolyLog() const
       {
           return PolyLog;
       }

    void updateGTHDMParameters();

private:
    const GeneralTHDM & myGTHDM;

    gslpp::matrix<gslpp::complex> myCKM;
    WilsonCoefficient mcdbs2, mcbtaunu, mccleptonnu, mcBMll, mcbsg, mcgminus2mu, mcbsmm;

    double GF, mMU;
    gslpp::complex CWbsgArrayLO[8], CWbsgArrayNLO[8], CWbsgArrayNNLO[8];
    double mtbsg, mhpbsg, mubsg; // caching
    gslpp::complex su, sd, sl; // caching

    const Polylogarithms PolyLog;
};

#endif	/* GENERALTHDMMATCHING_H */

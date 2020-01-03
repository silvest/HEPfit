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
     * @brief Calculates amplitudes for \f$ (g-2)_{\mu} \f$ at one loop from \cite Broggio:2014mna.
     * @return 
     */    
    virtual double gminus2muLO();

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
    virtual double gminus2muNLOF();
    
     /** Calculates the bosonic NLO contribution to the muon g-2**/
    /**
     * @brief Calculates amplitudes for \f$ (g-2)_{\mu} \f$ at approximate two-loop from \cite Cherchiglia:2016eui.
     * @return 
     */    
    virtual double gminus2muNLOB();

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
    virtual gslpp::complex TF(double m1, double m2, double m3);
    
    /**
     * @return GeneralTHDM Wilson coefficients for \f$ B_s \to \bar{B_s}\f$ according to @cite Geng:1988bq, @cite Deschamps:2009rh
     */
    virtual  std::vector<WilsonCoefficient>& CMdbs2();
//    virtual  std::vector<WilsonCoefficient>& CMdbsp2();

    /**
     * @return GeneralTHDM Wilson coefficient for \f$ B \to \tau \nu \f$ from @cite Hou:1992sy
     */
    virtual  std::vector<WilsonCoefficient>& CMbtaunu(QCD::meson meson_i);
    
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
    WilsonCoefficient mcdbs2, mcbtaunu, mcBMll, mcbsg, mcgminus2mu, mcbsmm;

    double GF, mMU;
    gslpp::complex CWbsgArrayLO[8], CWbsgArrayNLO[8], CWbsgArrayNNLO[8];
    double mtbsg, mhpbsg, mubsg; // caching
    gslpp::complex su, sd, sl; // caching

    const Polylogarithms PolyLog;
};

#endif	/* GENERALTHDMMATCHING_H */

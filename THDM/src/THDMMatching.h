/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMMATCHING_H
#define	THDMMATCHING_H

#include "gslpp.h"
#include "StandardModelMatching.h"

class THDM;

/**
 * @class THDMMatching
 * @ingroup THDM
 * @brief A class for the Wilson coefficients in the %THDM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details At the moment, this includes only the @f$B_s@f$ mass difference and the decay @f$B\to \tau \nu@f$.
 */
class THDMMatching : public StandardModelMatching {
public:
    THDMMatching(const THDM & THDM_i);

    /**
     * @return THDM Wilson coefficients for \f$ B_s \to \bar{B_s}\f$ according to @cite Geng:1988bq, @cite Deschamps:2009rh
     */
    virtual  std::vector<WilsonCoefficient>& CMdbs2();

    /**
     * @return THDM Wilson coefficient for \f$ B \to \tau \nu \f$ from @cite Hou:1992sy
     */
    virtual  std::vector<WilsonCoefficient>& CMbtaunu(QCD::meson meson_i);
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return THDM Wilson coefficients, Misiak basis, for \f$ B \rightarrow X_{s} \gamma, l^{+} l^{-} \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMbsg() ;
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow X_{s} \gamma, l^{+} l^{-} \f$
     */
    virtual   std::vector<WilsonCoefficient>& CMprimebsg() ;
    
    /**
     * 
     * @param i int, flag for the caching
     * @param tan \f$ \tan(\beta) \f$
     * @param mt top mass
     * @param mhp charged Higgs mass
     * @param mu matching scale
     * @param order
     * @return return the value of the wilson coefficients for \f$ B \rightarrow X_{s} \gamma, l^{+} l^{-} \f$
     */
    double setWCbsg (int i, double tan, double mt, double mhp, double mu, orders order);

    /** Calculates the muon g-2 at LO**/
    /**
     * @brief Calculates amplitudes for \f$ (g-2)_{\mu} \f$ at one loop from \cite Broggio:2014mna.
     * @return 
     */    
    virtual double gminus2muLO();

    /** Calculates the NLO contribution to the muon g-2**/
    /**
     * @brief Calculates amplitudes for \f$ (g-2)_{\mu} \f$ at approximate two-loop from \cite Broggio:2014mna.
     * @return \f$ (g-2)_{\mu} \f$ at LO
     */    
    virtual double gminus2muNLO();

    /**
     * @brief Wilson coefficient for \f$ (g-2)_{\mu} \f$.
     * @return \f$ (g-2)_{\mu} \f$ at NLO
     */
    virtual std::vector<WilsonCoefficient>& CMgminus2mu();

private:
    const THDM & myTHDM;
    gslpp::matrix<gslpp::complex> myCKM;
    WilsonCoefficient mcdbs2, mcbtaunu, mcbsg, mcprimebsg, mcgminus2mu;
    
    double CWbsgArrayLO[8], CWbsgArrayNLO[8], CWbsgArrayNNLO[8];
    double tanbsg, mtbsg, mhpbsg, mubsg; // caching
//    double gscalar(double r);
//    double gpseudoscalar(double r);

};

#endif	/* THDMMATCHING_H */

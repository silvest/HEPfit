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
    virtual  std::vector<WilsonCoefficient>& CMbtaunu();
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return THDM Wilson coefficients, Misiak basis, for \f$ B \rightarrow X_{s} \gamma, l^{+} l{-} \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMbsg() ;
    
    /** 
     * 
     * @brief operator basis: current current; qcd penguins; 
     * magnetic and chromomagnetic penguins; semileptonic 
     * @return Wilson coefficients, Misiak basis, for \f$ B \rightarrow X_{s} \gamma, l^{+} l{-} \f$
     */
    virtual   std::vector<WilsonCoefficient>& CMprimebsg() ;
    
    /**
     * 
     * @param i int, flag for the caching
     * @param x the square ratio between top mass and charged Higgs mass
     * @param tan \f$ \tan(\beta) \f$
     * @param mu matching scale
     * @param order
     * @return return the value of the wilson coefficients for \f$ B \rightarrow X_{s} \gamma, l^{+} l{-} \f$
     */
    double setWCbsg (int i, double x, double tan, double mu, orders order);

private:
    const THDM & myTHDM;
    gslpp::matrix<gslpp::complex> myCKM;
    WilsonCoefficient mcdbs2, mcbtaunu, mcbsg, mcprimebsg;
    
    double CWbsgArrayLO[8], CWbsgArrayNLO[8];
    double tanbsg, xbsg, mubsg; // caching

};

#endif	/* THDMMATCHING_H */

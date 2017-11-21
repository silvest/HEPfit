/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GENERALTHDMMATCHING_H
#define	GENERALTHDMMATCHING_H

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

    /**
     * @return GeneralTHDM Wilson coefficients for \f$ B_s \to \bar{B_s}\f$ according to @cite Geng:1988bq, @cite Deschamps:2009rh
     */
    virtual  std::vector<WilsonCoefficient>& CMdbs2();

    /**
     * @return GeneralTHDM Wilson coefficient for \f$ B \to \tau \nu \f$ from @cite Hou:1992sy
     */
    virtual  std::vector<WilsonCoefficient>& CMbtaunu();
    
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

    /**
     *
     * @brief Updates to new GTHDM parameter set.
     * @return
     */
    void updateGTHDMParameters();

private:
    const GeneralTHDM & myGTHDM;

    gslpp::matrix<gslpp::complex> myCKM;
    WilsonCoefficient mcdbs2, mcbtaunu, mcbsg, mcgminus2mu;

    double GF, mMU;
    gslpp::complex CWbsgArrayLO[8], CWbsgArrayNLO[8], CWbsgArrayNNLO[8];
    double mtbsg, mhpbsg, mubsg; // caching
    gslpp::complex su, sd; // caching

};

#endif	/* GENERALTHDMMATCHING_H */

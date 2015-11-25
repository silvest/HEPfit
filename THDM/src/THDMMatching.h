/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
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
 * @brief A class for the Wilson coefficients in the THDM.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details At the moment, this includes only the @f$B_s@f$ mass difference and the decay @f$B\to \tau \nu@f$.
 */
class THDMMatching : public StandardModelMatching {
public:
    THDMMatching(const THDM & THDM_i);

    /**
     * @return THDM Wilson coefficients for \f$ B_s \to \bar{B_s}\f$
     */
    virtual  std::vector<WilsonCoefficient>& CMdbs2();

    /**
     * @return THDM Wilson coefficient for \f$ B \to \tau \nu \f$
     */
    virtual  std::vector<WilsonCoefficient>& CMbtaunu();

private:
    const THDM & myTHDM;
    gslpp::matrix<gslpp::complex> myCKM;

    WilsonCoefficient mcdbs2, mcbtaunu;

};

#endif	/* THDMMATCHING_H */

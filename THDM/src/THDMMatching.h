/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef THDMMATCHING_H
#define	THDMMATCHING_H

#include <gslpp.h>
#include "StandardModelMatching.h"

class THDM;

/**
 * @class THDMMatching
 * @ingroup THDM
 * @brief A class for the matching in the THDM. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class THDMMatching : public StandardModelMatching {
public:
    THDMMatching(const THDM & THDM_i);
    
    /**
     *
     * @brief Updates to new THDM parameter sets.
     * @return
     */
    
    void updateTHDMParameters();

private:
    const THDM & myTHDM;
    gslpp::matrix<gslpp::complex> myCKM;

    double tanb;
    double v;
    double v1;
    double v2;
    double gW;
};

#endif	/* THDMMATCHING_H */


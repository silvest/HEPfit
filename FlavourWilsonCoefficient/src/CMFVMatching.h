/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CMFVMATCHING_H
#define	CMFVMATCHING_H

#include "gslpp.h"
#include "StandardModelMatching.h"

class CMFV;

/**
 * @class CMFVMatching
 * @ingroup CMFV
 * @brief A class for the matching in the CMFV. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class CMFVMatching : public StandardModelMatching {
public:
    CMFVMatching(const CMFV & CMFV_i);
    
    double S0(double, double) const;
    
private:
    const CMFV & myCMFV;
};

#endif	/* CMFVMATCHING_H */


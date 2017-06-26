/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MODELMATCHING_H
#define	MODELMATCHING_H

#include"WilsonCoefficient.h"
#include "QCD.h"
#include <vector>
#include <string>

/**
 * @class ModelMatching
 * @ingroup StandardModel
 * @brief A class for a template of model matching. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class ModelMatching {
public:
    
    virtual ~ModelMatching() {};
    
    virtual std::vector<WilsonCoefficient> CMDF1(std::string blocks, unsigned int nops, schemes scheme, orders order) = 0;
    
    virtual std::vector<WilsonCoefficient>& CMdbd2() = 0;
    virtual std::vector<WilsonCoefficient>& CMdbs2() = 0;
    virtual std::vector<WilsonCoefficient>& CMdd2() = 0;
};

#endif	/* MODELMATCHING_H */

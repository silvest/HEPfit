/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MODELMATCHING_H
#define	MODELMATCHING_H

#include"WilsonCoefficient.h"
#include <vector>

/**
 * @class ModelMatching
 * @ingroup StandardModel
 * @brief A class for a template of model matching. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class ModelMatching {
public:
    
    virtual std::vector<WilsonCoefficient>& CMdbd2() = 0;
    virtual std::vector<WilsonCoefficient>& CMdbs2() = 0;
    virtual std::vector<WilsonCoefficient>& CMdd2() = 0;

    virtual std::vector<WilsonCoefficient>& CMbnlep( const int a) = 0;
    virtual std::vector<WilsonCoefficient>& CMbnlepCC( const int a) = 0;
    virtual std::vector<WilsonCoefficient>& CMbsg()= 0;
    virtual std::vector<WilsonCoefficient>& CMd1()= 0;
    virtual std::vector<WilsonCoefficient>& CMd1Buras()= 0;
    
};

#endif	/* MODELMATCHING_H */


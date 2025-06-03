/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MODELMATCHING_H
#define	MODELMATCHING_H

#include"WilsonCoefficient.h"
#include"WilsonCoefficientNew.h"
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
    
    virtual std::vector<WilsonCoefficientNew>& CMDF1(std::string blocks, unsigned int nops) = 0;
    
    virtual std::vector<WilsonCoefficient>& CMdbd2() = 0;
    virtual std::vector<WilsonCoefficient>& CMdbs2() = 0;
    virtual std::vector<WilsonCoefficient>& CMdd2() = 0;

    virtual std::vector<WilsonCoefficient>& CMbnlep( int a) = 0;
    virtual std::vector<WilsonCoefficient>& CMbnlepCC( const int a) = 0;
    virtual std::vector<WilsonCoefficient>& CMbsg() = 0;
    virtual std::vector<WilsonCoefficient>& CMprimebsg() = 0;
    virtual std::vector<WilsonCoefficient>& CMBMll(QCD::lepton lepton) = 0;
    virtual std::vector<WilsonCoefficient>& CMprimeBMll(QCD::lepton lepton) = 0;
    virtual std::vector<WilsonCoefficient>& CMBXsnn(QCD::lepton lepton) = 0;
    virtual std::vector<WilsonCoefficient>& CMd1() = 0;
    virtual std::vector<WilsonCoefficient>& CMd1Buras() = 0;
    
};

#endif	/* MODELMATCHING_H */

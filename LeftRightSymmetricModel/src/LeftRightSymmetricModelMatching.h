/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEFTRIGHTSYMMETRICMODELMATCHING_H
#define	LEFTRIGHTSYMMETRICMODELMATCHING_H

#include "gslpp.h"
#include "StandardModelMatching.h"

class LeftRightSymmetricModel;

/**
 * @class LeftRightSymmetricModelMatching
 * @ingroup LeftRightSymmetricModel
 * @brief A class for the matching in the LeftRightSymmetricModel. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class LeftRightSymmetricModelMatching : public StandardModelMatching {
public:
    LeftRightSymmetricModelMatching(const LeftRightSymmetricModel & LeftRightSymmetricModel_i);
    
    /**
     *
     * @brief Updates to new LeftRightSymmetricModel parameter sets.
     * @return
     */
    
    void updateLeftRightSymmetricModelParameters();

    gslpp::complex setWCbsg (int i, double mu, orders order);

    virtual std::vector<WilsonCoefficient>& CMbsg();
    
    std::vector<WilsonCoefficient>& CMprimebsg();
    
    std::vector<WilsonCoefficient>& CMBMll(QCD::lepton lepton);
    
    std::vector<WilsonCoefficient>& CMprimeBMll(QCD::lepton lepton);

private:
    const LeftRightSymmetricModel & myLeftRightSymmetricModel;
    gslpp::matrix<gslpp::complex> myCKM;
    gslpp::matrix<gslpp::complex> myCKMR;
    
    WilsonCoefficient mcbsg, mcprimebsg;
    WilsonCoefficient mcBMll, mcprimeBMll;
    std::vector<WilsonCoefficient> vmcbsg, vmcprimebsg;
    std::vector<WilsonCoefficient> vmcBMll, vmcprimeBMll;
    gslpp::complex CWbsgArrayLO[8];
    double Muw,Mut,mW,mtop,mbottom,vev,gW;
    double mWR,mH2psq,xi,alpha;
};

#endif	/* LEFTRIGHTSYMMETRICMODELMATCHING_H */


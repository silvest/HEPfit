/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */



#include "HeffDF1_Plepnu.h"
#include "StandardModel.h"


HeffDF1_Plepnu::HeffDF1_Plepnu(const StandardModel & SM) 
: model(SM), coeffuleptonnu(4, NDR, LO)
{
  
    
}

HeffDF1_Plepnu::~HeffDF1_Plepnu() {
}


 /******************************************************************************
 * evolution Wilson Coefficien \pi -> lepton nu  
 * The Evolution is not yet implemented. In the SM is zero 
 * LEFT basis. The WC are written in the LEFT basis of arxiv:1709.04486
 *the expressions can be found in arxiv:1706.00410 and arxiv:1605.07114
 *in a similar basis
                                                             *
 ******************************************************************************/
gslpp::vector<gslpp::complex>** HeffDF1_Plepnu::ComputeCoeffuleptonnu(QCD::meson meson_i, QCD::lepton lepton_i) 
{
    const std::vector<WilsonCoefficient>& mcu = model.getMatching().CMuleptonnu(meson_i, lepton_i);
    coeffuleptonnu.resetCoefficient();
    orders ordDF1 = coeffuleptonnu.getOrder();
    for (unsigned int i = 0; i < mcu.size(); i++){
        for (int j = LO; j <= ordDF1; j++){
            coeffuleptonnu.setCoeff(*coeffuleptonnu.getCoeff(orders(j))
                                     + *mcu[i].getCoeff(orders(j)), orders(j));
        }
    }
     return coeffuleptonnu.getCoeff(); 
}
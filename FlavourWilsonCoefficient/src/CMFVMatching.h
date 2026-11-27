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
 * @ingroup FlavourWilsonCoefficent
 * @brief A class for the matching in the CMFV.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details The shift of S enters through S0(x_t, x_t). The shifts
 * of X, Y, Z and C7eff enter as one additional leading-order Wilson-coefficient block per
 * Delta F = 1 sector, with the normalisation of the corresponding SM matching.
 */
class CMFVMatching : public StandardModelMatching {
public:
    CMFVMatching(const CMFV & CMFV_i);

    const double S0(double, double) const;

    // Delta F = 1 sectors
    std::vector<WilsonCoefficient>& CMkpnn();
    std::vector<WilsonCoefficient>& CMBXsnn(QCD::lepton lepton);
    std::vector<WilsonCoefficient>& CMBXdnn();
    std::vector<WilsonCoefficient>& CMbsmm();
    std::vector<WilsonCoefficient>& CMbdmm();
    std::vector<WilsonCoefficient>& CMbsg();
    std::vector<WilsonCoefficient>& CMBMll(QCD::lepton lepton);
    std::vector<WilsonCoefficientNew>& CMDF1(std::string blocks, unsigned int nops);

private:
    const CMFV & myCMFV;

    // 
    void zeroAll(WilsonCoefficient& w, unsigned int n, bool qed);
    WilsonCoefficient npkpnn, npbsnn, npbdnn, npbsmm, npbdmm, npbsg, npBMll;
    std::vector<WilsonCoefficient> vkpnn, vbsnn, vbdnn, vbsmm, vbdmm, vbsg, vBMll;
    std::vector<WilsonCoefficientNew> vDF1;
};

#endif	/* CMFVMATCHING_H */


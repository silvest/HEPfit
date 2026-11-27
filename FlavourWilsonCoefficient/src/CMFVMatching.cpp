/*
 * Copyright (C) 2016 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#include "CMFVMatching.h"
#include "CMFV.h"
#include <stdexcept>

CMFVMatching::CMFVMatching(const CMFV & CMFV_i) :
    StandardModelMatching(CMFV_i),
    myCMFV(CMFV_i),
    npkpnn(2, NDR, NLO, NLO_QED11),
    npbsnn(2, NDR, NLO, NLO_QED11), npbdnn(2, NDR, NLO, NLO_QED11),
    npbsmm(8, NDR, NNLO, NLO_QED22), npbdmm(8, NDR, NNLO, NLO_QED22),
    npbsg(8, NDR, NNLO), npBMll(13, NDR, NLO) {};

const double CMFVMatching::S0(double x1, double x2) const {
    if (x1 > .02 && x2 > .02)
        return (myCMFV.getFtt() + StandardModelMatching::S0(x1, x2));
    return StandardModelMatching::S0(x1,x2);
}

/*******************************************************************************
 * Delta F = 1 master functions X, Y, Z and C7eff.
 * Each method returns the SM blocks plus one leading-order NP block that carries
 * the CMFV shift with the normalisation of the SM matching.
 ******************************************************************************/

void CMFVMatching::zeroAll(WilsonCoefficient& w, unsigned int n, bool qed)
{
    for (unsigned int j = 0; j < n; j++) {
        for (int o = LO; o <= w.getOrder(); o++) w.setCoeff(j, 0., orders(o));
        if (qed) {
            w.setCoeff(j, 0., LO_QED);
            if (w.getOrder_qed() >= NLO_QED11) w.setCoeff(j, 0., NLO_QED11);
            if (w.getOrder_qed() >= NLO_QED22) {
                w.setCoeff(j, 0., NLO_QED02);
                w.setCoeff(j, 0., NLO_QED21);
                w.setCoeff(j, 0., NLO_QED12);
                w.setCoeff(j, 0., NLO_QED22);
            }
        }
    }
}

// K -> pi nu nu. SM: 4 G_F/sqrt2 * alpha(M_Z)/(2 pi s_W^2) * lambda_t * X_0(x_t)
std::vector<WilsonCoefficient>& CMFVMatching::CMkpnn()
{
    vkpnn.clear();
    for (const WilsonCoefficient& w : StandardModelMatching::CMkpnn()) vkpnn.push_back(w);
    zeroAll(npkpnn, 2, true);
    double co = 4. * myCMFV.getGF() / sqrt(2.) * myCMFV.alphaMz() / 2. / M_PI / myCMFV.sW2_ND();
    npkpnn.setMu(myCMFV.getMut());
    npkpnn.setCoeff(0, co * myCMFV.getCKM().computelamt() * myCMFV.getdX(), LO);
    vkpnn.push_back(npkpnn);
    return (vkpnn);
}

// b -> s nu nu, b -> d nu nu. SM: -X_0(x_t)/s_W^2 (CKM factors in the observable classes)
std::vector<WilsonCoefficient>& CMFVMatching::CMBXsnn(QCD::lepton lepton)
{
    vbsnn.clear();
    for (const WilsonCoefficient& w : StandardModelMatching::CMBXsnn(lepton)) vbsnn.push_back(w);
    zeroAll(npbsnn, 2, true);
    npbsnn.setMu(myCMFV.getMut());
    npbsnn.setCoeff(0, -myCMFV.getdX() / myCMFV.sW2_ND(), LO);
    vbsnn.push_back(npbsnn);
    return (vbsnn);
}

std::vector<WilsonCoefficient>& CMFVMatching::CMBXdnn()
{
    vbdnn.clear();
    for (const WilsonCoefficient& w : StandardModelMatching::CMBXdnn()) vbdnn.push_back(w);
    zeroAll(npbdnn, 2, true);
    npbdnn.setMu(myCMFV.getMut());
    npbdnn.setCoeff(0, -myCMFV.getdX() / myCMFV.sW2_ND(), LO);
    vbdnn.push_back(npbdnn);
    return (vbdnn);
}

// B_q -> mu mu. In the basis of StandardModelMatching::setWCBsmm the O(alpha_e) coefficient [7] = -Y_0
// V_tb^* V_tq is C10. The improved calculation used by Bsmumu/Bdmumu (Bobeth, Gorbahn, Stamou: NNLO QCD, NLO EW)
// is valid for new physics in C10 only, so the CMFV shift enters only there: [7] += -dY V_tb^* V_tq.
// (dZ would enter the vector coefficient [6], which reaches the rate only through O(alpha_e) mixing into C10,
// a ~1% effect for dZ = 1; it is deliberately not included.)
std::vector<WilsonCoefficient>& CMFVMatching::CMbsmm()
{
    vbsmm.clear();
    for (const WilsonCoefficient& w : StandardModelMatching::CMbsmm()) vbsmm.push_back(w);
    zeroAll(npbsmm, 8, true);
    gslpp::complex lam = myCMFV.getCKM().getV_tb().conjugate() * myCMFV.getCKM().getV_ts();
    npbsmm.setMu(myCMFV.getMuw());
    npbsmm.setCoeff(7, lam * (-myCMFV.getdY()), NLO_QED11); // C10 only, see above
    vbsmm.push_back(npbsmm);
    return (vbsmm);
}

std::vector<WilsonCoefficient>& CMFVMatching::CMbdmm()
{
    vbdmm.clear();
    for (const WilsonCoefficient& w : StandardModelMatching::CMbdmm()) vbdmm.push_back(w);
    zeroAll(npbdmm, 8, true);
    gslpp::complex lam = myCMFV.getCKM().getV_tb().conjugate() * myCMFV.getCKM().getV_td();
    npbdmm.setMu(myCMFV.getMuw());
    npbdmm.setCoeff(7, lam * (-myCMFV.getdY()), NLO_QED11); // C10 only, see above
    vbdmm.push_back(npbdmm);
    return (vbdmm);
}

// b -> s gamma: shift of C7eff at muNP
std::vector<WilsonCoefficient>& CMFVMatching::CMbsg()
{
    vbsg.clear();
    for (const WilsonCoefficient& w : StandardModelMatching::CMbsg()) vbsg.push_back(w);
    zeroAll(npbsg, 8, false);
    npbsg.setMu(myCMFV.getmuNP());
    npbsg.setCoeff(6, myCMFV.getdC7(), LO);
    vbsg.push_back(npbsg);
    return (vbsg);
}

// b -> s l l: shifts of C7eff, C9 = dY/s_W^2 - 4 dZ, C10 = -dY/s_W^2 at muNP
std::vector<WilsonCoefficient>& CMFVMatching::CMBMll(QCD::lepton lepton)
{
    vBMll.clear();
    for (const WilsonCoefficient& w : StandardModelMatching::CMBMll(lepton)) vBMll.push_back(w);
    zeroAll(npBMll, 13, false);
    double sw2 = myCMFV.sW2_ND();
    npBMll.setMu(myCMFV.getmuNP());
    npBMll.setCoeff(6, myCMFV.getdC7(), LO);
    npBMll.setCoeff(8, myCMFV.getdY() / sw2 - 4. * myCMFV.getdZ(), LO);
    npBMll.setCoeff(9, -myCMFV.getdY() / sw2, LO);
    vBMll.push_back(npBMll);
    return (vBMll);
}

// Delta F = 1 in the Huber et al. basis used by HeffDF1 (e.g. B -> Xs l l, blocks "CPMLQB").
// Block sizes: C 2, P 4, M 2 (C7, C8), L 2 (C9, C10), Q 4, B 1. As in StandardModelMatching::mc_L,
// C9 and C10 are stored at order (QCD1, QED1) with a factor alpha_e(Muw)/(4 pi) and the on-shell s_W^2;
// C7 at (QCD0, QED0) in the standard normalisation. The C9/C10 shifts are matched at Muw like the SM,
// the C7 shift at muNP (as in CMbsg and CMBMll).
std::vector<WilsonCoefficientNew>& CMFVMatching::CMDF1(std::string blocks, unsigned int nops)
{
    vDF1.clear();
    for (const WilsonCoefficientNew& w : StandardModelMatching::CMDF1(blocks, nops)) vDF1.push_back(w);
    /* Copied, not referenced: the push_back calls below can reallocate vDF1. */
    const schemes scheme = vDF1.front().getScheme();
    const qcd_orders orderQCD = vDF1.front().getOrder_QCD();
    const qed_orders orderQED = vDF1.front().getOrder_QED();

    const std::string names = "CPMLQB";
    const unsigned int size[6] = {2, 4, 2, 2, 4, 1};
    int iM = -1, iL = -1;
    unsigned int tot = 0;
    for (unsigned int b = 0; b < 6; b++)
        if (blocks.find(names[b]) != std::string::npos) {
            if (names[b] == 'M') iM = tot;
            if (names[b] == 'L') iL = tot;
            tot += size[b];
        }
    if (tot != nops) throw std::runtime_error("CMFVMatching::CMDF1(): unexpected block layout");

    if (iL >= 0) {
        double muw = myCMFV.getMuw();
        double ale = myCMFV.Ale(muw, FULLNLO) / 4. / M_PI;
        double sw2 = myCMFV.sW2();
        WilsonCoefficientNew np(nops, scheme, orderQCD, orderQED);
        np.setMu(muw);
        np.setCoeff(iL, ale * (myCMFV.getdY() / sw2 - 4. * myCMFV.getdZ()), QCD1, QED1);
        np.setCoeff(iL + 1, ale * (-myCMFV.getdY() / sw2), QCD1, QED1);
        vDF1.push_back(np);
    }
    if (iM >= 0) {
        WilsonCoefficientNew np(nops, scheme, orderQCD, orderQED);
        np.setMu(myCMFV.getmuNP());
        np.setCoeff(iM, myCMFV.getdC7(), QCD0, QED0);
        vDF1.push_back(np);
    }
    return (vDF1);
}

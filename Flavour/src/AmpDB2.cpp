/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AmpDB2.h"
#include "EvolDF2.h"
#include "HeffDF2.h"

AmpDB2::AmpDB2(const StandardModel& SM_i)
: mySM(SM_i)
{
    mySM.initializeBParameter("BBs");
    mySM.initializeBParameter("BBd");
}

gslpp::complex AmpDB2::RBs(orders order)
{
    mySM.getFlavour().getHDF2().getCoeffBs().getOrder();
    if (mySM.getFlavour().getHDF2().getCoeffBs().getOrder() < order % 3)
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed");

    gslpp::vector<gslpp::complex> ** allcoeff_SM = mySM.getFlavour().ComputeCoeffBs(
            mySM.getBBs().getMu(),
            mySM.getBBs().getScheme(), true);
    
    C_1_SM = ((*(allcoeff_SM[LO]))(0) + (*(allcoeff_SM[NLO]))(0));
    
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffBs(
            mySM.getBBs().getMu(),
            mySM.getBBs().getScheme());

    gslpp::vector<double> me(mySM.getBBs().getBpars());
    double MBs = mySM.getMesons(QCD::B_S).getMass();
    double Mb = mySM.getQuarks(QCD::BOTTOM).getMass();
    double Ms = mySM.Mrun(mySM.getBBs().getMu(),
                mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double KBs = MBs/(Mb+Ms)*MBs/(Mb+Ms);
    double Fbs = mySM.getMesons(QCD::B_S).getDecayconst();
    me(0) *= 1./3.*MBs*Fbs*Fbs;
    me(1) *= -5./24.*KBs*MBs*Fbs*Fbs;
    me(2) *= 1./24.*KBs*MBs*Fbs*Fbs;
    me(3) *= 1./4.*KBs*MBs*Fbs*Fbs;
    me(4) *= 1./12.*KBs*MBs*Fbs*Fbs;

    /*std::cout << "low scale :" << std::endl << std::endl;
    std::cout << "C1_SM :" << C_1_SM << std::endl << std::endl;
    
    std::cout << "C1 :" << ((*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0)) << std::endl;
    std::cout << "C2 :" << ((*(allcoeff[LO]))(1) + (*(allcoeff[NLO]))(1)) << std::endl;
    std::cout << "C3 :" << ((*(allcoeff[LO]))(2) + (*(allcoeff[NLO]))(2)) << std::endl;
    std::cout << "C4 :" << ((*(allcoeff[LO]))(3) + (*(allcoeff[NLO]))(3)) << std::endl;
    std::cout << "C5 :" << ((*(allcoeff[LO]))(4) + (*(allcoeff[NLO]))(4)) << std::endl << std::endl;*/

    switch (order) {
        case FULLNLO:
            return (*(allcoeff[LO]) + *(allcoeff[NLO])) * me /
                    (C_1_SM * me(0));
        case LO:
            return ((*(allcoeff[LO])) * me / HCUT);
        default:
            throw std::runtime_error("RBs::RBs(): order not implemented");
    }
}

gslpp::complex AmpDB2::M12_Bd(orders order) {
    if (mySM.getFlavour().getHDF2().getCoeffBd().getOrder() < order % 3)
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed");

    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffBd(
            mySM.getBBd().getMu(),
            mySM.getBBd().getScheme());

    gslpp::vector<double> me(mySM.getBBd().getBpars());
    double MBd = mySM.getMesons(QCD::B_D).getMass();
    double Mb = mySM.Mrun(mySM.getBBd().getMu(),
            mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
            mySM.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    double Md = mySM.Mrun(mySM.getBBd().getMu(),
            mySM.getQuarks(QCD::DOWN).getMass_scale(),
            mySM.getQuarks(QCD::DOWN).getMass(), FULLNNLO);
    double KBd = MBd / (Mb + Md) * MBd / (Mb + Md);
    double Fb = mySM.getMesons(QCD::B_D).getDecayconst();
    me(0) *= 1. / 3. * MBd * Fb*Fb;
    me(1) *= -5. / 24. * KBd * MBd * Fb*Fb;
    me(2) *= 1. / 24. * KBd * MBd * Fb*Fb;
    me(3) *= 1. / 4. * KBd * MBd * Fb*Fb;
    me(4) *= 1. / 12. * KBd * MBd * Fb*Fb;

#if SUSYFIT_DEBUG & 1
    std::cout << "Bd: me(0) = " << me(0) << std::endl;
#endif
#if SUSYFIT_DEBUG & 2
    std::cout << "coefficient Bd: " << (*(allcoeff[LO]) + *(allcoeff[NLO]))(0) << std::endl;
    std::cout << "M: " << me << std::endl;
    std::cout << "mu : " << mySM.getBBd().getMu() << ", mut: " << mySM.getMut() << ", scheme: " << mySM.getBBd().getScheme() << ", B par.: " << mySM.getBBd().getBpars()(0) << std::endl;
    std::cout << "U (mut): " << (mySM.getFlavour().getHDF2().getUDF2().Df2Evol(mySM.getBBd().getMu(), mySM.getMut(), LO)(0, 0) +
            mySM.getFlavour().getHDF2().getUDF2().Df2Evol(mySM.getBBd().getMu(), mySM.getMut(), NLO)(0, 0)) << std::endl;
#endif

    switch (order) {
        case FULLNLO:
            return ((*(allcoeff[LO]) + *(allcoeff[NLO])) * me / HCUT);
        case LO:
            return ((*(allcoeff[LO])) * me / HCUT);
        default:
            throw std::runtime_error("AmpDB2::AmpBd(): order not implemented");
    }
}

gslpp::complex AmpDB2::M12_Bs(orders order) {
    if (mySM.getFlavour().getHDF2().getCoeffBs().getOrder() < order % 3)
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed");

    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffBs(
            mySM.getBBs().getMu(),
            mySM.getBBs().getScheme());

    gslpp::vector<double> me(mySM.getBBs().getBpars());
    double MBs = mySM.getMesons(QCD::B_S).getMass();
    double Mb = mySM.getQuarks(QCD::BOTTOM).getMass();
    double Ms = mySM.Mrun(mySM.getBBs().getMu(),
            mySM.getQuarks(QCD::STRANGE).getMass_scale(),
            mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double KBs = MBs / (Mb + Ms) * MBs / (Mb + Ms);
    double Fbs = mySM.getMesons(QCD::B_S).getDecayconst();
    me(0) *= 1. / 3. * MBs * Fbs*Fbs;
    me(1) *= -5. / 24. * KBs * MBs * Fbs*Fbs;
    me(2) *= 1. / 24. * KBs * MBs * Fbs*Fbs;
    me(3) *= 1. / 4. * KBs * MBs * Fbs*Fbs;
    me(4) *= 1. / 12. * KBs * MBs * Fbs*Fbs;
#if SUSYFIT_DEBUG & 1
    std::cout << "Bs: me(0) = " << me(0) << std::endl;
#endif

    switch (order) {
        case FULLNLO:
            return ((*(allcoeff[LO]) + *(allcoeff[NLO])) * me / HCUT);
        case LO:
            return ((*(allcoeff[LO])) * me / HCUT);
        default:
            throw std::runtime_error("AmpDB2::AmpBs(): order not implemented");
    }
}

//equation 18
gslpp::vector<gslpp::complex> AmpDB2::c(quark q) {
    gslpp::vector< complex > c(2, 0.);
    switch (q) {
        case d:
            for (int i = 0; i <= 1; i++) {
                c.assign(i,
                        VtbVtd2 * D(uu, i) + 2. * VcbVcd * VtbVtd * (D(uu, i) - D(cu, i))
                        + VcbVcd2 * (D(cc, i) + D(uu, i) - 2. * D(cu, i))
                        );
            }
            break;
        case s:
            for (int i = 0; i <= 1; i++) {
                c.assign(i,
                        VtbVts2 * D(uu, i) + 2. * VcbVcs * VtbVts * (D(uu, i) - D(cu, i))
                        + VcbVcs2 * (D(cc, i) + D(uu, i) - 2. * D(cu, i))
                        );
            }
            break;
        default:
            throw std::runtime_error("AmpDB2::c(quark q, double mu_2): invalid quark index: ");
    }
    return c;
}

void AmpDB2::compute_matrixelements(quark q){
    double Mq;
    double MBq2;
    double KBq;
    double FBq2;
    switch (q) {
        case d:
            Mq = Md;
            MBq2 = MB2;
            FBq2 = mySM.getMesons(QCD::B_D).getDecayconst() * mySM.getMesons(QCD::B_D).getDecayconst();
            break;
        case s:
            Mq = Ms;
            MBq2 = MB_s * MB_s;
            FBq2 = mySM.getMesons(QCD::B_S).getDecayconst() * mySM.getMesons(QCD::B_S).getDecayconst();            
            break;
        default:
            throw std::runtime_error("AmpDB2::compute_matrixelements(quark q): invalid quark index: ");
    }
    KBq = MBq2 / ((Mb + Mq) * (Mb + Mq));
    
//    me.assign(0, 0.);
//    me.assign(1, 0.);
//    me.assign(2, 0.);
//    me.assign(3, 0.);
//    me.assign(4, 0.);
//
//    me_R.assign(0, 0.);
//    me_R.assign(1, 0.); 
//    me_R.assign(2, 0.);
//    me_R.assign(3, 0.);
    
    me = mySM.getBBs().getBpars();
    me(0) *=  8. / 3. * MBq2 * FBq2;
    me(1) *= -5. / 3. * KBq * MBq2 * FBq2;
    me(2) *=  1. / 3. * KBq * MBq2 * FBq2;
    me(3) *=       2. * KBq * MBq2 * FBq2;
    me(4) *=  2. / 3. * KBq * MBq2 * FBq2;
    
    me_R(0) += Mq/Mb * me(0);
    //new lattice results?
    me_R(1) += 1.;
    me_R(2) += 1.;
    me_R(3) += me(2) + 0.5 * me(0) + me(1) - 2. * Mq/Mb * me(4) + me_R(1);
    return;
}

//equation 18
gslpp::complex AmpDB2::delta_1overm(quark q) {
    switch (q) {
        case d:
            compute_deltas_1overm(q);
            return VtbVtd2 * deltas_1overm(uu, d)
                    + 2. * VcbVcd * VtbVtd * (deltas_1overm(uu, d) - deltas_1overm(cu, d))
                    + VcbVcd2 * (deltas_1overm(cc, d) + deltas_1overm(uu, d) - 2. * deltas_1overm(cu, d));
        case s:
            compute_deltas_1overm(q);
            return VtbVts2 * deltas_1overm(uu, s)
                    + 2. * VcbVcs * VtbVts * (deltas_1overm(uu, s) - deltas_1overm(cu, s))
                    + VcbVcs2 * (deltas_1overm(cc, s) + deltas_1overm(uu, s) - 2. * deltas_1overm(cu, s));
        default:
            throw std::runtime_error("AmpDB2::delta_1overm(quark q): invalid quark index: ");
    }
}

//equation 20
void AmpDB2::compute_deltas_1overm(quark q) {
    cache_deltas_1overm[index_deltas(cc, q)] = sqrt1minus4z * ((1 + 2. * z) * (K_2 * (me_R(1) + 2. * me_R(3)) - 2. * K_1 * (me_R(0) + me_R(1)))
                    - 12. * z2 / (1. - 4. * z) * (K_1 * (me_R(1) + 2. * me_R(2)) + 2. * K_2 * me_R(2)));
    cache_deltas_1overm[index_deltas(cu, q)] = oneminusz2 * ((1. + 2. * z) * (K_2 * (me_R(1) + 2. * me_R(3)) - 2. * K_1 * (me_R(0) + me_R(1)))
                    - 6. * z2 / (1. - z) * (K_1 * (me_R(1) + 2. * me_R(2)) + 2. * K_2 * me_R(2)));
    cache_deltas_1overm[index_deltas(uu, q)] =  K_2 * (me_R(1) + 2. * me_R(3)) - 2. * K_1 * (me_R(0) + me_R(1));
    return;
}

double AmpDB2::F0(quarks qq, int k, int i, int j) {
    return cacheF0[indexF(qq, k, i, j)];
}

double AmpDB2::F1(quarks qq, int k, int i, int j) {
    return cacheF1[indexF(qq, k, i, j)];
}

double AmpDB2::F(quarks qq, int k, int i, int j) {
    return F0(qq, k, i, j) + as_4pi * F1(qq, k, i, j);
}

double AmpDB2::P(quarks qq, int k, int i, int j) {
    return cacheP[indexP(qq, k, i, j)];
}

gslpp::complex AmpDB2::D(quarks qq, int k) {
    return cacheD[indexD(qq, k)];
}

gslpp::complex AmpDB2::deltas_1overm(quarks qq, quark q) {
    return cache_deltas_1overm[index_deltas(qq, q)];
}

//indizes for F to save them in an array
int AmpDB2::indexF(quarks qq, int k, int i, int j) {
    return qq * 8 + (k - 1) * 4 + (i - 1) * 2 + (j - 1);
}
//indizes for P to save them in an array: c := cc and u := uu
int AmpDB2::indexP(quarks qq, int k, int i, int j) {
    return qq * 28 + (k - 1) * 14 + (i - 1) * 7 + (j - 2);
}

//indizes for D to save them in arrays
int AmpDB2::indexD(quarks qq, int k) {
    return qq * 2 + (k - 1);
}

//indizes for deltas_1overm to save them in an array
int AmpDB2::index_deltas(quarks qq, quark q) {
    return qq * 2 + (q - 1);
}

//equations 41-42: F0 = A
void AmpDB2::computeF0() {
    cacheF0[indexF(cu, 1, 1, 1)] = 1.5 * (2. - 3. * z + z3);
    cacheF0[indexF(cu, 1, 1, 2)] = 0.5 * (2. - 3. * z + z3);
    cacheF0[indexF(cu, 1, 2, 2)] = 0.5 * oneminusz2 * (1. - z);
    cacheF0[indexF(cu, 2, 1, 1)] = 3. * oneminusz2 * (1. + 2. * z);
    cacheF0[indexF(cu, 2, 1, 2)] = oneminusz2 * (1. + 2. * z);
    cacheF0[indexF(cu, 2, 2, 2)] = -oneminusz2 * (1. + 2. * z);
    cacheF0[indexF(cc, 1, 1, 1)] = 3. * sqrt1minus4z * (1. - z);
    cacheF0[indexF(cc, 1, 1, 2)] = sqrt1minus4z * (1. - z);
    cacheF0[indexF(cc, 1, 2, 2)] = 0.5 * sqrt1minus4z * sqrt1minus4z * sqrt1minus4z;
    cacheF0[indexF(cc, 2, 1, 1)] = 3. * sqrt1minus4z * (1. + 2. * z);
    cacheF0[indexF(cc, 2, 1, 2)] = sqrt1minus4z * (1. + 2. * z);
    cacheF0[indexF(cc, 2, 2, 2)] = -sqrt1minus4z * (1. + 2. * z);
    cacheF0[indexF(uu, 1, 1, 1)] = 3.;
    cacheF0[indexF(uu, 1, 1, 2)] = 1.;
    cacheF0[indexF(uu, 1, 2, 2)] = 0.5;
    cacheF0[indexF(uu, 2, 1, 1)] = 3.;
    cacheF0[indexF(uu, 2, 1, 2)] = 1.;
    cacheF0[indexF(uu, 2, 2, 2)] = -1.;
    for (int k = 1; k < 3; k++) {
        for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1)) {
            std::cout << "enum " << qq << "\n";
            cacheF0[indexF(qq, k, 2, 1)] = cacheF0[indexF(qq, k, 1, 2)];
        }
    }
    //std::cout << "F0 " << indexF(cu, 1, 1, 1) << F0(cu, 1, 1, 1) << "\n";
    return;
}

//equations 43-48: F1 = B
void AmpDB2::computeF1() {
    cacheF1[indexF(cu, 1, 1, 1)] = 109./6. - 37. * z + 1.5 * z2 + 52./3. * z3 + 2. * oneminusz2 * (5. + z) * logx_2 - 4. * oneminusz2 * (5. + 7. * z) * log1minusz -
            2. * z * (10. + 14. * z - 15. * z2) * logz + 8. * (2. - 3. * z + z3) * log1minusz * logz + 16. * (2. - 3. * z + z3) * Dilogz;
    cacheF1[indexF(cu, 2, 1, 1)] = -4./3. * (10. - 33. * z + 54. * z2 - 31. * z3) - 8. * oneminusz2 * (4. + 14. * z - 3. * z2) * log1minusz +
            8. * z * (2. - 23. * z + 21. * z2 - 3. * z3) * logz -
            16. * oneminusz2 * (1. + 2. * z) * (2. * logx_2 - log1minusz * logz - 2. * Dilogz);
    cacheF1[indexF(cu, 1, 1, 2)] = 2. * (
           (-502. + 912. * z - 387. * z2 - 23. * z3) / 36. - oneminusz2 * (17. + 4. * z) * logx_1 + 2./3. * oneminusz2 * (5. + z) * logx_2 -
           oneminusz2 / (12. * z) * (2. + 33. * z + 94. * z2) * log1minusz - z/12. * (80. + 69. * z - 126. * z2) * logz +
           8./3. * (2. - 3. * z + z3) * (log1minusz * logz + 2. * Dilogz)
           );
    cacheF1[indexF(cu, 2, 1, 2)] = 2. * (
           (-130. + 93. * z + 144. * z2 - 107. * z3) / 9. - 2./3. * oneminusz2 / z * (1. + 15. * z + 47. * z2 - 12. * z3) * log1minusz +
           2./3. * z * (8. - 93. * z + 87. * z2 - 12. * z3) * logz -
           8./3. * oneminusz2 * (1. + 2. * z) * (3. * logx_1 + 4. * logx_2 - 2. * log1minusz * logz - 4. * Dilogz)
           );
    cacheF1[indexF(cu, 1, 2, 2)] = -M_PI2/3. * (1. - 5. * z + 4. * z2) + (-136. - 159. * z + 738. * z2 - 443. * z3) / 18. - 2 * oneminusz2 * (5. + 4. * z) * logx_1 +
           2./3. * oneminusz2 * (4. - z) * logx_2 + oneminusz2 / (6. * z) * (7. + 32. * z2 + 3. * z3) * log1minusz -
           z/6. * (62. + 39. * z - 30. * z2 + 3. * z3) * logz + (5. - 3. * z - 18. * z2 + 16. * z3) / 3. * (log1minusz * logz + 2. * Dilogz);
    cacheF1[indexF(cu, 2, 2, 2)] = 8./3. * M_PI2 * (1. + z - 2. * z2) - 28./9. * (5. + 3. * z - 27. * z2 + 19. * z3) - 16. * oneminusz2 * (1. + 2. * z) * logx_1 +
           32./3. * oneminusz2 * (1. + 2. * z) * logx_2 - 4./3. * oneminusz2 / z * (1. - 12. * z - 16. * z2 - 3. * z3) * log1minusz +
           4./3. * z * (2. -  3. * z + 18. * z2 - 3. * z3) * logz + 8./3. * (1. - 3. * z - 6. * z2 + 8. * z3) * (log1minusz * logz + 2. * Dilogz);
    cacheF1[indexF(cc, 1, 1, 1)] = sqrt1minus4z * (109. - 226. * z + 168. * z2) / 6. - (52. - 104. * z - 16. * z2 + 56. * z3) * logsigma +
           2. * (5. - 8. * z) * sqrt1minus4z * logx_2 - 12. * sqrt1minus4z * (3. - 2. * z) * log1minus4z + 4. * (13. - 10. * z) * sqrt1minus4z * logz +
           16. * (1. - 3. * z + 2. * z2) * (3. * log2sigma + 2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2);
    cacheF1[indexF(cc, 2, 1, 1)] = -8./3. * sqrt1minus4z * (5. - 23. * z - 42. * z2) - 16. * (4. - 2. * z - 7. * z2 + 14. * z3) * logsigma - 32. * sqrt1minus4z * (1. + 2. * z) * logx_2 -
           48. * sqrt1minus4z * (1. + 2. * z) * log1minus4z + 64. * sqrt1minus4z * (1. + 2. * z) * logz +
           16. * (1. - 4. * z2) * (3. * log2sigma + 2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2);
    cacheF1[indexF(cc, 1, 1, 2)] = -sqrt1minus4z * (127. - 199. * z - 75. * z2) / 9. + (2. - 259. * z + 662. * z2 - 76. * z3 - 200. * z4) * logsigma / (12. * z) -
           (17. - 26. * z) * sqrt1minus4z * logx_1 + 2./3. * (5. - 8. * z) * sqrt1minus4z * logx_2 - 4. * sqrt1minus4z * (3. - 2. * z) * log1minus4z -
           sqrt1minus4z * (2. - 255. * z + 316. * z2) * logz / (12. * z) +
           16./3. * (1. - 3. * z + 2. * z2) * (3. * log2sigma + 2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2);
    cacheF1[indexF(cc, 2, 1, 2)] = -2. * sqrt1minus4z * (68. + 49. * z - 150. * z2) / 9. + 2./3. * (1. - 35. * z + 4. * z2 + 76. * z3 - 100. * z4) * logsigma/z +
           (16. - 64. * z2) * log2sigma - 8. * sqrt1minus4z * (1. + 2. * z) * logx_1 - 32./3. * sqrt1minus4z * (1. + 2. * z) * logx_2 -
           16. * sqrt1minus4z * (1. + 2. * z) * log1minus4z - 2./3. * sqrt1minus4z * (1. - 33. * z - 76. * z2) * logz/z +
           16./3. * (1. - 4. * z2) * (2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2);
    cacheF1[indexF(cc, 1, 2, 2)] = -M_PI2/3. * (1. - 10. * z) - sqrt1minus4z * (115. + 632. * z + 96. * z2) / 18. - (7. + 13. * z - 194. * z2 + 304. * z3 - 64. * z4) * logsigma / (6. * z) -
           2. * sqrt1minus4z * (5. - 2. * z) * logx_1 + 4./3. * (2. - 5. * z) * sqrt1minus4z * logx_2 - 4. * (1. - 6. * z) * sqrt1minus4z * log1minus4z +
           (13. - 54. * z + 8. * z2) * logsigma * log1minus4z / 3. + sqrt1minus4z * (7. + 27. * z - 250. * z2) * logz / (6. * z) +
           (7. - 32. * z + 4. * z2) * (log2sigma - logsigma * logz) + 4./3. * (5. - 12. * z + 4. * z2) * Dilogsigma + 4./3. * (4. - 21. * z + 2. * z2) * Dilogsigma2;
    cacheF1[indexF(cc, 2, 2, 2)] = 8./3. * M_PI2 * (1. + 2. * z) - 8./9. * sqrt1minus4z * (19. + 53. * z + 24. * z2) + 4./3. * (1. + 7. * z + 10. * z2 - 68. * z3 + 32. * z4) * logsigma/z -
            8. * (1. + 2. * z) * (1. + 2. * z) * log2sigma - 16. * sqrt1minus4z * (1. + 2. * z) * logx_1 + 32./3. * sqrt1minus4z * (1. + 2. * z) * logx_2 +
            16. * sqrt1minus4z * (1. + 2. * z) * log1minus4z - 8./3. * (1. + 6. * z + 8. * z2) * logsigma * log1minus4z -
            4./3. * sqrt1minus4z * (1. + 9. * z + 26. * z2) * logz / z + 8. * (1. + 2. * z) * (1. + 2. * z) * logsigma * logz +
            32./3. * (1. - 4. * z2) * Dilogsigma - 32./3. * (1. + 3. * z + 2. * z2) * Dilogsigma2;

    //differs dgamma.c
    cacheF1[indexF(uu, 1, 1, 1)] = 109./6. + 10. * logx_2;
    cacheF1[indexF(uu, 2, 1, 1)] = -40./3. - 32. * logx_2;
    cacheF1[indexF(uu, 1, 1, 2)] = -127./9. + 4./12. - 17. * logx_1 + 10./3. * logx_2;
    cacheF1[indexF(uu, 2, 1, 2)] = -(136. + 12.)/9. - 8. * logx_1 - 32./3. * logx_2;
    cacheF1[indexF(uu, 1, 2, 2)] = -M_PI2/3. - (115. + 42.)/18. - 10. * logx_1 + 8./3. * logx_2;
    cacheF1[indexF(uu, 2, 2, 2)] = 8. * M_PI2/3. - 8./9. * (19. - 3.) - 16. * logx_1 + 32./3. * logx_2;
    for (int k = 1; k < 3; k++) {
        for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1)) {
            std::cout << "enum " << qq << "\n";
            cacheF1[indexF(qq, k, 2, 1)] = cacheF1[indexF(qq, k, 1, 2)];
        }        
    }
   
//    std::cout << "F1 " << cacheF1[indexF(cc, 1, 1, 1)] << "\n";
//    std::cout << "F1 " << cacheF1[indexF(cc, 2, 1, 1)] << "\n";
//    std::cout << "F1 " << cacheF1[indexF(cc, 1, 1, 2)] << "\n";
//    std::cout << "F1 " << cacheF1[indexF(cc, 2, 1, 2)] << "\n";
//    std::cout << "F1 " << cacheF1[indexF(cc, 1, 2, 2)] << "\n";
//    std::cout << "F1 " << cacheF1[indexF(cc, 2, 2, 2)] << "\n";
    
    return;
}

void AmpDB2::computeP() {
    cacheP[indexP(cu, 1, 2, 2)] = -1./27. + 2./9. * z - logx_1/9. - sqrt1minus4z * (1. + 2. * z) * (2. + 3. * logsigma + 6. * logx_1) / 54. + logz / 18.;
    cacheP[indexP(cu, 2, 2, 2)] = 8./27. + 16./9. * z + 8./9. * logx_1 + 4./27. * sqrt1minus4z * (1. + 2. * z) * (2. + 3. * logsigma + 6. * logx_1) - 4. * logz / 9.;
    cacheP[indexP(cc, 1, 2, 2)] = -2./27. * sqrt1minus4z * (1. + 8. * z + 12. * z2) - logsigma / 9. + 4./3. * z2 * logsigma + 16./9. * z3 * logsigma -
           sqrt1minus4z * (1. + 2. * z) * (2. * logx_1 - logz) / 9.;
    cacheP[indexP(cc, 2, 2, 2)] = 16./27. * sqrt1minus4z * (1. + 8. * z + 12. * z2) + 8./9. * logsigma - 32./3. * z2 * logsigma - 128./9. * z3 * logsigma +
            8./9. * sqrt1minus4z * (1. + 2. * z) * (2. * logx_1 - logz);

    //differs to dgamma.c
    cacheP[indexP(uu, 1, 2, 2)] = -2./27. - 2./9. * logx_1;
    cacheP[indexP(uu, 2, 2, 2)] = 16./27. + 16./9. * logx_1;

    cacheP[indexP(cc, 1, 1, 3)] = 3. * sqrt1minus4z * (1. - z);
    cacheP[indexP(cc, 2, 1, 3)] = 3. * sqrt1minus4z * (1. + 2. * z);
    cacheP[indexP(cc, 1, 2, 3)] = sqrt1minus4z * (1. - z);
    cacheP[indexP(cc, 2, 2, 3)] = sqrt1minus4z * (1. + 2. * z);
    cacheP[indexP(cc, 1, 1, 4)] = sqrt1minus4z * (1. - z);
    cacheP[indexP(cc, 2, 1, 4)] = sqrt1minus4z * (1. + 2. * z);
    cacheP[indexP(cc, 1, 2, 4)] = 0.5 * sqrt1minus4z * sqrt1minus4z * sqrt1minus4z;
    cacheP[indexP(cc, 2, 2, 4)] = -sqrt1minus4z * (1. + 2. * z);
    cacheP[indexP(cc, 1, 1, 5)] = 9. * z * sqrt1minus4z;
    cacheP[indexP(cc, 2, 1, 5)] = 0.;
    cacheP[indexP(cc, 1, 2, 5)] = 3. * z * sqrt1minus4z;
    cacheP[indexP(cc, 2, 2, 5)] = 0.;
    cacheP[indexP(cc, 1, 1, 6)] = 3. * z * sqrt1minus4z;
    cacheP[indexP(cc, 2, 1, 6)] = 0.;
    cacheP[indexP(cc, 1, 2, 6)] = 3. * z * sqrt1minus4z;
    cacheP[indexP(cc, 2, 2, 6)] = 0.;
    cacheP[indexP(cc, 1, 2, 8)] = -1./6. * sqrt1minus4z * (1. + 2. * z);
    cacheP[indexP(cc, 2, 2, 8)] = 4./3. * sqrt1minus4z * (1. + 2. * z);
    
    cacheP[indexP(uu, 1, 1, 3)] = 3.;
    cacheP[indexP(uu, 2, 1, 3)] = 3.;
    cacheP[indexP(uu, 1, 2, 3)] = 1.;
    cacheP[indexP(uu, 2, 2, 3)] = 1.;
    cacheP[indexP(uu, 1, 1, 4)] = 1.;
    cacheP[indexP(uu, 2, 1, 4)] = 1.;
    cacheP[indexP(uu, 1, 2, 4)] = 0.5;
    cacheP[indexP(uu, 2, 2, 4)] = -1.;
    cacheP[indexP(uu, 1, 1, 5)] = 0.;
    cacheP[indexP(uu, 2, 1, 5)] = 0.;
    cacheP[indexP(uu, 1, 2, 5)] = 0.;
    cacheP[indexP(uu, 2, 2, 5)] = 0.;
    cacheP[indexP(uu, 1, 1, 6)] = 0.;
    cacheP[indexP(uu, 2, 1, 6)] = 0.;
    cacheP[indexP(uu, 1, 2, 6)] = 0.;
    cacheP[indexP(uu, 2, 2, 6)] = 0.;
    cacheP[indexP(uu, 1, 2, 8)] = -1./6.;
    cacheP[indexP(uu, 2, 2, 8)] = 4./3.;
    return;
}

//equation 39:
void AmpDB2::computeD() {
    //qq = uu and cc
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2)) {
        for (int k = 0; k <= 1; k++) {
            gslpp::complex result = 0.;
            for (int i = 1; i <= 2; i++) {
                for (int j = 1; j <= 2; j++) {
                    result += C(i) * C(j) * F(qq, k, i, j);
                }
            }
            result += + as_4pi * C(2) * C(2) * P(qq, k, 2, 2)
                    + 2. * as_4pi * C(2) * C_8G * P(qq, k, 2, 8);
            for (int i = 1; i <= 2; i++) {
                for (int r = 3; r <= 6; r++) {
                    result += 2. * C(i) * C(r) * P(qq, k, i, r);
                }                    
            }
            cacheD[indexD(qq, k)] = result;
        }
    }
    //qq = cu
    for (int k = 0; k <= 1; k++) {
        gslpp::complex result = 0.;
        for (int i = 1; i <= 2; i++) {
            for (int j = 1; j <= 2; j++) {
                result += C(i) * C(j) * F(cu, k, i, j);
            }
        }
        result += + as_4pi * C(2) * C(2) * P(cu, k, 2, 2)
                + as_4pi * C(2) * C_8G * (P(cc, k, 2, 8) + P(uu, k, 2, 8));
        for (int i = 1; i <= 2; i++) {
            for (int r = 3; r <= 6; r++) {
                result += C(i) * C(r) * (P(cc, k, i, r) + P(uu, k, i, r));
            }
        }
        cacheD[indexD(cu, k)] = result;
    }
}

void AmpDB2::computeCKMelements() {
    VtbVtd = mySM.getCKM().getV_tb().conjugate() * mySM.getCKM().getV_td();
    VtbVts = mySM.getCKM().getV_tb().conjugate() * mySM.getCKM().getV_ts();
    VtbVtd2 = VtbVtd * VtbVtd;
    VtbVts2 = VtbVts * VtbVts;
    VcbVcd = mySM.getCKM().getV_cb().conjugate() * mySM.getCKM().getV_cd();
    VcbVcs = mySM.getCKM().getV_cb().conjugate() * mySM.getCKM().getV_cs();
    VcbVcd2 = VcbVcd * VcbVcd;
    VcbVcs2 = VcbVcs * VcbVcs;
//    std::cout << F0(cu, 1, 2, 1)<< "\n";
//    std::cout << F0(cu, 1, 1, 2)<< "\n";
//    std::cout << F1(cc, 2, 2, 1)<< "\n";
//    std::cout << F1(cc, 2, 1, 2)<< "\n";
    
//    std::cout << F0(cu, 2, 1, 2)<< "\n";
//    std::cout << F0(cu, 2, 2, 2)<< "\n";
//    std::cout << F0(cc, 1, 1, 1)<< "\n";
//    std::cout << F0(cc, 1, 1, 2)<< "\n";
//    std::cout << F0(cc, 1, 2, 2)<< "\n";
//    std::cout << F0(cc, 2, 1, 1)<< "\n";
//    std::cout << F0(cc, 2, 1, 2)<< "\n";
//    std::cout << F0(cc, 2, 2, 2)<< "\n";
//    std::cout << F0(uu, 1, 1, 1)<< "\n";
//    std::cout << F0(uu, 1, 1, 2)<< "\n";
//    std::cout << F0(uu, 1, 2, 2)<< "\n";
//    std::cout << F0(uu, 2, 1, 1)<< "\n";
//    std::cout << F0(uu, 2, 1, 2)<< "\n";
//    std::cout << F0(uu, 2, 2, 2)<< "\n";
    
    std::cout << "--- P ---" << "\n";
//    std::cout << P(cu, 1, 2, 2)<< "\n";
//    std::cout << P(cu, 2, 2, 2)<< "\n";
//    std::cout << P(cc, 1, 2, 2)<< "\n";
//    std::cout << P(cc, 2, 2, 2)<< "\n";
//    std::cout << P(uu, 1, 2, 2)<< "\n";
//    std::cout << P(uu, 2, 2, 2)<< "\n";
//    
//    std::cout << P(uu, 1, 1, 3)<< "\n";
//    std::cout << P(uu, 2, 1, 3)<< "\n";
//    std::cout << P(uu, 1, 2, 3)<< "\n";
//    std::cout << P(uu, 2, 2, 3)<< "\n";
//    std::cout << P(uu, 1, 1, 4)<< "\n";
//    std::cout << P(uu, 2, 1, 4)<< "\n";
//    std::cout << P(uu, 1, 2, 4)<< "\n";
//    std::cout << P(uu, 2, 2, 4)<< "\n";
//    std::cout << P(uu, 1, 1, 5)<< "\n";
//    std::cout << P(uu, 2, 1, 5)<< "\n";
//    std::cout << P(uu, 1, 2, 5)<< "\n";
//    std::cout << P(uu, 2, 2, 5)<< "\n";
//    std::cout << P(uu, 1, 1, 6)<< "\n";
//    std::cout << P(uu, 2, 1, 6)<< "\n";
//    std::cout << P(uu, 1, 2, 6)<< "\n";
//    std::cout << P(uu, 2, 2, 6)<< "\n";
//    std::cout << P(uu, 1, 2, 8)<< "\n";
//    std::cout << P(uu, 2, 2, 8)<< "\n";    
//    
//    std::cout << P(cc, 1, 1, 3)<< "\n";
//    std::cout << P(cc, 2, 1, 3)<< "\n";
//    std::cout << P(cc, 1, 2, 3)<< "\n";
//    std::cout << P(cc, 2, 2, 3)<< "\n";
//    std::cout << P(cc, 1, 1, 4)<< "\n";
//    std::cout << P(cc, 2, 1, 4)<< "\n";
//    std::cout << P(cc, 1, 2, 4)<< "\n";
//    std::cout << P(cc, 2, 2, 4)<< "\n";
//    std::cout << P(cc, 1, 1, 5)<< "\n";
//    std::cout << P(cc, 2, 1, 5)<< "\n";
//    std::cout << P(cc, 1, 2, 5)<< "\n";
//    std::cout << P(cc, 2, 2, 5)<< "\n";
//    std::cout << P(cc, 1, 1, 6)<< "\n";
//    std::cout << P(cc, 2, 1, 6)<< "\n";
//    std::cout << P(cc, 1, 2, 6)<< "\n";
//    std::cout << P(cc, 2, 2, 6)<< "\n";
//    std::cout << P(cc, 1, 2, 8)<< "\n";
//    std::cout << P(cc, 2, 2, 8)<< "\n";
    
    //update others
    Gf2 = mySM.getGF() * mySM.getGF();
    z2 = z * z;
    z3 = z2 * z;
    z4 = z3 * z;
    logz = log(z);
    log1minusz = log(1 - z);
    log1minus4z = log(1. - 4. * z);
    oneminusz2 = (1 - z) * (1 - z);
    sqrt1minus4z = sqrt(1 - 4 * z);
    sigma = (1 - sqrt1minus4z)/(1 + sqrt1minus4z);
    logsigma = log(sigma);
    log2sigma = logsigma * logsigma;
    logx_1 = log(x_1);
    logx_2 = log(x_2);
    Dilogz = gslpp_special_functions::dilog(z);
    Dilogsigma = gslpp_special_functions::dilog(sigma);
    Dilogsigma2 = gslpp_special_functions::dilog(sigma * sigma);
    
    mu_1 = Mb;
    mu_2 = Mb;
    x_1 = mu_1/Mb;
    x_2 = mu_2/Mb;
    
    //check mu
    as_4pi = mySM.getAlsM() / (4. * M_PI);
    return;
}

void AmpDB2::computeWilsonCoeffs(QCD::lepton lep){
    gslpp::vector<gslpp::complex> ** WilsonCoeffs = mySM.getFlavour().ComputeCoeffBMll_Buras(Mb, lep);
    for (int i = 0; i < 6; i++) {
        cacheC[i] = (*(WilsonCoeffs[LO]))(i) + (*(WilsonCoeffs[NLO]))(i);
    }
    C_8G = 0.;

    std::cout.precision(4);
    std::cout << "C_1" << C(1) << "\n";
    std::cout << "C_2" << C(2) << "\n";
    std::cout << "C_3" << C(3) << "\n";
    std::cout << "C_4" << C(4) << "\n";
    std::cout << "C_5" << C(5) << "\n";
    std::cout << "C_6" << C(6) << "\n";

    K_1 = 3. * C(1) * C(1) + 2. * C(1) * C(2);
    K_2 = C(2) * C(2);
    K_1prime = 2. * (3. * C(1) * C(3) + C(1) * C(4) + C(2) * C(3));
    K_2prime = 2. * C(2) * C(4);
    K_3prime = 2 * (3 * C(1) * C(5) + C(1) * C(6) + C(2) * C(5) + C(2) * C(6));
    K12 = K_1 + K_2;
}

gslpp::complex AmpDB2::Gamma12_Bd(orders order) {
    if (order != FULLNLO) throw std::runtime_error("AmpDB2::Gamma12overM12_Bd(): order not implemented");
    else {
        //hep-ph/0308029v2
        
        //double mu_2 = 0;        

        //computeCKMelements();
 
        //equation 16
        gslpp::complex Gamma12_Bd = -Gf2 * Mb2 / (24 * M_PI * MB) *
                (c(d)(0) * me(0) + c(d)(1) * me(1) + delta_1overm(d));
        return Gamma12_Bd;
    }
}

gslpp::complex AmpDB2::Gamma12_Bs(orders order) {
    if (order != FULLNLO) throw std::runtime_error("AmpDB2::Gamma12overM12_Bs(): order not implemented");
    else {
        //hep-ph/0308029v2
        
        //double mu_2 = 0;        

        //computeCKMelements();
 
        //equation 16
        gslpp::complex Gamma12_Bs = -Gf2 * Mb2 / (24 * M_PI * MB_s) *
                (c(s)(0) * me(0) + c(s)(1) * me(1) + delta_1overm(s));
        return Gamma12_Bs;
    }
}

gslpp::complex AmpDB2::PBd() {
    double mbpole = mySM.Mbar2Mp(mySM.getQuarks(QCD::BOTTOM).getMass());
    double Mw = mySM.Mw();
    double kappa = -2. * M_PI * mbpole * mbpole /
            (3. * Mw * Mw * mySM.getFlavour().getHDF2().getUDF2().etabS0(mySM.getBBd().getMu()));

    double n[13] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    n[0] = 0.1797;
    n[1] = 0.1391;
    n[5] = 1.0116;
    n[6] = 0.0455;
    n[8] = -0.0714;
    n[10] = -0.3331;

    double B1 = mySM.getBBd().getBpars()(0);
    double B2 = mySM.getBBd().getBpars()(1);

    gslpp::complex PBd = -2. * kappa / mySM.getCBd() *
            (gslpp::complex(1, 2. * mySM.getPhiBd(), true) * (n[0] + (n[5] * B2 + n[10]) / B1)
            - gslpp::complex(1. / mySM.getCKM().computeRt(), mySM.getCKM().computeBeta() + 2. * mySM.getPhiBd(), true)
            * (n[1] + (n[6] * B2 + n[11]) / B1)
            + gslpp::complex(1. / mySM.getCKM().computeRt() / mySM.getCKM().computeRt(), 2. * (mySM.getCKM().computeBeta() + mySM.getPhiBd()), true)
            * (n[2] + (n[7] * B2 + n[12]) / B1));

    return PBd;
}

gslpp::complex AmpDB2::PBs() {
    double mbpole = mySM.Mbar2Mp(mySM.getQuarks(QCD::BOTTOM).getMass());
    double Mw = mySM.Mw();
    double kappa = -2. * M_PI * mbpole * mbpole /
            (3. * Mw * Mw * mySM.getFlavour().getHDF2().getUDF2().etabS0(mySM.getBBs().getMu()));

    double n[13] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    n[0] = 0.1797;
    n[1] = 0.1391;
    n[5] = 1.0116;
    n[6] = 0.0455;
    n[8] = -0.0714;
    n[10] = -0.3331;

    double B1 = mySM.getBBs().getBpars()(0);
    double B2 = mySM.getBBs().getBpars()(1);

    gslpp::complex PBs = -2. * kappa / mySM.getCBs() *
            (gslpp::complex(1, 2. * mySM.getPhiBs(), true) * (n[0] + (n[5] * B2 + n[10]) / B1)
            - gslpp::complex(1. / mySM.getCKM().computeRts(), -mySM.getCKM().computeBetas() + 2. * mySM.getPhiBs(), true)
            * (n[1] + (n[6] * B2 + n[11]) / B1)
            + gslpp::complex(1. / mySM.getCKM().computeRts() / mySM.getCKM().computeRts(), 2. * (-mySM.getCKM().computeBetas() + mySM.getPhiBs()), true)
            * (n[2] + (n[7] * B2 + n[12]) / B1));

    return PBs;
}

gslpp::complex AmpDB2::C(int i){
    return cacheC[i - 1];
}

double AmpDB2::Asl(orders order, QCD::lepton lep) {
    rhob = mySM.getCKM().getRhoBar();
    etab = mySM.getCKM().getEtaBar();
    Md = mySM.getQuarks(QCD::DOWN).getMass();
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass();
    Mc = mySM.getQuarks(QCD::CHARM).getMass();
    Ms = mySM.getQuarks(QCD::STRANGE).getMass();
    Mt = mySM.getQuarks(QCD::TOP).getMass();
    MW = mySM.Mw();
    MB = mySM.getMesons(QCD::B_D).getMass();
    MB_s = mySM.getMesons(QCD::B_S).getMass();

    //masses?
    std::cout << "Ms" << mySM.Mrun(mySM.getBBs().getMu(),
            mySM.getQuarks(QCD::STRANGE).getMass_scale(),
            mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO) << "\n";
    std::cout << "Ms" << mySM.getQuarks(QCD::STRANGE).getMass() << "\n";
    std::cout << "Mb_pole" << mySM.Mbar2Mp(Mb) << "\n";
    std::cout << "Mb" << Mb << "\n";

    std::cout << "test: " << mySM.getFlavour().getHDF2().getUDF2().etabS0(mySM.getBBd().getMu()) << "\n";

    Mb2 = Mb * Mb;
    Mt2 = Mt * Mt;
    MB2 = MB * MB;
    MW2 = MW * MW;

    z = Mc * Mc / Mb2;

    computeWilsonCoeffs(lep);
    computeCKMelements();
    computeF0();
    computeF1();
    computeP();
    computeD();
    compute_matrixelements(d);
    compute_deltas_1overm(d);
    std::cout << "Gamma12_Bd: "  << Gamma12_Bd(FULLNLO) << "\n";
    std::cout << "Gamma12oderM12_Bd: "  << Gamma12_Bd(FULLNLO)/M12_Bd(FULLNLO) << "\n";

    //mySM.getBBd().setFlagCsi(true);
    std::cout << "getBpars() " << mySM.getBBd().getBpars() << "\n";

    B1 = mySM.getBBd().getBpars()(0);
    B2 = mySM.getBBd().getBpars()(1);
    B_sprime = MB2 / ((Mb + Md)*(Mb + Md)) * B2;

    eta_B = 0.55;
    eta_Bb = eta_B * pow(mySM.Als(Mb), -6. / 23.) * (1. + mySM.Als(Mb) / (4. * M_PI) * 5165. / 3174.);
    S_0 = StandardModelMatching(mySM).S0(Mt2 / MW2, Mt2 / MW2);
    //std::cout << "eta_bar*S0" << eta_Bb*S_0 << "\n";

    kappa = 4. * M_PI * Mb2 / MW2 * K12 / (eta_Bb * S_0) * z;

    //StandardModel CP asymmetry in semileptonic B decay from hep-ph/0202010v2
    return (-kappa * etab / ((1. - rhob) * (1. - rhob) + etab * etab)
            *
            (
            1. + z * (5. / 4. * B_sprime / B1 * (K_2 - K_1) / K12 - K_2 / K12)
            + z * z * (K_2 - K_1) / K12 * (1. / 3. - 5. / 6. * B_sprime / B1)
            + 2. * z * z * (1. - etab) / ((1. - rhob) * (1. - rhob) + etab * etab) * (K_2 - K_1) / K12 * (5. / 2. * B_sprime / B1 - 1.)
            + z * (7. * K_1 + 3 * K_2) / (2. * K12 * B1) * (MB2 - Mb2) / Mb2
            + (K_1prime + K_2prime - K_3prime) / K12
            )
            ).real();
}


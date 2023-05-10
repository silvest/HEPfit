/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AmpDB2.h"
#include "EvolDF2.h"
#include "HeffDF2.h"
#include <chrono>

AmpDB2::AmpDB2(const StandardModel& SM_i)
: mySM(SM_i), meMStoRI(5, 0.), coeffsMStoRI(3, 0.)
{
    mySM.initializeBParameter("BBs");
    mySM.initializeBParameter("BBd");
    
    double log2 = log(2);
    double meMStoRI0[5] = {-3. - 5./3.+8.*log2, 0., 0., 0., 0.},
    meMStoRI1[5] = {0., -13./3. + 67./9.+44./9.*log2, -1./3. - 1./9.+28./9.*log2, 0., 0.},
    meMStoRI2[5] = {0., -29./6 - 28./9.+28./9.*log2, 7./6. - 68./9.+44./9.*log2, 0., 0.},
    meMStoRI3[5] = {0., 0., 0., -5./3. + 13.-2./3.*log2, -3. + 1.+2.*log2},
    meMStoRI4[5] = {0., 0., 0., -7./2. + 11./2.+2.*log2, -1./6. - 1./2.-2./3.*log2};
    meMStoRI.assign(0, meMStoRI0);
    meMStoRI.assign(1, meMStoRI1);
    meMStoRI.assign(2, meMStoRI2);
    meMStoRI.assign(3, meMStoRI3);
    meMStoRI.assign(4, meMStoRI4);
    for (int i=0; i<=2; i++ ) {
        for (int j=0; j<=2; j++ ) {
            coeffsMStoRI.assign(i, j, meMStoRI(i,j));
        }
    }
    //std::cout << meMStoRI << "\n";
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
    //better?
//    double Mb = mySM.Mrun(mySM.getBBs().getMu(),
//            mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
//            mySM.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
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

gslpp::complex AmpDB2::M21_Bd(orders order) {
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

gslpp::complex AmpDB2::M21_Bs(orders order) {
    if (mySM.getFlavour().getHDF2().getCoeffBs().getOrder() < order % 3)
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed");

    //Wilson coefficients in same mass scale and scheme as B parameters
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffBs(
            mySM.getBBs().getMu(),
            mySM.getBBs().getScheme());

    gslpp::vector<double> me(mySM.getBBs().getBpars());
    double MBs = mySM.getMesons(QCD::B_S).getMass();
    double Mb = mySM.getQuarks(QCD::BOTTOM).getMass();
    //better?
//    double Mb = mySM.Mrun(mySM.getBBs().getMu(),
//            mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
//            mySM.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    double Ms = mySM.Mrun(mySM.getBBs().getMu(),
            mySM.getQuarks(QCD::STRANGE).getMass_scale(),
            mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double KBs = MBs / (Mb + Ms) * MBs / (Mb + Ms);
    double Fbs = mySM.getMesons(QCD::B_S).getDecayconst();
    
    //matrix elements as in hep-ph/9604387v2 with KBs independent terms in 4 and 5
    //differ by factor of 2 to arXiv:1602.03560v2 etc.
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
    gslpp::vector< complex > c(3, 0.);
    switch (q) {
        case d:
            for (int i = 0; i <= 1; i++) {
                c.assign(i,
                        VtbVtd2 * D(uu, i) + 2. * VcbVcd * VtbVtd * (D(uu, i) - D(cu, i))
                        + VcbVcd2 * (D(cc, i) + D(uu, i) - 2. * D(cu, i))
                        );
            }
            //std::cout << "D " << D(uu, 1) << " " << D(uu,1)-D(cu,1) << " " << D(cc,1)+D(uu,1)-2.*D(cu,1) << "\n";
            //std::cout << "D " << D(uu, 2) << " " << D(uu,2)-D(cu,2) << " " << D(cc,2)+D(uu,2)-2.*D(cu,2) << "\n";            
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
    //change to RI
    //result += as_4pi * coeffsMStoRI.transpose() * result;
    return c;
}

void AmpDB2::compute_matrixelements(quark q){
    double Mq;
    double MBq;
    double MBq2;
    double KBq;
    double FBq2;
    switch (q) {
        case d:
            me = mySM.getBBd().getBpars();
            me1tilde = 1.2; //mySM.getOptionalParameter("BBd2tilde");
            Mq = Md;
            MBq2 = MB2;
            FBq2 = mySM.getMesons(QCD::B_D).getDecayconst() * mySM.getMesons(QCD::B_D).getDecayconst();
            break;
        case s:
            me = mySM.getBBs().getBpars();
            me1tilde = 1.31; // mySM.getOptionalParameter("BBs2tilde");
            Mq = Ms;
            MBq2 = MB_s * MB_s;
            FBq2 = mySM.getMesons(QCD::B_S).getDecayconst() * mySM.getMesons(QCD::B_S).getDecayconst();            
            break;
        default:
            throw std::runtime_error("AmpDB2::compute_matrixelements(quark q): invalid quark index: ");
    }
    KBq = MBq2 / ((Mb + Mq) * (Mb + Mq));
    
    //equation (26)
    me(0) *=  8. / 3. * MBq2 * FBq2;
    me(1) *= -5. / 3. * KBq * MBq2 * FBq2;
    me(2) *=  1. / 3. * KBq * MBq2 * FBq2;
    me(3) *=       2. * KBq * MBq2 * FBq2;
    me(4) *=  2. / 3. * KBq * MBq2 * FBq2;
    
    //switch to RI
    //me = me - as_4pi * meMStoRI;
    
    //equation (7.5)
    //should be just me(2)
    me1tilde *= 1. / 3. * MBq2 * FBq2;
    //std::cout << "me " << me1tilde << " " << me(2) << "\n";
    //me1tilde = me(2);
    
    me_R(0) = Mq/Mb * me(3);
    //new lattice results?
    me_R(1) = -2./3. * FBq2 * MBq2 * (MBq2 / Mb2 -1.);
    me_R(2) =  7./6. * FBq2 * MBq2 * (MBq2 / Mb2 -1.);
    me_R(3) = me(2) + 0.5 * me(0) + me(1) - 2. * Mq/Mb * me(4) + me_R(1);
    
     
    
    //as for M21 (different from equation (26))
    switch (q) {
        case d:
            me = mySM.getBBd().getBpars();
            MBq = MB;
            Mq = mySM.Mrun(mySM.getBBd().getMu(),
            mySM.getQuarks(QCD::DOWN).getMass_scale(),
            mySM.getQuarks(QCD::DOWN).getMass(), FULLNNLO);
            break;
        case s:
            me = mySM.getBBs().getBpars();
            MBq = MB_s;
            Mq = mySM.Mrun(mySM.getBBs().getMu(),
                mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO); 
            break;
    }
    KBq = MBq2 / ((Mb + Mq) * (Mb + Mq));
    me(0) *=  1. /  3. * MBq * FBq2;
    me(1) *= -5. / 24. * KBq * MBq * FBq2;
    me(2) *=  1. / 24. * KBq * MBq * FBq2;
    me(3) *=  1. /  4. * KBq * MBq * FBq2;
    me(4) *=  1. / 12. * KBq * MBq * FBq2;    

    std::cout << "me" << me << "\n";
    std::cout << "me_R " << me_R << "\n";
    return;
}

//equation 18
gslpp::complex AmpDB2::delta_1overm(quark q) {
    //return 0.;
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
            cacheF0[indexF(qq, k, 2, 1)] = cacheF0[indexF(qq, k, 1, 2)];
        }
    }
    return;
}

//equations 43-48: F1 = B
//with resummation of Mc via hep-ph/0307344 eq.(23) for that F0 has be computed before
void AmpDB2::computeF1() {
    cacheF1[indexF(cu, 1, 1, 1)] = 109./6. - 37. * z + 1.5 * z2 + 52./3. * z3 + 2. * oneminusz2 * (5. + z) * logx_2 - 4. * oneminusz2 * (5. + 7. * z) * log1minusz -
            2. * z * (10. + 14. * z - 15. * z2) * logz + 8. * (2. - 3. * z + z3) * log1minusz * logz + 16. * (2. - 3. * z + z3) * Dilogz
            + flag_resumz * (32./3. * F0(cu, 1, 1, 1) - 8. * z * logz * 1.5 * (-3. + 3. * z2));
    cacheF1[indexF(cu, 2, 1, 1)] = -4./3. * (10. - 33. * z + 54. * z2 - 31. * z3) - 8. * oneminusz2 * (4. + 14. * z - 3. * z2) * log1minusz +
            8. * z * (2. - 23. * z + 21. * z2 - 3. * z3) * logz -
            16. * oneminusz2 * (1. + 2. * z) * (2. * logx_2 - log1minusz * logz - 2. * Dilogz)
            + flag_resumz * (32./3. * F0(cu, 2, 1, 1) - 8. * z * logz * 3. * 6. * (z - 1.) * z);
    cacheF1[indexF(cu, 1, 1, 2)] = 2. * (
           (-502. + 912. * z - 387. * z2 - 23. * z3) / 36. - oneminusz2 * (17. + 4. * z) * logx_1 + 2./3. * oneminusz2 * (5. + z) * logx_2 -
           oneminusz2 / (12. * z) * (2. + 33. * z + 94. * z2) * log1minusz - z/12. * (80. + 69. * z - 126. * z2) * logz +
           8./3. * (2. - 3. * z + z3) * (log1minusz * logz + 2. * Dilogz)
           ) + flag_resumz * (32./3. * F0(cu, 1, 1, 2) - 8. * z * logz * 0.5 * (-3. + 3.* z2));
    cacheF1[indexF(cu, 2, 1, 2)] = 2. * (
           (-130. + 93. * z + 144. * z2 - 107. * z3) / 9. - 2./3. * oneminusz2 / z * (1. + 15. * z + 47. * z2 - 12. * z3) * log1minusz +
           2./3. * z * (8. - 93. * z + 87. * z2 - 12. * z3) * logz -
           8./3. * oneminusz2 * (1. + 2. * z) * (3. * logx_1 + 4. * logx_2 - 2. * log1minusz * logz - 4. * Dilogz)
           ) + flag_resumz * (32./3. * F0(cu, 2, 1, 2) - 8. * z * logz * 6. * (z - 1.) * z);
    cacheF1[indexF(cu, 1, 2, 2)] = -M_PI2/3. * (1. - 5. * z + 4. * z2) + (-136. - 159. * z + 738. * z2 - 443. * z3) / 18. - 2 * oneminusz2 * (5. + 4. * z) * logx_1 +
           2./3. * oneminusz2 * (4. - z) * logx_2 + oneminusz2 / (6. * z) * (7. + 32. * z2 + 3. * z3) * log1minusz -
           z/6. * (62. + 39. * z - 30. * z2 + 3. * z3) * logz + (5. - 3. * z - 18. * z2 + 16. * z3) / 3. * (log1minusz * logz + 2. * Dilogz)
            + flag_resumz * (32./3. * F0(cu, 1, 2, 2) - 8. * z * logz * -1.5 * oneminusz2);
    cacheF1[indexF(cu, 2, 2, 2)] = 8./3. * M_PI2 * (1. + z - 2. * z2) - 28./9. * (5. + 3. * z - 27. * z2 + 19. * z3) - 16. * oneminusz2 * (1. + 2. * z) * logx_1 +
           32./3. * oneminusz2 * (1. + 2. * z) * logx_2 - 4./3. * oneminusz2 / z * (1. - 12. * z - 16. * z2 - 3. * z3) * log1minusz +
           4./3. * z * (2. -  3. * z + 18. * z2 - 3. * z3) * logz + 8./3. * (1. - 3. * z - 6. * z2 + 8. * z3) * (log1minusz * logz + 2. * Dilogz)
            + flag_resumz * (32./3. * F0(cu, 2, 2, 2) - 8. * z * logz * -6. * (z - 1.) * z);

    cacheF1[indexF(cc, 1, 1, 1)] = sqrt1minus4z * (109. - 226. * z + 168. * z2) / 6. - (52. - 104. * z - 16. * z2 + 56. * z3) * logsigma +
           2. * (5. - 8. * z) * sqrt1minus4z * logx_2 - 12. * sqrt1minus4z * (3. - 2. * z) * log1minus4z + 4. * (13. - 10. * z) * sqrt1minus4z * logz +
           16. * (1. - 3. * z + 2. * z2) * (3. * log2sigma + 2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2)
            + flag_resumz * (32./3. * F0(cc, 1, 1, 1) - 8. * z * logz * 3. * (6. * z - 3.) / sqrt1minus4z);
    cacheF1[indexF(cc, 2, 1, 1)] = -8./3. * sqrt1minus4z * (5. - 23. * z - 42. * z2) - 16. * (4. - 2. * z - 7. * z2 + 14. * z3) * logsigma - 32. * sqrt1minus4z * (1. + 2. * z) * logx_2 -
           48. * sqrt1minus4z * (1. + 2. * z) * log1minus4z + 64. * sqrt1minus4z * (1. + 2. * z) * logz +
           16. * (1. - 4. * z2) * (3. * log2sigma + 2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2)
            + flag_resumz * (32./3. * F0(cc, 2, 1, 1) - 8. * z * logz * 3. * -12. * z / sqrt1minus4z);
    cacheF1[indexF(cc, 1, 1, 2)] = -sqrt1minus4z * (127. - 199. * z - 75. * z2) / 9. + (2. - 259. * z + 662. * z2 - 76. * z3 - 200. * z4) * logsigma / (12. * z) -
           (17. - 26. * z) * sqrt1minus4z * logx_1 + 2./3. * (5. - 8. * z) * sqrt1minus4z * logx_2 - 4. * sqrt1minus4z * (3. - 2. * z) * log1minus4z -
           sqrt1minus4z * (2. - 255. * z + 316. * z2) * logz / (12. * z) +
           16./3. * (1. - 3. * z + 2. * z2) * (3. * log2sigma + 2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2)
            + flag_resumz * (32./3. * F0(cc, 1, 1, 2) - 8. * z * logz * (6. * z - 3.) / sqrt1minus4z);
    cacheF1[indexF(cc, 2, 1, 2)] = -2. * sqrt1minus4z * (68. + 49. * z - 150. * z2) / 9. + 2./3. * (1. - 35. * z + 4. * z2 + 76. * z3 - 100. * z4) * logsigma/z +
           (16. - 64. * z2) * log2sigma - 8. * sqrt1minus4z * (1. + 2. * z) * logx_1 - 32./3. * sqrt1minus4z * (1. + 2. * z) * logx_2 -
           16. * sqrt1minus4z * (1. + 2. * z) * log1minus4z - 2./3. * sqrt1minus4z * (1. - 33. * z - 76. * z2) * logz/z +
           16./3. * (1. - 4. * z2) * (2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2)
            + flag_resumz * (32./3. * F0(cc, 2, 1, 2) - 8. * z * logz * -12. * z / sqrt1minus4z);
    cacheF1[indexF(cc, 1, 2, 2)] = -M_PI2/3. * (1. - 10. * z) - sqrt1minus4z * (115. + 632. * z + 96. * z2) / 18. - (7. + 13. * z - 194. * z2 + 304. * z3 - 64. * z4) * logsigma / (6. * z) -
           2. * sqrt1minus4z * (5. - 2. * z) * logx_1 + 4./3. * (2. - 5. * z) * sqrt1minus4z * logx_2 - 4. * (1. - 6. * z) * sqrt1minus4z * log1minus4z +
           (13. - 54. * z + 8. * z2) * logsigma * log1minus4z / 3. + sqrt1minus4z * (7. + 27. * z - 250. * z2) * logz / (6. * z) +
           (7. - 32. * z + 4. * z2) * (log2sigma - logsigma * logz) + 4./3. * (5. - 12. * z + 4. * z2) * Dilogsigma + 4./3. * (4. - 21. * z + 2. * z2) * Dilogsigma2
            + flag_resumz * (32./3. * F0(cc, 1, 2, 2) - 8. * z * logz * 0.5 * -6. * sqrt1minus4z);
    cacheF1[indexF(cc, 2, 2, 2)] = 8./3. * M_PI2 * (1. + 2. * z) - 8./9. * sqrt1minus4z * (19. + 53. * z + 24. * z2) + 4./3. * (1. + 7. * z + 10. * z2 - 68. * z3 + 32. * z4) * logsigma/z -
            8. * (1. + 2. * z) * (1. + 2. * z) * log2sigma - 16. * sqrt1minus4z * (1. + 2. * z) * logx_1 + 32./3. * sqrt1minus4z * (1. + 2. * z) * logx_2 +
            16. * sqrt1minus4z * (1. + 2. * z) * log1minus4z - 8./3. * (1. + 6. * z + 8. * z2) * logsigma * log1minus4z -
            4./3. * sqrt1minus4z * (1. + 9. * z + 26. * z2) * logz / z + 8. * (1. + 2. * z) * (1. + 2. * z) * logsigma * logz +
            32./3. * (1. - 4. * z2) * Dilogsigma - 32./3. * (1. + 3. * z + 2. * z2) * Dilogsigma2
            + flag_resumz * (32./3. * F0(cc, 2, 2, 2) - 8. * z * logz * 12. * z / sqrt1minus4z);

    cacheF1[indexF(uu, 1, 1, 1)] = 109./6. + 10. * logx_2
            + flag_resumz * (32./3. * F0(uu, 1, 1, 1));
    cacheF1[indexF(uu, 2, 1, 1)] = -40./3. - 32. * logx_2
            + flag_resumz * (32./3. * F0(uu, 2, 1, 1));
    cacheF1[indexF(uu, 1, 1, 2)] = -127./9. + 4./12. - 17. * logx_1 + 10./3. * logx_2
            + flag_resumz * (32./3. * F0(uu, 1, 1, 2));            
    cacheF1[indexF(uu, 2, 1, 2)] = -(136. + 12.)/9. - 8. * logx_1 - 32./3. * logx_2
            + flag_resumz * (32./3. * F0(uu, 2, 1, 2));            
    cacheF1[indexF(uu, 1, 2, 2)] = -M_PI2/3. - (115. + 42.)/18. - 10. * logx_1 + 8./3. * logx_2
            + flag_resumz * (32./3. * F0(uu, 1, 2, 2));            
    cacheF1[indexF(uu, 2, 2, 2)] = 8. * M_PI2/3. - 8./9. * (19. - 3.) - 16. * logx_1 + 32./3. * logx_2
            + flag_resumz * (32./3. * F0(uu, 2, 2, 2));
           
    for (int k = 1; k < 3; k++) {
        for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1)) {
            cacheF1[indexF(qq, k, 2, 1)] = cacheF1[indexF(qq, k, 1, 2)];
        }
    }
    return;
}

void AmpDB2::computeP() {
    cacheP[indexP(cu, 1, 2, 2)] = -1./27. + 2./9. * z - logx_1/9. - sqrt1minus4z * (1. + 2. * z) * (2. + 3. * logsigma + 6. * logx_1) / 54. + logz / 18.;
    cacheP[indexP(cu, 2, 2, 2)] = 8./27. + 16./9. * z + 8./9. * logx_1 + 4./27. * sqrt1minus4z * (1. + 2. * z) * (2. + 3. * logsigma + 6. * logx_1) - 4. * logz / 9.;
    cacheP[indexP(cc, 1, 2, 2)] = -2./27. * sqrt1minus4z * (1. + 8. * z + 12. * z2) - logsigma / 9. + 4./3. * z2 * logsigma + 16./9. * z3 * logsigma -
           sqrt1minus4z * (1. + 2. * z) * (2. * logx_1 - logz) / 9.;
    cacheP[indexP(cc, 2, 2, 2)] = 16./27. * sqrt1minus4z * (1. + 8. * z + 12. * z2) + 8./9. * logsigma - 32./3. * z2 * logsigma - 128./9. * z3 * logsigma +
            8./9. * sqrt1minus4z * (1. + 2. * z) * (2. * logx_1 - logz);
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
    return;
}

void AmpDB2::computeCKMandMasses(orders order) {
    if (order != NLO and order != NNLO)
        throw(std::runtime_error("computeCKMandMasses() order not present"));
    
    VtbVtd = mySM.getCKM().getV_tb().conjugate() * mySM.getCKM().getV_td();
    VtbVts = mySM.getCKM().getV_tb().conjugate() * mySM.getCKM().getV_ts();
    VtbVtd2 = VtbVtd * VtbVtd;
    VtbVts2 = VtbVts * VtbVts;
    VcbVcd = mySM.getCKM().getV_cb().conjugate() * mySM.getCKM().getV_cd();
    VcbVcs = mySM.getCKM().getV_cb().conjugate() * mySM.getCKM().getV_cs();
    VcbVcd2 = VcbVcd * VcbVcd;
    VcbVcs2 = VcbVcs * VcbVcs;
    
    Gf2 = mySM.getGF() * mySM.getGF();
    Md = mySM.getQuarks(QCD::DOWN).getMass();
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass();
    Ms = mySM.getQuarks(QCD::STRANGE).getMass();
    MB = mySM.getMesons(QCD::B_D).getMass();
    MB_s = mySM.getMesons(QCD::B_S).getMass();
    Mb2 = Mb * Mb;
    MB2 = MB * MB;
    
    mu_1 = mySM.getMub();
    mu_2 = Mb;
    as_4pi = mySM.Alstilde5(mySM.getMub());
    
    if (order == NLO){
        Mc = mySM.Mrun(mySM.getBBd().getMu(),
            mySM.getQuarks(QCD::CHARM).getMass_scale(),
            mySM.getQuarks(QCD::CHARM).getMass(), FULLNNLO);
    
//        flag_resumz = false;
//        Mc = mySM.getQuarks(QCD::CHARM).getMass();
        
        x_1 = mu_1/Mb;
        x_2 = mu_2/Mb;
        logx_1 = log(x_1);
        logx_2 = log(x_2);
    }
    if (order == NNLO){
        Mc = mySM.Mrun(mySM.getBBd().getMu(),
            mySM.getQuarks(QCD::CHARM).getMass_scale(),
            mySM.getQuarks(QCD::CHARM).getMass(), FULLNNLO);
    }
    
    z = Mc * Mc / Mb2;
    z2 = z * z;
    logz = log(z);
    oneminusz2 = (1. - z) * (1. - z);
    sqrt1minus4z = sqrt(1. - 4. * z);
    
    if (order == NLO){
        z3 = z2 * z;
        z4 = z3 * z;

        log1minusz = log(1. - z);
        log1minus4z = log(1. - 4. * z);

        sigma = (1. - sqrt1minus4z)/(1. + sqrt1minus4z);
        logsigma = log(sigma);
        log2sigma = logsigma * logsigma;

        Dilogz = gslpp_special_functions::dilog(z);
        Dilogsigma = gslpp_special_functions::dilog(sigma);
        Dilogsigma2 = gslpp_special_functions::dilog(sigma * sigma);
    }
    return;
}

void AmpDB2::computeWilsonCoeffs(){
    double currentInput_computeWilsonCoeffs = mu_1;
    if (lastInput_computeWilsonCoeffs == currentInput_computeWilsonCoeffs) return;
   gslpp::vector<gslpp::complex> ** WilsonCoeffs = mySM.getFlavour().ComputeCoeffBMll_Buras(mu_1, QCD::NOLEPTON);
    for (int i = 0; i < 6; i++) {
        cacheC[i] = (*(WilsonCoeffs[LO]))(i) + (*(WilsonCoeffs[NLO]))(i);
    }
    C_8G = (*(WilsonCoeffs[LO]))(7) + (*(WilsonCoeffs[NLO]))(7);
    
//    for(int i=0; i<=7; i++){
//        if(i==6) i++;
//        std::cout << "C_" << i << " "
//                << (*(WilsonCoeffs[LO]))(i).gslpp::complex::real() << " "
//                << (*(WilsonCoeffs[NLO]))(i).gslpp::complex::real() << "\n";       
//    }
//    std::cout << "--------\n" ;

    K_1 = 3. * C(1) * C(1) + 2. * C(1) * C(2);
    K_2 = C(2) * C(2);
    lastInput_computeWilsonCoeffs = mu_1;
}

void AmpDB2::computeWilsonCoeffsDB1bsg(){
    double currentInput_computeWilsonCoeffsDB1bsg = mu_1;
    if (lastInput_computeWilsonCoeffsDB1bsg == currentInput_computeWilsonCoeffsDB1bsg) return;
 
    gslpp::vector<gslpp::complex> ** WilsonCoeffsDB1bsg = mySM.getFlavour().ComputeCoeffsgamma_Buras(mu_1);
    for (int i = 0; i < 6; i++) {
        cacheC[i] = (*(WilsonCoeffsDB1bsg[LO]))(i) + (*(WilsonCoeffsDB1bsg[NLO]))(i) + (*(WilsonCoeffsDB1bsg[NNLO]))(i);
    }
    C_8G = (*(WilsonCoeffsDB1bsg[LO]))(7) + (*(WilsonCoeffsDB1bsg[NLO]))(7) + (*(WilsonCoeffsDB1bsg[NNLO]))(7);

    for(int i=0; i<=7; i++){
        if(i==6) i++;
        std::cout << "C_" << i << " "
                << (*(WilsonCoeffsDB1bsg[LO]))(i).gslpp::complex::real() << " "
                << (*(WilsonCoeffsDB1bsg[NLO]))(i).gslpp::complex::real() << " "
                << (*(WilsonCoeffsDB1bsg[NNLO]))(i).gslpp::complex::real() << "\n";       
    }
    std::cout << "--------\n" ;
    K_1 = 3. * C(1) * C(1) + 2. * C(1) * C(2);
    K_2 = C(2) * C(2);
    lastInput_computeWilsonCoeffsDB1bsg = currentInput_computeWilsonCoeffsDB1bsg;
    return;
}

gslpp::complex AmpDB2::Gamma21overM21_Bd(orders order) {
    std::cout.precision(4);
    //Marvin Gerlach 2022, Meson width differences and asymmetries 
    orderofp[2] = true;
    if (order == FULLNLO) {
        orderofp[2] = false;
        order = FULLNNLO;
    }
    if (order != FULLNNLO) throw std::runtime_error("AmpDB2::Gamma21overM21_Bd(): order not implemented");
//    if (mySM.getFlavour().getHDF2().getCoeffBd().getOrder() < order % 3)
//        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed");
    
    //get M21 without matrix element
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffBd(
            mySM.getBBd().getMu(),
            mySM.getBBd().getScheme());
    gslpp::complex M21overme0 = ((*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0));
    
    computeCKMandMasses(NNLO);
    computeWilsonCoeffsMisiak();
    lambda_c = mySM.getCKM().getV_cd().conjugate() * mySM.getCKM().getV_cb();
    lambda_u = mySM.getCKM().getV_ud().conjugate() * mySM.getCKM().getV_ub();

    //auto t2 = std::chrono::high_resolution_clock::now();
    compute_pp_s();
    //auto t3 = std::chrono::high_resolution_clock::now();

    compute_matrixelements(d); //compute_matrixelements(s);
    compute_deltas_1overm(d); //compute_deltas_1overm(s);
    
   //equation (6.1)
    gslpp::complex Gamma21overM21_Bd = Gf2 * Mb2 / (24 * M_PI * MB) / M21overme0 *
            ((H(0) + H(1) * me(1)/me(0) + H(2) * me1tilde/me(0)).conjugate() + delta_1overm(d)/me(0));

    //std::cout << "Computes_d" << (std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2)).count() << "\n";
    return Gamma21overM21_Bd;        
}


gslpp::complex AmpDB2::Gamma21overM21_Bs(orders order) {
    //Marvin Gerlach 2022, Meson width differences and asymmetries 
    std::cout.precision(4);    
    orderofp[2] = true;
    if (order == FULLNLO) {
        orderofp[2] = false;
        order = FULLNNLO;
    }
    if (order != FULLNNLO) throw std::runtime_error("AmpDB2::Gamma21overM21_Bs(): order not implemented");
//    if (mySM.getFlavour().getHDF2().getCoeffBs().getOrder() < order % 3)
//        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed");
    
    //get M21 without matrix element
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffBs(
            mySM.getBBs().getMu(),
            mySM.getBBs().getScheme());
    gslpp::complex M21overme0 = ((*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0));
    
    computeCKMandMasses(NNLO);
    computeWilsonCoeffsMisiak();
    lambda_c = mySM.getCKM().getV_cs().conjugate() * mySM.getCKM().getV_cb();
    lambda_u = mySM.getCKM().getV_us().conjugate() * mySM.getCKM().getV_ub();

    compute_pp_s();

    compute_matrixelements(s);
    compute_deltas_1overm(s);
    
   //equation (6.1)
    gslpp::complex Gamma21overM21_Bs = Gf2 * Mb2 / (24 * M_PI * MB_s) / M21overme0 *
            ((H(0) + H(1) * me(1)/me(0) + H(2) * me1tilde/me(0)).conjugate() + delta_1overm(s)/me(0));

    //std::cout << "Computes_d" << (std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2)).count() << "\n";
    return Gamma21overM21_Bs;
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


gslpp::complex AmpDB2::Gamma21overM21_BdFULLNLO1(){
    //hep-ph/0308029v2
    std::cout.precision(4);

    //get M21 without matrix element
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffBd(
            mySM.getBBd().getMu(),
            mySM.getBBd().getScheme());
    gslpp::complex M21overme0 = ((*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0));
         
    computeCKMandMasses(NLO);
    computeWilsonCoeffsDB1bsg();  
    computeF0();
    computeF1();
    computeP();
    computeD();
    compute_matrixelements(d);
    compute_deltas_1overm(d);
    
    //std::cout << "Test " << M21overme0 / (M21_Bd(FULLNLO)/me(0)) << "\n";
    
    //compare with hep-ph/0307344
//    gslpp::complex VudVub = mySM.getCKM().getV_ud().conjugate() * mySM.getCKM().getV_ub();
//
//    std::cout << "a " <<
//            -Gf2 * Mb2 / (24 * M_PI * MB) / M21overme0 *
//                (VtbVtd2 * (D(cc, 0) +  D(cc, 1) * me(1)/me(0) + deltas_1overm(cc, d)/me(0))).conjugate() << "\n";
//    std::cout << "b " <<
//            -Gf2 * Mb2 / (24 * M_PI * MB) / M21overme0 *
//                (2. * VudVub * VtbVtd * (D(cu, 0)-D(cc, 0) +  (D(cu, 1)-D(cc, 1)) * me(1)/me(0) + (deltas_1overm(cu, d)-deltas_1overm(cc, d))/me(0))).conjugate() << "\n";
//    std::cout << "c " <<
//            -Gf2 * Mb2 / (24 * M_PI * MB) / M21overme0 *
//                (VudVub*VudVub * (2.*D(cu, 0)-D(cc, 0)-D(uu, 0)) +  (2.*D(cu, 1)-D(cc, 1)-D(uu, 1)) * me(1)/me(0) + (2.*deltas_1overm(cu,d)-deltas_1overm(cc,d)-deltas_1overm(uu,d))/me(0)).conjugate() << "\n";
        
    //equation 16 divided by 12
    gslpp::complex Gamma21overM21_Bd = -Gf2 * Mb2 / (24 * M_PI * MB) / M21overme0 *
                (c(d)(0) + c(d)(1) * me(1)/me(0) + c(d)(2) * me(2)/me(0) + delta_1overm(d)/me(0));
    return Gamma21overM21_Bd;
}

gslpp::complex AmpDB2::Gamma21overM21_BsFULLNLO1(){
    //hep-ph/0308029v2
    std::cout.precision(4);

    //get M21 without matrix element
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffBs(
            mySM.getBBs().getMu(),
            mySM.getBBs().getScheme());
    gslpp::complex M21overme0 = ((*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0));
         
    computeCKMandMasses(NLO);
    computeWilsonCoeffsDB1bsg();    
    computeF0();
    computeF1();
    computeP();
    computeD();
    compute_matrixelements(s);
    compute_deltas_1overm(s);

    //equation 16 divided by 12
    gslpp::complex Gamma21overM21_Bs = -Gf2 * Mb2 / (24 * M_PI * MB_s) / M21overme0 *
                (c(s)(0) + c(s)(1) * me(1)/me(0) + c(s)(2) * me(2)/me(0) + delta_1overm(s)/me(0));
    return Gamma21overM21_Bs;
}

void AmpDB2::computeWilsonCoeffsMisiak(){
    double currentInput_computeWilsonCoeffsMisiak = mu_1;
    if (lastInput_computeWilsonCoeffsMisiak == currentInput_computeWilsonCoeffsMisiak) return;
    gslpp::vector<gslpp::complex> ** WilsonCoeffsMisiak = mySM.getFlavour().ComputeCoeffsgamma(mu_1);
    for (int i = 0; i < 6; i++) {
        cacheC[i] = (*(WilsonCoeffsMisiak[LO]))(i) + (*(WilsonCoeffsMisiak[NLO]))(i) + (*(WilsonCoeffsMisiak[NNLO]))(i);
    }
    C_8G = (*(WilsonCoeffsMisiak[LO]))(7) + (*(WilsonCoeffsMisiak[NLO]))(7) + (*(WilsonCoeffsMisiak[NNLO]))(7);

    for(int i=0; i<=7; i++){
        if(i==6) i++;
        std::cout << "C_" << i << " "
                << (*(WilsonCoeffsMisiak[LO]))(i).gslpp::complex::real() << " "
                << (*(WilsonCoeffsMisiak[NLO]))(i).gslpp::complex::real() << " "
                << (*(WilsonCoeffsMisiak[NNLO]))(i).gslpp::complex::real() << "\n";       
    }
    std::cout << "--------\n" ;
    K_1 = 3. * C(1) * C(1) + 2. * C(1) * C(2);
    K_2 = C(2) * C(2);
    lastInput_computeWilsonCoeffsMisiak = currentInput_computeWilsonCoeffsMisiak;
}


//equation (6.2)
gslpp::complex AmpDB2::H(){
    gslpp::vector< complex > H(3, 0.);
    H.assign(0, -lambda_c*lambda_c * H(cc) - 2. * lambda_c*lambda_u * H(cu) - lambda_u*lambda_u * H(uu));
    H.assign(2, -lambda_c*lambda_c * H_s(cc) - 2. * lambda_c*lambda_u * H_s(cu) - lambda_u*lambda_u * H_s(uu));
    //change to RI
    //H += as_4pi * coeffsMStoRI.transpose() * H;
    return H;
}

//equation (6.4)
gslpp::complex AmpDB2::H(quarks qq){
    gslpp::complex result = 0.;
    for (int i=1; i<=8; i++){
        if (i==7) i++;
        for (int j=i; j<=8; j++){
            if(j==7) j++;
            result += C(i) * C(j) * p(qq, i, j);
        }
    }
    return result;
}
gslpp::complex AmpDB2::H_s(quarks qq){
    gslpp::complex result = 0.;
    for (int i=1; i<=8; i++){
        if (i==7) i++;
        for (int j=i; j<=8; j++){
            if(j==7) j++;
            result += C(i) * C(j) * p_s(qq, i, j);
        }
    }
    return result;
}

//equation (6.5)
double AmpDB2::p(quarks qq, int i, int j){
    return p(qq, i, j, 0) * orderofp[0] + 
            as_4pi * p(qq, i, j, 1) * orderofp[1] +
            as_4pi * as_4pi * p(qq, i, j, 2) * orderofp[2];
            //+ as_4pi * as_4pi * as_4pi * p(qq, i, j, 3);
}
double AmpDB2::p_s(quarks qq, int i, int j){
    return p_s(qq, i, j, 0) * orderofp[0] + 
            as_4pi * p_s(qq, i, j, 1) * orderofp[1] +
            as_4pi * as_4pi * p_s(qq, i, j, 2) * orderofp[2];
            //+ as_4pi * as_4pi * as_4pi * p_s(qq, i, j, 3);
}

double AmpDB2::p(quarks qq, int i, int j, int n){
    return cache_p[index_p(qq, i, j, n)];
}
double AmpDB2::p_s(quarks qq, int i, int j, int n){
    return cache_ps[index_p(qq, i, j, n)];    
}

int AmpDB2::index_p(quarks qq, int i, int j, int n){
    return n * 192 + qq * 64 + (i - 1) * 8 + (j - 1);
}


void AmpDB2::compute_pp_s(){
    //input didn't change nothing to compute
    double currentInput_compute_pp_s[3] = {z, mu_1, mu_2};
    if (lastInput_compute_pp_s == currentInput_compute_pp_s) return;
    
    //remember value of z after setting it to 0 for calculation of uu coefficients
    double cache_z = z;
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2)) {    
        //equation (6.6)
        cache_p[index_p(qq, 1, 1, 0)]= 23./72. - 11./6. * z;
        cache_p[index_p(qq, 1, 2, 0)]= 1./6. - 2. * z;
        cache_p[index_p(qq, 2, 2, 0)]= 1. - 3. * z;
        cache_ps[index_p(qq, 1, 1, 0)]= -5./9.;
        cache_ps[index_p(qq, 1, 2, 0)]= -4./3.;
        cache_ps[index_p(qq, 2, 2, 0)]= 1.;

        //equation (6.7)
        cache_p[index_p(qq, 1, 3, 0)]= 4./3.;
        cache_p[index_p(qq, 1, 4, 0)]= -5./36.;
        cache_p[index_p(qq, 1, 5, 0)]= 64./3. - 96. * z;
        cache_p[index_p(qq, 1, 6, 0)]= 4. * z - 20./9.;
        cache_p[index_p(qq, 2, 3, 0)]= 1.;
        cache_p[index_p(qq, 2, 4, 0)]= 5./6.;
        cache_p[index_p(qq, 2, 5, 0)]= 16. - 72. * z;
        cache_p[index_p(qq, 2, 6, 0)]= 40./3. - 24. * z;

        cache_ps[index_p(qq, 1, 3, 0)]= -8./3.;
        cache_ps[index_p(qq, 1, 4, 0)]= -2./9.;
        cache_ps[index_p(qq, 1, 5, 0)]= -128./3.;
        cache_ps[index_p(qq, 1, 6, 0)]= -32./9.;
        cache_ps[index_p(qq, 2, 3, 0)]= -2.;
        cache_ps[index_p(qq, 2, 4, 0)]= 4./3.;
        cache_ps[index_p(qq, 2, 5, 0)]= -32.;
        cache_ps[index_p(qq, 2, 6, 0)]= 64./3.;
        
        //equation (6.9)
        z = 0;
    }
    z = cache_z;
    
    //equation (6.8)
    double n_l = 3.;
    double n_v = 1.;
    cache_p[index_p(cc, 3, 3, 0)]= 3. * (n_l + n_v) + 2.;
    cache_p[index_p(cc, 3, 4, 0)]= 7./3.;
    cache_p[index_p(cc, 3, 5, 0)]= 60. * (n_l + n_v) + 64.;
    cache_p[index_p(cc, 3, 6, 0)]= 112./3;
    cache_p[index_p(cc, 4, 4, 0)]= 5. * (n_l + n_v) / 12. + 13./72.;
    cache_p[index_p(cc, 4, 5, 0)]= 112./3.;
    cache_p[index_p(cc, 4, 6, 0)]= 25./3. * (n_l + n_v) + 52./9.;
    cache_p[index_p(cc, 5, 5, 0)]= -1296. * n_v * z + 408. * (n_l + n_v) + 512.;
    cache_p[index_p(cc, 5, 6, 0)]= 1792./3.;
    cache_p[index_p(cc, 6, 6, 0)]= -72. * n_v * z + 170./3. * (n_l + n_v) + 416./9.;
    
    cache_ps[index_p(cc, 3, 3, 0)]= -6. * (n_l + n_v) - 1.;
    cache_ps[index_p(cc, 3, 4, 0)]= -8./3;
    cache_ps[index_p(cc, 3, 5, 0)]= -120. * (n_l + n_v) - 32.;
    cache_ps[index_p(cc, 3, 6, 0)]= -128./3.;
    cache_ps[index_p(cc, 4, 4, 0)]= 2./3. * (n_l + n_v) - 7./9.;
    cache_ps[index_p(cc, 4, 5, 0)]= -128./3;
    cache_ps[index_p(cc, 4, 6, 0)]= 40./3. * (n_l + n_v) - 224./9.;
    cache_ps[index_p(cc, 5, 5, 0)]= -816. * (n_l + n_v) - 256.;
    cache_ps[index_p(cc, 5, 6, 0)]= -2048./3.;
    cache_ps[index_p(cc, 6, 6, 0)]= 272./3. * (n_l + n_v) - 1792./9.;
    
    //equation (6.9)    
    for (int i=3; i<=6; i++){
        for (int j=i; j<=6; j++){
            cache_p[index_p(uu, i, j, 0)] = cache_p[index_p(cc, i, j, 0)];
            cache_ps[index_p(uu, i, j, 0)] = cache_ps[index_p(cc, i, j, 0)];            
        }
    }
    
    //equation (6.10)
    for (int i=1; i<=6; i++){
        for (int j=i; j<=6; j++){
            cache_p[index_p(cu, i, j, 0)] =
                    0.5 * (cache_p[index_p(cc, i, j, 0)] +cache_p[index_p(uu, i, j, 0)]);
            cache_ps[index_p(cu, i, j, 0)] =
                    0.5 * (cache_ps[index_p(cc, i, j, 0)] +cache_ps[index_p(uu, i, j, 0)]);
        }
    }
    
    double L_1 = 2. * log(mu_1/Mb);
    double L_2 = 2. * log(mu_2/Mb);
    double sqrt3 = sqrt(3);
    
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2)) {    
        //equation (6.11)
        cache_p[index_p(qq, 1, 1, 1)]= z * (-14./3. * L_1 - 11./3. * L_2 - 44./3. * logz + M_PI2/54. - 4133./216.)
                + 337./324. * L_1 + 149./108. * L_2 - 5./108. * M_PI2 + 1789./486.;
        cache_p[index_p(qq, 1, 2, 1)]= z * (26. * L_1 - 4. * L_2 - 16. * logz - 2./9. * M_PI2 + 1199./18.)
                - 323./54. * L_1 + 19./9. * L_2 - 5./9. * M_PI2 + 1346./81.;
        cache_p[index_p(qq, 2, 2, 1)]= z * (12. * L_1 - 6. * L_2 - 24. * logz + 2./3. * M_PI2 + 115./6.)
                - 14./9. * L_1 + 2./3. * L_2 - 5./3. * M_PI2 + 91./54.;
        cache_ps[index_p(qq, 1, 1, 1)]= -38./81. * L_1 - 40./27. * L_2 + (-1159./27. - 4./27. * M_PI2) * z - 2./27. * M_PI2 + 2./243.;
        cache_ps[index_p(qq, 1, 2, 1)]= 44./27. * L_1 - 32./9. * L_2 + (16./9. * M_PI2 - 656./9.) * z - 8./9. * M_PI2 + 280./81.;
        cache_ps[index_p(qq, 2, 2, 1)]= 64./9. * L_1 + 8./3. * L_2 + (116./3. - 16./3. * M_PI2) * z - 8./3. * M_PI2 + 728./27.;
        //equation (6.14)
        cache_p[index_p(qq, 1, 3, 1)]= (320./9. - 4. * L_1) * z + 47./18. * L_1 + 56./9. * L_2 - 5./(18. * sqrt3) * M_PI + 1523./108.;
        cache_p[index_p(qq, 1, 4, 1)]= (59./3. * L_1 + 5./9. * M_PI2 + 4565./108.) * z - 281./108. * L_1 + L_2 / 54. + 5./18. * M_PI2
                - 25./(108. * sqrt3) * M_PI - 712./81.;
        cache_p[index_p(qq, 1, 5, 1)]= z * (-136. * L_1 - 192. * L_2 - 768. * logz - 16408./9.)
                + 376./9. * L_1 + 896./9. * L_2 - 40./(9. * sqrt3) * M_PI + 318.;
        cache_p[index_p(qq, 1, 6, 1)]= z * (764./3. * L_1 + 8. * L_2 + 32. * logz + 8./9. * M_PI2 + 22850./27.)
                - 1259./27. * L_1 + 8./27. * L_2 + 40./9. * M_PI2 - 55./(27. * sqrt3) * M_PI - 4243./27.;
        cache_p[index_p(qq, 2, 3, 1)]= (24. * L_1 + 170./3.) * z - 47./3. * L_1 + 14./3. * L_2 + 5./(3. * sqrt3) * M_PI - 677./18;
        cache_p[index_p(qq, 2, 4, 1)]= (26. * L_1 - 10./3. * M_PI2 + 1429./18.) * z - 35./9. * L_1 - L_2 / 9. - 5./3. * M_PI2
                + 25./(18. * sqrt3) * M_PI - 88./27.;
        cache_p[index_p(qq, 2, 5, 1)]= z * (816. * L_1 - 144. * L_2 - 576 * logz + 3656./3.)
                - 752./3. * L_1 + 224./3. * L_2 + 80./(3. * sqrt3) * M_PI - 580.;
        cache_p[index_p(qq, 2, 6, 1)]= z * (128. * L_1 - 48. * L_2 - 192. * logz - 16./3. * M_PI2 + 6140./9.)
                - 290./9. * L_1 - 16./9. * L_2 - 80./3. * M_PI2 + 110./(9. * sqrt3) * M_PI - 442./9.;

        cache_ps[index_p(qq, 1, 3, 1)]= -4./3. * L_1 - 64./9. * L_2 - 1720./9. * z - 4./(9. * sqrt3) * M_PI - 130./27.;
        cache_ps[index_p(qq, 1, 4, 1)]= 2. * L_1 - 16./27. * L_2 + (80./27. + 8./9. * M_PI2) * z + 4./9. * M_PI2 - 10./(27. * sqrt3) * M_PI + 404./81.;
        cache_ps[index_p(qq, 1, 5, 1)]= -64./3. * L_1 - 1024./9. * L_2 - 27952./9. * z - 64./(9 * sqrt3) * M_PI - 2128./9.;
        cache_ps[index_p(qq, 1, 6, 1)]= 24. * L_1 - 256./27. * L_2 + (128./9. * M_PI2 - 520./27.) * z - 64./9. * M_PI2
                - 88./(27. * sqrt3) * M_PI + 2824./27.;
        cache_ps[index_p(qq, 2, 3, 1)]= 8. * L_1 - 16./3. * L_2 - 448./3. * z + 8./(3. * sqrt3) * M_PI + 116./9.;
        cache_ps[index_p(qq, 2, 4, 1)]= 32./9. * L_2 + (488./9. - 16./3. * M_PI2) * z - 8./3. * M_PI2 + 20./(9. * sqrt3) * M_PI + 272./27.;
        cache_ps[index_p(qq, 2, 5, 1)]= 128. * L_1 - 256./3. * L_2 - 6304./3. * z + 128./(3. * sqrt3) * M_PI + 32./3.;
        cache_ps[index_p(qq, 2, 6, 1)]= 48. * L_1 + 512./9. * L_2 + (7520./9. - 256./3. * M_PI2) * z
                - 128./3. * M_PI2 + 176./(9. * sqrt3) * M_PI + 1840./9.;

        //equation (6.13) and equation (6.15)
        z = 0;
    }
    z = cache_z;
    //equation (6.15)
    cache_p[index_p(uu, 1, 4, 1)] += 5./9. * z;
    cache_p[index_p(uu, 1, 6, 1)] += 50./9. * z;
    cache_p[index_p(uu, 2, 4, 1)] += - 10./3. * z;
    cache_p[index_p(uu, 2, 6, 1)] += - 100./3. * z;

    cache_ps[index_p(uu, 1, 4, 1)] += 8./9. * z;
    cache_ps[index_p(uu, 1, 6, 1)] += 80./9. * z;
    cache_ps[index_p(uu, 2, 4, 1)] += - 16./3. * z;
    cache_ps[index_p(uu, 2, 6, 1)] += - 160./3. * z; 
    
    //equation (6.13) and (6.16)
    for (int i=1; i<=2; i++){
        for (int j=i; j<=6; j++){
            cache_p[index_p(cu, i, j, 1)] =
                    0.5 * (cache_p[index_p(cc, i, j, 1)] + cache_p[index_p(uu, i, j, 1)]);
            cache_ps[index_p(cu, i, j, 1)] =
                    0.5 * (cache_ps[index_p(cc, i, j, 1)] + cache_ps[index_p(uu, i, j, 1)]);
        }
    }
    
    //equation (6.18)
    cache_p[index_p(cc, 3, 3, 1)] = -154./9. * L_1 + 184./3. * L_2 + 90. * z - 5./3. * M_PI2 + 5./(3. * sqrt3) * M_PI + 1390./27.;
    cache_p[index_p(cc, 3, 4, 1)] = -811./54. * L_1 + 74./9. * L_2 - 10./3. * z - 10./9. * M_PI2 + 70./(9. * sqrt3) * M_PI - 27991./324.;
    cache_p[index_p(cc, 3, 5, 1)] = -4928./9. * L_1 + 3872./3. * L_2 + 1800. * z - 160./3. * M_PI2 + 160./(3. * sqrt3) * M_PI
            + 16880./27.;
    cache_p[index_p(cc, 3, 6, 1)] = (144. * L_1 + 440./3.) * z - 12932./27. * L_1 + 1184./9. * L_2 - 160./9. * M_PI2
            + 670./(9. * sqrt3) * M_PI - 131410./81.;
    cache_p[index_p(cc, 4, 4, 1)] = 181./162. * L_1 + 127./108. * L_2 + (323./36. - 5./3. * M_PI2) * z - 335./108. * M_PI2
            + 575/(108. * sqrt3) * M_PI + 779./486.;
    cache_p[index_p(cc, 4, 5, 1)] = (576. * L_1 + 3836./3.) * z - 14912./27. * L_1 + 1184./9. * L_2 - 160./9. * M_PI2
            + 1120./(9. * sqrt3) * M_PI - 127990./81.;
    cache_p[index_p(cc, 4, 6, 1)] = (60. * L_1 - 100./3. * M_PI2 + 2455./9.) * z - 8759./81. * L_1 + 1088./27. * L_2
            - 1600./27. * M_PI2 + 2665./(27. * sqrt3) * M_PI - 50083./243.;
    cache_p[index_p(cc, 5, 5, 1)] = z * (-2592. * L_2 - 10368. * logz - 33120.) - 39424./9. * L_1
            + 26944./3. * L_2 - 1280./3. * M_PI2 + 1280./(3. * sqrt3) * M_PI + 347104./27.;
    cache_p[index_p(cc, 5, 6, 1)] = (7200. * L_1 + 74000./3.) * z - 240608./27. * L_1 + 18944./9. * L_2
            - 2560./9. * M_PI2 + 10720./(9. * sqrt3) * M_PI - 2253568./81.;
    cache_p[index_p(cc, 6, 6, 1)] = z * (-48. * L_1 - 144. * L_2 - 576. * logz - 248./3. * M_PI2 + 12290./9.)
            - 59632./81. * L_1 + 8848./27. * L_2 - 10640./27. * M_PI2 + 12320./(27. * sqrt3) * M_PI - 662144./243.;
    
    cache_ps[index_p(cc, 3, 3, 1)] = 176./9. * L_1 - 200./3. * L_2 - 432. * z - 8./3. * M_PI2 + 8./(3. * sqrt3) * M_PI - 620./27.;
    cache_ps[index_p(cc, 3, 4, 1)] = 268./27. * L_1 - 64./9. * L_2 - 16./3. * z - 16./9. * M_PI2 + 112./(9. * sqrt3) * M_PI + 3506./81.;
    cache_ps[index_p(cc, 3, 5, 1)] = 5632./9. * L_1 - 4096./3. * L_2 - 8640. * z - 256./3. * M_PI2 + 256./(3. * sqrt3) * M_PI + 9728./27.;
    cache_ps[index_p(cc, 3, 6, 1)] = 9184./27. * L_1 - 1024./9. * L_2 - 160./3. * z - 256./9. * M_PI2 + 1072./(9. * sqrt3) * M_PI + 88688./81.;
    cache_ps[index_p(cc, 4, 4, 1)] = 1028./81. * L_1 + 136./27. * L_2 - 8./3. * M_PI2 * z + 230./9. * z - 134./27 * M_PI2 + 230./(27. * sqrt3) * M_PI
            + 6214./243.;
    cache_ps[index_p(cc, 4, 5, 1)] = 9472./27. * L_1 - 1024./9. * L_2 + 608./3. * z - 256./9. * M_PI2 + 1792./(9. * sqrt3) + 64784./81.;
    cache_ps[index_p(cc, 4, 6, 1)] = 10792./81. * L_1 + 2048./27. * L_2 - 160./3. * M_PI2 * z + 3568./9. * z - 2560./27. * M_PI2
            + 4264./(27. * sqrt3) * M_PI + 123080./243.;
    cache_ps[index_p(cc, 5, 5, 1)] = 45056./9. * L_1 - 28160./3. * L_2 - 58752. * z - 2048./3. * M_PI2 - 2048./(3. * sqrt3) * M_PI
            - 349184./27.;
    cache_ps[index_p(cc, 5, 6, 1)] = 167680./27. * L_1 - 16384./9. * L_2 + 6080./3. * z - 4096./9. * M_PI2 + 17152./(9. * sqrt3) * M_PI
            + 1502720./81.;
    cache_ps[index_p(cc, 6, 6, 1)] = 75392./81. * L_1 + 11776./27. * L_2 - 1088./3. * M_PI2 * z + 23696./9. * z - 17024./27. * M_PI2
            + 19712./(27. * sqrt3) * M_PI + 717184./243.;
    
    //equation (6.17)
    for (int i=3; i<=6; i++){
        for (int j=i; j<=6; j++){
            cache_p[index_p(cu, i, j, 1)] = cache_p[index_p(cc, i, j, 1)];
            cache_p[index_p(uu, i, j, 1)] = cache_p[index_p(cc, i, j, 1)];
            cache_ps[index_p(cu, i, j, 1)] = cache_ps[index_p(cc, i, j, 1)];
            cache_ps[index_p(uu, i, j, 1)] = cache_ps[index_p(cc, i, j, 1)];
        }
    }
    
    //equation (6.19)
    cache_p[index_p(cc, 1, 8, 1)] = 5./18.;
    cache_p[index_p(cc, 2, 8, 1)] = -5./3.;
    cache_ps[index_p(cc, 1, 8, 1)] = 4./9.;
    cache_ps[index_p(cc, 2, 8, 1)] = -8./3.;
    
    //equation (6.21)
    cache_p[index_p(cc, 3, 8, 1)] = -32./3.;
    cache_p[index_p(cc, 4, 8, 1)] = -169./18.;
    cache_p[index_p(cc, 5, 8, 1)] = -512./3.;
    cache_p[index_p(cc, 6, 8, 1)] = -992./9.;
    cache_ps[index_p(cc, 3, 8, 1)] = 64./3.;
    cache_ps[index_p(cc, 4, 8, 1)] = -20./9.;
    cache_ps[index_p(cc, 5, 8, 1)] = 1024./3.;
    cache_ps[index_p(cc, 6, 8, 1)] = 256./9.;
    
    //equation (6.20)
    for (int i=1; i<=6; i++){
        cache_p[index_p(cu, i, 8, 1)] = cache_p[index_p(cc, i, 8, 1)];
        cache_p[index_p(uu, i, 8, 1)] = cache_p[index_p(cc, i, 8, 1)];
        cache_ps[index_p(cu, i, 8, 1)] = cache_ps[index_p(cc, i, 8, 1)];
        cache_ps[index_p(uu, i, 8, 1)] = cache_ps[index_p(cc, i, 8, 1)];
    }
    
    double L_12 = L_1 * L_1;
    double L_22 = L_2 * L_2;
    
    //equations (3.102, 6.23, 6.25)
    double zeta3 = 1.20206;
    double t_2 = -0.389012;
    double Cl2PI3 = 1.014941;
    
    double sqrt5 = sqrt(5);
    double log2z = logz * logz;
    double log2 = log(2);
    double log3 = log(3);
    double log12sqrt52 = log(0.5 + sqrt5/2.);
    double sqrtz = sqrt(z);
    
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2)) {
        //equation (6.22)
        cache_p[index_p(qq, 1, 1, 2)] = z * (-1348./9. * L_1 * logz - 88./3. * L_2  * logz - 2347./54. * L_12 + 187./18. * L_22
                + 31./54. * M_PI2 * L_1 - 722039./1944. * L_1 - 337./9. * L_1 * L_2 + 19./81. * M_PI2 * L_2 + 1891./81. * L_2
                + 22./9. * log2z + 4./27. * M_PI2 * logz - 1591./3. * logz + 128581./216. * zeta3
                - 13637./116640. * M_PI4 + 203./81. * sqrt5 * M_PI2 + 235469./3888. * M_PI2 - 25./(162. * sqrt3) * M_PI
                - 601385353./583200. + 176./27. * M_PI2 * log2 - 4321./324. * M_PI2 * log2
                - 68./27. * M_PI2 * log12sqrt52) + 3211./324. * L_12 + 12911./972. * L_2 * L_1 - 5./4. * M_PI2 * L_1
                - 25./(972. * sqrt3) * M_PI * L_1 + 1320817./17496. * L_1 - 311./216. * L_22 + M_PI2 * L_2 /162. + 259603./11664. * L_2
                + 5./(162. * sqrt3) * t_2 + 28333./486. * zeta3 + 23./4860. * M_PI4 - 2197./972 * sqrt5 * M_PI2 - 216641./69984. * M_PI2
                - 25./(1458. * sqrt3) * M_PI + 814589597./4199040. - 71./972. * M_PI2 * log2 - 5./(1944. * sqrt3) * M_PI * log3
                - 169./81. * M_PI2 * log12sqrt52 * 56./(243. * sqrt3) * Cl2PI3
                + n_v * (0.0617284 * L_1 * z + 4.60031 * z + 4.75206 * sqrtz - 5.55556 * z * logz);
        cache_p[index_p(qq, 1, 2, 2)] = z * (256./3. * L_1 * logz - 32. * L_2 * logz + 1193./9. * L_12 + 34./3. * L_22
                - 44./9. * M_PI2 * L_1 + 117563./162. * L_1 + 64./3. * L_1 * L_2 - 76./27. * M_PI2 * L_2 + 6350./27. * L_2
                + 8./3. * log2z - 16./9. * M_PI2 * logz + 364./3. * logz + 85027./90. * zeta3
                + 20833./4860. * M_PI4 + 548./(27. * sqrt5) * M_PI2 - 11245./162. * M_PI2 + 50./(27. * sqrt3) * M_PI + 12685151./9720.
                + 64./9. * M_PI2 * log2 - 1361./27. * M_PI2 * log2 - 176./45. * M_PI2 * log12sqrt52)
                - 1751./54. * L_12 + 166./81. * L_2 * L_1 + 10. * M_PI2 * L_1 + 25./(81. * sqrt3) * M_PI * L_1 - 1026907./5832. * L_1
                - L_22/18. - 2./27. * M_PI2 * L_2 - 619./972. * L_2 - 10./(27. * sqrt3) * t_2 + 10573./324. * zeta3 - 799./810. * M_PI4
                - 299./81. * sqrt5 * M_PI2 + 497221./11664. * M_PI2 + 50./(243. * sqrt3) * M_PI - 95740679./349920.
                + 596./81. * M_PI2 * log2 + 5./(162. * sqrt3) * M_PI * log3 - 92./27. * M_PI2 * log12sqrt52
                - 224./(81. * sqrt3) * Cl2PI3 + n_v * (-0.740741 * L_1 * z - 85.8705 * z
                + 48.2515 * sqrtz + 2.66667 * z * logz);
        cache_p[index_p(qq, 2, 2, 2)] = z * (-88. * L_1 * logz - 48. * L_2 * logz - 122./3. * L_12 + 17. * L_22 + 26./3. * M_PI2 * L_1
                - 16583./54. * L_1 - 22. * L_1 * L_2 + 76./9. * M_PI2 * L_2 - 109./18. * L_2 + 4. * log2z
                + 16./3. * M_PI2 * logz - 464. * logz + 2521./15. * zeta3 + 5203./3240 * M_PI4
                + 28./(9. * sqrt5) * M_PI2 + 7097./108. * M_PI2 - 50./(9. * sqrt3) * M_PI - 12332857./16200. + 32./3. * M_PI2 * log2
                - 274./9. * M_PI2 * log2 - 16./15. * M_PI2 * log12sqrt52) + 239./18. * L_12
                - 202./27. * L_2 * L_1 - 15. * M_PI2 * L_1 - 25./(27. * sqrt3) * M_PI * L_1 + 106199./972. * L_1 - 19./3. * L_22
                + 2./9. * M_PI2 * L_2 - 5117./81. * L_2 + 10./(9. * sqrt3) * t_2 - 3157./54. * zeta3 + 971./540. * M_PI4
                - 13./27. * sqrt5 * M_PI2 - 177247./3888. * M_PI2 - 50./(81. * sqrt3) * M_PI + 74041./14580. + 148./27. * M_PI2 * log2
                + 5./(54. * sqrt3) * M_PI * log3 - 4./9. * M_PI2 * log12sqrt52 + 224./(27. * sqrt3) * Cl2PI3
                + n_v * (2.22222 * L_1 * z + 70.6121 * z - 105.276 * sqrtz - 32. * z * logz);
        cache_ps[index_p(qq, 1, 1, 2)] = z * (-4. * M_PI2 * L_1 - 98023./243. * L_1 - 32./81. * M_PI2 * L_2 - 9272./81. * L_2
                - 32./27. * M_PI2 * logz - 9272./27. * logz + 29./3. * zeta3 - 27529./14580. * M_PI4
                - 344./81. * sqrt5 * M_PI2 - 7103./486 * M_PI2 - 20./(81. * sqrt3) * M_PI - 33198263./36450.
                + 1826./81. * M_PI2 * log2) - 902./243. * L_12 - 3064./243. * L_2 * L_1 - 2. * M_PI2 * L_1
                - 10./(243. * sqrt3) * M_PI * L_1 - 77617./2187. * L_1 + 260./27. * L_22 - 16./81. * M_PI2 * L_2 - 5504./729 * L_2
                + 4./(81. * sqrt3) * t_2 + 28528./243. * zeta3 + 449./1215. * M_PI4 + 1118./243. * sqrt5 * M_PI2 - 44209./8748. * M_PI2
                - 20./(729. * sqrt3) * M_PI - 67489177./262440. - 3506./243. * M_PI2 * log2 - M_PI * log3 / (243. * sqrt3)
                + 344./81. * M_PI2 * log12sqrt52 + 104./(243. * sqrt3) * Cl2PI3
                + n_v * (0.0987654 * L_1 * z - 26.8617 * z + 27.7812 * sqrtz);
        cache_ps[index_p(qq, 1, 2, 2)] = z * (32. * M_PI2 * L_1 - 23276./81. * L_1 + 128./27. * M_PI2 * L_2 - 5248./27. * L_2
                + 128./9. * M_PI2 * logz - 5248./9. * logz - 244. * zeta3 + 5692./1215. * M_PI4
                - 160./27. * sqrt5 * M_PI2 + 3238./81. * M_PI2 + 80./(27. * sqrt3) * M_PI + 5060009./12150.
                + 208./27. * M_PI2 * log2) + 44./81. * L_12 - 1856./81. * L_2 * L_1 + 16. * M_PI2 * L_1 + 40./(81. * sqrt3) * M_PI * L_1
                - 28733./729. * L_1 + 208./9. * L_22 + 64./27. * M_PI2 * L_2 - 2176./243. * L_2 - 16./(27. * sqrt3) * t_2
                + 13934./81. * zeta3 + 226./405. * M_PI4 + 520./81. * sqrt5 * M_PI2 + 39995./1458. * M_PI2 + 80./(243. * sqrt3) * M_PI
                - 1336127./2187. - 1624./81. * M_PI2 * log2 + 4./(81. * sqrt3) * M_PI * log3
                + 160./27. * M_PI2 * log12sqrt52 - 416./(81. * sqrt3) * Cl2PI3
                + n_v * (-1.18519 * L_1 * z - 72.3265 * z + 87.73 * sqrtz);
        cache_ps[index_p(qq, 2, 2, 2)] = z * (-48. * M_PI2 * L_1 + 18740./27. * L_1 - 128./9. * M_PI2 * L_2 + 928./9. * L_2
                - 128./3. * M_PI2 * logz + 928./3. * logz - 600. * zeta3 + 7991./405 * M_PI4
                - 32./9. * sqrt5 * M_PI2 - 8038./27. * M_PI2 - 80./(9. * sqrt3) * M_PI + 6836747./2025.
                + 272./9. * M_PI2 * log2) + 604./27. * L_12 + 1064./27. * L_2 * L_1 - 24. * M_PI2 * L_1 - 40./(27. * sqrt3) * M_PI * L_1
                + 40370./243. * L_1 - 52./3. * L_22 - 64./9. * M_PI2 * L_2 + 6928./81. * L_2 + 16./(9. * sqrt3) * t_2
                - 4388./27. * zeta3 + 398./135. * M_PI4 + 104./27. * sqrt5 * M_PI2 - 41279./486. * M_PI2 - 80./(81. * sqrt3) * M_PI
                + 27476329./58320. - 656./27. * M_PI2 * log2 - 4./(27. * sqrt3) * M_PI * log3
                + 32./9. * M_PI2 * log12sqrt52 + 416./(27. * sqrt3) * Cl2PI3
                + n_v * (3.55556 * L_1 * z + 176.979 * z -105.276 * sqrtz);
    
        //equation (6.26)
        z = 0;
    }
    z = cache_z;
    
    //equation (6.26)
    cache_p[index_p(uu, 1, 1, 2)] += n_v * (-5.55556 * z * logz + 0.0617284 * L_1 * z + 4.60031 * z + 4.75206 * sqrtz);
    cache_p[index_p(uu, 1, 2, 2)] += n_v * (2.66667 * z * logz - 0.740741 * L_1 * z - 85.8705 * z + 48.2515 * sqrtz);
    cache_p[index_p(uu, 2, 2, 2)] += n_v * (-32. * z * logz + 2.22222 * L_1 * z + 70.6121 * z - 105.276 * sqrtz);
    cache_ps[index_p(uu, 1, 1, 2)] += n_v * (0.0987654 * L_1 * z - 26.8617 * z + 27.7812 * sqrtz);
    cache_ps[index_p(uu, 1, 2, 2)] += n_v * (-1.18519 * L_1 * z - 72.3265 * z + 87.73 * sqrtz);
    cache_ps[index_p(uu, 2, 2, 2)] += n_v * (3.55556 * L_1 * z + 176.979 * z - 105.276 * sqrtz);
    
    //equation (6.27)
    for (int i=1; i<=2; i++){
        for (int j=i; j<=2; j++){
            cache_p[index_p(cu, i, j, 2)] = 0.5 * (cache_p[index_p(cc, i, j, 2)] + cache_p[index_p(uu, i, j, 2)]);
            cache_ps[index_p(cu, i, j, 2)] = 0.5 * (cache_ps[index_p(cc, i, j, 2)] + cache_ps[index_p(uu, i, j, 2)]);
        }
    }
    
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2)) {
        //equation (6.28)
        cache_p[index_p(qq, 1, 8, 2)] = 208./81. * L_1 - L_2/27. + (2615./54. - 10./9. * M_PI2) * z - 5./9. * M_PI2 + 25./(54. * sqrt3) * M_PI
                - 115./486.;
        cache_p[index_p(qq, 2, 8, 2)] = -11./27. * L_1 + 2./9. * L_2 + (20./3. * M_PI2 - 833./9.) * z + 10./3. * M_PI2 - 25./(9. * sqrt3) * M_PI
                - 3125./81.;
        cache_ps[index_p(qq, 1, 8, 2)] = 448./81. * L_1 + 32./27. * L_2 + (1192./27. - 16./9. * M_PI2) * z - 8./9. * M_PI2 + 20./(27. * sqrt3) * M_PI
                + 3580./243.;
        cache_ps[index_p(qq, 2, 8, 2)] = -248./27. * L_1 - 64./9. * L_2 + (32./3. * M_PI2 - 1088./9.) * z + 16./3. * M_PI2 - 40./(9. * sqrt3) * M_PI
                - 4568./81.;
        
        //equation (6.30)
        cache_p[index_p(qq, 3, 8, 2)] = -85./27. * L_1 - 448./9. * L_2 - 196./3. * z + 25./6. * M_PI2 - 107./(18. * sqrt3) * M_PI - 17201./81.;
        cache_p[index_p(qq, 4, 8, 2)] = -3269./162. * L_1 - 427./27. * L_2 + (20./3. * M_PI2 - 404./3.) * z + 169./12. * M_PI2
                -514./(27. * sqrt3) * M_PI - 43016./243.;
        cache_p[index_p(qq, 5, 8, 2)] = 5120./27. * L_1 - 7168./9. * L_2 - 760./3. * z + 770./9. * M_PI2 - 28./(9. * sqrt3) * M_PI - 430238./81.;
        cache_p[index_p(qq, 6, 8, 2)] = -8962./81. * L_1 - 6976./27. * L_2 + (200./3. * M_PI2 - 4222./3.) * z + 3761./27. * M_PI2
                - 3220./(27. * sqrt3) * M_PI - 474656./243.;
        cache_ps[index_p(qq, 3, 8, 2)] = 440./27. * L_1 + 512./9. * L_2 + 608./3. * z - 596./27. * M_PI2 - 52./(9. * sqrt3) * M_PI + 22504./81.;
        cache_ps[index_p(qq, 4, 8, 2)] = -4804./81. * L_1 - 160./27. * L_2 + (32./3. * M_PI2 - 128./3.) * z + 1090./81. * M_PI2
                - 1120./(27. * sqrt3) * M_PI - 46988./243.;
        cache_ps[index_p(qq, 5, 8, 2)] = 17408./27. * L_1 + 8192./9. * L_2 + 11456./3. * z - 8912./27. * M_PI2 - 1984./(9. * sqrt3) * M_PI
                + 420304./81.;
        cache_ps[index_p(qq, 6, 8, 2)] = -28624./81. * L_1 + 2048./27. * L_2 + (320./3. * M_PI2 - 160./3.) * z + 5608./81. * M_PI2
                - 8416./(27. * sqrt3) * M_PI - 423440./243.;
        
        //equation (6.31) and equation (6.29)
        z = 0;
    }
    z = cache_z;
    
    //equation (6.29)
    cache_p[index_p(uu, 1, 8, 2)] += -10./9. * z;
    cache_p[index_p(uu, 2, 8, 2)] += 20./3. * z;
    cache_ps[index_p(uu, 1, 8, 2)] += -16./9. * z;
    cache_ps[index_p(uu, 1, 8, 2)] += 32./3. * z;
    
    //equation (6.31)
    cache_p[index_p(uu, 3, 8, 2)] += -196./3. * z;
    cache_p[index_p(uu, 4, 8, 2)] += (-404./3. + 20./3. * M_PI2) * z;
    cache_p[index_p(uu, 5, 8, 2)] += -760./3. * z;
    cache_p[index_p(uu, 6, 8, 2)] += (-4222./3. + 200./3. * M_PI2) * z;
    cache_ps[index_p(uu, 3, 8, 2)] += 608./3.* z;
    cache_ps[index_p(uu, 4, 8, 2)] += (-128./3. + 32./3. * M_PI2) * z;
    cache_ps[index_p(uu, 5, 8, 2)] += 11456./3.* z;
    cache_ps[index_p(uu, 6, 8, 2)] += (-160./3. + 320./3. * M_PI2) * z;
    
    //as in equation (6.20)
    for (int i=1; i<=6; i++){
        cache_p[index_p(cu, i, 8, 2)] =
                0.5 * (cache_p[index_p(cc, i, 8, 2)] + cache_p[index_p(uu, i, 8, 2)]);
        cache_ps[index_p(cu, i, 8, 2)] =
                0.5 * (cache_ps[index_p(cc, i, 8, 2)] + cache_ps[index_p(uu, i, 8, 2)]);
    }
    
    //equation (6.32)
    cache_p[index_p(cc, 8, 8, 2)] = -13./18.;
    cache_p[index_p(cu, 8, 8, 2)] = -13./18.;
    cache_p[index_p(uu, 8, 8, 2)] = -13./18.;
    cache_ps[index_p(cc, 8, 8, 2)] = -68./9.;
    cache_ps[index_p(cu, 8, 8, 2)] = -68./9.;
    cache_ps[index_p(uu, 8, 8, 2)] = -68./9.;

    lastInput_compute_pp_s[0] = z;
    lastInput_compute_pp_s[1] = mu_1;
    lastInput_compute_pp_s[2] = mu_2;
    return;
}

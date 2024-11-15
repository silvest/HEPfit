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

AmpDB2::AmpDB2(const StandardModel &SM_i, int BMeson_i, bool flag_fixmub, bool flag_RI)
    : mySM(SM_i), meMStoRI(5, 0.), coeffsMStoRI(3, 0.)
{
    BMeson = BMeson_i;
    this->flag_resumz = true;
    this->flag_fixmub = flag_fixmub;
    this->flag_RI = flag_RI;
    switch (BMeson)
    {
    case 0:
        mySM.initializeBParameter("BBd");
        mySM.initializeBParameter("BBd_subleading");   
        mySM.initializeMeson(QCD::B_D);
        break;
    case 1:
        mySM.initializeBParameter("BBs");
        mySM.initializeBParameter("BBs_subleading");
        mySM.initializeMeson(QCD::B_S);
        break;
    default:
        throw std::runtime_error("AmpDB2::AmpDB2(): invalid B meson index");
        break;
    }

    // hep-ph/0606197 eq. 4.7 - 4.10
    double meMStoRI0[5] = {-3. - 5. / 3. + 8. * log2, 0., 0., 0., 0.},
           meMStoRI1[5] = {0., 61. / 9. + 44. / 9. * log2, -7. / 9. + 28. / 9. * log2, 0., 0.},
           meMStoRI2[5] = {0., -25. / 9. + 28. / 9. * log2, -29. / 9. + 44. / 9. * log2, 0., 0.},
           meMStoRI3[5] = {0., 0., 0., -5. / 3. + 13. - 2. / 3. * log2, -3. + 1. + 2. * log2},
           meMStoRI4[5] = {0., 0., 0., -7. / 2. + 11. / 2. + 2. * log2, -1. / 6. - 1. / 2. - 2. / 3. * log2};
    meMStoRI.assign(0, meMStoRI0);
    meMStoRI.assign(1, meMStoRI1);
    meMStoRI.assign(2, meMStoRI2);
    meMStoRI.assign(3, meMStoRI3);
    meMStoRI.assign(4, meMStoRI4);
    for (int i = 0; i <= 2; i++)
    {
        for (int j = 0; j <= 2; j++)
        {
            coeffsMStoRI.assign(i, j, meMStoRI(i, j));
        }
    }
}

gslpp::complex AmpDB2::RBs(orders order)
{
    mySM.getFlavour().getHDF2().getCoeffBs().getOrder();
    if (mySM.getFlavour().getHDF2().getCoeffBs().getOrder() < getHighest(order))
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed");

    gslpp::vector<gslpp::complex> **allcoeff_SM = mySM.getFlavour().ComputeCoeffBs(
        mySM.getBBs().getMu(),
        mySM.getBBs().getScheme(), true);

    C_1_SM = ((*(allcoeff_SM[LO]))(0) + (*(allcoeff_SM[NLO]))(0));

    gslpp::vector<gslpp::complex> **allcoeff = mySM.getFlavour().ComputeCoeffBs(
        mySM.getBBs().getMu(),
        mySM.getBBs().getScheme());

    gslpp::vector<double> me(mySM.getBBs().getBpars());
    double MBs = mySM.getMesons(QCD::B_S).getMass();
    double Mb = mySM.Mrun(mySM.getBBs().getMu(),
                          mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
                          mySM.getQuarks(QCD::BOTTOM).getMass(), QCD::BOTTOM, FULLNNLO);
    double Ms = mySM.Mrun(mySM.getBBs().getMu(),
                          mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                          mySM.getQuarks(QCD::STRANGE).getMass(), QCD::STRANGE,FULLNNLO);
    double KBs = MBs / (Mb + Ms) * MBs / (Mb + Ms);
    double Fbs = mySM.getMesons(QCD::B_S).getDecayconst();
    me(0) *= 1. / 3. * MBs * Fbs * Fbs;
    me(1) *= -5. / 24. * KBs * MBs * Fbs * Fbs;
    me(2) *= 1. / 24. * KBs * MBs * Fbs * Fbs;
    me(3) *= 1. / 4. * KBs * MBs * Fbs * Fbs;
    me(4) *= 1. / 12. * KBs * MBs * Fbs * Fbs;

    /*std::cout << "low scale :" << std::endl << std::endl;
    std::cout << "C1_SM :" << C_1_SM << std::endl << std::endl;

    std::cout << "C1 :" << ((*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0)) << std::endl;
    std::cout << "C2 :" << ((*(allcoeff[LO]))(1) + (*(allcoeff[NLO]))(1)) << std::endl;
    std::cout << "C3 :" << ((*(allcoeff[LO]))(2) + (*(allcoeff[NLO]))(2)) << std::endl;
    std::cout << "C4 :" << ((*(allcoeff[LO]))(3) + (*(allcoeff[NLO]))(3)) << std::endl;
    std::cout << "C5 :" << ((*(allcoeff[LO]))(4) + (*(allcoeff[NLO]))(4)) << std::endl << std::endl;*/

    switch (order)
    {
    case FULLNLO:
        return (*(allcoeff[LO]) + *(allcoeff[NLO])) * me /
               (C_1_SM * me(0));
    case LO:
        return ((*(allcoeff[LO])) * me /(*(allcoeff_SM[LO]))(0) * me(0));
    default:
        throw std::runtime_error("RBs::RBs(): order not implemented");
    }
}

gslpp::complex AmpDB2::RBd(orders order)
{
    mySM.getFlavour().getHDF2().getCoeffBd().getOrder();
    if (mySM.getFlavour().getHDF2().getCoeffBd().getOrder() < getHighest(order))
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed");

    gslpp::vector<gslpp::complex> **allcoeff_SM = mySM.getFlavour().ComputeCoeffBd(
        mySM.getBBd().getMu(),
        mySM.getBBd().getScheme(), true);

    C_1_SM = ((*(allcoeff_SM[LO]))(0) + (*(allcoeff_SM[NLO]))(0));

    gslpp::vector<gslpp::complex> **allcoeff = mySM.getFlavour().ComputeCoeffBd(
        mySM.getBBd().getMu(),
        mySM.getBBd().getScheme());

    gslpp::vector<double> me(mySM.getBBd().getBpars());
    double MBd = mySM.getMesons(QCD::B_D).getMass();
    double Mb = mySM.Mrun(mySM.getBBd().getMu(),
                          mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
                          mySM.getQuarks(QCD::BOTTOM).getMass(), QCD::BOTTOM, FULLNNLO);
    double Md = mySM.Mrun(mySM.getBBd().getMu(),
                          mySM.getQuarks(QCD::DOWN).getMass_scale(),
                          mySM.getQuarks(QCD::DOWN).getMass(), QCD::DOWN, FULLNNLO);
    double KBd = MBd / (Mb + Md) * MBd / (Mb + Md);
    double Fbd = mySM.getMesons(QCD::B_D).getDecayconst();
    me(0) *= 1. / 3. * MBd * Fbd * Fbd;
    me(1) *= -5. / 24. * KBd * MBd * Fbd * Fbd;
    me(2) *= 1. / 24. * KBd * MBd * Fbd * Fbd;
    me(3) *= 1. / 4. * KBd * MBd * Fbd * Fbd;
    me(4) *= 1. / 12. * KBd * MBd * Fbd * Fbd;

    /*std::cout << "low scale :" << std::endl << std::endl;
    std::cout << "C1_SM :" << C_1_SM << std::endl << std::endl;

    std::cout << "C1 :" << ((*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0)) << std::endl;
    std::cout << "C2 :" << ((*(allcoeff[LO]))(1) + (*(allcoeff[NLO]))(1)) << std::endl;
    std::cout << "C3 :" << ((*(allcoeff[LO]))(2) + (*(allcoeff[NLO]))(2)) << std::endl;
    std::cout << "C4 :" << ((*(allcoeff[LO]))(3) + (*(allcoeff[NLO]))(3)) << std::endl;
    std::cout << "C5 :" << ((*(allcoeff[LO]))(4) + (*(allcoeff[NLO]))(4)) << std::endl << std::endl;*/

    switch (order)
    {
    case FULLNLO:
        return (*(allcoeff[LO]) + *(allcoeff[NLO])) * me /
               (C_1_SM * me(0));
    case LO:
        return ((*(allcoeff[LO])) * me / ((*(allcoeff_SM[LO]))(0) * me(0)));
    default:
        throw std::runtime_error("RBd::RBd(): order not implemented");
    }
}

gslpp::complex AmpDB2::M21_Bd(orders order)
{
    if (mySM.getFlavour().getHDF2().getCoeffBd().getOrder() < getHighest(order))
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed");

    gslpp::vector<gslpp::complex> **allcoeff = mySM.getFlavour().ComputeCoeffBd(
        mySM.getBBd().getMu(),
        mySM.getBBd().getScheme());

    gslpp::vector<double> me(mySM.getBBd().getBpars());
    double MBd = mySM.getMesons(QCD::B_D).getMass();
    double Mb = mySM.Mrun(mySM.getBBd().getMu(),
                          mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
                          mySM.getQuarks(QCD::BOTTOM).getMass(), QCD::BOTTOM, FULLNNLO);
    double Md = mySM.Mrun(mySM.getBBd().getMu(),
                          mySM.getQuarks(QCD::DOWN).getMass_scale(),
                          mySM.getQuarks(QCD::DOWN).getMass(), QCD::DOWN, FULLNNLO);
    double KBd = MBd / (Mb + Md) * MBd / (Mb + Md);
    double Fb = mySM.getMesons(QCD::B_D).getDecayconst();
    me(0) *= 1. / 3. * MBd * Fb * Fb;
    me(1) *= -5. / 24. * KBd * MBd * Fb * Fb;
    me(2) *= 1. / 24. * KBd * MBd * Fb * Fb;
    me(3) *= 1. / 4. * KBd * MBd * Fb * Fb;
    me(4) *= 1. / 12. * KBd * MBd * Fb * Fb;

#if SUSYFIT_DEBUG & 1
    std::cout << "Bd: me(0) = " << me(0) << std::endl;
#endif
#if SUSYFIT_DEBUG & 2
    std::cout << "coefficient Bd: " << (*(allcoeff[LO]) + *(allcoeff[NLO]))(0) << std::endl;
    std::cout << "M: " << me << std::endl;
    std::cout << "mu : " << mySM.getBBd().getMu() << ", mut: " << mySM.getMut() << ", scheme: " << mySM.getBBd().getScheme() << ", B par.: " << mySM.getBBd().getBpars()(0) << std::endl;
    std::cout << "U (mut): " << (mySM.getFlavour().getHDF2().getUDF2().Df2Evol(mySM.getBBd().getMu(), mySM.getMut(), LO)(0, 0) + mySM.getFlavour().getHDF2().getUDF2().Df2Evol(mySM.getBBd().getMu(), mySM.getMut(), NLO)(0, 0)) << std::endl;
#endif

    switch (order)
    {
    case FULLNLO:
        return ((*(allcoeff[LO]) + *(allcoeff[NLO])) * me / HCUT);
    case LO:
        return ((*(allcoeff[LO])) * me / HCUT);
    default:
        throw std::runtime_error("AmpDB2::AmpBd(): order not implemented");
    }
}

gslpp::complex AmpDB2::M21_Bs(orders order)
{
    if (mySM.getFlavour().getHDF2().getCoeffBs().getOrder() < getHighest(order))
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed");

    // Wilson coefficients in same mass scale and scheme as B parameters
    gslpp::vector<gslpp::complex> **allcoeff = mySM.getFlavour().ComputeCoeffBs(
        mySM.getBBs().getMu(),
        mySM.getBBs().getScheme());

    gslpp::vector<double> me(mySM.getBBs().getBpars());
    double MBs = mySM.getMesons(QCD::B_S).getMass();
    double Mb = mySM.Mrun(mySM.getBBs().getMu(),
                          mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
                          mySM.getQuarks(QCD::BOTTOM).getMass(), QCD::BOTTOM, FULLNNLO);
    double Ms = mySM.Mrun(mySM.getBBs().getMu(),
                          mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                          mySM.getQuarks(QCD::STRANGE).getMass(), QCD::STRANGE, FULLNNLO);
    double KBs = MBs / (Mb + Ms) * MBs / (Mb + Ms);
    double Fbs = mySM.getMesons(QCD::B_S).getDecayconst();

    // matrix elements as in hep-ph/9604387v2 with KBs independent terms in 4 and 5
    // differ by factor of 2 to arXiv:1602.03560v2 etc.
    me(0) *= 1. / 3. * MBs * Fbs * Fbs;
    me(1) *= -5. / 24. * KBs * MBs * Fbs * Fbs;
    me(2) *= 1. / 24. * KBs * MBs * Fbs * Fbs;
    me(3) *= 1. / 4. * KBs * MBs * Fbs * Fbs;
    me(4) *= 1. / 12. * KBs * MBs * Fbs * Fbs;
#if SUSYFIT_DEBUG & 1
    std::cout << "Bs: me(0) = " << me(0) << std::endl;
#endif

    switch (order)
    {
    case FULLNLO:
        return ((*(allcoeff[LO]) + *(allcoeff[NLO])) * me / HCUT);
    case LO:
        return ((*(allcoeff[LO])) * me / HCUT);
    default:
        throw std::runtime_error("AmpDB2::AmpBs(): order not implemented");
    }
}

/*******************************************************************************
 *  @f$\Gamma_{21}@f$ in NLO from Ciuchini (hep-ph/0308029v2)                   *
 * ****************************************************************************/

gslpp::complex AmpDB2::Gamma21overM21_tradBasis(orders order, quark q)
{
    computeCKMandMasses(NLO);
   // calculate DB=2 matrix elements for usage of "me" and "delta_1overm_tradBasis(quark)"
    compute_matrixelements(q, FULLNLO);

    // calculate M_21_q / <O_1>
    gslpp::complex M21overme0;
    double MBq;
    switch (q)
    {
    case d:
        M21overme0 = M21_Bd(FULLNLO) * HCUT / me(0); 
        MBq = MB;
        break;
    case s:
        M21overme0 = M21_Bs(FULLNLO) * HCUT / me(0); 
        MBq = MB_s;
        break;
    default:
        throw std::runtime_error("AmpDB2::Gamma21overM21_tradBasis(orders order, quark q): invalid quark index: ");
    }
    
    // calculate DB=1 Wilson coefficients
    computeWilsonCoeffsBuras();

    // calculate DB=2 coefficients for usage of "c(quark)"
    computeF0();
    computeF1();
    computeP();

    // staying in the old basis: hep-ph/0308029v2: eq. 16 divided by M_21
    /*
    computeD(order);
    gslpp::complex Gamma21overM21 = -Gf2 / (24 * M_PI * MBq) / M21overme0 *
                (Mb2_prefactor * (c(q)(0) + c(q)(1) * me(1)/me(0) + c(q)(2) * me(2)/me(0)) +
            Mb2_prefactor_1overm * delta_1overm_tradBasis(q)/me(0));
    */

    // switching to the new basis
    gslpp::complex Gamma21overM21;
    switch (order)
    {
    case FULLNLO:
        computeD(FULLNLO);
        Gamma21overM21 = Mb2_prefactor * (c(q, LO) * meoverme0) + Mb2_prefactor_1overm * (-delta_1overm(q)) / me(0);
        computeD(LO);
        Gamma21overM21 += Mb2_prefactor * (c(q, NLO) * meoverme0);
        Gamma21overM21 *= -Gf2 / (24. * M_PI * MBq) / M21overme0;
        break;
    case LO:
        computeD(LO);
        Gamma21overM21 = -Gf2 / (24. * M_PI * MBq) / M21overme0 *
                         (Mb2_prefactor * (c(q, LO) * meoverme0) + Mb2_prefactor_1overm * (-delta_1overm(q)) / me(0));
        break;
    default:
        throw std::runtime_error("AmpDB2::Gamma21overM21_tradBasis(orders order, quark q): order not implemented");
    }
    return Gamma21overM21;
}

void AmpDB2::computeCKMandMasses(orders order, mass_schemes mass_scheme)
{
    if (order != NLO and order != NNLO)
        throw(std::runtime_error("computeCKMandMasses() order not present"));

    VtbVtd = mySM.getCKM().computelamt_d();
    VtbVts = mySM.getCKM().computelamt_s();
    VtbVtd2 = VtbVtd * VtbVtd;
    VtbVts2 = VtbVts * VtbVts;
    VcbVcd = mySM.getCKM().computelamc_d();
    VcbVcs = mySM.getCKM().computelamc_s();
    VcbVcd2 = VcbVcd * VcbVcd;
    VcbVcs2 = VcbVcs * VcbVcs;
    lambda_c_d = VcbVcd.conjugate();
    lambda_u_d = mySM.getCKM().computelamu_d().conjugate();
    lambda_c_s = VcbVcs.conjugate();
    lambda_u_s = mySM.getCKM().computelamu_s().conjugate();

    // pole mass of bottom quark
    Mb_pole = mySM.Mbar2Mp(mySM.getQuarks(QCD::BOTTOM).getMass(), QCD::BOTTOM);

    // DB=1 matching scales (arxiv: 2205.07907 Results. or Gerlach thesis eq. 7.7) varied by "getMub()" or fixed to 4.2
    mu_1 = mySM.getMub();
    // mass scale mu_b=mu_c
    mu_b = mu_1;
    if (flag_fixmub)
        mu_b = 4.2;
    // mu_1_overm = Mb_pole; // Gerlach thesis
    mu_1_overm = mu_1;

    // DB=2 matching scale mu_2
    mu_2 = mySM.getBBs().getMu();

    // MSbar bottom quark mass Mb(Mb)
    Mb_Mb = mySM.getQuarks(QCD::BOTTOM).getMass();

    // MSbar bottom quark mass Mb(mu_b)
    double Mb_mub = mySM.Mrun(mu_b,
                              mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
                              mySM.getQuarks(QCD::BOTTOM).getMass(), QCD::BOTTOM, FULLNNLO);

    // MSbar charm quark mass Mc(Mc)
    double Mc_Mc = mySM.getQuarks(QCD::CHARM).getMass();

    // MSbar charm quark mass Mc(mu_b)
    double Mc_mub = mySM.Mrun(mu_b,
                              mySM.getQuarks(QCD::CHARM).getMass_scale(),
                              mySM.getQuarks(QCD::CHARM).getMass(), QCD::CHARM, FULLNNLO);

    // strong coupling constant divided by 4*Pi
    as_4pi_mu1 = mySM.Als(mu_1, FULLNNNLO, true, true) / (4. * M_PI);
    as_4pi_mu2 = mySM.Als(mu_2, FULLNNNLO, true, true) / (4. * M_PI);
    as_4pi = mySM.Als(Mb_Mb, FULLNNNLO, true, true) / (4. * M_PI);

    // resummed z always in MSbar
    z = Mc_mub * Mc_mub / (Mb_mub * Mb_mub);
    sqrtz = sqrt(z);

    // PS mass of bottom quark: NNLO evaluation from hep-ph/9804241v2 eq. (25)
    PoletoMS_as2 = 307. / 32. + M_PI2 / 3. + M_PI2 / 9. * log2 - 1. / 6. * zeta3 - (71. / 144. + M_PI2 / 18.) * nl + 4. / 3. * M_PI2 / 8. * sqrtz - z;
    Mb_PS = Mb_Mb * (1. + 16. / 3. * as_4pi * (1. - mu_f / Mb_Mb) + 16. * as_4pi * as_4pi * (PoletoMS_as2 - mu_f / (3. * Mb_Mb) * (a1 - b0 * (2. * log(mu_f / Mb_Mb) - 2.))));

    // adapt "Mb2_prefactor" to the used mass scheme; Mb, Mc, always in MSbar
    // explained in Gerlach thesis chapter 7.0 and arxiv:2205.07907 Results.
    if (order == NNLO)
    {
        Mb = Mb_mub; // if Mb changes independently of mub it has to be added to currentInput_compute_pp_s
        Mc = Mc_mub;
        this->flag_resumz = true;
        switch (mass_scheme)
        {
        case pole:
            Mb2_prefactor = Mb_pole * Mb_pole;
            break;
        case MSbar_partialNNLO:
        case MSbar_partialN3LO:
        case MSbar:
            Mb2_prefactor = Mb_mub * Mb_mub;
            break;
        case PS_partialNNLO:
        case PS_partialN3LO:
        case PS:
            Mb2_prefactor = Mb_PS * Mb_PS;
            break;
        case only1overmb:
            Mb2_prefactor = 0.;
            break;
        default:
            throw(std::runtime_error("mass_scheme not implemented"));
        }
    }
    Mb2_prefactor_1overm = Mb2_prefactor;
    if (mass_scheme == only1overmb) 
        Mb2_prefactor_1overm = Mb_PS * Mb_PS;

    // adapt to MSbar mass scheme for the traditional basis
    if (order == NLO)
    {
        Mb2_prefactor = Mb_mub * Mb_mub;
        Mb = Mb_pole;
        x_1 = mu_1 / Mb;
        x_2 = mu_2 / Mb;
        logx_1 = log(x_1);
        logx_2 = log(x_2);
        double Mc = Mc_mub;
        if (!flag_resumz)
            Mc = Mc_Mc;
        z = Mc * Mc / (Mb_mub * Mb_mub);
    }

    Gf2 = mySM.getGF() * mySM.getGF();
    MB = mySM.getMesons(QCD::B_D).getMass();
    MB_s = mySM.getMesons(QCD::B_S).getMass();

    z2 = z * z;
    logz = log(z);
    log2z = logz * logz;
    oneminusz2 = (1. - z) * (1. - z);
    sqrt1minus4z = sqrt(1. - 4. * z);

    // calculate "z" values for 1/m_b contributions
    // z_1overm = Mc_Mc * Mc_Mc / (Mb_Mb * Mb_Mb); //hep-ph/0612167 eq. 27
    // MSbar charm quark mass Mb(Mb)
    // double Mc_Mb = mySM.Mrun(mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
    //     mySM.getQuarks(QCD::CHARM).getMass_scale(),
    //     mySM.getQuarks(QCD::CHARM).getMass(), FULLNNLO);
    // z_1overm = Mc_Mb * Mc_Mb / (Mb_Mb * Mb_Mb); //arxiv:1910.00970 eq. 11
    // if z for the 1/m_b contributions shall be varied with "mu_b"
    z_1overm = Mc_mub * Mc_mub / (Mb_mub * Mb_mub);

    z_1overm2 = z_1overm * z_1overm;
    oneminusz_1overm2 = (1. - z_1overm) * (1. - z_1overm);
    sqrt1minus4z_1overm = sqrt(1. - 4. * z_1overm);

    z3 = z2 * z;
    z4 = z3 * z;

    log1minusz = log(1. - z);
    log1minus4z = log(1. - 4. * z);

    sigma = (1. - sqrt1minus4z) / (1. + sqrt1minus4z);
    logsigma = log(sigma);
    log2sigma = logsigma * logsigma;

    Dilogz = gslpp_special_functions::dilog(z);
    Dilogsigma = gslpp_special_functions::dilog(sigma);
    Dilogsigma2 = gslpp_special_functions::dilog(sigma * sigma);
    return;
}

void AmpDB2::computeWilsonCoeffsBuras()
{
    // NNLO DB=1 Wilson coefficients
    gslpp::vector<gslpp::complex> **WilsonCoeffsBuras = mySM.getFlavour().ComputeCoeffsgamma_Buras(mu_1);
    for (int i = 0; i < 8; i++)
    {
        if (i == 6)
            i = 7;
        C_Buras_LO[i] = (*(WilsonCoeffsBuras[LO]))(i);
        C_Buras_NLO[i] = (*(WilsonCoeffsBuras[NLO]))(i);
        C_Buras_NNLO[i] = (*(WilsonCoeffsBuras[NNLO]))(i);
    }

    // LO DB=1 Wilson coefficients for 1/mb corrections
    gslpp::vector<gslpp::complex> **WilsonCoeffs_1overm = mySM.getFlavour().ComputeCoeffsgamma_Buras(mu_1_overm);
    C1_LO_1overm = (*(WilsonCoeffs_1overm[LO]))(0);
    C2_LO_1overm = (*(WilsonCoeffs_1overm[LO]))(1);
    K_1 = 3. * C1_LO_1overm * C1_LO_1overm + 2. * C1_LO_1overm * C2_LO_1overm;
    K_2 = C2_LO_1overm * C2_LO_1overm;
    return;
}

// hep-ph/0308029v2: eq.: 41-42: F0 = A
void AmpDB2::computeF0()
{
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
    for (int k = 1; k < 3; k++)
    {
        for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1))
        {
            cacheF0[indexF(qq, k, 2, 1)] = cacheF0[indexF(qq, k, 1, 2)];
        }
    }
    return;
}

// hep-ph/0308029v2: F1 = B
// F0 has to be computed before if flag_resumz is enabled to resum z via hep-ph/0307344 eq.(23)
// see also 2106.05979 eq. (33): prefactor Mb2 in Gamma21 in MSbar scheme
void AmpDB2::computeF1()
{
    double log_muM = 2. * log(mu_b / Mb_pole);
    cacheF1[indexF(cu, 1, 1, 1)] = 109. / 6. - 37. * z + 1.5 * z2 + 52. / 3. * z3 + 2. * oneminusz2 * (5. + z) * logx_2 - 4. * oneminusz2 * (5. + 7. * z) * log1minusz -
                                   2. * z * (10. + 14. * z - 15. * z2) * logz + 8. * (2. - 3. * z + z3) * log1minusz * logz + 16. * (2. - 3. * z + z3) * Dilogz + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(cu, 1, 1, 1) - 8. * z * logz * 1.5 * (-3. + 3. * z2));
    cacheF1[indexF(cu, 2, 1, 1)] = -4. / 3. * (10. - 33. * z + 54. * z2 - 31. * z3) - 8. * oneminusz2 * (4. + 14. * z - 3. * z2) * log1minusz +
                                   8. * z * (2. - 23. * z + 21. * z2 - 3. * z3) * logz -
                                   16. * oneminusz2 * (1. + 2. * z) * (2. * logx_2 - log1minusz * logz - 2. * Dilogz) + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(cu, 2, 1, 1) - 8. * z * logz * 3. * 6. * (z - 1.) * z);
    cacheF1[indexF(cu, 1, 1, 2)] = 1. * ((-502. + 912. * z - 387. * z2 - 23. * z3) / 36. - oneminusz2 * (17. + 4. * z) * logx_1 + 2. / 3. * oneminusz2 * (5. + z) * logx_2 -
                                         oneminusz2 / (12. * z) * (2. + 33. * z + 94. * z2) * log1minusz - z / 12. * (80. + 69. * z - 126. * z2) * logz +
                                         8. / 3. * (2. - 3. * z + z3) * (log1minusz * logz + 2. * Dilogz)) +
                                   flag_resumz * ((32. / 3. + 8. * log_muM) * F0(cu, 1, 1, 2) - 8. * z * logz * 0.5 * (-3. + 3. * z2));
    cacheF1[indexF(cu, 2, 1, 2)] = 1. * ((-130. + 93. * z + 144. * z2 - 107. * z3) / 9. - 2. / 3. * oneminusz2 / z * (1. + 15. * z + 47. * z2 - 12. * z3) * log1minusz +
                                         2. / 3. * z * (8. - 93. * z + 87. * z2 - 12. * z3) * logz -
                                         8. / 3. * oneminusz2 * (1. + 2. * z) * (3. * logx_1 + 4. * logx_2 - 2. * log1minusz * logz - 4. * Dilogz)) +
                                   flag_resumz * ((32. / 3. + 8. * log_muM) * F0(cu, 2, 1, 2) - 8. * z * logz * 6. * (z - 1.) * z);
    cacheF1[indexF(cu, 1, 2, 2)] = -M_PI2 / 3. * (1. - 5. * z + 4. * z2) + (-136. - 159. * z + 738. * z2 - 443. * z3) / 18. - 2 * oneminusz2 * (5. + 4. * z) * logx_1 +
                                   2. / 3. * oneminusz2 * (4. - z) * logx_2 + oneminusz2 / (6. * z) * (7. + 32. * z2 + 3. * z3) * log1minusz -
                                   z / 6. * (62. + 39. * z - 30. * z2 + 3. * z3) * logz + (5. - 3. * z - 18. * z2 + 16. * z3) / 3. * (log1minusz * logz + 2. * Dilogz) + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(cu, 1, 2, 2) - 8. * z * logz * -1.5 * oneminusz2);
    cacheF1[indexF(cu, 2, 2, 2)] = 8. / 3. * M_PI2 * (1. + z - 2. * z2) - 28. / 9. * (5. + 3. * z - 27. * z2 + 19. * z3) - 16. * oneminusz2 * (1. + 2. * z) * logx_1 +
                                   32. / 3. * oneminusz2 * (1. + 2. * z) * logx_2 - 4. / 3. * oneminusz2 / z * (1. - 12. * z - 16. * z2 - 3. * z3) * log1minusz +
                                   4. / 3. * z * (2. - 3. * z + 18. * z2 - 3. * z3) * logz + 8. / 3. * (1. - 3. * z - 6. * z2 + 8. * z3) * (log1minusz * logz + 2. * Dilogz) + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(cu, 2, 2, 2) - 8. * z * logz * -6. * (z - 1.) * z);

    cacheF1[indexF(cc, 1, 1, 1)] = sqrt1minus4z * (109. - 226. * z + 168. * z2) / 6. - (52. - 104. * z - 16. * z2 + 56. * z3) * logsigma +
                                   2. * (5. - 8. * z) * sqrt1minus4z * logx_2 - 12. * sqrt1minus4z * (3. - 2. * z) * log1minus4z + 4. * (13. - 10. * z) * sqrt1minus4z * logz +
                                   16. * (1. - 3. * z + 2. * z2) * (3. * log2sigma + 2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2) + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(cc, 1, 1, 1) - 8. * z * logz * 3. * (6. * z - 3.) / sqrt1minus4z);
    cacheF1[indexF(cc, 2, 1, 1)] = -8. / 3. * sqrt1minus4z * (5. - 23. * z - 42. * z2) - 16. * (4. - 2. * z - 7. * z2 + 14. * z3) * logsigma - 32. * sqrt1minus4z * (1. + 2. * z) * logx_2 -
                                   48. * sqrt1minus4z * (1. + 2. * z) * log1minus4z + 64. * sqrt1minus4z * (1. + 2. * z) * logz +
                                   16. * (1. - 4. * z2) * (3. * log2sigma + 2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2) + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(cc, 2, 1, 1) - 8. * z * logz * 3. * -12. * z / sqrt1minus4z);
    cacheF1[indexF(cc, 1, 1, 2)] = -sqrt1minus4z * (127. - 199. * z - 75. * z2) / 9. + (2. - 259. * z + 662. * z2 - 76. * z3 - 200. * z4) * logsigma / (12. * z) -
                                   (17. - 26. * z) * sqrt1minus4z * logx_1 + 2. / 3. * (5. - 8. * z) * sqrt1minus4z * logx_2 - 4. * sqrt1minus4z * (3. - 2. * z) * log1minus4z -
                                   sqrt1minus4z * (2. - 255. * z + 316. * z2) * logz / (12. * z) +
                                   16. / 3. * (1. - 3. * z + 2. * z2) * (3. * log2sigma + 2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2) + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(cc, 1, 1, 2) - 8. * z * logz * (6. * z - 3.) / sqrt1minus4z);
    cacheF1[indexF(cc, 2, 1, 2)] = -2. * sqrt1minus4z * (68. + 49. * z - 150. * z2) / 9. + 2. / 3. * (1. - 35. * z + 4. * z2 + 76. * z3 - 100. * z4) * logsigma / z +
                                   (16. - 64. * z2) * log2sigma - 8. * sqrt1minus4z * (1. + 2. * z) * logx_1 - 32. / 3. * sqrt1minus4z * (1. + 2. * z) * logx_2 -
                                   16. * sqrt1minus4z * (1. + 2. * z) * log1minus4z - 2. / 3. * sqrt1minus4z * (1. - 33. * z - 76. * z2) * logz / z +
                                   16. / 3. * (1. - 4. * z2) * (2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2) + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(cc, 2, 1, 2) - 8. * z * logz * -12. * z / sqrt1minus4z);
    cacheF1[indexF(cc, 1, 2, 2)] = -M_PI2 / 3. * (1. - 10. * z) - sqrt1minus4z * (115. + 632. * z + 96. * z2) / 18. - (7. + 13. * z - 194. * z2 + 304. * z3 - 64. * z4) * logsigma / (6. * z) -
                                   2. * sqrt1minus4z * (5. - 2. * z) * logx_1 + 4. / 3. * (2. - 5. * z) * sqrt1minus4z * logx_2 - 4. * (1. - 6. * z) * sqrt1minus4z * log1minus4z +
                                   (13. - 54. * z + 8. * z2) * logsigma * log1minus4z / 3. + sqrt1minus4z * (7. + 27. * z - 250. * z2) * logz / (6. * z) +
                                   (7. - 32. * z + 4. * z2) * (log2sigma - logsigma * logz) + 4. / 3. * (5. - 12. * z + 4. * z2) * Dilogsigma + 4. / 3. * (4. - 21. * z + 2. * z2) * Dilogsigma2 + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(cc, 1, 2, 2) - 8. * z * logz * 0.5 * -6. * sqrt1minus4z);
    cacheF1[indexF(cc, 2, 2, 2)] = 8. / 3. * M_PI2 * (1. + 2. * z) - 8. / 9. * sqrt1minus4z * (19. + 53. * z + 24. * z2) + 4. / 3. * (1. + 7. * z + 10. * z2 - 68. * z3 + 32. * z4) * logsigma / z -
                                   8. * (1. + 2. * z) * (1. + 2. * z) * log2sigma - 16. * sqrt1minus4z * (1. + 2. * z) * logx_1 + 32. / 3. * sqrt1minus4z * (1. + 2. * z) * logx_2 +
                                   16. * sqrt1minus4z * (1. + 2. * z) * log1minus4z - 8. / 3. * (1. + 6. * z + 8. * z2) * logsigma * log1minus4z -
                                   4. / 3. * sqrt1minus4z * (1. + 9. * z + 26. * z2) * logz / z + 8. * (1. + 2. * z) * (1. + 2. * z) * logsigma * logz +
                                   32. / 3. * (1. - 4. * z2) * Dilogsigma - 32. / 3. * (1. + 3. * z + 2. * z2) * Dilogsigma2 + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(cc, 2, 2, 2) - 8. * z * logz * 12. * z / sqrt1minus4z);

    cacheF1[indexF(uu, 1, 1, 1)] = 109. / 6. + 10. * logx_2 + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(uu, 1, 1, 1));
    cacheF1[indexF(uu, 2, 1, 1)] = -40. / 3. - 32. * logx_2 + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(uu, 2, 1, 1));
    cacheF1[indexF(uu, 1, 1, 2)] = -127. / 9. + 4. / 12. - 17. * logx_1 + 10. / 3. * logx_2 + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(uu, 1, 1, 2));
    cacheF1[indexF(uu, 2, 1, 2)] = -(136. + 12.) / 9. - 8. * logx_1 - 32. / 3. * logx_2 + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(uu, 2, 1, 2));
    cacheF1[indexF(uu, 1, 2, 2)] = -M_PI2 / 3. - (115. + 42.) / 18. - 10. * logx_1 + 8. / 3. * logx_2 + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(uu, 1, 2, 2));
    cacheF1[indexF(uu, 2, 2, 2)] = 8. * M_PI2 / 3. - 8. / 9. * (19. - 3.) - 16. * logx_1 + 32. / 3. * logx_2 + flag_resumz * ((32. / 3. + 8. * log_muM) * F0(uu, 2, 2, 2));

    for (int k = 1; k < 3; k++)
    {
        for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1))
        {
            cacheF1[indexF(qq, k, 2, 1)] = cacheF1[indexF(qq, k, 1, 2)];
        }
    }
    return;
}

// hep-ph/0308029v2: eq. 50-54
void AmpDB2::computeP()
{
    cacheP[indexP(cu, 1, 2, 2)] = -1. / 27. + 2. / 9. * z - logx_1 / 9. - sqrt1minus4z * (1. + 2. * z) * (2. + 3. * logsigma + 6. * logx_1) / 54. + logz / 18.;
    cacheP[indexP(cu, 2, 2, 2)] = 8. / 27. + 16. / 9. * z + 8. / 9. * logx_1 + 4. / 27. * sqrt1minus4z * (1. + 2. * z) * (2. + 3. * logsigma + 6. * logx_1) - 4. * logz / 9.;
    cacheP[indexP(cc, 1, 2, 2)] = -2. / 27. * sqrt1minus4z * (1. + 8. * z + 12. * z2) - logsigma / 9. + 4. / 3. * z2 * logsigma + 16. / 9. * z3 * logsigma -
                                  sqrt1minus4z * (1. + 2. * z) * (2. * logx_1 - logz) / 9.;
    cacheP[indexP(cc, 2, 2, 2)] = 16. / 27. * sqrt1minus4z * (1. + 8. * z + 12. * z2) + 8. / 9. * logsigma - 32. / 3. * z2 * logsigma - 128. / 9. * z3 * logsigma +
                                  8. / 9. * sqrt1minus4z * (1. + 2. * z) * (2. * logx_1 - logz);
    cacheP[indexP(uu, 1, 2, 2)] = -2. / 27. - 2. / 9. * logx_1;
    cacheP[indexP(uu, 2, 2, 2)] = 16. / 27. + 16. / 9. * logx_1;

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
    cacheP[indexP(cc, 1, 2, 8)] = -1. / 6. * sqrt1minus4z * (1. + 2. * z);
    cacheP[indexP(cc, 2, 2, 8)] = 4. / 3. * sqrt1minus4z * (1. + 2. * z);

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
    cacheP[indexP(uu, 1, 2, 8)] = -1. / 6.;
    cacheP[indexP(uu, 2, 2, 8)] = 4. / 3.;
    return;
}

// hep-ph/0308029v2: 39:
void AmpDB2::computeD(orders order)
{
    if (order == FULLNLO)
    {
        // qq = uu and cc
        for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2))
        {
            for (int k = 1; k <= 2; k++)
            {
                gslpp::complex result = 0.;
                for (int i = 1; i <= 2; i++)
                {
                    for (int j = 1; j <= 2; j++)
                    {
                        result += C_Buras_LO[i - 1] * C_Buras_LO[j - 1] * F(qq, k, i, j) + (C_Buras_NLO[i - 1] * C_Buras_LO[j - 1] + C_Buras_LO[i - 1] * C_Buras_NLO[j - 1]) * F0(qq, k, i, j);
                    }
                }
                result += +as_4pi_mu1 * C_Buras_LO[2 - 1] * C_Buras_LO[2 - 1] * P(qq, k, 2, 2) + 2. * as_4pi_mu1 * C_Buras_LO[2 - 1] * C_Buras_LO[8 - 1] * P(qq, k, 2, 8);
                for (int i = 1; i <= 2; i++)
                {
                    for (int r = 3; r <= 6; r++)
                    {
                        result += 2. * C_Buras_LO[i - 1] * C_Buras_LO[r - 1] * P(qq, k, i, r);
                    }
                }
                cacheD[indexD(qq, k)] = result;
            }
        }
        // qq = cu
        for (int k = 1; k <= 2; k++)
        {
            gslpp::complex result = 0.;
            for (int i = 1; i <= 2; i++)
            {
                for (int j = 1; j <= 2; j++)
                {
                    result += C_Buras_LO[i - 1] * C_Buras_LO[j - 1] * F(cu, k, i, j) + (C_Buras_NLO[i - 1] * C_Buras_LO[j - 1] + C_Buras_LO[i - 1] * C_Buras_NLO[j - 1]) * F0(cu, k, i, j);
                }
            }
            result += +as_4pi_mu1 * C_Buras_LO[2 - 1] * C_Buras_LO[2 - 1] * P(cu, k, 2, 2) + as_4pi_mu1 * C_Buras_LO[2 - 1] * C_Buras_LO[8 - 1] * (P(cc, k, 2, 8) + P(uu, k, 2, 8));
            for (int i = 1; i <= 2; i++)
            {
                for (int r = 3; r <= 6; r++)
                {
                    result += C_Buras_LO[i - 1] * C_Buras_LO[r - 1] * (P(cc, k, i, r) + P(uu, k, i, r));
                }
            }
            cacheD[indexD(cu, k)] = result;
        }
    }
    else if (order == LO)
    {
        // qq = uu and cc
        for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2))
        {
            for (int k = 1; k <= 2; k++)
            {
                gslpp::complex result = 0.;
                for (int i = 1; i <= 2; i++)
                {
                    for (int j = 1; j <= 2; j++)
                    {
                        result += C_Buras_LO[i - 1] * C_Buras_LO[j - 1] * F0(qq, k, i, j);
                    }
                }
                for (int i = 1; i <= 2; i++)
                {
                    for (int r = 3; r <= 6; r++)
                    {
                        result += 2. * C_Buras_LO[i - 1] * C_Buras_LO[r - 1] * P(qq, k, i, r);
                    }
                }
                cacheD[indexD(qq, k)] = result;
            }
        }
        // qq = cu
        for (int k = 1; k <= 2; k++)
        {
            gslpp::complex result = 0.;
            for (int i = 1; i <= 2; i++)
            {
                for (int j = 1; j <= 2; j++)
                {
                    result += C_Buras_LO[i - 1] * C_Buras_LO[j - 1] * F0(cu, k, i, j);
                }
            }
            for (int i = 1; i <= 2; i++)
            {
                for (int r = 3; r <= 6; r++)
                {
                    result += C_Buras_LO[i - 1] * C_Buras_LO[r - 1] * (P(cc, k, i, r) + P(uu, k, i, r));
                }
            }
            cacheD[indexD(cu, k)] = result;
        }
    }
    else
    {
        throw std::runtime_error("AmpDB2::computeD(orders order): order not implemented");
    }
    return;
}

double AmpDB2::F0(quarks qq, int k, int i, int j)
{
    return cacheF0[indexF(qq, k, i, j)];
}

double AmpDB2::F1(quarks qq, int k, int i, int j)
{
    return cacheF1[indexF(qq, k, i, j)];
}

double AmpDB2::F(quarks qq, int k, int i, int j)
{
    return F0(qq, k, i, j) + as_4pi_mu1 * F1(qq, k, i, j);
}

double AmpDB2::P(quarks qq, int k, int i, int j)
{
    return cacheP[indexP(qq, k, i, j)];
}

// indizes for F to save them in an array
int AmpDB2::indexF(quarks qq, int k, int i, int j)
{
    return qq * 8 + (k - 1) * 4 + (i - 1) * 2 + (j - 1);
}
// indizes for P to save them in an array: c := cc and u := uu
int AmpDB2::indexP(quarks qq, int k, int i, int j)
{
    return qq * 28 + (k - 1) * 14 + (i - 1) * 7 + (j - 2);
}

void AmpDB2::compute_matrixelements(quark q, orders order)
{
    double Mq;
    double Mb_mu2 = mySM.Mrun(mu_2,
                              mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
                              mySM.getQuarks(QCD::BOTTOM).getMass(), QCD::BOTTOM, FULLNNLO);
    double MBq2;
    double KBq;
    double FBq2;
    switch (q)
    {
    case d:
        me = mySM.getBBd().getBpars();
        Mq = mySM.Mrun(mu_2,
                       mySM.getQuarks(QCD::DOWN).getMass_scale(),
                       mySM.getQuarks(QCD::DOWN).getMass(), QCD::DOWN, FULLNNLO);
        MBq2 = MB * MB;
        FBq2 = mySM.getMesons(QCD::B_D).getDecayconst() * mySM.getMesons(QCD::B_D).getDecayconst();
        break;
    case s:
        me = mySM.getBBs().getBpars();
        Mq = mySM.Mrun(mu_2,
                       mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                       mySM.getQuarks(QCD::STRANGE).getMass(), QCD::STRANGE, FULLNNLO);
        MBq2 = MB_s * MB_s;
        FBq2 = mySM.getMesons(QCD::B_S).getDecayconst() * mySM.getMesons(QCD::B_S).getDecayconst();
        break;
    default:
        throw std::runtime_error("AmpDB2::compute_matrixelements(quark q): invalid quark index: ");
    }
    // pole mass for mbpow like in: hep-ph/0612167 eq. 28
    // double Mbpow = mySM.Mbar2Mp(mySM.getQuarks(QCD::BOTTOM).getMass());
    // double Mbpow2 = Mbpow * Mbpow;
    // KBq = MBq2 / ((Mbpow + Mq) * (Mbpow + Mq));

    // arXiv:1907.01025v2 equation (4)
    KBq = MBq2 / ((Mb_mu2 + Mq) * (Mb_mu2 + Mq));
    me(0) *= 8. / 3. * MBq2 * FBq2;
    me(1) *= -5. / 3. * KBq * MBq2 * FBq2;
    me(2) *= 1. / 3. * KBq * MBq2 * FBq2;
    me(3) *= 2. * (KBq + 1. / 6.) * MBq2 * FBq2;
    me(4) *= 2. / 3. * (KBq + 3. / 2.) * MBq2 * FBq2;

    // old parameterization from hep-ph/0612167 eq. 28 (or hep-ph/0308029 eq.26)
    /*
    me_R(0) = me(1) + me(2) + 0.5 * me(0);
    me_R(1) = Mq/Mb * me(3);
    me_R(2) = -2./3. * FBq2 * MBq2 * (MBq2 / Mbpow2 - 1.);
    me_R(3) =  7./6. * FBq2 * MBq2 * (MBq2 / Mbpow2 - 1.);
    me_R(4) = 0.5 * (me(2) + 0.5 * me(0) + me(1) - 2. * Mq/Mbpow * me(4) + me_R(2));
    */

    // Gerlach thesis eq.7.5, 7.6
    switch (q)
    {
    case d:
        me_R(2) = mySM.getBBd_subleading().getBpars()(0);
        me_R(3) = mySM.getBBd_subleading().getBpars()(1);
        break;
    case s:
        me_R(2) = mySM.getBBs_subleading().getBpars()(0);
        me_R(3) = mySM.getBBs_subleading().getBpars()(1);
        break;
    }
    me_R(2) *= MBq2 * FBq2;
    me_R(3) *= MBq2 * FBq2;
    me_R(4) = 0.5 * (me(2) + 0.5 * me(0) + me(1) - 2. * Mq / Mb_mu2 * me(4) + me_R(2));

    me_Rtilde(1) = -me_R(2);
    me_Rtilde(2) = me_R(3) + 0.5 * me_R(2);

    // subleading matrix elements that are not independent
    // Gerlach thesis eq. (3.64)
    double L_2 = 2. * log(mu_2 / Mb_mu2);
    double as1_me0 = 0., as1_me2 = 0.;
    if (order == FULLNLO)
    {
        as1_me0 = 4. * L_2 + 26. / 3.;
        as1_me2 = 8. * L_2 + 8.;
    }
    me_R(0) = 0.5 * (1. + as1_me0 * as_4pi_mu2) * me(0) + me(1) + (1. + as1_me2 * as_4pi_mu2) * me(2);
    // arxiv/1907.01025 eq.26
    me_R(1) = Mq / Mb_mu2 * me(3);
    me_Rtilde(0) = Mq / Mb_mu2 * me(4);

    // switch matrix elements to RI scheme
    if (flag_RI)
    {
        me_R(0) = 0.5 * me(0) + me(1) + me(2);
        if (order == FULLNLO)
            me += -as_4pi_mu2 * meMStoRI * me;
    }

    // matrix elements needed for the leading term of Gamma21overM21
    meoverme0(0) = 1.;
    meoverme0(1) = me(1) / me(0);
    meoverme0(2) = me(2) / me(0);
    return;
}

// hep-ph/0308029v2 eq. 18
gslpp::vector<gslpp::complex> AmpDB2::c(quark q, orders order)
{
    gslpp::vector<gslpp::complex> result = gslpp::vector<gslpp::complex>(3, 0.);
    switch (q)
    {
    case d:
        for (int i = 1; i <= 2; i++)
        {
            result.assign(i - 1,
                          VtbVtd2 * D(uu, i) + 2. * VcbVcd * VtbVtd * (D(uu, i) - D(cu, i)) + VcbVcd2 * (D(cc, i) + D(uu, i) - 2. * D(cu, i)));
        }
        break;
    case s:
        for (int i = 1; i <= 2; i++)
        {
            result.assign(i - 1,
                          VtbVts2 * D(uu, i) + 2. * VcbVcs * VtbVts * (D(uu, i) - D(cu, i)) + VcbVcs2 * (D(cc, i) + D(uu, i) - 2. * D(cu, i)));
        }
        break;
    default:
        throw std::runtime_error("AmpDB2::c(quark q, double mu_2): invalid quark index: ");
    }

    // transformation to the new basis (hep-ph/0612167 eq. 18 + 21)
    if (order == LO)
    {
        result.assign(0, result(0) - 0.5 * result(1));
        result.assign(2, -result(1));
        result.assign(1, 0.);
    }
    else if (order == NLO)
    {
        // save NLO contribution from transformation of the tradBasis to RI scheme
        gslpp::vector<gslpp::complex> RI_NLO = gslpp::vector<gslpp::complex>(3, 0.);
        if (flag_RI)
            RI_NLO = as_4pi_mu2 * coeffsMStoRI.transpose() * result;
        double alpha1 = as_4pi_mu2 * 4. / 3. * (12. * log(mu_2 / Mb_pole) + 6.);
        double alpha2 = as_4pi_mu2 * 4. / 3. * (6. * log(mu_2 / Mb_pole) + 13. / 2.);
        result.assign(0, -alpha2 / 2. * result(1));
        result.assign(2, -alpha1 * result(1));
        if (flag_RI)
            result += RI_NLO;
        result.assign(1, 0.);
    }
    else
    {
        throw std::runtime_error("AmpDB2::c(quark q, orders order): order not implemented");
    }
    return result;
}

gslpp::complex AmpDB2::D(quarks qq, int k)
{
    return cacheD[indexD(qq, k)];
}

// indizes for D to save them in arrays
int AmpDB2::indexD(quarks qq, int k)
{
    if (k != 1 and k != 2)
    {
        throw std::runtime_error("AmpDB2::indexD(quarks qq, int k): invalid k");
    }
    return qq * 2 + (k - 1);
}

/*******************************************************************************
 *  1/mb corrections of @f$\Gamma_{21}@f$                                       *
 * ****************************************************************************/

gslpp::complex AmpDB2::delta_1overm_tradBasis(quark q)
{
    // hep-ph/0308029: equation 18
    compute_deltas_1overm_tradBasis(q);
    switch (q)
    {
    case d:
        return VtbVtd2 * deltas_1overm_tradBasis(uu, d) + 2. * VcbVcd * VtbVtd * (deltas_1overm_tradBasis(uu, d) - deltas_1overm_tradBasis(cu, d)) + VcbVcd2 * (deltas_1overm_tradBasis(cc, d) + deltas_1overm_tradBasis(uu, d) - 2. * deltas_1overm_tradBasis(cu, d));
    case s:
        return VtbVts2 * deltas_1overm_tradBasis(uu, s) + 2. * VcbVcs * VtbVts * (deltas_1overm_tradBasis(uu, s) - deltas_1overm_tradBasis(cu, s)) + VcbVcs2 * (deltas_1overm_tradBasis(cc, s) + deltas_1overm_tradBasis(uu, s) - 2. * deltas_1overm_tradBasis(cu, s));
    default:
        throw std::runtime_error("AmpDB2::delta_1overm_tradBasis(quark q): invalid quark index: ");
    }
}

void AmpDB2::compute_deltas_1overm_tradBasis(quark q)
{
    // hep-ph/0308029v2 eq.20
    cache_deltas_1overm_tradBasis[index_deltas(cc, q)] = sqrt1minus4z_1overm * ((1 + 2. * z_1overm) * (K_2 * (me_R(2) + 2. * me_R(4)) - 2. * K_1 * (me_R(1) + me_R(2))) - 12. * z_1overm2 / (1. - 4. * z_1overm) * (K_1 * (me_R(2) + 2. * me_R(3)) + 2. * K_2 * me_R(3)));
    cache_deltas_1overm_tradBasis[index_deltas(cu, q)] = oneminusz_1overm2 * ((1. + 2. * z_1overm) * (K_2 * (me_R(2) + 2. * me_R(4)) - 2. * K_1 * (me_R(1) + me_R(2))) - 6. * z_1overm2 / (1. - z_1overm) * (K_1 * (me_R(2) + 2. * me_R(3)) + 2. * K_2 * me_R(3)));
    cache_deltas_1overm_tradBasis[index_deltas(uu, q)] = K_2 * (me_R(2) + 2. * me_R(4)) - 2. * K_1 * (me_R(1) + me_R(2));
    return;
}

gslpp::complex AmpDB2::deltas_1overm_tradBasis(quarks qq, quark q)
{
    return cache_deltas_1overm_tradBasis[index_deltas(qq, q)];
}

gslpp::complex AmpDB2::deltas_1overm(quarks qq, quark q)
{
    return cache_deltas_1overm[index_deltas(qq, q)];
}

// indizes for deltas_1overm to save them in an array
int AmpDB2::index_deltas(quarks qq, quark q)
{
    return qq * 2 + q;
}

gslpp::complex AmpDB2::delta_1overm(quark q)
{
    // hep-ph/0612167 equation (10)
    gslpp::complex lambda_c, lambda_u;
    switch (q)
    {
    case d:
        lambda_c = lambda_c_d;
        lambda_u = lambda_u_d;
        break;
    case s:
        lambda_c = lambda_c_s;
        lambda_u = lambda_u_s;
        break;
    default:
        throw std::runtime_error("AmpDB2::delta_1overm(quark q): invalid quark index: ");
    }
    compute_deltas_1overm(q);
    return -(lambda_c * lambda_c * deltas_1overm(cc, q) + 2. * lambda_c * lambda_u * deltas_1overm(cu, q) + lambda_u * lambda_u * deltas_1overm(uu, q)).conjugate();
}

void AmpDB2::compute_deltas_1overm(quark q)
{
    // hep-ph/0612167 eq.24
    compute_g();
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1))
    {
        cache_deltas_1overm[index_deltas(qq, q)] = 0.;
        for (int i = 0; i <= 3; i++)
        {
            cache_deltas_1overm[index_deltas(qq, q)] += g(qq, i) * me_R(i);
        }
        for (int i = 0; i <= 2; i++)
        {
            cache_deltas_1overm[index_deltas(qq, q)] += gtilde(qq, i) * me_Rtilde(i);
        }
    }
    return;
}

void AmpDB2::compute_g()
{
    // hep-ph/0612167
    // equation (25)
    double cache_z = z_1overm;
    double cache_z2 = z_1overm2;
    double cache_oneminusz_1overm2 = oneminusz_1overm2;
    double cache_sqrt1minus4z = sqrt1minus4z_1overm;
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2))
    {
        cacheg[indexg(qq, 0)] = sqrt1minus4z_1overm * (1. + 2. * z_1overm) * K_1;
        cacheg[indexg(qq, 1)] = -2. * sqrt1minus4z_1overm * (1. + 2. * z_1overm) * K_1;
        cacheg[indexg(qq, 2)] = -2. * (1. - 2. * z_1overm - 2. * z_1overm2) / sqrt1minus4z_1overm * K_1;
        cacheg[indexg(qq, 3)] = -24. * z_1overm2 / sqrt1minus4z_1overm * K_1;
        cachegtilde[indexg(qq, 0)] = -2. * sqrt1minus4z_1overm * (1. + 2. * z_1overm) * K_2;
        cachegtilde[indexg(qq, 1)] = -2. * (1. - 2. * z_1overm - 2. * z_1overm2) / sqrt1minus4z_1overm * K_2;
        cachegtilde[indexg(qq, 2)] = -24. * z_1overm2 / sqrt1minus4z_1overm * K_2;
        z_1overm = 0.;
        z_1overm2 = 0.;
        oneminusz_1overm2 = 1.;
        sqrt1minus4z_1overm = 1.;
    }
    // equation (26)
    z_1overm = cache_z;
    z_1overm2 = cache_z2;
    oneminusz_1overm2 = cache_oneminusz_1overm2;
    sqrt1minus4z_1overm = cache_sqrt1minus4z;
    cacheg[indexg(cu, 0)] = oneminusz_1overm2 * (1. + 2. * z_1overm) * K_1;
    cacheg[indexg(cu, 1)] = -2. * oneminusz_1overm2 * (1. + 2. * z_1overm) * K_1;
    cacheg[indexg(cu, 2)] = -2. * (1. - z_1overm) * (1. + z_1overm + z_1overm2) * K_1;
    ;
    cacheg[indexg(cu, 3)] = -12. * (1. - z_1overm) * z_1overm2 * K_1;
    cachegtilde[indexg(cu, 0)] = -2. * oneminusz_1overm2 * (1. + 2. * z_1overm) * K_2;
    cachegtilde[indexg(cu, 1)] = -2. * (1. - z_1overm) * (1. + z_1overm + z_1overm2) * K_2;
    cachegtilde[indexg(cu, 2)] = -12. * (1. - z_1overm) * z_1overm2 * K_2;
    return;
}

gslpp::complex AmpDB2::g(quarks qq, int i)
{
    return cacheg[indexg(qq, i)];
}

gslpp::complex AmpDB2::gtilde(quarks qq, int i)
{
    return cachegtilde[indexg(qq, i)];
}

int AmpDB2::indexg(quarks qq, int i)
{
    return qq * 4 + i;
}

/*******************************************************************************
 *  @f$\Gamma_{21}@f$ in NNLO from Marvin Gerlach (2205.07907 and thesis)       *
 * ****************************************************************************/

gslpp::complex AmpDB2::Gamma21overM21(orders order, mass_schemes mass_scheme, int BMeson)
{
    computeCKMandMasses(NNLO, mass_scheme);

    // calculate M_21_q / <O_1>
    gslpp::complex M21overme0;
    double MBq;
    quark q;
    switch (BMeson)
    {
    case 0:
        q = d;
        // calculate DB=2 matrix elements for usage of "me" and "delta_1overm_(quark)"
        compute_matrixelements(q, FULLNLO);
        M21overme0 = M21_Bd(FULLNLO) * HCUT / me(0);
        MBq = MB;
        break;
    case 1:
        q = s;
        // calculate DB=2 matrix elements for usage of "me" and "delta_1overm_(quark)"
        compute_matrixelements(q, FULLNLO);
        M21overme0 = M21_Bs(FULLNLO) * HCUT / me(0);
        MBq = MB_s;
        break;
    default:
        throw std::runtime_error("AmpDB2::Gamma21overM21(orders order, int BMeson): invalid B meson index: ");
    }

    // calculate DB=1 Wilson coefficients
    computeWilsonCoeffsMisiak();

    // calculate DB=2 coefficients in pole scheme for usage of "c_H()"
    compute_pp_s();

    // switch to another scheme if needed
    if (mass_scheme == MSbar)
    {
        poletoMSbar_pp_s();
    }
    else if (mass_scheme == PS)
    {
        poletoPS_pp_s();
    }
    else if (mass_scheme == pole or mass_scheme == only1overmb)
    {
    }
    else if (mass_scheme == MSbar_partialNNLO)
    {
        compute_partialNNLO();
        poletoMSbar_pp_s();
    }
    else if (mass_scheme == PS_partialNNLO)
    {
        compute_partialNNLO();
        poletoPS_pp_s();
    }
    else if (mass_scheme == MSbar_partialN3LO)
    {
        poletoMSbar_pp_s_partialN3LO();
    }
    else if (mass_scheme == PS_partialN3LO)
    {
        poletoPS_pp_s_partialN3LO();
    }
    else
    {
        std::cerr << "WARNING: mass_scheme might no be implemented.\n";
    }

    /* partial contribution without 1overm
    gslpp::vector<gslpp::complex> Gamma21overM21_partial = Mb2_prefactor * (c_H_partial(q, 0) + c_H_partial(q, 1) * me(1)/me(0) + c_H_partial(q, 2) * me(2)/me(0)).conjugate()
            * Gf2 / (24 * M_PI * MBq) / M21overme0; */

    // Gerlach thesis eq. 6.1 divided by M_21
    gslpp::complex Gamma21overM21 = 0.;
    if (flag_RI)
    {
        Gamma21overM21 = Mb2_prefactor * (c_H(q, LO) * meoverme0).conjugate();
        // compute_matrixelements(q, LO); //minimal difference
        Gamma21overM21 += Mb2_prefactor * (c_H(q, NLO) * meoverme0).conjugate();
    }
    else
    {
        Gamma21overM21 = Mb2_prefactor * (c_H(q, order) * meoverme0).conjugate();
    }
    computeWilsonCoeffsBuras(); /*calculate DB=1 Wilson coefficients in the basis for "delta_1overm" */
    Gamma21overM21 += Mb2_prefactor_1overm * delta_1overm(q) / me(0);
    Gamma21overM21 *= Gf2 / (24. * M_PI * MBq) / M21overme0;
    return Gamma21overM21;
}

void AmpDB2::computeWilsonCoeffsMisiak()
{
    // NLO DB=1 Wilson coefficients C_i, i=1-6,8
    gslpp::vector<gslpp::complex> **WilsonCoeffsMisiak = mySM.getFlavour().ComputeCoeffsgamma(mu_1);
    for (int i = 0; i < 8; i++)
    {
        if (i == 6)
            i = 7;
        C_Misiak_LO[i] = (*(WilsonCoeffsMisiak[LO]))(i);
        C_Misiak_NLO[i] = (*(WilsonCoeffsMisiak[NLO]))(i);
        C_Misiak_NNLO[i] = (*(WilsonCoeffsMisiak[NNLO]))(i);
    }
}

// Gerlach thesis eq. (6.2)
gslpp::vector<gslpp::complex> AmpDB2::c_H(quark q, orders order)
{
    gslpp::vector<gslpp::complex> result = gslpp::vector<gslpp::complex>(3, 0.);
    gslpp::complex lambda_c, lambda_u;
    switch (q)
    {
    case d:
        lambda_c = lambda_c_d;
        lambda_u = lambda_u_d;
        break;
    case s:
        lambda_c = lambda_c_s;
        lambda_u = lambda_u_s;
        break;
    default:
        throw std::runtime_error("AmpDB2::c_H(quark q): invalid quark index: ");
    }
    if (flag_RI and order == NLO)
    {
        // LO transformed to the traditional basis
        result.assign(0, -lambda_c * lambda_c * H(cc, LO) - 2. * lambda_c * lambda_u * H(cu, LO) - lambda_u * lambda_u * H(uu, LO));
        result.assign(2, -lambda_c * lambda_c * H_s(cc, LO) - 2. * lambda_c * lambda_u * H_s(cu, LO) - lambda_u * lambda_u * H_s(uu, LO));
        result.assign(0, result(0) - 0.5 * result(2));
        result.assign(1, -result(2));
        result.assign(2, 0.);
        // only NLO change to RI
        result = as_4pi_mu2 * coeffsMStoRI.transpose() * result;
        // return to tradBasis + NLO
        result.assign(0, result(0) - 0.5 * result(1) +
                             -lambda_c * lambda_c * H(cc, NLO) - 2. * lambda_c * lambda_u * H(cu, NLO) - lambda_u * lambda_u * H(uu, NLO));
        result.assign(2, -result(1) +
                             -lambda_c * lambda_c * H_s(cc, NLO) - 2. * lambda_c * lambda_u * H_s(cu, NLO) - lambda_u * lambda_u * H_s(uu, NLO));
        result.assign(1, 0.);
        return result;
    }
    if (flag_RI and order != LO)
        throw std::runtime_error("AmpDB2::c_H(quark q, orders order) order for RI not implemented");

    result.assign(0, -lambda_c * lambda_c * H(cc, order) - 2. * lambda_c * lambda_u * H(cu, order) - lambda_u * lambda_u * H(uu, order));
    result.assign(2, -lambda_c * lambda_c * H_s(cc, order) - 2. * lambda_c * lambda_u * H_s(cu, order) - lambda_u * lambda_u * H_s(uu, order));
    return result;
}

// Gerlach thesis eq. (6.4)
gslpp::complex AmpDB2::H(quarks qq, orders order)
{
    if (order == LO)
        return H_partial(qq, 1, 8, 1, 8, 0);
    if (order == NLO)
        return H_partial(qq, 1, 8, 1, 8, 1);
    if (order == NNLO)
        return H_partial(qq, 1, 8, 1, 8, 2);
    if (order == FULLNLO)
        return H_partial(qq, 1, 8, 1, 8, 0) + H_partial(qq, 1, 8, 1, 8, 1);
    if (order == FULLNNLO)
        return H_partial(qq, 1, 8, 1, 8, 0) + H_partial(qq, 1, 8, 1, 8, 1) + H_partial(qq, 1, 8, 1, 8, 2);
    if (order == FULLNNNLO)
        return H_partial(qq, 1, 8, 1, 8, 0) + H_partial(qq, 1, 8, 1, 8, 1) + H_partial(qq, 1, 8, 1, 8, 2) + H_partial(qq, 1, 8, 1, 8, 3);
    throw std::runtime_error("AmpDB2::H(quarks qq, orders order) order not implemented");
}
gslpp::complex AmpDB2::H_s(quarks qq, orders order)
{
    if (order == LO)
        return H_s_partial(qq, 1, 8, 1, 8, 0);
    if (order == NLO)
        return H_s_partial(qq, 1, 8, 1, 8, 1);
    if (order == NNLO)
        return H_s_partial(qq, 1, 8, 1, 8, 2);
    if (order == FULLNLO)
        return H_s_partial(qq, 1, 8, 1, 8, 0) + H_s_partial(qq, 1, 8, 1, 8, 1);
    if (order == FULLNNLO)
        return H_s_partial(qq, 1, 8, 1, 8, 0) + H_s_partial(qq, 1, 8, 1, 8, 1) + H_s_partial(qq, 1, 8, 1, 8, 2);
    if (order == FULLNNNLO)
        return H_s_partial(qq, 1, 8, 1, 8, 0) + H_s_partial(qq, 1, 8, 1, 8, 1) + H_s_partial(qq, 1, 8, 1, 8, 2) + H_s_partial(qq, 1, 8, 1, 8, 3);
    throw std::runtime_error("AmpDB2::H_s(quarks qq, orders order) order not implemented");
}

gslpp::vector<gslpp::complex> AmpDB2::c_H_partial(quark q, int i)
{
    gslpp::complex lambda_c, lambda_u;
    switch (q)
    {
    case d:
        lambda_c = lambda_c_d;
        lambda_u = lambda_u_d;
        break;
    case s:
        lambda_c = lambda_c_s;
        lambda_u = lambda_u_s;
        break;
    default:
        throw std::runtime_error("AmpDB2::c_H_partial(quark q, int i): invalid quark index: ");
    }
    gslpp::vector<gslpp::complex> zeros(12, 0.);
    switch (i)
    {
    case 0:
        return -lambda_c * lambda_c * H_allpartial(cc) - 2. * lambda_c * lambda_u * H_allpartial(cu) - lambda_u * lambda_u * H_allpartial(uu);
    case 1:
        return zeros;
    case 2:
        return -lambda_c * lambda_c * H_s_allpartial(cc) - 2. * lambda_c * lambda_u * H_s_allpartial(cu) - lambda_u * lambda_u * H_s_allpartial(uu);
    default:
        throw std::runtime_error("AmpDB2::c_H_partial(quark q, int i): invalid index i");
    }
}

gslpp::vector<gslpp::complex> AmpDB2::H_allpartial(quarks qq)
{
    gslpp::vector<gslpp::complex> result(12, 0.);
    result.assign(0, H_partial(qq, 1, 2, 1, 2, 0));
    result.assign(1, H_partial(qq, 1, 2, 1, 2, 1));
    result.assign(2, H_partial(qq, 1, 2, 1, 2, 2));
    result.assign(3, H_partial(qq, 1, 2, 3, 6, 0));
    result.assign(4, H_partial(qq, 1, 2, 3, 6, 1));
    result.assign(5, H_partial(qq, 3, 6, 3, 6, 0));
    result.assign(6, H_partial(qq, 3, 6, 3, 6, 1));
    result.assign(7, H_partial(qq, 1, 2, 8, 8, 1));
    result.assign(8, H_partial(qq, 1, 2, 8, 8, 2));
    result.assign(9, H_partial(qq, 3, 6, 8, 8, 1));
    result.assign(10, H_partial(qq, 3, 6, 8, 8, 2));
    result.assign(11, H_partial(qq, 8, 8, 8, 8, 2));
    return result;
}

gslpp::vector<gslpp::complex> AmpDB2::H_s_allpartial(quarks qq)
{
    gslpp::vector<gslpp::complex> result(12, 0.);
    result.assign(0, H_s_partial(qq, 1, 2, 1, 2, 0));
    result.assign(1, H_s_partial(qq, 1, 2, 1, 2, 1));
    result.assign(2, H_s_partial(qq, 1, 2, 1, 2, 2));
    result.assign(3, H_s_partial(qq, 1, 2, 3, 6, 0));
    result.assign(4, H_s_partial(qq, 1, 2, 3, 6, 1));
    result.assign(5, H_s_partial(qq, 3, 6, 3, 6, 0));
    result.assign(6, H_s_partial(qq, 3, 6, 3, 6, 1));
    result.assign(7, H_s_partial(qq, 1, 2, 8, 8, 1));
    result.assign(8, H_s_partial(qq, 1, 2, 8, 8, 2));
    result.assign(9, H_s_partial(qq, 3, 6, 8, 8, 1));
    result.assign(10, H_s_partial(qq, 3, 6, 8, 8, 2));
    result.assign(11, H_s_partial(qq, 8, 8, 8, 8, 2));
    return result;
}

gslpp::complex AmpDB2::H_partial(quarks qq, int i_start, int i_end, int j_start, int j_end, int n)
{
    gslpp::complex result = 0.;
    for (int i = i_start; i <= i_end; i++)
    {
        if (i == 7)
            i++;
        for (int j = j_start; j <= j_end; j++)
        {
            if (j == 7)
                j++;
            if (j < i)
                continue;
            if (n == 0)
            {
                result += C_Misiak_LO[i - 1] * C_Misiak_LO[j - 1] * p(qq, i, j, 0);
            }
            else if (n == 1)
            {
                if (j == 1 or j == 2 or j == 8)
                {
                    result += as_4pi_mu1 * C_Misiak_LO[i - 1] * C_Misiak_LO[j - 1] * p(qq, i, j, 1) + (C_Misiak_NLO[i - 1] * C_Misiak_LO[j - 1] + C_Misiak_LO[i - 1] * C_Misiak_NLO[j - 1]) * p(qq, i, j, 0);
                }
                else if (3 <= j and j <= 6)
                {
                    result += as_4pi_mu1 * C_Misiak_LO[i - 1] * C_Misiak_LO[j - 1] * p(qq, i, j, 1) + (C_Misiak_NLO[i - 1] * C_Misiak_LO[j - 1] + C_Misiak_LO[i - 1] * C_Misiak_NLO[j - 1]) * p(qq, i, j, 0, flag_LOz);
                }
            }
            else if (n == 2)
            {
                if (j == 1 or j == 2 or j == 8)
                {
                    result += as_4pi_mu1 * as_4pi_mu1 * C_Misiak_LO[i - 1] * C_Misiak_LO[j - 1] * p(qq, i, j, 2) + as_4pi_mu1 * (C_Misiak_NLO[i - 1] * C_Misiak_LO[j - 1] + C_Misiak_LO[i - 1] * C_Misiak_NLO[j - 1]) * p(qq, i, j, 1, flag_LOz) + (C_Misiak_NNLO[i - 1] * C_Misiak_LO[j - 1] + C_Misiak_NLO[i - 1] * C_Misiak_NLO[j - 1] + C_Misiak_LO[i - 1] * C_Misiak_NNLO[j - 1]) * p(qq, i, j, 0, flag_LOz);
                }
            }
            else if (n == 3)
            {
                result += as_4pi_mu1 * as_4pi_mu1 * as_4pi_mu1 * C_Misiak_LO[i - 1] * C_Misiak_LO[j - 1] * p(qq, i, j, 3) + as_4pi_mu1 * as_4pi_mu1 * (C_Misiak_NLO[i - 1] * C_Misiak_LO[j - 1] + C_Misiak_LO[i - 1] * C_Misiak_NLO[j - 1]) * p(qq, i, j, 2) + as_4pi_mu1 * (C_Misiak_NNLO[i - 1] * C_Misiak_LO[j - 1] + C_Misiak_NLO[i - 1] * C_Misiak_NLO[j - 1] + C_Misiak_LO[i - 1] * C_Misiak_NNLO[j - 1]) * p(qq, i, j, 1) + (C_Misiak_NNLO[i - 1] * C_Misiak_NLO[j - 1] + C_Misiak_NLO[i - 1] * C_Misiak_NNLO[j - 1]) * p(qq, i, j, 0);
            }
            else
            {
                throw(std::runtime_error("AmpDB2::H_partial order not implemented"));
            }
        }
    }
    return result;
}

gslpp::complex AmpDB2::H_s_partial(quarks qq, int i_start, int i_end, int j_start, int j_end, int n)
{
    gslpp::complex result = 0.;
    for (int i = i_start; i <= i_end; i++)
    {
        if (i == 7)
            i++;
        for (int j = j_start; j <= j_end; j++)
        {
            if (j == 7)
                j++;
            if (j < i)
                continue;
            if (n == 0)
            {
                result += C_Misiak_LO[i - 1] * C_Misiak_LO[j - 1] * p_s(qq, i, j, 0);
            }
            else if (n == 1)
            {
                if (j == 1 or j == 2 or j == 8)
                {
                    result += as_4pi_mu1 * C_Misiak_LO[i - 1] * C_Misiak_LO[j - 1] * p_s(qq, i, j, 1) + (C_Misiak_NLO[i - 1] * C_Misiak_LO[j - 1] + C_Misiak_LO[i - 1] * C_Misiak_NLO[j - 1]) * p_s(qq, i, j, 0);
                }
                else if (3 <= j and j <= 6)
                {
                    result += as_4pi_mu1 * C_Misiak_LO[i - 1] * C_Misiak_LO[j - 1] * p_s(qq, i, j, 1) + (C_Misiak_NLO[i - 1] * C_Misiak_LO[j - 1] + C_Misiak_LO[i - 1] * C_Misiak_NLO[j - 1]) * p_s(qq, i, j, 0, flag_LOz);
                }
            }
            else if (n == 2)
            {
                if (j == 1 or j == 2 or j == 8)
                {
                    result += as_4pi_mu1 * as_4pi_mu1 * C_Misiak_LO[i - 1] * C_Misiak_LO[j - 1] * p_s(qq, i, j, 2) + as_4pi_mu1 * (C_Misiak_NLO[i - 1] * C_Misiak_LO[j - 1] + C_Misiak_LO[i - 1] * C_Misiak_NLO[j - 1]) * p_s(qq, i, j, 1, flag_LOz) + (C_Misiak_NNLO[i - 1] * C_Misiak_LO[j - 1] + C_Misiak_NLO[i - 1] * C_Misiak_NLO[j - 1] + C_Misiak_LO[i - 1] * C_Misiak_NNLO[j - 1]) * p_s(qq, i, j, 0, flag_LOz);
                }
            }
            else if (n == 3)
            {
                result += as_4pi_mu1 * as_4pi_mu1 * as_4pi_mu1 * C_Misiak_LO[i - 1] * C_Misiak_LO[j - 1] * p_s(qq, i, j, 3) + as_4pi_mu1 * as_4pi_mu1 * (C_Misiak_NLO[i - 1] * C_Misiak_LO[j - 1] + C_Misiak_LO[i - 1] * C_Misiak_NLO[j - 1]) * p_s(qq, i, j, 2) + as_4pi_mu1 * (C_Misiak_NNLO[i - 1] * C_Misiak_LO[j - 1] + C_Misiak_NLO[i - 1] * C_Misiak_NLO[j - 1] + C_Misiak_LO[i - 1] * C_Misiak_NNLO[j - 1]) * p_s(qq, i, j, 1) + (C_Misiak_NNLO[i - 1] * C_Misiak_NLO[j - 1] + C_Misiak_NLO[i - 1] * C_Misiak_NNLO[j - 1]) * p_s(qq, i, j, 0);
            }
            else
            {
                throw(std::runtime_error("AmpDB2::H_s_partial order not implemented"));
            }
        }
    }
    return result;
}

double AmpDB2::p(quarks qq, int i, int j, int n, bool flag_LOz)
{
    if (n < 0 or n >= 768)
        throw(std::runtime_error("AmpDB2::p(quarks qq, int i, int j, int n, bool flag_LOz) out of index"));
    if (flag_LOz)
    {
        if (n >= 576)
            throw(std::runtime_error("AmpDB2::p(quarks qq, int i, int j, int n, bool flag_LOz) out of index"));
        return cache_p_LO[index_p(qq, i, j, n)];
    }
    return cache_p[index_p(qq, i, j, n)];
}
double AmpDB2::p_s(quarks qq, int i, int j, int n, bool flag_LOz)
{
    if (n < 0 or n >= 768)
        throw(std::runtime_error("AmpDB2::p_s(quarks qq, int i, int j, int n, bool flag_LOz) out of index"));
    if (flag_LOz)
    {
        if (n >= 576)
            throw(std::runtime_error("AmpDB2::p_s(quarks qq, int i, int j, int n, bool flag_LOz) out of index"));
        return cache_ps_LO[index_p(qq, i, j, n)];
    }
    return cache_ps[index_p(qq, i, j, n)];
}

int AmpDB2::index_p(quarks qq, int i, int j, int n)
{
    return n * 192 + qq * 64 + (i - 1) * 8 + (j - 1);
}

void AmpDB2::compute_pp_s()
{
    // input didn't change nothing to compute
    // double currentInput_compute_pp_s[4] = {z, mu_1, mu_2, mu_b}; // Removed unused variable
    // if (lastInput_compute_pp_s == currentInput_compute_pp_s) return;

    double cache_logz = logz;
    if (!flag_resumz)
    {
        throw std::runtime_error("AmpDB2::compute_pp_s(): only implemented for resummed z");
    }

    // remember value of z after setting it to 0 for calculation of uu coefficients
    const double cache_z = z;
    const double cache_sqrt1minus4z = sqrt1minus4z;
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2))
    {
        // Gerlach thesis eq. (6.6)
        cache_p[index_p(qq, 1, 1, 0)] = (23. * sqrt1minus4z) / 72. - (43. * sqrt1minus4z * z) / 36.;
        cache_p[index_p(qq, 1, 2, 0)] = sqrt1minus4z / 6. - (5. * sqrt1minus4z * z) / 3.;
        cache_p[index_p(qq, 2, 2, 0)] = sqrt1minus4z - sqrt1minus4z * z;
        cache_ps[index_p(qq, 1, 1, 0)] = (-5. * sqrt1minus4z) / 9. - (10. * sqrt1minus4z * z) / 9.;
        cache_ps[index_p(qq, 1, 2, 0)] = (-4. * sqrt1minus4z) / 3. - (8. * sqrt1minus4z * z) / 3.;
        cache_ps[index_p(qq, 2, 2, 0)] = sqrt1minus4z + 2. * sqrt1minus4z * z;

        cache_p[index_p(qq, 1, 3, 0)] = 4. / 3.;
        cache_p[index_p(qq, 1, 4, 0)] = -5. / 36.;
        cache_p[index_p(qq, 1, 5, 0)] = (64. * sqrt1minus4z) / 3. - (160. * sqrt1minus4z * z) / 3.;
        cache_p[index_p(qq, 1, 6, 0)] = (-20. * sqrt1minus4z) / 9. - (4. * sqrt1minus4z * z) / 9.;
        cache_p[index_p(qq, 2, 3, 0)] = 1.;
        cache_p[index_p(qq, 2, 4, 0)] = 5. / 6.;
        cache_p[index_p(qq, 2, 5, 0)] = 16. * sqrt1minus4z - 40. * sqrt1minus4z * z;
        cache_p[index_p(qq, 2, 6, 0)] = (40. * sqrt1minus4z) / 3. + (8. * sqrt1minus4z * z) / 3.;

        cache_ps[index_p(qq, 1, 3, 0)] = -8. / 3.;
        cache_ps[index_p(qq, 1, 4, 0)] = -2. / 9.;
        cache_ps[index_p(qq, 1, 5, 0)] = -128. / 3.;
        cache_ps[index_p(qq, 1, 6, 0)] = -32. / 9.;
        cache_ps[index_p(qq, 2, 3, 0)] = -2.;
        cache_ps[index_p(qq, 2, 4, 0)] = 4. / 3.;
        cache_ps[index_p(qq, 2, 5, 0)] = -32.;
        cache_ps[index_p(qq, 2, 6, 0)] = 64. / 3.;

        // Gerlach thesis eq. (6.9)
        z = 0.;
        sqrt1minus4z = 1.;
    }
    z = cache_z;
    sqrt1minus4z = cache_sqrt1minus4z;

    // Gerlach thesis eq. (6.8)
    double n_l = 3.;
    double n_v = 1.;
    cache_p[index_p(cc, 3, 3, 0)] = 3. * (n_l + n_v) + 2.;
    cache_p[index_p(cc, 3, 4, 0)] = 7. / 3.;
    cache_p[index_p(cc, 3, 5, 0)] = 60. * (n_l + n_v) + 64.;
    cache_p[index_p(cc, 3, 6, 0)] = 112. / 3.;
    cache_p[index_p(cc, 4, 4, 0)] = 5. * (n_l + n_v) / 12. + 13. / 72.;
    cache_p[index_p(cc, 4, 5, 0)] = 112. / 3.;
    cache_p[index_p(cc, 4, 6, 0)] = 25. / 3. * (n_l + n_v) + 52. / 9.;
    cache_p[index_p(cc, 5, 5, 0)] = 512. + 408. * n_l + 408. * n_v * sqrt1minus4z - 480. * n_v * sqrt1minus4z * z;
    cache_p[index_p(cc, 5, 6, 0)] = 1792. / 3.;
    cache_p[index_p(cc, 6, 6, 0)] = 416. / 9. + (170. * n_l) / 3. + (170. * n_v * sqrt1minus4z) / 3. + (124. * n_v * sqrt1minus4z * z) / 3.;

    cache_ps[index_p(cc, 3, 3, 0)] = -6. * (n_l + n_v) - 1.;
    cache_ps[index_p(cc, 3, 4, 0)] = -8. / 3.;
    cache_ps[index_p(cc, 3, 5, 0)] = -120. * (n_l + n_v) - 32.;
    cache_ps[index_p(cc, 3, 6, 0)] = -128. / 3.;
    cache_ps[index_p(cc, 4, 4, 0)] = 2. / 3. * (n_l + n_v) - 7. / 9.;
    cache_ps[index_p(cc, 4, 5, 0)] = -128. / 3.;
    cache_ps[index_p(cc, 4, 6, 0)] = 40. / 3. * (n_l + n_v) - 224. / 9.;
    cache_ps[index_p(cc, 5, 5, 0)] = -816. * (n_l + n_v) - 256.;
    cache_ps[index_p(cc, 5, 6, 0)] = -2048. / 3.;
    cache_ps[index_p(cc, 6, 6, 0)] = 272. / 3. * (n_l + n_v) - 1792. / 9.;

    // Gerlach thesis eq. (6.9)
    for (int i = 3; i <= 6; i++)
    {
        for (int j = i; j <= 6; j++)
        {
            cache_p[index_p(uu, i, j, 0)] = cache_p[index_p(cc, i, j, 0)];
            cache_ps[index_p(uu, i, j, 0)] = cache_ps[index_p(cc, i, j, 0)];
        }
    }

    // Gerlach thesis eq. (6.10)
    for (int i = 1; i <= 6; i++)
    {
        for (int j = i; j <= 6; j++)
        {
            cache_p[index_p(cu, i, j, 0)] =
                0.5 * (cache_p[index_p(cc, i, j, 0)] + cache_p[index_p(uu, i, j, 0)]);
            cache_ps[index_p(cu, i, j, 0)] =
                0.5 * (cache_ps[index_p(cc, i, j, 0)] + cache_ps[index_p(uu, i, j, 0)]);
        }
    }
    double L_1 = 2. * log(mu_1 / Mb);
    double L_2 = 2. * log(mu_2 / Mb);
    double log_some_z = log((2 * sqrt1minus4z) / (1 + sqrt1minus4z));

    cache_p[index_p(cc, 1, 1, 1)] = (169. * Dilogsigma) / 27. + (92. * Dilogsigma2) / 27. + (92. * log2sigma) / 27. - (191. * logsigma) / 54. +
                                    (46. * log1minus4z * logsigma) / 27. + (169. * logsigma * log_some_z) / 54. - (92. * logsigma * logz) / 27. +
                                    (1211. * sqrt1minus4z) / 324. + (19. * L_1 * sqrt1minus4z) / 18. + (149. * L_2 * sqrt1minus4z) / 108. -
                                    (8. * log1minus4z * sqrt1minus4z) / 3. + (43. * logz * sqrt1minus4z) / 12. - (5. * logsigma) / (216. * z) +
                                    (5. * logz * sqrt1minus4z) / (216. * z) - (340. * Dilogsigma * z) / 9. - 19. * Dilogsigma2 * z -
                                    19. * log2sigma * z + (371. * logsigma * z) / 27. - (19. * log1minus4z * logsigma * z) / 2. -
                                    (170. * logsigma * log_some_z * z) / 9. + 19. * logsigma * logz * z - (6917. * sqrt1minus4z * z) / 324. -
                                    (23. * L_1 * sqrt1minus4z * z) / 9. - (49. * L_2 * sqrt1minus4z * z) / 54. +
                                    (128. * log1minus4z * sqrt1minus4z * z) / 9. - (1951. * logz * sqrt1minus4z * z) / 108. +
                                    (1364. * Dilogsigma * z2) / 27. + (682. * Dilogsigma2 * z2) / 27. + (682. * log2sigma * z2) / 27. -
                                    (263. * logsigma * z2) / 54. + (341. * log1minus4z * logsigma * z2) / 27. +
                                    (682. * logsigma * log_some_z * z2) / 27. - (682. * logsigma * logz * z2) / 27. -
                                    (295. * sqrt1minus4z * z2) / 54. + (295. * logsigma * z3) / 27. - (5. * (M_PI2)) / 108. + (z * (M_PI2)) / 54.;
    cache_p[index_p(cc, 1, 2, 1)] = (92. * Dilogsigma) / 9. + (16. * Dilogsigma2) / 9. + (16. * log2sigma) / 9. - (275. * logsigma) / 36. +
                                    (8. * log1minus4z * logsigma) / 9. + (46. * logsigma * log_some_z) / 9. - (16. * logsigma * logz) / 9. -
                                    (476. * sqrt1minus4z) / 27. - (37. * L_1 * sqrt1minus4z) / 6. + (19. * L_2 * sqrt1minus4z) / 9. +
                                    (27. * logz * sqrt1minus4z) / 4. + (4. * logsigma) / (9. * z) - (4. * logz * sqrt1minus4z) / (9. * z) -
                                    (176. * Dilogsigma * z) / 3. - 28. * Dilogsigma2 * z - 28. * log2sigma * z + (815. * logsigma * z) / 18. -
                                    14. * log1minus4z * logsigma * z - (88. * logsigma * log_some_z * z) / 3. + 28. * logsigma * logz * z +
                                    (917. * sqrt1minus4z * z) / 27. + (41. * L_1 * sqrt1minus4z * z) / 3. + (2. * L_2 * sqrt1minus4z * z) / 9. +
                                    (64. * log1minus4z * sqrt1minus4z * z) / 3. - (392. * logz * sqrt1minus4z * z) / 9. + (688. * Dilogsigma * z2) / 9. +
                                    (344. * Dilogsigma2 * z2) / 9. + (344. * log2sigma * z2) / 9. - (269. * logsigma * z2) / 9. +
                                    (172. * log1minus4z * logsigma * z2) / 9. + (344. * logsigma * log_some_z * z2) / 9. - (344. * logsigma * logz * z2) / 9. -
                                    (91. * sqrt1minus4z * z2) / 9. + (182. * logsigma * z3) / 9. + (5. * (M_PI2)) / 9. - (2. * z * (M_PI2)) / 9.;
    cache_p[index_p(cc, 2, 2, 1)] = (4. * Dilogsigma) / 3. + (32. * Dilogsigma2) / 3. + (32. * log2sigma) / 3. - (41. * logsigma) / 6. +
                                    (16. * log1minus4z * logsigma) / 3. + (2. * logsigma * log_some_z) / 3. - (32. * logsigma * logz) / 3. +
                                    (103. * sqrt1minus4z) / 18. - L_1 * sqrt1minus4z + (2. * L_2 * sqrt1minus4z) / 3. -
                                    12. * log1minus4z * sqrt1minus4z + (21. * logz * sqrt1minus4z) / 2. - (11. * logsigma) / (6. * z) +
                                    (11. * logz * sqrt1minus4z) / (6. * z) - 16. * Dilogsigma * z - 12. * Dilogsigma2 * z - 12. * log2sigma * z +
                                    (77. * logsigma * z) / 3. - 6. * log1minus4z * logsigma * z - 8. * logsigma * log_some_z * z + 12. * logsigma * logz * z +
                                    (34. * sqrt1minus4z * z) / 9. + 10. * L_1 * sqrt1minus4z * z - (14. * L_2 * sqrt1minus4z * z) / 3. +
                                    8. * log1minus4z * sqrt1minus4z * z - (73. * logz * sqrt1minus4z * z) / 3. + (80. * Dilogsigma * z2) / 3. +
                                    (40. * Dilogsigma2 * z2) / 3. + (40. * log2sigma * z2) / 3. - (16. * logsigma * z2) / 3. +
                                    (20. * log1minus4z * logsigma * z2) / 3. + (40. * logsigma * log_some_z * z2) / 3. - (40. * logsigma * logz * z2) / 3. +
                                    (16. * sqrt1minus4z * z2) / 3. - (32. * logsigma * z3) / 3. - (5. * (M_PI2)) / 3. + (2. * z * (M_PI2)) / 3.;
    cache_ps[index_p(cc, 1, 1, 1)] = (-344. * Dilogsigma) / 27. - (160. * Dilogsigma2) / 27. - (160. * log2sigma) / 27. + (320. * logsigma) / 27. -
                                     (80. * log1minus4z * logsigma) / 27. - (172. * logsigma * log_some_z) / 27. + (160. * logsigma * logz) / 27. -
                                     (10. * sqrt1minus4z) / 81. - (4. * L_1 * sqrt1minus4z) / 9. - (40. * L_2 * sqrt1minus4z) / 27. +
                                     (80. * log1minus4z * sqrt1minus4z) / 9. - 12. * logz * sqrt1minus4z + (2. * logsigma) / (27. * z) -
                                     (2. * logz * sqrt1minus4z) / (27. * z) + (8. * Dilogsigma2 * z) / 9. + (8. * log2sigma * z) / 9. -
                                     (214. * logsigma * z) / 27. + (4. * log1minus4z * logsigma * z) / 9. - (8. * logsigma * logz * z) / 9. -
                                     (1511. * sqrt1minus4z * z) / 81. - (8. * L_1 * sqrt1minus4z * z) / 9. - (80. * L_2 * sqrt1minus4z * z) / 27. +
                                     (160. * log1minus4z * sqrt1minus4z * z) / 9. - (610. * logz * sqrt1minus4z * z) / 27. +
                                     (1376. * Dilogsigma * z2) / 27. + (688. * Dilogsigma2 * z2) / 27. + (688. * log2sigma * z2) / 27. -
                                     (460. * logsigma * z2) / 27. + (344. * log1minus4z * logsigma * z2) / 27. + (688. * logsigma * log_some_z * z2) / 27. -
                                     (688. * logsigma * logz * z2) / 27. - (590. * sqrt1minus4z * z2) / 27. + (1180. * logsigma * z3) / 27. -
                                     (2. * (M_PI2)) / 27. - (4. * z * (M_PI2)) / 27.;
    cache_ps[index_p(cc, 1, 2, 1)] = (-160. * Dilogsigma) / 9. - (128. * Dilogsigma2) / 9. - (128. * log2sigma) / 9. + (238. * logsigma) / 9. -
                                     (64. * log1minus4z * logsigma) / 9. - (80. * logsigma * log_some_z) / 9. + (128. * logsigma * logz) / 9. +
                                     (100. * sqrt1minus4z) / 27. + (4. * L_1 * sqrt1minus4z) / 3. - (32. * L_2 * sqrt1minus4z) / 9. +
                                     (64. * log1minus4z * sqrt1minus4z) / 3. - 26. * logz * sqrt1minus4z - (2. * logsigma) / (9. * z) +
                                     (2. * logz * sqrt1minus4z) / (9. * z) - (32. * Dilogsigma2 * z) / 3. - (32. * log2sigma * z) / 3. +
                                     (16. * logsigma * z) / 9. - (16. * log1minus4z * logsigma * z) / 3. + (32. * logsigma * logz * z) / 3. -
                                     (442. * sqrt1minus4z * z) / 27. + (8. * L_1 * sqrt1minus4z * z) / 3. - (64. * L_2 * sqrt1minus4z * z) / 9. +
                                     (128. * log1minus4z * sqrt1minus4z * z) / 3. - (560. * logz * sqrt1minus4z * z) / 9. + (640. * Dilogsigma * z2) / 9. +
                                     (320. * Dilogsigma2 * z2) / 9. + (320. * log2sigma * z2) / 9. - (728. * logsigma * z2) / 9. +
                                     (160. * log1minus4z * logsigma * z2) / 9. + (320. * logsigma * log_some_z * z2) / 9. - (320. * logsigma * logz * z2) / 9. -
                                     (364. * sqrt1minus4z * z2) / 9. + (728. * logsigma * z3) / 9. + (8. * (M_PI2)) / 9. + (16. * z * (M_PI2)) / 9.;
    cache_ps[index_p(cc, 2, 2, 1)] = (-32. * Dilogsigma) / 3. + (32. * Dilogsigma2) / 3. + (32. * log2sigma) / 3. - (28. * logsigma) / 3. +
                                     (16. * log1minus4z * logsigma) / 3. - (16. * logsigma * log_some_z) / 3. - (32. * logsigma * logz) / 3. +
                                     (272. * sqrt1minus4z) / 9. + 8. * L_1 * sqrt1minus4z + (8. * L_2 * sqrt1minus4z) / 3. -
                                     16. * log1minus4z * sqrt1minus4z + 12. * logz * sqrt1minus4z - (4. * logsigma) / (3. * z) +
                                     (4. * logz * sqrt1minus4z) / (3. * z) + 32. * Dilogsigma2 * z + 32. * log2sigma * z - (40. * logsigma * z) / 3. +
                                     16. * log1minus4z * logsigma * z - 32. * logsigma * logz * z + (664. * sqrt1minus4z * z) / 9. +
                                     16. * L_1 * sqrt1minus4z * z + (16. * L_2 * sqrt1minus4z * z) / 3. - 32. * log1minus4z * sqrt1minus4z * z +
                                     (104. * logz * sqrt1minus4z * z) / 3. + (128. * Dilogsigma * z2) / 3. + (64. * Dilogsigma2 * z2) / 3. +
                                     (64. * log2sigma * z2) / 3. + (272. * logsigma * z2) / 3. + (32. * log1minus4z * logsigma * z2) / 3. +
                                     (64. * logsigma * log_some_z * z2) / 3. - (64. * logsigma * logz * z2) / 3. + (64. * sqrt1minus4z * z2) / 3. -
                                     (128. * logsigma * z3) / 3. - (8. * (M_PI2)) / 3. - (16. * z * (M_PI2)) / 3.;

    cache_p[index_p(uu, 1, 1, 1)] = 299. / 81. + (19. * L_1) / 18. + (149. * L_2) / 108. - (5. * (M_PI2)) / 108.;
    cache_p[index_p(uu, 1, 2, 1)] = -452. / 27. - (37. * L_1) / 6. + (19. * L_2) / 9. + (5. * (M_PI2)) / 9.;
    cache_p[index_p(uu, 2, 2, 1)] = 37. / 18. - L_1 + (2. * L_2) / 3. - (5. * (M_PI2)) / 3.;
    cache_ps[index_p(uu, 1, 1, 1)] = 2. / 81. - (4. * L_1) / 9. - (40. * L_2) / 27. - (2. * (M_PI2)) / 27.;
    cache_ps[index_p(uu, 1, 2, 1)] = 88. / 27. + (4. * L_1) / 3. - (32. * L_2) / 9. + (8. * (M_PI2)) / 9.;
    cache_ps[index_p(uu, 2, 2, 1)] = 248. / 9. + 8. * L_1 + (8. * L_2) / 3. - (8. * (M_PI2)) / 3.;

    // pengFlag terms from P12-36 here in full z dependence
    cache_p[index_p(cc, 1, 1, 1)] += (-5. * logsigma) / 324. - (5. * sqrt1minus4z) / 486. -
                                     (5. * L_1 * sqrt1minus4z) / 324. + (5. * logz * sqrt1minus4z) / 324. -
                                     (20. * sqrt1minus4z * z) / 243. - (5. * L_1 * sqrt1minus4z * z) / 162. +
                                     (5. * logz * sqrt1minus4z * z) / 162. + (5. * logsigma * z2) / 27. -
                                     (10. * sqrt1minus4z * z2) / 81. + (20. * logsigma * z3) / 81.;
    cache_p[index_p(cc, 1, 2, 1)] += (5. * logsigma) / 27. + (10. * sqrt1minus4z) / 81. + (5. * L_1 * sqrt1minus4z) / 27. -
                                     (5. * logz * sqrt1minus4z) / 27. + (80. * sqrt1minus4z * z) / 81. +
                                     (10. * L_1 * sqrt1minus4z * z) / 27. - (10. * logz * sqrt1minus4z * z) / 27. -
                                     (20. * logsigma * z2) / 9. + (40. * sqrt1minus4z * z2) / 27. -
                                     (80. * logsigma * z3) / 27.;
    cache_p[index_p(cc, 2, 2, 1)] += (-5. * logsigma) / 9. - (10. * sqrt1minus4z) / 27. - (5. * L_1 * sqrt1minus4z) / 9. +
                                     (5. * logz * sqrt1minus4z) / 9. - (80. * sqrt1minus4z * z) / 27. -
                                     (10. * L_1 * sqrt1minus4z * z) / 9. + (10. * logz * sqrt1minus4z * z) / 9. +
                                     (20. * logsigma * z2) / 3. - (40. * sqrt1minus4z * z2) / 9. + (80. * logsigma * z3) / 9.;
    cache_ps[index_p(cc, 1, 1, 1)] += (-2. * logsigma) / 81. - (4. * sqrt1minus4z) / 243. -
                                      (2. * L_1 * sqrt1minus4z) / 81. + (2. * logz * sqrt1minus4z) / 81. -
                                      (32. * sqrt1minus4z * z) / 243. - (4. * L_1 * sqrt1minus4z * z) / 81. +
                                      (4. * logz * sqrt1minus4z * z) / 81. + (8. * logsigma * z2) / 27. -
                                      (16. * sqrt1minus4z * z2) / 81. + (32. * logsigma * z3) / 81.;
    cache_ps[index_p(cc, 1, 2, 1)] += (8. * logsigma) / 27. + (16. * sqrt1minus4z) / 81. + (8. * L_1 * sqrt1minus4z) / 27. -
                                      (8. * logz * sqrt1minus4z) / 27. + (128. * sqrt1minus4z * z) / 81. +
                                      (16. * L_1 * sqrt1minus4z * z) / 27. - (16. * logz * sqrt1minus4z * z) / 27. -
                                      (32. * logsigma * z2) / 9. + (64. * sqrt1minus4z * z2) / 27. -
                                      (128. * logsigma * z3) / 27.;
    cache_ps[index_p(cc, 2, 2, 1)] += (-8. * logsigma) / 9. - (16. * sqrt1minus4z) / 27. - (8. * L_1 * sqrt1minus4z) / 9. +
                                      (8. * logz * sqrt1minus4z) / 9. - (128. * sqrt1minus4z * z) / 27. -
                                      (16. * L_1 * sqrt1minus4z * z) / 9. + (16. * logz * sqrt1minus4z * z) / 9. +
                                      (32. * logsigma * z2) / 3. - (64. * sqrt1minus4z * z2) / 9. + (128. * logsigma * z3) / 9.;

    cache_p[index_p(uu, 1, 1, 1)] += -5. / 486. - (5. * L_1) / 324.;
    cache_p[index_p(uu, 1, 2, 1)] += 10. / 81. + (5. * L_1) / 27.;
    cache_p[index_p(uu, 2, 2, 1)] += -10. / 27. - (5. * L_1) / 9.;
    cache_ps[index_p(uu, 1, 1, 1)] += -4. / 243. - (2. * L_1) / 81.;
    cache_ps[index_p(uu, 1, 2, 1)] += 16. / 81. + (8. * L_1) / 27.;
    cache_ps[index_p(uu, 2, 2, 1)] += -16. / 27. - (8. * L_1) / 9.;

    // Corrections from LO coefficients when resummation logz
    double z_replace = -6. * 4. / 3. * logz * z;
    cache_p[index_p(cc, 1, 1, 1)] += (-11. / 6. + 43. / 6. * z) / sqrt1minus4z * z_replace;
    cache_p[index_p(cc, 1, 2, 1)] += (-2. + 10. * z) / sqrt1minus4z * z_replace;
    cache_p[index_p(cc, 2, 2, 1)] += (-3. + 6. * z) / sqrt1minus4z * z_replace;
    cache_ps[index_p(cc, 1, 1, 1)] += (20. / 3. * z) / sqrt1minus4z * z_replace;
    cache_ps[index_p(cc, 1, 2, 1)] += (16. * z) / sqrt1minus4z * z_replace;
    cache_ps[index_p(cc, 2, 2, 1)] += (-12. * z) / sqrt1minus4z * z_replace;

    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2))
    {
        // equation (6.14)
        cache_p[index_p(qq, 1, 3, 1)] = (320. / 9. - 4. * L_1) * z + 47. / 18. * L_1 + 56. / 9. * L_2 - 5. / (18. * sqrt3) * M_PI + 1523. / 108.;
        cache_p[index_p(qq, 1, 4, 1)] = (59. / 3. * L_1 + 5. / 9. * M_PI2 + 4565. / 108.) * z - 281. / 108. * L_1 + L_2 / 54. + 5. / 18. * M_PI2 - 25. / (108. * sqrt3) * M_PI - 712. / 81.;
        cache_p[index_p(qq, 1, 5, 1)] = z * (-136. * L_1 - 192. * L_2 - 768. * logz - 16408. / 9.) + 376. / 9. * L_1 + 896. / 9. * L_2 - 40. / (9. * sqrt3) * M_PI + 318.;
        cache_p[index_p(qq, 1, 6, 1)] = z * (764. / 3. * L_1 + 8. * L_2 + 32. * logz + 8. / 9. * M_PI2 + 22850. / 27.) - 1259. / 27. * L_1 + 8. / 27. * L_2 + 40. / 9. * M_PI2 - 55. / (27. * sqrt3) * M_PI - 4243. / 27.;
        cache_p[index_p(qq, 2, 3, 1)] = (24. * L_1 + 170. / 3.) * z - 47. / 3. * L_1 + 14. / 3. * L_2 + 5. / (3. * sqrt3) * M_PI - 677. / 18.;
        cache_p[index_p(qq, 2, 4, 1)] = (26. * L_1 - 10. / 3. * M_PI2 + 1429. / 18.) * z - 35. / 9. * L_1 - L_2 / 9. - 5. / 3. * M_PI2 + 25. / (18. * sqrt3) * M_PI - 88. / 27.;
        cache_p[index_p(qq, 2, 5, 1)] = z * (816. * L_1 - 144. * L_2 - 576. * logz + 3656. / 3.) - 752. / 3. * L_1 + 224. / 3. * L_2 + 80. / (3. * sqrt3) * M_PI - 580.;
        cache_p[index_p(qq, 2, 6, 1)] = z * (128. * L_1 - 48. * L_2 - 192. * logz - 16. / 3. * M_PI2 + 6140. / 9.) - 290. / 9. * L_1 - 16. / 9. * L_2 - 80. / 3. * M_PI2 + 110. / (9. * sqrt3) * M_PI - 442. / 9.;

        cache_ps[index_p(qq, 1, 3, 1)] = -4. / 3. * L_1 - 64. / 9. * L_2 - 1720. / 9. * z - 4. / (9. * sqrt3) * M_PI - 130. / 27.;
        cache_ps[index_p(qq, 1, 4, 1)] = 2. * L_1 - 16. / 27. * L_2 + (80. / 27. + 8. / 9. * M_PI2) * z + 4. / 9. * M_PI2 - 10. / (27. * sqrt3) * M_PI + 404. / 81.;
        cache_ps[index_p(qq, 1, 5, 1)] = -64. / 3. * L_1 - 1024. / 9. * L_2 - 27952. / 9. * z - 64. / (9. * sqrt3) * M_PI - 2128. / 9.;
        cache_ps[index_p(qq, 1, 6, 1)] = 24. * L_1 - 256. / 27. * L_2 + (128. / 9. * M_PI2 - 520. / 27.) * z - 64. / 9. * M_PI2 - 88. / (27. * sqrt3) * M_PI + 2824. / 27.;
        cache_ps[index_p(qq, 2, 3, 1)] = 8. * L_1 - 16. / 3. * L_2 - 448. / 3. * z + 8. / (3. * sqrt3) * M_PI + 116. / 9.;
        cache_ps[index_p(qq, 2, 4, 1)] = 32. / 9. * L_2 + (488. / 9. - 16. / 3. * M_PI2) * z - 8. / 3. * M_PI2 + 20. / (9. * sqrt3) * M_PI + 272. / 27.;
        cache_ps[index_p(qq, 2, 5, 1)] = 128. * L_1 - 256. / 3. * L_2 - 6304. / 3. * z + 128. / (3. * sqrt3) * M_PI + 32. / 3.;
        cache_ps[index_p(qq, 2, 6, 1)] = 48. * L_1 + 512. / 9. * L_2 + (7520. / 9. - 256. / 3. * M_PI2) * z - 128. / 3. * M_PI2 + 176. / (9. * sqrt3) * M_PI + 1840. / 9.;

        // Gerlach thesis eq. (6.13) and (6.15)
        z = 0.;
    }
    z = cache_z;
    // Gerlach thesis eq. (6.15)
    cache_p[index_p(uu, 1, 4, 1)] += 5. / 9. * z;
    cache_p[index_p(uu, 1, 6, 1)] += 50. / 9. * z;
    cache_p[index_p(uu, 2, 4, 1)] += -10. / 3. * z;
    cache_p[index_p(uu, 2, 6, 1)] += -100. / 3. * z;

    cache_ps[index_p(uu, 1, 4, 1)] += 8. / 9. * z;
    cache_ps[index_p(uu, 1, 6, 1)] += 80. / 9. * z;
    cache_ps[index_p(uu, 2, 4, 1)] += -16. / 3. * z;
    cache_ps[index_p(uu, 2, 6, 1)] += -160. / 3. * z;

    // Corrections from LO coefficients when resummation logz
    cache_p[index_p(cc, 1, 5, 1)] += -96. * z_replace;
    cache_p[index_p(cc, 1, 6, 1)] += 4. * z_replace;
    cache_p[index_p(cc, 2, 5, 1)] += -72. * z_replace;
    cache_p[index_p(cc, 2, 6, 1)] += -24. * z_replace;

    // Gerlach thesis eq. (6.13) and (6.16)
    for (int i = 1; i <= 2; i++)
    {
        for (int j = i; j <= 6; j++)
        {
            cache_p[index_p(cu, i, j, 1)] =
                0.5 * (cache_p[index_p(cc, i, j, 1)] + cache_p[index_p(uu, i, j, 1)]);
            cache_ps[index_p(cu, i, j, 1)] =
                0.5 * (cache_ps[index_p(cc, i, j, 1)] + cache_ps[index_p(uu, i, j, 1)]);
        }
    }

    cache_p[index_p(cu, 1, 1, 1)] = (15. * log1minusz + 2407. * z + 894. * L_2 * z - 1065. * log1minusz * z + 1014. * log1minusz * logz * z - 1539. * log1minusz * z4 * z +
                                     1539. * logz * z4 * z - 36. * L_1 * oneminusz2 * z * (-19. + 4. * z) - 9255. * z2 - 1188. * L_2 * z2 + 5679. * log1minusz * z2 - 3738. * logz * z2 -
                                     2970. * log1minusz * logz * z2 + 10008. * z3 - 306. * L_2 * z3 - 9762. * log1minusz * z3 + 7794. * logz * z3 + 3060. * log1minusz * logz * z3 -
                                     12. * Dilogz * z * (-169. + 495. * z - 510. * z2 + 184. * z3) - 3160. * z4 + 600. * L_2 * z4 + 6672. * log1minusz * z4 - 6876. * logz * z4 -
                                     1104. * log1minusz * logz * z4 - 30. * z * (M_PI2) + 6. * z2 * (M_PI2) + 24. * z3 * (M_PI2)) /
                                    (648. * z);
    cache_p[index_p(cu, 1, 2, 1)] = -1. / 108. * (48. * log1minusz + 1856. * z - 228. * L_2 * z - 627. * log1minusz * z - 276. * log1minusz * logz * z + 378. * log1minusz * z4 * z - 378. * logz * z4 * z - 18. * L_1 * oneminusz2 * z * (-37. + 4. * z) - 3588. * z2 + 216. * L_2 * z2 + 72. * log1minusz * z2 + 588. * logz * z2 + 972. * log1minusz * logz * z2 + 1143. * z3 + 252. * L_2 * z3 + 1923. * log1minusz * z3 - 2421. * logz * z3 - 792. * log1minusz * logz * z3 + 24. * Dilogz * z * (-23. + 81. * z - 66. * z2 + 8. * z3) + 589. * z4 - 240. * L_2 * z4 - 1794. * log1minusz * z4 + 1746. * logz * z4 + 96. * log1minusz * logz * z4 - 60. * z * (M_PI2) + 12. * z2 * (M_PI2) + 48. * z3 * (M_PI2)) / z;
    cache_p[index_p(cu, 2, 2, 1)] = (70. + 12. * L_2 - 210. * log1minusz + 6. * log1minusz * logz + (33. * log1minusz) / z - 3. * z - 54. * L_2 * z + 225. * log1minusz * z - 210. * logz * z +
                                     54. * log1minusz * logz * z + 18. * L_1 * oneminusz2 * (-1. + 4. * z) - 360. * z2 + 72. * L_2 * z2 + 21. * log1minusz * z2 + 153. * logz * z2 +
                                     36. * log1minusz * logz * z2 + 293. * z3 - 30. * L_2 * z3 - 42. * log1minusz * z3 - 126. * logz * z3 - 96. * log1minusz * logz * z3 -
                                     12. * Dilogz * (-1. - 9. * z - 6. * z2 + 16. * z3) - 27. * log1minusz * z4 + 27. * logz * z4 - 30. * (M_PI2) + 6. * z * (M_PI2) + 24. * z2 * (M_PI2)) /
                                    18.;
    cache_ps[index_p(cu, 1, 1, 1)] = -1. / 162. * (12. * log1minusz + 8. * z + 240. * L_2 * z - 978. * log1minusz * z + 516. * log1minusz * logz * z + 774. * log1minusz * z4 * z - 774. * logz * z4 * z + 72. * L_1 * oneminusz2 * z * (1. + 2. * z) + 1461. * z2 - 1674. * log1minusz * z2 + 516. * logz * z2 - 36. * log1minusz * logz * z2 - 3654. * z3 - 720. * L_2 * z3 + 7008. * log1minusz * z3 - 5796. * logz * z3 - 1584. * log1minusz * logz * z3 + 24. * Dilogz * z * (43. - 3. * z - 132. * z2 + 92. * z3) + 2185. * z4 + 480. * L_2 * z4 - 5142. * log1minusz * z4 + 5346. * logz * z4 + 1104. * log1minusz * logz * z4 + 12. * z * (M_PI2) + 12. * z2 * (M_PI2)-24. * z3 * (M_PI2)) / z;
    cache_ps[index_p(cu, 1, 2, 1)] = (6. * log1minusz + 94. * z - 96. * L_2 * z + 402. * log1minusz * z - 120. * log1minusz * logz * z - 180. * log1minusz * z4 * z +
                                      180. * logz * z4 * z + 36. * L_1 * oneminusz2 * z * (1. + 2. * z) - 363. * z2 + 216. * log1minusz * z2 - 120. * logz * z2 -
                                      72. * log1minusz * logz * z2 + 792. * z3 + 288. * L_2 * z3 - 1842. * log1minusz * z3 + 1638. * logz * z3 +
                                      288. * log1minusz * logz * z3 - 48. * Dilogz * z * (5. + 3. * z - 12. * z2 + 4. * z3) - 523. * z4 - 192. * L_2 * z4 +
                                      1398. * log1minusz * z4 - 1350. * logz * z4 - 96. * log1minusz * logz * z4 + 24. * z * (M_PI2) + 24. * z2 * (M_PI2)-48. * z3 * (M_PI2)) /
                                     (27. * z);
    cache_ps[index_p(cu, 2, 2, 1)] = (4. * (12. * oneminusz2 * (1. + 2. * z) + 18. * L_1 * oneminusz2 * (1. + 2. * z) - 12. * L_2 * oneminusz2 * (1. + 2. * z) +
                                            18. * (1. + L_2) * oneminusz2 * (1. + 2. * z) - 6. * (2. * Dilogz + log1minusz * logz) * (-1. + z) * (1. + 2. * z) * (-1. + 4. * z) - (3. * log1minusz * oneminusz2 * (1. + z) * (-1. + 13. * z + 3. * z2)) / z +
                                            7. * (-1. + z) * (-5. - 8. * z + 19. * z2) + 3. * logz * z * (-2. + 3. * z - 18. * z2 + 3. * z3) +
                                            6. * (-1. + z) * (1. + 2. * z) * (M_PI2))) /
                                     9.;

    // pengFlag terms from P12-36 here in full z dependence
    cache_p[index_p(cu, 1, 1, 1)] += (5. * (-2. - 3. * L_1 + 3. * logz - 12. * z - 2. * sqrt1minus4z * (1. + 2. * z) -
                                            3. * L_1 * sqrt1minus4z * (1. + 2. * z) - 3. * logsigma * sqrt1minus4z * (1. + 2. * z))) /
                                     1944.;
    cache_p[index_p(cu, 1, 2, 1)] += (5. * (2. + 3. * L_1 - 3. * logz + 12. * z + 2. * sqrt1minus4z * (1. + 2. * z) +
                                            3. * L_1 * sqrt1minus4z * (1. + 2. * z) + 3. * logsigma * sqrt1minus4z * (1. + 2. * z))) /
                                     162.;
    cache_p[index_p(cu, 2, 2, 1)] += (5. * (-2. - 3. * L_1 + 3. * logz - 12. * z - 2. * sqrt1minus4z * (1. + 2. * z) -
                                            3. * L_1 * sqrt1minus4z * (1. + 2. * z) - 3. * logsigma * sqrt1minus4z * (1. + 2. * z))) /
                                     54.;
    cache_ps[index_p(cu, 1, 1, 1)] += (-2. - 3. * L_1 + 3. * logz - 12. * z - 2. * sqrt1minus4z * (1. + 2. * z) -
                                       3. * L_1 * sqrt1minus4z * (1. + 2. * z) - 3. * logsigma * sqrt1minus4z * (1. + 2. * z)) /
                                      243.;
    cache_ps[index_p(cu, 1, 2, 1)] += (4. * (2. + 3. * L_1 - 3. * logz + 12. * z + 2. * sqrt1minus4z * (1. + 2. * z) +
                                             3. * L_1 * sqrt1minus4z * (1. + 2. * z) + 3. * logsigma * sqrt1minus4z * (1. + 2. * z))) /
                                      81.;
    cache_ps[index_p(cu, 2, 2, 1)] += (4. * (-2. - 3. * L_1 + 3. * logz - 12. * z - 2. * sqrt1minus4z * (1. + 2. * z) -
                                             3. * L_1 * sqrt1minus4z * (1. + 2. * z) - 3. * logsigma * sqrt1minus4z * (1. + 2. * z))) /
                                      27.;

    // Corrections from LO coefficients when resummation logz
    cache_p[index_p(cu, 1, 1, 1)] += (-11. / 12. + 7. / 4. * z + 5. / 6. * z2) * z_replace;
    cache_p[index_p(cu, 1, 2, 1)] += (-1. + 3. * z - 2. * z2) * z_replace;
    cache_p[index_p(cu, 2, 2, 1)] += (-3. / 2. + 3. / 2. * z2) * z_replace;
    cache_ps[index_p(cu, 1, 1, 1)] += (10. / 3. * z - 10. / 3. * z2) * z_replace;
    cache_ps[index_p(cu, 1, 2, 1)] += (8. * z - 8. * z2) * z_replace;
    cache_ps[index_p(cu, 2, 2, 1)] += (-6. * z + 6. * z2) * z_replace;

    // Gerlach thesis eq. (6.18)
    cache_p[index_p(cc, 3, 3, 1)] = -154. / 9. * L_1 + 184. / 3. * L_2 + 90. * z - 5. / 3. * M_PI2 + 5. / (3. * sqrt3) * M_PI + 1390. / 27.;
    cache_p[index_p(cc, 3, 4, 1)] = -811. / 54. * L_1 + 74. / 9. * L_2 - 10. / 3. * z - 10. / 9. * M_PI2 + 70. / (9. * sqrt3) * M_PI - 27991. / 324.;
    cache_p[index_p(cc, 3, 5, 1)] = -4928. / 9. * L_1 + 3872. / 3. * L_2 + 1800. * z - 160. / 3. * M_PI2 + 160. / (3. * sqrt3) * M_PI + 16880. / 27.;
    cache_p[index_p(cc, 3, 6, 1)] = (144. * L_1 + 440. / 3.) * z - 12932. / 27. * L_1 + 1184. / 9. * L_2 - 160. / 9. * M_PI2 + 670. / (9. * sqrt3) * M_PI - 131410. / 81.;
    cache_p[index_p(cc, 4, 4, 1)] = 181. / 162. * L_1 + 127. / 108. * L_2 + (323. / 36. - 5. / 3. * M_PI2) * z - 335. / 108. * M_PI2 + 575. / (108. * sqrt3) * M_PI + 779. / 486.;
    cache_p[index_p(cc, 4, 5, 1)] = (576. * L_1 + 3836. / 3.) * z - 14912. / 27. * L_1 + 1184. / 9. * L_2 - 160. / 9. * M_PI2 + 1120. / (9. * sqrt3) * M_PI - 127990. / 81.;
    cache_p[index_p(cc, 4, 6, 1)] = (60. * L_1 - 100. / 3. * M_PI2 + 2455. / 9.) * z - 8759. / 81. * L_1 + 1088. / 27. * L_2 - 1600. / 27. * M_PI2 + 2665. / (27. * sqrt3) * M_PI - 50083. / 243.;
    cache_p[index_p(cc, 5, 5, 1)] = z * (-2592. * L_2 - 10368. * logz - 33120.) - 39424. / 9. * L_1 + 26944. / 3. * L_2 - 1280. / 3. * M_PI2 + 1280. / (3. * sqrt3) * M_PI + 347104. / 27.;
    cache_p[index_p(cc, 5, 6, 1)] = (7200. * L_1 + 74000. / 3.) * z - 240608. / 27. * L_1 + 18944. / 9. * L_2 - 2560. / 9. * M_PI2 + 10720. / (9. * sqrt3) * M_PI - 2253568. / 81.;
    cache_p[index_p(cc, 6, 6, 1)] = z * (-48. * L_1 - 144. * L_2 - 576. * logz - 248. / 3. * M_PI2 + 12290. / 9.) - 59632. / 81. * L_1 + 8848. / 27. * L_2 - 10640. / 27. * M_PI2 + 12320. / (27. * sqrt3) * M_PI - 662144. / 243.;

    cache_ps[index_p(cc, 3, 3, 1)] = 176. / 9. * L_1 - 200. / 3. * L_2 - 432. * z - 8. / 3. * M_PI2 + 8. / (3. * sqrt3) * M_PI - 620. / 27.;
    cache_ps[index_p(cc, 3, 4, 1)] = 268. / 27. * L_1 - 64. / 9. * L_2 - 16. / 3. * z - 16. / 9. * M_PI2 + 112. / (9. * sqrt3) * M_PI + 3506. / 81.;
    cache_ps[index_p(cc, 3, 5, 1)] = 5632. / 9. * L_1 - 4096. / 3. * L_2 - 8640. * z - 256. / 3. * M_PI2 + 256. / (3. * sqrt3) * M_PI + 9728. / 27.;
    cache_ps[index_p(cc, 3, 6, 1)] = 9184. / 27. * L_1 - 1024. / 9. * L_2 - 160. / 3. * z - 256. / 9. * M_PI2 + 1072. / (9. * sqrt3) * M_PI + 88688. / 81.;
    cache_ps[index_p(cc, 4, 4, 1)] = 1028. / 81. * L_1 + 136. / 27. * L_2 - 8. / 3. * M_PI2 * z + 230. / 9. * z - 134. / 27. * M_PI2 + 230. / (27. * sqrt3) * M_PI + 6214. / 243.;
    cache_ps[index_p(cc, 4, 5, 1)] = 9472. / 27. * L_1 - 1024. / 9. * L_2 + 608. / 3. * z - 256. / 9. * M_PI2 + 1792. / (9. * sqrt3) + 64784. / 81.;
    cache_ps[index_p(cc, 4, 6, 1)] = 10792. / 81. * L_1 + 2048. / 27. * L_2 - 160. / 3. * M_PI2 * z + 3568. / 9. * z - 2560. / 27. * M_PI2 + 4264. / (27. * sqrt3) * M_PI + 123080. / 243.;
    cache_ps[index_p(cc, 5, 5, 1)] = 45056. / 9. * L_1 - 28160. / 3. * L_2 - 58752. * z - 2048. / 3. * M_PI2 - 2048. / (3. * sqrt3) * M_PI - 349184. / 27.;
    cache_ps[index_p(cc, 5, 6, 1)] = 167680. / 27. * L_1 - 16384. / 9. * L_2 + 6080. / 3. * z - 4096. / 9. * M_PI2 + 17152. / (9. * sqrt3) * M_PI + 1502720. / 81.;
    cache_ps[index_p(cc, 6, 6, 1)] = 75392. / 81. * L_1 + 11776. / 27. * L_2 - 1088. / 3. * M_PI2 * z + 23696. / 9. * z - 17024. / 27. * M_PI2 + 19712. / (27. * sqrt3) * M_PI + 717184. / 243.;

    // Corrections from LO coefficients when resummation logz
    cache_p[index_p(cc, 5, 5, 1)] += -1296. * n_v * z_replace;
    cache_p[index_p(cc, 6, 6, 1)] += -72. * n_v * z_replace;

    // Gerlach thesis eq. (6.17)
    for (int i = 3; i <= 6; i++)
    {
        for (int j = i; j <= 6; j++)
        {
            cache_p[index_p(cu, i, j, 1)] = cache_p[index_p(cc, i, j, 1)];
            cache_p[index_p(uu, i, j, 1)] = cache_p[index_p(cc, i, j, 1)];
            cache_ps[index_p(cu, i, j, 1)] = cache_ps[index_p(cc, i, j, 1)];
            cache_ps[index_p(uu, i, j, 1)] = cache_ps[index_p(cc, i, j, 1)];
        }
    }

    // Gerlach thesis eq. (6.19)
    cache_p[index_p(cc, 1, 8, 1)] = 5. / 18. * sqrt1minus4z + 5. / 9. * z * sqrt1minus4z;
    cache_p[index_p(cc, 2, 8, 1)] = -5. / 3. * sqrt1minus4z - 10. / 3. * z * sqrt1minus4z;
    cache_ps[index_p(cc, 1, 8, 1)] = 4. / 9. * sqrt1minus4z + 8. / 9. * z * sqrt1minus4z;
    cache_ps[index_p(cc, 2, 8, 1)] = -8. / 3. * sqrt1minus4z - 16. / 3. * z * sqrt1minus4z;
    cache_p[index_p(uu, 1, 8, 1)] = 5. / 18.;
    cache_p[index_p(uu, 2, 8, 1)] = -5. / 3.;
    cache_ps[index_p(uu, 1, 8, 1)] = 4. / 9.;
    cache_ps[index_p(uu, 2, 8, 1)] = -8. / 3.;

    // Gerlach thesis eq. (6.27)
    for (int i = 1; i <= 2; i++)
    {
        cache_p[index_p(cu, i, 8, 1)] = 0.5 * (cache_p[index_p(cc, i, 8, 1)] + cache_p[index_p(uu, i, 8, 1)]);
        cache_ps[index_p(cu, i, 8, 1)] = 0.5 * (cache_ps[index_p(cc, i, 8, 1)] + cache_ps[index_p(uu, i, 8, 1)]);
    }

    // Gerlach thesis eq. (6.21)
    cache_p[index_p(cc, 3, 8, 1)] = -32. / 3.;
    cache_p[index_p(cc, 4, 8, 1)] = (-139. - 30. * sqrt1minus4z - 60. * sqrt1minus4z * z) / 18.;
    cache_p[index_p(cc, 5, 8, 1)] = -512. / 3.;
    cache_p[index_p(cc, 6, 8, 1)] = (-1684. - 300. * sqrt1minus4z - 600. * sqrt1minus4z * z) / 18.;
    cache_ps[index_p(cc, 3, 8, 1)] = 64. / 3.;
    cache_ps[index_p(cc, 4, 8, 1)] = (4. * (1. - 6. * sqrt1minus4z - 12. * sqrt1minus4z * z)) / 9.;
    cache_ps[index_p(cc, 5, 8, 1)] = 1024. / 3.;
    cache_ps[index_p(cc, 6, 8, 1)] = (4 * (124. - 60. * sqrt1minus4z - 120. * sqrt1minus4z * z)) / 9.;

    // Gerlach thesis eq. (6.20)
    for (int i = 3; i <= 6; i++)
    {
        cache_p[index_p(cu, i, 8, 1)] = cache_p[index_p(cc, i, 8, 1)];
        cache_p[index_p(uu, i, 8, 1)] = cache_p[index_p(cc, i, 8, 1)];
        cache_ps[index_p(cu, i, 8, 1)] = cache_ps[index_p(cc, i, 8, 1)];
        cache_ps[index_p(uu, i, 8, 1)] = cache_ps[index_p(cc, i, 8, 1)];
    }

    double L_12 = L_1 * L_1;
    double L_22 = L_2 * L_2;

    gslpp::complex log1M1OverSqrtz = log(abs(1. - 1. / sqrtz));
    if (1. / sqrtz > 1)
        log1M1OverSqrtz += M_PI * complex::i();
    gslpp::complex logM1OverSqrtz = log(abs(-1. / sqrtz)) + M_PI * complex::i();
    gslpp::complex logM1OverSqrtz2 = logM1OverSqrtz * logM1OverSqrtz;

    // LO in z
    cache_p[index_p(cc, 1, 1, 2)] = (814589597. / 4199040. + (1320817. * L_1) / 17496. + (3211. * L_12) / 324. + (259603. * L_2) / 11664. + (12911. * L_1 * L_2) / 972. - (311. * L_22) / 216. + (56. * Cl2PI3 * sqrt3) / 729. -
                                     ((5. * complex::i()) / 3888.) * log3 * log3 * sqrt3 + (855371. * sqrtz) / 180000. + (22. * log2z * sqrtz) / 9. + (88. * logM1OverSqrtz2 * sqrtz) / 9. + (88. * logM1OverSqrtz * logz * sqrtz) / 9. -
                                     ((5. * -1.) / 486.) * sqrt3 * t_2 - (22. * log2z) / (9. * z) - (88. * logM1OverSqrtz2) / (9. * z) - (88. * logM1OverSqrtz * logz) / (9. * z) + (22. * log2z) / (9. * sqrtz) +
                                     (88. * logM1OverSqrtz2) / (9. * sqrtz) + (88. * logM1OverSqrtz * logz) / (9. * sqrtz) - (3547371902315179. * z) / 3455518320000. - (721919. * L_1 * z) / 1944. - (2347. * L_12 * z) / 54. +
                                     (1891. * L_2 * z) / 81. - (337. * L_1 * L_2 * z) / 9. + (187. * L_22 * z) / 18. - (88. * logM1OverSqrtz2 * z) / 9. - (4823. * logz * z) / 9. - (1348. * L_1 * logz * z) / 9. - (88. * L_2 * logz * z) / 3. -
                                     (88. * logM1OverSqrtz * logz * z) / 9. - (28333. * zeta3) / 486. + (128581. * z * zeta3) / 216. - (5. * sqrt3 * (20. + 30. * L_1 + 3. * log3 + 180. * z) * (M_PI)) / 17496. +
                                     ((-684288. + 684288. * sqrtz - (216641. + 87480. * L_1 - 432. * L_2 + 146016. * log12sqrt52 + 5112. * log2 - (50. * complex::i()) * sqrt3 + 158184. * sqrt5 - 684288. * sqrtz) * z +
                                       18. * (197453. + 2232. * L_1 + 912. * L_2 - 9792. * log12sqrt52 - 26508. * log2 + 576. * logz + 9744. * sqrt5) * z2) *
                                      (M_PI2)) /
                                         (69984. * z) +
                                     (23. * (M_PI4)) / 4860. - (13637. * z * (M_PI4)) / 116640.)
                                        .real();
    cache_p[index_p(cc, 1, 2, 2)] = (-95740679. / 349920. - (1026907. * L_1) / 5832. - (1751. * L_12) / 54. - (619. * L_2) / 972. + (166. * L_1 * L_2) / 81. - L_22 / 18. - (224. * Cl2PI3 * sqrt3) / 243. + ((5. * complex::i()) / 324.) * log3 * log3 * sqrt3 +
                                     (2895091. * sqrtz) / 60000. + (8. * log2z * sqrtz) / 3. + (32. * logM1OverSqrtz2 * sqrtz) / 3. + (32. * logM1OverSqrtz * logz * sqrtz) / 3. + ((10. * -1.) / 81.) * sqrt3 * t_2 - (8. * log2z) / (3. * z) -
                                     (32. * logM1OverSqrtz2) / (3. * z) - (32. * logM1OverSqrtz * logz) / (3. * z) + (8. * log2z) / (3. * sqrtz) + (32. * logM1OverSqrtz2) / (3. * sqrtz) +
                                     (32. * logM1OverSqrtz * logz) / (3. * sqrtz) + (2370097879. * z) / 1944000. + (117443. * L_1 * z) / 162. + (1193. * L_12 * z) / 9. + (6350. * L_2 * z) / 27. + (64. * L_1 * L_2 * z) / 3. + (34. * L_22 * z) / 3. -
                                     (32. * logM1OverSqrtz2 * z) / 3. + 124. * logz * z + (256. * L_1 * logz * z) / 3. - 32. * L_2 * logz * z - (32. * logM1OverSqrtz * logz * z) / 3. - (10573. * zeta3) / 324. + (85027. * z * zeta3) / 90. +
                                     (5. * sqrt3 * (20. + 30. * L_1 + 3. * log3 + 180. * z) * (M_PI)) / 1458. -
                                     ((622080. - 622080. * sqrtz - 5. * (497221. + 116640. * L_1 - 864. * L_2 - 39744. * log12sqrt52 + 85824. * log2 - (100. * complex::i()) * sqrt3 - 43056. * sqrt5 + 124416. * sqrtz) * z +
                                       72. * (64865. + 3960. * L_1 + 2280. * L_2 + 3168. * log12sqrt52 + 35070. * log2 + 1440. * logz - 3288. * sqrt5) * z2) *
                                      (M_PI2)) /
                                         (58320. * z) -
                                     (799. * (M_PI4)) / 810. + (20833. * z * (M_PI4)) / 4860.)
                                        .real();
    cache_p[index_p(cc, 2, 2, 2)] = (74041. / 14580. + (106199. * L_1) / 972. + (239. * L_12) / 18. - (5117. * L_2) / 81. - (202. * L_1 * L_2) / 27. - (19. * L_22) / 3. + (224. * Cl2PI3 * sqrt3) / 81. - ((5. * complex::i()) / 108.) * log3 * log3 * sqrt3 -
                                     (1052759. * sqrtz) / 10000. + 4. * log2z * sqrtz + 16. * logM1OverSqrtz2 * sqrtz + 16. * logM1OverSqrtz * logz * sqrtz - ((10. * -1.) / 27.) * sqrt3 * t_2 - (4. * log2z) / z -
                                     (16. * logM1OverSqrtz2) / z - (16. * logM1OverSqrtz * logz) / z + (4. * log2z) / sqrtz + (16. * logM1OverSqrtz2) / sqrtz + (16. * logM1OverSqrtz * logz) / sqrtz -
                                     (69930883. * z) / 101250. - (16463. * L_1 * z) / 54. - (122. * L_12 * z) / 3. - (109. * L_2 * z) / 18. - 22. * L_1 * L_2 * z + 17. * L_22 * z - 16. * logM1OverSqrtz2 * z - 496. * logz * z - 88. * L_1 * logz * z -
                                     48. * L_2 * logz * z - 16. * logM1OverSqrtz * logz * z - (3157. * zeta3) / 54. + (2521. * z * zeta3) / 15. - (5. * sqrt3 * (20. + 30. * L_1 + 3. * log3 + 180. * z) * (M_PI)) / 486. +
                                     (-177247. / 3888. - (4. * log12sqrt52) / 9. + (148. * log2) / 27. + ((25. * complex::i()) / 972.) * sqrt3 - (13. * sqrt5) / 27. + 16. * sqrtz - 16. / z + 16. / sqrtz + (5369. * z) / 108. - (16. * log12sqrt52 * z) / 15. -
                                      (178. * log2 * z) / 9. + (16. * logz * z) / 3. + (28. * sqrt5 * z) / 45. + L_1 * (-15. + (26. * z) / 3.) + (2. * L_2 * (1. + 38. * z)) / 9.) *
                                         (M_PI2) +
                                     (971. * (M_PI4)) / 540. + (5203. * z * (M_PI4)) / 3240.)
                                        .real();

    cache_ps[index_p(cc, 1, 1, 2)] = (-67489177. / 262440. - (77617. * L_1) / 2187. - (902. * L_12) / 243. -
                                      (5504. * L_2) / 729. - (3064. * L_1 * L_2) / 243. + (260. * L_22) / 27. +
                                      (104. * Cl2PI3 * sqrt3) / 729. - (complex::i() / 486.) * log3 * log3 * sqrt3 +
                                      (166687. * sqrtz) / 6000. - ((4. * -1.) / 243.) * sqrt3 * t_2 -
                                      (27341897029. * z) / 29160000. - (97999. * L_1 * z) / 243. - (9272. * L_2 * z) / 81. -
                                      (9272. * logz * z) / 27. + (28528. * zeta3) / 243. + (29. * z * zeta3) / 3. -
                                      (sqrt3 * (20. + 30. * L_1 + 3. * log3 + 180. * z) * (M_PI)) / 2187. -
                                      ((44209. - 37152. * log12sqrt52 + 126216. * log2 - (10. * complex::i()) * sqrt3 -
                                        40248. * sqrt5 + 127854. * z - 197208. * log2 * z + 10368. * logz * z +
                                        37152. * sqrt5 * z + 17496. * L_1 * (1. + 2. * z) + 1728. * L_2 * (1. + 2. * z)) *
                                       (M_PI2)) /
                                          8748. +
                                      (449. * (M_PI4)) / 1215. - (27529. * z * (M_PI4)) / 14580.)
                                         .real();
    cache_ps[index_p(cc, 1, 2, 2)] = (-1336127. / 2187. - (28733. * L_1) / 729. + (44. * L_12) / 81. - (2176. * L_2) / 243. -
                                      (1856. * L_1 * L_2) / 81. + (208. * L_22) / 9. - (416. * Cl2PI3 * sqrt3) / 243. +
                                      ((2. * complex::i()) / 81.) * log3 * log3 * sqrt3 + (8773. * sqrtz) / 100. + ((16. * -1.) / 81.) * sqrt3 * t_2 +
                                      (1672496939. * z) / 4860000. - (23372. * L_1 * z) / 81. - (5248. * L_2 * z) / 27. -
                                      (5248. * logz * z) / 9. + (13934. * zeta3) / 81. - 244. * z * zeta3 +
                                      (4. * sqrt3 * (20. + 30. * L_1 + 3. * log3 + 180. * z) * (M_PI)) / 729. +
                                      ((39995. + 8640. * log12sqrt52 - 29232. * log2 - (20. * complex::i()) * sqrt3 + 9360. * sqrt5 +
                                        58284. * z + 11232. * log2 * z + 20736. * logz * z - 8640. * sqrt5 * z +
                                        23328. * L_1 * (1. + 2. * z) + 3456. * L_2 * (1. + 2. * z)) *
                                       (M_PI2)) /
                                          1458. +
                                      (226. * (M_PI4)) / 405. + (5692. * z * (M_PI4)) / 1215.)
                                         .real();
    cache_ps[index_p(cc, 2, 2, 2)] = (27476329. / 58320. + (40370. * L_1) / 243. + (604. * L_12) / 27. + (6928. * L_2) / 81. +
                                      (1064. * L_1 * L_2) / 27. - (52. * L_22) / 3. + (416. * Cl2PI3 * sqrt3) / 81. -
                                      ((2. * complex::i()) / 27.) * log3 * log3 * sqrt3 - (26319. * sqrtz) / 250. - ((16. * -1.) / 27.) * sqrt3 * t_2 +
                                      (287805209. * z) / 81000. + (18836. * L_1 * z) / 27. + (928. * L_2 * z) / 9. + (928. * logz * z) / 3. -
                                      (4388. * zeta3) / 27. - 600. * z * zeta3 -
                                      (4. * sqrt3 * (20. + 30. * L_1 + 3. * log3 + 180. * z) * (M_PI)) / 243. -
                                      ((41279. - 1728. * log12sqrt52 + 11808. * log2 - (20. * complex::i()) * sqrt3 - 1872. * sqrt5 +
                                        144684. * z - 14688. * log2 * z + 20736. * logz * z + 1728. * sqrt5 * z +
                                        11664. * L_1 * (1. + 2. * z) + 3456. * L_2 * (1. + 2. * z)) *
                                       (M_PI2)) /
                                          486. +
                                      (398. * (M_PI4)) / 135. + (7991. * z * (M_PI4)) / 405.)
                                         .real();

    cache_p[index_p(uu, 1, 1, 2)] = (814589597. / 4199040. + (1320817. * L_1) / 17496. + (3211. * L_12) / 324. +
                                     (259603. * L_2) / 11664. + (12911. * L_1 * L_2) / 972. - (311. * L_22) / 216. +
                                     (56. * Cl2PI3 * sqrt3) / 729. - ((5. * complex::i()) / 3888.) * log3 * log3 * sqrt3 +
                                     (855371. * sqrtz) / 180000. - ((5. * -1.) / 486.) * sqrt3 * t_2 +
                                     (5298817581707. * z) / 1151839440000. + (5. * L_1 * z) / 81. - (50. * logz * z) / 9. -
                                     (28333. * zeta3) / 486. - (5. * (20. + 30. * L_1 + 3. * log3) * sqrt3 * (M_PI)) / 17496. -
                                     ((216641. + 87480. * L_1 - 432. * L_2 + 146016. * log12sqrt52 + 5112. * log2 -
                                       (50. * complex::i()) * sqrt3 + 158184. * sqrt5) *
                                      (M_PI2)) /
                                         69984. +
                                     (23. * (M_PI4)) / 4860.)
                                        .real();
    cache_p[index_p(uu, 1, 2, 2)] = (-95740679. / 349920. - (1026907. * L_1) / 5832. - (1751. * L_12) / 54. - (619. * L_2) / 972. +
                                     (166. * L_1 * L_2) / 81. - L_22 / 18. - (224. * Cl2PI3 * sqrt3) / 243. +
                                     ((5. * complex::i()) / 324.) * log3 * log3 * sqrt3 + (2895091. * sqrtz) / 60000. +
                                     ((10. * -1.) / 81.) * sqrt3 * t_2 - (55644107. * z) / 648000. - (20. * L_1 * z) / 27. +
                                     (8. * logz * z) / 3. - (10573. * zeta3) / 324. +
                                     (5. * (20. + 30. * L_1 + 3. * log3) * sqrt3 * (M_PI)) / 1458. +
                                     ((497221. + 116640. * L_1 - 864. * L_2 - 39744. * log12sqrt52 + 85824. * log2 -
                                       (100. * complex::i()) * sqrt3 - 43056. * sqrt5) *
                                      (M_PI2)) /
                                         11664. -
                                     (799. * (M_PI4)) / 810.)
                                        .real();
    cache_p[index_p(uu, 2, 2, 2)] = (74041. / 14580. + (106199. * L_1) / 972. + (239. * L_12) / 18. - (5117. * L_2) / 81. -
                                     (202. * L_1 * L_2) / 27. - (19. * L_22) / 3. + (224. * Cl2PI3 * sqrt3) / 81. -
                                     ((5. * complex::i()) / 108.) * log3 * log3 * sqrt3 - (1052759. * sqrtz) / 10000. -
                                     ((10. * -1.) / 27.) * sqrt3 * t_2 + (9532631. * z) / 135000. + (20. * L_1 * z) / 9. - 32. * logz * z -
                                     (3157. * zeta3) / 54. - (5. * (20. + 30. * L_1 + 3. * log3) * sqrt3 * (M_PI)) / 486. -
                                     ((177247. + 58320. * L_1 - 864. * L_2 + 1728. * log12sqrt52 - 21312. * log2 -
                                       (100. * complex::i()) * sqrt3 + 1872. * sqrt5) *
                                      (M_PI2)) /
                                         3888. +
                                     (971. * (M_PI4)) / 540.)
                                        .real();
    cache_ps[index_p(uu, 1, 1, 2)] = (-67489177. / 262440. - (77617. * L_1) / 2187. - (902. * L_12) / 243. - (5504. * L_2) / 729. -
                                      (3064. * L_1 * L_2) / 243. + (260. * L_22) / 27. + (104. * Cl2PI3 * sqrt3) / 729. -
                                      (complex::i() / 486.) * log3 * log3 * sqrt3 + (166687. * sqrtz) / 6000. - ((4. * -1.) / 243.) * sqrt3 * t_2 -
                                      (261095543. * z) / 9720000. + (8. * L_1 * z) / 81. + (28528. * zeta3) / 243. -
                                      ((20. + 30. * L_1 + 3. * log3) * sqrt3 * (M_PI)) / 2187. -
                                      ((44209. + 17496. * L_1 + 1728. * L_2 - 37152. * log12sqrt52 + 126216. * log2 -
                                        (10. * complex::i()) * sqrt3 - 40248. * sqrt5) *
                                       (M_PI2)) /
                                          8748. +
                                      (449. * (M_PI4)) / 1215.)
                                         .real();
    cache_ps[index_p(uu, 1, 2, 2)] = (-1336127. / 2187. - (28733. * L_1) / 729. + (44. * L_12) / 81. - (2176. * L_2) / 243. -
                                      (1856. * L_1 * L_2) / 81. + (208. * L_22) / 9. - (416. * Cl2PI3 * sqrt3) / 243. +
                                      ((2. * complex::i()) / 81.) * log3 * log3 * sqrt3 + (8773. * sqrtz) / 100. + ((16. * -1.) / 81.) * sqrt3 * t_2 -
                                      (117168887. * z) / 1620000. - (32. * L_1 * z) / 27. + (13934. * zeta3) / 81. +
                                      (4. * (20. + 30. * L_1 + 3. * log3) * sqrt3 * (M_PI)) / 729. +
                                      ((39995. + 23328. * L_1 + 3456. * L_2 + 8640. * log12sqrt52 - 29232. * log2 -
                                        (20. * complex::i()) * sqrt3 + 9360. * sqrt5) *
                                       (M_PI2)) /
                                          1458. +
                                      (226. * (M_PI4)) / 405.)
                                         .real();
    cache_ps[index_p(uu, 2, 2, 2)] = (27476329. / 58320. + (40370. * L_1) / 243. + (604. * L_12) / 27. + (6928. * L_2) / 81. +
                                      (1064. * L_1 * L_2) / 27. - (52. * L_22) / 3. + (416. * Cl2PI3 * sqrt3) / 81. -
                                      ((2. * complex::i()) / 27.) * log3 * log3 * sqrt3 - (26319. * sqrtz) / 250. - ((16. * -1.) / 27.) * sqrt3 * t_2 +
                                      (4778443. * z) / 27000. + (32. * L_1 * z) / 9. - (4388. * zeta3) / 27. -
                                      (4. * (20. + 30. * L_1 + 3. * log3) * sqrt3 * (M_PI)) / 243. -
                                      ((41279. + 11664. * L_1 + 3456. * L_2 - 1728. * log12sqrt52 + 11808. * log2 -
                                        (20. * complex::i()) * sqrt3 - 1872. * sqrt5) *
                                       (M_PI2)) /
                                          486. +
                                      (398. * (M_PI4)) / 135.)
                                         .real();

    // Corrections from LO and NLO coefficients (incl. pengFlag)
    // resummation of logz
    cache_p[index_p(cc, 1, 1, 2)] += (-2. * (1661. + 1113. * log2z + 132. * M_PI2) * z + logz * (-30. + (39727. + 12132. * L_1 + 2376. * L_2 - 12. * M_PI2) * z)) / 81.;
    cache_p[index_p(cc, 1, 2, 2)] += (-8. * ((151 + 147. * log2z + 12. * M_PI2) * z + 2 * logz * (-12. + (448. + 144. * L_1 - 54. * L_2 - 3. * M_PI2) * z))) / 27.;
    cache_p[index_p(cc, 2, 2, 2)] += (-4. * (66. * logz - 2. * logz * (518. + 99. * L_1 + 54. * L_2 - 6. * M_PI2) * z + (151. + 21. * log2z + 12. * M_PI2) * z)) / 9.;
    cache_ps[index_p(cc, 1, 1, 2)] += (8. * logz * (4. + (507. + 172. * logz + 4. * M_PI2) * z)) / 27.;
    cache_ps[index_p(cc, 1, 2, 2)] += (32. * logz * (-1. + 4. * (12. + 5. * logz - M_PI2) * z)) / 9.;
    cache_ps[index_p(cc, 2, 2, 2)] += (32. * logz * (-2. + (-9. + 4. * logz + 4. * M_PI2) * z)) / 3.;

    // resummation affecting arguments of logs and Dilogs
    cache_p[index_p(cc, 1, 1, 2)] += 2. / 81. * logz * (15. + 1615. * z + 1014. * z * logz);
    cache_p[index_p(cc, 1, 2, 2)] += 4. / 27. * logz * (-48. + 973. * z + 276. * z * logz);
    cache_p[index_p(cc, 2, 2, 2)] += 8. / 9. * logz * (33. + 4. * z + 6. * z * logz);
    cache_ps[index_p(cc, 1, 1, 2)] += -32. / 27. * logz + 5216. / 27. * logz * z - 1376. / 27. * log2z * z;
    cache_ps[index_p(cc, 1, 2, 2)] += 32. / 9. * logz + 3712. / 9. * logz * z - 640. / 9. * log2z * z;
    cache_ps[index_p(cc, 2, 2, 2)] += 64. / 3. * logz - 640. / 3. * logz * z - 128. / 3. * log2z * z;

    double L_b = 2. * log(mu_b / Mb);

    // Transforming mb in NLO coefficients from Pole scheme to MSbar
    cache_p[index_p(cc, 1, 1, 2)] += (8. * (4. + 3. * L_b) * (-196. + 675. * z)) / 243.;
    cache_p[index_p(cc, 1, 2, 2)] += (-44. * (4. + 3. * L_b) * (-19. + 108. * z)) / 81.;
    cache_p[index_p(cc, 2, 2, 2)] += (-16. * (4. + 3. * L_b) * (-4. + 27. * z)) / 27.;
    cache_ps[index_p(cc, 1, 1, 2)] += (1264. * (4. + 3. * L_b)) / 243.;
    cache_ps[index_p(cc, 1, 2, 2)] += (416. * (4. + 3. * L_b)) / 81.;
    cache_ps[index_p(cc, 2, 2, 2)] += (-704. * (4. + 3. * L_b)) / 27.;

    cache_p[index_p(uu, 1, 1, 2)] += (-1568. * (4. + 3. * L_1)) / 243.;
    cache_p[index_p(uu, 1, 2, 2)] += (836. * (4. + 3. * L_1)) / 81.;
    cache_p[index_p(uu, 2, 2, 2)] += (64. * (4. + 3. * L_1)) / 27.;
    cache_ps[index_p(uu, 1, 1, 2)] += (1264. * (4. + 3. * L_b)) / 243.;
    cache_ps[index_p(uu, 1, 2, 2)] += (416. * (4. + 3. * L_b)) / 81.;
    cache_ps[index_p(uu, 2, 2, 2)] += (-704. * (4. + 3. * L_b)) / 27.;

    // Gerlach thesis eq. (6.27)
    for (int i = 1; i <= 2; i++)
    {
        for (int j = i; j <= 2; j++)
        {
            cache_p[index_p(cu, i, j, 2)] = 0.5 * (cache_p[index_p(cc, i, j, 2)] + cache_p[index_p(uu, i, j, 2)]);
            cache_ps[index_p(cu, i, j, 2)] = 0.5 * (cache_ps[index_p(cc, i, j, 2)] + cache_ps[index_p(uu, i, j, 2)]);
        }
    }

    // LO in z
    cache_p[index_p(cu, 1, 1, 2)] = (814589597. / 4199040. + (56. * Cl2PI3) / (243. * sqrt3) + (1320817. * L_1) / 17496. + (3211. * L_12) / 324. +
                                     (259603. * L_2) / 11664. + (12911. * L_1 * L_2) / 972. - (311. * L_22) / 216. - (((5. * complex::i()) / 1296.) * log3 * log3) / sqrt3 +
                                     (855371. * sqrtz) / 180000. + (11. * log2z * sqrtz) / 9. + (44. * logM1OverSqrtz2 * sqrtz) / 9. +
                                     (44. * logM1OverSqrtz * logz * sqrtz) / 9. - (((5. * -1.) / 162.) * t_2) / sqrt3 - (11. * log2z) / (9. * z) -
                                     (44. * logM1OverSqrtz2) / (9. * z) - (44. * logM1OverSqrtz * logz) / (9. * z) + (11. * log2z) / (9. * sqrtz) +
                                     (44. * logM1OverSqrtz2) / (9. * sqrtz) + (44. * logM1OverSqrtz * logz) / (9. * sqrtz) -
                                     (1765737724785029. * z) / 3455518320000. - (721799. * L_1 * z) / 3888. - (2347. * L_12 * z) / 108. + (1891. * L_2 * z) / 162. -
                                     (337. * L_1 * L_2 * z) / 18. + (187. * L_22 * z) / 36. - (44. * logM1OverSqrtz2 * z) / 9. - (4873. * logz * z) / 18. -
                                     (674. * L_1 * logz * z) / 9. - (44. * L_2 * logz * z) / 3. - (44. * logM1OverSqrtz * logz * z) / 9. - (28333. * zeta3) / 486. +
                                     (128581. * z * zeta3) / 432. - (5. * (3. * sqrt3 * log3 + 30. * L_1 * sqrt3 + 10. * sqrt3 * (2. + 9. * z)) * (M_PI)) / 17496. +
                                     ((-342144. + 342144. * sqrtz + (-216641. + (50. * complex::i()) * sqrt3 - 87480. * L_1 + 432. * L_2 - 5112. * log2 - 158184. * sqrt5 + 342144. * sqrtz) * z + 9. * (197453. + 2232. * L_1 + 912. * L_2 - 26508. * log2 + 576. * logz + 9744. * sqrt5) * z2 - 864. * z * (169. + 102. * z) * log12sqrt52) * (M_PI2)) / (69984. * z) +
                                     (23. * (M_PI4)) / 4860. - (13637. * z * (M_PI4)) / 233280.)
                                        .real();
    cache_p[index_p(cu, 1, 2, 2)] = (-95740679. / 349920. - (1026907. * L_1) / 5832. - (1751. * L_12) / 54. - (619. * L_2) / 972. + (166. * L_1 * L_2) / 81. - L_22 / 18. -
                                     (224. * Cl2PI3 * sqrt3) / 243. + ((5. * complex::i()) / 324.) * log3 * log3 * sqrt3 + (2895091. * sqrtz) / 60000. + (4. * log2z * sqrtz) / 3. +
                                     (16. * logM1OverSqrtz2 * sqrtz) / 3. + (16. * logM1OverSqrtz * logz * sqrtz) / 3. + ((10. * -1.) / 81.) * sqrt3 * t_2 -
                                     (4. * log2z) / (3. * z) - (16. * logM1OverSqrtz2) / (3. * z) - (16. * logM1OverSqrtz * logz) / (3. * z) +
                                     (4. * log2z) / (3. * sqrtz) + (16. * logM1OverSqrtz2) / (3. * sqrtz) +
                                     (16. * logM1OverSqrtz * logz) / (3. * sqrtz) + (1101582779. * z) / 1944000. + (117323. * L_1 * z) / 324. +
                                     (1193. * L_12 * z) / 18. + (3175. * L_2 * z) / 27. + (32. * L_1 * L_2 * z) / 3. + (17. * L_22 * z) / 3. - (16. * logM1OverSqrtz2 * z) / 3. +
                                     (190. * logz * z) / 3. + (128. * L_1 * logz * z) / 3. - 16. * L_2 * logz * z - (16. * logM1OverSqrtz * logz * z) / 3. -
                                     (10573. * zeta3) / 324. + (85027. * z * zeta3) / 180. + (5. * sqrt3 * (20. + 30. * L_1 + 3. * log3 + 90. * z) * (M_PI)) / 1458. -
                                     ((311040. - 311040. * sqrtz - 5. * (497221. + 116640. * L_1 - 864. * L_2 - 39744. * log12sqrt52 + 85824. * log2 - (100. * complex::i()) * sqrt3 - 43056. * sqrt5 + 62208. * sqrtz) * z +
                                       36. * (64865. + 3960. * L_1 + 2280. * L_2 + 3168. * log12sqrt52 + 35070. * log2 + 1440. * logz - 3288. * sqrt5) * z2) *
                                      (M_PI2)) /
                                         (58320. * z) -
                                     (799. * (M_PI4)) / 810. + (20833. * z * (M_PI4)) / 9720.)
                                        .real();
    cache_p[index_p(cu, 2, 2, 2)] = (74041. / 14580. + (106199. * L_1) / 972. + (239. * L_12) / 18. - (5117. * L_2) / 81. - (202. * L_1 * L_2) / 27. - (19. * L_22) / 3. +
                                     (224. * Cl2PI3 * sqrt3) / 81. - ((5. * complex::i()) / 108.) * log3 * log3 * sqrt3 - (1052759. * sqrtz) / 10000. + 2. * log2z * sqrtz +
                                     8. * logM1OverSqrtz2 * sqrtz + 8. * logM1OverSqrtz * logz * sqrtz - ((10. * -1.) / 27.) * sqrt3 * t_2 - (2. * log2z) / z -
                                     (8. * logM1OverSqrtz2) / z - (8. * logM1OverSqrtz * logz) / z + (2. * log2z) / sqrtz +
                                     (8. * logM1OverSqrtz2) / sqrtz + (8. * logM1OverSqrtz * logz) / sqrtz - (251125639. * z) / 810000. -
                                     (16343. * L_1 * z) / 108. - (61. * L_12 * z) / 3. - (109. * L_2 * z) / 36. - 11. * L_1 * L_2 * z + (17. * L_22 * z) / 2. -
                                     8. * logM1OverSqrtz2 * z - 264. * logz * z - 44. * L_1 * logz * z - 24. * L_2 * logz * z - 8. * logM1OverSqrtz * logz * z -
                                     (3157. * zeta3) / 54. + (2521. * z * zeta3) / 30. - (5. * sqrt3 * (20. + 30. * L_1 + 3. * log3 + 90. * z) * (M_PI)) / 486. +
                                     (-177247. / 3888. - (4. * log12sqrt52) / 9. + (148. * log2) / 27. + ((25. * complex::i()) / 972.) * sqrt3 - (13. * sqrt5) / 27. +
                                      8. * sqrtz - 8. / z + 8. / sqrtz + (5369. * z) / 216. - (8. * log12sqrt52 * z) / 15. - (89. * log2 * z) / 9. +
                                      (8. * logz * z) / 3. + (14. * sqrt5 * z) / 45. + L_1 * (-15. + (13. * z) / 3.) + (2. * L_2 * (1. + 19. * z)) / 9.) *
                                         (M_PI2) +
                                     (971. * (M_PI4)) / 540. + (5203. * z * (M_PI4)) / 6480.)
                                        .real();
    cache_ps[index_p(cu, 1, 1, 2)] = (-67489177. / 262440. + (104. * Cl2PI3) / (243. * sqrt3) - (77617. * L_1) / 2187. -
                                      (902. * L_12) / 243. - (5504. * L_2) / 729. - (3064. * L_1 * L_2) / 243. + (260. * L_22) / 27. -
                                      ((complex::i() / 162.) * log3 * log3) / sqrt3 + (166687. * sqrtz) / 6000. - (((4. * -1.) / 81.) * t_2) / sqrt3 -
                                      (14062591829. * z) / 29160000. - (97975. * L_1 * z) / 486. - (4636. * L_2 * z) / 81. -
                                      (4636. * logz * z) / 27. + (28528. * zeta3) / 243. + (29. * z * zeta3) / 6. -
                                      ((3. * sqrt3 * log3 + 30. * L_1 * sqrt3 + 10. * sqrt3 * (2. + 9. * z)) * (M_PI)) / 2187. +
                                      ((-44209. + (10. * complex::i()) * sqrt3 - 1728. * L_2 - 163368. * log2 + 40248. * sqrt5 - 63927. * z -
                                        1728. * L_2 * z + 98604. * log2 * z - 5184. * logz * z - 18576. * sqrt5 * z - 17496. * L_1 * (1. + z) +
                                        37152. * (log12sqrt52 + log2)) *
                                       (M_PI2)) /
                                          8748. +
                                      (449. * (M_PI4)) / 1215. -
                                      (27529. * z * (M_PI4)) / 29160.)
                                         .real();
    cache_ps[index_p(cu, 1, 2, 2)] = (-1336127. / 2187. - (28733. * L_1) / 729. + (44. * L_12) / 81. - (2176. * L_2) / 243. - (1856. * L_1 * L_2) / 81. +
                                      (208. * L_22) / 9. - (416. * Cl2PI3 * sqrt3) / 243. + ((2. * complex::i()) / 81.) * log3 * log3 * sqrt3 +
                                      (8773. * sqrtz) / 100. + ((16. * -1.) / 81.) * sqrt3 * t_2 + (660495139. * z) / 4860000. -
                                      (11734. * L_1 * z) / 81. - (2624. * L_2 * z) / 27. - (2624. * logz * z) / 9. + (13934. * zeta3) / 81. -
                                      122. * z * zeta3 + (4. * sqrt3 * (20. + 30. * L_1 + 3. * log3 + 90. * z) * (M_PI)) / 729. +
                                      ((39995. + 8640. * log12sqrt52 - 29232. * log2 - (20. * complex::i()) * sqrt3 + 9360. * sqrt5 + 29142. * z +
                                        5616. * log2 * z + 10368. * logz * z - 4320. * sqrt5 * z + 23328. * L_1 * (1. + z) +
                                        3456. * L_2 * (1. + z)) *
                                       (M_PI2)) /
                                          1458. +
                                      (226. * (M_PI4)) / 405. + (2846. * z * (M_PI4)) / 1215.)
                                         .real();
    cache_ps[index_p(cu, 2, 2, 2)] = (27476329. / 58320. + (40370. * L_1) / 243. + (604. * L_12) / 27. + (6928. * L_2) / 81. +
                                      (1064. * L_1 * L_2) / 27. - (52. * L_22) / 3. + (416. * Cl2PI3 * sqrt3) / 81. -
                                      ((2. * complex::i()) / 27.) * log3 * log3 * sqrt3 - (26319. * sqrtz) / 250. - ((16. * -1.) / 27.) * sqrt3 * t_2 +
                                      (151070269. * z) / 81000. + (9466. * L_1 * z) / 27. + (464. * L_2 * z) / 9. + (464. * logz * z) / 3. -
                                      (4388. * zeta3) / 27. - 300. * z * zeta3 - (4. * sqrt3 * (20. + 30. * L_1 + 3. * log3 + 90. * z) * (M_PI)) / 243. - ((41279. - 1728. * log12sqrt52 + 11808. * log2 - (20. * complex::i()) * sqrt3 - 1872. * sqrt5 + 72342. * z - 7344. * log2 * z + 10368. * logz * z + 864. * sqrt5 * z + 11664. * L_1 * (1. + z) + 3456. * L_2 * (1. + z)) * (M_PI2)) / 486. + (398. * (M_PI4)) / 135. + (7991. * z * (M_PI4)) / 810.)
                                         .real();

    // Corrections from LO and NLO coefficients (incl. pengFlag)
    // resummation logz in NLO coefficients
    cache_p[index_p(cu, 1, 1, 2)] += ((logz * (33433. + 12132. * L_1 + 2376. * L_2 - 12. * M_PI2) - 22. * (151. + 9. * log2z + 12. * M_PI2)) * z) / 162.;
    cache_p[index_p(cu, 1, 2, 2)] += (-2. * (302. + 18. * log2z + logz * (1663. + 576. * L_1 - 216. * L_2 - 12. * M_PI2) + 24. * M_PI2) * z) / 27.;
    cache_p[index_p(cu, 2, 2, 2)] += (2. * (-151. - 9. * log2z + 2. * logz * (296. + 99. * L_1 + 54. * L_2 - 6. * M_PI2) - 12. * M_PI2) * z) / 9.;
    cache_ps[index_p(cu, 1, 1, 2)] += (4. * logz * (3473. + 12. * M_PI2) * z) / 81.;
    cache_ps[index_p(cu, 1, 2, 2)] += (-64. * logz * (-124. + 3. * M_PI2) * z) / 27.;
    cache_ps[index_p(cu, 2, 2, 2)] += (16. * logz * (-91. + 12. * M_PI2) * z) / 9.;

    // resummation affecting arguments of logs and Dilogs
    cache_p[index_p(cu, 1, 1, 2)] += 4762. / 81. * logz * z;
    cache_p[index_p(cu, 1, 2, 2)] += 1688. / 27. * logz * z;
    cache_p[index_p(cu, 2, 2, 2)] += 904. / 9. * logz * z;
    cache_ps[index_p(cu, 1, 1, 2)] += 16. / 81. * logz * z;
    cache_ps[index_p(cu, 1, 2, 2)] += -64. / 27. * logz * z;
    cache_ps[index_p(cu, 2, 2, 2)] += 64. / 9. * logz * z;

    // Transforming mb in NLO coefficients from Pole scheme to MSbar
    cache_p[index_p(cu, 1, 1, 2)] += (4. * (4. + 3. * L_b) * (-392. + 675. * z)) / 243.;
    cache_p[index_p(cu, 1, 2, 2)] += (-44. * (4. + 3. * L_b) * (-19. + 54. * z)) / 81.;
    cache_p[index_p(cu, 2, 2, 2)] += (-8. * (4. + 3. * L_b) * (-8. + 27. * z)) / 27.;
    cache_ps[index_p(cu, 1, 1, 2)] += (1264. * (4. + 3. * L_b)) / 243.;
    cache_ps[index_p(cu, 1, 2, 2)] += (416. * (4. + 3. * L_b)) / 81.;
    cache_ps[index_p(cu, 2, 2, 2)] += (-704. * (4. + 3. * L_b)) / 27.;

    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2))
    {
        // Gerlach thesis eq. (6.28)
        cache_p[index_p(qq, 1, 8, 2)] = 208. / 81. * L_1 - L_2 / 27. + (2615. / 54. - 10. / 9. * M_PI2) * z - 5. / 9. * M_PI2 + 25. / (54. * sqrt3) * M_PI - 115. / 486.;
        cache_p[index_p(qq, 2, 8, 2)] = -11. / 27. * L_1 + 2. / 9. * L_2 + (20. / 3. * M_PI2 - 833. / 9.) * z + 10. / 3. * M_PI2 - 25. / (9. * sqrt3) * M_PI - 3125. / 81.;
        cache_ps[index_p(qq, 1, 8, 2)] = 448. / 81. * L_1 + 32. / 27. * L_2 + (1192. / 27. - 16. / 9. * M_PI2) * z - 8. / 9. * M_PI2 + 20. / (27. * sqrt3) * M_PI + 3580. / 243.;
        cache_ps[index_p(qq, 2, 8, 2)] = -248. / 27. * L_1 - 64. / 9. * L_2 + (32. / 3. * M_PI2 - 1088. / 9.) * z + 16. / 3. * M_PI2 - 40. / (9. * sqrt3) * M_PI - 4568. / 81.;

        // equation (6.30)
        cache_p[index_p(qq, 3, 8, 2)] = -85. / 27. * L_1 - 448. / 9. * L_2 - 196. / 3. * z + 25. / 6. * M_PI2 - 107. / (18. * sqrt3) * M_PI - 17201. / 81.;
        cache_p[index_p(qq, 4, 8, 2)] = -3269. / 162. * L_1 - 427. / 27. * L_2 + (20. / 3. * M_PI2 - 404. / 3.) * z + 169. / 12. * M_PI2 - 514. / (27. * sqrt3) * M_PI - 43016. / 243.;
        cache_p[index_p(qq, 5, 8, 2)] = 5120. / 27. * L_1 - 7168. / 9. * L_2 - 760. / 3. * z + 770. / 9. * M_PI2 - 28. / (9. * sqrt3) * M_PI - 430238. / 81.;
        cache_p[index_p(qq, 6, 8, 2)] = -8962. / 81. * L_1 - 6976. / 27. * L_2 + (200. / 3. * M_PI2 - 4222. / 3.) * z + 3761. / 27. * M_PI2 - 3220. / (27. * sqrt3) * M_PI - 474656. / 243.;
        cache_ps[index_p(qq, 3, 8, 2)] = 440. / 27. * L_1 + 512. / 9. * L_2 + 608. / 3. * z - 596. / 27. * M_PI2 - 52. / (9. * sqrt3) * M_PI + 22504. / 81.;
        cache_ps[index_p(qq, 4, 8, 2)] = -4804. / 81. * L_1 - 160. / 27. * L_2 + (32. / 3. * M_PI2 - 128. / 3.) * z + 1090. / 81. * M_PI2 - 1120. / (27. * sqrt3) * M_PI - 46988. / 243.;
        cache_ps[index_p(qq, 5, 8, 2)] = 17408. / 27. * L_1 + 8192. / 9. * L_2 + 11456. / 3. * z - 8912. / 27. * M_PI2 - 1984. / (9. * sqrt3) * M_PI + 420304. / 81.;
        cache_ps[index_p(qq, 6, 8, 2)] = -28624. / 81. * L_1 + 2048. / 27. * L_2 + (320. / 3. * M_PI2 - 160. / 3.) * z + 5608. / 81. * M_PI2 - 8416. / (27. * sqrt3) * M_PI - 423440. / 243.;

        // Gerlach thesis eq. (6.31) and (6.29)
        z = 0.;
    }
    z = cache_z;

    // Gerlach thesis eq. (6.29)
    cache_p[index_p(uu, 1, 8, 2)] += -10. / 9. * z;
    cache_p[index_p(uu, 2, 8, 2)] += 20. / 3. * z;
    cache_ps[index_p(uu, 1, 8, 2)] += -16. / 9. * z;
    cache_ps[index_p(uu, 1, 8, 2)] += 32. / 3. * z;

    // Gerlach thesis eq. (6.31)
    cache_p[index_p(uu, 3, 8, 2)] += -196. / 3. * z;
    cache_p[index_p(uu, 4, 8, 2)] += (-404. / 3. + 20. / 3. * M_PI2) * z;
    cache_p[index_p(uu, 5, 8, 2)] += -760. / 3. * z;
    cache_p[index_p(uu, 6, 8, 2)] += (-4222. / 3. + 200. / 3. * M_PI2) * z;
    cache_ps[index_p(uu, 3, 8, 2)] += 608. / 3. * z;
    cache_ps[index_p(uu, 4, 8, 2)] += (-128. / 3. + 32. / 3. * M_PI2) * z;
    cache_ps[index_p(uu, 5, 8, 2)] += 11456. / 3. * z;
    cache_ps[index_p(uu, 6, 8, 2)] += (-160. / 3. + 320. / 3. * M_PI2) * z;

    // as in Gerlach thesis eq. (6.20)
    for (int i = 1; i <= 6; i++)
    {
        cache_p[index_p(cu, i, 8, 2)] =
            0.5 * (cache_p[index_p(cc, i, 8, 2)] + cache_p[index_p(uu, i, 8, 2)]);
        cache_ps[index_p(cu, i, 8, 2)] =
            0.5 * (cache_ps[index_p(cc, i, 8, 2)] + cache_ps[index_p(uu, i, 8, 2)]);
    }

    // Gerlach thesis eq. (6.32)
    cache_p[index_p(cc, 8, 8, 2)] = -13. / 18.;
    cache_p[index_p(cu, 8, 8, 2)] = -13. / 18.;
    cache_p[index_p(uu, 8, 8, 2)] = -13. / 18.;
    cache_ps[index_p(cc, 8, 8, 2)] = -68. / 9.;
    cache_ps[index_p(cu, 8, 8, 2)] = -68. / 9.;
    cache_ps[index_p(uu, 8, 8, 2)] = -68. / 9.;

    // compute flag_LOz terms
    for (int i = 0; i < 576; i++)
    {
        cache_p_LO[i] = cache_p[i];
        cache_ps_LO[i] = cache_ps[i];
    }

    cache_p_LO[index_p(cc, 1, 1, 0)] = 23. / 72. - 11. / 6. * z;
    cache_p_LO[index_p(cc, 1, 2, 0)] = 1. / 6. - 2. * z;
    cache_p_LO[index_p(cc, 2, 2, 0)] = 1. - 3. * z;
    cache_ps_LO[index_p(cc, 1, 1, 0)] = -5. / 9.;
    cache_ps_LO[index_p(cc, 1, 2, 0)] = -4. / 3.;
    cache_ps_LO[index_p(cc, 2, 2, 0)] = 1.;

    cache_p_LO[index_p(cc, 1, 5, 0)] = 64. / 3. - 96. * z;
    cache_p_LO[index_p(cc, 1, 6, 0)] = 4. * z - 20. / 9.;
    cache_p_LO[index_p(cc, 2, 5, 0)] = 16. - 72. * z;
    cache_p_LO[index_p(cc, 2, 6, 0)] = 40. / 3. - 24. * z;

    cache_p_LO[index_p(cc, 5, 5, 0)] = -1296. * n_v * z + 408. * (n_l + n_v) + 512.;
    cache_p_LO[index_p(cc, 6, 6, 0)] = -72. * n_v * z + 170. / 3. * (n_l + n_v) + 416. / 9.;

    for (int i = 3; i <= 6; i++)
    {
        for (int j = i; j <= 6; j++)
        {
            cache_p_LO[index_p(uu, i, j, 0)] = cache_p_LO[index_p(cc, i, j, 0)];
            cache_ps_LO[index_p(uu, i, j, 0)] = cache_ps_LO[index_p(cc, i, j, 0)];
        }
    }

    // Gerlach thesis eq. (6.10)
    for (int i = 1; i <= 6; i++)
    {
        for (int j = i; j <= 6; j++)
        {
            cache_p_LO[index_p(cu, i, j, 0)] =
                0.5 * (cache_p_LO[index_p(cc, i, j, 0)] + cache_p_LO[index_p(uu, i, j, 0)]);
            cache_ps_LO[index_p(cu, i, j, 0)] =
                0.5 * (cache_ps_LO[index_p(cc, i, j, 0)] + cache_ps_LO[index_p(uu, i, j, 0)]);
        }
    }

    cache_p_LO[index_p(cc, 1, 1, 1)] = (1196. + 342. * L_1 + 447. * L_2 - 15. * (M_PI2)) / 324. + (z * (-4113. - 1008. * L_1 - 792. * L_2 - 3168. * logz + 4. * (M_PI2))) / 216.;
    cache_p_LO[index_p(cc, 1, 2, 1)] = (z * (1179. + 468. * L_1 - 72. * L_2 - 288. * logz - 4. * (M_PI2))) / 18. + (-452. - (333. * L_1) / 2. + 57. * L_2 + 15. * (M_PI2)) / 27.;
    cache_p_LO[index_p(cc, 2, 2, 1)] = (37. - 18. * L_1 + 12. * L_2 - 30. * (M_PI2)) / 18. + (z * (135. + 72. * L_1 - 36. * L_2 - 144. * logz + 4. * (M_PI2))) / 6.;
    cache_ps_LO[index_p(cc, 1, 1, 1)] = z * (-385. / 9. - (4. * (M_PI2)) / 27.) - (2. * (-1. + 18. * L_1 + 60. * L_2 + 3. * (M_PI2))) / 81.;
    cache_ps_LO[index_p(cc, 1, 2, 1)] = (16. * z * (-42. + (M_PI2))) / 9. + (8. * (11. + (9. * L_1) / 2. - 12. * L_2 + 3. * (M_PI2))) / 27.;
    cache_ps_LO[index_p(cc, 2, 2, 1)] = z * (44. - (16. * (M_PI2)) / 3.) - (8. * (-31. - 9. * L_1 - 3. * L_2 + 3. * (M_PI2))) / 9.;

    cache_p_LO[index_p(cc, 1, 1, 1)] += -5. / 486. - (5. * L_1) / 324. - (5. * z) / 54.;
    cache_p_LO[index_p(cc, 1, 2, 1)] += 10. / 81. + (5. * L_1) / 27. + (10. * z) / 9.;
    cache_p_LO[index_p(cc, 2, 2, 1)] += -10. / 27. - (5. * L_1) / 9. - (10. * z) / 3.;
    cache_ps_LO[index_p(cc, 1, 1, 1)] += -4. / 243. - (2. * L_1) / 81. - (4. * z) / 27.;
    cache_ps_LO[index_p(cc, 1, 2, 1)] += 16. / 81. + (8. * L_1) / 27. + (16. * z) / 9.;
    cache_ps_LO[index_p(cc, 2, 2, 1)] += -16. / 27. - (8. * L_1) / 9. - (16. * z) / 3.;

    cache_p_LO[index_p(cc, 1, 1, 1)] += -11. / 6. * z_replace;
    cache_p_LO[index_p(cc, 1, 2, 1)] += -2. * z_replace;
    cache_p_LO[index_p(cc, 2, 2, 1)] += -3. * z_replace;

    // Gerlach thesis eq. (6.13) and (6.16)
    for (int i = 1; i <= 2; i++)
    {
        for (int j = i; j <= 2; j++)
        {
            cache_p_LO[index_p(cu, i, j, 1)] =
                0.5 * (cache_p_LO[index_p(cc, i, j, 1)] + cache_p_LO[index_p(uu, i, j, 1)]);
            cache_ps_LO[index_p(cu, i, j, 1)] =
                0.5 * (cache_ps_LO[index_p(cc, i, j, 1)] + cache_ps_LO[index_p(uu, i, j, 1)]);
        }
    }

    // Gerlach thesis eq. (6.19)
    cache_p_LO[index_p(cc, 1, 8, 1)] = 5. / 18.;
    cache_p_LO[index_p(cc, 2, 8, 1)] = -5. / 3.;
    cache_ps_LO[index_p(cc, 1, 8, 1)] = 4. / 9.;
    cache_ps_LO[index_p(cc, 2, 8, 1)] = -8. / 3.;

    // Gerlach thesis eq. (6.21)
    cache_p_LO[index_p(cc, 3, 8, 1)] = -32. / 3.;
    cache_p_LO[index_p(cc, 4, 8, 1)] = -169. / 18.;
    cache_p_LO[index_p(cc, 5, 8, 1)] = -512. / 3.;
    cache_p_LO[index_p(cc, 6, 8, 1)] = -992. / 9.;
    cache_ps_LO[index_p(cc, 3, 8, 1)] = 64. / 3.;
    cache_ps_LO[index_p(cc, 4, 8, 1)] = -20. / 9.;
    cache_ps_LO[index_p(cc, 5, 8, 1)] = 1024. / 3.;
    cache_ps_LO[index_p(cc, 6, 8, 1)] = 256. / 9.;

    // Gerlach thesis eq. (6.20)
    for (int i = 1; i <= 6; i++)
    {
        cache_p_LO[index_p(cu, i, 8, 1)] = cache_p_LO[index_p(cc, i, 8, 1)];
        cache_p_LO[index_p(uu, i, 8, 1)] = cache_p_LO[index_p(cc, i, 8, 1)];
        cache_ps_LO[index_p(cu, i, 8, 1)] = cache_ps_LO[index_p(cc, i, 8, 1)];
        cache_ps_LO[index_p(uu, i, 8, 1)] = cache_ps_LO[index_p(cc, i, 8, 1)];
    }

    logz = cache_logz;
    //lastInput_compute_pp_s[0] = z;
    //lastInput_compute_pp_s[1] = mu_1;
    //lastInput_compute_pp_s[2] = mu_2;
    //lastInput_compute_pp_s[3] = mu_b;
    return;
}

void AmpDB2::poletoMSbar_pp_s()
{
    // arxiv:2106.05979 eq. (33)
    double log_mub_Mb = 2. * log(mu_b / Mb);
    PoletoMS_as1 = 4. / 3. + log_mub_Mb;
    PoletoMS_as2_z0 = (307. / 32. + M_PI2 / 3. + M_PI2 / 9. * log2 - 1. / 6. * zeta3 - (71. / 144. + M_PI2 / 18.) * nl + log_mub_Mb * (493. / 72. - nl * 13. / 36. + log_mub_Mb * (43. / 24. - nl / 12.)));
    PoletoMS_as2_z1 = 4. / 3. * M_PI2 / 8. * Mc / Mb - z;
    // a_s(mu_b) to a_s(mu_1)
    double beta0 = 23. / 3.;
    double log_mu1_mub = 2. * log(mu_1 / mu_b);
    double as_running1 = beta0 * log_mu1_mub;
    bool flag_LOz = true;
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1))
    {
        for (int i = 1; i <= 6; i++)
        {
            // not all terms used for n=2
            for (int j = i; j <= 8; j++)
            {
                if (j == 3)
                    j = 8;
                if (i >= 3)
                    j = 8;
                cache_p[index_p(qq, i, j, 2)] += 8. * PoletoMS_as1 * p(qq, i, j, 1, flag_LOz) + (32. * PoletoMS_as2_z0 + 16. * PoletoMS_as1 * PoletoMS_as1 + 8. * PoletoMS_as1 * as_running1) * p(qq, i, j, 0, flag_LOz) + 16. * 2. * PoletoMS_as2_z1 * p(uu, i, j, 0, flag_LOz);
                cache_ps[index_p(qq, i, j, 2)] += 8. * PoletoMS_as1 * p_s(qq, i, j, 1, flag_LOz) + (32. * PoletoMS_as2_z0 + 16. * PoletoMS_as1 * PoletoMS_as1 + 8. * PoletoMS_as1 * as_running1) * p_s(qq, i, j, 0, flag_LOz) + 16. * 2. * PoletoMS_as2_z1 * p_s(uu, i, j, 0, flag_LOz);
                cache_p_LO[index_p(qq, i, j, 2)] = cache_p[index_p(qq, i, j, 2)];
                cache_ps_LO[index_p(qq, i, j, 2)] = cache_ps[index_p(qq, i, j, 2)];
            }
            for (int j = i; j <= 8; j++)
            {
                cache_p_LO[index_p(qq, i, j, 1)] += 8. * PoletoMS_as1 * p(qq, i, j, 0, flag_LOz);
                cache_ps_LO[index_p(qq, i, j, 1)] += 8. * PoletoMS_as1 * p_s(qq, i, j, 0, flag_LOz);
                if (j == 1 or j == 2 or j == 8)
                {
                    cache_p[index_p(qq, i, j, 1)] += 8. * PoletoMS_as1 * p(qq, i, j, 0);
                    cache_ps[index_p(qq, i, j, 1)] += 8. * PoletoMS_as1 * p_s(qq, i, j, 0);
                }
                else if (3 <= j and j <= 6)
                {
                    cache_p[index_p(qq, i, j, 1)] += 8. * PoletoMS_as1 * p(qq, i, j, 0, flag_LOz);
                    cache_ps[index_p(qq, i, j, 1)] += 8. * PoletoMS_as1 * p_s(qq, i, j, 0, flag_LOz);
                }
            }
        }
    }
    return;
}

void AmpDB2::poletoPS_pp_s()
{
    // analog to arxiv:2106.05979 eq. (33) for PS mass
    double log_mu1_Mb = 2. * log(mu_1 / Mb);
    double log_mub_Mb = 2. * log(mu_b / Mb);
    // hep-ph/9804241v2 eq. (21)
    PoletoPS_as1 = 4. / 3. * mu_f / Mb;
    PoletoPS_as2 = 4. / 3. * mu_f / (Mb * 4.) * (a1 - b0 * (2. * log(mu_f / Mb) - 2.));
    bool flag_LOz = true;
    double beta0 = 23. / 3.;
    double as_running1 = beta0 * log_mu1_Mb; // a_s(Mb) to a_s(mu_1)
    PoletoMS_as1 = 4. / 3. + log_mub_Mb;
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1))
    {
        for (int i = 1; i <= 6; i++)
        {
            // not all terms used for n=2
            for (int j = i; j <= 8; j++)
            {
                if (j == 3)
                    j = 8;
                if (i >= 3)
                    j = 8;
                cache_p[index_p(qq, i, j, 2)] += 8. * PoletoPS_as1 * p(qq, i, j, 1, flag_LOz) + (32. * PoletoPS_as2 + 16. * PoletoPS_as1 * PoletoPS_as1 +
                                                                                                 8. * PoletoPS_as1 * as_running1 -
                                                                                                 8. * PoletoPS_as1 * 4. * (-PoletoPS_as1 + PoletoMS_as1) // internal Mb_PS to Mb(mu_b)
                                                                                                 ) * p(qq, i, j, 0, flag_LOz);
                cache_ps[index_p(qq, i, j, 2)] += 8. * PoletoPS_as1 * p_s(qq, i, j, 1, flag_LOz) + (32. * PoletoPS_as2 + 16. * PoletoPS_as1 * PoletoPS_as1 +
                                                                                                    8. * PoletoPS_as1 * as_running1 -
                                                                                                    8. * PoletoPS_as1 * 4. * (-PoletoPS_as1 + PoletoMS_as1)) *
                                                                                                       p_s(qq, i, j, 0, flag_LOz);
                cache_p_LO[index_p(qq, i, j, 2)] = cache_p[index_p(qq, i, j, 2)];
                cache_ps_LO[index_p(qq, i, j, 2)] = cache_ps[index_p(qq, i, j, 2)];
            }
            for (int j = i; j <= 8; j++)
            {
                cache_p_LO[index_p(qq, i, j, 1)] += 8. * PoletoPS_as1 * p(qq, i, j, 0, flag_LOz);
                cache_ps_LO[index_p(qq, i, j, 1)] += 8. * PoletoPS_as1 * p_s(qq, i, j, 0, flag_LOz);
                if (j == 1 or j == 2 or j == 8)
                {
                    cache_p[index_p(qq, i, j, 1)] += 8. * PoletoPS_as1 * p(qq, i, j, 0);
                    cache_ps[index_p(qq, i, j, 1)] += 8. * PoletoPS_as1 * p_s(qq, i, j, 0);
                }
                else if (3 <= j and j <= 6)
                {
                    cache_p[index_p(qq, i, j, 1)] += 8. * PoletoPS_as1 * p(qq, i, j, 0, flag_LOz);
                    cache_ps[index_p(qq, i, j, 1)] += 8. * PoletoPS_as1 * p_s(qq, i, j, 0, flag_LOz);
                }
            }
        }
    }
    return;
}

void AmpDB2::poletoMSbar_pp_s_partialN3LO()
{
    // arxiv:2106.05979 eq. (33)
    double log_mub_Mb = 2. * log(mu_b / Mb);
    double log_mub_Mb_2 = log_mub_Mb * log_mub_Mb;
    double log_mub_Mb_3 = log_mub_Mb_2 * log_mub_Mb;
    double nl2 = nl * nl;
    double Dilogsqrtz = gslpp_special_functions::dilog(sqrtz);
    double Dilogminussqrtz = gslpp_special_functions::dilog(-sqrtz);
    double oneminussqrtz2 = (1. - sqrtz) * (1. - sqrtz);
    double log1minussqrtz = log(1. - sqrtz);
    double log1plussqrtz = log(1. + sqrtz);
    double sqrtz3 = sqrtz * sqrtz * sqrtz;
    double log_mub_Mc = 2. * log(mu_b / Mc);
    PoletoMS_as1 = 4. / 3. + log_mub_Mb;
    PoletoMS_as2 = (307. / 32. + M_PI2 / 3. + M_PI2 / 9. * log2 - 1. / 6. * zeta3 - (71. / 144. + M_PI2 / 18.) * nl + log_mub_Mb * (493. / 72. - nl * 13. / 36. + log_mub_Mb * (43. / 24. - nl / 12.))) + 4. / 3. * M_PI2 / 8. * sqrtz - z;
    PoletoMS_as3 = (8481925. / 93312. + (2. * Dilogminussqrtz) / 9. + (2. * Dilogsqrtz) / 9. - (55. * log2_4) / 162. + (652841. * M_PI2) / 38880. - (575. * log2 * M_PI2) / 162. - (22. * log2_2 * M_PI2) / 81. - (695. * M_PI4) / 7776. + (137. * nl) / 216. - (M_PI2 * nl) / 27. +
                    (log_mub_Mb_3 * (1591. - 160. * nl + 4. * nl2)) / 432. + (log_mub_Mb_2 * (19315. - 2206. * nl + 52. * nl2)) / 864. - (220. * polylog4_12) / 27. + 22.968067156 * sqrtz + (2. * Dilogminussqrtz * sqrtz) / 9. - (2. * Dilogsqrtz * sqrtz) / 9. - 0.982666666 * nl * sqrtz +
                    1.022111111 * sqrtz3 + (2. * Dilogminussqrtz * sqrtz3) / 9. - (2. * Dilogsqrtz * sqrtz3) / 9. + (M_PI2 * sqrtz3) / 9. + 4.014666667 * z - 0.300333333 * nl * z - (log_mub_Mc * sqrtz * (3. * log1minussqrtz * logz - 3. * log1plussqrtz * logz - 3. * M_PI2 + 24. * sqrtz + 6. * logz * sqrtz + 6. * log2z * sqrtz3 - 12. * log1minussqrtz * logz * sqrtz3 - 12. * log1plussqrtz * logz * sqrtz3 + 4. * M_PI2 * sqrtz3 + 9. * log1minussqrtz * logz * z - 9. * log1plussqrtz * logz * z - 9. * M_PI2 * z - 6. * Dilogminussqrtz * (1. + 4. * sqrtz3 + 3. * z) + Dilogsqrtz * (6. - 24. * sqrtz3 + 18. * z))) / 18. + 0.0493333333 * zeta2 + (2. * Dilogminussqrtz * zeta2) / 9. + (2. * Dilogsqrtz * zeta2) / 9. - (log2z * zeta2) / 18. - (M_PI2 * zeta2) / 27. +
                    (log_mub_Mb * (12. * Dilogminussqrtz + 12. * Dilogsqrtz + 118.990000008 * sqrtz + 6. * Dilogminussqrtz * sqrtz - 6. * Dilogsqrtz * sqrtz + 3. * M_PI2 * sqrtz - 9.624000006 * nl * sqrtz + 13.037999994 * sqrtz3 - 6. * Dilogminussqrtz * sqrtz3 +
                                   6. * Dilogsqrtz * sqrtz3 - 3. * M_PI2 * sqrtz3 - 1.206 * nl * sqrtz3 - 38.372000004 * z + 3.96 * nl * z - 3. * logz * (-1. + z) * (log1minussqrtz * (2. - sqrtz + 2. * z) + log1plussqrtz * (2. + sqrtz + 2. * z)) - 12. * Dilogminussqrtz * zeta2 -
                                   12. * Dilogsqrtz * zeta2 + 3. * log2z * zeta2 + 2. * M_PI2 * zeta2)) /
                        18. +
                    (logz * (log1minussqrtz * oneminussqrtz2 * (1. + sqrtz + z) + 1. * (sqrtz * (-76.2645000015 + 4.9559999985 * nl - 13.544000001 * sqrtz + 0.1544999985 * z) +
                                                                                        log1plussqrtz * (1. + sqrtz + sqrtz3 + zeta2)))) /
                        9. +
                    (nl * (-246643. + 288. * log2_4 + 36. * (-967. - 88. * log2 + 16. * log2_2) * M_PI2 + 732. * M_PI4 + 6912. * polylog4_12 - 78084. * zeta3)) / 23328. + (58. * zeta3) / 27. - (1439. * M_PI2 * zeta3) / 432. +
                    (nl2 * (2353. + 936. * M_PI2 + 3024. * zeta3)) / 23328. + (log_mub_Mb * (177305. + 10656. * (3. + log2) * zeta2 + 4. * nl2 * (89. + 72. * zeta2) - 4824. * zeta3 - 2. * nl * (72. * (49. + 4. * log2) * zeta2 + 7. * (1447. + 144. * zeta3)))) / 2592. + (1975. * zeta5) / 216.) /
                   M_PI3;
    // a_s(mu_b) to a_s(mu_1)
    double beta0 = 23. / 3.;
    double beta1 = 51. / 8. - 19. / 24. * 5.;
    double log_mu1_mub = 2. * log(mu_1 / mu_b);
    double as_running1 = beta0 * log_mu1_mub;
    double as_running2 = -beta0 * beta0 * log_mu1_mub * log_mu1_mub + beta1 * log_mu1_mub;

    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1))
    {
        for (int i = 1; i <= 8; i++)
        {
            if (i == 7)
                i++;
            for (int j = i; j <= 8; j++)
            {
                if (j == 7)
                    j++;
                if (j < i)
                    continue;
                cache_p[index_p(qq, i, j, 3)] = 8. * PoletoMS_as1 * p(qq, i, j, 2) + (32. * PoletoMS_as2 + 16. * PoletoMS_as1 * PoletoMS_as1 + 8. * PoletoMS_as1 * as_running1) * p(qq, i, j, 1) + (128. * (PoletoMS_as3 + PoletoMS_as1 * PoletoMS_as2) + 8. * PoletoMS_as1 * as_running2 + 32. * PoletoMS_as2 * as_running1) * p(qq, i, j, 0);
                cache_ps[index_p(qq, i, j, 3)] = 8. * PoletoMS_as1 * p_s(qq, i, j, 2) + (32. * PoletoMS_as2 + 16. * PoletoMS_as1 * PoletoMS_as1 + 8. * PoletoMS_as1 * as_running1) * p_s(qq, i, j, 1) + (128. * (PoletoMS_as3 + PoletoMS_as1 * PoletoMS_as2) + 8. * PoletoMS_as1 * as_running2 + 32. * PoletoMS_as2 * as_running1) * p_s(qq, i, j, 0);
                if (j >= 3 and j <= 6)
                {
                    cache_p[index_p(qq, i, j, 2)] = 0.;
                    cache_ps[index_p(qq, i, j, 2)] = 0.;
                }
                cache_p[index_p(qq, i, j, 2)] += 8. * PoletoMS_as1 * p(qq, i, j, 1) + (32. * PoletoMS_as2 + 16. * PoletoMS_as1 * PoletoMS_as1 +
                                                                                       8. * PoletoMS_as1 * as_running1) *
                                                                                          p(qq, i, j, 0);
                cache_ps[index_p(qq, i, j, 2)] += 8. * PoletoMS_as1 * p_s(qq, i, j, 1) + (32. * PoletoMS_as2 + 16. * PoletoMS_as1 * PoletoMS_as1 +
                                                                                          8. * PoletoMS_as1 * as_running1) *
                                                                                             p_s(qq, i, j, 0);
                cache_p_LO[index_p(qq, i, j, 2)] = cache_p[index_p(qq, i, j, 2)];
                cache_ps_LO[index_p(qq, i, j, 2)] = cache_ps[index_p(qq, i, j, 2)];

                cache_p[index_p(qq, i, j, 1)] += 8. * PoletoMS_as1 * p(qq, i, j, 0);
                cache_ps[index_p(qq, i, j, 1)] += 8. * PoletoMS_as1 * p_s(qq, i, j, 0);
                cache_p_LO[index_p(qq, i, j, 1)] = cache_p[index_p(qq, i, j, 1)];
                cache_ps_LO[index_p(qq, i, j, 1)] = cache_ps[index_p(qq, i, j, 1)];
                cache_p_LO[index_p(qq, i, j, 0)] = cache_p[index_p(qq, i, j, 0)];
                cache_ps_LO[index_p(qq, i, j, 0)] = cache_ps[index_p(qq, i, j, 0)];
            }
        }
    }
    return;
}

void AmpDB2::poletoPS_pp_s_partialN3LO()
{
    // analog to arxiv:2106.05979 eq. (33) for PS mass
    double log_mu1_Mb = 2. * log(mu_1 / Mb);
    double log_mub_Mb = 2. * log(mu_b / Mb);
    // hep-ph/9804241v2 eq. (21)
    double log_muf_Mb = 2. * log(mu_f / Mb);
    PoletoPS_as1 = 4. / 3. * mu_f / Mb;
    double PoletoPS_as1_2 = PoletoPS_as1 * PoletoPS_as1;
    PoletoPS_as2 = 4. / 3. * mu_f / (Mb * 4.) * (a1 - b0 * (log_muf_Mb - 2.));
    PoletoPS_as3 = 4. / 3. * mu_f / (Mb * 16.) * (a2 - (2. * a1 * b0 + b1) * (log_muf_Mb - 2.) + b0_2 * (log_muf_Mb * log_muf_Mb - 4. * log_muf_Mb + 8.));
    PoletoMS_as1 = 4. / 3. + log_mub_Mb;
    PoletoMS_as2 = (307. / 32. + M_PI2 / 3. + M_PI2 / 9. * log2 - 1. / 6. * zeta3 - (71. / 144. + M_PI2 / 18.) * nl + log_mub_Mb * (493. / 72. - nl * 13. / 36. + log_mub_Mb * (43. / 24. - nl / 12.))) + 4. / 3. * M_PI2 / 8. * sqrtz - z;
    double beta0 = 23. / 3.;
    double beta1 = 51. / 8. - 19. / 24. * 5.;
    double as_running1 = beta0 * log_mu1_Mb;
    double as_running1_mub = beta0 * 2. * log(mu_1 / mu_b);
    double as_running2 = -beta0 * beta0 * log_mu1_Mb * log_mu1_Mb + beta1 * log_mu1_Mb;
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1))
    {
        for (int i = 1; i <= 8; i++)
        {
            if (i == 7)
                i++;
            for (int j = i; j <= 8; j++)
            {
                if (j == 7)
                    j++;
                if (j < i)
                    continue;
                cache_p[index_p(qq, i, j, 3)] = 8. * PoletoPS_as1 * p(qq, i, j, 2) + (32. * PoletoPS_as2 + 16. * PoletoPS_as1_2 + 8. * PoletoPS_as1 * as_running1 - 8. * PoletoPS_as1 * 4. * (-PoletoPS_as1 + PoletoMS_as1)) * p(qq, i, j, 1) + (128. * (PoletoPS_as3 + PoletoPS_as1 * PoletoPS_as2) - (-PoletoPS_as1 + PoletoMS_as1) * (128. * (PoletoPS_as1_2 + PoletoPS_as2) + 32. * PoletoPS_as1 * (as_running1 + as_running1_mub)) - (-PoletoPS_as2 + PoletoMS_as2) * 128. * PoletoPS_as1 + (32. * PoletoPS_as1_2 + 64. * PoletoPS_as2) * as_running1 + 8. * PoletoPS_as1 * as_running2) * p(qq, i, j, 0);
                cache_ps[index_p(qq, i, j, 3)] = 8. * PoletoPS_as1 * p_s(qq, i, j, 2) + (32. * PoletoPS_as2 + 16. * PoletoPS_as1_2 + 8. * PoletoPS_as1 * as_running1 - 8. * PoletoPS_as1 * 4. * (-PoletoPS_as1 + PoletoMS_as1)) * p_s(qq, i, j, 1) + (128. * (PoletoPS_as3 + PoletoPS_as1 * PoletoPS_as2) - (-PoletoPS_as1 + PoletoMS_as1) * (128. * (PoletoPS_as1_2 + PoletoPS_as2) + 32. * PoletoPS_as1 * (as_running1 + as_running1_mub)) - (-PoletoPS_as2 + PoletoMS_as2) * 128. * PoletoPS_as1 + (32. * PoletoPS_as1_2 + 64. * PoletoPS_as2) * as_running1 + 8. * PoletoPS_as1 * as_running2) * p_s(qq, i, j, 0);
                if (j >= 3 and j <= 6)
                {
                    cache_p[index_p(qq, i, j, 2)] = 0.;
                    cache_ps[index_p(qq, i, j, 2)] = 0.;
                }
                cache_p[index_p(qq, i, j, 2)] += 8. * PoletoPS_as1 * p(qq, i, j, 1) + (32. * PoletoPS_as2 + 16. * PoletoPS_as1_2 +
                                                                                       8. * PoletoPS_as1 * as_running1 -
                                                                                       8. * PoletoPS_as1 * 4. * (-PoletoPS_as1 + PoletoMS_as1)) *
                                                                                          p(qq, i, j, 0);
                cache_ps[index_p(qq, i, j, 2)] += 8. * PoletoPS_as1 * p_s(qq, i, j, 1) + (32. * PoletoPS_as2 + 16. * PoletoPS_as1_2 +
                                                                                          8. * PoletoPS_as1 * as_running1 -
                                                                                          8. * PoletoPS_as1 * 4. * (-PoletoPS_as1 + PoletoMS_as1)) *
                                                                                             p_s(qq, i, j, 0);
                cache_p_LO[index_p(qq, i, j, 2)] = cache_p[index_p(qq, i, j, 2)];
                cache_ps_LO[index_p(qq, i, j, 2)] = cache_ps[index_p(qq, i, j, 2)];

                cache_p[index_p(qq, i, j, 1)] += 8. * PoletoPS_as1 * p(qq, i, j, 0);
                cache_ps[index_p(qq, i, j, 1)] += 8. * PoletoPS_as1 * p_s(qq, i, j, 0);
                cache_p_LO[index_p(qq, i, j, 1)] = cache_p[index_p(qq, i, j, 1)];
                cache_ps_LO[index_p(qq, i, j, 1)] = cache_ps[index_p(qq, i, j, 1)];
                cache_p_LO[index_p(qq, i, j, 0)] = cache_p[index_p(qq, i, j, 0)];
                cache_ps_LO[index_p(qq, i, j, 0)] = cache_ps[index_p(qq, i, j, 0)];
            }
        }
    }
    return;
}

void AmpDB2::compute_partialNNLO()
{
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1))
    {
        for (int i = 1; i <= 6; i++)
        {
            // not all terms used for n=2
            for (int j = i; j <= 8; j++)
            {
                if (j == 3)
                    j = 8;
                if (i >= 3)
                    j = 8;
                cache_p[index_p(qq, i, j, 2)] = 0.;
                cache_ps[index_p(qq, i, j, 2)] = 0.;
            }
        }
    }
    for (int i = 0; i < 8; i++)
    {
        if (i == 6)
            i = 7;
        C_Misiak_NNLO[i] = 0.;
    }
    return;
}

gslpp::complex AmpDB2::PBd()
{
    double mbpole = mySM.Mbar2Mp(mySM.getQuarks(QCD::BOTTOM).getMass(), QCD::BOTTOM);
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
                         (gslpp::complex(1, 2. * mySM.getPhiBd(), true) * (n[0] + (n[5] * B2 + n[10]) / B1) - gslpp::complex(1. / mySM.getCKM().computeRt(), mySM.getCKM().computeBeta() + 2. * mySM.getPhiBd(), true) * (n[1] + (n[6] * B2 + n[11]) / B1) + gslpp::complex(1. / mySM.getCKM().computeRt() / mySM.getCKM().computeRt(), 2. * (mySM.getCKM().computeBeta() + mySM.getPhiBd()), true) * (n[2] + (n[7] * B2 + n[12]) / B1));

    return PBd;
}

gslpp::complex AmpDB2::PBs()
{
    double mbpole = mySM.Mbar2Mp(mySM.getQuarks(QCD::BOTTOM).getMass(), QCD::BOTTOM);

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
                         (gslpp::complex(1, 2. * mySM.getPhiBs(), true) * (n[0] + (n[5] * B2 + n[10]) / B1) - gslpp::complex(1. / mySM.getCKM().computeRts(), -mySM.getCKM().computeBetas() + 2. * mySM.getPhiBs(), true) * (n[1] + (n[6] * B2 + n[11]) / B1) + gslpp::complex(1. / mySM.getCKM().computeRts() / mySM.getCKM().computeRts(), 2. * (-mySM.getCKM().computeBetas() + mySM.getPhiBs()), true) * (n[2] + (n[7] * B2 + n[12]) / B1));

    return PBs;
}
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

AmpDB2::AmpDB2(const StandardModel& SM_i, bool flag_RI)
: mySM(SM_i), meMStoRI(5, 0.), coeffsMStoRI(3, 0.)
{
    mySM.initializeBParameter("BBs");
    mySM.initializeBParameter("BBd");
    mySM.initializeBParameter("BBs_subleading");
    mySM.initializeBParameter("BBd_subleading");    
    this->flag_resumz = true;    
    this->flag_RI = flag_RI; 
    
    //hep-ph/0606197 eq. 4.7 - 4.10
    double meMStoRI0[5] = {-3. - 5./3.+8.*log2, 0., 0., 0., 0.},
    meMStoRI1[5] = {0., 61./9.+44./9.*log2, -7./9.+28./9.*log2, 0., 0.},
    meMStoRI2[5] = {0., -25./9.+28./9.*log2, -29./9.+44./9.*log2, 0., 0.},
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
}

gslpp::complex AmpDB2::RBs(orders order)
{
    mySM.getFlavour().getHDF2().getCoeffBs().getOrder();
    if (mySM.getFlavour().getHDF2().getCoeffBs().getOrder() < getHighest(order))
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
    double Mb = mySM.Mrun(mySM.getBBs().getMu(),
            mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
            mySM.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
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
    if (mySM.getFlavour().getHDF2().getCoeffBd().getOrder() < getHighest(order))
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
    if (mySM.getFlavour().getHDF2().getCoeffBs().getOrder() < getHighest(order))
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed");

    //Wilson coefficients in same mass scale and scheme as B parameters
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffBs(
            mySM.getBBs().getMu(),
            mySM.getBBs().getScheme());

    gslpp::vector<double> me(mySM.getBBs().getBpars());
    double MBs = mySM.getMesons(QCD::B_S).getMass();
    double Mb = mySM.Mrun(mySM.getBBs().getMu(),
            mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
            mySM.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
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

 /*******************************************************************************
 *  @f$\Gamma_{21}@f$ in NLO from Ciuchini (hep-ph/0308029v2)                   * 
 * ****************************************************************************/  

gslpp::complex AmpDB2::Gamma21overM21_BdFULLNLO_tradBasis(){
    //source: hep-ph/0308029v2
    std::cout.precision(4);

    computeCKMandMasses(NLO);
    //calculate M_21 / <O_1>
    gslpp::vector<gslpp::complex> ** M21overme0_times_8MB = mySM.getFlavour().ComputeCoeffBd(
            mySM.getBBd().getMu(),
            mySM.getBBd().getScheme());
    gslpp::complex M21overme0 = ((*(M21overme0_times_8MB[LO]))(0) + (*(M21overme0_times_8MB[NLO]))(0)) / (8. * MB);
    
    //calculate DB=1 Wilson coefficients
    computeWilsonCoeffsDB1bsg(); 
    
    //calculate DB=2 coefficients for usage of "c(quark)"
    computeF0();
    computeF1();
    computeP();
    computeD();

    //calculate DB=2 matrix elements for usage of "me" and "delta_1overm_tradBasis(quark)"
    compute_matrixelements(d);
        
    //hep-ph/0308029v2: eq. 16 divided by M_21
    gslpp::complex Gamma21overM21_Bd = -Gf2 / (24 * M_PI * MB) / M21overme0 *
                (Mb2_prefactor * (c(d)(0) + c(d)(1) * me(1)/me(0) + c(d)(2) * me(2)/me(0)) + 
            Mb_PS * Mb_PS * delta_1overm_tradBasis(d)/me(0));
    return Gamma21overM21_Bd;
}

gslpp::complex AmpDB2::Gamma21overM21_BsFULLNLO_tradBasis(){
    //source: hep-ph/0308029v2
    std::cout.precision(4);

    computeCKMandMasses(NLO);
    
    //calculate M_21 / <O_1>
    gslpp::vector<gslpp::complex> ** M21overme0_times_8MB = mySM.getFlavour().ComputeCoeffBs(
            mySM.getBBs().getMu(),
            mySM.getBBs().getScheme());
    gslpp::complex M21overme0 = ((*(M21overme0_times_8MB[LO]))(0) + (*(M21overme0_times_8MB[NLO]))(0)) / (8. * MB_s);
    
    //calculate DB=1 Wilson coefficients
    computeWilsonCoeffsDB1bsg(); 
    
    //calculate DB=2 coefficients for usage of "c(quark)"
    computeF0();
    computeF1();
    computeP();
    computeD();
    
    //calculate DB=2 matrix elements for usage of "me" and "delta_1overm_tradBasis(quark)"
    compute_matrixelements(s);

    //hep-ph/0308029v2: eq. 16 divided by M_21
    gslpp::complex Gamma21overM21_Bs = -Gf2 / (24 * M_PI * MB_s) / M21overme0 *
                (Mb2_prefactor * (c(s)(0) + c(s)(1) * me(1)/me(0) + c(s)(2) * me(2)/me(0)) + 
            Mb_PS * Mb_PS * delta_1overm_tradBasis(s)/me(0));
    return Gamma21overM21_Bs;
}

gslpp::complex AmpDB2::Gamma21overM21_BsLO_tradBasis(){
    //source: hep-ph/0308029v2
    std::cout.precision(4);

    computeCKMandMasses(NLO);
    
    //calculate M_21 / <O_1>
    gslpp::vector<gslpp::complex> ** M21overme0_times_8MB = mySM.getFlavour().ComputeCoeffBs(
            mySM.getBBs().getMu(),
            mySM.getBBs().getScheme());
    gslpp::complex M21overme0 = ((*(M21overme0_times_8MB[LO]))(0) + (*(M21overme0_times_8MB[NLO]))(0)) / (8. * MB_s);
    
    //calculate DB=1 Wilson coefficients
    computeWilsonCoeffsDB1bsg(); 
    
    //calculate DB=2 coefficients for usage of "c(quark)"
    computeF0();
    computeF1();
    computeP();
    computeD_LO();
    
    //calculate DB=2 matrix elements for usage of "me" and "delta_1overm_tradBasis(quark)"
    compute_matrixelements(s);

    //hep-ph/0308029v2: eq. 16 divided by M_21
    gslpp::complex Gamma21overM21_Bs = -Gf2 / (24 * M_PI * MB_s) / M21overme0 *
                (Mb2_prefactor * (c(s)(0) + c(s)(1) * me(1)/me(0) + c(s)(2) * me(2)/me(0)) + 
            Mb_PS * Mb_PS * delta_1overm_tradBasis(s)/me(0));
    return Gamma21overM21_Bs;
}

void AmpDB2::computeCKMandMasses(orders order, mass_schemes mass_scheme) {
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
    
    //DB=1 matching scales (arxiv: 2205.07907 Results. or Gerlach thesis eq. 7.7) varied by "getMub()" or fixed to 4.2
    mu_1 = mySM.getMub();
    mu_b = mu_1;
    mu_1_overm = 4.2;       
    
    //MSbar bottom quark mass Mb(Mb)
    Mb_Mb = mySM.getQuarks(QCD::BOTTOM).getMass();
    
    //MSbar bottom quark mass Mb(mu_b)
    double Mb_mub = mySM.Mrun(mu_b,
        mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
        mySM.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    
    //MSbar charm quark mass Mc(Mc)
    Mc_Mc = mySM.getQuarks(QCD::CHARM).getMass();
    
    //MSbar charm quark mass Mb(mu_b)    
    double Mc_mub = mySM.Mrun(mu_b,
        mySM.getQuarks(QCD::CHARM).getMass_scale(),
        mySM.getQuarks(QCD::CHARM).getMass(), FULLNNLO);
    //MSbar charm quark mass Mb(Mb)    
    double Mc_Mb = mySM.Mrun(mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
        mySM.getQuarks(QCD::CHARM).getMass_scale(),
        mySM.getQuarks(QCD::CHARM).getMass(), FULLNNLO);
    std::cout << Mc_mub << "\n";
    //double Mc_pole = 1.67;
    
    //pole mass of bottom:quark
    Mb_pole = 4.757; //fixed to compare with (Gerlach thesis)
    //Mb_pole = mySM.Mbar2Mp(mySM.getQuarks(QCD::BOTTOM).getMass
    
    //PS mass of bottom quark
    Mb_PS = 4.479; //fixed to compare with (Gerlach thesis)
    //NNLO evaluation from hep-ph/9804241v2 eq. (25)
    //Mb_PS = Mb_Mb * (1. + 16./3. * as_4pi * (1. - mu_f/Mb_Mb) + 16 * as_4pi * as_4pi * (K - mu_f/(3. * Mb_Mb) * (a1 - b0 * (2. * log(mu_f/Mb_Mb) - 2))));    
    
    //hep-ph/9804241v2 eq. (21)
    PoletoPS_as1 = 4./3. * mu_f / Mb_pole;
    PoletoPS_as2 = 4./3. * mu_f / (Mb_pole * 4.) * (a1 - b0 * (2. * log(mu_f/Mb_Mb) - 2.));
    
    //DB=2 matching scale mu_2
    mu_2 = 4.757;//mySM.getBBs().getMu();
    
    //strong coupling constant divided by 4*Pi
    as_4pi_mu1 = mySM.Als(mu_1, FULLNNNLO, true)/(4.*M_PI);//mySM.Alstilde5(mu_1);
    as_4pi_mu2 = mySM.Als(mu_2, FULLNNNLO, true)/(4.*M_PI);;//mySM.Alstilde5(mu_2);
    as_4pi = mySM.Als(Mb_Mb, FULLNNNLO, true)/(4.*M_PI);;//mySM.Alstilde5(Mb_Mb);

    //adapt "Mb2_prefactor", "Mb" and "z" to the used mass scheme
    //explained in Gerlach thesis chapter 7.0 and arxiv:2205.07907 Results.
    if(order == NNLO){
        switch (mass_scheme) {
            case pole:
                Mb2_prefactor = Mb_pole * Mb_pole;
                Mb = Mb_pole;
                z = Mc_mub * Mc_mub / (Mb_mub * Mb_mub);
                this->flag_resumz = true;
                break;
            case MSbar:
                Mb2_prefactor = Mb_mub * Mb_mub;
                Mb = Mb_pole;
                z = Mc_mub * Mc_mub / (Mb_mub * Mb_mub);
                this->flag_resumz = true;            
                break;
            case PS:
                Mb2_prefactor = Mb_PS * Mb_PS;
                Mb = Mb_pole;
                //arxiv:2106.05979 eq. 31
                z = Mc_mub * Mc_mub / (Mb_mub * Mb_mub);// * (1. - 2. * (16./3. * as_4pi * (1. - mu_f/Mb_Mb) + 16 * as_4pi * as_4pi * (K - mu_f/(3. * Mb_Mb) * (a1 - b0 * (2. * log(mu_f/Mb_Mb) - 2)))));
                this->flag_resumz = true;            
                break;
            case only1overmb:
                Mb2_prefactor = 0.;
                Mb = Mb_pole;
                break;                
            default:
                throw(std::runtime_error("mass_scheme not implemented"));
        }
    }
    //adapt to MSbar mass scheme for the traditional basis
    if (order == NLO){
        Mb2_prefactor = Mb_mub * Mb_mub;
        Mb = Mb_mub;
        x_1 = mu_1/Mb;
        x_2 = mu_2/Mb;
        logx_1 = log(x_1);
        logx_2 = log(x_2);
        double Mc = Mc_mub;
        if (!flag_resumz) Mc = Mc_Mc;
        z = Mc * Mc / (Mb_mub * Mb_mub);
    }

    Gf2 = mySM.getGF() * mySM.getGF();
    Md = mySM.getQuarks(QCD::DOWN).getMass();
    Ms = mySM.Mrun(mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
                mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    MB = mySM.getMesons(QCD::B_D).getMass();
    MB_s = mySM.getMesons(QCD::B_S).getMass();

    z2 = z * z;
    logz = log(z);
    oneminusz2 = (1. - z) * (1. - z);
    sqrt1minus4z = sqrt(1. - 4. * z);
    
    //calculate "z" values for 1/m_b contributions
    z_1overm = Mc_Mc * Mc_Mc / (Mb_Mb * Mb_Mb); //hep-ph/0612167 eq. 27
    z_1overm = Mc_Mb * Mc_Mb / (Mb_mub * Mb_mub); //arxiv:1910.00970 eq. 11
    //if z for the 1/m_b contributions shall be varied with "mu_b"
    //z_1overm = Mc_mub * Mc_mub / (Mb_mub * Mb_mub);
    
    z_1overm2 = z_1overm * z_1overm;
    oneminusz_1overm2 = (1. - z_1overm) * (1. - z_1overm);
    sqrt1minus4z_1overm = sqrt(1. - 4. * z_1overm);
    
    //values needed only for the traditional basis
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

void AmpDB2::computeWilsonCoeffsDB1bsg(){
    //NNLO DB=1 Wilson coefficients
    gslpp::vector<gslpp::complex> ** WilsonCoeffsDB1bsg = mySM.getFlavour().ComputeCoeffsgamma_Buras(mu_1);
    for (int i = 0; i < 8; i++) {
        if (i==6) i=7;
        cacheC[i] = (*(WilsonCoeffsDB1bsg[LO]))(i) + (*(WilsonCoeffsDB1bsg[NLO]))(i) + (*(WilsonCoeffsDB1bsg[NNLO]))(i);
        cacheC_LO[i] = (*(WilsonCoeffsDB1bsg[LO]))(i);
        cacheC_NLO[i] = (*(WilsonCoeffsDB1bsg[NLO]))(i);
        cacheC_NNLO[i] = (*(WilsonCoeffsDB1bsg[NNLO]))(i);
    } 
//    for(int i=0; i<=7; i++){
//        if(i==6) i++;
//        std::cout << "C_" << i << " "
//                << cacheC[i].gslpp::complex::real() << " "
//                << cacheC_LO[i].gslpp::complex::real() << " "
//                << cacheC_NLO[i].gslpp::complex::real() << " "
//                << cacheC_NNLO[i].gslpp::complex::real() << "\n";
//    }
//    std::cout << "--------\n";
    
    //LO DB=1 Wilson coefficients for 1/mb corrections
    WilsonCoeffsDB1bsg = mySM.getFlavour().ComputeCoeffsgamma_Buras(mu_1_overm);    
    C_1LO = (*(WilsonCoeffsDB1bsg[LO]))(0);
    C_2LO = (*(WilsonCoeffsDB1bsg[LO]))(1);
    K_1 = 3. * C_1LO * C_1LO + 2. * C_1LO * C_2LO;
    K_2 = C_2LO * C_2LO;
    return;
}

void AmpDB2::computeWilsonCoeffs(){
    //NLO DB=1 Wilson coefficients C_i, i=1-6,8
   gslpp::vector<gslpp::complex> ** WilsonCoeffs = mySM.getFlavour().ComputeCoeffBMll_Buras(mu_1, QCD::NOLEPTON);
    for (int i = 0; i < 8; i++) {
        if (i==6) i=7;
        cacheC[i] = (*(WilsonCoeffs[LO]))(i) + (*(WilsonCoeffs[NLO]))(i);
        cacheC_LO[i] = (*(WilsonCoeffs[LO]))(i);
        cacheC_NLO[i] = (*(WilsonCoeffs[NLO]))(i);        
    }  
//    for(int i=0; i<=7; i++){
//        if(i==6) i++;
//        std::cout << "C_" << i << " "
//                << cacheC_LO[i].gslpp::complex::real() << " "
//                << cacheC_NLO[i].gslpp::complex::real() << "\n";     
//    }
//    std::cout << "--------\n" ;
   
    //LO DB=1 Wilson coefficients for 1/mb corrections    
    WilsonCoeffs = mySM.getFlavour().ComputeCoeffBMll_Buras(mu_1_overm, QCD::NOLEPTON);
    C_1LO = (*(WilsonCoeffs[LO]))(0);
    C_2LO = (*(WilsonCoeffs[LO]))(1);
    K_1 = 3. * C_1LO * C_1LO + 2. * C_1LO * C_2LO;
    K_2 = C_2LO * C_2LO;
}

gslpp::complex AmpDB2::C(int i){
    if (i>=1 and (i<=6 or i==8)) return cacheC[i - 1];
    throw std::runtime_error("Wilson cofficient out of order");
}

//hep-ph/0308029v2: eq.: 41-42: F0 = A
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

//hep-ph/0308029v2: F1 = B
//F0 has to be computed before if flag_resumz is enabled to resum z via hep-ph/0307344 eq.(23)
//see also 2106.05979 eq. (33): prefactor Mb2 in Gamma21 in MSbar scheme
void AmpDB2::computeF1() {
    double log_muM = 2. * log(mu_b/Mb_pole);
    cacheF1[indexF(cu, 1, 1, 1)] = 109./6. - 37. * z + 1.5 * z2 + 52./3. * z3 + 2. * oneminusz2 * (5. + z) * logx_2 - 4. * oneminusz2 * (5. + 7. * z) * log1minusz -
            2. * z * (10. + 14. * z - 15. * z2) * logz + 8. * (2. - 3. * z + z3) * log1minusz * logz + 16. * (2. - 3. * z + z3) * Dilogz
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(cu, 1, 1, 1) - 8. * z * logz * 1.5 * (-3. + 3. * z2));
    cacheF1[indexF(cu, 2, 1, 1)] = -4./3. * (10. - 33. * z + 54. * z2 - 31. * z3) - 8. * oneminusz2 * (4. + 14. * z - 3. * z2) * log1minusz +
            8. * z * (2. - 23. * z + 21. * z2 - 3. * z3) * logz -
            16. * oneminusz2 * (1. + 2. * z) * (2. * logx_2 - log1minusz * logz - 2. * Dilogz)
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(cu, 2, 1, 1) - 8. * z * logz * 3. * 6. * (z - 1.) * z);
    cacheF1[indexF(cu, 1, 1, 2)] = 2. * (
           (-502. + 912. * z - 387. * z2 - 23. * z3) / 36. - oneminusz2 * (17. + 4. * z) * logx_1 + 2./3. * oneminusz2 * (5. + z) * logx_2 -
           oneminusz2 / (12. * z) * (2. + 33. * z + 94. * z2) * log1minusz - z/12. * (80. + 69. * z - 126. * z2) * logz +
           8./3. * (2. - 3. * z + z3) * (log1minusz * logz + 2. * Dilogz)
           ) + flag_resumz * ((32./3. + 8. * log_muM) * F0(cu, 1, 1, 2) - 8. * z * logz * 0.5 * (-3. + 3.* z2));
    cacheF1[indexF(cu, 2, 1, 2)] = 2. * (
           (-130. + 93. * z + 144. * z2 - 107. * z3) / 9. - 2./3. * oneminusz2 / z * (1. + 15. * z + 47. * z2 - 12. * z3) * log1minusz +
           2./3. * z * (8. - 93. * z + 87. * z2 - 12. * z3) * logz -
           8./3. * oneminusz2 * (1. + 2. * z) * (3. * logx_1 + 4. * logx_2 - 2. * log1minusz * logz - 4. * Dilogz)
           ) + flag_resumz * ((32./3. + 8. * log_muM) * F0(cu, 2, 1, 2) - 8. * z * logz * 6. * (z - 1.) * z);
    cacheF1[indexF(cu, 1, 2, 2)] = -M_PI2/3. * (1. - 5. * z + 4. * z2) + (-136. - 159. * z + 738. * z2 - 443. * z3) / 18. - 2 * oneminusz2 * (5. + 4. * z) * logx_1 +
           2./3. * oneminusz2 * (4. - z) * logx_2 + oneminusz2 / (6. * z) * (7. + 32. * z2 + 3. * z3) * log1minusz -
           z/6. * (62. + 39. * z - 30. * z2 + 3. * z3) * logz + (5. - 3. * z - 18. * z2 + 16. * z3) / 3. * (log1minusz * logz + 2. * Dilogz)
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(cu, 1, 2, 2) - 8. * z * logz * -1.5 * oneminusz2);
    cacheF1[indexF(cu, 2, 2, 2)] = 8./3. * M_PI2 * (1. + z - 2. * z2) - 28./9. * (5. + 3. * z - 27. * z2 + 19. * z3) - 16. * oneminusz2 * (1. + 2. * z) * logx_1 +
           32./3. * oneminusz2 * (1. + 2. * z) * logx_2 - 4./3. * oneminusz2 / z * (1. - 12. * z - 16. * z2 - 3. * z3) * log1minusz +
           4./3. * z * (2. -  3. * z + 18. * z2 - 3. * z3) * logz + 8./3. * (1. - 3. * z - 6. * z2 + 8. * z3) * (log1minusz * logz + 2. * Dilogz)
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(cu, 2, 2, 2) - 8. * z * logz * -6. * (z - 1.) * z);

    cacheF1[indexF(cc, 1, 1, 1)] = sqrt1minus4z * (109. - 226. * z + 168. * z2) / 6. - (52. - 104. * z - 16. * z2 + 56. * z3) * logsigma +
           2. * (5. - 8. * z) * sqrt1minus4z * logx_2 - 12. * sqrt1minus4z * (3. - 2. * z) * log1minus4z + 4. * (13. - 10. * z) * sqrt1minus4z * logz +
           16. * (1. - 3. * z + 2. * z2) * (3. * log2sigma + 2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2)
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(cc, 1, 1, 1) - 8. * z * logz * 3. * (6. * z - 3.) / sqrt1minus4z);
    cacheF1[indexF(cc, 2, 1, 1)] = -8./3. * sqrt1minus4z * (5. - 23. * z - 42. * z2) - 16. * (4. - 2. * z - 7. * z2 + 14. * z3) * logsigma - 32. * sqrt1minus4z * (1. + 2. * z) * logx_2 -
           48. * sqrt1minus4z * (1. + 2. * z) * log1minus4z + 64. * sqrt1minus4z * (1. + 2. * z) * logz +
           16. * (1. - 4. * z2) * (3. * log2sigma + 2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2)
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(cc, 2, 1, 1) - 8. * z * logz * 3. * -12. * z / sqrt1minus4z);
    cacheF1[indexF(cc, 1, 1, 2)] = -sqrt1minus4z * (127. - 199. * z - 75. * z2) / 9. + (2. - 259. * z + 662. * z2 - 76. * z3 - 200. * z4) * logsigma / (12. * z) -
           (17. - 26. * z) * sqrt1minus4z * logx_1 + 2./3. * (5. - 8. * z) * sqrt1minus4z * logx_2 - 4. * sqrt1minus4z * (3. - 2. * z) * log1minus4z -
           sqrt1minus4z * (2. - 255. * z + 316. * z2) * logz / (12. * z) +
           16./3. * (1. - 3. * z + 2. * z2) * (3. * log2sigma + 2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2)
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(cc, 1, 1, 2) - 8. * z * logz * (6. * z - 3.) / sqrt1minus4z);
    cacheF1[indexF(cc, 2, 1, 2)] = -2. * sqrt1minus4z * (68. + 49. * z - 150. * z2) / 9. + 2./3. * (1. - 35. * z + 4. * z2 + 76. * z3 - 100. * z4) * logsigma/z +
           (16. - 64. * z2) * log2sigma - 8. * sqrt1minus4z * (1. + 2. * z) * logx_1 - 32./3. * sqrt1minus4z * (1. + 2. * z) * logx_2 -
           16. * sqrt1minus4z * (1. + 2. * z) * log1minus4z - 2./3. * sqrt1minus4z * (1. - 33. * z - 76. * z2) * logz/z +
           16./3. * (1. - 4. * z2) * (2. * logsigma * log1minus4z - 3. * logsigma * logz + 4. * Dilogsigma + 2. * Dilogsigma2)
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(cc, 2, 1, 2) - 8. * z * logz * -12. * z / sqrt1minus4z);
    cacheF1[indexF(cc, 1, 2, 2)] = -M_PI2/3. * (1. - 10. * z) - sqrt1minus4z * (115. + 632. * z + 96. * z2) / 18. - (7. + 13. * z - 194. * z2 + 304. * z3 - 64. * z4) * logsigma / (6. * z) -
           2. * sqrt1minus4z * (5. - 2. * z) * logx_1 + 4./3. * (2. - 5. * z) * sqrt1minus4z * logx_2 - 4. * (1. - 6. * z) * sqrt1minus4z * log1minus4z +
           (13. - 54. * z + 8. * z2) * logsigma * log1minus4z / 3. + sqrt1minus4z * (7. + 27. * z - 250. * z2) * logz / (6. * z) +
           (7. - 32. * z + 4. * z2) * (log2sigma - logsigma * logz) + 4./3. * (5. - 12. * z + 4. * z2) * Dilogsigma + 4./3. * (4. - 21. * z + 2. * z2) * Dilogsigma2
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(cc, 1, 2, 2) - 8. * z * logz * 0.5 * -6. * sqrt1minus4z);
    cacheF1[indexF(cc, 2, 2, 2)] = 8./3. * M_PI2 * (1. + 2. * z) - 8./9. * sqrt1minus4z * (19. + 53. * z + 24. * z2) + 4./3. * (1. + 7. * z + 10. * z2 - 68. * z3 + 32. * z4) * logsigma/z -
            8. * (1. + 2. * z) * (1. + 2. * z) * log2sigma - 16. * sqrt1minus4z * (1. + 2. * z) * logx_1 + 32./3. * sqrt1minus4z * (1. + 2. * z) * logx_2 +
            16. * sqrt1minus4z * (1. + 2. * z) * log1minus4z - 8./3. * (1. + 6. * z + 8. * z2) * logsigma * log1minus4z -
            4./3. * sqrt1minus4z * (1. + 9. * z + 26. * z2) * logz / z + 8. * (1. + 2. * z) * (1. + 2. * z) * logsigma * logz +
            32./3. * (1. - 4. * z2) * Dilogsigma - 32./3. * (1. + 3. * z + 2. * z2) * Dilogsigma2
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(cc, 2, 2, 2) - 8. * z * logz * 12. * z / sqrt1minus4z);

    cacheF1[indexF(uu, 1, 1, 1)] = 109./6. + 10. * logx_2
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(uu, 1, 1, 1));
    cacheF1[indexF(uu, 2, 1, 1)] = -40./3. - 32. * logx_2
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(uu, 2, 1, 1));
    cacheF1[indexF(uu, 1, 1, 2)] = -127./9. + 4./12. - 17. * logx_1 + 10./3. * logx_2
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(uu, 1, 1, 2));            
    cacheF1[indexF(uu, 2, 1, 2)] = -(136. + 12.)/9. - 8. * logx_1 - 32./3. * logx_2
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(uu, 2, 1, 2));            
    cacheF1[indexF(uu, 1, 2, 2)] = -M_PI2/3. - (115. + 42.)/18. - 10. * logx_1 + 8./3. * logx_2
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(uu, 1, 2, 2));            
    cacheF1[indexF(uu, 2, 2, 2)] = 8. * M_PI2/3. - 8./9. * (19. - 3.) - 16. * logx_1 + 32./3. * logx_2
            + flag_resumz * ((32./3. + 8. * log_muM) * F0(uu, 2, 2, 2));
           
    for (int k = 1; k < 3; k++) {
        for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1)) {
            cacheF1[indexF(qq, k, 2, 1)] = cacheF1[indexF(qq, k, 1, 2)];
        }
    }
    return;
}

//hep-ph/0308029v2: eq. 50-54
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
////using the FULLNLO DB=1 Wilson coefficients
////hep-ph/0308029v2: eq. 39:
//void AmpDB2::computeD() {
//    //qq = uu and cc
//    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2)) {
//        for (int k = 1; k <= 2; k++) {
//            gslpp::complex result = 0.;
//            for (int i = 1; i <= 2; i++) {
//                for (int j = 1; j <= 2; j++) {
//                    result += C(i) * C(j) * F(qq, k, i, j);
//                }
//            }
//            result += + as_4pi_mu1 * C(2) * C(2) * P(qq, k, 2, 2)
//                    + 2. * as_4pi_mu1 * C(2) * C(8) * P(qq, k, 2, 8);
//            for (int i = 1; i <= 2; i++) {
//                for (int r = 3; r <= 6; r++) {
//                    result += 2. * C(i) * C(r) * P(qq, k, i, r);
//                }                    
//            }
//            cacheD[indexD(qq, k)] = result;
//        }
//    }
//    //qq = cu
//    for (int k = 1; k <= 2; k++) {
//        gslpp::complex result = 0.;
//        for (int i = 1; i <= 2; i++) {
//            for (int j = 1; j <= 2; j++) {
//                result += C(i) * C(j) * F(cu, k, i, j);
//            }
//        }
//        result += + as_4pi_mu1 * C(2) * C(2) * P(cu, k, 2, 2)
//                + as_4pi_mu1 * C(2) * C(8) * (P(cc, k, 2, 8) + P(uu, k, 2, 8));
//        for (int i = 1; i <= 2; i++) {
//            for (int r = 3; r <= 6; r++) {
//                result += C(i) * C(r) * (P(cc, k, i, r) + P(uu, k, i, r));
//            }
//        }
//        cacheD[indexD(cu, k)] = result;
//    }         
//    return;
//}

//hep-ph/0308029v2: 39:
void AmpDB2::computeD() {
    //qq = uu and cc
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2)) {
        for (int k = 1; k <= 2; k++) {
            gslpp::complex result = 0.;
            for (int i = 1; i <= 2; i++) {
                for (int j = 1; j <= 2; j++) {
                    result += cacheC_LO[i-1] * cacheC_LO[j-1] * F(qq, k, i, j) + (cacheC_NLO[i-1] * cacheC_LO[j-1] + cacheC_LO[i-1] * cacheC_NLO[j-1]) * F0(qq, k, i, j);
                }
            }
            result += + as_4pi_mu1 * cacheC_LO[2-1] * cacheC_LO[2-1] * P(qq, k, 2, 2)
                    + 2. * as_4pi_mu1 * cacheC_LO[2-1] * cacheC_LO[8-1] * P(qq, k, 2, 8);
            for (int i = 1; i <= 2; i++) {
                for (int r = 3; r <= 6; r++) {
                    result += 2. * cacheC_LO[i-1] * cacheC_LO[r-1] * P(qq, k, i, r);
                }                    
            }
            cacheD[indexD(qq, k)] = result;
        }
    }
    //qq = cu
    for (int k = 1; k <= 2; k++) {
        gslpp::complex result = 0.;
        for (int i = 1; i <= 2; i++) {
            for (int j = 1; j <= 2; j++) {
                result += cacheC_LO[i-1] * cacheC_LO[j-1] * F(cu, k, i, j) + (cacheC_NLO[i-1] * cacheC_LO[j-1] + cacheC_LO[i-1] * cacheC_NLO[j-1]) * F0(cu, k, i, j);
            }
        }
        result += + as_4pi_mu1 * cacheC_LO[2-1] * cacheC_LO[2-1] * P(cu, k, 2, 2)
                + as_4pi_mu1 * cacheC_LO[2-1] * cacheC_LO[8-1] * (P(cc, k, 2, 8) + P(uu, k, 2, 8));
        for (int i = 1; i <= 2; i++) {
            for (int r = 3; r <= 6; r++) {
                result += cacheC_LO[i-1] * cacheC_LO[r-1] * (P(cc, k, i, r) + P(uu, k, i, r));
            }
        }
        cacheD[indexD(cu, k)] = result;
    }         
    return;
}

void AmpDB2::computeD_LO() {
    //qq = uu and cc
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2)) {
        for (int k = 1; k <= 2; k++) {
            gslpp::complex result = 0.;
            for (int i = 1; i <= 2; i++) {
                for (int j = 1; j <= 2; j++) {
                    result += cacheC_LO[i-1] * cacheC_LO[j-1] * F0(qq, k, i, j);
                }
            }
            for (int i = 1; i <= 2; i++) {
                for (int r = 3; r <= 6; r++) {
                    result += 2. * cacheC_LO[i-1] * cacheC_LO[r-1] * P(qq, k, i, r);
                }                    
            }
            cacheD[indexD(qq, k)] = result;
        }
    }
    //qq = cu
    for (int k = 1; k <= 2; k++) {
        gslpp::complex result = 0.;
        for (int i = 1; i <= 2; i++) {
            for (int j = 1; j <= 2; j++) {
                result += cacheC_LO[i-1] * cacheC_LO[j-1] * F0(cu, k, i, j);
            }
        }
        for (int i = 1; i <= 2; i++) {
            for (int r = 3; r <= 6; r++) {
                result += cacheC_LO[i-1] * cacheC_LO[r-1] * (P(cc, k, i, r) + P(uu, k, i, r));
            }
        }
        cacheD[indexD(cu, k)] = result;
    }         
    return;
}

double AmpDB2::F0(quarks qq, int k, int i, int j) {
    return cacheF0[indexF(qq, k, i, j)];
}

double AmpDB2::F1(quarks qq, int k, int i, int j) {
    return cacheF1[indexF(qq, k, i, j)];
}

double AmpDB2::F(quarks qq, int k, int i, int j) {
    return F0(qq, k, i, j) + as_4pi_mu1 * F1(qq, k, i, j);
}

double AmpDB2::P(quarks qq, int k, int i, int j) {
    return cacheP[indexP(qq, k, i, j)];
} 

//indizes for F to save them in an array
int AmpDB2::indexF(quarks qq, int k, int i, int j) {
    return qq * 8 + (k - 1) * 4 + (i - 1) * 2 + (j - 1);
}
//indizes for P to save them in an array: c := cc and u := uu
int AmpDB2::indexP(quarks qq, int k, int i, int j) {
    return qq * 28 + (k - 1) * 14 + (i - 1) * 7 + (j - 2);
}


void AmpDB2::compute_matrixelements(quark q){
    double Mq;
    double Mb_mu;
    double MBq2;
    double KBq;
    double FBq2;
    switch (q) {
        case d:
            me = mySM.getBBd().getBpars();
            Mq = Md;
            Mb_mu = mySM.getBBd().getMu();
            MBq2 = MB * MB;
            FBq2 = mySM.getMesons(QCD::B_D).getDecayconst() * mySM.getMesons(QCD::B_D).getDecayconst();
            break;
        case s:
            me = mySM.getBBs().getBpars();
            Mq = Ms;
            Mb_mu = mySM.getBBs().getMu();
            MBq2 = MB_s * MB_s;
            FBq2 = mySM.getMesons(QCD::B_S).getDecayconst() * mySM.getMesons(QCD::B_S).getDecayconst();            
            break;
        default:
            throw std::runtime_error("AmpDB2::compute_matrixelements(quark q): invalid quark index: ");
    }
    //pole mass for mbpow like in: hep-ph/0612167 eq. 28
    //double Mbpow = mySM.Mbar2Mp(mySM.getQuarks(QCD::BOTTOM).getMass());
    //double Mbpow2 = Mbpow * Mbpow;
    //KBq = MBq2 / ((Mbpow + Mq) * (Mbpow + Mq));
    
    
    //arXiv:1907.01025v2 equation (4)
    KBq = MBq2 / ((Mb_Mb + Mq) * (Mb_Mb + Mq));
    me(0) *=  8. / 3. * MBq2 * FBq2;
    me(1) *= -5. / 3. * KBq * MBq2 * FBq2;
    me(2) *=  1. / 3. * KBq * MBq2 * FBq2;
    me(3) *=       2. * (KBq + 1./6.) * MBq2 * FBq2;
    me(4) *=  2. / 3. * (KBq + 2./3.) * MBq2 * FBq2;
      
    //switch matrix elements to RI scheme
    if (flag_RI) me += -as_4pi_mu2 * meMStoRI * me;

    //old parameterization from hep-ph/0612167 eq. 28 (or hep-ph/0308029 eq.26)
    /*
    me_R(0) = me(1) + me(2) + 0.5 * me(0);
    me_R(1) = Mq/Mb * me(3);
    me_R(2) = -2./3. * FBq2 * MBq2 * (MBq2 / Mbpow2 - 1.);
    me_R(3) =  7./6. * FBq2 * MBq2 * (MBq2 / Mbpow2 - 1.);
    me_R(4) = 0.5 * (me(2) + 0.5 * me(0) + me(1) - 2. * Mq/Mbpow * me(4) + me_R(2));
    */
    
    //Gerlach thesis eq.7.5, 7.6
    switch (q) {
        case d:
            //me_R(0) = -0.35; //value in Gerlach thesis
            me_Rtilde(0) = mySM.getBBd_subleading().getBpars()(0);
            me_R(1) = mySM.getBBd_subleading().getBpars()(1);       
            me_R(2) = mySM.getBBd_subleading().getBpars()(2);                    
            me_R(3) = mySM.getBBd_subleading().getBpars()(3);            
            break;
        case s:
            //me_R(0) = -0.43; //value in Gerlach thesis
            me_Rtilde(0) = mySM.getBBs_subleading().getBpars()(0);
            me_R(1) = mySM.getBBs_subleading().getBpars()(1);       
            me_R(2) = mySM.getBBs_subleading().getBpars()(2);                    
            me_R(3) = mySM.getBBs_subleading().getBpars()(3);            
            break;
    } 
    me_R(0) *= MBq2 * FBq2;
    me_R(1) *= MBq2 * FBq2;
    me_R(2) *= MBq2 * FBq2;
    me_R(3) *= MBq2 * FBq2;
    me_R(4) = 0.5 * (me(2) + 0.5 * me(0) + me(1) - 2. * Mq/Mb_mu * me(4) + me_R(2));
    
    me_Rtilde(0) *= MBq2 * FBq2;        
    me_Rtilde(1) = -me_R(2);
    me_Rtilde(2) = me_R(3) + 0.5 * me_R(2);
    
    double n_l = 3.; //number of massless quark flavors
    double n_h = 1.; //number of quarks with mass of mb
    double L = 2. * log(mu_2/Mb_pole);
    double L2 = L * L;
    
    //Gerlach thesis eq. (3.84)
    double as1_me0 = 4. * L + 26./3.;
    double as1_me2 = 8. * L + 8.;
    
    //Gerlach thesis eq. (3.104, 3.105)
    double as2_me0 = (n_l + n_h) * (-4./3. * L2 - 52./9. * L - 8./9. * M_PI2 - 218./27.) + n_h * (8./3. * M_PI2 - 8.)
                        + 58./3. * L2 + 649./6. * L + 17./3. * M_PI2 + 11183./48. + 16./3. * M_PI2 * log2 - 8. * zeta3;
    double as2_me2 = (n_l + n_h) * (-8./3. * L2 - 104./9. * L - 16./9. * M_PI2 - 422./27.) + n_h * (16./3. * M_PI2 - 16.)
                        + 188./3. * L2 + 220. * L + 320./27. * M_PI2 + 326047./720. + 32./3. * M_PI2 * log2 - 16. * zeta3;
    //std::cout << "me_R" << me_R << "\n";
    me_R(0) = 0.5 * (1. + as1_me0 * as_4pi_mu2 + as2_me0 * as_4pi_mu2 * as_4pi_mu2) * me(0) + me(1) + (1. + as1_me2 * as_4pi_mu2 + as2_me2 * as_4pi_mu2 * as_4pi_mu2) * me(2);
    
    //std::cout << "me" << me << "\n";
    //std::cout << "me_R" << me_R(0) << " meR(0): " << 0.5 * (1. + 26./3. * as_4pi_mu2) * me(0) + me(1) + (1. + 8. * as_4pi_mu2) * me(2) << "\n";
    
    //fix leading matrix elements "me" to obtain only the uncertainties  from the subleading matrix elements "me_R"
//    if (q==s) {me(0) = 0.813; me(1) = 0.817; me(2) = 0.816;}
//    if (q==d) {me(0) = 0.806; me(1) = 0.769; me(2) = 0.747;}
//    me(0) *=  8. / 3. * MBq2 * FBq2;
//    me(1) *= -5. / 3. * KBq * MBq2 * FBq2;
//    me(2) *=  1. / 3. * KBq * MBq2 * FBq2;
    return;
}

//hep-ph/0308029v2 eq. 18
gslpp::vector<gslpp::complex> AmpDB2::c(quark q) {
    gslpp::vector< complex > c(3, 0.);
    switch (q) {
        case d:
            for (int i = 1; i <= 2; i++) {
                c.assign(i-1,
                        VtbVtd2 * D(uu, i) + 2. * VcbVcd * VtbVtd * (D(uu, i) - D(cu, i))
                        + VcbVcd2 * (D(cc, i) + D(uu, i) - 2. * D(cu, i))
                        );
            }
            //std::cout << "D " << D(uu, 1) << " " << D(uu,1)-D(cu,1) << " " << D(cc,1)+D(uu,1)-2.*D(cu,1) << "\n";            
            //std::cout << "D " << D(uu, 2) << " " << D(uu,2)-D(cu,2) << " " << D(cc,2)+D(uu,2)-2.*D(cu,2) << "\n";            
            break;
        case s:
            for (int i = 1; i <= 2; i++) {              
                c.assign(i-1,
                        VtbVts2 * D(uu, i) + 2. * VcbVcs * VtbVts * (D(uu, i) - D(cu, i))
                        + VcbVcs2 * (D(cc, i) + D(uu, i) - 2. * D(cu, i))
                        );
            }
            break;
        default:
            throw std::runtime_error("AmpDB2::c(quark q, double mu_2): invalid quark index: ");
    }
    //switch Wilson coefficients to RI scheme
    if (flag_RI) c += as_4pi_mu2 * coeffsMStoRI.transpose() * c;
    return c;
}

gslpp::complex AmpDB2::D(quarks qq, int k) {
    return cacheD[indexD(qq, k)];
}

//indizes for D to save them in arrays
int AmpDB2::indexD(quarks qq, int k) {
    if (k!=1 and k!=2){
        throw std::runtime_error("AmpDB2::indexD(quarks qq, int k): invalid k");       
    }
    return qq * 2 + (k - 1);
}

 /*******************************************************************************
 *  1/mb corrections of @f$\Gamma_{21}@f$                                       * 
 * ****************************************************************************/

gslpp::complex AmpDB2::delta_1overm_tradBasis(quark q) {
    //hep-ph/0308029: equation 18
    compute_deltas_1overm_tradBasis(q);
    switch (q) {
        case d:
            return VtbVtd2 * deltas_1overm_tradBasis(uu, d)
                    + 2. * VcbVcd * VtbVtd * (deltas_1overm_tradBasis(uu, d) - deltas_1overm_tradBasis(cu, d))
                    + VcbVcd2 * (deltas_1overm_tradBasis(cc, d) + deltas_1overm_tradBasis(uu, d) - 2. * deltas_1overm_tradBasis(cu, d));
        case s:
            return VtbVts2 * deltas_1overm_tradBasis(uu, s)
                    + 2. * VcbVcs * VtbVts * (deltas_1overm_tradBasis(uu, s) - deltas_1overm_tradBasis(cu, s))
                    + VcbVcs2 * (deltas_1overm_tradBasis(cc, s) + deltas_1overm_tradBasis(uu, s) - 2. * deltas_1overm_tradBasis(cu, s));
        default:
            throw std::runtime_error("AmpDB2::delta_1overm(quark q): invalid quark index: ");
    }
}


void AmpDB2::compute_deltas_1overm_tradBasis(quark q) {
    //hep-ph/0308029v2 eq.20
    cache_deltas_1overm_tradBasis[index_deltas(cc, q)] = sqrt1minus4z_1overm * ((1 + 2. * z_1overm) * (K_2 * (me_R(2) + 2. * me_R(4)) - 2. * K_1 * (me_R(1) + me_R(2)))
                    - 12. * z_1overm2 / (1. - 4. * z_1overm) * (K_1 * (me_R(2) + 2. * me_R(3)) + 2. * K_2 * me_R(3)));
    cache_deltas_1overm_tradBasis[index_deltas(cu, q)] = oneminusz_1overm2 * ((1. + 2. * z_1overm) * (K_2 * (me_R(2) + 2. * me_R(4)) - 2. * K_1 * (me_R(1) + me_R(2)))
                    - 6. * z_1overm2 / (1. - z_1overm) * (K_1 * (me_R(2) + 2. * me_R(3)) + 2. * K_2 * me_R(3)));
    cache_deltas_1overm_tradBasis[index_deltas(uu, q)] = K_2 * (me_R(2) + 2. * me_R(4)) - 2. * K_1 * (me_R(1) + me_R(2));
    return;
}


gslpp::complex AmpDB2::deltas_1overm_tradBasis(quarks qq, quark q) {
    return cache_deltas_1overm_tradBasis[index_deltas(qq, q)];
}

gslpp::complex AmpDB2::deltas_1overm(quarks qq, quark q) {
    return cache_deltas_1overm[index_deltas(qq, q)];
}

//indizes for deltas_1overm to save them in an array
int AmpDB2::index_deltas(quarks qq, quark q) {
    return qq * 2 + q;
}


gslpp::complex AmpDB2::delta_1overm(quark q) {
    //hep-ph/0612167 equation (10)
    switch (q) {
        case d:
            lambda_c = mySM.getCKM().getV_cd().conjugate() * mySM.getCKM().getV_cb();
            lambda_u = mySM.getCKM().getV_ud().conjugate() * mySM.getCKM().getV_ub();
            break;
        case s:
            lambda_c = mySM.getCKM().getV_cs().conjugate() * mySM.getCKM().getV_cb();
            lambda_u = mySM.getCKM().getV_us().conjugate() * mySM.getCKM().getV_ub();
            break;
        default:
            throw std::runtime_error("AmpDB2::delta_1overm(quark q): invalid quark index: ");
    }
    compute_deltas_1overm(q);
    return -(lambda_c * lambda_c * deltas_1overm(cc, q) + 2. * lambda_c * lambda_u * deltas_1overm(cu, q)
            + lambda_u * lambda_u * deltas_1overm(uu, q)).conjugate();
}


void AmpDB2::compute_deltas_1overm(quark q) {
    //hep-ph/0612167 eq.24
    compute_g();
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1)) {
        cache_deltas_1overm[index_deltas(qq, q)] = 0.;
        for (int i=0; i<=3; i++) {
            cache_deltas_1overm[index_deltas(qq, q)] += g(qq, i) * me_R(i);
        }
        for (int i=0; i<=2; i++){
            cache_deltas_1overm[index_deltas(qq, q)] += gtilde(qq, i) * me_Rtilde(i);
        }
    }
    return;
}


void AmpDB2::compute_g(){
    //hep-ph/0612167
    //equation (25)
    double cache_zc = z_1overm;
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2)) {
        cacheg[indexg(qq, 0)] = sqrt1minus4z_1overm * (1. + 2. * z_1overm) * K_1;
        cacheg[indexg(qq, 1)] = -2. * sqrt1minus4z_1overm * (1. + 2. * z_1overm) * K_1;
        cacheg[indexg(qq, 2)] = -2. * (1. - 2. * z_1overm - 2. * z_1overm2) / sqrt1minus4z_1overm * K_1;
        cacheg[indexg(qq, 3)] = -24. * z_1overm2 / sqrt1minus4z_1overm * K_1;
        cachegtilde[indexg(qq, 0)] = -2. * sqrt1minus4z_1overm * (1. + 2. * z_1overm) * K_2;
        cachegtilde[indexg(qq, 1)] = -2. * (1. - 2. * z_1overm - 2. * z_1overm2) / sqrt1minus4z_1overm * K_2;
        cachegtilde[indexg(qq, 2)] = -24. * z_1overm2 / sqrt1minus4z_1overm * K_2;
        z_1overm = 0;
    }
    //equation (26)
    z_1overm = cache_zc;
    cacheg[indexg(cu, 0)] = oneminusz_1overm2 * (1. + 2. * z_1overm) * K_1;
    cacheg[indexg(cu, 1)] = -2. * oneminusz_1overm2 * (1. + 2. * z_1overm) * K_1;
    cacheg[indexg(cu, 2)] = -2. * (1. - z_1overm) * (1. + z_1overm + z_1overm2) * K_1;;
    cacheg[indexg(cu, 3)] = -12. * (1. - z_1overm) * z_1overm2 * K_1;
    cachegtilde[indexg(cu, 0)] = -2. * oneminusz_1overm2 * (1. + 2. * z_1overm) * K_2;
    cachegtilde[indexg(cu, 1)] = -2. * (1. - z_1overm) * (1. + z_1overm + z_1overm2) * K_2;
    cachegtilde[indexg(cu, 2)] = -12. * (1. - z_1overm) * z_1overm2 * K_2;
    return;
}


gslpp::complex AmpDB2::g(quarks qq, int i){
    return cacheg[indexg(qq, i)];
}

gslpp::complex AmpDB2::gtilde(quarks qq, int i){
    return cachegtilde[indexg(qq, i)];
}

int AmpDB2::indexg(quarks qq, int i){
    return qq * 4 + i;
}

 /*******************************************************************************
 *  @f$\Gamma_{21}@f$ in NNLO from Marvin Gerlach (2205.07907 and thesis)       * 
 * ****************************************************************************/

gslpp::complex AmpDB2::Gamma21overM21_Bd(orders order, mass_schemes mass_scheme) {
    std::cout.precision(4);
    
    //save the order that shall be computed
    orderofp[1] = true;
    orderofp[2] = true;
    if (order == LO) {
        orderofp[1] = false;
        orderofp[2] = false;        
        order = FULLNNLO;
    }
    else if (order == FULLNLO) {
        orderofp[2] = false;
        order = FULLNNLO;
    }
    if (order != FULLNNLO) throw std::runtime_error("AmpDB2::Gamma21overM21_Bd(): order not implemented");
 
    computeCKMandMasses(NNLO, mass_scheme);
    
    //calculate M_21 / <O_1>
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffBd(
            mySM.getBBd().getMu(),
            mySM.getBBd().getScheme());
    gslpp::complex M21overme0 = ((*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0)) /( 8. * MB);

    //calculate DB=1 Wilson coefficients    
    computeWilsonCoeffsMisiak();
    lambda_c = mySM.getCKM().getV_cd().conjugate() * mySM.getCKM().getV_cb();
    lambda_u = mySM.getCKM().getV_ud().conjugate() * mySM.getCKM().getV_ub();
    //std::cout << "lambda_u/t: " << lambda_u/(VtbVtd.conjugate()) << "\n";

    //calculate DB=2 coefficients in pole scheme for usage of "c_H()"
    compute_pp_s();
    
    //switch to another scheme if needed
    if (mass_scheme == MSbar) poletoMSbar_pp_s();
    if (mass_scheme == PS) poletoPS_pp_s();    

    //calculate DB=2 matrix elements for usage of "me" and "delta_1overm_tradBasis(quark)"    
    compute_matrixelements(d);
    
   //Gerlach thesis eq. 6.1 divided by M_21
    gslpp::complex Gamma21overM21_Bd = Mb2_prefactor * (c_H()(0) + c_H()(1) * me(1)/me(0) + c_H()(2) * me(2)/me(0)).conjugate();
    computeWilsonCoeffsDB1bsg(); /*calculate DB=1 Wilson coefficients in the basis for "delta_1overm" */  
    Gamma21overM21_Bd += Mb_PS * Mb_PS * delta_1overm(d)/me(0);
    Gamma21overM21_Bd *= Gf2 / (24 * M_PI * MB) / M21overme0;
    return Gamma21overM21_Bd;  

    //parameterization from Gerlach gives same result    
    /*
    gslpp::complex lambda_uovert = 0.0122 - 0.4203 * gslpp::complex::i();
    gslpp::complex lambda_uovert2 = lambda_uovert * lambda_uovert;   
    Gamma21overM21_Bd =-(
            (H(cc) + 2. * lambda_uovert * (H(cc) - H(cu)) + lambda_uovert2 * (H(uu) - 2. * H(cu) + H(cc)))
            + (H_s(cc) + 2. * lambda_uovert * (H_s(cc) - H_s(cu)) + lambda_uovert2 * (H_s(uu) - 2. * H_s(cu) + H_s(cc))) * me(2)/me(0)).conjugate();
    computeWilsonCoeffsDB1bsg();
    compute_deltas_1overm(d); 
    Gamma21overM21_Bd += -(deltas_1overm(cc, d) + 2. * lambda_uovert * (deltas_1overm(cc, d) - deltas_1overm(cu, d))
            + lambda_uovert2 * (deltas_1overm(uu, d) - 2. * deltas_1overm(cu, d) + deltas_1overm(cc, d))).conjugate() /me(0);
    Gamma21overM21_Bd *= Gf2 * Mb2_prefactor / (24 * M_PI * MB) / (M21overme0 / (mySM.getCKM().computelamt_d() * mySM.getCKM().computelamt_d()));
    return Gamma21overM21_Bd;
     */
}


gslpp::complex AmpDB2::Gamma21overM21_Bs(orders order, mass_schemes mass_scheme) {
    std::cout.precision(4); 
    
    //save the order that shall be computed    
    orderofp[1] = true;
    orderofp[2] = true;
    if (order == LO) {
        orderofp[1] = false;
        orderofp[2] = false;        
        order = FULLNNLO;
    }
    else if (order == FULLNLO) {
        orderofp[2] = false;
        order = FULLNNLO;
    }
    if (order != FULLNNLO) throw std::runtime_error("AmpDB2::Gamma21overM21_Bs(): order not implemented");
    
    computeCKMandMasses(NNLO, mass_scheme);

    //calculate M_21 / <O_1>
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getFlavour().ComputeCoeffBs(
            mySM.getBBs().getMu(),
            mySM.getBBs().getScheme());
    gslpp::complex M21overme0 = ((*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0)) / (8. * MB_s);

    //calculate DB=1 Wilson coefficients        
    computeWilsonCoeffsMisiak();
    lambda_c = mySM.getCKM().getV_cs().conjugate() * mySM.getCKM().getV_cb();
    lambda_u = mySM.getCKM().getV_us().conjugate() * mySM.getCKM().getV_ub();   
//    std::cout <<  "lambda_us/t: " << lambda_u/(VtbVts.conjugate()) << "\n";
//    std::cout << "delta " << mySM.getCKM().getdelta() << "\n";

    //calculate DB=2 coefficients in pole scheme for usage of "c_H()"    
    compute_pp_s();
    
    //switch to another scheme if needed    
    if (mass_scheme == MSbar) poletoMSbar_pp_s();
    if (mass_scheme == PS) poletoPS_pp_s();

    //calculate DB=2 matrix elements for usage of "me" and "delta_1overm_tradBasis(quark)"    
    compute_matrixelements(s);
    
   //Gerlach thesis eq. 6.1 divided by M_21
    gslpp::complex Gamma21overM21_Bs = Mb2_prefactor * (c_H()(0) + c_H()(1) * me(1)/me(0) + c_H()(2) * me(2)/me(0)).conjugate();
    computeWilsonCoeffsDB1bsg(); /*calculate DB=1 Wilson coefficients in the basis for "delta_1overm" */ 
    Gamma21overM21_Bs += Mb_PS * Mb_PS * delta_1overm(s)/me(0);
    Gamma21overM21_Bs *= Gf2 / (24 * M_PI * MB_s) / M21overme0;
    return Gamma21overM21_Bs;
    
    //parameterization from Gerlach gives same result
    /*
    computeWilsonCoeffsMisiak();    
    gslpp::complex lambda_uovert = -0.00865 + 0.01832 * gslpp::complex::i();
    gslpp::complex lambda_uovert2 = lambda_uovert * lambda_uovert;
    Gamma21overM21_Bs =-((H(cc) + 2. * lambda_uovert * (H(cc) - H(cu)) + lambda_uovert2 * (H(uu) - 2. * H(cu) + H(cc)))
            + (H_s(cc) + 2. * lambda_uovert * (H_s(cc) - H_s(cu)) + lambda_uovert2 * (H_s(uu) - 2. * H_s(cu) + H_s(cc))) * me(2)/me(0)).conjugate();
    compute_deltas_1overm(s); 
    Gamma21overM21_Bs += -(deltas_1overm(cc, s) + 2. * lambda_uovert * (deltas_1overm(cc, s) - deltas_1overm(cu, s))
            + lambda_uovert2 * (deltas_1overm(uu, s) - 2. * deltas_1overm(cu, s) + deltas_1overm(cc, s))).conjugate()/me(0);    
    Gamma21overM21_Bs *= Gf2 * Mb2_prefactor / (24 * M_PI * MB_s) / (M21overme0 / (mySM.getCKM().computelamt_s()*mySM.getCKM().computelamt_s()));
    */
    return Gamma21overM21_Bs;
}


void AmpDB2::computeWilsonCoeffsMisiak(){
    //NLO DB=1 Wilson coefficients C_i, i=1-6,8
    gslpp::vector<gslpp::complex> ** WilsonCoeffsMisiak = mySM.getFlavour().ComputeCoeffsgamma(mu_1);
    for (int i = 0; i < 8; i++) {
       if (i==6) i=7;
        cacheC[i] = (*(WilsonCoeffsMisiak[LO]))(i) + (*(WilsonCoeffsMisiak[NLO]))(i) + (*(WilsonCoeffsMisiak[NNLO]))(i);
        cacheC_LO[i] = (*(WilsonCoeffsMisiak[LO]))(i);
        cacheC_NLO[i] = (*(WilsonCoeffsMisiak[NLO]))(i);
        cacheC_NNLO[i] = (*(WilsonCoeffsMisiak[NNLO]))(i);
    }

    //Wilson coefficients stated in Gerlach (2022)
    /*
    cacheC[0] = -0.6367 + 0.2986 + 0.0455;
    cacheC[1] =  1.0389 - 0.0322 + 0.0026;
    cacheC[2] = -0.0078 + 0.0023 - 0.0005;
    cacheC[3] = -0.0898 - 0.0013 + 0.0042;
    cacheC[4] =  0.0007 - 0.0004 + 0.00005;
    cacheC[5] =  0.0016 - 0.0005 + 0.000008;
    cacheC[7] = -0.1580 - 0.0104 + 0.0057; 
    */                  
    /*
    for(int i=0; i<=7; i++){
        if(i==6) i++;
        std::cout << "C_" << i << " "
                << cacheC_LO[i].gslpp::complex::real() << " "
                << cacheC_NLO[i].gslpp::complex::real() << " "
                << cacheC_NNLO[i].gslpp::complex::real() << "\n";
    }
    std::cout << "--------\n" ;
    */
}


//Gerlach thesis eq. (6.2)
gslpp::vector<gslpp::complex> AmpDB2::c_H(){
    gslpp::vector< gslpp::complex > result = gslpp::vector< gslpp::complex >(3, 0.);    
    gslpp::vector< gslpp::complex > result_LO = gslpp::vector< gslpp::complex >(3, 0.);
    gslpp::vector< gslpp::complex > result_NLO = gslpp::vector< gslpp::complex >(3, 0.);
    gslpp::vector< gslpp::complex > result_NNLO = gslpp::vector< gslpp::complex >(3, 0.);
    
    
    result_LO.assign(0, -lambda_c*lambda_c * H(cc, LO) - 2. * lambda_c*lambda_u * H(cu, LO) - lambda_u*lambda_u * H(uu, LO));
    result_LO.assign(2, -lambda_c*lambda_c * H_s(cc, LO) - 2. * lambda_c*lambda_u * H_s(cu, LO) - lambda_u*lambda_u * H_s(uu, LO));
    result_NLO.assign(0, -lambda_c*lambda_c * H(cc, NLO) - 2. * lambda_c*lambda_u * H(cu, NLO) - lambda_u*lambda_u * H(uu, NLO));
    result_NLO.assign(2, -lambda_c*lambda_c * H_s(cc, NLO) - 2. * lambda_c*lambda_u * H_s(cu, NLO) - lambda_u*lambda_u * H_s(uu, NLO));
    result_NNLO.assign(0, -lambda_c*lambda_c * H(cc, NNLO) - 2. * lambda_c*lambda_u * H(cu, NNLO) - lambda_u*lambda_u * H(uu, NNLO));
    result_NNLO.assign(2, -lambda_c*lambda_c * H_s(cc, NNLO) - 2. * lambda_c*lambda_u * H_s(cu, NNLO) - lambda_u*lambda_u * H_s(uu, NNLO));
 
    if (orderofp[2]) {
        result = transformation(result_LO, NNLO) + transformation(result_NLO, NLO) + transformation(result_NNLO, NNLO);
    } else if (orderofp[1]) {
        result = transformation(result_LO, NLO) + transformation(result_NLO, LO);
    } else {
        result = transformation(result_LO, LO);
    }
    //change to RI
    if (flag_RI) result += as_4pi_mu2 * coeffsMStoRI.transpose() * result;
    return result;
}

gslpp::vector<gslpp::complex> AmpDB2::transformation(gslpp::vector< gslpp::complex > result, orders order){ 
    bool orderofH[3] = {false, false, false};    
    if (order == LO) orderofH[0] = true;
    if (order == NLO) {orderofH[0] = true; orderofH[1] = true;}
    if (order == NNLO) {orderofH[0] = true; orderofH[1] = true; orderofH[2] = true;}
    //transformation to adapt to DB=2 matrix elements at scale Mb_pole
    //Gerlach thesis 3.60 - 3.62: mu2=mu1 , mu1=4.2
    double logmu1mu2 = log(mu_2) - log(4.2);
    double nf = 5.;
    double e21 = -4.;
    double e31 = -4.;
    double e421 = 8., e521 = 0., e411 = 0., e511 = 8.;
    gslpp::complex transformedresult0 = result(0) * (1. - 4. * orderofH[1] * as_4pi * logmu1mu2
             + orderofH[2] * as_4pi * as_4pi * (logmu1mu2 * (-2./3. * nf * e21 + 2./9. * nf * e31 + 11. * e21 - 11./3. * e31 - 20./9. * nf + 109./3.)
                           + (52. - 8./3. * nf) * logmu1mu2 * logmu1mu2));
    gslpp::complex transformedresult1 = result(1) * (1. + 28./3. * orderofH[1] * as_4pi * logmu1mu2
            + orderofH[2] * as_4pi * as_4pi * (logmu1mu2 * (2./3. * nf * e421 - 2./9. * nf * e521 - 8./3. * e411 - 92./9. * e421 + 8./9. * e511 + 4. * e521 - 232./27. * nf + 484./3.)
                           + (56./9. * nf - 500./9.) * logmu1mu2 * logmu1mu2))
        + result(2) * (-4./3. * orderofH[1] * as_4pi * logmu1mu2
            + orderofH[2] * as_4pi * as_4pi * (logmu1mu2 * (2./3. * nf * e411 - 2./9. * nf * e511 - 182./9. * e411 - 2./3. * e421 + 22./3. * e511 + 2./9. * e521 + 40./27. * nf - 116./3.)
                           + (140./9. - 8./9. * nf) * logmu1mu2 * logmu1mu2));
    gslpp::complex transformedresult2 = result(1) * (-16./3. * orderofH[1] * as_4pi * logmu1mu2
            + orderofH[2] * as_4pi * as_4pi * (logmu1mu2 * (7./9. * nf * e421 + 1./3. * nf * e521 - 28./9. * e411 + 3./2. * e421 - 4./3. * e511 - 25./18. * e521 - 92./27. * nf - 82.)
                           + (560./9. - 32./9. * nf) * logmu1mu2 * logmu1mu2))
        + result(2) * (1. - 32./3. * orderofH[1] * as_4pi * logmu1mu2
            + orderofH[2] * as_4pi * as_4pi * (logmu1mu2 * (7./9. * nf * e411 + 1./3. * nf * e511 - 61./6. * e411 - 7./9. * e421 - 115./18. * e511 - 1./3. * e521 + 260./27. * nf - 422./3.)
                           + (1600./9. - 64./9. * nf) * logmu1mu2 * logmu1mu2));
    result.assign(0, transformedresult0);
    result.assign(1, transformedresult1);
    result.assign(2, transformedresult2);
    return result;
}

////using the FULLNLO DB=1 Wilson coefficients
////Gerlach thesis eq. (6.4)
//gslpp::complex AmpDB2::H(quarks qq){
//    gslpp::complex result = 0.;
//    for (int i=1; i<=8; i++){
//        if (i==7) i++;
//        for (int j=i; j<=8; j++){
//            if(j==7) j++;
//            result += C(i) * C(j) * p(qq, i, j);
//        }
//    }
//    return result;
//}
//gslpp::complex AmpDB2::H_s(quarks qq){
//    gslpp::complex result = 0.;
//    for (int i=1; i<=8; i++){
//        if (i==7) i++;
//        for (int j=i; j<=8; j++){
//            if(j==7) j++;
//            result += C(i) * C(j) * p_s(qq, i, j);
//        }
//    }
//    return result;
//}

//Gerlach thesis eq. (6.4)
gslpp::complex AmpDB2::H(quarks qq, orders order){
    gslpp::complex result = 0.;
    bool orderofH[3] = {false, false, false};    
    if (order == LO) orderofH[0] = true;
    if (order == NLO) orderofH[1] = true;
    if (order == NNLO) orderofH[2] = true;       
    for (int i=1; i<=8; i++){
        if (i==7) i++;
        for (int j=i; j<=8; j++){
            if(j==7) j++;
            result += orderofH[0] * cacheC_LO[i-1] * cacheC_LO[j-1] * p(qq, i, j, 0)
                    + orderofH[1] * as_4pi_mu1 * (cacheC_LO[i-1] * cacheC_LO[j-1] * p(qq, i, j, 1)
                        + (cacheC_NLO[i-1] * cacheC_LO[j-1] + cacheC_LO[i-1] * cacheC_NLO[j-1]) * p(qq, i, j, 0))
                    + orderofH[2] * as_4pi_mu1 * as_4pi_mu1 * (cacheC_LO[i-1] * cacheC_LO[j-1] * p(qq, i, j, 2)
                        + (cacheC_NLO[i-1] * cacheC_LO[j-1] + cacheC_LO[i-1] * cacheC_NLO[j-1]) * p(qq, i, j, 1)
                        + (cacheC_NNLO[i-1] * cacheC_LO[j-1] + cacheC_NLO[i-1] * cacheC_NLO[j-1] + cacheC_LO[i-1] * cacheC_NNLO[j-1]) * p(qq, i, j, 0));
        }
    }
    return result;
}
gslpp::complex AmpDB2::H_s(quarks qq, orders order){
    gslpp::complex result = 0.;
    bool orderofH[3] = {false, false, false};    
    if (order == LO) orderofH[0] = true;
    if (order == NLO) orderofH[1] = true;
    if (order == NNLO) orderofH[2] = true;  
    for (int i=1; i<=8; i++){
        if (i==7) i++;
        for (int j=i; j<=8; j++){
            if(j==7) j++;
            result += orderofH[0] * cacheC_LO[i-1] * cacheC_LO[j-1] * p_s(qq, i, j, 0)
                    + orderofH[1] * as_4pi_mu1 * (cacheC_LO[i-1] * cacheC_LO[j-1] * p_s(qq, i, j, 1)
                        + (cacheC_NLO[i-1] * cacheC_LO[j-1] + cacheC_LO[i-1] * cacheC_NLO[j-1]) * p_s(qq, i, j, 0))
                    + orderofH[2] * as_4pi_mu1 * as_4pi_mu1 * (cacheC_LO[i-1] * cacheC_LO[j-1] * p_s(qq, i, j, 2)
                        + (cacheC_NLO[i-1] * cacheC_LO[j-1] + cacheC_LO[i-1] * cacheC_NLO[j-1]) * p_s(qq, i, j, 1)
                        + (cacheC_NNLO[i-1] * cacheC_LO[j-1] + cacheC_NLO[i-1] * cacheC_NLO[j-1] + cacheC_LO[i-1] * cacheC_NNLO[j-1]) * p_s(qq, i, j, 0));

        }
    }
    return result;
}

//Gerlach thesis eq. (6.5)
double AmpDB2::p(quarks qq, int i, int j){
    return p(qq, i, j, 0) * orderofp[0] + 
            as_4pi_mu1 * p(qq, i, j, 1) * orderofp[1] +
            as_4pi_mu1 * as_4pi_mu1 * p(qq, i, j, 2) * orderofp[2];
            //+ as_4pi_mu1 * as_4pi_mu1 * as_4pi_mu1 * p(qq, i, j, 3);
}
double AmpDB2::p_s(quarks qq, int i, int j){
    return p_s(qq, i, j, 0) * orderofp[0] + 
            as_4pi_mu1 * p_s(qq, i, j, 1) * orderofp[1] +
            as_4pi_mu1 * as_4pi_mu1 * p_s(qq, i, j, 2) * orderofp[2];
            //+ as_4pi_mu1 * as_4pi_mu1 * as_4pi_mu1 * p_s(qq, i, j, 3);
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
    
    double cache_logz = logz;
    if (flag_resumz){
        logz = 0.;        
    }
    
    //remember value of z after setting it to 0 for calculation of uu coefficients
    double cache_z = z;
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2)) {    
        //Gerlach thesis eq. (6.6)
        cache_p[index_p(qq, 1, 1, 0)]= 23./72. - 11./6. * z;
        cache_p[index_p(qq, 1, 2, 0)]= 1./6. - 2. * z;
        cache_p[index_p(qq, 2, 2, 0)]= 1. - 3. * z;
        cache_ps[index_p(qq, 1, 1, 0)]= -5./9.;
        cache_ps[index_p(qq, 1, 2, 0)]= -4./3.;
        cache_ps[index_p(qq, 2, 2, 0)]= 1.;

        //Gerlach thesis eq. (6.7)
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
        
        //Gerlach thesis eq. (6.9)
        z = 0;
    }
    z = cache_z;
    
    //Gerlach thesis eq. (6.8)
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
    
    //Gerlach thesis eq. (6.9)    
    for (int i=3; i<=6; i++){
        for (int j=i; j<=6; j++){
            cache_p[index_p(uu, i, j, 0)] = cache_p[index_p(cc, i, j, 0)];
            cache_ps[index_p(uu, i, j, 0)] = cache_ps[index_p(cc, i, j, 0)];            
        }
    }
    
    //Gerlach thesis eq. (6.10)
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
    
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2)) {    
        //Gerlach thesis eq. (6.11)
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

        //Gerlach thesis eq. (6.13) and (6.15)
        z = 0;
    }
    z = cache_z;
    //Gerlach thesis eq. (6.15)
    cache_p[index_p(uu, 1, 4, 1)] += 5./9. * z;
    cache_p[index_p(uu, 1, 6, 1)] += 50./9. * z;
    cache_p[index_p(uu, 2, 4, 1)] += - 10./3. * z;
    cache_p[index_p(uu, 2, 6, 1)] += - 100./3. * z;

    cache_ps[index_p(uu, 1, 4, 1)] += 8./9. * z;
    cache_ps[index_p(uu, 1, 6, 1)] += 80./9. * z;
    cache_ps[index_p(uu, 2, 4, 1)] += - 16./3. * z;
    cache_ps[index_p(uu, 2, 6, 1)] += - 160./3. * z; 
    
    //Gerlach thesis eq. (6.13) and (6.16)
    for (int i=1; i<=2; i++){
        for (int j=i; j<=6; j++){
            cache_p[index_p(cu, i, j, 1)] =
                    0.5 * (cache_p[index_p(cc, i, j, 1)] + cache_p[index_p(uu, i, j, 1)]);
            cache_ps[index_p(cu, i, j, 1)] =
                    0.5 * (cache_ps[index_p(cc, i, j, 1)] + cache_ps[index_p(uu, i, j, 1)]);
        }
    }
    
    //Gerlach thesis eq. (6.18)
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
    
    //Gerlach thesis eq. (6.17)
    for (int i=3; i<=6; i++){
        for (int j=i; j<=6; j++){
            cache_p[index_p(cu, i, j, 1)] = cache_p[index_p(cc, i, j, 1)];
            cache_p[index_p(uu, i, j, 1)] = cache_p[index_p(cc, i, j, 1)];
            cache_ps[index_p(cu, i, j, 1)] = cache_ps[index_p(cc, i, j, 1)];
            cache_ps[index_p(uu, i, j, 1)] = cache_ps[index_p(cc, i, j, 1)];
        }
    }
    
    //Gerlach thesis eq. (6.19)
    cache_p[index_p(cc, 1, 8, 1)] = 5./18.;
    cache_p[index_p(cc, 2, 8, 1)] = -5./3.;
    cache_ps[index_p(cc, 1, 8, 1)] = 4./9.;
    cache_ps[index_p(cc, 2, 8, 1)] = -8./3.;
    
    //Gerlach thesis eq. (6.21)
    cache_p[index_p(cc, 3, 8, 1)] = -32./3.;
    cache_p[index_p(cc, 4, 8, 1)] = -169./18.;
    cache_p[index_p(cc, 5, 8, 1)] = -512./3.;
    cache_p[index_p(cc, 6, 8, 1)] = -992./9.;
    cache_ps[index_p(cc, 3, 8, 1)] = 64./3.;
    cache_ps[index_p(cc, 4, 8, 1)] = -20./9.;
    cache_ps[index_p(cc, 5, 8, 1)] = 1024./3.;
    cache_ps[index_p(cc, 6, 8, 1)] = 256./9.;
    
    //Gerlach thesis eq. (6.20)
    for (int i=1; i<=6; i++){
        cache_p[index_p(cu, i, 8, 1)] = cache_p[index_p(cc, i, 8, 1)];
        cache_p[index_p(uu, i, 8, 1)] = cache_p[index_p(cc, i, 8, 1)];
        cache_ps[index_p(cu, i, 8, 1)] = cache_ps[index_p(cc, i, 8, 1)];
        cache_ps[index_p(uu, i, 8, 1)] = cache_ps[index_p(cc, i, 8, 1)];
    }
    
    double L_12 = L_1 * L_1;
    double L_22 = L_2 * L_2;

    double log2z = logz * logz;
    double sqrtz = sqrt(z);
    
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2)) {
        //Gerlach thesis eq. (6.22)
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
    
        //Gerlach thesis eq. (6.26)
        z = 0;
    }
    z = cache_z;
    
    //Gerlach thesis eq. (6.26)
    cache_p[index_p(uu, 1, 1, 2)] += n_v * (-5.55556 * z * logz + 0.0617284 * L_1 * z + 4.60031 * z + 4.75206 * sqrtz);
    cache_p[index_p(uu, 1, 2, 2)] += n_v * (2.66667 * z * logz - 0.740741 * L_1 * z - 85.8705 * z + 48.2515 * sqrtz);
    cache_p[index_p(uu, 2, 2, 2)] += n_v * (-32. * z * logz + 2.22222 * L_1 * z + 70.6121 * z - 105.276 * sqrtz);
    cache_ps[index_p(uu, 1, 1, 2)] += n_v * (0.0987654 * L_1 * z - 26.8617 * z + 27.7812 * sqrtz);
    cache_ps[index_p(uu, 1, 2, 2)] += n_v * (-1.18519 * L_1 * z - 72.3265 * z + 87.73 * sqrtz);
    cache_ps[index_p(uu, 2, 2, 2)] += n_v * (3.55556 * L_1 * z + 176.979 * z - 105.276 * sqrtz);
    
    //Gerlach thesis eq. (6.27)
    for (int i=1; i<=2; i++){
        for (int j=i; j<=2; j++){
            cache_p[index_p(cu, i, j, 2)] = 0.5 * (cache_p[index_p(cc, i, j, 2)] + cache_p[index_p(uu, i, j, 2)]);
            cache_ps[index_p(cu, i, j, 2)] = 0.5 * (cache_ps[index_p(cc, i, j, 2)] + cache_ps[index_p(uu, i, j, 2)]);
        }
    }
    
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 2)) {
        //Gerlach thesis eq. (6.28)
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
        
        //Gerlach thesis eq. (6.31) and (6.29)
        z = 0;
    }
    z = cache_z;
    
    //Gerlach thesis eq. (6.29)
    cache_p[index_p(uu, 1, 8, 2)] += -10./9. * z;
    cache_p[index_p(uu, 2, 8, 2)] += 20./3. * z;
    cache_ps[index_p(uu, 1, 8, 2)] += -16./9. * z;
    cache_ps[index_p(uu, 1, 8, 2)] += 32./3. * z;
    
    //Gerlach thesis eq. (6.31)
    cache_p[index_p(uu, 3, 8, 2)] += -196./3. * z;
    cache_p[index_p(uu, 4, 8, 2)] += (-404./3. + 20./3. * M_PI2) * z;
    cache_p[index_p(uu, 5, 8, 2)] += -760./3. * z;
    cache_p[index_p(uu, 6, 8, 2)] += (-4222./3. + 200./3. * M_PI2) * z;
    cache_ps[index_p(uu, 3, 8, 2)] += 608./3.* z;
    cache_ps[index_p(uu, 4, 8, 2)] += (-128./3. + 32./3. * M_PI2) * z;
    cache_ps[index_p(uu, 5, 8, 2)] += 11456./3.* z;
    cache_ps[index_p(uu, 6, 8, 2)] += (-160./3. + 320./3. * M_PI2) * z;
    
    //as in Gerlach thesis eq. (6.20)
    for (int i=1; i<=6; i++){
        cache_p[index_p(cu, i, 8, 2)] =
                0.5 * (cache_p[index_p(cc, i, 8, 2)] + cache_p[index_p(uu, i, 8, 2)]);
        cache_ps[index_p(cu, i, 8, 2)] =
                0.5 * (cache_ps[index_p(cc, i, 8, 2)] + cache_ps[index_p(uu, i, 8, 2)]);
    }
    
    //Gerlach thesis eq. (6.32)
    cache_p[index_p(cc, 8, 8, 2)] = -13./18.;
    cache_p[index_p(cu, 8, 8, 2)] = -13./18.;
    cache_p[index_p(uu, 8, 8, 2)] = -13./18.;
    cache_ps[index_p(cc, 8, 8, 2)] = -68./9.;
    cache_ps[index_p(cu, 8, 8, 2)] = -68./9.;
    cache_ps[index_p(uu, 8, 8, 2)] = -68./9.;

    logz = cache_logz;
    lastInput_compute_pp_s[0] = z;
    lastInput_compute_pp_s[1] = mu_1;
    lastInput_compute_pp_s[2] = mu_2;
    return;
}

void AmpDB2::poletoMSbar_pp_s(){
    //arxiv:2106.05979 eq. (33)
    double logmuM = 2. * log(mu_b/Mb_pole);
    double n_l = 4.;
    PoletoMS_as1 = 4./3. + logmuM;
    PoletoMS_as2 = -(- 3019./288. -  M_PI2/3. - M_PI2/9. * log2 + 1./6. * zeta3 + (71./144. + M_PI2/18.) * n_l + logmuM * (-445./72. + n_l * 13./36. + logmuM * (-19./24. + n_l/12.)) - 4./3. * M_PI2/8. * Mc_Mc/Mb);
    double as_mub_mu1 = mySM.Alstilde5(mu_b)/as_4pi_mu1;
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1)) {
        for (int i=1; i<=6; i++){
            //not all terms used for n=2
            for (int j=i; j<=8; j++){
                if(j==3) j=8;
                if(i>=3) j=8;
                cache_p[index_p(qq, i, j, 2)] += as_mub_mu1 * as_mub_mu1 * (8. * PoletoMS_as1 * p(qq, i, j, 1) + 16.* (2. * PoletoMS_as2 + PoletoMS_as1 * PoletoMS_as1) * p(qq, i, j, 0));
                cache_ps[index_p(qq, i, j, 2)] += as_mub_mu1 * as_mub_mu1 * (8. * PoletoMS_as1 * p_s(qq, i, j, 1) + 16.* (2. * PoletoMS_as2 + PoletoMS_as1 * PoletoMS_as1) * p_s(qq, i, j, 0));                
            }
            for (int j=i; j<=8; j++){
                if(j==7) j++;
                cache_p[index_p(qq, i, j, 1)] += as_mub_mu1 * 8. * PoletoMS_as1 * p(qq, i, j, 0);
                cache_ps[index_p(qq, i, j, 1)] += as_mub_mu1 * 8. * PoletoMS_as1 * p_s(qq, i, j, 0);
            }            
        }
    }
    return;
}

void AmpDB2::poletoPS_pp_s(){
    //analog to arxiv:2106.05979 eq. (33) for PS mass
    double as_Mb_mu1 = as_4pi/as_4pi_mu1;
    for (quarks qq = cc; qq <= uu; qq = quarks(qq + 1)) {
        for (int i=1; i<=6; i++){
            //not all terms used for n=2
            for (int j=i; j<=8; j++){
                if(j==3) j=8;
                if(i>=3) j=8;
                cache_p[index_p(qq, i, j, 2)] += as_Mb_mu1 * as_Mb_mu1 * (8. * PoletoPS_as1 * p(qq, i, j, 1) + 16.* (2. * PoletoPS_as2 + PoletoPS_as1 * PoletoPS_as1) * p(qq, i, j, 0));
                cache_ps[index_p(qq, i, j, 2)] += as_Mb_mu1 * as_Mb_mu1 * (8. * PoletoPS_as1 * p_s(qq, i, j, 1) + 16.* (2. * PoletoPS_as2 + PoletoPS_as1 * PoletoPS_as1) * p_s(qq, i, j, 0));
            }            
            for (int j=i; j<=8; j++){
                if(j==7) j++;
                cache_p[index_p(qq, i, j, 1)] += as_Mb_mu1 * 8. * PoletoPS_as1 * p(qq, i, j, 0);
                cache_ps[index_p(qq, i, j, 1)] += as_Mb_mu1 * 8. * PoletoPS_as1 * p_s(qq, i, j, 0);
            }
        }
    }
    return;
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
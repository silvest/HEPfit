/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <gsl/gsl_sf.h>
#include <gslpp.h>
#include <boost/bind.hpp>
#include "std_make_vector.h"
//#include <limits>
#include "BXqll.h"
#include "StandardModel.h"
#include "gslpp_function_adapter.h"

#define MPI2 (M_PI * M_PI)

BXqll::BXqll(const StandardModel& SM_i, QCD::quark quark_i, QCD::lepton lep_i)
: mySM(SM_i), myF_1(), myF_2(), myHeff("CPMLQB", SM_i, QCD2, QED2), WC(15, 9, 0.)
{    
    lep = lep_i;
    quark = quark_i;    
    w_H = gsl_integration_cquad_workspace_alloc(100);
    CF = mySM.getCF();
    phi1 = 50./3. - 8. * MPI2 / 3.; // functions for Phi_u at nh = 2 and nl = 3 
    phi2 = 2. * (- 2048. / 9. * gslpp_special_functions::zeta(3) + 16987. / 54. -340. / 81. * MPI2)
            + 3. * (256. / 9. * gslpp_special_functions::zeta(3) -1009. / 27. + 308. / 81. * MPI2)
            - 41848. / 81. * gslpp_special_functions::zeta(3) + 578. / 81. * MPI2 * MPI2 
            - 104480. / 729. * MPI2 + 1571095. / 1458. - 848. / 27. * MPI2 * log(2.);

}

BXqll::~BXqll() 
{
}

std::vector<std::string> BXqll::initializeBXqllParameters()
{
    BXqllParameters = make_vector<std::string>() << "BI_lambda1" << "BI_lambda2";
    
    return (BXqllParameters);
}

void BXqll::updateParameters() 
{
    GF = mySM.getGF();
    Mlep = mySM.getLeptons(lep).getMass();
//    mu_b = mySM.getMub();
    
    //RYOUTARO'S VALUES
    mu_b = 5.0;
    alstilde = 0.0170686027;
    ale = 0.00757373;
    aletilde = ale / 4. / M_PI;
    kappa = aletilde / alstilde;
    Mb_pole = 4.8;
    Mc_pole = 2.04545;
    Lbl = 2. * log(Mb_pole / Mlep);
    
    mu_c = mySM.getMuc();
    Mb = mySM.getQuarks(QCD::BOTTOM).getMass(); // add the PS b mass
    Mc = mySM.getQuarks(QCD::CHARM).getMass();
    Ms = mySM.getQuarks(QCD::STRANGE).getMass();
    MW = mySM.Mw();
    abslambdat_over_Vcb = mySM.getCKM().computelamt_s().abs()/mySM.getCKM().getV_cb().abs();
    Vts_over_Vcb = mySM.getCKM().getV_ts().abs()/mySM.getCKM().getV_cb().abs();
    alsmu = mySM.Als(mu_b, FULLNNNLO, true);
    alsmuc = mySM.Als(mu_c, FULLNNNLO, true);
//    ale = mySM.Ale(mu_b, FULLNLO);
//    alstilde = alsmu / 4. / M_PI;
//    aletilde = ale / 4. / M_PI;
//    kappa = ale / alsmu;
    Mtau = mySM.getLeptons(QCD::TAU).getMass(); // pole mass?
//    Mb_pole = mySM.Mbar2Mp(Mb, FULLNNLO);
    //Mc_pole = mySM.Mbar2Mp(Mc, FULLNNLO); //*** Mbar2Mp does not receive Mc ***/
//    Mc_pole = Mc*(1.+4.*alsmuc/3./M_PI+alsmuc*alsmuc/M_PI/M_PI*(-1.0414*(1.-4.*Ms/3.*Mc)+13.4434));
    muh = mu_b/Mb_pole; // log(muh) uses the pole mass as stated in hep-ph/9910220
    z = Mc_pole*Mc_pole/Mb_pole/Mb_pole; //****** Must be pole masses ****/
//    Lbl = 2.*log(Mb/Mlep);
    
    BR_BXcenu = 0.1051; // Branching ratio of B -> Xc e nu
    C_ratio = 0.574; // Ratio of branching ratios as defined by Gambino, Misiak, arXiv:hep-ph/0104034
//    pre = BR_BXcenu * abslambdat_over_Vcb * abslambdat_over_Vcb * 4. / C_ratio;
    pre = BR_BXcenu * 0.9621 * 4. / C_ratio;
    
    lambda_1 = -0.362; //mySM.getOptionalParameter("BI_lambda1");
    lambda_2 = 0.12; //mySM.getOptionalParameter("BI_lambda2");
    
    //ISIDORI VALUES
//    z = 0.29*0.29;
//    mu_b = 5.0;
//    Mb = 4.9;
//    Mtau = 1.77;
//    muh = mu_b/Mb;
//    ale = 0.0078125;
//    abslambdat_over_Vcb = 0.97;
//    Vts_over_Vcb = 0.97;
//    alsmu = 0.215;

//    allcoeff_smm = mySM.getFlavour().ComputeCoeffsmumu(mu_b, NDR);
    allcoeffDF1 = myHeff.ComputeCoeff(mu_b); // TO BE CHANGED TO ALLCOEFF ONLY

    double alpha_kappa;
    
    for(unsigned int qed_ord = QED0; qed_ord <= QED2; qed_ord++)
        for(unsigned int qcd_ord = QCD0; qcd_ord <= QCD2; qcd_ord++)
        {
            alpha_kappa = pow(alstilde, qcd_ord) * pow(kappa, qed_ord);
            for (unsigned int i = 0; i < 15; i++)
                WC.assign(i, qcd_ord + 3*qed_ord, alpha_kappa * (myHeff.LowScaleCoeff((qcd_orders) qcd_ord, (qed_orders) qed_ord))(i));
        }    
}

double BXqll::getR_LOWQ2(double sh)
{
    updateParameters();
    
    //To test HeffDF1 Wilson coefficients and Expanded multiplications
    Test_WC_DF1();
    return 0.;
    
//    computeMi(sh);
//    return H_A(sh);
}

gslpp::complex BXqll::F19(double sh) 
{
    gslpp::complex i = gslpp::complex::i();
    
    return (F_19re(muh, z, sh) + i*F_19im(muh, z, sh));
}

gslpp::complex BXqll::F29(double sh) 
{
    gslpp::complex i = gslpp::complex::i();
    
    return (F_29re(muh, z, sh) + i*F_29im(muh, z, sh));
}

gslpp::complex BXqll::F17(double sh) 
{
    gslpp::complex i = gslpp::complex::i();
    
    return (F_17re(muh, z, sh) + i*F_17im(muh, z, sh));
}

gslpp::complex BXqll::F27(double sh) 
{
    gslpp::complex i = gslpp::complex::i();
    
    return (F_27re(muh, z, sh) + i*F_27im(muh, z, sh));
}

gslpp::complex BXqll::F87(double sh) 
{
    gslpp::complex i = gslpp::complex::i();
    double ash = asin(sqrt(sh)/2.);
    double umsh = 1.-sh;

    return (4.*M_PI*M_PI/27.*(2.+sh)/umsh/umsh/umsh/umsh-4./9.*(11.-16.*sh+8.*sh*sh)/umsh/umsh-
            8./9.*sqrt(sh*(4.-sh))/umsh/umsh/umsh*(9.-5.*sh+2.*sh*sh)*ash-16./3.*(2+sh)/
            umsh/umsh/umsh/umsh*ash*ash-8./9.*sh/umsh*log(sh)-32./9.*log(muh)-i*8./9.*M_PI);
}

double BXqll::F89(double sh) 
{
    double ash = asin(sqrt(sh)/2.);
    double umsh = 1.-sh;

    return (-8.*M_PI*M_PI/27.*(4.-sh)/umsh/umsh/umsh/umsh+8./9.*(5.-2.*sh)/umsh/umsh+
            16./9.*sqrt((4.-sh)/sh)/umsh/umsh/umsh*(4.+3.*sh-sh*sh)*ash+32./3.*(4.-sh)/
            umsh/umsh/umsh/umsh*ash*ash+16./9./umsh*log(sh));
}

double BXqll::F_17re(double muh, double z, double sh, int maxpow)
{
    return myF_1.F_17re(muh, z, sh, maxpow);
};

double BXqll::F_17im(double muh, double z, double sh, int maxpow)
{
    return myF_1.F_17im(muh, z, sh, maxpow);
};

double BXqll::F_19re(double muh, double z, double sh, int maxpow)
{
    return myF_1.F_19re(muh, z, sh, maxpow);
};

double BXqll::F_19im(double muh, double z, double sh, int maxpow)
{
    return myF_1.F_19im(muh, z, sh, maxpow);
};

double BXqll::F_27re(double muh, double z, double sh, int maxpow)
{
    return myF_2.F_27re(muh, z, sh, maxpow);
};

double BXqll::F_27im(double muh, double z, double sh, int maxpow)
{
    return myF_2.F_27im(muh, z, sh, maxpow);
};

double BXqll::F_29re(double muh, double z, double sh, int maxpow)
{
    return myF_2.F_29re(muh, z, sh, maxpow);
};

double BXqll::F_29im(double muh, double z, double sh, int maxpow)
{
    return myF_2.F_29im(muh, z, sh, maxpow);
};

double BXqll::DeltaF_19re(double muh, double z, double sh, int maxpow)
{
    return myF_1.DeltaF_19re(muh, z, sh, maxpow);
};

double BXqll::DeltaF_19im(double muh, double z, double sh, int maxpow)
{
    return myF_1.DeltaF_19im(muh, z, sh, maxpow);
};

double BXqll::DeltaF_29re(double muh, double z, double sh, int maxpow)
{
    return myF_2.DeltaF_29re(muh, z, sh, maxpow);
};

double BXqll::DeltaF_29im(double muh, double z, double sh, int maxpow)
{
    return myF_2.DeltaF_29im(muh, z, sh, maxpow);
};


/*********************************************************
 * Implementation of the notation of @cite Huber:2015sra *
 *********************************************************/

double BXqll::integrateH(std::string obs, double q_min, double q_max)
{
    updateParameters();

    old_handler = gsl_set_error_handler_off();
    
    double sh_min = q_min/Mb_pole/Mb_pole, sh_max = q_max/Mb_pole/Mb_pole; // pole mass as explicitly stated in hep-ph/051206
    
    FH = convertToGslFunction(boost::bind(&BXqll::getH, &(*this), obs, _1));
    
    if (gsl_integration_cquad(&FH, sh_min, sh_max, 1.e-5, 1.e-4, w_H, &aveH, &errH, NULL) != 0)
        return std::numeric_limits<double>::quiet_NaN();
    return aveH;
        
    gsl_set_error_handler(old_handler);
}

double BXqll::getH(std::string obs, double sh)
{
    updateParameters();
    computeMi(sh);
    
    if (obs == "T")
        return H_T(sh);
    else if (obs == "L")
        return H_L(sh);
    else if (obs == "A")
        return H_A(sh);
    else if (obs == "TL")
        return (H_T(sh) + H_L(sh));
    else
        throw std::runtime_error("BXqll::getH: Angular observable not implemented");
}

double BXqll::H_T(double sh)
{
    computeHij_T(sh);
    
    double Phi_ll = CCH_multiplication(Hij_T);
    
    // 1/mc^2 corrections
    for(unsigned int j = 0; j < 15; j++)
        for(unsigned int i = 0; i <= j; i++)
            Phi_ll += (WC(i, LO).conjugate() * WC(j, LO) * cij_T(i, j, sh) ).real();
    
    // 1/mc^2 corrections associated with C9
    gslpp::complex C9 = WC(8, int_qed(LO_QED)) + WC(8, int_qed(NLO_QED11));
    
    Phi_ll += (C9 * (WC(0, LO).conjugate() * cij_T(0, 8, sh) + WC(1, LO).conjugate() * cij_T(1, 8, sh)) ).real();
    
    // log-enhanced electromagnetic corrections
    for (unsigned int i = 0; i < 7; i++)
    {
        for (unsigned int j = i; j < 7; j++)
            Phi_ll += (WC(i, LO).conjugate() * WC(j, LO) * eij_T(i, j, sh) ).real();
        
        Phi_ll += (WC(i, LO).conjugate() * C9 * eij_T(i, 8, sh) ).real();
    }
    
    Phi_ll += (C9.abs2() * eij_T(8, 8, sh) ).real();
    Phi_ll += (WC(9, int_qed(NLO_QED11)).abs2() * eij_T(9, 9, sh) ).real();

    return pre * Phi_ll / Phi_u(FULLNLO_QED);
}

double BXqll::H_L(double sh)
{
    computeHij_L(sh);
    
    double Phi_ll = CCH_multiplication(Hij_L);
    
    // 1/mc^2 corrections
    for(unsigned int j = 0; j < 15; j++)
        for(unsigned int i = 0; i <= j; i++)
            Phi_ll += (WC(i, LO).conjugate() * WC(j, LO) * cij_L(i, j, sh) ).real();
    
    // 1/mc^2 corrections associated with C9
    gslpp::complex C9 = WC(8, int_qed(LO_QED)) + WC(8, int_qed(NLO_QED11));
    
    Phi_ll += (C9 * (WC(0, LO).conjugate() * cij_L(0, 8, sh) + WC(1, LO).conjugate() * cij_L(1, 8, sh)) ).real();
    
    // log-enhanced electromagnetic corrections
    for (unsigned int i = 0; i < 7; i++)
    {
        for (unsigned int j = i; j < 7; j++)
            Phi_ll += (WC(i, LO).conjugate() * WC(j, LO) * eij_L(i, j, sh) ).real();
        
        Phi_ll += (WC(i, LO).conjugate() * C9 * eij_L(i, 8, sh) ).real();
    }
    
    Phi_ll += (C9.abs2() * eij_L(8, 8, sh) ).real();
    Phi_ll += (WC(9, int_qed(NLO_QED11)).abs2() * eij_L(9, 9, sh) ).real();
    
    return pre * Phi_ll / Phi_u(FULLNLO_QED);
}

double BXqll::H_A(double sh)
{
    computeHij_A(sh);
    
    double Phi_ll = CCH_multiplication(Hij_A);
//    double Phi_ll = FULLCCH_multiplication(Hij_A);

    // 1/mc^2 corrections associated with C10
    Phi_ll += (WC(0, LO).conjugate() * WC(9, int_qed(NLO_QED11)) * cij_A(0, 9, sh) +
               WC(1, LO).conjugate() * WC(9, int_qed(NLO_QED11)) * cij_A(1, 9, sh) ).real();
    
    // log-enhanced electromagnetic corrections
    Phi_ll += (WC(0, LO).conjugate() * WC(9, int_qed(NLO_QED11)) * eij_A(0, 9, sh) ).real();
    Phi_ll += (WC(1, LO).conjugate() * WC(9, int_qed(NLO_QED11)) * eij_A(1, 9, sh) ).real();
    Phi_ll += (WC(6, LO).conjugate() * WC(9, int_qed(NLO_QED11)) * eij_A(6, 9, sh) ).real();
    Phi_ll += ((WC(8, int_qed(LO_QED)) + WC(8, int_qed(NLO_QED11))) * WC(9, int_qed(NLO_QED11)) * eij_A(8, 9, sh) ).real();
    
    return pre * Phi_ll / Phi_u(FULLNLO_QED);
}

void BXqll::computeHij_T(double sh)
{
    // Clears the vector for a new value of sh
    Hij_T.clear();
    
    gslpp::matrix<gslpp::complex> Hij(15, 15, 0.);

    // LO
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
        {
            if (i == j)
                Hij.assign(i, j,
                (M_9[LO])(i).abs2() * S99_T(sh, LO) +
                (M_10[LO])(i).abs2() * S1010_T(sh, LO) );

            else
                Hij.assign(i, j,
                2. * (M_9[LO])(i) * (M_9[LO])(j).conjugate() * S99_T(sh, LO) );
        }
    }
    
    Hij_T.push_back(Hij);
    
    
    // NLO
    Hij.reset(); // To clear Hij
    
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
        {
            if (i == j)
                Hij.assign(i, j,
                (M_9[LO])(i).abs2() * S99_T(sh, NLO) +
                (M_10[LO])(i).abs2() * S1010_T(sh, NLO) );

            else
                Hij.assign(i, j,
                2. * (M_9[LO])(i) * (M_9[LO])(j).conjugate() * S99_T(sh, NLO) );
        }
    }
    
    Hij_T.push_back(Hij);
    
    
    // NNLO -> LO_QED
    Hij.reset(); // To clear Hij
    Hij_T.push_back(Hij);
    Hij_T.push_back(Hij);
    
    
    // NLO_QED11
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
        {
            if (i == j)
                Hij.assign(i, j,
                ((M_9[LO])(i) * (M_9[int_qed(NLO_QED11)])(i).conjugate()).real() * 2. * S99_T(sh, LO) +
                ((M_7[LO])(i) * (M_9[int_qed(NLO_QED11)])(i).conjugate()).real() * S79_T(sh, LO) );

            else
                Hij.assign(i, j,
                2. * ((M_9[LO])(i) * (M_9[int_qed(NLO_QED11)])(j).conjugate() +
                      (M_9[int_qed(NLO_QED11)])(i) * (M_9[LO])(j).conjugate()) * S99_T(sh, LO) +
                ((M_7[int_qed(NLO_QED11)])(i) * (M_9[LO])(j).conjugate() +
                 (M_9[LO])(i) * (M_7[int_qed(NLO_QED11)])(j).conjugate()) * S79_T(sh, LO) );
        }
    }
    
    Hij_T.push_back(Hij);
    
    
    // NLO_QED21
    Hij.reset(); // To clear Hij
    
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
        {
            if (i == j)
                Hij.assign(i, j,
                ((M_9[LO])(i) * (M_9[int_qed(NLO_QED21)])(i).conjugate()).real() * 2. * S99_T(sh, LO) +
                ((M_9[LO])(i) * (M_9[int_qed(NLO_QED11)])(i).conjugate()).real() * 2. * S99_T(sh, NLO) +
                ((M_7[int_qed(NLO_QED21)])(i) * (M_9[LO])(j).conjugate() +
                 (M_9[LO])(i) * (M_7[int_qed(NLO_QED21)])(j).conjugate()).real() * S79_T(sh, LO) +
                ((M_7[int_qed(NLO_QED11)])(i) * (M_9[LO])(j).conjugate() +
                 (M_9[LO])(i) * (M_7[int_qed(NLO_QED11)])(j).conjugate()).real() * S79_T(sh, NLO) );

            else
                Hij.assign(i, j,
                2. * ((M_9[LO])(i) * (M_9[int_qed(NLO_QED21)])(j).conjugate() +
                      (M_9[int_qed(NLO_QED21)])(i) * (M_9[LO])(j).conjugate()) * S99_T(sh, LO) +
                2. * ((M_9[LO])(i) * (M_9[int_qed(NLO_QED11)])(j).conjugate() +
                      (M_9[int_qed(NLO_QED11)])(i) * (M_9[LO])(j).conjugate()) * S99_T(sh, NLO) +
                ((M_7[int_qed(NLO_QED21)])(i) * (M_9[LO])(j).conjugate() +
                 (M_9[LO])(i) * (M_7[int_qed(NLO_QED21)])(j).conjugate()) * S79_T(sh, LO) +
                ((M_7[int_qed(NLO_QED11)])(i) * (M_9[LO])(j).conjugate() +
                 (M_9[LO])(i) * (M_7[int_qed(NLO_QED11)])(j).conjugate()) * S79_T(sh, NLO) );
        }
    }
    
    Hij_T.push_back(Hij);
    
    
    // NLO_QED02 -> NLO_QED12
    Hij.reset(); // To clear Hij
    Hij_T.push_back(Hij);
    Hij_T.push_back(Hij);
    
    
    // NLO_QED22
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
        {
            if (i == j)
                Hij.assign(i, j,
                (M_7[int_qed(NLO_QED11)])(i).abs2() * S77_T(sh, LO) +
                (M_9[int_qed(NLO_QED11)])(i).abs2() * S99_T(sh, LO) +
                ((M_7[int_qed(NLO_QED11)])(i) * (M_9[int_qed(NLO_QED11)])(i).conjugate()).real() * S79_T(sh, LO) );

            else
                Hij.assign(i, j,
                2. * (M_7[int_qed(NLO_QED11)])(i) * (M_7[int_qed(NLO_QED11)])(j).conjugate() * S77_T(sh, LO) +
                2. * (M_9[int_qed(NLO_QED11)])(i) * (M_9[int_qed(NLO_QED11)])(j).conjugate() * S99_T(sh, LO) +
                ((M_7[int_qed(NLO_QED11)])(i) * (M_9[int_qed(NLO_QED11)])(j).conjugate() +
                 (M_9[int_qed(NLO_QED11)])(i) * (M_7[int_qed(NLO_QED11)])(j).conjugate()) * S79_T(sh, LO) );
        }
    }
    
    Hij_T.push_back(Hij);
}

void BXqll::computeHij_L(double sh)
{
    // Clears the vector for a new value of sh
    Hij_L.clear();
    
    gslpp::matrix<gslpp::complex> Hij(15, 15, 0.);

    // LO
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
        {
            if (i == j)
                Hij.assign(i, j,
                (M_9[LO])(i).abs2() * S99_L(sh, LO) +
                (M_10[LO])(i).abs2() * S1010_L(sh, LO) );

            else
                Hij.assign(i, j,
                2. * (M_9[LO])(i) * (M_9[LO])(j).conjugate() * S99_L(sh, LO) );
        }
    }
    
    Hij_L.push_back(Hij);
    
    
    // NLO
    Hij.reset(); // To clear Hij
    
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
        {
            if (i == j)
                Hij.assign(i, j,
                (M_9[LO])(i).abs2() * S99_L(sh, NLO) +
                (M_10[LO])(i).abs2() * S1010_L(sh, NLO) );

            else
                Hij.assign(i, j,
                2. * (M_9[LO])(i) * (M_9[LO])(j).conjugate() * S99_L(sh, NLO) );
        }
    }
    
    Hij_L.push_back(Hij);
    
    
    // NNLO -> LO_QED
    Hij.reset(); // To clear Hij
    Hij_L.push_back(Hij);
    Hij_L.push_back(Hij);
    
    
    // NLO_QED11
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
        {
            if (i == j)
                Hij.assign(i, j,
                ((M_9[LO])(i) * (M_9[int_qed(NLO_QED11)])(i).conjugate()).real() * 2. * S99_L(sh, LO) +
                ((M_7[LO])(i) * (M_9[int_qed(NLO_QED11)])(i).conjugate()).real() * S79_L(sh, LO) );

            else
                Hij.assign(i, j,
                2. * ((M_9[LO])(i) * (M_9[int_qed(NLO_QED11)])(j).conjugate() +
                      (M_9[int_qed(NLO_QED11)])(i) * (M_9[LO])(j).conjugate()) * S99_L(sh, LO) +
                ((M_7[int_qed(NLO_QED11)])(i) * (M_9[LO])(j).conjugate() +
                 (M_9[LO])(i) * (M_7[int_qed(NLO_QED11)])(j).conjugate()) * S79_L(sh, LO) );
        }
    }
    
    Hij_L.push_back(Hij);
    
    
    // NLO_QED21
    Hij.reset(); // To clear Hij
    
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
        {
            if (i == j)
                Hij.assign(i, j,
                ((M_9[LO])(i) * (M_9[int_qed(NLO_QED21)])(i).conjugate()).real() * 2. * S99_L(sh, LO) +
                ((M_9[LO])(i) * (M_9[int_qed(NLO_QED11)])(i).conjugate()).real() * 2. * S99_L(sh, NLO) +
                ((M_7[int_qed(NLO_QED21)])(i) * (M_9[LO])(j).conjugate() +
                 (M_9[LO])(i) * (M_7[int_qed(NLO_QED21)])(j).conjugate()).real() * S79_L(sh, LO) +
                ((M_7[int_qed(NLO_QED11)])(i) * (M_9[LO])(j).conjugate() +
                 (M_9[LO])(i) * (M_7[int_qed(NLO_QED11)])(j).conjugate()).real() * S79_L(sh, NLO) );

            else
                Hij.assign(i, j,
                2. * ((M_9[LO])(i) * (M_9[int_qed(NLO_QED21)])(j).conjugate() +
                      (M_9[int_qed(NLO_QED21)])(i) * (M_9[LO])(j).conjugate()) * S99_L(sh, LO) +
                2. * ((M_9[LO])(i) * (M_9[int_qed(NLO_QED11)])(j).conjugate() +
                      (M_9[int_qed(NLO_QED11)])(i) * (M_9[LO])(j).conjugate()) * S99_L(sh, NLO) +
                ((M_7[int_qed(NLO_QED21)])(i) * (M_9[LO])(j).conjugate() +
                 (M_9[LO])(i) * (M_7[int_qed(NLO_QED21)])(j).conjugate()) * S79_L(sh, LO) +
                ((M_7[int_qed(NLO_QED11)])(i) * (M_9[LO])(j).conjugate() +
                 (M_9[LO])(i) * (M_7[int_qed(NLO_QED11)])(j).conjugate()) * S79_L(sh, NLO) );
        }
    }
    
    Hij_L.push_back(Hij);
    
    
    // NLO_QED02 -> NLO_QED12
    Hij.reset(); // To clear Hij
    Hij_L.push_back(Hij);
    Hij_L.push_back(Hij);
    
    
    // NLO_QED22
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
        {
            if (i == j)
                Hij.assign(i, j,
                (M_7[int_qed(NLO_QED11)])(i).abs2() * S77_L(sh, LO) +
                (M_9[int_qed(NLO_QED11)])(i).abs2() * S99_L(sh, LO) +
                ((M_7[int_qed(NLO_QED11)])(i) * (M_9[int_qed(NLO_QED11)])(i).conjugate()).real() * S79_L(sh, LO) );

            else
                Hij.assign(i, j,
                2. * (M_7[int_qed(NLO_QED11)])(i) * (M_7[int_qed(NLO_QED11)])(j).conjugate() * S77_L(sh, LO) +
                2. * (M_9[int_qed(NLO_QED11)])(i) * (M_9[int_qed(NLO_QED11)])(j).conjugate() * S99_L(sh, LO) +
                ((M_7[int_qed(NLO_QED11)])(i) * (M_9[int_qed(NLO_QED11)])(j).conjugate() +
                 (M_9[int_qed(NLO_QED11)])(i) * (M_7[int_qed(NLO_QED11)])(j).conjugate()) * S79_L(sh, LO) );
        }
    }
    
    Hij_L.push_back(Hij);
}

void BXqll::computeHij_A(double sh)
{
    // Clears the vector for a new value of sh
    Hij_A.clear();
    
    gslpp::matrix<gslpp::complex> Hij(15, 15, 0.);

    // LO
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
            if (i != j)
                Hij.assign(i, j,
                ((M_9[LO])(i) * (M_10[LO])(j).conjugate() +
                 (M_10[LO])(i) * (M_9[LO])(j).conjugate()) * S910_A(sh, LO) );
    }
    
    Hij_A.push_back(Hij);
    
    
    // NLO
    Hij.reset(); // To clear Hij
    
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
            if (i != j)
                Hij.assign(i, j,
                ((M_9[LO])(i) * (M_10[LO])(j).conjugate() +
                 (M_10[LO])(i) * (M_9[LO])(j).conjugate()) * S910_A(sh, NLO) );
    }
    
    Hij_A.push_back(Hij);
    
    
    // NNLO -> LO_QED
    Hij.reset(); // To clear Hij
    Hij_A.push_back(Hij);
    Hij_A.push_back(Hij);
    
    
    // NLO_QED11
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
            if (i != j)
                Hij.assign(i, j,
                ((M_7[int_qed(NLO_QED11)])(i) * (M_10[LO])(j).conjugate() +
                 (M_10[LO])(i) * (M_7[int_qed(NLO_QED11)])(j).conjugate()) * S710_A(sh, LO) +
                ((M_9[int_qed(NLO_QED11)])(i) * (M_10[LO])(j).conjugate() +
                 (M_10[LO])(i) * (M_9[int_qed(NLO_QED11)])(j).conjugate()) * S910_A(sh, LO) );
    }
    
    Hij_A.push_back(Hij);
    
    
    // NLO_QED21
    Hij.reset(); // To clear Hij
    
    for (unsigned int j = 0; j < 15; j++)
    {
        for (unsigned int i = 0; i <= j; i++)
            if (i != j)
                Hij.assign(i, j,
                ((M_7[int_qed(NLO_QED21)])(i) * (M_10[LO])(j).conjugate() +
                 (M_10[LO])(i) * (M_7[int_qed(NLO_QED21)])(j).conjugate()) * S710_A(sh, LO) +
                ((M_9[int_qed(NLO_QED21)])(i) * (M_10[LO])(j).conjugate() +
                 (M_10[LO])(i) * (M_9[int_qed(NLO_QED21)])(j).conjugate()) * S910_A(sh, LO) +
                ((M_7[int_qed(NLO_QED11)])(i) * (M_10[LO])(j).conjugate() +
                 (M_10[LO])(i) * (M_7[int_qed(NLO_QED11)])(j).conjugate()) * S710_A(sh, NLO) +
                ((M_9[int_qed(NLO_QED11)])(i) * (M_10[LO])(j).conjugate() +
                 (M_10[LO])(i) * (M_9[int_qed(NLO_QED11)])(j).conjugate()) * S910_A(sh, NLO) );
    }
    
    Hij_A.push_back(Hij);
    
    
//    // NLO_QED02 -> NLO_QED22
    Hij.reset(); // To clear Hij
    Hij_A.push_back(Hij);
    Hij_A.push_back(Hij);    
    Hij_A.push_back(Hij);
}

void BXqll::computeMi(double sh)
{
    // Clears the vectors for a new value of sh
    M_7.clear();
    M_9.clear();
    M_10.clear();
    
    gslpp::vector<gslpp::complex> M7i(15, 0.);
    gslpp::vector<gslpp::complex> M9i(15, 0.);
    gslpp::vector<gslpp::complex> M10i(15, 0.);
    
    // M_7: LO -> LO_QED
    for (unsigned int i = LO; i <= int_qed(LO_QED); i++)
        M_7.push_back(M7i);
    
    // M_7: NLO_QED11
    M7i.assign(6, aletilde);
    M_7.push_back(M7i);
    
    // M_7: NLO_QED21
    M7i.reset(); // To clear M7i
    M7i.assign(0, -alstilde * aletilde * F17(sh));
    M7i.assign(1, -alstilde * aletilde * F27(sh));
    M7i.assign(7, -alstilde * aletilde * F87(sh));
    M_7.push_back(M7i);
    
    // M_9: LO
    M9i.assign(8, 1.);
    M_9.push_back(M9i);
    
    // M_9: NLO -> LO_QED
    M9i.reset(); // To clear M9i
    for (unsigned int i = NLO; i <= int_qed(LO_QED); i++)
        M_9.push_back(M9i);
    
    // M_9: NLO_QED11
    M9i.assign(0,  aletilde * f_Huber(sh,   -32./27., 4./3.,      0.,        0.,   -16./27.));
    M9i.assign(1,  aletilde * f_Huber(sh,     -8./9.,    1.,      0.,        0.,     -4./9.));
    M9i.assign(2,  aletilde * f_Huber(sh,    -16./9.,    6.,  -7./2.,     2./9.,     2./27.));
    M9i.assign(3,  aletilde * f_Huber(sh,    32./27.,    0.,  -2./3.,    8./27.,     8./81.));
    M9i.assign(4,  aletilde * f_Huber(sh,   -112./9.,   60.,    -38.,    32./9.,  -136./27.));
    M9i.assign(5,  aletilde * f_Huber(sh,   512./27.,    0., -32./3.,  128./27.,   320./81.));
    M9i.assign(10, aletilde * f_Huber(sh,  -272./27.,    4.,   7./6.,  -74./27.,   358./81.));
    M9i.assign(11, aletilde * f_Huber(sh,   -32./81.,    0.,   2./9.,   -8./81.,   -8./243.));
    M9i.assign(12, aletilde * f_Huber(sh, -2768./27.,   40.,  38./3., -752./27.,  1144./81.));
    M9i.assign(13, aletilde * f_Huber(sh,  -512./81.,    0.,  32./9., -128./81., -320./243.));
    M9i.assign(14, aletilde * f_Huber(sh,     16./9.,    0.,     -2.,        0.,    26./27.));
    M9i.assign(8,  aletilde * f9pen_Huber(sh));
    M_9.push_back(M9i);

    // M_9: NLO_QED21
    M9i.reset(); // To clear M9i
    M9i.assign(0, -alstilde * aletilde * F19(sh));
    M9i.assign(1, -alstilde * aletilde * F29(sh));
    M9i.assign(7, -alstilde * aletilde * F89(sh));
    M_9.push_back(M9i);

    // M_10: LO
    M10i.assign(9, 1.);
    M_10.push_back(M10i);
    
    // M_10: NLO -> NLO_QED21
    M10i.reset(); // To clear M7i
    for (unsigned int i = NLO; i <= int_qed(NLO_QED21); i++)
        M_10.push_back(M10i);
}

double BXqll::S77_T(double sh, orders order)
{
    double umsh = 1. - sh;
    double sigma = 8.*umsh*umsh/sh;
    double chi_1 = 4.*umsh*(5.*sh + 3.)/3./sh; 
    double chi_2 = 4.*(3.*sh*sh + 2.*sh - 9.)/sh;
    double deltaMb2 = (lambda_1*chi_1 + lambda_2*chi_2)/Mb_pole/Mb_pole; // pole mass as stated in hep-ph/9801456
    
    switch(order)
    {
        case LO:
            return sigma + deltaMb2;
        case NLO:
            return sigma*8.*alstilde*omega77_T(sh);
        default:
            throw std::runtime_error("BXqll::S77_T: order not implemented");
    }
}

double BXqll::S79_T(double sh, orders order)
{
    double umsh = 1. - sh;
    double sigma = 8.*umsh*umsh;
    double chi_1 = 4.*umsh*umsh;
    double chi_2 = 4.*(9.*sh*sh - 6.*sh - 7.);
    double deltaMb2 = (lambda_1*chi_1 + lambda_2*chi_2)/Mb_pole/Mb_pole;
    
    switch(order)
    {
        case LO:
            return sigma + deltaMb2;
        case NLO:
            return sigma*8.*alstilde*omega79_T(sh);
        default:
            throw std::runtime_error("BXqll::S79_T: order not implemented");
    }
}

double BXqll::S99_T(double sh, orders order)
{
    double umsh = 1. - sh;
    double sigma = 2.*sh*umsh*umsh;
    double chi_1 = -sh*umsh*(3.*sh + 5.)/3.;
    double chi_2 = sh*(15.*sh*sh - 14.*sh - 5.);
    double deltaMb2 = (lambda_1*chi_1 + lambda_2*chi_2)/Mb_pole/Mb_pole;
    double omega_2 = mySM.Beta0(5.)*log(muh)*omega99_T(sh) + 54.919*umsh*umsh*umsh*umsh -
                    136.374*umsh*umsh*umsh + 119.344*umsh*umsh - 15.6175*umsh - 31.1706;
    
    switch(order)
    {
        case LO:
            return sigma + deltaMb2;
        case NLO:
            return sigma*8.*alstilde*omega99_T(sh);
        case NNLO:
            return sigma*16.*alstilde*alstilde*omega_2;
        default:
            throw std::runtime_error("BXqll::S99_T: order not implemented");
    }
}

double BXqll::S1010_T(double sh, orders order)
{
    return S99_T(sh,order);
}

double BXqll::S77_L(double sh, orders order)
{
    double umsh = 1. - sh;
    double sigma = 4.*umsh*umsh;
    double chi_1 = -2.*umsh*(3.*sh + 13.)/3.;
    double chi_2 = 2.*(15.*sh*sh - 6.*sh - 13.);
    double deltaMb2 = (lambda_1*chi_1 + lambda_2*chi_2)/Mb_pole/Mb_pole;
    
    switch(order)
    {
        case LO:
            return sigma + deltaMb2;
        case NLO:
            return sigma*8.*alstilde*omega77_L(sh);
        default:
            throw std::runtime_error("BXqll::S77_L: order not implemented");
    }
}

double BXqll::S79_L(double sh, orders order)
{
    double umsh = 1. - sh;
    double sigma = 4.*umsh*umsh;
    double chi_1 = 2.*umsh*umsh;
    double chi_2 = 2.*(3.*sh*sh - 6.*sh - 1.);
    double deltaMb2 = (lambda_1*chi_1 + lambda_2*chi_2)/Mb_pole/Mb_pole;
    
    switch(order)
    {
        case LO:
            return sigma + deltaMb2;
        case NLO:
            return sigma*8.*alstilde*omega79_L(sh);
        default:
            throw std::runtime_error("BXqll::S79_L: order not implemented");
    }
}

double BXqll::S99_L(double sh, orders order)
{
    double umsh = 1. - sh;
    double sigma = umsh*umsh;
    double chi_1 = umsh*(13.*sh + 3.)/6.;
    double chi_2 = (-17.*sh*sh + 10.*sh + 3.)/2.;
    double deltaMb2 = (lambda_1*chi_1 + lambda_2*chi_2)/Mb_pole/Mb_pole;
    double omega_2 = mySM.Beta0(5.)*log(muh)*omega99_L(sh) - 5.95974*umsh*umsh*umsh +
                    11.7493*umsh*umsh + 12.2293*umsh - 38.6457;
    
    switch(order)
    {
        case LO:
            return sigma + deltaMb2;
        case NLO:
            return sigma*8.*alstilde*omega99_L(sh);
        case NNLO:
            return sigma*16.*alstilde*alstilde*omega_2;
        default:
            throw std::runtime_error("BXqll::S99_L: order not implemented");
    }
}

double BXqll::S1010_L(double sh, orders order)
{
    return S99_L(sh,order);
}

double BXqll::S710_A(double sh, orders order)
{
    double umsh = 1. - sh;
    double sigma = -8.*umsh*umsh;
    double chi_1 = -4.*(3.*sh*sh + 2.*sh + 3.)/3.;
    double chi_2 = -4.*(9.*sh*sh - 10.*sh - 7.);
    double deltaMb2 = (lambda_1*chi_1 + lambda_2*chi_2)/Mb_pole/Mb_pole;
    
    switch(order)
    {
        case LO:
            return sigma + deltaMb2;
        case NLO:
            return sigma*8.*alstilde*omega710_A(sh);
        default:
            throw std::runtime_error("BXqll::S710_A: order not implemented");
    }
}

double BXqll::S910_A(double sh, orders order)
{
    double umsh = 1. - sh;
    double sigma = -4.*umsh*umsh;
    double chi_1 = -2.*sh*(3.*sh*sh + 2.*sh + 3.)/3.;
    double chi_2 = -2.*sh*(15.*sh*sh - 14.*sh - 9.);
    double deltaMb2 = (lambda_1*chi_1 + lambda_2*chi_2)/Mb_pole/Mb_pole;
    double omega_2 = mySM.Beta0(5.)*log(muh)*omega910_A(sh) + 74.3717*umsh*umsh*umsh*umsh -
                    183.885*umsh*umsh*umsh + 158.739*umsh*umsh - 29.0124*umsh - 30.8056;
    
    switch(order)
    {
        case LO:
            return sigma + deltaMb2;
        case NLO:
            return sigma*8.*alstilde*omega910_A(sh);
        case NNLO:
            return sigma*16.*alstilde*alstilde*omega_2;
        default:
            throw std::runtime_error("BXqll::S910_A: order not implemented");
    }
}

gslpp::complex BXqll::cij_T(unsigned int i, unsigned int j, double sh)
{
    unsigned int ij = 10*(i + 1) + (j + 1);
    double umsh = (1. - sh), uptsh = 1. + 3.*sh;
//    double r = sh * Mb_pole * Mb_pole / 4. / Mc / Mc;
    double r = sh * Mb_pole * Mb_pole / 4. / Mc_pole / Mc_pole; //FOR TESTING STAGE
    gslpp::complex Mj7, Mj9;
    
//    switch(j + 1)
//    {
//        case 9:
//            Mj7 = 0.;
//            Mj9 = (M_9[LO])(j);
//            break;
//        case 7:
//            Mj7 = (M_7[int_qed(NLO_QED11)])(j);
//            Mj9 = 0.;
//            break;
//        default:
//            Mj7 = 0.;
//            Mj9 = (M_9[int_qed(NLO_QED11)])(j);
//            break;
//    }
    
    //FOR TESTING STAGE
    Mj7 = M_7[LO](j) + M_7[NLO](j) + M_7[NNLO](j) + M_7[int_qed(LO_QED)](j) +
                M_7[int_qed(NLO_QED11)](j) + M_7[int_qed(NLO_QED21)](j);
    Mj9 = M_9[LO](j) + M_9[NLO](j) + M_9[NNLO](j) + M_9[int_qed(LO_QED)](j) +
                M_9[int_qed(NLO_QED11)](j) + M_9[int_qed(NLO_QED21)](j);

    gslpp::complex F_M7c_M9c = F_BIR(r) * (Mj7.conjugate() / sh + Mj9.conjugate() / 2.);
    
    if (i + 1 == 2 && j + 1 > 1)
//        return (-aletilde*8.*lambda_2/9./Mc/Mc*umsh*umsh*uptsh * F_M7c_M9c);
        return (-aletilde*8.*lambda_2/9./Mc_pole/Mc_pole*umsh*umsh*uptsh * F_M7c_M9c); //FOR TESTING STAGE
    else if (ij == 11)
//        return (aletilde*4.*lambda_2/27./Mc/Mc*umsh*umsh*uptsh * F_M7c_M9c);
        return (aletilde*4.*lambda_2/27./Mc_pole/Mc_pole*umsh*umsh*uptsh * F_M7c_M9c); //FOR TESTING STAGE
    else if (ij == 12)
    {
//        return (-aletilde*8.*lambda_2/9./Mc/Mc*umsh*umsh*uptsh *
//                (F_BIR(r).conjugate() * (M_9[int_qed(NLO_QED11)])(0) / 2. - F_M7c_M9c / 6.));
        
        //FOR TESTING STAGE
        Mj7 = M_7[LO](0) + M_7[NLO](0) + M_7[NNLO](0) + M_7[int_qed(LO_QED)](0) +
                M_7[int_qed(NLO_QED11)](0) + M_7[int_qed(NLO_QED21)](0);
        Mj9 = M_9[LO](0) + M_9[NLO](0) + M_9[NNLO](0) + M_9[int_qed(LO_QED)](0) +
                M_9[int_qed(NLO_QED11)](0) + M_9[int_qed(NLO_QED21)](0);
        
        return (-aletilde*8.*lambda_2/9./Mc_pole/Mc_pole*umsh*umsh*uptsh *
                (F_BIR(r).conjugate() * (Mj7 / sh + Mj9 / 2.) - F_M7c_M9c / 6.));
    }
    else if (i + 1 == 1)
//        return (aletilde*8.*lambda_2/54./Mc/Mc*umsh*umsh*uptsh * F_M7c_M9c);
        return (aletilde*8.*lambda_2/54./Mc_pole/Mc_pole*umsh*umsh*uptsh * F_M7c_M9c); //FOR TESTING STAGE
    else
        return (0.);
}

gslpp::complex BXqll::cij_L(unsigned int i, unsigned int j, double sh)
{
    unsigned int ij = 10*(i + 1) + (j + 1);
    double umsh = (1. - sh), tmsh = 3. - sh;
//    double r = sh * Mb_pole * Mb_pole / 4. / Mc / Mc;
    double r = sh * Mb_pole * Mb_pole / 4. / Mc_pole / Mc_pole; //FOR TESTING STAGE
    gslpp::complex Mj7, Mj9;
    
//    switch(j + 1)
//    {
//        case 9:
//            Mj7 = 0.;
//            Mj9 = (M_9[LO])(j);
//            break;
//        case 7:
//            Mj7 = (M_7[int_qed(NLO_QED11)])(j);
//            Mj9 = 0.;
//            break;
//        default:
//            Mj7 = 0.;
//            Mj9 = (M_9[int_qed(NLO_QED11)])(j);
//            break;
//    }
    
    //FOR TESTING STAGE
    Mj7 = M_7[LO](j) + M_7[NLO](j) + M_7[NNLO](j) + M_7[int_qed(LO_QED)](j) +
                M_7[int_qed(NLO_QED11)](j) + M_7[int_qed(NLO_QED21)](j);
    Mj9 = M_9[LO](j) + M_9[NLO](j) + M_9[NNLO](j) + M_9[int_qed(LO_QED)](j) +
                M_9[int_qed(NLO_QED11)](j) + M_9[int_qed(NLO_QED21)](j);
    
    gslpp::complex F_M7c_M9c = F_BIR(r) * (Mj7.conjugate() + Mj9.conjugate() / 2.);
    
    if (i + 1 == 2 && j + 1 > 1)
//        return (-aletilde*8.*lambda_2/9./Mc/Mc*umsh*umsh*tmsh * F_M7c_M9c);
        return (-aletilde*8.*lambda_2/9./Mc_pole/Mc_pole*umsh*umsh*tmsh * F_M7c_M9c); //FOR TESTING STAGE
    else if (ij == 11)
//        return (aletilde*4.*lambda_2/27./Mc/Mc*umsh*umsh*tmsh * F_M7c_M9c);        
        return (aletilde*4.*lambda_2/27./Mc_pole/Mc_pole*umsh*umsh*tmsh * F_M7c_M9c); //FOR TESTING STAGE
    else if (ij == 12)
    {
//        return (-aletilde*8.*lambda_2/9./Mc/Mc*umsh*umsh*tmsh *
//                (F_BIR(r).conjugate() * (M_9[int_qed(NLO_QED11)])(0) / 2. - F_M7c_M9c / 6.));
        
        //FOR TESTING STAGE
        Mj7 = M_7[LO](0) + M_7[NLO](0) + M_7[NNLO](0) + M_7[int_qed(LO_QED)](0) +
                M_7[int_qed(NLO_QED11)](0) + M_7[int_qed(NLO_QED21)](0);
        Mj9 = M_9[LO](0) + M_9[NLO](0) + M_9[NNLO](0) + M_9[int_qed(LO_QED)](0) +
                M_9[int_qed(NLO_QED11)](0) + M_9[int_qed(NLO_QED21)](0);
        
        return (-aletilde*8.*lambda_2/9./Mc_pole/Mc_pole*umsh*umsh*tmsh *
                (F_BIR(r).conjugate() * (Mj7 + Mj9 / 2.) - F_M7c_M9c / 6.));
    }
    else if (i + 1 == 1)
//        return (aletilde*8.*lambda_2/54./Mc/Mc*umsh*umsh*tmsh * F_M7c_M9c);
        return (aletilde*8.*lambda_2/54./Mc_pole/Mc_pole*umsh*umsh*tmsh * F_M7c_M9c); //FOR TESTING STAGE
    else
        return (0.);
}

gslpp::complex BXqll::cij_A(unsigned int i, unsigned int j, double sh)
{
    unsigned int ij = 100*(i + 1) + (j + 1);
    double umsh = (1. - sh), uptsh = 1. + 3.*sh;
//    double r = sh * Mb_pole * Mb_pole / 4. / Mc / Mc;
    double r = sh * Mb_pole * Mb_pole / 4. / Mc_pole / Mc_pole; //FOR TESTING STAGE
    
    switch(ij)
    {
        case 110:
//            return (aletilde*4.*lambda_2/9./Mc/Mc*umsh*umsh*uptsh * F_BIR(r));
            return (aletilde*4.*lambda_2/9./Mc_pole/Mc_pole*umsh*umsh*uptsh * F_BIR(r)); //FOR TESTING STAGE
        case 210:
//            return (-aletilde*4.*lambda_2/54./Mc/Mc*umsh*umsh*uptsh * F_BIR(r));
            return (-aletilde*4.*lambda_2/54./Mc_pole/Mc_pole*umsh*umsh*uptsh * F_BIR(r)); //FOR TESTING STAGE
        default:
            return (0.);
    }
}

gslpp::complex BXqll::eij_T(unsigned int i, unsigned int j, double sh)
{
    unsigned int ij;
    double umsh = (1. - sh);
    double sigma77_T = 8.*umsh*umsh/sh, sigma79_T = 8.*umsh*umsh, sigma99_T = 2.*sh*umsh*umsh;
    
    if (j == 9)
        ij = 100*(i + 1) + (j + 1);
    else
        ij = 10*(i + 1) + (j + 1);
    
    switch(ij)
    {
        case 11:
            return (128./9. * aletilde*aletilde*aletilde * sigma99_T*omega22em_T(sh));
        case 12:
            return (64./3. * aletilde*aletilde*aletilde * sigma99_T*omega22em_T(sh));
        case 17:
            return (32./3. * aletilde*aletilde*aletilde * sigma79_T*omega27em_T(sh));
        case 19:
            return (32./3. * aletilde*aletilde * sigma99_T*omega29em_T(sh));
        case 22:
            return (8. * aletilde*aletilde*aletilde * sigma99_T*omega22em_T(sh));
        case 27:
            return (8. * aletilde*aletilde*aletilde * sigma79_T*omega27em_T(sh));
        case 29:
            return (8. * aletilde*aletilde * sigma99_T*omega29em_T(sh));
        case 77:
            return (8. * aletilde*aletilde*aletilde * sigma77_T*omega77em_T(sh));
        case 79:
            return (8. * aletilde*aletilde * sigma79_T*omega79em_T(sh));
        case 99:
            return (8. * aletilde * sigma99_T*omega99em_T(sh));
        case 1010:
            return (8. * aletilde * sigma99_T*omega99em_T(sh));
        default:
            return (0.);
    }
}

gslpp::complex BXqll::eij_L(unsigned int i, unsigned int j, double sh)
{
    unsigned int ij;
    double umsh = (1. - sh);
    double sigma77_L = 4.*umsh*umsh, sigma79_L = 4.*umsh*umsh, sigma99_L = umsh*umsh;
    
    if (j == 9)
        ij = 100* (i + 1) + (j + 1);
    else
        ij = 10*(i + 1) + (j + 1);
    
    switch(ij)
    {
        case 11:
            return (128./9. * aletilde*aletilde*aletilde * sigma99_L*omega22em_L(sh));
        case 12:
            return (64./3. * aletilde*aletilde*aletilde * sigma99_L*omega22em_L(sh));
        case 17:
            return (32./3. * aletilde*aletilde*aletilde * sigma79_L*omega27em_L(sh));
        case 19:
            return (32./3. * aletilde*aletilde * sigma99_L*omega29em_L(sh));
        case 22:
            return (8. * aletilde*aletilde*aletilde * sigma99_L*omega22em_L(sh));
        case 27:
            return (8. * aletilde*aletilde*aletilde * sigma79_L*omega27em_L(sh));
        case 29:
            return (8. * aletilde*aletilde * sigma99_L*omega29em_L(sh));
        case 77:
            return (8. * aletilde*aletilde*aletilde * sigma77_L*omega77em_L(sh));
        case 79:
            return (8. * aletilde*aletilde * sigma79_L*omega79em_L(sh));
        case 99:
            return (8. * aletilde * sigma99_L*omega99em_L(sh));
        case 1010:
            return (8. * aletilde * sigma99_L*omega99em_L(sh));
        default:
            return (0.);
    }
}

gslpp::complex BXqll::eij_A(unsigned int i, unsigned int j, double sh)
{
    unsigned int ij;
    double umsh = (1. - sh);
    double sigma710_A = -8.*umsh*umsh, sigma910_A = -4.*sh*umsh*umsh;
    
    if (j == 9)
        ij = 100*(i + 1) + (j + 1);
    else
        ij = 10*(i + 1) + (j + 1);
    
    switch(ij)
    {
        case 110:
            return (32./3. * aletilde*aletilde * sigma910_A*omega210em_A(sh));
        case 210:
            return (8. * aletilde*aletilde * sigma910_A*omega210em_A(sh));
        case 710:
            return (8. * aletilde*aletilde * sigma710_A*omega710em_A(sh));
        case 910:
            return (8. * aletilde * sigma910_A*omega910em_A(sh));        
        default:
            return (0.);
    }
}

double BXqll::omega77_T(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    double dilog_umsh = gslpp_special_functions::dilog(umsh);
    double dilog_umsqrt = gslpp_special_functions::dilog(umsqrt);
    
    return (-8./3.*log(muh) - (sqrt(sh)+1.)*(sqrt(sh)+1.)*(pow(sh,1.5)-10.*sh+13.*sqrt(sh)-8.)*
            dilog_umsh/6./umsh/umsh + 2.*sqrt(sh)*(sh*sh-6.*sh-3.)*dilog_umsqrt/3./umsh/umsh -
            M_PI*M_PI*(3.*pow(sh,1.5)+22.*sh+23.*sqrt(sh)+16.)*umsqrt*umsqrt/36./umsh/umsh +
            (5.*sh*sh*sh-54.*sh*sh+57.*sh-8.)/18./umsh/umsh - log(umsh) +
            sh*(5.*sh+1.)*log(sh)/3./umsh/umsh + 2./3.*log(umsh)*log(sh));
}

double BXqll::omega79_T(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    double dilog_umsh = gslpp_special_functions::dilog(umsh);
    double dilog_umsqrt = gslpp_special_functions::dilog(umsqrt);
    
    return (-4./3.*log(muh) - 2.*sqrt(sh)*(sh+3.)*dilog_umsqrt/3./umsh/umsh - M_PI*M_PI*
            (16.*sh+29.*sqrt(sh)+19.)*umsqrt*umsqrt/36./umsh/umsh + (sh*sh-6.*sh+5.)/6./umsh/umsh +
            (sqrt(sh)+1.)*(sqrt(sh)+1.)*(8.*sh-15.*sqrt(sh)+9.)*dilog_umsh/6./umsh/umsh - (5.*sh+1.)*
            log(umsh)/6./sh + sh*(3.*sh+1.)*log(sh)/6./umsh/umsh + 2./3.*log(umsh)*log(sh));
}

double BXqll::omega99_T(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    double dilog_umsh = gslpp_special_functions::dilog(umsh);
    double dilog_umsqrt = gslpp_special_functions::dilog(umsqrt);
    
    return ((sqrt(sh)+1.)*(sqrt(sh)+1.)*(8.*pow(sh,1.5)-15.*sh+4.*sqrt(sh)-5.)*dilog_umsh/
            6./umsh/umsh/sqrt(sh) - 2.*(sh*sh-12.*sh-5.)*dilog_umsqrt/3./umsh/umsh/sqrt(sh) -
            M_PI*M_PI*(16.*pow(sh,1.5)+29.*sh+4.*sqrt(sh)+15.)*umsqrt*umsqrt/36./umsh/umsh/sqrt(sh) +
            (2.*sh*sh-7.*sh-5.)*log(sh)/3./umsh/umsh + (sh*sh+18.*sh-19.)/6./umsh/umsh -
            (2.*sh+1)*log(umsh)/3./sh + 2./3.*log(umsh)*log(sh));
}

double BXqll::omega77_L(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    double dilog_umsh = gslpp_special_functions::dilog(umsh);
    double dilog_umsqrt = gslpp_special_functions::dilog(umsqrt);
    
    return (-8./3.*log(muh) + (sqrt(sh)+1.)*(sqrt(sh)+1.)*(4.*pow(sh,1.5)-7.*sh+2.*sqrt(sh)-3.)*
            dilog_umsh/3./umsh/umsh/sqrt(sh) - (9.*sh*sh-38.*sh+29.)/6./umsh/umsh -
            4.*(sh*sh-6.*sh-3.)*dilog_umsqrt/3./umsh/umsh/sqrt(sh) - M_PI*M_PI*
            (8.*pow(sh,1.5)+13.*sh+2.*sqrt(sh)+9.)*umsqrt*umsqrt/18./umsh/umsh/sqrt(sh) - (sh*sh*sh-3.*sh+2.)*
            log(umsh)/3./umsh/umsh/sh + 2.*(sh*sh-3.*sh-3.)*log(sh)/3./umsh/umsh + 2./3.*log(umsh)*log(sh));
}

double BXqll::omega79_L(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    double dilog_umsh = gslpp_special_functions::dilog(umsh);
    double dilog_umsqrt = gslpp_special_functions::dilog(umsqrt);
    
    return (-4./3.*log(muh) + 4.*sqrt(sh)*(sh+3.)*dilog_umsqrt/3./umsh/umsh + (sqrt(sh)+1.)*
            (sqrt(sh)+1.)*(4.*sh-9.*sqrt(sh)+3.)*dilog_umsh/3./umsh/umsh + (7.*sh*sh-2.*sh-5.)/
            6./umsh/umsh - M_PI*M_PI*(8.*sh+19.*sqrt(sh)+5.)*umsqrt*umsqrt/18./umsh/umsh -
            (2.*sh+1.)*log(umsh)/3./sh + (sh-7.)*sh*log(sh)/3./umsh/umsh + 2./3.*log(umsh)*log(sh));
}

double BXqll::omega99_L(double sh)
{
    double umsh = 1.-sh;
    double umsqrt = 1.-sqrt(sh);
    double dilog_umsh = gslpp_special_functions::dilog(umsh);
    double dilog_umsqrt = gslpp_special_functions::dilog(umsqrt);
    
    return (-(sqrt(sh)+1.)*(sqrt(sh)+1.)*(pow(sh,1.5)-8.*sh+3.*sqrt(sh)-4.)*dilog_umsh/
            3./umsh/umsh + 4.*sqrt(sh)*(sh*sh-12.*sh-5.)*dilog_umsqrt/3./umsh/umsh -
            M_PI*M_PI*(3.*pow(sh,1.5)+20.*sh+sqrt(sh)+8.)*umsqrt*umsqrt/18./umsh/umsh +
            (4.*sh*sh*sh-51.*sh*sh+42.*sh+5.)/6./umsh/umsh - log(umsh) +
            8.*sh*(2.*sh+1.)*log(sh)/3./umsh/umsh + + 2./3.*log(umsh)*log(sh));
}

double BXqll::omega710_A(double sh)
{
    double umsh = 1.-sh;
    double num = 3.*umsh*umsh;
    double umsqrt = 1.-sqrt(sh);
    double dilog_umsh = gslpp_special_functions::dilog(umsh);
    double dilog_umsqrt = gslpp_special_functions::dilog(umsqrt);
    
    return (-4./3.*log(muh) + 2.*(4.*sh*sh-13.*sh-1.)*dilog_umsqrt/num - (2.*sh*sh-9.*sh-3.)*dilog_umsh/num -
            (3.*sh*sh-16.*sh+13.)*log(umsqrt)/num + (4.*sh*sh-13.*sh-1.)*log(umsqrt)*log(sh)/num -
            (2.*sh*sh-9.*sh-3.)*log(umsh)*log(sh)/num + (sh*sh*sh-23.*sh*sh+23.*sh-1.)*log(umsh)/2./num/sh +
            (sh-20.*sqrt(sh)+5.)*umsqrt*umsqrt/2./num - M_PI*M_PI/3.);
}

double BXqll::omega910_A(double sh)
{
    double umsh = 1.-sh;
    double num = 3.*umsh*umsh;
    double umsqrt = 1.-sqrt(sh);
    double dilog_umsh = gslpp_special_functions::dilog(umsh);
    double dilog_umsqrt = gslpp_special_functions::dilog(umsqrt);
    
    return (-2.*(sh*sh-3.*sh-1.)*dilog_umsh/num - 4.*(5.-2.*sh)*sh*dilog_umsqrt/num -
            (4.*sqrt(sh)-3.)*umsqrt*umsqrt/num - 2.*(2.*sh*sh-7.*sh+5.)*log(umsqrt)/num -
            2.*(sh*sh-3.*sh-1.)*log(umsh)*log(sh)/num + (2.*sh*sh*sh-11.*sh*sh+10.*sh-1.)*
            log(umsh)/num/sh + 2.*sh*(2.*sh-5.)*log(umsqrt)*log(sh)/num - M_PI*M_PI/3.);
}

double BXqll::omega77em_T(double sh)
{
    double umsh = 1. - sh;
    
    return (Lbl*(1.54986 - 1703.72*sh*sh*sh*sh*sh + 1653.38*sh*sh*sh*sh - 683.608*sh*sh*sh +
            179.279*sh*sh - 35.5047*sh)/8./umsh/umsh);
}

double BXqll::omega79em_T(double sh)
{
    double umsh = 1. - sh;
    
    return (Lbl*(19.063 + 2158.03*sh*sh*sh*sh - 2062.92*sh*sh*sh + 830.53*sh*sh -
            186.12*sh + 0.324236/sh)/8./umsh/umsh);
}

double BXqll::omega99em_T(double sh)
{
    double umsh = 1. - sh;
    
    return (Lbl*(2.2596 + 157.984*sh*sh*sh*sh - 141.281*sh*sh*sh + 52.8914*sh*sh -
            13.5377*sh + 0.0284049/sh)/2./sh/umsh/umsh);
}

double BXqll::omega22em_T(double sh)
{
    double umsh = 1. - sh;
    double Lmub = log(mu_b/5.);
    
    return (Lbl*((2.84257 + 269.974*sh*sh*sh*sh - 194.443*sh*sh*sh + 48.4535*sh*sh -
            8.24929*sh + 0.0111118/sh)/2./sh/umsh/umsh +
            Lmub*4.*(4.54727 + 330.182*sh*sh*sh*sh - 258.194*sh*sh*sh + 79.8713*sh*sh -
            19.6855*sh +  0.0371348/sh)/9./sh/umsh/umsh) +
            64./81.*omega99em_T(sh)*Lmub*Lmub);
}

gslpp::complex BXqll::omega27em_T(double sh)
{
    double umsh = 1. - sh;
    double Lmub = log(mu_b/5.);
    gslpp::complex i = gslpp::complex::i();
    
    return (Lbl*((21.5291 + 3044.94*sh*sh*sh*sh - 2563.05*sh*sh*sh + 874.074*sh*sh -
            175.874*sh + 0.121398/sh)/8./umsh/umsh +
            i*(2.49475 + 598.376*sh*sh*sh*sh - 456.831*sh*sh*sh + 117.683*sh*sh -
            9.90525*sh - 0.0116501/sh)/8./umsh/umsh) +
            8./9.*omega79em_T(sh)*Lmub);
}

gslpp::complex BXqll::omega29em_T(double sh)
{
    double umsh = 1. - sh;
    double Lmub = log(mu_b/5.);
    gslpp::complex i = gslpp::complex::i();
    
    return (Lbl*((4.54727 + 330.182*sh*sh*sh*sh - 258.194*sh*sh*sh + 79.8713*sh*sh -
            19.6855*sh + 0.0371348/sh)/2./sh/umsh/umsh +
            i*(73.9149*sh*sh*sh*sh - 61.1338*sh*sh*sh + 14.6517*sh*sh - 0.102331*sh +
            0.710037)/2./sh/umsh/umsh) +
            16./9.*omega99em_T(sh)*Lmub);
}

double BXqll::omega77em_L(double sh)
{
    double umsh = 1. - sh;
    
    return (Lbl*(9.73761 + 647.747*sh*sh*sh*sh - 642.637*sh*sh*sh + 276.839*sh*sh -
            68.3562*sh - 1.6755/sh)/4./umsh/umsh);
}

double BXqll::omega79em_L(double sh)
{
    double umsh = 1. - sh;
    
    return (Lbl*(-6.03641 - 896.643*sh*sh*sh*sh + 807.349*sh*sh*sh - 278.559*sh*sh +
            47.6636*sh - 0.190701/sh)/4./umsh/umsh);
}

double BXqll::omega99em_L(double sh)
{
    double umsh = 1. - sh;
    
    return (Lbl*(-0.768521 - 80.8068*sh*sh*sh*sh + 70.0821*sh*sh*sh - 21.2787*sh*sh +
            2.9335*sh - 0.0180809/sh)/umsh/umsh);
}

double BXqll::omega22em_L(double sh)
{
    double umsh = 1. - sh;
    double Lmub = log(mu_b/5.);
    
    return (Lbl*((-1.71832 - 234.11*sh*sh*sh*sh + 162.126*sh*sh*sh - 37.2361*sh*sh +
            6.29949*sh - 0.00810233/sh)/umsh/umsh +
            Lmub*8.*(224.662*sh*sh*sh - 2.27221 - 298.369*sh*sh*sh*sh - 65.1375*sh*sh +
            11.5686*sh - 0.0233098/sh)/9./umsh/umsh) +
            64./81.*omega99em_L(sh)*Lmub*Lmub);
}

gslpp::complex BXqll::omega27em_L(double sh)
{
    double umsh = 1. - sh;
    double Lmub = log(mu_b/5.);
    gslpp::complex i = gslpp::complex::i();
    
    return (Lbl*((-8.01684 - 1121.13*sh*sh*sh*sh + 882.711*sh*sh*sh - 280.866*sh*sh +
            54.1943*sh - 0.128988/sh)/4./umsh/umsh +
            i*(-2.14058 - 588.771*sh*sh*sh*sh + 483.997*sh*sh*sh - 124.579*sh*sh +
            12.3282*sh + 0.0145059/sh)/4./umsh/umsh) +
            8./9.*omega79em_L(sh)*Lmub);
}

gslpp::complex BXqll::omega29em_L(double sh)
{
    double umsh = 1. - sh;
    double Lmub = log(mu_b/5.);
    gslpp::complex i = gslpp::complex::i();
    
    return (Lbl*((-2.27221 - 298.369*sh*sh*sh*sh + 224.662*sh*sh*sh - 65.1375*sh*sh +
            11.5686*sh - 0.0233098/sh)/umsh/umsh +
            i*(-0.666157 - 120.303*sh*sh*sh*sh + 109.315*sh*sh*sh - 28.2734*sh*sh +
            2.44527*sh + 0.00279781/sh)/umsh/umsh) +
            16./9.*omega99em_L(sh)*Lmub);
}

double BXqll::omega710em_A(double sh)
{
    double umsh = 1. - sh;
    
    return (Lbl*((7. - 16.*sqrt(sh) + 9.*sh)/4./umsh + log(1. - sqrt(sh)) +
            (1. + 3.*sh)/umsh*log((1. + sqrt(sh))/2.) - sh*log(sh)/umsh));
}

double BXqll::omega910em_A(double sh)
{
    double umsh = 1. - sh;
    
    return (Lbl*(log(1. - sqrt(sh)) - (5. - 16.*sqrt(sh) + 11.*sh)/4./umsh +
            (1. - 5.*sh)/umsh*log((1. + sqrt(sh))/2.) - (1. - 3.*sh)*log(sh)/umsh));
}

gslpp::complex BXqll::omega210em_A(double sh)
{
    double umsh = 1. - sh;
    double a = 16.*Mc*Mc*Mc*Mc/Mb/Mb/Mb/Mb;
    double shma = sh - a;
    double Lmub = log(mu_b/5.);
    gslpp::complex omega, i = gslpp::complex::i();
    
    omega = Lbl*((-351.322*sh*sh*sh*sh + 378.173*sh*sh*sh - 160.158*sh*sh + 24.2096*sh +
            0.305176)/24./sh/umsh/umsh) +
            8./9.*omega910em_A(sh)*Lmub;
    if(sh > a)
        omega += Lbl*i*(7.98625 + 238.507*shma - 766.869*shma*shma)*shma/24./sh/umsh/umsh;
    
    return omega;
}

gslpp::complex BXqll::f_Huber(double sh, double gamma_9, double rho_c, double rho_b, double rho_0, double rho_num)
{
    gslpp::complex i = gslpp::complex::i();
    
//    return (gamma_9 * log(Mb/mu_b) + rho_c * (g_Huber(4.*Mc_pole*Mc_pole/Mb_pole/Mb_pole/sh) +
//            8./9. * log(Mb/Mc)) + rho_b * g_Huber(4./sh) +
//            rho_0 * (log(sh) - i*M_PI) + rho_num);
    
    return (-gamma_9 * log(muh) + rho_c * (g_Huber(4.*z/sh) - 4./9. * log(z) + KS_cc(sh)) +
            rho_b * g_Huber(4./sh) + rho_0 * (log(sh) - i*M_PI) + rho_num);
}

gslpp::complex BXqll::f9pen_Huber(double sh)
{
    gslpp::complex i = gslpp::complex::i();
    
//    return (8. * log(Mb/mu_b) - 3. * g_Huber(4.*Mtau*Mtau/Mb_pole/Mb_pole/sh) - 8./3. * log(Mb/Mtau) +
//            8./3. * (log(sh) - i*M_PI) - 40./9.);
    
    return (-8. * log(muh) - 3. * g_Huber(4.*Mtau*Mtau/Mb_pole/Mb_pole/sh) - 8./3. * log(Mb_pole/Mtau) +
            8./3. * (log(sh) - i*M_PI) - 40./9.);
}

gslpp::complex BXqll::g_Huber(double y)
{
    gslpp::complex i = gslpp::complex::i();
    gslpp::complex g_y;
    
    g_y = -2./9.*(2.+y)*sqrt(fabs(1.-y));
    if(y < 1.)
        g_y *= log(fabs((1.+sqrt(1.-y))/(1.-sqrt(1.-y))))-i*M_PI;
    else
        g_y *= 2.*atan(1./sqrt(y-1.));
        
    g_y += 20./27. + 4./9.*y;
    
    return (g_y);
}

gslpp::complex BXqll::KS_cc(double sh)
{
    double pre = 3. * sh / ale / ale;
    gslpp::complex i = gslpp::complex::i();
    gslpp::complex kscc;
    
    // Contributions of J/Psi, Psi(2S)
    // Masses, total decay widths, and branching fractions taken from Kruger:1996cv
    kscc  = KS_aux(sh, 3.097, 0.088e-3, 6.0e-2, 0.877);
    kscc += KS_aux(sh, 3.686, 0.277e-3, 8.3e-3, 0.9785);
    
    kscc  = pre * kscc;
    
    // Contributions from the continuum for the large-q^2 region
    // 0.6 < sh < 0.69 region divergent, not implemented
    if (sh > 0.69)
        kscc += (i * M_PI * 1.02 - log(1. - sh / 4. / z)) / 3.;
    
    return (kscc);
}

gslpp::complex BXqll::KS_aux(double sh, double m, double Gamma, double Br_ll, double Br_had)
{
    gslpp::complex i = gslpp::complex::i();
    double Gamma_had = Br_had * Gamma / Mb_pole;
    double a = m * m / Mb_pole / Mb_pole;
    double b = a * Gamma * Gamma / Mb_pole / Mb_pole;
    double sqrtb = sqrt(b);
    double amsh = a - sh;
    double am4z = a - 4. * z;
    
    return (Br_ll * Gamma / Mb_pole * Gamma_had * ((M_PI * amsh + 2. * amsh * atan(am4z / sqrtb) +
            sqrtb * (log(b + am4z * am4z) - 2. * log(4. * z - sh))) /
            (2. * sqrtb * (b + amsh * amsh)) + i * M_PI / (amsh * amsh + b)));
}

gslpp::complex BXqll::F_BIR(double r)
{
    gslpp::complex i = gslpp::complex::i();

    if(r > 0. && r < 1.)
        return (1.5 / r * (atan(sqrt(r / (1. - r))) / sqrt(r - r*r) - 1.));
    else if (r > 1.)
        return (1.5 / r * ((log((1. - sqrt(1. - 1./r)) / (1. + sqrt(1. - 1./r))) + i * M_PI) /
                2. / sqrt(r*r - r) - 1.));
    else
        throw std::runtime_error("BXqll::F_BIR(): 1/mc^2 corrections diverge at q^2 = 4*mc^2");
}

double BXqll::Phi_u(orders ord)
{
    switch(ord)
    {
        case LO:
            return(1.);
        case NLO:
            return(alstilde * phi1);
        case NNLO:
            return(alstilde * alstilde * (phi2 + 2. * mySM.Beta0(5) * phi1 * log(muh))
                   + (lambda_1 / 2. - 9. / 2. * lambda_2) / Mb / Mb );
        case FULLNNLO:
            return (1. + alstilde * phi1 +
                    alstilde * alstilde * (phi2 + 2. * mySM.Beta0(5) * phi1 * log(muh)) +
                    (lambda_1 / 2. - 9. / 2. * lambda_2) / Mb / Mb);
        default:
            throw std::runtime_error("BXqll::Phi_u(): order not implemented.");
     }
}

double BXqll::Phi_u(orders_qed ord_qed)
{
    switch(ord_qed)
    {
        case LO_QED:
            return(kappa * (12. / 23. * (1. - alsmu / mySM.Als(mySM.getMuw(), FULLNNNLO, true))));
        case FULLNLO_QED:
            return (Phi_u(FULLNNLO) +
                    kappa * (12. / 23. * (1. - alsmu / mySM.Als(mySM.getMuw(), FULLNNNLO, true))));
        default:
            throw std::runtime_error("BXqll::Phi_u(): order not implemented.");
     }
}

unsigned int BXqll::int_qed(orders_qed order_qed)
{
    // For LO_QED to come right after NNLO
    return (order_qed - NO_QED + NNLO);
}

double BXqll::CCH_multiplication(std::vector< gslpp::matrix<gslpp::complex> >& Hij)
{
    double Phi_ll = 0.;
    
    for(unsigned int j = 0; j < 15; j++)
    {
        for(unsigned int i = 0; i <= j; i++)
        {
            // 00, with 1/mc^2 corrections
            Phi_ll += (WC(i, LO).conjugate() * WC(j, LO) * Hij[LO](i,j) ).real();
            
            // 10
            Phi_ll += (WC(i, LO).conjugate()  * WC(j, NLO) * Hij[LO](i,j) +
                       WC(i, NLO).conjugate() * WC(j, LO)  * Hij[LO](i,j) +
                       WC(i, LO).conjugate()  * WC(j, LO)  * Hij[NLO](i,j) ).real();
            
            // 20
            Phi_ll += (WC(i, LO).conjugate()   * WC(j, NNLO) * Hij[LO](i,j) +
                       WC(i, NLO).conjugate()  * WC(j, NLO)  * Hij[LO](i,j) +
                       WC(i, NNLO).conjugate() * WC(j, LO)   * Hij[LO](i,j) +
                       WC(i, LO).conjugate()   * WC(j, NLO)  * Hij[NLO](i,j) +
                       WC(i, NLO).conjugate()  * WC(j, LO)   * Hij[NLO](i,j) ).real();
            
            // 01
            Phi_ll += (WC(i, LO).conjugate()              * WC(j, int_qed(LO_QED)) * Hij[LO](i,j) +
                       WC(i, int_qed(LO_QED)).conjugate() * WC(j, LO)              * Hij[LO](i,j) ).real();
            
            // 11
            Phi_ll += (WC(i, LO).conjugate()                 * WC(j, int_qed(NLO_QED11)) * Hij[LO](i,j) +
                       WC(i, NLO).conjugate()                * WC(j, int_qed(LO_QED))    * Hij[LO](i,j) +
                       WC(i, int_qed(LO_QED)).conjugate()    * WC(j, NLO)                * Hij[LO](i,j) +
                       WC(i, int_qed(NLO_QED11)).conjugate() * WC(j, LO)                 * Hij[LO](i,j) +
                       WC(i, LO).conjugate()                 * WC(j, int_qed(LO_QED))    * Hij[NLO](i,j) +
                       WC(i, int_qed(LO_QED)).conjugate()    * WC(j, LO)                 * Hij[NLO](i,j) +
                       WC(i, LO).conjugate()                 * WC(j, LO)                 * Hij[int_qed(NLO_QED11)](i,j) ).real();
            
            // 21
            Phi_ll += (WC(i, LO).conjugate()                 * WC(j, int_qed(NLO_QED21)) * Hij[LO](i,j) +
                       WC(i, NLO).conjugate()                * WC(j, int_qed(NLO_QED11)) * Hij[LO](i,j) +
                       WC(i, NNLO).conjugate()               * WC(j, int_qed(LO_QED))    * Hij[LO](i,j) +
                       WC(i, int_qed(LO_QED)).conjugate()    * WC(j, NNLO)               * Hij[LO](i,j) +
                       WC(i, int_qed(NLO_QED11)).conjugate() * WC(j, NLO)                * Hij[LO](i,j) +
                       WC(i, int_qed(NLO_QED21)).conjugate() * WC(j, LO)                 * Hij[LO](i,j) +
                       WC(i, LO).conjugate()                 * WC(j, int_qed(NLO_QED11)) * Hij[NLO](i,j) +                 
                       WC(i, NLO).conjugate()                * WC(j, int_qed(LO_QED))    * Hij[NLO](i,j) +
                       WC(i, int_qed(LO_QED)).conjugate()    * WC(j, NLO)                * Hij[NLO](i,j) +
                       WC(i, int_qed(NLO_QED11)).conjugate() * WC(j, LO)                 * Hij[NLO](i,j) +
                       WC(i, LO).conjugate()                 * WC(j, NLO)                * Hij[int_qed(NLO_QED11)](i,j) +
                       WC(i, NLO).conjugate()                * WC(j, LO)                 * Hij[int_qed(NLO_QED11)](i,j) +
                       WC(i, LO).conjugate()                 * WC(j, LO)                 * Hij[int_qed(NLO_QED21)](i,j) ).real();
            
            // 02
            Phi_ll += (WC(i, LO).conjugate()                 * WC(j, int_qed(NLO_QED02)) * Hij[LO](i,j) +
                       WC(i, int_qed(LO_QED)).conjugate()    * WC(j, int_qed(LO_QED))    * Hij[LO](i,j) +
                       WC(i, int_qed(NLO_QED02)).conjugate() * WC(j, LO)                 * Hij[LO](i,j) ).real();
            
            // 12
            Phi_ll += (WC(i, LO).conjugate()                 * WC(j, int_qed(NLO_QED12)) * Hij[LO](i,j) +
                       WC(i, NLO).conjugate()                * WC(j, int_qed(NLO_QED02)) * Hij[LO](i,j) +
                       WC(i, int_qed(LO_QED)).conjugate()    * WC(j, int_qed(NLO_QED11)) * Hij[LO](i,j) +
                       WC(i, int_qed(NLO_QED11)).conjugate() * WC(j, int_qed(LO_QED))    * Hij[LO](i,j) +
                       WC(i, int_qed(NLO_QED02)).conjugate() * WC(j, NLO)                * Hij[LO](i,j) +
                       WC(i, int_qed(NLO_QED12)).conjugate() * WC(j, LO)                 * Hij[LO](i,j) +
                       WC(i, LO).conjugate()                 * WC(j, int_qed(NLO_QED02)) * Hij[NLO](i,j) +
                       WC(i, int_qed(LO_QED)).conjugate()    * WC(j, int_qed(LO_QED))    * Hij[NLO](i,j) +
                       WC(i, int_qed(NLO_QED02)).conjugate() * WC(j, LO)                 * Hij[NLO](i,j) +
                       WC(i, LO).conjugate()                 * WC(j, int_qed(LO_QED))    * Hij[int_qed(NLO_QED11)](i,j) +
                       WC(i, int_qed(LO_QED)).conjugate()    * WC(j, LO)                 * Hij[int_qed(NLO_QED11)](i,j) ).real();
            
            // 22
            Phi_ll += (WC(i, LO).conjugate()                 * WC(j, int_qed(NLO_QED22)) * Hij[LO](i,j) +
                       WC(i, NLO).conjugate()                * WC(j, int_qed(NLO_QED12)) * Hij[LO](i,j) +
                       WC(i, NNLO).conjugate()               * WC(j, int_qed(NLO_QED02)) * Hij[LO](i,j) +
                       WC(i, int_qed(LO_QED)).conjugate()    * WC(j, int_qed(NLO_QED21)) * Hij[LO](i,j) +
                       WC(i, int_qed(NLO_QED11)).conjugate() * WC(j, int_qed(NLO_QED11)) * Hij[LO](i,j) +
                       WC(i, int_qed(NLO_QED21)).conjugate() * WC(j, int_qed(LO_QED))    * Hij[LO](i,j) +
                       WC(i, int_qed(NLO_QED02)).conjugate() * WC(j, NNLO)               * Hij[LO](i,j) +
                       WC(i, int_qed(NLO_QED12)).conjugate() * WC(j, NLO)                * Hij[LO](i,j) +
                       WC(i, int_qed(NLO_QED22)).conjugate() * WC(j, LO)                 * Hij[LO](i,j) +
                       WC(i, LO).conjugate()                 * WC(j, int_qed(NLO_QED12)) * Hij[NLO](i,j) +
                       WC(i, NLO).conjugate()                * WC(j, int_qed(NLO_QED02)) * Hij[NLO](i,j) +
                       WC(i, int_qed(LO_QED)).conjugate()    * WC(j, int_qed(NLO_QED11)) * Hij[NLO](i,j) +
                       WC(i, int_qed(NLO_QED11)).conjugate() * WC(j, int_qed(LO_QED))    * Hij[NLO](i,j) +
                       WC(i, int_qed(NLO_QED02)).conjugate() * WC(j, NLO)                * Hij[NLO](i,j) +
                       WC(i, int_qed(NLO_QED12)).conjugate() * WC(j, LO)                 * Hij[NLO](i,j) +
                       WC(i, LO).conjugate()                 * WC(j, int_qed(NLO_QED11)) * Hij[int_qed(NLO_QED11)](i,j) +
                       WC(i, NLO).conjugate()                * WC(j, int_qed(LO_QED))    * Hij[int_qed(NLO_QED11)](i,j) +
                       WC(i, int_qed(LO_QED)).conjugate()    * WC(j, NLO)                * Hij[int_qed(NLO_QED11)](i,j) +
                       WC(i, int_qed(NLO_QED11)).conjugate() * WC(j, LO)                 * Hij[int_qed(NLO_QED11)](i,j) +
                       WC(i, LO).conjugate()                 * WC(j, int_qed(LO_QED))    * Hij[int_qed(NLO_QED21)](i,j) +
                       WC(i, int_qed(LO_QED)).conjugate()    * WC(j, LO)                 * Hij[int_qed(NLO_QED21)](i,j) +
                       WC(i, LO).conjugate()                 * WC(j, LO)                 * Hij[int_qed(NLO_QED22)](i,j) ).real();
        }
    }
    
    // orders above 22 (only available for C9 and C10)
    // 32
    Phi_ll += (WC(9, int_qed(NLO_QED11)).conjugate() * WC(9, int_qed(NLO_QED21)) * Hij[LO](9,9) +
               WC(9, int_qed(NLO_QED21)).conjugate() * WC(9, int_qed(NLO_QED11)) * Hij[LO](9,9) +
               WC(9, int_qed(NLO_QED11)).conjugate() * WC(9, int_qed(NLO_QED11)) * Hij[NLO](9,9) ).real();
    
    // 03
    Phi_ll += (WC(8, int_qed(LO_QED)).conjugate()    * WC(8, int_qed(NLO_QED02)) * Hij[LO](8,8) +
               WC(8, int_qed(NLO_QED02)).conjugate() * WC(8, int_qed(LO_QED))    * Hij[LO](8,8) ).real();
    
    for(unsigned int i = 8; i < 10; i++)
    {
        // 13
        Phi_ll += (WC(i, int_qed(LO_QED)).conjugate()    * WC(i, int_qed(NLO_QED12)) * Hij[LO](i,i) +
                   WC(i, int_qed(NLO_QED11)).conjugate() * WC(i, int_qed(NLO_QED02)) * Hij[LO](i,i) +
                   WC(i, int_qed(NLO_QED02)).conjugate() * WC(i, int_qed(NLO_QED11)) * Hij[LO](i,i) +
                   WC(i, int_qed(NLO_QED12)).conjugate() * WC(i, int_qed(LO_QED))    * Hij[LO](i,i) +
                   WC(i, int_qed(LO_QED)).conjugate()    * WC(i, int_qed(NLO_QED02)) * Hij[NLO](i,i) +
                   WC(i, int_qed(NLO_QED02)).conjugate() * WC(i, int_qed(LO_QED))    * Hij[NLO](i,i) +
                   WC(i, int_qed(LO_QED)).conjugate()    * WC(i, int_qed(LO_QED))    * Hij[int_qed(NLO_QED11)](i,i) ).real();
        
        // 23, MISSING: Ci02 Cj01 Hij20 + Ci01 Cj02 Hij20
        Phi_ll += (WC(i, int_qed(LO_QED)).conjugate()    * WC(i, int_qed(NLO_QED22)) * Hij[LO](i,i) +
                   WC(i, int_qed(NLO_QED11)).conjugate() * WC(i, int_qed(NLO_QED12)) * Hij[LO](i,i) +
                   WC(i, int_qed(NLO_QED21)).conjugate() * WC(i, int_qed(NLO_QED02)) * Hij[LO](i,i) +
                   WC(i, int_qed(NLO_QED02)).conjugate() * WC(i, int_qed(NLO_QED21)) * Hij[LO](i,i) +
                   WC(i, int_qed(NLO_QED12)).conjugate() * WC(i, int_qed(NLO_QED11)) * Hij[LO](i,i) +
                   WC(i, int_qed(NLO_QED22)).conjugate() * WC(i, int_qed(LO_QED))    * Hij[LO](i,i) +
                   WC(i, int_qed(LO_QED)).conjugate()    * WC(i, int_qed(NLO_QED12)) * Hij[NLO](i,i) +
                   WC(i, int_qed(NLO_QED11)).conjugate() * WC(i, int_qed(NLO_QED02)) * Hij[NLO](i,i) +
                   WC(i, int_qed(NLO_QED02)).conjugate() * WC(i, int_qed(NLO_QED11)) * Hij[NLO](i,i) +
                   WC(i, int_qed(NLO_QED12)).conjugate() * WC(i, int_qed(LO_QED))    * Hij[NLO](i,i) +
                   WC(i, int_qed(LO_QED)).conjugate()    * WC(i, int_qed(NLO_QED11)) * Hij[int_qed(NLO_QED11)](i,i) +
                   WC(i, int_qed(NLO_QED11)).conjugate() * WC(i, int_qed(LO_QED))    * Hij[int_qed(NLO_QED11)](i,i) +
                   WC(i, int_qed(LO_QED)).conjugate()    * WC(i, int_qed(LO_QED))    * Hij[int_qed(NLO_QED21)](i,i) ).real();   
    }
    
    // 33, MISSING: Ci11 Cj02 Hij20 + Ci02 Cj11 Hij20
    Phi_ll += (WC(9, int_qed(NLO_QED11)).conjugate() * WC(9, int_qed(NLO_QED22)) * Hij[LO](9,9) +
               WC(9, int_qed(NLO_QED21)).conjugate() * WC(9, int_qed(NLO_QED12)) * Hij[LO](9,9) +
               WC(9, int_qed(NLO_QED12)).conjugate() * WC(9, int_qed(NLO_QED21)) * Hij[LO](9,9) +
               WC(9, int_qed(NLO_QED22)).conjugate() * WC(9, int_qed(NLO_QED11)) * Hij[LO](9,9) +
               WC(9, int_qed(NLO_QED11)).conjugate() * WC(9, int_qed(NLO_QED12)) * Hij[NLO](9,9) +
               WC(9, int_qed(NLO_QED21)).conjugate() * WC(9, int_qed(NLO_QED02)) * Hij[NLO](9,9) +
               WC(9, int_qed(NLO_QED02)).conjugate() * WC(9, int_qed(NLO_QED21)) * Hij[NLO](9,9) +
               WC(9, int_qed(NLO_QED12)).conjugate() * WC(9, int_qed(NLO_QED11)) * Hij[NLO](9,9) +
               WC(9, int_qed(NLO_QED11)).conjugate() * WC(9, int_qed(NLO_QED11)) * Hij[int_qed(NLO_QED11)](9,9) ).real();
    
    return Phi_ll;
}

double BXqll::FULLCCH_multiplication(std::vector< gslpp::matrix<gslpp::complex> >& Hij)
{
    double Phi_ll = 0.;
    gslpp::vector<gslpp::complex> FULLWC(15, 0.);
    gslpp::matrix<gslpp::complex> FULLHij(15, 15, 0.);
    
    for(unsigned int j = 0; j < 15; j++)
    {
        FULLWC.assign(j, WC(j, LO) + WC(j, NLO) + WC(j, NNLO) +
                WC(j, int_qed(LO_QED)) + WC(j, int_qed(NLO_QED11)) + WC(j, int_qed(NLO_QED21)) +
                WC(j, int_qed(NLO_QED02)) + WC(j, int_qed(NLO_QED12)) + WC(j, int_qed(NLO_QED22)) );
        for(unsigned int i = 0; i <= j; i++)
        {
            FULLHij.assign(i,j,
                    Hij[LO](i,j) + Hij[NLO](i,j) + Hij[NNLO](i,j) +
                    Hij[int_qed(LO_QED)](i,j) + Hij[int_qed(NLO_QED11)](i,j) + Hij[int_qed(NLO_QED21)](i,j) +
                    Hij[int_qed(NLO_QED02)](i,j) + Hij[int_qed(NLO_QED12)](i,j) + Hij[int_qed(NLO_QED22)](i,j) );
            
            Phi_ll += (FULLWC(i) * FULLWC (j).conjugate() * FULLHij(i,j) ).real();
        }
    }
    
    return (Phi_ll);
}

void BXqll::Test_WC_DF1()
{
//    double muw = mySM.getMuw();
//    double alstilde_2 = alstilde * alstilde;
//    double kappa_2 = kappa * kappa;
//    double GF_2 = GF * GF;
//    double MW_2 = MW * MW;
////    double sw2 = 0.2312;
//    double MZ = mySM.getMz();
//    double sw2 = 1. - (MW_2 / MZ / MZ);
//    gslpp::complex Vtb = mySM.getCKM().getV_tb();
//    gslpp::complex Vts = mySM.getCKM().getV_ts();
//    gslpp::complex lambda_t = Vtb.conjugate() * Vts;
//    double xt = mySM.getMatching().x_t(muw);
//    double alemuw = mySM.Ale(muw, FULLNLO);
//    double aletildemuw = alemuw / 4. / M_PI;
//    
//    std::cout << "MW         = " << MW << std::endl;
//    std::cout << "1/ale(muw) = " << 1./alemuw << std::endl;
//    std::cout << "1/ale(mub) = " << 1./ale << std::endl;
//    std::cout << "xt         = " << xt << std::endl; 
//     
//    gslpp::complex c7_df1 = 0.;
//    gslpp::complex c9_df1 = 0.;
//    gslpp::complex c10_df1 = 0.;
//    
//    for (int i = 0; i < 5; i++)
//    {
//        c7_df1 += WC(6, i);
//        c9_df1 += WC(8, i);
//        c10_df1 += WC(9, i);
//    }
//    for (int i = 5; i < 9; i++)
//    {
//        c9_df1 += WC(8, i);
//        c10_df1 += WC(9, i);
//    }
//    
//    std::cout << std::endl;
//    
//    std::cout << "00:      " << WC(6,0) << std::endl;
//    std::cout << "10:      " << WC(6,1) / alstilde << std::endl;
//    std::cout << "20:      " << WC(6,2) / alstilde / alstilde << std::endl;
//    std::cout << "01:      " << WC(6,3) / kappa << std::endl;
//    std::cout << "11:      " << WC(6,4) / alstilde / kappa << std::endl;
//    std::cout << "C7:      " << c7_df1 << std::endl << std::endl;
//    
//    std::cout << "01:      " << WC(8,3) / kappa << std::endl;
//    std::cout << "11:      " << WC(8,4) / alstilde / kappa << std::endl;
//    std::cout << "21:      " << WC(8,5) / alstilde_2 / kappa << std::endl;
//    std::cout << "02:      " << WC(8,6) / kappa_2 << std::endl;
//    std::cout << "12:      " << WC(8,7) / alstilde / kappa_2 << std::endl;
//    std::cout << "22:      " << WC(8,8) / alstilde_2 / kappa_2 << std::endl;
//    std::cout << "C9:      " << c9_df1 / aletilde << std::endl << std::endl;
//    
//    std::cout << "11:      " << WC(9,4) / alstilde / kappa << std::endl;
//    std::cout << "21:      " << WC(9,5) / alstilde_2 / kappa << std::endl;
//    std::cout << "02:      " << WC(9,6) / kappa_2 << std::endl;
//    std::cout << "12:      " << WC(9,7) / alstilde / kappa_2 << std::endl;
//    std::cout << "22:      " << WC(9,8) / alstilde_2 / kappa_2 << std::endl;
//    std::cout << "C10:     " << c10_df1 / aletilde << std::endl;
//    
//    std::cout << "C9 (27): " << (c9_df1 - WC(8, int_qed(NLO_QED22)) + alstilde_2*kappa_2*27.32) / aletilde << std::endl;
//    
//    
//    std::cout << "C10(-36): " << (c10_df1 - WC(9, int_qed(NLO_QED22)) +
//                                  alstilde_2 * kappa_2 * (-36.09) ) / aletilde << std::endl;
//
//    gslpp::complex c10_stu = ( (*(allcoeff_smm[NLO_QED11]))(7) +
//                            alstilde * (*(allcoeff_smm[NLO_QED21]))(7) +
//                            kappa / alstilde * (*(allcoeff_smm[NLO_QED02]))(7) +
//                            kappa * (*(allcoeff_smm[NLO_QED12]))(7) +
//                            alstilde * kappa * (*(allcoeff_smm[NLO_QED22]))(7) ) / lambda_t;
//    
//    std::cout << "C10_stu:  " << c10_stu / sw2 << std::endl;
//    std::cout << "Should be zero: " << (*(allcoeff_smm[LO]))(7) + (*(allcoeff_smm[NLO]))(7) +
//                                    (*(allcoeff_smm[NNLO]))(7) + (*(allcoeff_smm[LO_QED]))(7) << std::endl;
//    
//    std::cout << std::endl;
//    std::cout << "4 GF / sqrt(2) * C10   = " << 4. * GF * c10_df1 / sqrt(2.) << std::endl;
//    std::cout << "GF^2 MW^2 / pi^2 * C10 = " << GF_2 * MW_2 * c10_stu / M_PI / M_PI << std::endl;
//
//    std::cout << std::endl;
//    std::cout << "C10_OS1: " << mySM.getMatching().C10_OS1(mt_w * mt_w / MW_2, muw) << std::endl;
//    std::cout << "Rest: " << mySM.getMatching().Rest(mt_w * mt_w / MW_2, muw) << std::endl;
//
//    const std::vector<WilsonCoefficientNew>& match_df1 = mySM.getMatching().CMDF1("CPMLQB",15);
//    const std::vector<WilsonCoefficient>& match_stu = mySM.getMatching().CMbsmm();
//    gslpp::complex c10_11_df1 = match_df1[0].getCoeff(QCD1, QED1)(9);
//    gslpp::complex c10_22_df1 = match_df1[0].getCoeff(QCD2, QED2)(9);
//    gslpp::complex c10_11_stu = (*(match_stu[0].getCoeff(orders_qed(NLO_QED11))))(7) / lambda_t;
//    gslpp::complex c10_22_stu = (*(match_stu[0].getCoeff(orders_qed(NLO_QED22))))(7) / lambda_t;
//    
//    std::cout << std::endl;
//    std::cout << "c10_11_df1:     " << c10_11_df1/aletildemuw << std::endl;
//    std::cout << "c10_11_stu:     " << c10_11_stu << std::endl;
//    std::cout << "c10_11_stu/sW2: " << c10_11_stu/sw2 << std::endl;
//    std::cout << "c10_22_df1:     " << c10_22_df1/aletildemuw/aletildemuw << std::endl;
//    std::cout << "c10_22_stu:     " << c10_22_stu << std::endl;
//    std::cout << "c10_22_stu/sW2: " << c10_22_stu/sw2 << std::endl;
//    std::cout << "GF matching:    " << (4. * GF / sqrt(2.)) * (c10_11_df1 + c10_22_df1) << std::endl;
//    std::cout << "GF^2 matching:  " << (GF_2 * MW_2 / M_PI / M_PI) * (c10_11_stu +
//                                    (mySM.Ale(muw,FULLNLO) / 4. / M_PI) * c10_22_stu) << std::endl;
        
    //ATTEMPT AT A WORKING EXPANDED * EXPANDED
//    gslpp::vector<double> vec2(2,0.);
//    std::vector<gslpp::vector<double> > vtmp;
//    std::vector<std::vector<gslpp::vector<double> > > vtmp12;
//    std::vector<std::vector<gslpp::vector<double> > > vtmp910;
//
//    for (int j = 0; j <= QCD2; j++)
//        vtmp.push_back(vec2);
//    for (int i = 0; i <= QED1; i++)
//        vtmp12.push_back(vtmp);
//    for (int i = 0; i <= QED2; i++)
//        vtmp910.push_back(vtmp);
//    
//    Expanded<gslpp::vector<double> > wilson12(vtmp12);
//    Expanded<gslpp::vector<double> > wilson910(vtmp910);
//    
//    wilson12.setVectorElement(QCD0, QED0, 0, 1.);
//    wilson12.setVectorElement(QCD1, QED0, 0, 2.);
//    wilson12.setVectorElement(QCD2, QED0, 0, 3.);
//    wilson12.setVectorElement(QCD0, QED1, 0, 4.);
//    wilson12.setVectorElement(QCD0, QED0, 1, 1.);
//    wilson12.setVectorElement(QCD1, QED0, 1, 2.);
//    wilson12.setVectorElement(QCD2, QED0, 1, 3.);
//    wilson12.setVectorElement(QCD0, QED1, 1, 4.);
//    
//    wilson910.setVectorElement(QCD1, QED1, 0, 5.);
//    wilson910.setVectorElement(QCD2, QED2, 0, 6.);
//    wilson910.setVectorElement(QCD1, QED1, 1, 5.);
//    wilson910.setVectorElement(QCD2, QED2, 1, 6.);
//
//    
//    std::cout << "00 " << wilson12.getOrd(QCD0,QED0) << std::endl;
//    std::cout << "10 " << wilson12.getOrd(QCD1,QED0) << std::endl;
//    std::cout << "20 " << wilson12.getOrd(QCD2,QED0) << std::endl;
//    std::cout << "01 " << wilson12.getOrd(QCD0,QED1) << std::endl;
//    
//    std::cout << "11 " << wilson910.getOrd(QCD1,QED1) << std::endl;
//    std::cout << "22 " << wilson910.getOrd(QCD2,QED2) << std::endl;
//    
//    std::cout << std::endl;
//    std::cout << wilson12*wilson12 << std::endl;
//    std::cout << std::endl;
//    std::cout << wilson12*wilson910 << std::endl;
//    std::cout << std::endl;
//    std::cout << wilson910*wilson910 << std::endl;
    
    //PRINT OF AUXILIARY FUNCTIONS
    std::cout << "alstilde = " << alstilde << std::endl;
    std::cout << "aletilde = " << aletilde << std::endl;
    std::cout << "Mb_pole  = " << Mb_pole << std::endl;
    std::cout << "Mc_pole  = " << Mc_pole << std::endl;
    std::cout << std::endl;
    
    double s = 0.15;
    computeMi(s);
    
    std::cout << "H_T = " << getH("T", s) << std::endl;
    std::cout << "H_L = " << getH("L", s) << std::endl;
    
    std::cout << std::endl;    
    std::cout << "KS_cc = " << KS_cc(s) << std::endl;
    
    std::cout << std::endl;
    for (unsigned int i = 0; i < 10; i++)
    {
        std::cout << "M_" << i+1 << "^7 = " << M_7[LO](i) + M_7[NLO](i) + M_7[NNLO](i) + M_7[int_qed(LO_QED)](i) +
                M_7[int_qed(NLO_QED11)](i) + M_7[int_qed(NLO_QED21)](i) << std::endl;
        std::cout << "M_" << i+1 << "^9 = " << M_9[LO](i) + M_9[NLO](i) + M_9[NNLO](i) + M_9[int_qed(LO_QED)](i) +
                M_9[int_qed(NLO_QED11)](i) + M_9[int_qed(NLO_QED21)](i) << std::endl;
    }
    
    for (unsigned int i = 10; i < 14; i++)
    {
        std::cout << "M_" << i-7 << "Q^7 = " << M_7[LO](i) + M_7[NLO](i) + M_7[NNLO](i) + M_7[int_qed(LO_QED)](i) +
                M_7[int_qed(NLO_QED11)](i) + M_7[int_qed(NLO_QED21)](i) << std::endl;
        std::cout << "M_" << i-7 << "Q^9 = " << M_9[LO](i) + M_9[NLO](i) + M_9[NNLO](i) + M_9[int_qed(LO_QED)](i) +
                M_9[int_qed(NLO_QED11)](i) + M_9[int_qed(NLO_QED21)](i) << std::endl;
    }
    
    std::cout << "M_b^7 = " << M_7[LO](14) + M_7[NLO](14) + M_7[NNLO](14) + M_7[int_qed(LO_QED)](14) +
                M_7[int_qed(NLO_QED11)](14) + M_7[int_qed(NLO_QED21)](14) << std::endl;
    std::cout << "M_b^9 = " << M_9[LO](14) + M_9[NLO](14) + M_9[NNLO](14) + M_9[int_qed(LO_QED)](14) +
                M_9[int_qed(NLO_QED11)](14) + M_9[int_qed(NLO_QED21)](14) << std::endl;
    
//    std::cout << std::endl;
//    std::cout << "S_77^T   = " << S77_T(s, LO) + S77_T(s, NLO) << std::endl;
//    std::cout << "S_99^T   = " << S99_T(s, LO) + S99_T(s, NLO) + S99_T(s, NNLO) << std::endl;
//    std::cout << "S_79^T   = " << S79_T(s, LO) + S79_T(s, NLO) << std::endl;
//    std::cout << "S_1010^T = " << S1010_T(s, LO) + S1010_T(s, NLO) + S1010_T(s, NNLO) << std::endl;
//    
//    std::cout << std::endl;
//    std::cout << "S_77^L   = " << S77_L(s, LO) + S77_L(s, NLO) << std::endl;
//    std::cout << "S_99^L   = " << S99_L(s, LO) + S99_L(s, NLO) + S99_L(s, NNLO) << std::endl;
//    std::cout << "S_79^L   = " << S79_L(s, LO) + S79_L(s, NLO) << std::endl;
//    std::cout << "S_1010^L = " << S1010_L(s, LO) + S1010_L(s, NLO) + S1010_L(s, NNLO) << std::endl;
//    
//    std::cout << std::endl;
//    std::cout << "S_710^A   = " << S710_A(s, LO) + S710_A(s, NLO) << std::endl;
//    std::cout << "S_910^A   = " << S910_A(s, LO) + S910_A(s, NLO) + S910_A(s, NNLO) << std::endl;
//    
    std::cout << std::endl;
    std::cout << "f_1  = " << f_Huber(s,   -32./27., 4./3.,      0.,        0.,   -16./27.) << std::endl;
    std::cout << "f_2  = " << f_Huber(s,     -8./9.,    1.,      0.,        0.,     -4./9.) << std::endl;
    std::cout << "f_3  = " << f_Huber(s,    -16./9.,    6.,  -7./2.,     2./9.,     2./27.) << std::endl;
    std::cout << "f_4  = " << f_Huber(s,    32./27.,    0.,  -2./3.,    8./27.,     8./81.) << std::endl;
    std::cout << "f_5  = " << f_Huber(s,   -112./9.,   60.,    -38.,    32./9.,  -136./27.) << std::endl;
    std::cout << "f_6  = " << f_Huber(s,   512./27.,    0., -32./3.,  128./27.,   320./81.) << std::endl;
    std::cout << "f_9  = " << f9pen_Huber(s) << std::endl;
    std::cout << "f_3Q = " << f_Huber(s,  -272./27.,    4.,   7./6.,  -74./27.,   358./81.) << std::endl;
    std::cout << "f_4Q = " << f_Huber(s,   -32./81.,    0.,   2./9.,   -8./81.,   -8./243.) << std::endl;
    std::cout << "f_5Q = " << f_Huber(s, -2768./27.,   40.,  38./3., -752./27.,  1144./81.) << std::endl;
    std::cout << "f_6Q = " << f_Huber(s,  -512./81.,    0.,  32./9., -128./81., -320./243.) << std::endl;
    std::cout << "f_b  = " << f_Huber(s,     16./9.,    0.,     -2.,        0.,    26./27.) << std::endl;
//    
//    std::cout << std::endl;
//    std::cout << "F_1^7 = " << F17(s) << std::endl;
//    std::cout << "F_2^7 = " << F27(s) << std::endl;
//    std::cout << "F_8^7 = " << F87(s) << std::endl;
//    std::cout << "F_1^9 = " << F19(s) << std::endl;
//    std::cout << "F_2^9 = " << F29(s) << std::endl;
//    std::cout << "F_8^9 = " << F89(s) << std::endl;
//    
//    double umsh = 1. - s;
//    std::cout << std::endl;
//    std::cout << "sigma_77^T = " << 8.*umsh*umsh/s << std::endl;
//    std::cout << "sigma_99^T = " << 2.*s*umsh*umsh << std::endl;
//    std::cout << "sigma_79^T = " << 8.*umsh*umsh << std::endl;
//    std::cout << "sigma_77^L = " << 4.*umsh*umsh << std::endl;
//    std::cout << "sigma_99^L = " << umsh*umsh << std::endl;
//    std::cout << "sigma_79^L = " << 4.*umsh*umsh << std::endl;
//    
//    std::cout << std::endl;
//    std::cout << "w1_77^T = " << omega77_T(s) << std::endl;
//    std::cout << "w1_99^T = " << omega99_T(s) << std::endl;
//    std::cout << "w1_79^T = " << omega79_T(s) << std::endl;
//    std::cout << "w1_77^L = " << omega77_L(s) << std::endl;
//    std::cout << "w1_99^L = " << omega99_L(s) << std::endl;
//    std::cout << "w1_79^L = " << omega79_L(s) << std::endl;
//    
//    std::cout << std::endl;
//    std::cout << "chi1_77^T = " << 4.*umsh*(5.*s + 3.)/3./s << std::endl;
//    std::cout << "chi1_99^T = " << -s*umsh*(3.*s + 5.)/3. << std::endl;
//    std::cout << "chi1_79^T = " << 4.*umsh*umsh << std::endl;
//    std::cout << "chi1_77^L = " << -2.*umsh*(3.*s + 13.)/3. << std::endl;
//    std::cout << "chi1_99^L = " << umsh*(13.*s + 3.)/6. << std::endl;
//    std::cout << "chi1_79^L = " << 2.*umsh*umsh << std::endl;
//    
//    std::cout << std::endl;
//    std::cout << "chi2_77^T = " << 4.*(3.*s*s + 2.*s - 9.)/s << std::endl;
//    std::cout << "chi2_99^T = " << s*(15.*s*s - 14.*s - 5.) << std::endl;
//    std::cout << "chi2_79^T = " << 4.*(9.*s*s - 6.*s - 7.) << std::endl;
//    std::cout << "chi2_77^L = " << 2.*(15.*s*s - 6.*s - 13.) << std::endl;
//    std::cout << "chi2_99^L = " << (-17.*s*s + 10.*s + 3.)/2. << std::endl;
//    std::cout << "chi2_79^L = " << 2.*(3.*s*s - 6.*s - 1.) << std::endl;
//    
//    std::cout << std::endl;
//    for (unsigned int i = 0; i < 2; i++)
//    {
//        for (unsigned int j = 0; j < 15; j++)
//            std::cout << "c^T_" << i+1 << "," << j+1 << " = " << cij_T(i, j, s) << std::endl;
//        std::cout << std::endl;
//    }
//    
//    for (unsigned int i = 0; i < 2; i++)
//    {
//        for (unsigned int j = 0; j < 15; j++)
//            std::cout << "c^L_" << i+1 << "," << j+1 << " = " << cij_L(i, j, s) << std::endl;
//        std::cout << std::endl;
//    }
//    
//    std::cout << "c^A_1,10 = " << cij_A(0, 9, s) << std::endl;
//    std::cout << "c^A_2,10 = " << cij_A(1, 9, s) << std::endl;
//    
//    std::cout << std::endl;  
//    std::cout << "e^T_1,1 = " << eij_T(0, 0, s) << std::endl;
//    std::cout << "e^T_1,2 = " << eij_T(0, 1, s) << std::endl;
//    std::cout << "e^T_1,7 = " << eij_T(0, 6, s) << std::endl;
//    std::cout << "e^T_1,9 = " << eij_T(0, 8, s) << std::endl;
//    std::cout << "e^T_2,2 = " << eij_T(1, 1, s) << std::endl;
//    std::cout << "e^T_2,7 = " << eij_T(1, 6, s) << std::endl;
//    std::cout << "e^T_2,9 = " << eij_T(1, 8, s) << std::endl;
//    std::cout << "e^T_7,7 = " << eij_T(6, 6, s) << std::endl;
//    std::cout << "e^T_7,9 = " << eij_T(6, 8, s) << std::endl;
//    std::cout << "e^T_9,9 = " << eij_T(8, 8, s) << std::endl;
//    std::cout << "e^T_10,10 = " << eij_T(9, 9, s) << std::endl;
//    
//    std::cout << std::endl;
//    std::cout << "e^L_1,1 = " << eij_L(0, 0, s) << std::endl;
//    std::cout << "e^L_1,2 = " << eij_L(0, 1, s) << std::endl;
//    std::cout << "e^L_1,7 = " << eij_L(0, 6, s) << std::endl;
//    std::cout << "e^L_1,9 = " << eij_L(0, 8, s) << std::endl;
//    std::cout << "e^L_2,2 = " << eij_L(1, 1, s) << std::endl;
//    std::cout << "e^L_2,7 = " << eij_L(1, 6, s) << std::endl;
//    std::cout << "e^L_2,9 = " << eij_L(1, 8, s) << std::endl;
//    std::cout << "e^L_7,7 = " << eij_L(6, 6, s) << std::endl;
//    std::cout << "e^L_7,9 = " << eij_L(6, 8, s) << std::endl;
//    std::cout << "e^L_9,9 = " << eij_L(8, 8, s) << std::endl;
//    std::cout << "e^L_10,10 = " << eij_L(9, 9, s) << std::endl;
//    
//    std::cout << std::endl;
//    std::cout << "e^A_1,10 = " << eij_A(0, 9, s) << std::endl;
//    std::cout << "e^A_2,10 = " << eij_A(1, 9, s) << std::endl;
//    std::cout << "e^A_7,10 = " << eij_A(6, 9, s) << std::endl;
//    std::cout << "e^A_9,10 = " << eij_A(8, 9, s) << std::endl;
}


/************************************************************************
 * Finite bremsstrahlung corrections in the notation of Asatryan:2002iy *
 ************************************************************************/

double BXqll::R_bremsstrahlung(double sh, q2regions q2region)
{
    gslpp::complex c78, c89, c17, c27, c18, c28, c19, c29;
    double c88, c11, c12, c22;
    double Brem_a;
    double Brem_b;
    double ctau1 = (3. * 3. - 1.) / 8. / 3. / 3. / 3.;
    double ctau2 = - (3. * 3. - 1.) / 4. / 3. / 3.;
   
    c78 = CF * WC(6, LO) * WC(7, LO).conjugate();
    c89 = CF * WC(7, LO) * C9mod(sh).conjugate();
    c88 = CF * WC(7,LO).abs2();
    
    c11 = ctau1 * WC(0, LO).abs2();
    c12 = ctau2 * 2. * (WC(0, LO) * WC(1, LO).conjugate()).real();
    c22 = CF * WC(1, LO).abs2();
    
    c17 = ctau2 * WC(0, LO) * WC(6, LO).conjugate();
    c18 = ctau2 * WC(0, LO) * WC(7, LO).conjugate();
    c19 = ctau2 * WC(0, LO) * C9mod(sh).conjugate();
    c27 = CF * WC(1, LO) * WC(6, LO).conjugate();
    c28 = CF * WC(1, LO) * WC(7, LO).conjugate();
    c29 = CF * WC(1, LO) * C9mod(sh).conjugate();
    
    Brem_a = 2.*(c78*tau78(sh)+c89*tau89(sh)).real() + c88*tau88(sh);
    
    if(q2region <= HIGHQ2)
        Brem_b = (c11 + c12 + c22)*tau22fit(sh,q2region)+2.*(c17 + c27).real()*tau27fit_Re(sh,q2region)-
            2.*(c17 + c27).imag()*tau27fit_Im(sh,q2region)+2.*(c18 + c28).real()*tau28fit_Re(sh,q2region)-
            2.*(c18 + c28).imag()*tau28fit_Im(sh,q2region)+2.*(c19 + c29).real()*tau29fit_Re(sh,q2region)-
            2.*(c19 + c29).imag()*tau29fit_Im(sh,q2region);
    else
        throw std::runtime_error("BXqll::R_bremsstrahlung: q2 region not implemented");
    
    return (Brem_a + Brem_b);
}

gslpp::complex BXqll::C9mod(double sh)
{    
    gslpp::complex T9 = 4./3. * WC(0, LO) + WC(1, LO) + 6. * WC(2, LO) + 60. * WC(4, LO);
    gslpp::complex U9 = -7./2. * WC(2, LO) - 2./3. * WC(3, LO) - 38. * WC(4, LO) - 32./3. * WC(5, LO);
    gslpp::complex W9 = -1./2. * WC(2, LO) - 2./3. * WC(3, LO) -  8. * WC(4, LO) - 32./3. * WC(5, LO);
    gslpp::complex A9 = -32./27. * WC(0, LO) -   8./9. * WC(1, LO) -   16./9. * WC(2, LO) +
                         32./27. * WC(3, LO) - 112./9. * WC(4, LO) + 512./27. * WC(5, LO);
    A9 *= log(1. / muh);
    A9 += (WC(8, int_qed(LO_QED)) + WC(8, int_qed(NLO_QED11))) / aletilde;
    A9 += 4./3. * WC(2, LO) + 64./9. * WC(4, LO) + 64./27. * WC(5, LO);
    
    return (A9 + T9 * h_z(z, sh) + U9 * h_z(1., sh) + W9 * h_z(0., sh));
}

gslpp::complex BXqll::h_z(double zed, double sh)
{
    gslpp::complex i = gslpp::complex::i();
    gslpp::complex h_z;
    
    if (zed == 0.)
    {
        h_z = 8./27.-4./9.*(log(sh)-i*M_PI);
    }
    else
    {
        h_z = -2./9.*(2.+4./sh*zed)*sqrt(std::abs(4.*zed-sh)/sh);
        if(sh > 4.*zed)
            h_z *= log((sqrt(sh)+sqrt(sh-4.*zed))/(sqrt(sh)-sqrt(sh-4.*zed)))-i*M_PI;
        else
            h_z *= 2.*atan(sqrt(sh/(4.*zed-sh)));
        
        h_z += -4./9.*log(zed)+8./27.+16./9.*zed/sh;
   }   
     
    return(h_z);
}

double BXqll::tau78(double sh)
{
    gslpp::complex i = gslpp::complex::i();
    double dmsh = 2.-sh;
    double qmsh = 4.-sh;
    
    return (8./9./sh*(25.-2.*M_PI*M_PI-27.*sh+3.*sh*sh-sh*sh*sh+12.*(sh+sh*sh)*log(sh)+
            6.*(M_PI/2.-atan((2.-4.*sh+sh*sh)/dmsh/sqrt(sh)/sqrt(qmsh)))*
            (M_PI/2.-atan((2.-4.*sh+sh*sh)/dmsh/sqrt(sh)/sqrt(qmsh)))-
            24.*(dilog((sh-i*sqrt(sh)*sqrt(qmsh))/2.)).real()-12.*((1.-sh)*sqrt(sh)*sqrt(qmsh)-
            2.*atan(sqrt(sh)*sqrt(qmsh)/dmsh))*(atan(sqrt(qmsh/sh))-atan(sqrt(sh)*sqrt(qmsh)/dmsh))));
}

double BXqll::tau89(double sh)
{
    gslpp::complex i = gslpp::complex::i();
    double dmsh = 2.-sh;
    double qmsh = 4.-sh;
    
    return (2./3.*(sh*qmsh-3.-4.*log(sh)*(1.-sh-sh*sh)-8.*(dilog(sh/2.+i*sqrt(sh)*sqrt(qmsh)/2.)-
            dilog((sh*qmsh-2.)/2.+i*dmsh*sqrt(sh)*sqrt(qmsh)/2.)).real()+
            4.*(sh*sh*sqrt(qmsh/sh)+2.*atan(sqrt(sh)*sqrt(qmsh)/dmsh))*(atan(sqrt(qmsh/sh))-
            atan(sqrt(sh)*sqrt(qmsh)/dmsh))));
}

double BXqll::tau88(double sh)
{   
    gslpp::complex i = gslpp::complex::i();
    double umsh = 1.-sh;
    double qmsh = 4.-sh;
    
    return (4./27./sh*(-8.*M_PI*M_PI+umsh*(77.-sh-4.*sh*sh)-24.*dilog((gslpp::complex) umsh).real()+
            3.*(10.-4.*sh-9.*sh*sh+8.*log(sqrt(sh)/umsh))*log(sh)+48.*(dilog((3.-sh)/2.+
            i*umsh*sqrt(qmsh)/2./sqrt(sh))).real()-6.*((20.*sh+10.*sh*sh-3.*sh*sh*sh)/sqrt(sh)/sqrt(qmsh)
            -8.*M_PI+8.*atan(sqrt(qmsh/sh)))*(atan(sqrt(qmsh/sh))-atan(sqrt(sh)*sqrt(qmsh)/(2.-sh)))));
}

double BXqll::tau22fit(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = -186.96738127 + 1313.45792139*sh - 8975.40399683*sh*sh + 47018.56440838*sh*sh*sh -
                    159603.3217871*sh*sh*sh*sh + 309228.13379963*sh*sh*sh*sh*sh - 258317.14719949*sh*sh*sh*sh*sh*sh -
                    51.2467544*log(sh);
            break;
        case HIGHQ2:
            fit = -322.73989723 + 4.75813656/sh/sh - 80.36414222/sh + 687.70415138*sh - 491.08241967*sh*sh +
                    303.28125994*sh*sh*sh - 132.82124268*sh*sh*sh*sh + 35.63127394*sh*sh*sh*sh*sh -
                    4.36712003*sh*sh*sh*sh*sh*sh - 306.899641*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau22fit: region of q^2 not implemented");
    }
    
    return (fit);
}

double BXqll::tau27fit_Re(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = -45.40905903+334.92509492*sh-2404.69181358*sh*sh+12847.93973401*sh*sh*sh-
                    44421.35127703*sh*sh*sh*sh+87786.54536182*sh*sh*sh*sh*sh-75574.96266083*sh*sh*sh*sh*sh*sh-
                    13.79251644*log(sh);
            break;
        case HIGHQ2:
            fit = 87.43391175-196.67646862*sh+219.51106756*sh*sh-184.44868587*sh*sh*sh+
                    103.59892754*sh*sh*sh*sh-34.56056777*sh*sh*sh*sh*sh+5.14181565*sh*sh*sh*sh*sh*sh+
                    38.55667004*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau27fit_Re: region of q^2 not implemented");
    }
    
    return (fit);
}

double BXqll::tau27fit_Im(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = -189.61083508+1349.85607262*sh-9198.62227938*sh*sh+48104.40980548*sh*sh*sh-
                    162998.75872037*sh*sh*sh*sh+315224.375522*sh*sh*sh*sh*sh-262649.64320483*sh*sh*sh*sh*sh*sh-
                    52.52183304*log(sh);
            break;
        case HIGHQ2:
            fit = 523.76263422+49.97156504/sh-1120.42920341*sh+1024.46949612*sh*sh-767.28958612*sh*sh*sh+
                    393.62561539*sh*sh*sh*sh-120.74162898*sh*sh*sh*sh*sh+16.63110789*sh*sh*sh*sh*sh*sh+
                    352.74960196*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau27fit_Im: region of q^2 not implemented");
    }
    
    return (fit);
}

double BXqll::tau28fit_Re(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = 8.67757227-85.91172547*sh+666.57779497*sh*sh-3619.65234448*sh*sh*sh+
                    12475.74303361*sh*sh*sh*sh-24365.45545631*sh*sh*sh*sh*sh+20446.33269814*sh*sh*sh*sh*sh*sh+
                    1.54278226*log(sh);
            break;
        case HIGHQ2:
            fit = -4.11234905-0.52681762/sh+8.21844628*sh-6.04601107*sh*sh+3.67099354*sh*sh*sh-
                    1.57120958*sh*sh*sh*sh+0.41975346*sh*sh*sh*sh*sh-0.05280596*sh*sh*sh*sh*sh*sh-
                    3.16331567*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau28fit_Re: region of q^2 not implemented");
    }
    
    return (fit);
}

double BXqll::tau28fit_Im(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = 57.88258299-430.77957254*sh+3002.9999511*sh*sh-15808.63980887*sh*sh*sh+
                    53787.08410769*sh*sh*sh*sh-104360.60205475*sh*sh*sh*sh*sh+87294.84251167*sh*sh*sh*sh*sh*sh+
                    14.61062129*log(sh);
            break;
        case HIGHQ2:
            fit = -24.92802842+0.3842418/sh/sh-6.38294139/sh+53.15600599*sh-37.59024636*sh*sh+
                    23.04316804*sh*sh*sh-10.03556518*sh*sh*sh*sh+2.68088049*sh*sh*sh*sh*sh-
                    0.32751495*sh*sh*sh*sh*sh*sh-24.01652729*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau28fit_Im: region of q^2 not implemented");
    }
    
    return (fit);
}

double BXqll::tau29fit_Re(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = 0.53834924+0.47775224*sh-16.20063604*sh*sh+101.36668267*sh*sh*sh-
                    466.94537092*sh*sh*sh*sh+1224.77742613*sh*sh*sh*sh*sh-1469.41817323*sh*sh*sh*sh*sh*sh-
                    0.01686348*log(sh);
            break; 
        case HIGHQ2:
            fit = 4.46985355-6.16130742*sh+0.84917331*sh*sh+1.7696124*sh*sh*sh-1.14453916*sh*sh*sh*sh+
                    0.24261178*sh*sh*sh*sh*sh-0.02540446*sh*sh*sh*sh*sh*sh+2.67164817*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau29fit_Re: region of q^2 not implemented");
    }
    
    return (fit);
}

double BXqll::tau29fit_Im(double sh, q2regions q2region)
{
    double fit;
    
    switch(q2region)
    {
        case LOWQ2:
            fit = 0.7688748-0.21680402*sh-1.16934757*sh*sh+8.31833871*sh*sh*sh-4.81289468*sh*sh*sh*sh-
                    51.53765482*sh*sh*sh*sh*sh+158.06040784*sh*sh*sh*sh*sh*sh-0.00485643*log(sh);
            break;
        case HIGHQ2:
            fit = -38.80905455+95.60697233*sh-124.04368889*sh*sh+118.64599185*sh*sh*sh-
                    73.76081228*sh*sh*sh*sh+26.55080999*sh*sh*sh*sh*sh-4.19021877*sh*sh*sh*sh*sh*sh-
                    16.02711369*log(sh);
            break;
        default:
            throw std::runtime_error("BXqll::tau29fit_Im: region of q^2 not implemented");
    }
    
    return (fit);
}

/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cstring>
#include "EvolDF2.h"

EvolDF2::EvolDF2(unsigned int dim_i, schemes scheme, orders order, const StandardModel& model) 
:   RGEvolutor(dim_i, scheme, order),
    model(model),
    dim(dim_i)
{
    //double Nc = model.getNc();
    int basis = 0; //0: Gabbiani, 1: Buras
    gslpp::matrix<double> g0t(AnomalousDimension(LO, 3, basis).transpose());

    gslpp::matrix<gslpp::complex>vv(dim, dim, 0.);
    gslpp::vector<gslpp::complex>ee(dim, 0.);

    g0t.eigensystem(vv, ee);

    gslpp::matrix<double> v(vv.real());
    gslpp::vector<double> e(ee.real());

    gslpp::matrix<double> vi = v.inverse();
    for (unsigned int k = 0; k < dim; k++) {
        a[k] = e(k);
//        std::cout << "a[" << k << "] = " << a[k] << std::endl;
        for (unsigned int j = 0; j < dim; j++) {
            for (unsigned int i = 0; i < dim; i++) {
                b[i][j][k] = v(i, k) * vi(k, j);
//                if (b[i][j][k] != 0.)
//                    std::cout << "b[" << i << "][" << j << "][" << k << "] = " << b[i][j][k] << std::endl;
            }
        }
    }

    gslpp::matrix<double> h(dim, dim, 0.);
    for (int l = 0; l < 3; l++) {
        gslpp::matrix<double> gg = vi * (AnomalousDimension(NLO, 6 - l, basis).transpose()) * v;
        double b0 = model.Beta0(6 - l);
        for (unsigned int i = 0; i < dim; i++)
            for (unsigned int j = 0; j < dim; j++)
                h(i, j) = (i == j) * e(i) * model.Beta1(6 - l) / (2. * b0 * b0) -
                gg(i, j) / (2. * b0 + e(i) - e(j));
        gslpp::matrix<double> j = v * h * vi;
        gslpp::matrix<double> jv = j*v;
        gslpp::matrix<double> vij = vi*j;
        for (unsigned int i = 0; i < dim; i++)
            for (unsigned int j = 0; j < dim; j++)
                for (unsigned int k = 0; k < dim; k++) {
                    c[l][i][j][k] = jv(i, k) * vi(k, j);
//                    if (c[l][i][j][k] != 0.)
//                        std::cout << "c[" << l << "][" << i << "][" << j << "][" << k << "] = " << c[l][i][j][k] << std::endl;
                    d[l][i][j][k] = -v(i, k) * vij(k, j);
//                    if (d[l][i][j][k] != 0.)
//                        std::cout << "d[" << l << "][" << i << "][" << j << "][" << k << "] = " << d[l][i][j][k] << std::endl;
                }
    }
}

EvolDF2::~EvolDF2()
{}

gslpp::matrix<double> EvolDF2::AnomalousDimension(orders order, unsigned int nf, int basis) const
{
    gslpp::matrix<double> ad(dim, dim, 0.);
    double Nc = model.getNc();
    switch (basis) {
        case 0:
            switch (order) {
                case LO:
                    ad(0, 0) = (6. * Nc - 6.) / Nc;
                    ad(1, 1) = -(6. * Nc * Nc - 8. * Nc - 2.) / Nc;
                    ad(1, 2) = (4. * Nc - 8.) / Nc;
                    ad(2, 1) = (4. * Nc * Nc - 4. * Nc - 8.) / Nc;
                    ad(2, 2) = (2. * Nc * Nc + 4. * Nc + 2.) / Nc;
                    ad(3, 3) = -(6. * Nc * Nc - 6.) / Nc;
                    ad(4, 3) = -6.;
                    ad(4, 4) = 6. / Nc;
                    break;
                case NLO:

                    // MSbar-NDR scheme with evanescent operators of Buras, Misiak & Urban

                    if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6))
                        throw std::runtime_error("EvolDF2::AnomalousDimension(): wrong number of flavours");
                    ad(0, 0) = -(-1. + Nc)*(-171. + 19. * Nc * Nc + Nc * (63. - 4. * nf)) / (6. * Nc * Nc);
                    ad(1, 1) = (-1251. - 609. * Nc * Nc * Nc * Nc + Nc * (432. - 52. * nf) - 8. * Nc * Nc * (-71 + 2. * nf) +
                            20. * Nc * Nc * Nc * (32. + 3. * nf)) / (18. * Nc * Nc);
                    ad(1, 2) = -(2. * (-2. + Nc)*(-72. + Nc * Nc + 2. * Nc * (63. + nf))) / (9. * Nc * Nc);
                    ad(2, 1) = 2. * (119. * Nc * Nc * Nc + 8. * (-9. + 5. * nf) + 2. * Nc * (-125. + 7. * nf) -
                            Nc * Nc * (101. + 14. * nf)) / (9. * Nc);
                    ad(2, 2) = (477. + 343. * Nc * Nc * Nc * Nc + Nc * Nc * Nc * (380. - 52. * nf) - 4. * Nc * (-36. + nf) -
                            8. * Nc * Nc * (16. + 13. * nf)) / (18. * Nc * Nc);
                    ad(3, 3) = (45. + 479. * Nc * Nc - 203. * Nc * Nc * Nc * Nc - 44. * Nc * nf + 20. * Nc * Nc * Nc * nf) /
                            (6. * Nc * Nc);
                    ad(3, 4) = -(18. / Nc)-(71. * Nc) / 2. + 4. * nf;
                    ad(4, 3) = 3. / Nc - (100. * Nc) / 3. + (22. * nf) / 3.;
                    ad(4, 4) = (45. + 137. * Nc * Nc - 44. * Nc * nf) / (6. * Nc * Nc);
                    break;
                default:
                    throw std::runtime_error("EvolDF2::AnomalousDimension(): order not implemented");
            }
            break;
        case 1:
            switch (order) {
                case LO:
                    ad(0, 0) = 6. - 6. / Nc;
                    ad(1, 1) = 6. / Nc;
                    ad(1, 2) = 12.;
                    ad(2, 2) = 6. / Nc - 6. * Nc;
                    ad(3, 3) = 6. + 6. / Nc - 6. * Nc;
                    ad(3, 4) = 1. / 2. - 1. / Nc;
                    ad(4, 3) = -24. - 48. / Nc;
                    ad(4, 4) = 6. - 2. / Nc + 2. * Nc;
                    break;
                case NLO:

                    // MSbar-NDR scheme with evanescent operators of Buras, Misiak & Urban

                    if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6))
                        throw std::runtime_error("EvolDF2::AnomalousDimension(): wrong number of flavours");
                    ad(0, 0) = -22. / 3. - 57. / (2. * Nc * Nc) + 39. / Nc - (19. * Nc) / 6. + (2. * nf) / 3. - (2. * nf) / 3. / Nc;
                    ad(1, 1) = 137. / 6. + 15. / (2 * Nc * Nc) - (22 * nf) / (3 * Nc);
                    ad(1, 2) = -(6. / Nc) + (200. * Nc) / 3. - (44 * nf) / 3.;
                    ad(2, 1) = 9. / Nc + (71. * Nc) / 4. - 2. * nf;
                    ad(2, 2) = 479. / 6. + 15. / (2. * Nc * Nc) - (203. * Nc * Nc) / 6. - (22. * nf) / (3. * Nc) + (10. * Nc * nf) / 3.;
                    ad(3, 3) = 136. / 3. - 107. / (2. * Nc * Nc) - 12. / Nc + (107. * Nc) / 3. - (203. * Nc * Nc) / 6. - (2. * nf) / 3. - (10. * nf) / (3. * Nc) + (10. * Nc * nf) / 3.;
                    ad(3, 4) = -31. / 9. - 4. / (Nc * Nc) + 9. / Nc - Nc / 36. - nf / 18. + nf / (9 * Nc);
                    ad(4, 3) = -704. / 3. - 320. / (Nc * Nc) - 208. / Nc - (364. * Nc) / 3. + (136. * nf) / 3. + (176. * nf) / (3. * Nc);
                    ad(4, 4) = -188. / 9. + 21. / (2. * Nc * Nc) + 44. / Nc + 21. * Nc + (343. * Nc * Nc) / 18. - 6. * nf + (2. * nf) / (9. * Nc) - (26. * Nc * nf) / 9.;
                    break;
                default:
                    throw std::runtime_error("EvolDF2::AnomalousDimension(): order not implemented");
            }
            break;
        default:
            throw std::runtime_error("EvolDF2::AnomalousDimension(): basis not implemented");
    }
    return (ad);
}

gslpp::matrix<double>& EvolDF2::Df2Evol(double mu, double M, orders order, schemes scheme)
{
    switch (scheme) {
        case NDR:
            break;
        case LRI:
        case HV:
        default:
            std::stringstream out;
            out << scheme;
            throw std::runtime_error("EvolDF2::Df2Evol(): scheme " + out.str() + " not implemented ");
    }

    double alsMZ = model.getAlsMz();
    double Mz = model.getMz();
    if(alsMZ == alsMZ_cache && Mz == Mz_cache) {
        if (mu == this->mu && M == this->M && scheme == this->scheme)
            return (*Evol(order));        
    }
    alsMZ_cache = alsMZ;
    Mz_cache = Mz;

    if (M < mu) {
        std::stringstream out;
        out << "M = " << M << " < mu = " << mu;
        throw out.str();
    }

    setScales(mu, M); // also assign evol to identity

    double m_down = mu;
    double m_up = model.AboveTh(m_down);
    double nf = model.Nf(m_down);

    while (m_up < M) {
        Df2Evol(m_down, m_up, nf, scheme);
        m_down = m_up;
        m_up = model.AboveTh(m_down);
        nf += 1.;
    }
    Df2Evol(m_down, M, nf, scheme);
    return (*Evol(order));
}

void EvolDF2::Df2Evol(double mu, double M, double nf, schemes scheme)
{

    gslpp::matrix<double> resLO(dim, 0.), resNLO(dim, 0.), resNNLO(dim, 0.);

    int l = 6 - (int) nf;
    double alsM = model.Als(M, FULLNLO) / 4. / M_PI;
    double alsmu = model.Als(mu, FULLNLO) / 4. / M_PI;

    double eta = alsM / alsmu;

    for (unsigned int k = 0; k < dim; k++) {
        double etap = pow(eta, a[k] / 2. / model.Beta0(nf));
        for (unsigned int i = 0; i < dim; i++)
            for (unsigned int j = 0; j < dim; j++) {
                resNNLO(i, j) += 0.;
                resNLO(i, j) += c[l][i][j][k] * etap * alsmu;
                resNLO(i, j) += d[l][i][j][k] * etap * alsM;
                resLO(i, j) += b[i][j][k] * etap;
            }
    }
    switch (order) {
        case NNLO:
            *elem[NNLO] = 0.; // Marco can implement it if he wishes to!
        case NLO:
            *elem[NLO] = (*elem[LO]) * resNLO + (*elem[NLO]) * resLO;
        case LO:
            *elem[LO] = (*elem[LO]) * resLO;
            break;
        case FULLNNLO:
        case FULLNLO:
        default:
            throw std::runtime_error("Error in EvolDF2::Df2Evol()");
    }
}

double EvolDF2::etacc(double mu) const
{
    double N = model.getNc();
    double N2 = N * N;
    double gamma0[2] = {0.}, gamma1[2][4] = {
        {0.}
    }, d[2][4] = {
        {0.}
    },
    J[2][4] = {
        {0.}
    }, dd[2][2][4] = {
        {
            {0.}
        }
    }, JJ[2][2][4] = {
        {
            {0.}
        }
    },
    B[2] = {0.}; // 0 = + ; 1 = - ;
    double tau[2][2] = {
        {0.}
    }, K[2][2] = {
        {0.}
    }, beta[2][2] = {
        {0.}
    };

    gamma0[0] = 6. * (N - 1.) / N;
    gamma0[1] = -6. * (N + 1.) / N;

    gamma1[0][0] = ((N - 1.) / (2. * N)) * (-21. + 57. / N - 19. / 3. * N + 4. / 3. * 3.);
    gamma1[1][0] = ((N + 1.) / (2. * N)) * (-21. - 57. / N + 19. / 3. * N - 4. / 3. * 3.);
    gamma1[0][1] = ((N - 1.) / (2. * N)) * (-21. + 57. / N - 19. / 3. * N + 4. / 3. * 4.);
    gamma1[1][1] = ((N + 1.) / (2. * N)) * (-21. - 57. / N + 19. / 3. * N - 4. / 3. * 4.);
    gamma1[0][2] = ((N - 1.) / (2. * N)) * (-21. + 57. / N - 19. / 3. * N + 4. / 3. * 5.);
    gamma1[1][2] = ((N + 1.) / (2. * N)) * (-21. - 57. / N + 19. / 3. * N - 4. / 3. * 5.);
    gamma1[0][3] = ((N - 1.) / (2. * N)) * (-21. + 57. / N - 19. / 3. * N + 4. / 3. * 6.);
    gamma1[1][3] = ((N + 1.) / (2. * N)) * (-21. - 57. / N + 19. / 3. * N - 4. / 3. * 6.);

    for (int i = 0; i < 2; i++) { // 0 = + ; 1 = - ;
        for (int j = 0; j < 4; j++) { // j = nf - 3;
            d[i][j] = gamma0[i] / 2. / model.Beta0(j + 3);
            J[i][j] = d[i][j] / model.Beta0(j + 3) * model.Beta1(j + 3) -
                    gamma1[i][j] / 2. / model.Beta0(j + 3);
        }

    }

    for (int i = 0; i < 2; i++) { // 0 = + ; 1 = - ;
        for (int j = 0; j < 2; j++) { // 0 = + ; 1 = - ;
            for (int k = 0; k < 4; k++) { // k = nf - 3;
                dd[i][j][k] = d[i][k] + d[j][k];
                JJ[i][j][k] = J[i][k] + J[j][k];
            }
        }
    }

    tau[0][0] = (N + 3.) / 4.;
    tau[1][0] = -(N - 1.) / 4.;
    tau[0][1] = tau[1][0];
    tau[1][1] = (N - 1.) / 4.;

    K[0][0] = 3. * (N - 1.) * tau[0][0];
    K[1][0] = 3. * (N + 1.) * tau[0][1];
    K[0][1] = K[1][0];
    K[1][1] = 3. * (N + 3.) * tau[1][1];

    beta[0][0] = (1. - N) * (M_PI * M_PI * (N2 - 6.) / (12. * N) + 3. * (-N2 + 2. * N + 13.) / (4. * N));
    beta[1][0] = (1. - N) * (M_PI * M_PI * (-N2 + 2. * N - 2.) / (12. * N) + (3. * N2 + 13.) / (4. * N));
    beta[0][1] = beta[1][0];
    beta[1][1] = (1. - N) * (M_PI * M_PI * (N2 - 4. * N + 2) / (12. * N) - (3. * N2 + 10. * N + 13.) / (4. * N));

    B[0] = 11. * (N - 1.) / (2. * N);
    B[1] = -11. * (N + 1.) / (2. * N);

    double eta = 0.;
    double as_mb__as_muc = model.Als(model.getMub(), FULLNLO) / model.Als(model.getMuc(), FULLNLO);
    double as_muw__as_mb = model.Als(model.getMuw(), FULLNLO) / model.Als(model.getMub(), FULLNLO);
    double as_muc__4pi = model.Als(model.getMuc(), FULLNLO) / 4. / M_PI;
    double log_muc__mc = log(model.getMuc() / model.getQuarks(QCD::CHARM).getMass());
    double as_mb__4pi = model.Als(model.getMub(), FULLNLO) / 4. / M_PI;
    double as_muw__4pi = model.Als(model.getMuw(), FULLNLO) / 4. / M_PI;


    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            eta += pow(as_mb__as_muc, dd[i][j][1]) *
                    pow(as_muw__as_mb, dd[i][j][2]) *
                    (tau[i][j] + as_muc__4pi * (2. * K[i][j] *
                    log_muc__mc + tau[i][j] * (JJ[i][j][1] - J[0][0]) + beta[i][j]) +
                    tau[i][j]*(as_mb__4pi * (JJ[i][j][2] - JJ[i][j][1]) +
                    as_muw__4pi * (B[i] + B[j] - JJ[i][j][2])));
        }
    }

    eta *= pow(model.Als(model.getMuc()), d[0][0]);
    eta *= pow(model.Als(mu, FULLNLO), -2. / 9.);

    return (eta * (1. + model.Als(mu, FULLNLO) / 4. / M_PI * J[0][0]));
}

double EvolDF2::etact(double mu) const
{

    //temporary fix waiting for NNLO

    double K = model.Als4(model.getMuw()) / model.Als4(model.getMuc());
    double Kpp = pow(K, 12. / 25.);
    double Kpm = pow(K, -6. / 25.);
    double Kmm = pow(K, -24. / 25.);
    double K7 = pow(K, 1. / 5.);
    double xc = model.Mrun4(model.getMuc(), model.getQuarks(QCD::CHARM).getMass_scale(), model.getQuarks(QCD::CHARM).getMass()) / model.Mw();
    xc *= xc;
    double xt = model.Mrun4(model.getMuw(), model.getQuarks(QCD::TOP).getMass_scale(), model.getQuarks(QCD::TOP).getMass()) / model.Mw();
    xt *= xt;

    double J3 = 6. * (3. - 1.) / 3. / 2. / model.Beta0(3) / model.Beta0(3) * model.Beta1(3) -
            ((3. - 1.) / (2. * 3.)) * (-21. + 57. / 3. - 19. + 4.) / 2. / model.Beta0(3);


    double eta = 0.;
    double AlsC = model.Als4(model.getMuc());

    eta = (M_PI / AlsC * (-18. / 7. * Kpp - 12. / 11. * Kpm + 6. / 29. * Kmm + 7716. / 2233. * K7) *
            (1. - AlsC / (4. * M_PI) * 307. / 162.) + (log(model.getMuc() /
            model.getQuarks(QCD::CHARM).getMass()) - 0.25) * (3. * Kpp - 2. * Kpm + Kmm) +
            262497. / 35000. * Kpp - 123. / 625. * Kpm + 1108657. / 1305000. * Kmm - 277133. / 50750. * K7 +
            K * (-21093. / 8750. * Kpp + 13331. / 13750. * Kpm - 10181. / 18125. * Kmm - 1731104. / 2512125. * K7) +
            (log(xt) - (3. * xt) / (4. - 4. * xt) - log(xt) * (3. * xt * xt) / (4. * (1. - xt) * (1. - xt)) + 0.5) * K * K7) * xc / (model.getMyMatching()->S0(xc, xt)) * pow(AlsC, 2. / 9.);

    return (eta * (1. + model.Als(mu, FULLNLO) / 4. / M_PI * J3) * pow(model.Als(mu, FULLNLO), -2. / 9.));
}

double EvolDF2::etatt(double m) const
{
    double N = model.getNc();
    double x = model.getMyMatching()->x_t(model.getMut());
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double xm2 = (1 - x) * (1 - x);
    double xm3 = xm2 * (1 - x);
    //double xm4 = xm3 * (1 - x);
    double logx = log(x);
    //double Li2 = gsl_sf_dilog(1-x);

    double S0tt = (4. * x - 11. * x2 + x3) / 4. / xm2 -
            3. * x3 / (2. * xm3) * logx;

    double Bt = 5. * (N - 1.) / 2. / N + 3. * (N * N - 1.) / 2. / N;

    double gamma0 = 6. * (N - 1.) / N;

    double gamma1[4] = {0.}, J[4] = {0.};

    for (int i = 0; i < 4; ++i) {
        gamma1[i] = (N - 1.) / (2. * N) * (-21. + 57. / N - 19. / 3. * N + 4. / 3. * (i + 3.));
        J[i] = gamma0 * model.Beta1(3 + i) / 2. / model.Beta0(3 + i) / model.Beta0(3 + i)
                - gamma1[i] / 2. / model.Beta0(3 + i);
    }

    double b = (4. - 22. * x + 15. * x2 + 2. * x3 + x4 - 18. * x2 * logx)
            / ((-1. + x) * (-4. + 15. * x - 12. * x2 + x3 + 6. * x2 * logx));
    //double AlsT = model.Als(model.getMut());
    //double AlsB = model.Als(model.getMub());
    //double AlsC = model.Als(model.getMuc());

    double eta = pow(model.Als(model.getMuc()), 6. / 27.) *
            pow(model.Als(model.getMub()) / model.Als(model.getMuc()), 6. / 25.) *
            pow(model.Als(model.getMut()) / model.Als(model.getMub()), 6. / 23.) *
            (1. + model.Als(model.getMuc()) / 4. / M_PI * (J[1] - J[0]) +
            model.Als(model.getMub()) / 4. / M_PI * (J[2] - J[1])
            + model.Als(model.getMut()) / 4. / M_PI * (model.getMyMatching()->S1(x) / S0tt
            + Bt - J[2] + gamma0 * log(model.getMut() / model.getMuw())
            + 6 * (N * N - 1) / N * log(model.getMut() / model.getMuw()) * b));
    /* double J3 = 6. * (N - 1.) / N * (model.Beta1(3) / 2. / model.Beta0(3) / model.Beta0(3)) -
            (N - 1.) / (2. * N) * (-21. + 57. / N - 19./3. * N + 4.) / 2. / model.Beta0(3);*/

    return (eta * (1. + model.Als(m, FULLNLO) / 4. / M_PI * J[0]) * pow(model.Als(m, FULLNLO), -2. / 9.));
}

/*double EvolDF2::S1tt() const {
    double N = model.getNc();
    double x = model.getMyMatching()->x_t(model.getMut());
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double xm2 = (1 - x) * (1 - x);
    double xm3 = xm2 * (1 - x);
    double xm4 = xm3 * (1 - x);
    double logx = log(x);
    double Li2 = gsl_sf_dilog(1-x);
    
    double S0tt = (4. * x - 11. * x2 + x3) / 4. / xm2 -
                  3. * x3 / (2. * xm3) * logx;
   
    double Bt = 5. * (N - 1.) / 2. / N + 3. * (N * N - 1.) / 2. / N;
    
    double gamma0 = 6. * (N - 1.) / N;
    
    double gamma1[4] = {0.}, J[4] = {0.};
    
    for(int i = 0; i < 4; ++i){
        gamma1[i] = (N - 1.)/(2. * N) * (-21. + 57./N - 19./3. * N + 4./3. * (i + 3.));
        J[i] = gamma0 * model.Beta1(3 + i) / 2. / model.Beta0(3 + i) / model.Beta0(3 + i)
        - gamma1[i] / 2. / model.Beta0(3 + i);
    }
        
    double b = (4. - 22. * x + 15. * x2 + 2. * x3 + x4 - 18. * x2 * logx)
               / ((-1. + x) * (-4. + 15. * x - 12. * x2 + x3 + 6. * x2 * logx));
    double AlsT = model.Als(model.getMut());
    double AlsB = model.Als(model.getMub());
    double AlsC = model.Als(model.getMuc());
    
    return (pow(model.Als(model.getMuc()), 6./27.) *
            pow(model.Als(model.getMub())/model.Als(model.getMuc()), 6./25.) *
            pow(model.Als(model.getMut())/model.Als(model.getMub()), 6./23.) *
            (1. + model.Als(model.getMuc())/4./M_PI * (J[1]-J[0]) +
            model.Als(model.getMub())/4./M_PI * (J[2]-J[1])
            + model.Als(model.getMut())/4./M_PI * (model.getMyMatching()->S1(x)/S0tt
            + Bt - J[2] + gamma0 * log(model.getMut() / model.getMuw())
            + 6 * (N * N - 1) / N * log(model.getMut() / model.getMuw()) * b)));
}*/

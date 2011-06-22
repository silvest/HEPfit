/* 
 * File:   EvolDF2.cpp
 * Author: marco
 * 
 * Created on May 20, 2011, 3:55 PM
 */

#include "EvolDF2.h"

EvolDF2::EvolDF2(unsigned int dim, schemes scheme, orders order,
        const StandardModel& model) : model(model), RGEvolutor(dim, scheme, order) {
    double Nc = model.getNc();
    matrix<double> v(5, 5, 0.);
    vector<double> e(5, 0.);

    // Magic numbers in the basis of Gabbiani et al.
    
    e(0) = 6. / Nc;
    e(1) = 6. * (-1. + Nc) / Nc;
    e(2) = -6. * (-1. + Nc * Nc) / Nc;
    e(3) = (2. * (Nc + 3. * Nc * Nc - Nc * Nc * Nc + sqrt(Nc * Nc * (16. -
            11. * Nc * Nc + 4. * Nc * Nc * Nc * Nc)))) / (Nc * Nc);
    e(4) = (-2. * (-Nc - 3. * Nc * Nc + Nc * Nc * Nc + sqrt(Nc * Nc * (16. -
            11. * Nc * Nc + 4. * Nc * Nc * Nc * Nc)))) / (Nc * Nc);
    v(0, 1) = 1.;
    v(1, 3) = -(-Nc * Nc + 2. * Nc * Nc * Nc - sqrt(Nc * Nc * (16. -
            11. * Nc * Nc + 4. * Nc * Nc * Nc * Nc))) / (2. * (-2. + Nc) * Nc);
    v(1, 4) = -(-Nc * Nc + 2. * Nc * Nc * Nc + sqrt(Nc * Nc * (16. -
            11. * Nc * Nc + 4. * Nc * Nc * Nc * Nc))) / (2. * (-2. + Nc) * Nc);
    v(2, 3) = 1.;
    v(2, 4) = 1.;
    v(3, 0) = -1. / Nc;
    v(3, 2) = 1.;
    v(4, 0) = 1.;

    matrix<double> vi = v.inverse();
    for (int i = 0; i < 5; i++){
        a[i] = e(i);
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < 5; k++)
                b[i][j][k] = v(i, k) * vi(k, j);
    }
    
    matrix<double> h(5, 5, 0.);
    for (int l = 0; l < 3; l++) {
        matrix<double> gg = vi * (AnomalousDimension(NLO, 6 - l).transpose()) * v;
        double b0 = model.Beta0(6 - l);
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                h(i, j) = (i == j) * e(i) * model.Beta1(6 - l) / (2. * b0 * b0) -
                gg(i, j) / (2. * b0 + e(i) - e(j));
        matrix<double> j = v * h * vi;
        matrix<double> jv = j*v;
        matrix<double> vij = vi*j;
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++)
                for (int k = 0; k < 5; k++) {
                    c[l][i][j][k] = jv(i, k) * vi(k, j);
                    d[l][i][j][k] = -v(i, k) * vij(k, j);
                }
    }

}

EvolDF2::~EvolDF2() {
}

matrix<double> EvolDF2::AnomalousDimension(orders order, unsigned int nf) const {
    matrix<double> ad(5, 5, 0.);
    double Nc = model.getNc();
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
                throw "EvolDF2::AnomalousDimension(): wrong number of flavours";
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
            throw "EvolDF2::AnomalousDimension(): order not implemented";
    }
    return (ad);
}

matrix<double>& EvolDF2::Df2Evol(double mu, double M, orders order, schemes scheme) {
    switch (scheme) {
        case NDR:
            break;
        case LRI:
        case HV:
        default:
            std::stringstream out;
            out << scheme;
            throw "EvolDF2::Df2Evol(): scheme " + out.str() + " not implemented ";
    }

    if (mu == this->mu && M == this->M && scheme == this->scheme)
        return (*Evol(order));

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

void EvolDF2::Df2Evol(double mu, double M, double nf, schemes scheme) {

    matrix<double> resLO(5, 0.), resNLO(5, 0.), resNNLO(5, 0.);

    int l = 6 - (int) nf;
    double alsM = model.Als(M) / 4. / M_PI;
    double alsmu = model.Als(mu) / 4. / M_PI;

    double eta = alsM / alsmu;

    for (int k = 0; k < 5; k++) {
        double etap = pow(eta, a[k] / 2. / model.Beta0(nf));
        for (int i = 0; i < 5; i++)
            for (int j = 0; j < 5; j++) {
                resNNLO(i, j) += 0.;
                resNLO(i, j) += c[l][i][j][k] * etap * alsmu;
                resNLO(i, j) += d[l][i][j][k] * etap * alsM;
                resLO(i, j) += b[i][j][k] * etap;
            }
    }
    switch(order) {
        case NNLO:
            *elem[NNLO] = 0.; // Marco can implement it if he wishes to!
        case NLO:
            *elem[NLO] = (*elem[LO]) * resNLO + (*elem[NLO]) * resLO;
        case LO:
            *elem[LO] = (*elem[LO]) * resLO;
    }
}

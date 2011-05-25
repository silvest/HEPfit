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
            throw "EvolDF2::AnomalousDimension(): scheme not implemented";
    }
    return (ad);
}

matrix<double>& EvolDF2::Df2Evol(double mu, double M, orders order, schemes scheme) {
    switch (scheme) {
        case LRI:
            break;
        case NDR:
        case HV:
        default:
            std::stringstream out;
            out << scheme;
            throw "HeffDF2::Df2Evol(): scheme " + out.str() + " not implemented ";
    }
    
    if(mu == this->mu && M == this->M && scheme == this->scheme)
        return(*Evol(order));
    
    double eta = model.Als(M) / model.Als(mu);

    setScales(mu, M); // also assign evol to zero

    matrix<double> ev(5, 5, 0.);
    for (int ord = LO; ord < MAXORDER; ord++) {
        ev = 0.;
        if (ord <= this->order) {
            if (mu >= model.BelowTh(M))
                ev = Df2Evol(mu, M, model.Nf(M), orders(ord), scheme);
            else {
                double nterm = model.Nf(M) - model.Nf(mu)- 2;
                for (int ord1 = LO; ord1 <= ord; ord1++)
                    for (int ord2 = LO; ord2 <= ord1; ord2++) {
                        ev += Df2Evol(model.BelowTh(M), M, model.Nf(M), orders(ord1), scheme);
                        for (int n = model.Nf(M) - 1; n > model.Nf(mu); n--)
                            ev = Df2Evol(model.Thresholds(7 - n), model.Thresholds(8 - n),
                                n, LO, scheme) * ev; //todo
                        ev += Df2Evol(mu, model.AboveTh(mu), model.Nf(mu), orders(ord - ord1 - ord2),
                                scheme) * ev;
                    }
            }
            setEvol(ev, orders(ord));
        }
    }
}

matrix<double> EvolDF2::Df2Evol(double mu, double M, double nf, orders order, 
        schemes scheme) {
        for (int i = 0; i < 5; i++) {
            double eta = model.Als(M) / model.Als(mu);
        double etap = pow(eta, a[i]);
        for (int s = 0; s < 5; s++) {
//            switch (getOrder()) {
//                case NNLO:
//                    setEvol = 0.;
//                    *(uDF2.Evol(NNLO))(1, s) += 0.;
//                    *(uDF2.Evol(NNLO))(2, s) += 0.;
//                    *(uDF2.Evol(NNLO))(3, s) += 0.;
//                    *(uDF2.Evol(NNLO))(4, s) += 0.;
//                    *(uDF2.Evol(NNLO))(5, s) += 0.;
//                    *(uDF2.Evol(NNLO))(6, s) += 0.;
//                    *(uDF2.Evol(NNLO))(7, s) += 0.;
//                case NLO:
//                    *(uDF2.Evol(NLO))(0, s) += eta * c1[s][i] * etap;
//                    *(uDF2.Evol(NLO))(1, s) += eta * c2[s][i] * etap;
//                    *(uDF2.Evol(NLO))(2, s) += eta * c3[s][i] * etap;
//                    *(uDF2.Evol(NLO))(3, s) += eta * c4[s][i] * etap;
//                    *(uDF2.Evol(NLO))(4, s) += eta * c5[s][i] * etap;
//                    *(uDF2.Evol(NLO))(5, s) += eta * c1[s][i] * etap;
//                    *(uDF2.Evol(NLO))(6, s) += eta * c2[s][i] * etap;
//                    *(uDF2.Evol(NLO))(7, s) += eta * c3[s][i] * etap;
//                case LO:
//                    *(uDF2.Evol(LO))(0, s) += b1[s][i] * etap;
//                    *(uDF2.Evol(LO))(1, s) += b2[s][i] * etap;
//                    *(uDF2.Evol(LO))(2, s) += b3[s][i] * etap;
//                    *(uDF2.Evol(LO))(3, s) += b4[s][i] * etap;
//                    *(uDF2.Evol(LO))(4, s) += b5[s][i] * etap;
//                    *(uDF2.Evol(LO))(5, s) += b1[s][i] * etap;
//                    *(uDF2.Evol(LO))(6, s) += b2[s][i] * etap;
//                    *(uDF2.Evol(LO))(7, s) += b3[s][i] * etap;
//            }

        }
    }

}


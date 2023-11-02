/* 
 * File:   EvolDF1_diujlknu.cpp
 * Author: silvest
 * 
 * Created on 31 ottobre 2023, 18:07
 */

#include "EvolDF1_diujlknu.h"
#include <cstring>
#include "StandardModel.h"


EvolDF1_diujlknu::EvolDF1_diujlknu(unsigned int dim_i, schemes scheme, orders order, const StandardModel& model) :   
RGEvolutor(dim_i, scheme, order),
    model(model),
    dim(dim_i){

    gslpp::matrix<double> g0t(AnomalousDimension(LO, 3).transpose());

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
    for (int l = 0; l < 2; l++) {
        gslpp::matrix<double> gg = vi * (AnomalousDimension(NLO, 5 - l).transpose()) * v;
        double b0 = model.Beta0(5 - l);
        for (unsigned int i = 0; i < dim; i++)
            for (unsigned int j = 0; j < dim; j++)
                h(i, j) = (i == j) * e(i) * model.Beta1(5 - l) / (2. * b0 * b0) -
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

EvolDF1_diujlknu::~EvolDF1_diujlknu() {
}

gslpp::matrix<double> EvolDF1_diujlknu::AnomalousDimension(orders order, unsigned int nf) const
{
    gslpp::matrix<double> ad(dim, dim, 0.);
    double Nc = model.getNc();
            switch (order) {
                case LO:
                    ad(2, 2) = -4.;
                    ad(3, 3) = -4.;
                    ad(4, 4) = 4./3.;
                    break;
                case NLO:
                    ad(2, 2) = 2./9.*(-303.+10.*nf);
                    ad(3, 3) = ad(2,2);
                    ad(4, 4) = 2./27.*(543.-26.*nf);
                    break;
                default:
                    throw std::runtime_error("EvolDF1_diujlknu::AnomalousDimension(): order not implemented");
            }
    return (ad);
}

gslpp::matrix<double> EvolDF1_diujlknu::AnomalousDimension(orders_qed order_qed, unsigned int nf) const
{
    gslpp::matrix<double> ad(dim, dim, 0.);
    double Nc = model.getNc();
            switch (order_qed) {
                case LO_QED:
                    ad(0, 0) = -2.;
                    ad(1, 1) = -1.;
                    ad(2, 2) = 2./3.;
                    ad(3, 3) = 2./3.;
                    ad(3, 4) = 1./12.;
                    ad(4, 3) = 4.;
                    ad(4, 4) = -20./9.;
                    break;
                default:
                    throw std::runtime_error("EvolDF1_diujlknu::AnomalousDimension(): order_qed not implemented");
            }
    return (ad);
}

gslpp::matrix<double>& EvolDF1_diujlknu::Df1diujlknuEvol(double mu, double M, orders order, schemes scheme)
{
    switch (scheme) {
        case NDR:
            break;
        case LRI:
        case HV:
        default:
            std::stringstream out;
            out << scheme;
            throw std::runtime_error("EvolDF1_diujlknu::Df1diujlknuEvol(): scheme " + out.str() + " not implemented ");
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
    if (M != mu) {
        double m_down = mu;
        double m_up = model.AboveTh(m_down);
        double nf = model.Nf(m_down);

        while (m_up < M) {
            Df1diujlknuEvol(m_down, m_up, nf, scheme);
            m_down = m_up;
            m_up = model.AboveTh(m_down);
            nf += 1.;
        }
        Df1diujlknuEvol(m_down, M, nf, scheme);
    }
    
    //add the leading QED log without resumming it for the moment, to be improved later
    *elem[LO] = (*elem[LO]) + AnomalousDimension(LO_QED,3).transpose()*model.getAle()/2./M_PI*log(mu/M);
    
    return (*Evol(order));
}

void EvolDF1_diujlknu::Df1diujlknuEvol(double mu, double M, double nf, schemes scheme)
{

    gslpp::matrix<double> resLO(dim, 0.), resNLO(dim, 0.), resNNLO(dim, 0.);

    int l = 5 - (int) nf;
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
            throw std::runtime_error("Error in EvolDF1_diujlknu::Df1diujlknuEvol()");
    }
}

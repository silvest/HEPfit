/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <gsl/gsl_sf_zeta.h>
#include "EvolDC1Buras.h"
#include <stdexcept>


EvolDC1Buras::EvolDC1Buras(unsigned int dim_i, schemes scheme, orders order, const StandardModel& model) 
:           RGEvolutor(dim_i, scheme, order), model(model),
            v(dim_i,0.), vi(dim_i,0.), js(dim_i,0.), h(dim_i,0.), gg(dim_i,0.), s_s(dim_i,0.),
            jssv(dim_i,0.), jss(dim_i,0.), jv(dim_i,0.), vij(dim_i,0.), e(dim_i,0.), dim(dim_i)  
{
    
    /*magic numbers a & b */
    
    for(int L=2; L>-1; L--){
        
    /* L=2 --> u,d,s,c (nf=4)  L=1 --> u,d,s,c,b (nf=5) L=0 --> u,d,s,c,b,t (nf=6)*/
        
    nu = L;  nd = L;
    if(L == 1){nd = 3; nu = 2;} 
    if(L == 0){nd = 3; nu = 3;}
    
    AnomalousDimension_DC1_Buras(LO,nu,nd).transpose().eigensystem(v,e);
    vi = v.inverse();
    for(unsigned int i = 0; i < dim; i++){
       a[L][i] = e(i).real();
       for (unsigned int j = 0; j < dim; j++) {
           for (unsigned int k = 0; k < dim; k++)  {
                b[L][i][j][k] = v(i, k).real() * vi(k, j).real();
               }
           }
       }
    
    // LO evolutor in the standard basis
    
    gg = vi * AnomalousDimension_DC1_Buras(NLO,nu,nd).transpose() * v;
    double b0 = model.Beta0(6-L);
    double b1 = model.Beta1(6-L);
    for (unsigned int i = 0; i < dim; i++){
        for (unsigned int j = 0; j < dim; j++){
            s_s.assign( i, j, (b1 / b0) * (i==j) * e(i).real() - gg(i,j));    
            if(fabs(e(i).real() - e(j).real() + 2. * b0)>0.00000000001){
                h.assign( i, j, s_s(i,j) / (2. * b0 + e(i) - e(j)));
                }
            }
        }
    js = v * h * vi;
    jv = js * v;
    vij = vi * js;
    jss = v * s_s * vi;
    jssv = jss * v;        
    for (unsigned int i = 0; i < dim; i++){
        for (unsigned int j = 0; j < dim; j++){
            if(fabs(e(i).real() - e(j).real() + 2. * b0) > 0.00000000001){
                for(unsigned int k = 0; k < dim; k++){
                        c[L][i][j][k] = jv(i, k).real() * vi(k, j).real();
                        d[L][i][j][k] = -v(i, k).real() * vij(k, j).real();
                        }
                    }
            else{    
                for(unsigned int k = 0; k < dim; k++){
                   c[L][i][j][k] = (1./(2. * b0)) * jssv(i, k).real() * vi(k, j).real();
                   d[L][i][j][k] = 0.;
                   }   
                }
            }
        }
    }
}
    
EvolDC1Buras::~EvolDC1Buras() 
{}

gslpp::matrix<double> EvolDC1Buras::AnomalousDimension_DC1_Buras(orders order, unsigned int n_u, unsigned int n_d) const
{
   
    /* anomalous dimension related to Delta F = 1 operators in Buras basis, hep-ph/9512380v1 */
    
    /* gamma(row, column) at the LO */
    
    unsigned int nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    gslpp::matrix<double> gammaDF1(dim, 0.);
    
    switch(order){
        
    case LO:
    
    gammaDF1(0,0) = -2.;
    gammaDF1(0,1) = 6. ;    
    
    gammaDF1(1,0) = 6.;
    gammaDF1(1,1) = -2.;
    
    if(nf<5){ 
    gammaDF1(1,2) = -2./9.;
    gammaDF1(1,3) = 2./3.;
    gammaDF1(1,4) = -2./9.;
    gammaDF1(1,5) = 2./3.;
    
    }
    
    gammaDF1(2,2) = -22./9.;
    gammaDF1(2,3) = 22./3.;
    gammaDF1(2,4) = -4./9.;
    gammaDF1(2,5) = 4./3.;
    
    gammaDF1(3,2) = 6.-2./9.*nf;
    gammaDF1(3,3) = -2.+2./3.*nf;
    gammaDF1(3,4) = -2./9.*nf;
    gammaDF1(3,5) = 2./3.*nf;
    
    gammaDF1(4,4) = 2.;
    gammaDF1(4,5) = -6.;
    
    gammaDF1(5,2) = -2./9.*nf;
    gammaDF1(5,3) = 2./3.*nf;
    gammaDF1(5,4) = -2./9.*nf;
    gammaDF1(5,5) = -16.+2./3.*nf;
   
    break;    
    
    case NLO:
        
    if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolDF1nlep::AnomalousDimension_B(): wrong number of flavours"); 
        }
    
    /* gamma(row, column) at the NLO */
    
    gammaDF1(0,0) = -21./2.-2./9.*nf;
    gammaDF1(0,1) = 7./2.+2./3.*nf;
    
    gammaDF1(1,0) = 7./2.+2./3.*nf;
    gammaDF1(1,1) = -21./2.-2./9.*nf;
    
    if(nf<5){    
        
    gammaDF1(0,2) = 79./9.;
    gammaDF1(0,3) = -7./3.;
    gammaDF1(0,4) = -65./9.;
    gammaDF1(0,5) = -7./3.;
    
    gammaDF1(1,2) = -202./243.;
    gammaDF1(1,3) = 1354./81.;
    gammaDF1(1,4) = -1192./243.;
    gammaDF1(1,5) = 904./81.;
    
    }
    
    gammaDF1(2,2) = -5911./486.+71./9.*nf;
    gammaDF1(2,3) = 5983./162.+1./3.*nf;
    gammaDF1(2,4) = -2384./243.-71./9.*nf;
    gammaDF1(2,5) = 1808./81.-1./3.*nf;
    
    
    gammaDF1(3,2) = 379./18.+56./243.*nf;
    gammaDF1(3,3) = -91./6.+808./81.*nf;
    gammaDF1(3,4) = -130./9.-502./243.*nf;
    gammaDF1(3,5) = -14./3.+646./81.*nf;
    
    gammaDF1(4,2) = -61./9.*nf;
    gammaDF1(4,3) = -11./3.*nf;
    gammaDF1(4,4) = 71./3.+61./9.*nf;
    gammaDF1(4,5) = -99.+11./3.*nf;
    
    gammaDF1(5,2) = -682./243.*nf;
    gammaDF1(5,3) = 106./81.*nf;
    gammaDF1(5,4) = -225./2.+1676./243.*nf;
    gammaDF1(5,5) = -1343./6.+1348./81.*nf;
        
    break;
    
    default:
            throw std::runtime_error("EvolDF1nlep::AnomalousDimensio_B_S(): order not implemented"); 
    }
   
    return (gammaDF1);
    
  }

gslpp::matrix<double>& EvolDC1Buras::DC1EvolBuras(double mu, double M, orders order, schemes scheme) 
{
    switch (scheme) {
        case NDR:
            break;
        case LRI:
        case HV:
        default:
            std::stringstream out;
            out << scheme;
            throw std::runtime_error("EvolDC1::Df1EvolDC1(): scheme " + out.str() + " not implemented "); 
    }

    double alsMZ = model.getAlsMz();
    double Mz = model.getMz();
    if(alsMZ == alsMZ_cache && Mz == Mz_cache) {
        if (mu == this-> mu && M == this->M && scheme == this->scheme)
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
        DC1EvolBuras(m_down, m_up, nf, scheme);
        DC1PenguinThresholds(m_up, order);
        m_down = m_up;
        m_up = model.AboveTh(m_down);
        nf += 1.;
    }
    DC1EvolBuras(m_down, M, nf, scheme);
    return (*Evol(order));
    
    }
    
void EvolDC1Buras::DC1EvolBuras(double mu, double M, double nf, schemes scheme) 
{

    gslpp::matrix<double> resLO(dim, 0.), resNLO(dim, 0.), resNNLO(dim, 0.);

    int L = 6 - (int) nf;
    double alsM = model.Als(M) / 4. / M_PI;
    double alsmu = model.Als(mu) / 4. / M_PI;
    
    double eta = alsM / alsmu;
    
    for (unsigned int k = 0; k < dim; k++) {
        double etap = pow(eta, a[L][k] / 2. / model.Beta0(nf));
        for (unsigned int i = 0; i < dim; i++){
            for (unsigned int j = 0; j < dim; j++) {
                resNNLO(i, j) += 0.;
                
                if(fabs(e(i).real() - e(j).real() + 2. * model.Beta0(nf))>0.000000000001)  {
                    resNLO(i, j) += c[L][i][j][k] * etap * alsmu;
                    resNLO(i, j) += d[L][i][j][k] * etap * alsM;
                }
                else{
                    resNLO(i, j) += - c[L][i][j][k] * etap * alsmu * log(eta);    
                }        
                resLO(i, j) += b[L][i][j][k] * etap;
            }
        }
    }
    switch(order) {
        case NNLO:
            *elem[NNLO] = 0.;
        case NLO:
            *elem[NLO] = (*elem[LO]) * resNLO + (*elem[NLO]) * resLO;
        case LO:
            *elem[LO] = (*elem[LO]) * resLO;
            break;
        case FULLNNLO:
        case FULLNLO:
        default:
            throw std::runtime_error("Error in EvolDC1Buras::DC1EvolBuras()");
    } 
    
  }
 
gslpp::matrix<double> EvolDC1Buras::StrongThresholds() const
{

// entries of the threshold matrix for the evolution at the NLO
    
gslpp::matrix<double> deltarsT(dim,0.);

deltarsT(2,3) = 5./27.;
deltarsT(2,5) = 5./27.;
deltarsT(3,3) = -5./9.;
deltarsT(3,5) = -5./9.;
deltarsT(4,3) = 5./27.;
deltarsT(4,5) = 5./27.;
deltarsT(5,3) = -5./9.;
deltarsT(5,5) = -5./9.;

return(deltarsT);

}

void EvolDC1Buras::DC1PenguinThresholds(double M, orders order) 
{

    double alsM = model.Als(M) / 4. / M_PI;
    gslpp::matrix<double> drsT(dim,0.);
    drsT = alsM * StrongThresholds();
    *elem[NLO] = (*elem[NLO]) + (*elem[LO]) * drsT;  
  }

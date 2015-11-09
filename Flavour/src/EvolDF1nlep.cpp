/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EvolDF1nlep.h"

EvolDF1nlep::EvolDF1nlep(unsigned int dim_i, schemes scheme, orders order, orders_ew order_ew, const StandardModel& model)
:   RGEvolutor(dim_i, scheme, order, order_ew), model(model), V(dim_i,0.), Vi(dim_i,0.),
    gs(dim_i,0.), Js(dim_i,0.), ge0(dim_i,0.), K0(dim_i,0.), ge11(dim_i,0.), K11(dim_i,0.),
    JsK0V(dim_i,0.), ViK0Js(dim_i,0.), Gamma_s0T(dim_i,0.), Gamma_s1T(dim_i,0.), 
    Gamma_eT(dim_i,0.), Gamma_seT(dim_i,0.), JsV(dim_i,0.), ViJs(dim_i,0.), K0V(dim_i,0.), 
    ViK0(dim_i,0.), K11V(dim_i,0.), ViK11(dim_i,0.), ge11sing(dim_i,0.), K11sing(dim_i,0.), 
    K11singV(dim_i,0.), e(dim_i,0.), dim(dim_i) 
{
    
    int nu = 0, nd = 0;
    double  b0 = 0., b1 = 0.;
    
    /* L=3 --> u,d,s,c (nf=3) L=2 --> u,d,s,c (nf=4)  L=1 --> u,d,s,c,b (nf=5) L=0 --> u,d,s,c,b,t (nf=6)*/
    for(int L=3; L>-1; L--){
        
        b0 = model.Beta0(6-L);
        b1 = model.Beta1(6-L);
        
	if(L == 3){nd = 2; nu = 1;} 
        if(L == 2){nd = 2; nu = 2;}
        if(L == 1){nd = 3; nu = 2;} 
        if(L == 0){nd = 3; nu = 3;}
        
        Gamma_s0T = AnomalousDimension_nlep_S(LO,nu,nd).transpose();
        Gamma_s1T = AnomalousDimension_nlep_S(NLO,nu,nd).transpose();
        Gamma_eT = AnomalousDimension_nlep_EM(LO,nu,nd).transpose();
        Gamma_seT = AnomalousDimension_nlep_EM(NLO,nu,nd).transpose();
        
        AnomalousDimension_nlep_S(LO,nu,nd).transpose().eigensystem(V,e);
        Vi = V.inverse();
        
        /* magic numbers of U0 */
        for(unsigned int i = 0; i < dim; i++){
            a[L][i] = e(i).real()/2./b0;
            for (unsigned int j = 0; j < dim; j++){
                for (unsigned int k = 0; k < dim; k++){
                    b[L][i][j][k] = V(i, k).real() * Vi(k, j).real();
                }
            }
        }
    
        gs = (b1/2./b0/b0) * Vi * Gamma_s0T * V - (1./2./b0) * Vi * Gamma_s1T * V;
        for(unsigned int i = 0; i<dim ; i++){
            for(unsigned int j = 0; j<dim ; j++){  
                gs.assign( i , j, gs(i,j)/(1. + a[L][i] - a[L][j]));
            }
        }
        Js = V * gs * Vi;
        
        /*magic numbers related to Js*/
        JsV = Js*V;
        ViJs = Vi * Js;
        for(unsigned int i = 0; i<dim; i++){
            for(unsigned int j = 0; j<dim; j++){
                for(unsigned int k = 0; k<dim; k++){
                    c[L][i][j][k] = JsV(i, k).real() * Vi(k, j).real();
                    d[L][i][j][k] = -V(i, k).real() * ViJs(k, j).real();
                }
            }
        }
        
        ge0 = (1./2./b0) *  Vi * Gamma_eT * V;
        for(unsigned int i = 0; i<dim ; i++){
            for(unsigned int j = 0; j<dim ; j++){
                ge0.assign( i , j, ge0(i,j)/(1. - a[L][i] + a[L][j]));
            }
        }
        K0 = V * ge0 * Vi;
        
        /*magic numbers related to K0*/
        K0V = K0*V;
        ViK0 = Vi * K0;
        for(unsigned int i = 0; i<dim; i++){
            for(unsigned int j = 0; j<dim; j++){
                for(unsigned int k = 0; k<dim; k++){
                    m[L][i][j][k] = K0V(i, k).real() * Vi(k, j).real();
                    n[L][i][j][k] = -V(i, k).real() * ViK0(k, j).real();
                }
            }
        }
        
        ge11 = Gamma_seT - (b1/b0) * Gamma_eT + Gamma_eT * Js - Js * Gamma_eT;
        ge11 = Vi * ge11;
        ge11 = ge11 * V;
        for(unsigned int i = 0; i<dim ; i++){
            for(unsigned int j = 0; j<dim ; j++){
                if(fabs(a[L][j]-a[L][i])> 0.00000000001){
                    ge11.assign( i , j, ge11(i,j)/( 2. * b0 * (a[L][j] - a[L][i])));
                }
                else{
                    ge11sing.assign( i, j, ge11(i,j)/2./b0);
                    ge11.assign( i , j, 0.);
                }
            }
        }
        K11 = V * ge11 * Vi;
        K11sing = V * ge11sing * Vi;
        /*magic numbers related to K11*/
        K11V = K11 * V;
        ViK11 = Vi * K11;
        K11singV = K11sing * V;
        if(L==1){
        }
        for(unsigned int i = 0; i<dim ; i++){
            for(unsigned int j = 0; j<dim ; j++){
                    for(unsigned int k = 0; k<dim ; k++){
                        o[L][i][j][k] = K11V(i, k).real() * Vi(k, j).real();
                        p[L][i][j][k] = -V(i, k).real() * ViK11(k, j).real();
                        u[L][i][j][k] = K11singV(i, k).real() * Vi(k, j).real();
                    }
                }    
            }
        
        /*magic numbers related to K12 and K13*/
        JsK0V = Js * K0 * V; 
        ViK0Js = Vi * K0 * Js;
        for(unsigned int i = 0; i<dim ; i++){
            for(unsigned int j = 0; j<dim ; j++){
                for(unsigned int k=0; k<dim; k++){
                    q[L][i][j][k] =  JsK0V(i, k).real() * Vi(k, j).real();
                    r[L][i][j][k] =  V(i, k).real() * ViK0Js(k, j).real();
                    s[L][i][j][k] = -JsV(i, k).real() * ViK0(k, j).real();
                    t[L][i][j][k] = -K0V(i, k).real() * ViJs(k, j).real();
                }
            }
        }
    }        
}
    
EvolDF1nlep::~EvolDF1nlep() 
{}

gslpp::matrix<double> EvolDF1nlep::AnomalousDimension_nlep_S(orders order, unsigned int n_u, unsigned int n_d) const
{
   
    /* anomalous dimension related to Delta F = 1 operators in Buras basis, hep-ph/9512380v1 */
    
    /*gamma(row, column) leading order*/
    
    unsigned int nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    gslpp::matrix<double> gammaDF1(dim, 0.);
   
    switch(order){
        
    case LO:
        
    gammaDF1(0,0) = -2.;
    gammaDF1(0,1) = 6. ;
    
    
    gammaDF1(1,0) = 6.;
    gammaDF1(1,1) = -2.;
    gammaDF1(1,2) = -2./9.;
    gammaDF1(1,3) = 2./3.;
    gammaDF1(1,4) = -2./9.;
    gammaDF1(1,5) = 2./3.;
    
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
    
    gammaDF1(6,6) = 2.;   
    gammaDF1(6,7) = -6.;
    
    gammaDF1(7,2) = -2./9.*(n_u-n_d/2.);
    gammaDF1(7,3) = 2./3.*(n_u-n_d/2.);
    gammaDF1(7,4) = -2./9.*(n_u-n_d/2.);
    gammaDF1(7,5) = 2./3.*(n_u-n_d/2.);
    gammaDF1(7,7) = -16.;
    
    gammaDF1(8,2) = 2./9.;
    gammaDF1(8,3) = -2./3.;
    gammaDF1(8,4) = 2./9.;
    gammaDF1(8,5) = -2./3.;
    gammaDF1(8,8) = -2.; 
    gammaDF1(8,9) = 6.;
    
    gammaDF1(9,2) = -2./9.*(n_u-n_d/2.);
    gammaDF1(9,3) = 2./3.*(n_u-n_d/2.);
    gammaDF1(9,4) = -2./9.*(n_u-n_d/2.);
    gammaDF1(9,5) = 2./3.*(n_u-n_d/2.);
    gammaDF1(9,8) = 6.;
    gammaDF1(9,9) = -2.; 
    
    break;    
    
    case NLO:
        
    if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolDF1nlep::AnomalousDimension_nlep_S("
                "orders order, unsigned int n_u, unsigned int n_d) " " wrong number of flavour"); 
        }
    
    /*gamma(riga, colonna) next to leading order*/
        
    gammaDF1(0,0) = -21./2.-2./9.*nf;
    gammaDF1(0,1) = 7./2.+2./3.*nf;
    gammaDF1(0,2) = 79./9.;
    gammaDF1(0,3) = -7./3.;
    gammaDF1(0,4) = -65./9.;
    gammaDF1(0,5) = -7./3.;
    
    
    gammaDF1(1,0) = 7./2.+2./3.*nf;
    gammaDF1(1,1) = -21./2.-2./9.*nf;
    gammaDF1(1,2) = -202./243.;
    gammaDF1(1,3) = 1354./81.;
    gammaDF1(1,4) = -1192./243.;
    gammaDF1(1,5) = 904./81.;
    
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
    
    gammaDF1(6,2) = -61./9.*(n_u-n_d/2.);
    gammaDF1(6,3) = -11./3.*(n_u-n_d/2.);
    gammaDF1(6,4) = 83./9.*(n_u-n_d/2.);
    gammaDF1(6,5) = -11./3.*(n_u-n_d/2.);
    gammaDF1(6,6) = 71./3.-22./9.*nf;   
    gammaDF1(6,7) = -99.+22./3.*nf;
    
    gammaDF1(7,2) = -682./243*(n_u-n_d/2.);
    gammaDF1(7,3) = 106./81.*(n_u-n_d/2.);
    gammaDF1(7,4) = 704./243.*(n_u-n_d/2.);
    gammaDF1(7,5) = 736./81.*(n_u-n_d/2.);
    gammaDF1(7,6) = -225./2.+4*nf;
    gammaDF1(7,7) = -1343./6.+68./9.*nf;
    
    gammaDF1(8,2) = 202./243.+73./9.*(n_u-n_d/2.);
    gammaDF1(8,3) = -1354./81.-1./3.*(n_u-n_d/2.);
    gammaDF1(8,4) = 1192./243.-71./9.*(n_u-n_d/2.);
    gammaDF1(8,5) = -904./81.-1./3.*(n_u-n_d/2.);
    gammaDF1(8,8) = -21./2.-2./9.*nf;
    gammaDF1(8,9) = 7./2.+2./3.*nf;
    
    gammaDF1(9,2) = -79./9.-106./243.*(n_u-n_d/2.);
    gammaDF1(9,3) = 7./3.+826./81.*(n_u-n_d/2.);
    gammaDF1(9,4) = 65./9.-502./243.*(n_u-n_d/2.);
    gammaDF1(9,5) = 7./3.+646./81.*(n_u-n_d/2.);
    gammaDF1(9,8) = 7./2.+2./3.*nf;
    gammaDF1(9,9) = -21./2.-2./9.*nf; 
    
    break;
    
    default:
        std::stringstream out;
        out << order;
        throw std::runtime_error("EvolDF1nlep::AnomalousDimension_nlep_S("
                "orders order, unsigned int n_u, unsigned int n_d) " 
                + out.str() + " not implemented"); 
        
    }
   
    return (gammaDF1);
    
  }

gslpp::matrix<double> EvolDF1nlep::AnomalousDimension_nlep_EM(orders order, unsigned int n_u, unsigned int n_d) const
{
   
    /* anomalous dimension related to Buras operators hep-ph/9512380v1 */
    /*gamma(riga, colonna) leading order*/
    unsigned int nf = n_u + n_d; /*n_u\d = active type up/down flavor d.o.f.*/
    gslpp::matrix<double> gammaDF1(dim, 0.);
   
    switch(order){
        
    case LO:
        
    gammaDF1(0,0) = -8./3.;
    gammaDF1(0,6) = 16./9.;
    gammaDF1(0,8) = 16./9.;
    
    gammaDF1(1,1) = -8./3.;
    gammaDF1(1,6) = 16./27.;
    gammaDF1(1,8) = 16./27.;
    
    gammaDF1(2,6) = -16./27.+16./9.*(n_u-n_d/2.);
    gammaDF1(2,8) = -88./27.+16./9.*(n_u-n_d/2.);
    
    gammaDF1(3,6) = -16./9.+16./27.*(n_u-n_d/2.);
    gammaDF1(3,8) = -16./9.+16./27.*(n_u-n_d/2.);
    gammaDF1(3,9) = -8./3.;
    
    gammaDF1(4,6) = 8./3.+16./9.*(n_u-n_d/2.);
    gammaDF1(4,8) = 16./9.*(n_u-n_d/2.);
    
    gammaDF1(5,6) = 16./27.*(n_u-n_d/2.);
    gammaDF1(5,7) = 8./3.;
    gammaDF1(5,8) = 16./27.*(n_u-n_d/2.);
    
    gammaDF1(6,4) = 4./3.;
    gammaDF1(6,6) = 4./3.+16./9.*(n_u+n_d/4.);   
    gammaDF1(6,8) = 16./9.*(n_u+n_d/4.);
    
    gammaDF1(7,5) = 4./3.;
    gammaDF1(7,6) = 16./27.*(n_u+n_d/4.);
    gammaDF1(7,7) = 4./3.;
    gammaDF1(7,8) = 16./27.*(n_u+n_d/4.);
    
    gammaDF1(8,2) = -4./3.;
    gammaDF1(8,6) = 8./27.+16./9.*(n_u+n_d/4.);
    gammaDF1(8,8) = -28./27.+16./9.*(n_u+n_d/4.); 
    
    gammaDF1(9,3) = -4./3.;
    gammaDF1(9,6) = 8./9.+16./27.*(n_u+n_d/4.);
    gammaDF1(9,8) = 8./9.+16./27.*(n_u+n_d/4.);
    gammaDF1(9,9) = -4./3.; 
    
    break;    
    
    case NLO:
        
     if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
               throw std::runtime_error("EvolDF1nlep::AnomalousDimension_nlep_EM("
                "orders order, unsigned int n_u, unsigned int n_d) " " wrong number of flavour"); 
      }
    
    /*gamma(riga, colonna) next to leading order*/
        
    gammaDF1(0,0) = 194./9.;
    gammaDF1(0,1) = -2./3.;
    gammaDF1(0,2) = -88./243.;
    gammaDF1(0,3) = 88./81.;
    gammaDF1(0,4) = -88./243.;
    gammaDF1(0,5) = 88./81.;
    gammaDF1(0,6) = 152./27.;
    gammaDF1(0,7) = 40./9.;
    gammaDF1(0,8) = 136./27.;
    gammaDF1(0,9) = 56./9.;
    
    gammaDF1(1,0) = 25./3.;
    gammaDF1(1,1) = -49./9.;
    gammaDF1(1,2) = -556./729.;
    gammaDF1(1,3) = 556./243.;
    gammaDF1(1,4) = -556./729.;
    gammaDF1(1,5) = 556./243.;
    gammaDF1(1,6) = -484./729.;
    gammaDF1(1,7) = -124./27.;
    gammaDF1(1,8) = -3148./729.;
    gammaDF1(1,9) = 172./27.;
    
    gammaDF1(2,2) = 1690./729.-136./243.*(n_u-n_d/2.);
    gammaDF1(2,3) = -1690./243.+136./81.*(n_u-n_d/2.);
    gammaDF1(2,4) = 232./729.-136./243.*(n_u-n_d/2.);
    gammaDF1(2,5) = -232./243.+136./81.*(n_u-n_d/2.);
    gammaDF1(2,6) = 3136./729.+104./27.*(n_u-n_d/2.);
    gammaDF1(2,7) = 64./27.+88./9.*(n_u-n_d/2.);
    gammaDF1(2,8) = 20272./729.+184./27.*(n_u-n_d/2.);
    gammaDF1(2,9) = -112./27.+8./9.*(n_u-n_d/2.);
    
    gammaDF1(3,2) = -641./243.-388./729.*n_u+32./729.*n_d;
    gammaDF1(3,3) = -655./81.+388./243.*n_u-32./243.*n_d;
    gammaDF1(3,4) = 88./243.-388./729*n_u+32./729.*n_d;
    gammaDF1(3,5) = -88./81.+388./243.*n_u-32./243.*n_d;
    gammaDF1(3,6) = -152./27.+3140./729.*n_u+656./729.*n_d;
    gammaDF1(3,7) = -40./9.-100./27.*n_u-16./27.*n_d;
    gammaDF1(3,8) = 170./27.+908./729.*n_u+1232./729.*n_d;
    gammaDF1(3,9) = -14./3.+148./27.*n_u-80./27*n_d;
    
    gammaDF1(4,2) = -136./243.*(n_u-n_d/2.);
    gammaDF1(4,3) = 136./81.*(n_u-n_d/2.);
    gammaDF1(4,4) = -2.-136./243.*(n_u-n_d/2.);
    gammaDF1(4,5) = 6.+136./81.*(n_u-n_d/2.);
    gammaDF1(4,6) = -232./9.+104./27.*(n_u-n_d/2.);
    gammaDF1(4,7) = 40./3.+88./9.*(n_u-n_d/2.);
    gammaDF1(4,8) = 184./27.*(n_u-n_d/2.);
    gammaDF1(4,9) = 8./9.*(n_u-n_d/2.);
    
    gammaDF1(5,2) = -748./729.*n_u+212./729.*n_d;
    gammaDF1(5,3) = 748./243.*n_u-212./243.*n_d;        
    gammaDF1(5,4) = 3.-748./729.*n_u+212./729.*n_d;
    gammaDF1(5,5) = 7.+748./243.*n_u-212./243.*n_d;
    gammaDF1(5,6) = -2.-5212./729.*n_u+4832./729.*n_d;
    gammaDF1(5,7) = 182./9.+188./27.*n_u-160./27.*n_d;
    gammaDF1(5,8) = -2260./729.*n_u+2816./729.*n_d;
    gammaDF1(5,9) = -140./27.*n_u+64./27.*n_d;
   
    gammaDF1(6,2) = -136./243.*(n_u+n_d/4.);
    gammaDF1(6,3) = 136./81.*(n_u+n_d/4.);
    gammaDF1(6,4) = -116./9.-136./243.*(n_u+n_d/4.);
    gammaDF1(6,5) = 20./3.+136./81.*(n_u+n_d/4.);
    gammaDF1(6,6) = -134./9.+104./27.*(n_u+n_d/4.);   
    gammaDF1(6,7) = 38./3.+88./9.*(n_u+n_d/4.);
    gammaDF1(6,8) = 184./27.*(n_u+n_d/4.);
    gammaDF1(6,9) = 8./9.*(n_u+n_d/4.);
    
    gammaDF1(7,2) = -748./729.*n_u-106./729.*n_d; 
    gammaDF1(7,3) = 748./243.*n_u+106./243.*n_d;
    gammaDF1(7,4) = -1.-748./729.*n_u-106./729.*n_d;
    gammaDF1(7,5) = 91./9.+748./243.*n_u+106./243.*n_d;
    gammaDF1(7,6) = 2.-5212./729.*n_u-2416./729.*n_d;
    gammaDF1(7,7) = 154./9.+188./27.*n_u+80./27.*n_d;
    gammaDF1(7,8) = -2260./729.*n_u-1408./729.*n_d;
    gammaDF1(7,9) = -140./27.*n_u-32./27.*n_d;
    
    gammaDF1(8,2) = 7012./729.-136./243.*(n_u+n_d/4.);
    gammaDF1(8,3) = 764./243.+136./81.*(n_u+n_d/4.);
    gammaDF1(8,4) = -116./729.-136./243.*(n_u+n_d/4.);
    gammaDF1(8,5) = 116./243.+136./81.*(n_u+n_d/4.);
    gammaDF1(8,6) = -1568./729.+104./27.*(n_u+n_d/4.);
    gammaDF1(8,7) = -32./27.+88./9.*(n_u+n_d/4.);
    gammaDF1(8,8) = 5578./729.+184./27.*(n_u+n_d/4.);
    gammaDF1(8,9) = 38./27.+8./9.*(n_u+n_d/4.);
    
    gammaDF1(9,2) = 1333./243.-388./729.*n_u-16./729.*n_d;
    gammaDF1(9,3) = 107./81.+388./243.*n_u+16./243.*n_d;
    gammaDF1(9,4) = -44./243.-388./729.*n_u-16./729.*n_d;
    gammaDF1(9,5) = 44./81.+388./243.*n_u+16./243.*n_d; 
    gammaDF1(9,6) = 76./27.+3140./729.*n_u-328./729.*n_d;
    gammaDF1(9,7) = 20./9.-100./27.*n_u+8./27.*n_d;
    gammaDF1(9,8) = 140./27.+908./729.*n_u-616./729.*n_d;
    gammaDF1(9,9) = -28./9.+148./27.*n_u+40./27.*n_d; 
    
    break;
    
    default:
        std::stringstream out;
        out << order;
        throw std::runtime_error("EvolDF1nlep::AnomalousDimension_nlep_EM("
                "orders order, unsigned int n_u, unsigned int n_d) " 
                + out.str() + " not implemented"); 
        
    }
   
    return (gammaDF1);
    
}


gslpp::matrix<double> EvolDF1nlep::Df1threshold_deltarsT(double nf) const 
{
    
    gslpp::matrix <double> delta_rsT(dim,0.); 
  
    delta_rsT(2,3) = 5./27.;
    delta_rsT(2,5) = 5./27.;
    delta_rsT(3,3) = -5./9.;
    delta_rsT(4,5) = -5./9.;
    delta_rsT(4,3) = 5./27.;
    delta_rsT(4,5) = 5./27.;
    delta_rsT(5,3) = -5./9.;
    delta_rsT(5,5) = -5./9.;
    
    if(nf == 3. || nf == 5.){
      
    delta_rsT(2,7) = -5./54.;
    delta_rsT(2,9) = -5./54.;    
    delta_rsT(3,7) = 5./18.;
    delta_rsT(3,9) = 5./18.;
    delta_rsT(4,7) = -5./54.; 
    delta_rsT(4,9) = -5./54.;    
    delta_rsT(5,7) = 5./18.;    
    delta_rsT(5,9) = 5./18.;
    }   

    
    
    else {
      
    delta_rsT(2,7) = 5./27.;
    delta_rsT(2,9) = 5./27.;
    delta_rsT(3,7) = -5./9.;
    delta_rsT(3,9) = -5./9.;
    delta_rsT(4,7) = 5./27.;
    delta_rsT(4,9) = 5./27.;
    delta_rsT(5,7) = -5./9.;   
    delta_rsT(5,9) = -5./9.;
   
    }  

    return(delta_rsT);

}

gslpp::matrix<double> EvolDF1nlep::Df1threshold_deltareT(double nf) const 
{
 
    gslpp::matrix<double> delta_reT(dim,0.);    
    
    if(nf == 3. || nf == 5.){
        
    delta_reT(6,2) = 20./27.;
    delta_reT(6,4) = 20./81.;
    delta_reT(6,4) = 20./27.;
    delta_reT(6,5) = 20./81.;
    delta_reT(6,6) = -10./27.;
    delta_reT(6,7) = -10./81.;
    delta_reT(6,8) = -10./27.;
    delta_reT(6,9) = -10./81.;
    delta_reT(8,2) = 20./27.;
    delta_reT(8,3) = 20./81.;
    delta_reT(8,4) = 20./27.;
    delta_reT(8,5) = 20./81.;
    delta_reT(8,6) = -10./27.;
    delta_reT(8,7) = -10./81.;
    delta_reT(8,8) = -10./27.;
    delta_reT(8,9) = -10./81.;
    
    }

    else {
    
    delta_reT(6,2) = -40./27.;
    delta_reT(6,3) = -40./81.;
    delta_reT(6,4) = -40./27.;
    delta_reT(6,5) = -40./81.;
    delta_reT(6,5) = -40./27.;
    delta_reT(6,6) = -40./81.;
    delta_reT(6,7) = -40./27.;
    delta_reT(6,8) = -40./81.;
    delta_reT(8,2) = -40./27.;
    delta_reT(8,3) = -40./81.;
    delta_reT(8,4) = -40./27.;
    delta_reT(8,5) = -40./81.;
    delta_reT(8,6) = -40./27.;
    delta_reT(8,7) = -40./81.;
    delta_reT(8,8) = -40./27.;
    delta_reT(8,9) = -40./81.;    
  
    }

    return(delta_reT);
    
}

gslpp::matrix<double>& EvolDF1nlep::Df1Evolnlep(double mu, double M, orders order, orders_ew order_ew, schemes scheme) 
{
    switch (scheme) {
        case NDR:
            break;
        case LRI:
        case HV:
        default:
            std::stringstream out;
            out << scheme;
            throw std::runtime_error("EvolDF1nlep::Df1Evolnlep_EM(): scheme " + out.str()
                    + " not implemented ");
    }
/* IMPORTANT!!: Please check cache for variation in AlsMZ and Ale. Ayan Paul*/
    if (mu == this->mu && M == this->M && scheme == this->scheme && order_ew == NULL_ew)
       return (*Evol(order));
    
    if (mu == this->mu && M == this->M && scheme == this->scheme &&  order_ew == NLO_ew)
       return (*Evol(order_ew));
        
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
        Df1Evolnlep(m_down, m_up, nf, scheme);
        Df1threshold_nlep(m_up, nf+1.);
        m_down = m_up;
        m_up = model.AboveTh(m_down);
        nf += 1.;
    }
    
    Df1Evolnlep(m_down, M, nf, scheme);
    
    if(order_ew != NULL_ew){
    
    return (*Evol(order_ew));
    }
    
    else { 
    return (*Evol(order)); 
    }
   
   }

void EvolDF1nlep::Df1Evolnlep(double mu, double M, double nf, schemes scheme) 
{

  gslpp::matrix<double> resLO(dim, 0.), resNLO(dim, 0.), resLO_ew(dim,0.), resNLO_ew(dim,0.);

    int L = 6 - (int) nf;
    double alsM = model.Als(M) / 4. / M_PI;
    double alsmu = model.Als(mu) / 4. / M_PI;
    double ale = model.getAle()/ 4. / M_PI ;
    
    double eta = alsM / alsmu;
    
    for (unsigned int k = 0; k < dim; k++) {
        double etap = pow(eta, a[L][k]);
        for (unsigned int i = 0; i < dim; i++) {
            for (unsigned int j = 0; j < dim; j++) {
                
                resLO(i,j) += b[L][i][j][k] * etap;
                
                resNLO(i,j) += c[L][i][j][k] * etap * alsmu;
                resNLO(i,j) += d[L][i][j][k] * etap * alsM;
                
                resLO_ew(i,j) +=  m[L][i][j][k] * etap * ale/alsmu;
                resLO_ew(i,j) +=  n[L][i][j][k] * etap * ale/alsM;
                
                resNLO_ew(i,j) += o[L][i][j][k] * etap * ale;
                resNLO_ew(i,j) += p[L][i][j][k] * etap * ale;
                resNLO_ew(i,j) += u[L][i][j][k] * etap * ale * log(eta);
                
                resNLO_ew(i,j) += q[L][i][j][k] * etap * ale;
                resNLO_ew(i,j) += r[L][i][j][k] * etap * ale;
                resNLO_ew(i,j) += s[L][i][j][k] * etap * ale / eta;
                resNLO_ew(i,j) += t[L][i][j][k] * etap * ale * eta;
            }
         }   
     }
 
    switch(order_ew) {
        case NLO_ew:
            *elem[NLO_ew] = (*elem[NLO]) * resLO_ew +
                            (*elem[NLO_ew]) * resLO + (*elem[LO]) *resNLO_ew;
        case LO_ew:
            *elem[LO_ew] =  (*elem[LO]) * resLO_ew;
            break;
        default:
            throw std::runtime_error("Error in EvolDF1nlep::Df1Evolnlep()");
    }   
   
    switch(order) {
        case NNLO:
            *elem[NNLO] = 0.;
        case NLO:
            *elem[NLO] = (*elem[LO]) * resNLO + (*elem[NLO]) * resLO;
        case LO:
            *elem[LO] = (*elem[LO]) * resLO;
            break;
        default:
            throw std::runtime_error("Error in EvolDF1nlep::Df1Evolnlep()");
    } 
}

void EvolDF1nlep::Df1threshold_nlep(double M, double nf){
 
    gslpp::matrix<double> drsT(dim,0.), dreT(dim,0.);
    
    double alsM = model.Als(M) / 4. / M_PI;
    double ale = model.getAle() / 4. / M_PI ;
    
    drsT = alsM * Df1threshold_deltarsT(nf);
    dreT = ale * Df1threshold_deltareT(nf);
    
     switch(order_ew){
         case NLO_ew:
             *elem[NLO_ew] += (*elem[LO])*dreT + (*elem[LO_ew]) * drsT ; 
             break;
         default:
             throw std::runtime_error("Error in EvolDF1nlep::Df1threshold_nlep()");
     }
    
    switch(order){
        case NNLO:
            *elem[NNLO] = 0.;
        case NLO:
            *elem[NLO] += (*elem[LO])* drsT ; 
            break;
        default:
            throw std::runtime_error("Error in EvolDF1nlep::Df1threshold_nlep()");
    }
    
   
}

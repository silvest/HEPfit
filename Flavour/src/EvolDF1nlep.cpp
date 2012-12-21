/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EvolDF1nlep.h"
#include <stdexcept>

EvolDF1nlep::EvolDF1nlep(unsigned int dim, schemes scheme, orders order, orders_ew
                         order_ew, const StandardModel& model) : model(model),
                         v(10,0.), vi(10,0.), js(10,0.), h(10,0.), gg(10,0.),
                         jv(10,0.), vij(10,0.), g_0(10,0.), k_0(10,0.), vk_0(10,0.),
                         k_0vi(10,0.), g_1(10,0.), k_11(10,0.),vk_11(10,0.), 
                         k_11vi(10,0.), k_12(10,0.), vh(10,0.),k_12vi(10,0.), 
                         k_13(10,0.), vk_13(10,0.), vg_1(10,0.), vg_0h(10,0.),
                         k12s(10,0.), s_svi(10,0.), k_12s(10,0.), k_12svi(10,0.),
                         k_13s(10,0.), vk_13s(10,0.), vs_s(10,0.), hvi(10,0.), 
                         vhg_0(10,0.), vg_0(10,0.), Gamma_T(10,0.), Gamma_1(10,0.),
                         Gamma_ew(10,0.), s_s(10,0.), jss(10,0.), jssv(10,0.), 
                         e(10,0.), RGEvolutor(dim, scheme, order, order_ew){
    
    double b0 = 0., b1 = 0.;
    
    for(int L=2; L>-1; L--){
    
    /* L=2 --> u,d,s,c (nf=4)  L=1 --> u,d,s,c,b (nf=5) L=0 --> u,d,s,c,b,t (nf=6)*/
     
    nu = L;  nd = L;
     
    b0 = model.Beta0(6-L);
    b1 = model.Beta1(6-L);
     
    if(L == 1){nd = 3; nu = 2;} 
    if(L == 0){nd = 3; nu = 3;}
    
   AnomalousDimension_nlep_S(LO,nu,nd).transpose().eigensystem(v,e);
   vi = v.inverse();
   for(int i = 0; i < 10; i++){
       a[L][i] = e(i).real();
       for (int j = 0; j < 10; j++){
            for (int k = 0; k < 10; k++) {
                b[L][i][j][k] = v(i, k).real() * vi(k, j).real();
            }
        }
    }
   
   /*e.m. part of the evolutor at the 1/ alpha_s order in alpha_s*/
    g_0 = vi * AnomalousDimension_nlep_EM(LO,nu,nd).transpose() * v;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++)  {
            if(fabs(e(i).real() - e(j).real() - 2. * b0)>0.000000000001){    
            k_0.assign( i , j , g_0(i,j)/( e(i).real() - e(j).real() - 2. * b0));
            } 
        }
    }
    vk_0 = v * k_0;
    k_0vi = k_0 * vi;    
    vg_0 = v * g_0;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++)  {
            if(fabs(e(i).real() - e(j).real() - 2. * b0)>0.000000000001){        
                    for (int k = 0; k < 10; k++) {
                        c[L][i][j][k] = - vk_0(i, k).real() * vi(k, j).real();
                        d[L][i][j][k] =  v(i, k).real() * k_0vi(k, j).real();
                }
            }
        
    else{
        for (int k = 0; k < 10; k++) {
                c[L][i][j][k] = (- 1. /(b0)) * vg_0(i, k).real() * vi(k, j).real();
                d[L][i][j][k] = 0.;
                }
            }
        }  
    }
    
      
    /*strong part of the evolutor at the NLO order*/
    gg = vi * AnomalousDimension_nlep_S(NLO,nu,nd).transpose() * v;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++)  {
            s_s.assign( i, j, (b1 / b0) * (i==j) * e(i).real() - gg(i,j));    
            if(fabs(e(i).real() - e(j).real() + 2. * b0)>0.000000000001){
                h.assign( i, j, s_s(i,j) / (2. * b0 + e(i) - e(j)));
            }
        }
    }
    js = v * h * vi;
    jv = js * v;
    vij = vi * js;
    jss = v * s_s * vi;
    jssv = jss * v;       
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++)   {
            if(fabs(e(i).real() - e(j).real() + 2. * b0)>0.000000000001)  {
                for (int k = 0; k < 10; k++)     {
                        m[L][i][j][k] = jv(i, k).real() * vi(k, j).real();
                        n[L][i][j][k] = -v(i, k).real() * vij(k, j).real();
                        }
                    }
            else{    
                for (int k = 0; k < 10; k++)     {
                        m[L][i][j][k] = (1./(4. * b0)) * jssv(i, k).real() * vi(k, j).real();
                        n[L][i][j][k] = 0.;
                        }   
                    }
                }
            }
     
    /*e.m. part of the evolutor at the NLO order*/
    Gamma_T = AnomalousDimension_nlep_EM(LO,nu,nd).transpose();
    Gamma_1 = AnomalousDimension_nlep_EM(NLO,nu,nd).transpose() -
             (b1/b0) * AnomalousDimension_nlep_EM(LO,nu,nd).transpose() ;
    Gamma_ew =  Gamma_1 + Gamma_T * js - js * Gamma_T;
    g_1 = vi * Gamma_ew * v;  
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++)   {
            if(fabs(e(i).real() - e(j).real()) > 0.000000000001){
                k_11.assign( i, j, g_1(i,j)/(e(i).real()-e(j).real()));
                }
            }
        }
    vk_11 = v * k_11;
    k_11vi = k_11 * vi; 
    k_12 = - k_0 * h; 
    k_12s = - k_0 * s_s;
    hvi = - h * vi;    
    s_svi = -s_s * vi;
    k_12vi = k_12 * vi;
    k_12svi = k_12s * vi;
    k_13 = h * k_0;
    k_13s = k_13 * s_s;
    vk_13 = v * k_13;
    vk_13s = v * k_13s;
    vh = v * h ;
    vs_s = v * s_s;
    vg_1 = v * g_1;
    vg_0h = v * g_0 * h;
    vhg_0 = v * h * g_0; 
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++)   {
            if(fabs(e(i).real() - e(j).real()) > 0.000000000001) {
                for (int k = 0; k < 10; k++)  {
                    o[L][i][j][k] = - vk_11(i,k).real() * vi(k,j).real();
                    p[L][i][j][k] = v(i,k).real() * k_11vi(k,j).real();
                }
            }   
            else{
                for (int k = 0; k < 10; k++)  {
                o[L][i][j][k] = - (1. / (2. * b0)) * vg_1(i,k).real() * vi(k,j).real();
                p[L][i][j][k] = 0.;
                }
        }
        if(fabs(e(i).real() - e(j).real() - 2. * b0) > 0.000000000001){
            if(fabs(e(i).real() - e(j).real() + 2. * b0) > 0.000000000001){
                for (int k = 0; k < 10; k++)  { 
                    q[L][i][j][k] = - vk_0(i, k).real() * hvi(k, j).real();
                    r[L][i][j][k] = v(i, k).real() * k_12vi(k, j).real();
                    s[L][i][j][k] = - vk_13(i, k).real() * vi(k, j).real();
                    t[L][i][j][k] = vh(i, k).real() * k_0vi(k, j).real();    
              }
            }
            else{
                for (int k = 0; k < 10; k++)  { 
                    q[L][i][j][k] = - (1./(2. * b0)) * vk_0(i, k).real() * s_svi(k, j).real();
                    r[L][i][j][k] = (1./(2. * b0)) * v(i, k).real() * k_12svi(k, j).real();
                    s[L][i][j][k] = - (1./(2. * b0)) * vk_13s(i, k).real() * vi(k, j).real();
                    t[L][i][j][k] = (1./(2. * b0)) * vs_s(i, k).real() * k_0vi(k, j).real();    
              }   
            }  
        }
        else{
            for (int k = 0; k < 10; k++)  {
                q[L][i][j][k] = (1. / (2. * b0)) * vg_0h(i,k).real() * vi(k,j).real();
                r[L][i][j][k] = 0.;
                s[L][i][j][k] = - (1. / (2. * b0)) * vhg_0(i,k).real() * vi(k,j).real();
                t[L][i][j][k] = 0.;
                    }
                }
            }
        } 
    }
}

EvolDF1nlep::~EvolDF1nlep() {
    
}

matrix<double> EvolDF1nlep::AnomalousDimension_nlep_S(orders order, unsigned int n_u,
        unsigned int n_d) const{
   
    /* anomalous dimension related to Delta F = 1 operators in Buras basis, hep-ph/9512380v1 */
    
    /*gamma(row, column) leading order*/
    
    unsigned int nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    matrix<double> gammaDF1(10, 0.);
   
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
                throw std::runtime_error("EvolDF1nlep::AnomalousDimension_B(): wrong number of flavours"); 
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
    
    gammaDF1(3,2) = 379./18.-56./243.*nf;
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
            throw std::runtime_error("EvolDF1nlep::AnomalousDimensio_B_S(): order not implemented"); 
    }
   
    return (gammaDF1);
    
  }

matrix<double> EvolDF1nlep::AnomalousDimension_nlep_EM(orders order, unsigned int n_u,
        unsigned int n_d) const{
   
    /* anomalous dimension related to Buras operators hep-ph/9512380v1 */
    /*gamma(riga, colonna) leading order*/
    unsigned int nf = n_u + n_d; /*n_u\d = active type up/down flavor d.o.f.*/
    matrix<double> gammaDF1(10, 0.);
   
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
                throw std::runtime_error("EvolDF1nlep::AnomalousDimension_EM(): wrong number of flavours"); 
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
    gammaDF1(4,4) = -2-136./243.*(n_u-n_d/2.);
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
    gammaDF1(7,6) = 2-5212./729.*n_u-2416./729.*n_d;
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
            throw std::runtime_error("EvolDF1nlep::AnomalousDimension_B_EM(): order not implemented"); 
    }
   
    return (gammaDF1);
    
}


matrix<double> EvolDF1nlep::Df1threshold_deltarsT(double nf) const {
    
    matrix <double> delta_rsT(10,0.); 
  
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

matrix<double> EvolDF1nlep::Df1threshold_deltareT(double nf) const {
 
    matrix<double> delta_reT(10,0.);    
    
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

matrix<double>& EvolDF1nlep::Df1Evolnlep(double mu, double M, orders order, orders_ew order_ew, schemes scheme) {
    switch (scheme) {
        case NDR:
            break;
        case LRI:
        case HV:
        default:
            std::stringstream out;
            out << scheme;
            throw std::runtime_error("EvolDF1nlep::Df1Evolnlep_EM(): scheme " + out.str() + " not implemented "); 
    }

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
        Df1threshold_nlep(m_up, nf);
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

void EvolDF1nlep::Df1Evolnlep(double mu, double M, double nf, schemes scheme) {

  matrix<double> resLO(10, 0.), resNLO(10, 0.), resLO_ew(10,0.), resNLO_ew(10,0.),
                 resNNLO(10, 0.);

    int L = 6 - (int) nf;
    double alsM = model.Als(M) / 4. / M_PI;
    double alsmu = model.Als(mu) / 4. / M_PI;
    double ale = model.getAle()/ 4. / M_PI ;
    double B_0 = model.Beta0(nf);
    
    double eta = alsM / alsmu;
    
    for (int k = 0; k < 10; k++) {
        double etap = pow(eta, a[L][k] / 2. / B_0);
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
               
               resNNLO(i, j) += 0.;
               if(fabs(e(i).real() - e(j).real() + 2. * B_0)>0.000000000001)  {
               resNLO(i, j) += m[L][i][j][k] * etap * alsmu;
               resNLO(i, j) += n[L][i][j][k] * etap * alsM;
               }
               else{
               resNLO(i, j) += - m[L][i][j][k] * etap * alsmu * log(eta);    
               }
               if(fabs(e(i).real() - e(j).real())>0.000000000001){
                 resNLO_ew(i, j) += o[L][i][j][k] * etap * ale;
                 resNLO_ew(i, j) += p[L][i][j][k] * etap * ale;
               } 
               else{ 
                 resNLO_ew(i, j) += - o[L][i][j][k] * etap * ale * log(eta); 
               }    
               
               if(fabs(e(i).real() - e(j).real() - 2. * B_0)>0.000000000001){
                 resLO_ew(i,j) += c[L][i][j][k] * etap * (1./alsmu) * ale;
                 resLO_ew(i,j) += d[L][i][j][k] * etap * (1./alsM) * ale;
                 if(fabs(e(i).real() - e(j).real() + 2. * B_0) > 0.000000000001){
                 resNLO_ew(i, j) += q[L][i][j][k] * etap *(alsM/alsmu) * ale;
                 resNLO_ew(i, j) += r[L][i][j][k] * etap * ale;
                 resNLO_ew(i, j) += s[L][i][j][k] * etap * ale;
                 resNLO_ew(i, j) += t[L][i][j][k] * etap * (alsmu/alsM) * ale;
                 }
                 else{
                 resNLO_ew(i, j) += - q[L][i][j][k] * etap *(alsM/alsmu) * ale * log(eta);
                 resNLO_ew(i, j) += - r[L][i][j][k] * etap * ale * log(eta);
                 resNLO_ew(i, j) += - s[L][i][j][k] * etap * ale * log(eta);
                 resNLO_ew(i, j) += - t[L][i][j][k] * etap * (alsmu/alsM) * ale * log(eta);    
                 }
               }
               else{
                 resLO_ew(i,j) += - c[L][i][j][k] * etap * (1./(2. * alsmu)) * ale * log(eta);  
                 resNLO_ew(i, j) += - q[L][i][j][k] * etap * (alsM/alsmu)* ale *log(eta) ; 
                 resNLO_ew(i, j) += - s[L][i][j][k] * etap *  ale * log(eta); 
               }          
                          
               resLO(i, j) += b[L][i][j][k] * etap;
              
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
       case NULL_ew:
       default:
           throw std::runtime_error("Error in EvolDF1nlep::Df1Evolnlep()");
   }   
    
    switch(order) {
   
        case NNLO:
            *elem[NNLO] = 0.;
            break;
        case NLO:
            *elem[NLO] = (*elem[LO]) * resNLO + (*elem[NLO]) * resLO;
        case LO:
            *elem[LO] = (*elem[LO]) * resLO;
            break;
        case FULLNNLO:
        case FULLNLO:
        default:
            throw std::runtime_error("Error in EvolDF1nlep::Df1Evolnlep()");

    } 
    
 
    }

void EvolDF1nlep::Df1threshold_nlep(double M, double nf){
 
    matrix<double> drsT(10,0.), dreT(10,0.);
    
    double alsM = model.Als(M) / 4. / M_PI;
    double ale = model.getAle() / 4. / M_PI ;
    
    drsT = alsM * Df1threshold_deltarsT(nf);
    dreT = ale * Df1threshold_deltareT(nf);
    
     switch(order_ew){
         case NLO_ew:
             *elem[NLO_ew] += (*elem[LO])*dreT + (*elem[LO_ew]) * drsT ; 
             break;
         case NULL_ew:
         case LO_ew:
         default:
             throw std::runtime_error("Error in EvolDF1nlep::Df1threshold_nlep()");
     }
    
    switch(order){
        case NNLO:
            *elem[NNLO] = 0.;
            break;
        case NLO:
            *elem[NLO] += (*elem[LO])* drsT ; 
            break;
        case LO:
        case FULLNNLO:
        case FULLNLO:
        default:
            throw std::runtime_error("Error in EvolDF1nlep::Df1threshold_nlep()");
    }
    
   
}

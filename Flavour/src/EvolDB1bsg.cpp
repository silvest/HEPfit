/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <gsl/gsl_sf_zeta.h>
#include "EvolDB1bsg.h"
#include "StandardModel.h"

EvolDB1bsg::EvolDB1bsg(unsigned int dim_i, schemes scheme, orders order, const StandardModel& model) 
:           RGEvolutor(dim_i, scheme, order), model(model),
            v(dim_i,0.), vi(dim_i,0.), js(dim_i,0.), h(dim_i,0.), gg(dim_i,0.), s_s(dim_i,0.),
            jssv(dim_i,0.), jss(dim_i,0.), jv(dim_i,0.), vij(dim_i,0.), gg2(dim_i,0.),
            s_s2(dim_i,0.),h2(dim_i,0.),js2(dim_i,0.), j2v(dim_i,0.), vijj(dim_i,0.),
            vij2(dim_i,0.), e(dim_i,0.), dim(dim_i) 
{
    if (dim != 8 ) throw std::runtime_error("ERROR: EvolDB1bsg can only be of dimension 8"); 
    
    /* magic numbers a & b */ 
    
    for(int L=2; L>-1; L--){
        
    /* L=2 --> u,d,s,c (nf=4)  L=1 --> u,d,s,c,b (nf=5) L=0 --> u,d,s,c,b,t (nf=6) */
        
    nu = L;  nd = L;
    if(L == 1){nd = 3; nu = 2;} 
    if(L == 0){nd = 3; nu = 3;}
    
    // LO evolutor of the effective Wilson coefficients in the Chetyrkin, Misiak and Munz basis
    
    (ToEffectiveBasis(ToRescaleBasis(LO,nu,nd))).transpose().eigensystem(v,e);
    vi = v.inverse();
    for(unsigned int i = 0; i < dim; i++){
       a[L][i] = e(i).real();
       for (unsigned int j = 0; j < dim; j++) {
           for (unsigned int k = 0; k < dim; k++)  {
                b[L][i][j][k] = v(i, k).real() * vi(k, j).real();
               }
           }
       }
    
    // NLO evolutor of the effective Wilson coefficients in the Chetyrkin, Misiak and Munz basis
    
    gg = vi * (ToEffectiveBasis(ToRescaleBasis(NLO,nu,nd))).transpose() * v;
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
    gg2 = vi * ( AnomalousDimension_M(NNLO,nu,nd).transpose() ) * v ;
    double b2 = model.Beta2(nu+nd);
    for (unsigned int i = 0; i < dim; i++){
        for (unsigned int j = 0; j < dim; j++){
            gslpp::complex s_s2_temp = 0.;
            for (unsigned int k = 0; k < dim; k++){
                s_s2_temp += (2. * b0 + e(i).real() - e(k).real()) * 
                ( h(i,k) * h(k,j) - b1/b0 * h(i,j) * (j==k) );
            }
            s_s2.assign( i, j, s_s2_temp + b2/b0 * (i==j) * e(i).real() - gg2(i,j));
            if(fabs(e(i).real() - e(j).real() + 4. * b0)>0.00000000001){
                h2.assign( i, j, s_s2(i,j) / (4. * b0 + e(i) - e(j)));
                }
            else{
                throw std::runtime_error("EvolDF1bsg::EvolDF1bsg(): singular case at NNLO not yet implemented!");
                }
            }
        }
    js2 = v * h2 * vi;
    j2v = js2 * v;
    vijj = vij * js;
    vij2 = vi * js2; 
    for (unsigned int i = 0; i < dim; i++){ 
        for (unsigned int j = 0; j < dim; j++){  
            for(unsigned int k = 0; k < dim; k++){ 
                c2[L][i][j][k] = j2v(i, k).real() * vi(k, j).real(); 
                d2[L][i][j][k] = -jv(i, k).real() * vij(k, j).real();  
                e2[L][i][j][k] = +v(i, k).real() * vijj(k, j).real(); 
                f2[L][i][j][k] = - v(i, k).real() * vij2(k, j).real();                
                }
            }
        }
    }
}
    
EvolDB1bsg::~EvolDB1bsg() 
{}

gslpp::matrix<double> EvolDB1bsg::AnomalousDimension_M(orders order, unsigned int n_u, unsigned int n_d) const
{
    
    /* Delta F = 1 anomalous dimension in Misiak basis, 
       ref.: M. Misiak, Nucl. Phys. B393 (1993) 23, B439 (1995) 461 (E),  
             A.J. Buras and M. Munz, Phys. Rev. D52 (1995) 186. */   
    
    /* gamma(row, coloumn) at the LO */
    
    unsigned int nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    double z3 = gsl_sf_zeta_int(3);
  
    gslpp::matrix<double> gammaDF1(dim, dim, 0.);
   
    switch(order){
        
    case LO:
        
    gammaDF1(0,0) = -4. ;
    gammaDF1(0,1) = 8./3. ;
    gammaDF1(0,3) = -2./9.;
    
    gammaDF1(1,0) = 12.;
    gammaDF1(1,3) = 4./3.;
    
    gammaDF1(2,3) = -52./3.;
    gammaDF1(2,5) = 2.;
    
    gammaDF1(3,2) = -40./9.;
    gammaDF1(3,3) = -160./9. + 4./3.*nf;
    gammaDF1(3,4) = 4./9.;
    gammaDF1(3,5) = 5./6.;
    
    gammaDF1(4,3) = -256./3.;
    gammaDF1(4,5) = 20.;
    
    gammaDF1(5,2) = -256./9.;
    gammaDF1(5,3) = -544./9. + (40./3.)*nf;
    gammaDF1(5,4) = 40./9.;
    gammaDF1(5,5) = -2./3.;
    
    gammaDF1(6,6) = 32./3. - 2.*model.Beta0(nf);   
    
    gammaDF1(7,6) = -32./9.;
    gammaDF1(7,7) = 28./3. - 2.*model.Beta0(nf);
    
    break;    
    case NLO:
        
    if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolDF1::AnomalousDimension_M(): wrong number of flavours"); 
    }
    
    /* gamma(row, coloumn) at the NLO */
    
    gammaDF1(0,0) = -145./3. + (16./9.)*nf;
    gammaDF1(0,1) = -26. + (40./27.)*nf;
    gammaDF1(0,2) = -1412./243.;
    gammaDF1(0,3) = -1369./243.;
    gammaDF1(0,4) = 134./243.;
    gammaDF1(0,5) = -35./162.;
    gammaDF1(0,6) = -232./243.;
    gammaDF1(0,7) = 167./162.;
            
    gammaDF1(1,0) = -45. + (20./3.)*nf; 
    gammaDF1(1,1) = -28./3.;
    gammaDF1(1,2) = -416./81.;
    gammaDF1(1,3) = 1280./81.;
    gammaDF1(1,4) = 56./81.;
    gammaDF1(1,5) = 35./27.;
    gammaDF1(1,6) = 464./81.;
    gammaDF1(1,7) = 76./27.;
    
    gammaDF1(2,2) = -4468./81.;
    gammaDF1(2,3) = -29129./81. - (52./9.)*nf;
    gammaDF1(2,4) = 400./81.;
    gammaDF1(2,5) = 3493./108. - (2./9.)*nf;
    gammaDF1(2,6) = 64./81.;
    gammaDF1(2,7) = 368./27.;
    
    gammaDF1(3,2) = -13678./243. + (368.*nf)/81.;       
    gammaDF1(3,3) = -79409./243. + (1334.*nf)/81.;
    gammaDF1(3,4) = 509./486. - (8.*nf)/81.;
    gammaDF1(3,5) = 13499./648. - (5.*nf)/27.;
    gammaDF1(3,6) = -680./243. + (32.*nf)/81;
    gammaDF1(3,7) = -427./81. - (37.*nf)/54.; 
            
    gammaDF1(4,2) = -244480./81. - (160./9.)*nf;
    gammaDF1(4,3) = -29648./81. - (2200./9.)*nf;
    gammaDF1(4,4) = 23116./81. + (16./9.)*nf;
    gammaDF1(4,5) = 3886./27. + (148./9.)*nf;
    gammaDF1(4,6) = -6464./81.;
    gammaDF1(4,7) = 8192./27. + 36.*nf;
            
    gammaDF1(5,2) = 77600./243. - (1264./81.)*nf;
    gammaDF1(5,3) = -28808./243. + (164./81.)*nf;
    gammaDF1(5,4) = -20324./243. + (400./81.)*nf;
    gammaDF1(5,5) = -21211./162.+ (622./27.)*nf;
    gammaDF1(5,6) = -20096./243. - (976.*n_d)/81. + (2912.*n_u)/81.;
    gammaDF1(5,7) = -6040./81. + (220./27.)*nf;
           
    gammaDF1(6,6) = 1936./9.-224./27.*nf-2*model.Beta1(nf); 
    
    gammaDF1(7,6) = -368./9.+224./81.*nf;
    gammaDF1(7,7) = 1456./9.-61./27.*nf-2*model.Beta1(nf); 

    break; 
    case NNLO:
        
    if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolDF1::AnomalousDimension_M(): wrong number of flavours"); 
    }
    
    // ADM -- already given in the effective and rescaled basis -- @ NNLO
    
    // See hep-ph/0411071 for the mixing among Q1-Q6 @ NNLO hep-ph/0504194
    
    gammaDF1(0,0) = -1927./2.+40./9.*(n_d*n_d+n_u*n_u)+224.*z3+(257./9.+160./3.*z3)*n_u+(257./9.+80./9.*n_u+160./3.*z3)*n_d;
    gammaDF1(0,1) = 475./9.-40./27.*n_d*n_d-40./27.*n_u*n_u+(362./27.-320./9.*z3)*n_u+(362./27.-80./27.*n_u-320./9.*z3)*n_d-896./3.*z3;
    gammaDF1(0,2) = 269107./13122.-2288./729.*nf-1360./81.*z3;
    gammaDF1(0,3) = -2425817./13122.+30815./4374.*nf-776./81.*z3;
    gammaDF1(0,4) = -343783./52488.+392./729.*nf+124./81.*z3;
    gammaDF1(0,5) = -37573./69984.+35./972.*nf+100./27.*z3;
            
    gammaDF1(1,0) = 307./2.-20./3.*(n_d*n_d+n_u*n_u)+(361./3.-160.*z3)*n_u+(361./3.-40./3.*n_u-160.*z3)*n_d-1344*z3;
    gammaDF1(1,1) = 1298./3.-76./3.*nf-224.*z3;
    gammaDF1(1,2) = 69797./2187.+904./243.*nf+2720./27.*z3;
    gammaDF1(1,3) = 1457549./8748.-22067./729.*nf-2768./27.*z3;
    gammaDF1(1,4) = -37889./8748.-28./243.*nf-248./27.*z3;
    gammaDF1(1,5) = 366919./11664.-35./162.*nf-110./9.*z3;
    
    gammaDF1(2,2) = -4203068./2187.+14012./243.*nf-608./27.*z3;
    gammaDF1(2,3) = -18422762./2187.+888605./2916. *nf+272./27.*nf*nf+(39824./27. + 160.*nf)*z3;
    gammaDF1(2,4) = 674281./4374.-1352./243.*nf-496./27.*z3;
    gammaDF1(2,5) = 9284531./11664.-26./27.*(n_d*n_d+n_u*n_u)+(-2798./81.-20.*z3)*n_u+(-2798./81.-52./27.*n_u-20.*z3)*n_d-1921./9.*z3;
    
    gammaDF1(3,2) = -5875184./6561.+217892./2187.*nf+472./81.*nf*nf+(27520./81.+1360./9.*nf)*z3;
    gammaDF1(3,3) = -70274587./13122.+8860733./17496.*nf-4010./729.*nf*nf+(16592./81.+2512./27.*nf)*z3;
    gammaDF1(3,4) = 2951809./52488.-52./81.*(n_d*n_d+n_u*n_u)+(-31175./8748.-136./9.*z3)*n_u+(-31175./8748.-104./81.*n_u-136./9.*z3)*n_d-3154./81.*z3;
    gammaDF1(3,5) = 3227801./8748.-65./54.*(n_d*n_d+n_u*n_u)+(-105293./11664.-220./9.*z3)*n_u+(-105293./11664.-65./27.*n_u-220./9.*z3)*n_d+200./27.*z3;
    
    gammaDF1(4,2) = -194951552./2187.+358672./81.*nf-2144./81.*nf*nf+87040./27.*z3;
    gammaDF1(4,3) =  -130500332./2187.+3088/27.*(n_d*n_d+n_u*n_u)+238016./27.*z3+(-2949616./729.+640.*z3)*n_u+(-2949616./729.+6176./27.*n_u+640.*z3)*n_d;
    gammaDF1(4,4) = 14732222./2187.+272./81.*(n_d*n_d+n_u*n_u)-27428./81.*n_u+(-27428./81.+544./81.*n_u)*n_d-13984./27.*z3;
    gammaDF1(4,5) = 16521659./2916.-316./27.*(n_d*n_d+n_u*n_u)+(8081./54.-200.*z3)*n_u+(8081./54.-632./27.*n_u-200.*z3)*n_d-22420./9.*z3;
    
    gammaDF1(5,2) = 162733912./6561.+17920./243.*(n_d*n_d+n_u*n_u)+174208./81.*z3+(-2535466./2187.+12160./9.*z3)*n_u+(-2535466./2187.+35840./243.*n_u+12160./9.*z3)*n_d;
    gammaDF1(5,3) = 13286236./6561.-159548./729.*(n_d*n_d+n_u*n_u)+(-1826023./4374.-9440./27.*z3)*n_u+(-1826023./4374.-319096./729.*n_u-9440./27.*z3)*n_d-24832./81.*z3;
    gammaDF1(5,4) = -22191107./13122.+395783./4374.*nf-1720./243.*nf*nf-(33832./81.+1360./9.*nf)*z3;
    gammaDF1(5,5) = -32043361./8748.+3353393./5832.*nf-533./81.*nf*nf+(9248./27.-1120./9.*nf)*z3;
            
    // See hep-ph/0612329 for the mixing of Q1-Q6 into Q7-Q8 @ NNLO
    
    gammaDF1(0,6) = 77506102./531441. - 875374./177147.*nf + 560./19683.*nf*nf - 9731./162.*(2./3.) + 11045./729.*nf*(2./3.) + 316./729.*nf*nf*(2./3.) + 3695./486.*(1./3.);
    gammaDF1(0,6) += z3*(-112216./6561. + 728./729.*nf + 25508./81.*(2./3.) - 64./81.*nf*(2./3.)-100./27.*(1./3.));
    gammaDF1(1,6) = -15463055./177147. + 242204./59049.*nf - 1120./6561.*nf*nf + 55748./27.*(2./3.) - 33970./243.*nf*(2./3.) - 632./243.*nf*nf*(2./3.) - 3695./81.*(1./3.);
    gammaDF1(1,6) += z3*(365696./2187. - 1168./243.*nf - 51232./27.*(2./3.) - 1024./27.*nf*(2./3.) + 200./9.*(1./3.));
    gammaDF1(2,6) =  102439553./177147. - 12273398./59049.*nf + 5824./6561.*nf*nf + 26639./81.*(1./3.) - 8./27.*nf*(1./3.) ;
    gammaDF1(2,6) += z3*(3508864./2187. - 1904./243.*nf - 1984./9.*(1./3.) -64./9.*nf*(1./3.));
    gammaDF1(3,6) = -2493414077./1062882. - 9901031./354294.*nf + 243872./59049.*nf*nf - 1184./6561.*nf*nf*nf - 49993./972.*(1./3.) + 305./27.*nf*(1./3.);
    gammaDF1(3,6) += z3*(-1922264./6561. + 308648./2187.*nf - 1280./243.*nf*nf + 1010./9.*(1./3.) - 200./27.*nf*(1./3.));
    gammaDF1(4,6) =  8808397748./177147. - 174839456./59049.*nf + 1600./729.*nf*nf - 669694./81.*(1./3.) + 10672./27.*nf*(1./3.);
    gammaDF1(4,6) += z3*(123543040./2187. - 207712./243.*nf + 128./27.*nf*nf - 24880./9.*(1./3.) - 640./9.*nf*(1./3.));
    gammaDF1(5,6) = 7684242746./531441. - 351775414./177147.*nf - 479776./59049.*nf*nf - 11456./6561.*nf*nf*nf + 3950201./243.*(1./3.) - 130538./81.*nf*(1./3.) - 592./81.*nf*nf*(1./3.);
    gammaDF1(5,6) += z3*(7699264./6561. + 2854976./2187.*nf - 12320./243.*nf*nf - 108584./9.*(1./3.) - 1136./27.*nf*(1./3.));
            
    gammaDF1(0,7) = -421272953./1417176. - 8210077./472392.*nf - 1955./6561.*nf*nf + z3*(-953042./2187. - 10381./486.*nf);
    gammaDF1(1,7) = 98548513./472392 - 5615165./78732.*nf - 2489./2187.*nf*nf + z3*(-607103./729. - 1679./81.*nf);
    gammaDF1(2,7) = 3205172129./472392. - 108963529./314928.*nf+58903./4374.*nf*nf + z3*(-1597588./729. + 13028./81.*nf - 20./9.*nf*nf);
    gammaDF1(3,7) = -6678822461./2834352. + 127999025./1889568.*nf + 1699073./157464.*nf*nf + 505./4374.*nf*nf*nf + z3*(2312684./2187. + 128347./729.*nf + 920./81.*nf*nf) ;
    gammaDF1(4,7) = 29013624461./118098. - 64260772./19683.*nf - 230962./243.*nf*nf - 148./27.*nf*nf*nf + z3*(-69359224./729. - 885356./81.*nf -5080./9.*nf*nf);
    gammaDF1(5,7) = -72810260309./708588. + 2545824851./472392.*nf - 33778271./78732.*nf*nf - 3988./2187.*nf*nf*nf + z3*(-61384768./2187. - 685472./729.*nf +350./81.*nf*nf);
    
    // See hep-ph/0504194 for the mixing among Q7-Q8 @ NNLO
    
    gammaDF1(6,6) = 307448./81.-23776./81.*nf-352./81.*nf*nf+(-1856./27.-1280./9.*nf)*z3;
    gammaDF1(7,6) = -164672./243.+17108./243.*nf+352./243.*nf*nf+(3776./81.+1280./27.*nf)*z3;
    gammaDF1(7,7) = 268807./81.-4343./27.*nf-461./81.*nf*nf+(-28624./27.-1312./9.*nf)*z3;
    
    break;
    default:
            throw std::runtime_error("EvolDF1bsg::AnomalousDimension_M(): order not implemented"); 
    }
    return (gammaDF1);
}
    
gslpp::matrix<double> EvolDB1bsg::ToRescaleBasis(orders order, unsigned int n_u, unsigned int n_d) const
{
    
    /* matrix entries for the anomalous dimension in the Chetyrkin, Misiak and Munz basis,
       ref. hep-ph/9711280v1, hep-ph/0306079 */ 

    gslpp::matrix<double> mat(dim, 0.);
    gslpp::matrix<double> mat1(dim, 0.);
    unsigned int nf = n_u + n_d;
            
    mat1(0,6) = - 13454./2187. + 44./2187.*nf;
    mat1(1,6) = 20644./729. - 88./729.*nf;
    mat1(2,6) = 119456./729. + 5440./729.*n_d -21776./729.*n_u;
    mat1(3,6) = - 202990./2187. + 32./729.*n_d*n_d + n_d*(16888./2187. + 64./729.*n_u) - 17132./2187.*n_u + 32./729.*n_u*n_u;
    mat1(4,6) = 530240./243. + 300928./729.*n_d - 461120./729.*n_u;
    mat1(5,6) = - 1112344./729. + 5432./729.*n_d*n_d + n_d*(419440./2187. - 2744./729.*n_u) + 143392./2187.*n_u - 8176./729.*n_u*n_u;
    
    mat1(0,7) = 25759./5832. + 431./5832.*nf;
    mat1(1,7) = 9733./486. - 917./972.*nf;
    mat1(2,7) = 82873./243. - 3361./243.*nf;
    mat1(3,7) = - 570773./2916. - 253./486.*n_d*n_d +n_d*(-40091./5832. - 253./243.*n_u) - 40091./5832.*n_u - 253./486.*n_u*n_u;
    mat1(4,7) = 838684./81. - 14.*n_d*n_d + n_d*(129074./243. - 28.*n_u) + 129074./243.*n_u - 14.*n_u*n_u;
    mat1(5,7) = - 923522./243. - 6031./486.*n_d*n_d + n_d*(-13247./1458. - 6031./243.*n_u) - 13247./1458.*n_u - 6031./486.*n_u*n_u;
    
    
    switch(order){
        case(NLO): 
            mat = AnomalousDimension_M(NLO, n_u, n_d);
            for (int i=0; i<6; i++){
                for (unsigned int j=6; j<dim; j++){
                    mat(i,j) = mat1(i,j);
                }
            }
            for (unsigned int i=6; i<dim; i++){
                for (unsigned int j=6; j<dim; j++){
                    mat(i,j) = mat(i,j) + 2. * (i==j) * model.Beta1(nf);
                }
            }
            return (mat);
        case(LO):
            mat = AnomalousDimension_M(LO, n_u, n_d);
            for (int i=0; i<6; i++){
                for (unsigned int j=6; j<dim; j++){
                    mat(i,j) = AnomalousDimension_M(NLO, n_u, n_d)(i,j);
                }
            }
            for (unsigned int i=6; i<dim; i++){
                for (unsigned int j=6; j<dim; j++){
                    mat(i,j) = mat(i,j) + 2. * (i==j) * model.Beta0(nf);
                }
            }
            return (mat);
    default:
            throw std::runtime_error("change to rescaled operator basis: order not implemented"); 
    }
    
}

gslpp::matrix<double> EvolDB1bsg::ToEffectiveBasis(gslpp::matrix<double> mat) const
{
    
    gslpp::matrix<double> y(dim, 0.);
    
    y(0,0) = 1.;
    y(1,1) = 1.;
    y(2,2) = 1.;
    y(3,3) = 1.;
    y(4,4) = 1.;
    y(5,5) = 1.;
    y(6,6) = 1.;
    y(7,7) = 1.;
    
    y(6,2) = -1./3.;
    y(6,3) = -4./9.;
    y(6,4) = -20./3.;
    y(6,5) = -80./9.;
    
    y(7,2) = 1.;
    y(7,3) = -1./6.;
    y(7,4) = 20.;
    y(7,5) = -10./3.;
    
    return( (y.inverse()).transpose() * mat * y.transpose() );
    
}

gslpp::matrix<double>& EvolDB1bsg::Df1Evolbsg(double mu, double M, orders order, schemes scheme) 
{
    
    switch (scheme) {
        case NDR:
            break;
        case LRI:
        case HV:
        default:
            std::stringstream out;
            out << scheme;
            throw std::runtime_error("EvolDF1bsg::Df1Evolbsg(): scheme " + out.str() + " not implemented "); 
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
        double nf = 5;//model.Nf(m_down); b to s gamma is always 5 flavour. This erroneously makes the evolutor cross thresholds.

//        while (m_up < M) {
//            Df1Evolbsg(m_down, m_up, nf, scheme);
//            m_down = m_up;
//            m_up = model.AboveTh(m_down);
//            nf += 1.;
//        } This code is commented out since it is not necessary unless thresholds are crossed.
        Df1Evolbsg(m_down, M, nf, scheme);
    }
    
    return (*Evol(order));
    
}
    
 void EvolDB1bsg::Df1Evolbsg(double mu, double M, double nf, schemes scheme) 
 {

    gslpp::matrix<double> resLO(dim, 0.), resNLO(dim, 0.), resNNLO(dim, 0.);

    int L = 6 - (int) nf;
    double alsM = model.Alstilde5(M); 
    double alsmu = model.Alstilde5(mu);
    
    double eta = alsM / alsmu;
    
    for (unsigned int k = 0; k < dim; k++) {
        double etap = pow(eta, a[L][k] / 2. / model.Beta0(nf));
        for (unsigned int i = 0; i < dim; i++){
            for (unsigned int j = 0; j < dim; j++) {
                resNNLO(i, j) += c2[L][i][j][k] * etap * alsmu * alsmu;
                resNNLO(i, j) += d2[L][i][j][k] * etap * alsmu * alsM;
                resNNLO(i, j) += e2[L][i][j][k] * etap * alsM * alsM;
                resNNLO(i, j) += f2[L][i][j][k] * etap * alsM * alsM;
                if(fabs(e(i).real() - e(j).real() + 2. * model.Beta0(nf))>0.000000000001)  {
                    resNLO(i, j) += c[L][i][j][k] * etap * alsmu;
                    resNLO(i, j) += d[L][i][j][k] * etap * alsM;
                }
                else{
                    resNLO(i, j) += - c[L][i][j][k] * etap * alsmu * log(eta);    
                }        
                resLO(i, j) += b[L][i][j][k] * etap;
                if (fabs(resLO(i, j)) <= 1.e-12) {resLO(i, j) = 0.;}
                if (fabs(resNLO(i, j)) <= 1.e-12) {resNLO(i, j) = 0.;}
                if (fabs(resNNLO(i, j)) <= 1.e-12) {resNNLO(i, j) = 0.;}
            }
        }
    }
    
    switch(order) {
        case NNLO:
            *elem[NNLO] = (*elem[LO]) * resNNLO + (*elem[NLO]) * resNLO + (*elem[NNLO]) * resLO;
        case NLO:
            *elem[NLO] = (*elem[LO]) * resNLO + (*elem[NLO]) * resLO;
        case LO:
            *elem[LO] = (*elem[LO]) * resLO;
            break;
        case FULLNNLO:
        case FULLNLO:
        default:
            throw std::runtime_error("Error in EvolDF1bsg::Df1Evolbsg()");
    } 
    
 }
 

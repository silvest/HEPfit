/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <gsl/gsl_sf_zeta.h>
#include "EvolDB1Mll.h"

EvolDB1Mll::EvolDB1Mll(unsigned int dim_i, schemes scheme, orders order, const StandardModel& model) 
:           RGEvolutor(dim_i, scheme, order), model(model),
            v(dim_i,0.), vi(dim_i,0.), js(dim_i,0.), h(dim_i,0.), gg(dim_i,0.), s_s(dim_i,0.),
            jssv(dim_i,0.), jss(dim_i,0.), jv(dim_i,0.), vij(dim_i,0.), e(dim_i,0.), dim(dim_i) 
{
    if (dim != 13 ) throw std::runtime_error("ERROR: EvolDB1Mll can only be of dimension 13"); 
    
    /* magic numbers a & b */ 
    
    for(int L=3; L>-1; L--){
        
    /* L=3 --> u,d,s (nf=3) L=2 --> u,d,s,c (nf=4)  L=1 --> u,d,s,c,b (nf=5) L=0 --> u,d,s,c,b,t (nf=6) */
        
    nu = L;  nd = L;
    if(L == 3){nd = 2; nu = 1;} 
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
    double b0 = model.Beta0(nu+nd);
    double b1 = model.Beta1(nu+nd);
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
    
EvolDB1Mll::~EvolDB1Mll() 
{}

gslpp::matrix<double> EvolDB1Mll::AnomalousDimension_M(orders order, unsigned int n_u, unsigned int n_d) const
{
    
    /* Delta F = 1 anomalous dimension in Misiak basis, 
       ref.: M. Misiak, Nucl. Phys. B393 (1993) 23, B439 (1995) 461 (E),  
             A.J. Buras and M. Munz, Phys. Rev. D52 (1995) 186. */   
    
    /* gamma(row, coloumn) at the LO */
    
    unsigned int nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
  
    gslpp::matrix<double> gammaDF1(dim, dim, 0.);
   
    switch(order){
        
    case LO:
        
    gammaDF1(0,0) = -4. ;
    gammaDF1(0,1) = 8./3. ;
    gammaDF1(0,3) = -2./9.;
    gammaDF1(0,8) = -32./27.;
    
    
    gammaDF1(1,0) = 12.;
    gammaDF1(1,3) = 4./3.;
    gammaDF1(1,8) = -8./9.;
    
    gammaDF1(2,3) = -52./3.;
    gammaDF1(2,5) = 2.;
    gammaDF1(2,8) = 8./9. + (8.*n_d)/3. - (16.*n_u)/3.;
    
    gammaDF1(3,2) = -40./9.;
    gammaDF1(3,3) = -160./9. + 4./3.*nf;
    gammaDF1(3,4) = 4./9.;
    gammaDF1(3,5) = 5./6.;
    gammaDF1(3,8) = 32./27.;
    
    gammaDF1(4,3) = -256./3.;
    gammaDF1(4,5) = 20.;
    gammaDF1(4,8) = 128./9.+(80.*n_d)/3. - (160.*n_u)/3.;
    
    gammaDF1(5,2) = -256./9.;
    gammaDF1(5,3) = -544./9. + (40./3.)*nf;
    gammaDF1(5,4) = 40./9.;
    gammaDF1(5,5) = -2./3.;
    gammaDF1(5,8) = 512./27.;
    
    gammaDF1(6,6) = 32./3. - 2.*model.Beta0(nf);   
    
    gammaDF1(7,6) = -32./9.;
    gammaDF1(7,7) = 28./3. - 2.*model.Beta0(nf);
    
    gammaDF1(8,8) = -2.*model.Beta0(nf); 
    
    gammaDF1(9,9) = -2.*model.Beta0(nf); 
    
    gammaDF1(10,10)= -2.*model.Beta0(nf);
    
    gammaDF1(11,11)= -2.*model.Beta0(nf);
    
    gammaDF1(12,12)= -2.*model.Beta0(nf);
    
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
    gammaDF1(0,7) = +167./162.;
    gammaDF1(0,8) = -2272./729.;
            
    gammaDF1(1,0) = -45. + (20./3.)*nf; 
    gammaDF1(1,1) = -28./3.;
    gammaDF1(1,2) = -416./81.;
    gammaDF1(1,3) = 1280./81.;
    gammaDF1(1,4) = 56./81.;
    gammaDF1(1,5) = 35./27.;
    gammaDF1(1,6) = 464./81.;
    gammaDF1(1,7) = 76./27.;
    gammaDF1(1,8) = 1952./243.;
    
    gammaDF1(2,2) = -4468./81.;
    gammaDF1(2,3) = -29129./81. - (52./9.)*nf;
    gammaDF1(2,4) = 400./81.;
    gammaDF1(2,5) = 3493./108. - (2./9.)*nf;
    gammaDF1(2,6) = 64./81.;
    gammaDF1(2,7) = 368./27.;
    gammaDF1(2,8) = -4160./243. + (32.*n_d)/3. - (64.*n_u)/3.;
    
    gammaDF1(3,2) = -13678./243. + (368.*nf)/81.;       
    gammaDF1(3,3) = -79409./243. + (1334.*nf)/81.;
    gammaDF1(3,4) = 509./486. - (8.*nf)/81.;
    gammaDF1(3,5) = 13499./648. - (5.*nf)/27.;
    gammaDF1(3,6) = -680./243. + (32.*nf)/81;
    gammaDF1(3,7) = -427./81. - (37.*nf)/54.; 
    gammaDF1(3,8) = 784./729. - (2272.*n_d)/243. + (2912.*n_u)/243.;
            
    gammaDF1(4,2) = -244480./81. - (160./9.)*nf;
    gammaDF1(4,3) = -29648./81. - (2200./9.)*nf;
    gammaDF1(4,4) = 23116./81. + (16./9.)*nf;
    gammaDF1(4,5) = 3886./27. + (148./9.)*nf;
    gammaDF1(4,6) = -6464./81.;
    gammaDF1(4,7) = 8192./27. + 36.*nf;
    gammaDF1(4,8) = -58112./243. + (320.*n_d)/3. - (640.*n_u)/3.;
            
    gammaDF1(5,2) = 77600./243. - (1264./81.)*nf;
    gammaDF1(5,3) = -28808./243. + (164./81.)*nf;
    gammaDF1(5,4) = -20324./243. + (400./81.)*nf;
    gammaDF1(5,5) = -21211./162.+ (622./27.)*nf;
    gammaDF1(5,6) = -20096./243. - (976.*n_d)/81. + (2912.*n_u)/81.;
    gammaDF1(5,7) = -6040./81. + (220./27.)*nf;
    gammaDF1(5,8) = -22784./729. - (20704.*n_d)/243. + (28544.*n_u)/243.;
           
    gammaDF1(6,6) = 1936./9.-224./27.*nf-2*model.Beta1(nf); 
    
    gammaDF1(7,6) = -368./9.+224./81.*nf;
    gammaDF1(7,7) = 1456./9.-61./27.*nf-2*model.Beta1(nf);  
    
    gammaDF1(8,8) = -2.*model.Beta1(nf); 
    
    gammaDF1(9,9) = -2.*model.Beta1(nf); 
    
    gammaDF1(10,10)= -2.*model.Beta1(nf);
    
    gammaDF1(11,11)= -2.*model.Beta1(nf);
    
    gammaDF1(12,12)= -2.*model.Beta1(nf);
   
    break;
    default:
            throw std::runtime_error("EvolDF1bsg::AnomalousDimension_M(): order not implemented"); 
    }
    return (gammaDF1);
}

gslpp::matrix<double> EvolDB1Mll::ToRescaleBasis(orders order, unsigned int n_u, unsigned int n_d) const
{
    
    /* matrix entries for the anomalous dimension in the Chetyrkin, Misiak and Munz basis,
       ref. hep-ph/9711280v1, hep-ph/0504194 */
    
    gslpp::matrix<double> mat(dim, 0.);
    gslpp::matrix<double> mat1(dim, 0.);
    unsigned int nf = n_u + n_d;
    double z3 = gsl_sf_zeta_int(3);
    
    mat1(0,6) = - 13454./2187. + 44./2187.*nf;
    mat1(1,6) = 20644./729. - 88./729.*nf;
    mat1(2,6) = 119456./729. + 5440./729.*n_d -21776./729.*n_u;
    mat1(3,6) = - 202990./2187. + 32./729.*n_d*n_d + n_d*(16888./2187. + 64./729.*n_u)
                - 17132./2187.*n_u + 32./729.*n_u*n_u;
    mat1(4,6) = 530240./243. + 300928./729.*n_d - 461120./729.*n_u;
    mat1(5,6) = - 1112344./729. + 5432./729.*n_d*n_d + n_d*(419440./2187. - 
                2744./729.*n_u) + 143392./2187.*n_u - 8176./729.*n_u*n_u;
    
    mat1(0,7) = 25759./5832. + 431./5832.*nf;
    mat1(1,7) = 9733./486. - 917./972.*nf;
    mat1(2,7) = 82873./243. - 3361./243.*nf;
    mat1(3,7) = - 570773./2916. - 253./486.*n_d*n_d +n_d*(-40091./5832. - 
                253./243.*n_u) - 40091./5832.*n_u - 253./486.*n_u*n_u;
    mat1(4,7) = 838684./81. - 14.*n_d*n_d + n_d*(129074./243. - 28.*n_u) + 
                129074./243.*n_u - 14.*n_u*n_u;
    mat1(5,7) = - 923522./243. - 6031./486.*n_d*n_d + n_d*(-13247./1458. - 6031./243.*n_u)
                -13247./1458.*n_u - 6031./486.*n_u*n_u;
    
    mat1(0,8) = - 2357278./19683. + 14440./6561.*n_d + 144688./6561.*n_u + 6976./243.*z3;
    mat1(1,8) = - 200848./6561. - 23696./2187.*n_d + 30736./2187.*n_u - 3584./81.*z3;
    mat1(2,8) = - 1524104./6561. - 176./27.*n_d*n_d + 352./27.*n_u*n_u +
                n_d*(257564./2187. + 176./27.*n_u - 128./3.*z3) - 256./81.*z3 + 
                n_u*(-382984./2187. + 256./3.*z3); 
    mat1(3,8) = 1535926./19683. + 1984./2187.*n_d*n_d - 5792./2187.*n_u*n_u +
                n_d*(-256901./6561. - 3808./2187.*n_u - 2720./81.*z3) - 
                5056./243.*z3 + n_u*(34942./6561. + 1600./81.*z3);
    mat1(4,8) = - 31433600./6561. - 2912./27.*n_d*n_d + 5824./27.*n_u*n_u +
                n_d*(- 3786616./2187. + 2912./27.*n_u - 1280./3.*z3) -
                4096./81.*z3 + n_u*(7525520./2187. + 2560./3.*z3);
    mat1(5,8) = 48510784./19683. -51296./2187.*n_d*n_d + 54976./2187.*n_u*n_u +
                n_u*(-11231648./6561. - 22016./81.*z3) + n_d*(340984./6561. + 
                3680./2187.*n_u - 8192./81.*z3) - 80896./243.*z3;
     
    
    switch(order){
        case(NLO): 
            mat = AnomalousDimension_M(NLO, n_u, n_d);
            for (unsigned int i=0; i<6; i++){
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
            for (unsigned int i=0; i<6; i++){
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

gslpp::matrix<double> EvolDB1Mll::ToEffectiveBasis(gslpp::matrix<double> mat) const
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
    y(8,8) = 1.;
    y(9,9) = 1.;
    y(10,10) = 1.;
    y(11,11) = 1.;
    y(12,12) = 1.;
    
    y(6,2) = -1./3.;
    y(6,3) = -4./9.;
    y(6,4) = -20./3.;
    y(6,5) = -80./9.;
    
    y(7,2) = 1.;
    y(7,3) = -1./6.;
    y(7,4) = 20.;
    y(7,5) = -10./3.;
    
    y(8,2) = 4./3.;
    y(8,4) = 64./9.;
    y(8,5) = 64./27.; // Add terms proportional to Log(mb/mub))
    
    return( (y.inverse()).transpose() * mat * y.transpose() );
    
}

gslpp::matrix<double>& EvolDB1Mll::Df1EvolMll(double mu, double M, orders order, schemes scheme) 
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

    double m_down = mu;
    double m_up = model.AboveTh(m_down);
    double nf = model.Nf(m_down);
    
    while (m_up < M) {
        Df1EvolMll(m_down, m_up, nf, scheme);
        m_down = m_up;
        m_up = model.AboveTh(m_down);
        nf += 1.;
    }
    Df1EvolMll(m_down, M, nf, scheme);
    
    return (*Evol(order));
    
    }
    
 void EvolDB1Mll::Df1EvolMll(double mu, double M, double nf, schemes scheme) 
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
                if (fabs(resLO(i, j)) < 1.e-12) {resLO(i, j) = 0.;}
                if (fabs(resNLO(i, j)) < 1.e-12) {resNLO(i, j) = 0.;}
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
            throw std::runtime_error("Error in EvolDF1bsg::Df1Evolbsg()");
    } 
    
  }
 


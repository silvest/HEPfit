/* 
 * File:   EvolBsmm.cpp
 * Author: claudio
 * 
 * Created on 30 settembre 2015, 18.44
 */

#include "EvolBsmm.h"
#include <gsl/gsl_sf.h>
#include <stdlib.h>
#include <gsl/gsl_types.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_inline.h>
#include <gsl/gsl_check_range.h>
#include <gsl/gsl_block_double.h>

EvolBsmm::EvolBsmm(unsigned int dim_i, schemes scheme, orders order, orders_ew order_ew, const StandardModel & model)
:   RGEvolutor(dim_i, scheme, order, order_ew), model(model), V(dim_i,0.), Vi(dim_i,0.),
    AA(dim_i,0.), BB(dim_i,0.), CC(dim_i,0.), DD(dim_i,0.), EE(dim_i,0.), FF(dim_i,0.),
    RR(dim_i,0.), e(dim_i,0.), vavi(0,0.), vbvi(0,0.), vcvi(0,0.), vdvi(0,0.),
        vevi(0,0.), vfvi(0,0.), vrvi(0,0.),vaevi(0,0.), vbbvi(0,0.), vbdvi(0,0.), vbevi(0,0.),
        vdbvi(0,0.), vdevi(0,0.), veavi(0,0.), vebvi(0,0.), vedvi(0,0.), veevi(0,0.), vbeevi(0,0.), 
        vebevi(0,0.), veebvi(0,0.), vbbevi(0,0.), vbebvi(0,0.), vebbvi(0,0.), dim(dim_i)
{
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int l = 0;
    unsigned int m = 0;
    unsigned int p = 0;
    unsigned int q = 0;
    int L = 1;
    double  b0 = 0., b1 = 0.;

  
    vavi.clear();
    vbvi.clear();
    vcvi.clear();
    vdvi.clear();
    vfvi.clear();
    vevi.clear();
    vrvi.clear();
    vaevi.clear();
    vbbvi.clear(); 
    vbdvi.clear();
    vbevi.clear();
    vdbvi.clear();
    vdevi.clear();
    veavi.clear();
    vebvi.clear();
    vedvi.clear();
    veevi.clear();
    vbeevi.clear();
    vebevi.clear(); 
    veebvi.clear();
    vbbevi.clear(); 
    vbebvi.clear(); 
    vebbvi.clear();
    
    
 
/* define L, nu, nd */
    if(L == 1){nd = 3; nu = 2;} 
    b0 = model.Beta0(6-L);
    b1 = model.Beta1(6-L);
	    
    AnomalousDimension(10,nu,nd).transpose().eigensystem(V,e); 
   
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
	    


    AA = BuiltB('A', nu, nd);  
    BB = BuiltB('B', nu, nd);
    CC = BuiltB('C', nu, nd);
    DD = BuiltB('D', nu, nd);
    EE = BuiltB('E', nu, nd);
    FF = BuiltB('F', nu, nd);
    RR = BuiltB('R', nu, nd);
        
    double cutoff = 0.000000001 ;
    
    for(l = 0; l < dim; l++){
        for(i = 0; i < dim; i++){
            for(j = 0; j < dim; j++){
            	for(m = 0; m < dim; m++){
                  
                    if(fabs(V(l, i).real() * AA(i, j).real() * Vi(j, m).real()) > cutoff){
                    
                        vavi.push_back(l);
                        vavi.push_back(i);
                        vavi.push_back(j);
                        vavi.push_back(m);
                        vavi.push_back(V(l, i).real() * AA(i, j).real() * Vi(j, m).real());
                   
                    }
                    if(fabs(V(l, i).real() * BB(i, j).real() * Vi(j, m).real()) > cutoff){
                    
                        vbvi.push_back(l);
                        vbvi.push_back(i);
                        vbvi.push_back(j);
                        vbvi.push_back(m);
                        vbvi.push_back(V(l, i).real() * BB(i, j).real() * Vi(j, m).real());
                    
                    }
                    if(fabs(V(l, i).real() * CC(i, j).real() * Vi(j, m).real()) > cutoff){
                    
                        vcvi.push_back(l);
                        vcvi.push_back(i);
                        vcvi.push_back(j);
                        vcvi.push_back(m);
                        vcvi.push_back(V(l, i).real() * CC(i, j).real() * Vi(j, m).real());
                    
                    }
                    if(fabs(V(l, i).real() * DD(i, j).real() * Vi(j, m).real()) > cutoff){
                    
                        vdvi.push_back(l);
                        vdvi.push_back(i);
                        vdvi.push_back(j);
                        vdvi.push_back(m);
                        vdvi.push_back(V(l, i).real() * DD(i, j).real() * Vi(j, m).real());
                    
                    }
                    if(fabs(V(l, i).real() * EE(i, j).real() * Vi(j, m).real()) > cutoff){
                    
                        vevi.push_back(l);
                        vevi.push_back(i);
                        vevi.push_back(j);
                        vevi.push_back(m);
                        vevi.push_back(V(l, i).real() * EE(i, j).real() * Vi(j, m).real());
                    
                    }
                    if(fabs(V(l, i).real() * FF(i, j).real() * Vi(j, m).real()) > cutoff){
                    
                        vfvi.push_back(l);
                        vfvi.push_back(i);
                        vfvi.push_back(j);
                        vfvi.push_back(m);
                        vfvi.push_back(V(l, i).real() * FF(i, j).real() * Vi(j, m).real());

                    }
                    if(fabs(V(l, i).real() * RR(i, j).real() * Vi(j, m).real()) > cutoff){
                    
                        vrvi.push_back(l);
                        vrvi.push_back(i);
                        vrvi.push_back(j);
                        vrvi.push_back(m);
                        vrvi.push_back(V(l, i).real() * RR(i, j).real() * Vi(j, m).real());

                    }
                   
                    for(p = 0; p < dim; p++){

                        if(fabs(V(l, i).real() * AA(i, p).real() * EE(p, j).real() * Vi(j, m).real()) > cutoff){
                    
                            vaevi.push_back(l);
                            vaevi.push_back(i);
                            vaevi.push_back(p);
                            vaevi.push_back(j);
                            vaevi.push_back(m);
                            vaevi.push_back(V(l, i).real() * AA(i, p).real() * EE(p, j).real() * Vi(j, m).real());
                   
                        }
                        if(fabs(V(l, i).real() * BB(i, p).real() * BB(p, j).real() * Vi(j, m).real()) > cutoff){
                    
                            vbbvi.push_back(l);
                            vbbvi.push_back(i);
                            vbbvi.push_back(p);
                            vbbvi.push_back(j);
                            vbbvi.push_back(m);
                            vbbvi.push_back(V(l, i).real() * BB(i, p).real() * BB(p, j).real() * Vi(j, m).real());
                    
                        }
                        if(fabs(V(l, i).real() * BB(i, p).real() * DD(p, j).real() * Vi(j, m).real()) > cutoff){
                    
                            vbdvi.push_back(l);
                            vbdvi.push_back(i);
                            vbdvi.push_back(p);
                            vbdvi.push_back(j);
                            vbdvi.push_back(m);
                            vbdvi.push_back(V(l, i).real() * BB(i, p).real() * DD(p, j).real() * Vi(j, m).real());
                   
                        }
                        if(fabs(V(l, i).real() * BB(i, p).real() * EE(p, j).real() * Vi(j, m).real()) > cutoff){
                    
                            vbevi.push_back(l);
                            vbevi.push_back(i);
                            vbevi.push_back(p);
                            vbevi.push_back(j);
                            vbevi.push_back(m);
                            vbevi.push_back(V(l, i).real() * BB(i, p).real() * EE(p, j).real() * Vi(j, m).real());
                   
                        }
                        if(fabs(V(l, i).real() * DD(i, p).real() * BB(p, j).real() * Vi(j, m).real()) > cutoff){
                    
                            vdbvi.push_back(l);
                            vdbvi.push_back(i);
                            vdbvi.push_back(p);
                            vdbvi.push_back(j);
                            vdbvi.push_back(m);
                            vdbvi.push_back(V(l, i).real() * DD(i, p).real() * BB(p, j).real() * Vi(j, m).real());
                   
                        }
                        if(fabs(V(l, i).real() * DD(i, p).real() * EE(p, j).real() * Vi(j, m).real()) > cutoff){
                    
                            vdevi.push_back(l);
                            vdevi.push_back(i);
                            vdevi.push_back(p);
                            vdevi.push_back(j);
                            vdevi.push_back(m);
                            vdevi.push_back(V(l, i).real() * DD(i, p).real() * EE(p, j).real() * Vi(j, m).real());
                   
                        }
                        if(fabs(V(l, i).real() * EE(i, p).real() * AA(p, j).real() * Vi(j, m).real()) > cutoff){
                    
                            veavi.push_back(l);
                            veavi.push_back(i);
                            veavi.push_back(p);
                            veavi.push_back(j);
                            veavi.push_back(m);
                            veavi.push_back(V(l, i).real() * EE(i, p).real() * AA(p, j).real() * Vi(j, m).real());
                   
                        }
                        if(fabs(V(l, i).real() * EE(i, p).real() * BB(p, j).real() * Vi(j, m).real()) > cutoff){
                    
                            vebvi.push_back(l);
                            vebvi.push_back(i);
                            vebvi.push_back(p);
                            vebvi.push_back(j);
                            vebvi.push_back(m);
                            vebvi.push_back(V(l, i).real() * EE(i, p).real() * BB(p, j).real() * Vi(j, m).real());
                   
                        }    
                        if(fabs(V(l, i).real() * EE(i, p).real() * DD(p, j).real() * Vi(j, m).real()) > cutoff){
                    
                            vedvi.push_back(l);
                            vedvi.push_back(i);
                            vedvi.push_back(p);
                            vedvi.push_back(j);
                            vedvi.push_back(m);
                            vedvi.push_back(V(l, i).real() * EE(i, p).real() * DD(p, j).real() * Vi(j, m).real());
                   
                        }    
                        if(fabs(V(l, i).real() * EE(i, p).real() * EE(p, j).real() * Vi(j, m).real()) > cutoff){
                    
                            veevi.push_back(l);
                            veevi.push_back(i);
                            veevi.push_back(p);
                            veevi.push_back(j);
                            veevi.push_back(m);
                            veevi.push_back(V(l, i).real() * EE(i, p).real() * EE(p, j).real() * Vi(j, m).real());
                   
                        }
                        
                        for(q = 0; q < dim; q++){
			
                            if(fabs(V(l, i).real() * BB(i, p).real() * EE(p, q).real() * EE(q, j).real() * Vi(j, m).real()) > cutoff){
                    
                                vbeevi.push_back(l);
                                vbeevi.push_back(i);
                                vbeevi.push_back(p);
                                vbeevi.push_back(q);
                                vbeevi.push_back(j);
                                vbeevi.push_back(m);
                                vbeevi.push_back(V(l, i).real() * BB(i, p).real() * EE(p, q).real() * EE(q, j).real() * Vi(j, m).real());
                   
                            }
                            
                            if (fabs(V(l, i).real() * EE(i, p).real() * BB(p, q).real() * EE(q, j).real() * Vi(j, m).real()) > cutoff) {

                                vebevi.push_back(l);
                                vebevi.push_back(i);
                                vebevi.push_back(p);
                                vebevi.push_back(q);
                                vebevi.push_back(j);
                                vebevi.push_back(m);
                                vebevi.push_back(V(l, i).real() * EE(i, p).real() * BB(p, q).real() * EE(q, j).real() * Vi(j, m).real());

                            }
                            if (fabs(V(l, i).real() * EE(i, p).real() * EE(p, q).real() * BB(q, j).real() * Vi(j, m).real()) > cutoff) {

                                veebvi.push_back(l);
                                veebvi.push_back(i);
                                veebvi.push_back(p);
                                veebvi.push_back(q);
                                veebvi.push_back(j);
                                veebvi.push_back(m);
                                veebvi.push_back(V(l, i).real() * EE(i, p).real() * EE(p, q).real() * BB(q, j).real() * Vi(j, m).real());

                            }
                            if (fabs(V(l, i).real() * BB(i, p).real() * BB(p, q).real() * EE(q, j).real() * Vi(j, m).real()) > cutoff) {

                                vbbevi.push_back(l);
                                vbbevi.push_back(i);
                                vbbevi.push_back(p);
                                vbbevi.push_back(q);
                                vbbevi.push_back(j);
                                vbbevi.push_back(m);
                                vbbevi.push_back(V(l, i).real() * BB(i, p).real() * BB(p, q).real() * EE(q, j).real() * Vi(j, m).real());

                            }
                            if (fabs(V(l, i).real() * BB(i, p).real() * EE(p, q).real() * BB(q, j).real() * Vi(j, m).real()) > cutoff) {

                                vbebvi.push_back(l);
                                vbebvi.push_back(i);
                                vbebvi.push_back(p);
                                vbebvi.push_back(q);
                                vbebvi.push_back(j);
                                vbebvi.push_back(m);
                                vbebvi.push_back(V(l, i).real() * BB(i, p).real() * EE(p, q).real() * BB(q, j).real() * Vi(j, m).real());

                            }
                            if (fabs(V(l, i).real() * EE(i, p).real() * BB(p, q).real() * BB(q, j).real() * Vi(j, m).real()) > cutoff) {

                                vebbvi.push_back(l);
                                vebbvi.push_back(i);
                                vebbvi.push_back(p);
                                vebbvi.push_back(q);
                                vebbvi.push_back(j);
                                vebbvi.push_back(m);
                                vebbvi.push_back(V(l, i).real() * EE(i, p).real() * BB(p, q).real() * BB(q, j).real() * Vi(j, m).real());

                            }
			}
                    }
		}
            }
	}
    }
            
}


EvolBsmm::~EvolBsmm()
{}

gslpp::matrix<double> EvolBsmm::AnomalousDimension(int gam, unsigned int n_u, unsigned int n_d) const

{
	       
	   /* Delta F = 1 anomalous dimension in Misiak basis, 
	        ref.: ref. hep-ph/1311.1348v2, hep-ph/0512066 */   
	        
	/* gamma(row, coloumn) at the LO */
	    
    unsigned int nf = n_u + n_d; /*n_u/d = active type up/down flavor d.o.f.*/
    double zeta3 = 0;
        
    gslpp::matrix<double> gammaDF1(dim, 0.);
	
    zeta3 = gsl_sf_zeta_int(3);
        
    switch(gam){
	       
        case 10:
            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }

            gammaDF1(0,0) = -4. ;
            gammaDF1(0,1) = 8./3. ;
            gammaDF1(0,3) = -2./9.;

            gammaDF1(1,0) = 12.;
            gammaDF1(1,3) = 4./3.;

            gammaDF1(2,3) = -52./3.;
            gammaDF1(2,5) = 2.;

            gammaDF1(3,2) = -40./9.;
            gammaDF1(3,3) = -100./9.;
            gammaDF1(3,4) = 4./9.;
            gammaDF1(3,5) = 5./6.;

            gammaDF1(4,3) = -256./3.;
            gammaDF1(4,5) = 20.;

            gammaDF1(5,2) = -256./9.;
            gammaDF1(5,3) = 56./9.;
            gammaDF1(5,4) = 40./9.;
            gammaDF1(5,5) = -2./3.;   
	  
	break;    
	case 20:
	   
            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }

            /* gamma(row, coloumn) at the NLO */

            gammaDF1(0,0) = -355./9.;
            gammaDF1(0,1) = -502./27.;
            gammaDF1(0,2) = -1412./243.;
            gammaDF1(0,3) = -1369./243.;
            gammaDF1(0,4) = 134./243.;
            gammaDF1(0,5) = -35./162.;

            gammaDF1(1,0) = -35./3.; 
            gammaDF1(1,1) = -28./3.;
            gammaDF1(1,2) = -416./81.;
            gammaDF1(1,3) = 1280./81.;
            gammaDF1(1,4) = 56./81.;
            gammaDF1(1,5) = 35./27.;

            gammaDF1(2,2) = -4468./81.;
            gammaDF1(2,3) = -31469./81.;
            gammaDF1(2,4) = 400./81.;
            gammaDF1(2,5) = 3373./108.;

            gammaDF1(3,2) = -8158./243.;       
            gammaDF1(3,3) = -59399./243.;
            gammaDF1(3,4) = 269./486.;
            gammaDF1(3,5) = 12899./648.;

            gammaDF1(4,2) = -251680./81.;
            gammaDF1(4,3) = -128648./81.;
            gammaDF1(4,4) = 23836./81.;
            gammaDF1(4,5) = 6106./27.;

            gammaDF1(5,2) = 58640./243.;
            gammaDF1(5,3) = -26348./243.;
            gammaDF1(5,4) = -14324./243.;
            gammaDF1(5,5) = -2551./162.;
	 
        break;

	case 30:
   
            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }

            /* gamma(row, coloumn) at the NNLO */

            gammaDF1(0,0) = -12773./18. +  zeta3 * 1472./3.;
            gammaDF1(0,1) = 745./9. - zeta3 *4288./9.;
            gammaDF1(0,2) = 63187./13122. - zeta3 * 1360./81.;
            gammaDF1(0,3) = -981796./6561. - zeta3 * 776./81.;
            gammaDF1(0,4) = -202663./52488. + zeta3 * 124./81.;
            gammaDF1(0,5) = -24973./69984. + zeta3 * 100./27.;

            gammaDF1(1,0) = 1177./2. - zeta3 * 2144.; 
            gammaDF1(1,1) = 306. - zeta3 * 224.;
            gammaDF1(1,2) = 110477./2187. + zeta3 * 2720./27.;
            gammaDF1(1,3) = 133529./8748. - zeta3 * 2768./27.;
            gammaDF1(1,4) = -42929./8748. - zeta3 * 248./27.;
            gammaDF1(1,5) = 354319./11664. - zeta3 * 110./9.;

            gammaDF1(2,2) = -3572528./2187. - zeta3 * 608./27.;
            gammaDF1(2,3) = -58158773./8748. + zeta3 * 61424./27.;
            gammaDF1(2,4) = 552601./4374. - zeta3 * 496./27.;
            gammaDF1(2,5) = 6989171./11664. - zeta3 * 2821./9.;

            gammaDF1(3,2) = -1651004./6561. + zeta3 * 88720./81.;       
            gammaDF1(3,3) = -155405353./52488 + zeta3 * 54272./81.;
            gammaDF1(3,4) = 1174159./52488. - zeta3 * 9274./81.;
            gammaDF1(3,5) = 10278809./34992. - zeta3 * 3100./27.;

            gammaDF1(4,2) = -147978032./2187. + zeta3 * 87040./27.;
            gammaDF1(4,3) = -168491372./2187. + zeta3 * 324416./27.;
            gammaDF1(4,4) = 11213042./2187. - zeta3 * 13984./27.;
            gammaDF1(4,5) = 17850329./2916. - zeta3 * 31420./9.;

            gammaDF1(5,2) = 136797922./6561. + zeta3 * 721408./81.;
            gammaDF1(5,3) = -72614473./13122. - zeta3 * 166432./81.;
            gammaDF1(5,4) = -9288181./6561. - zeta3 * 95032./81.;
            gammaDF1(5,5) = -16664027./17496. - zeta3 * 7552./27.;
	
        break;	
	case 01:

            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }    

            gammaDF1(0,0) = -8./3.;
            gammaDF1(0,6) = -32./27.;

            gammaDF1(1,1) = -8./3.;
            gammaDF1(1,6) = -8./9.;

            gammaDF1(2,6) = -16./9.;

            gammaDF1(3,6) = 32./27.;

            gammaDF1(4,6) = -112./9.;

            gammaDF1(5,6) = 512./27.;  

            gammaDF1(6,6) = 8.;
            gammaDF1(6,7) = -4.; 

            gammaDF1(7,6) = -4.; 
  	
        break;  
	case 11:
   
            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }


            gammaDF1(0,0) = 169./9.;
            gammaDF1(0,1) = 100./27.;
            gammaDF1(0,3) = 254./729.;
            gammaDF1(0,6) = -2272./729.;

            gammaDF1(1,0) = 50./3.; 
            gammaDF1(1,1) = -8./3.;
            gammaDF1(1,3) = 1076./243.;
            gammaDF1(1,6) = 1952./243.;

            gammaDF1(2,3) = 11116./243.;
            gammaDF1(2,5) = -14./3.;
            gammaDF1(2,6) = -6752./243.;

            gammaDF1(3,2) = 280./27.;       
            gammaDF1(3,3) = 18763./729.;
            gammaDF1(3,4) = -28./27.;
            gammaDF1(3,5) = -35./18.;
            gammaDF1(3,6) = -2192./729.;  

            gammaDF1(4,3) = 111136./243.;
            gammaDF1(4,5) = -140./3.;
            gammaDF1(4,6) = -84032./243.;

            gammaDF1(5,2) = 2944./27.;
            gammaDF1(5,3) = 193312./729.;
            gammaDF1(5,4) = -280./27.;
            gammaDF1(5,5) = -175./9.;
            gammaDF1(5,6) = -37856./729.;

            gammaDF1(6,7) = 16.;

            gammaDF1(7,6) = 16.; 

        break;
	case 02:
   
            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                    throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }

             gammaDF1(0,6)= -11680./2187.;
             gammaDF1(0,7)= -416./81.;

             gammaDF1(1,6)= -2920./729.;
             gammaDF1(1,7)= -104./27.;

             gammaDF1(2,6)= -39752./729.;
             gammaDF1(2,7)= -136./27.;

             gammaDF1(3,6)= 1024./2187.;
             gammaDF1(3,7)= -448./81.;

             gammaDF1(4,6)= -381344./729.;
             gammaDF1(4,7)= -15616./27.;

             gammaDF1(5,6)= 24832./2187.;
             gammaDF1(5,7)= -7936./81.;


        break;
	case 21:
   
            if (!(nf == 3 || nf == 4 || nf == 5 || nf == 6)){ 
                throw std::runtime_error("EvolBsmm::AnomalousDimension(): wrong number of flavours"); 
            }

            gammaDF1(0,6)= -1359190./19683. + zeta3 * 6976./243.;

            gammaDF1(1,6)= -229696./6561 - zeta3 * 3584./81.;

            gammaDF1(2,6)= -1290092./6561 + zeta3 * 3200./81.;

            gammaDF1(3,6)= -819971./19683. - zeta3 * 19936./243;

            gammaDF1(4,6)= -16821944./6561 + zeta3 * 30464./81.;

            gammaDF1(5,6)= -17787368./19683. - zeta3 * 286720./243.;


        break;

	default:
            throw std::runtime_error("EvolBsmm::AnomalousDimension(): order not implemented"); 
    }
    return (gammaDF1);
}


gslpp::matrix<double> EvolBsmm::BuiltB(char letter, unsigned int n_u, unsigned int n_d)

{

    unsigned int nf = 5; //all the class works for nf = 5

    double B00S = model.Beta0(nf), B10S = model.Beta1(nf), B20S = model.Beta2(nf), 
            B01S = -22./9., B11S = -308./27.; 

    double B00E = 80./9., B01E = 176./9.;

    double b1 = B10S/(2. * B00S * B00S), b2 = B20S/(4. * B00S * B00S * B00S) - b1 * b1 , 
            b3 = B01S/(2. * B00S * B00E), b4 = B11S /(4. * B00S * B00S * B00E) - 2 * b1 * b3, 
            b5 = B01E/(2. * B00S * B00E) - b1;



    gslpp::matrix<double> B(dim, 0.);
    gslpp::matrix<double> W10T(dim, 0.);
    gslpp::matrix<double> W20T(dim, 0.);
    gslpp::matrix<double> W30T(dim, 0.);
    gslpp::matrix<double> W01T(dim, 0.);
    gslpp::matrix<double> W02T(dim, 0.);
    gslpp::matrix<double> W11T(dim, 0.);
    gslpp::matrix<double> W21T(dim, 0.);


    W10T = (AnomalousDimension(10, n_u, n_d).transpose())/2./B00S;
    W20T = (AnomalousDimension(20, n_u, n_d).transpose())/4./B00S/B00S;
    W30T = (AnomalousDimension(30, n_u, n_d).transpose())/8./B00S/B00S/B00S;
    W01T = (AnomalousDimension(01, n_u, n_d).transpose())/2./B00E;
    W02T = (AnomalousDimension(02, n_u, n_d).transpose())/4./B00E/B00E;
    W11T = (AnomalousDimension(11, n_u, n_d).transpose())/4./B00S/B00E;
    W21T = (AnomalousDimension(21, n_u, n_d).transpose())/8./B00S/B00S/B00E;



    switch(letter){

        case 'A':

            B = Vi.real() * (W30T - b1 * W20T - b2 * W10T) * V.real();

        break;    
        case 'B':

            B = Vi.real() * (W20T - b1 * W10T) * V.real();

        break;
        case 'C':

            B = Vi.real() * (W21T - b1 * W11T - b2 * W01T - b3 * W20T - b4 * W10T) * V.real();

        break;
        case 'D':

            B = Vi.real() * (W11T - b1 * W01T - b3 * W10T) * V.real();

        break;
        case 'E':

            B = Vi.real() * (W01T) * V.real();

        break;
        case 'F':

            B = Vi.real() * (W02T + W11T - (b1 + b3) * W01T - b3 * W10T) * V.real();

        break;
        case 'R':

            B = b5 * Vi.real() * (W01T) * V.real();

        break;
        default:
            throw std::runtime_error("EvolBsmm::BuiltB(): order not implemented"); 
    }
    return (B);
}



gslpp::matrix< double > & EvolBsmm::Df1Evol(double mu, double M, orders order, orders_ew order_ew, schemes scheme) 

{
    switch (scheme) {             /*  complete this method */
        case NDR:
        break;
        case LRI:
        case HV:
        default:
            std::stringstream out;
            out << scheme;
            throw std::runtime_error("EvolDF1nlep::Df1Evol(): scheme " + out.str()  + " not implemented ");
    }

    if (mu == this->mu && M == this->M && scheme == this->scheme && order_ew == NULL_ew)
        return (*Evol(order));

    if (mu == this->mu && M == this->M && scheme == this->scheme &&  order_ew != NULL_ew)
        return (*Evol(order_ew));


    if (M < mu) {
        std::stringstream out;
        out << "M = " << M << " < mu = " << mu;
        throw out.str();
    }

    setScales(mu, M); // also assign evol to identity
    
    

    double m_down = mu;
    //double m_up = model.AboveTh(m_down);
    //double nf = model.Nf(m_down);
    double nf = 5;                     //all the process in implemented for nf = 5
    alsM = alphatilde_s(M);
    alsmu = alphatilde_s(mu);
    
    eta = alsM / alsmu;
    logeta = log(eta);

    /*while (m_up < M) {                         //there is no thresholds
            Df1Evol(m_down, m_up, nf, scheme);
            //Df1threshold_nlep(m_up, nf+1.);
            m_down = m_up;
            m_up = model.AboveTh(m_down);
            nf += 1.;
    } */ 

    Df1Evol(m_down, M, nf, scheme);

    if(order_ew != NULL_ew) return (*Evol(order_ew));
    else return (*Evol(order)); 
    
}	

void EvolBsmm::Df1Evol(double mu, double M, double nf, schemes scheme)

{ 
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int k = 0;

    gslpp::matrix<double> resLO(dim, 0.);
    gslpp::matrix<double> Ueos(dim, 0.), Ue(dim, 0.), Ues(dim,0.), Us(dim,0.), Ue2os(dim, 0.), Ue2os2(dim, 0.), Ue2(dim,0.), Us2(dim,0.);

    int L = 6 - (int) nf;

    double eta = alsM / alsmu;

    double B00S = model.Beta0(nf);// B10S = model.Beta1(nf); /*inizializza i B*/

    double B00E = 80./9.;// B01E = 176./9.; /*inizializza i B*/

    double omega = 2. * B00S , lambda = B00E  /B00S ;    /*   E ale dipende da mu?*/


    for (k = 0; k < dim; k++) {
        double etap = pow(eta, a[L][k]);
        for (i = 0; i < dim; i++) {
            for (j = 0; j < dim; j++) {
                resLO(i, j) += b[L][i][j][k] * etap;
                Ue2(i, j) = (i == j);
            }
        }
    }

    unsigned int ind = 0;
    double max = 0;
    gsl_vector * list = gsl_vector_alloc (23);
    gsl_vector_set (list, 0, vbeevi.size()/7.);
    gsl_vector_set (list, 1, vebevi.size()/7.);
    gsl_vector_set (list, 2, veebvi.size()/7.);
    gsl_vector_set (list, 3, vbbevi.size()/7.);
    gsl_vector_set (list, 4, vbebvi.size()/7.);
    gsl_vector_set (list, 5, vebbvi.size()/7.);
    gsl_vector_set (list, 6, vaevi.size()/6.);
    gsl_vector_set (list, 7, vbbvi.size()/6.);
    gsl_vector_set (list, 8, vbdvi.size()/6.);
    gsl_vector_set (list, 9, vbevi.size()/6.);
    gsl_vector_set (list, 10, vdbvi.size()/6.);
    gsl_vector_set (list, 11, vdevi.size()/6.);
    gsl_vector_set (list, 12, veavi.size()/6.);
    gsl_vector_set (list, 13, vebvi.size()/6.);
    gsl_vector_set (list, 14, vedvi.size()/6.);
    gsl_vector_set (list, 15, veevi.size()/6.);
    gsl_vector_set (list, 16, vavi.size()/5.);
    gsl_vector_set (list, 17, vbvi.size()/5.);
    gsl_vector_set (list, 18, vcvi.size()/5.);
    gsl_vector_set (list, 19, vdvi.size()/5.);
    gsl_vector_set (list, 20, vevi.size()/5.);
    gsl_vector_set (list, 21, vfvi.size()/5.);
    gsl_vector_set (list, 22, vrvi.size()/5.);

    max = gsl_vector_max (list);

    for (ind = 0; ind < max; ind++) {

        if (ind < gsl_vector_get(list, 0)) {

            Ue2os(vbeevi[7 * ind], vbeevi[7 * ind + 5]) += vbeevi[7 * ind + 6] * H(vbeevi[7 * ind + 1], vbeevi[7 * ind + 2], vbeevi[7 * ind + 3], vbeevi[7 * ind + 4], 2, 4, 4, mu, M, nf);
        }
        if (ind < gsl_vector_get(list, 1)) {

            Ue2os(vebevi[7 * ind], vebevi[7 * ind + 5]) += vebevi[7 * ind + 6] * H(vebevi[7 * ind + 1], vebevi[7 * ind + 2], vebevi[7 * ind + 3], vebevi[7 * ind + 4], 4, 2, 4, mu, M, nf);
        }
        if (ind < gsl_vector_get(list, 2)) {

            Ue2os(veebvi[7 * ind], veebvi[7 * ind + 5]) += veebvi[7 * ind + 6] * H(veebvi[7 * ind + 1], veebvi[7 * ind + 2], veebvi[7 * ind + 3], veebvi[7 * ind + 4], 4, 4, 2, mu, M, nf);
        }
        if (ind < gsl_vector_get(list, 3)) {

            Ues(vbbevi[7 * ind], vbbevi[7 * ind + 5]) += vbbevi[7 * ind + 6] * H(vbbevi[7 * ind + 1], vbbevi[7 * ind + 2], vbbevi[7 * ind + 3], vbbevi[7 * ind + 4], 2, 2, 4, mu, M, nf);
        }
        if (ind < gsl_vector_get(list, 4)) {

            Ues(vbebvi[7 * ind], vbebvi[7 * ind + 5]) += vbebvi[7 * ind + 6] * H(vbebvi[7 * ind + 1], vbebvi[7 * ind + 2], vbebvi[7 * ind + 3], vbebvi[7 * ind + 4], 2, 4, 2, mu, M, nf);
        }
        if (ind < gsl_vector_get(list, 5)) {

            Ues(vebbvi[7 * ind], vebbvi[7 * ind + 5]) += vebbvi[7 * ind + 6] * H(vebbvi[7 * ind + 1], vebbvi[7 * ind + 2], vebbvi[7 * ind + 3], vebbvi[7 * ind + 4], 4, 2, 2, mu, M, nf);
        }


        if (ind < gsl_vector_get(list, 6)) {

            Ues(vaevi[6 * ind], vaevi[6 * ind + 4]) += vaevi[6 * ind + 5] * G(vaevi[6 * ind + 1], vaevi[6 * ind + 2], vaevi[6 * ind + 3], 1, 4, mu, M, nf);
        }
        if (ind < gsl_vector_get(list, 7)) {

            Us2(vbbvi[6 * ind], vbbvi[6 * ind + 4]) += vbbvi[6 * ind + 5] * G(vbbvi[6 * ind + 1], vbbvi[6 * ind + 2], vbbvi[6 * ind + 3], 2, 2, mu, M, nf);
        }
        if (ind < gsl_vector_get(list, 8)) {

            Ues(vbdvi[6 * ind], vbdvi[6 * ind + 4]) += vbdvi[6 * ind + 5] * G(vbdvi[6 * ind + 1], vbdvi[6 * ind + 2], vbdvi[6 * ind + 3], 2, 3, mu, M, nf);
        }
        if (ind < gsl_vector_get(list, 9)) {

            Ue(vbevi[6 * ind], vbevi[6 * ind + 4]) += vbevi[6 * ind + 5] * G(vbevi[6 * ind + 1], vbevi[6 * ind + 2], vbevi[6 * ind + 3], 2, 4, mu, M, nf);
            Ue2os(vbevi[6 * ind], vbevi[6 * ind + 4]) += vbevi[6 * ind + 5] * (-G(vbevi[6 * ind + 1], vbevi[6 * ind + 2], vbevi[6 * ind + 3], 2, 4, mu, M, nf)
                    + G(vbevi[6 * ind + 1], vbevi[6 * ind + 2], vbevi[6 * ind + 3], 2, 5, mu, M, nf));
        }
        if (ind < gsl_vector_get(list, 10)) {

            Ues(vdbvi[6 * ind], vdbvi[6 * ind + 4]) += vdbvi[6 * ind + 5] * G(vdbvi[6 * ind + 1], vdbvi[6 * ind + 2], vdbvi[6 * ind + 3], 3, 2, mu, M, nf);
        }
        if (ind < gsl_vector_get(list, 11)) {

            Ue2os(vdevi[6 * ind], vdevi[6 * ind + 4]) += vdevi[6 * ind + 5] * G(vdevi[6 * ind + 1], vdevi[6 * ind + 2], vdevi[6 * ind + 3], 3, 4, mu, M, nf);
        }
        if (ind < gsl_vector_get(list, 12)) {

            Ues(veavi[6 * ind], veavi[6 * ind + 4]) += veavi[6 * ind + 5] * G(veavi[6 * ind + 1], veavi[6 * ind + 2], veavi[6 * ind + 3], 4, 1, mu, M, nf);

        }
        if (ind < gsl_vector_get(list, 13)) {

            Ue(vebvi[6 * ind], vebvi[6 * ind + 4]) += vebvi[6 * ind + 5] * G(vebvi[6 * ind + 1], vebvi[6 * ind + 2], vebvi[6 * ind + 3], 4, 2, mu, M, nf);
            Ue2os(vebvi[6 * ind], vebvi[6 * ind + 4]) += vebvi[6 * ind + 5] * (-G(vebvi[6 * ind + 1], vebvi[6 * ind + 2], vebvi[6 * ind + 3], 4, 2, mu, M, nf)
                    + G(vebvi[6 * ind + 1], vebvi[6 * ind + 2], vebvi[6 * ind + 3], 5, 2, mu, M, nf));
        }
        if (ind < gsl_vector_get(list, 14)) {

            Ue2os(vedvi[6 * ind], vedvi[6 * ind + 4]) += vedvi[6 * ind + 5] * G(vedvi[6 * ind + 1], vedvi[6 * ind + 2], vedvi[6 * ind + 3], 4, 3, mu, M, nf);

        }
        if (ind < gsl_vector_get(list, 15)) {

            Ue2os2(veevi[6 * ind], veevi[6 * ind + 4]) += (veevi[6 * ind + 5] * G(veevi[6 * ind + 1], veevi[6 * ind + 2], veevi[6 * ind + 3], 4, 4, mu, M, nf));

        }


        if (ind < gsl_vector_get(list, 16)) {

            Us2(vavi[5 * ind], vavi[5 * ind + 3]) += (vavi[5 * ind + 4] * F(vavi[5 * ind + 1], vavi[5 * ind + 2], 1, mu, M, nf));
        }


        if (ind < gsl_vector_get(list, 17)) {

            Us(vbvi[5 * ind], vbvi[5 * ind + 3]) += vbvi[5 * ind + 4] * F(vbvi[5 * ind + 1], vbvi[5 * ind + 2], 2, mu, M, nf);

        }

        if (ind < gsl_vector_get(list, 18)) {

            Ues(vcvi[5 * ind], vcvi[5 * ind + 3]) += vcvi[5 * ind + 4] * F(vcvi[5 * ind + 1], vcvi[5 * ind + 2], 2, mu, M, nf);
        }
        if (ind < gsl_vector_get(list, 19)) {

            Ue(vdvi[5 * ind], vdvi[5 * ind + 3]) += vdvi[5 * ind + 4] * F(vdvi[5 * ind + 1], vdvi[5 * ind + 2], 3, mu, M, nf);
            Ue2os(vdvi[5 * ind], vdvi[5 * ind + 3]) += (lambda * lambda * omega)
                    * (-vdvi[5 * ind + 4] * F(vdvi[5 * ind + 1], vdvi[5 * ind + 2], 3, mu, M, nf));
        }
        if (ind < gsl_vector_get(list, 20)) {

            Ueos(vevi[5 * ind], vevi[5 * ind + 3]) += vevi[5 * ind + 4] * F(vevi[5 * ind + 1], vevi[5 * ind + 2], 4, mu, M, nf);
            Ue2os2(vevi[5 * ind], vevi[5 * ind + 3]) += (vevi[5 * ind + 4] * (F(vevi[5 * ind + 1], vevi[5 * ind + 2], 5, mu, M, nf)
                    - F(vevi[5 * ind + 1], vevi[5 * ind + 2], 4, mu, M, nf)));
        }
        if (ind < gsl_vector_get(list, 21)) {

            Ue2os(vfvi[5 * ind], vfvi[5 * ind + 3]) += (vfvi[5 * ind + 4] * F(vfvi[5 * ind + 1], vfvi[5 * ind + 2], 4, mu, M, nf));
        }
        if (ind < gsl_vector_get(list, 22)) {

            Ue2os(vrvi[5 * ind], vrvi[5 * ind + 3]) += vrvi[5 * ind + 4] * R(vrvi[5 * ind + 1], vrvi[5 * ind + 2], 4, mu, M, nf);

        }

    } 


    Us = omega * Us;
    Us2 = omega * omega * Us2;
    Ueos = lambda * Ueos;
    Ue = lambda * omega * Ue;
    Ue2os2 = lambda * lambda * Ue2os2;
    Ues = (omega * omega * lambda) * Ues; 
    Ue2os = (omega * lambda * lambda) * Ue2os;

    switch(order_ew) {    


        case NLO_ewt4:

            *elem[NLO_ewt4] = (*elem[NLO_ewt4]) * resLO + (*elem[NLO_ew]) * Ue + (*elem[LO]) * Ue2 + 
                    (*elem[NLO_ewt2]) * Ueos + (*elem[NNLO]) * Ue2os2 + (*elem[NLO]) * Ue2os;

        case NLO_ewt3:

            *elem[NLO_ewt3] =(*elem[NLO_ew]) * Ueos + (*elem[NLO]) * Ue2os2 + (*elem[LO]) * Ue2os;

        case NLO_ewt2:    

            *elem[NLO_ewt2] = (*elem[NLO_ewt2]) * resLO + (*elem[NLO_ew]) * Us +
                        (*elem[NLO]) * Ue + (*elem[LO]) * Ues + (*elem[NNLO]) * Ueos;

        case NLO_ewt1:   

            *elem[NLO_ewt1] = (*elem[LO]) * Ue2os2;

        case NLO_ew:

            *elem[NLO_ew] = (*elem[NLO_ew]) * resLO + (*elem[LO]) * Ue + (*elem[NLO]) * Ueos;


        case LO_ew:

            *elem[LO_ew] = (*elem[LO]) * Ueos;
            break;
            default:
            throw std::runtime_error("Error in EvolBsmm::Df1Evol");
    }   

    switch(order) {
        case NNLO:

            *elem[NNLO] =  (*elem[LO]) * Us2 + (*elem[NNLO]) * resLO + (*elem[NLO]) * Us;
        case NLO:

            *elem[NLO] = (*elem[LO]) * Us + (*elem[NLO]) * resLO;
        case LO:

            *elem[LO] = (*elem[LO]) * resLO;	
        break;
            default:
            throw std::runtime_error("Error in EvolBsmm::Df1Evol");
    } 
}

double EvolBsmm::F(unsigned int i, unsigned int j, int x, double mu, double M, double nf)

{
    int value = 0;
    int L = 6 - (int) nf;

    double etai = pow(eta, a[L][i]);    
    double etajx = pow(eta, a[L][j]+ x - 3.);
        
    double cut = 0.000000001;
    double result = 0.;

    if(fabs(a[L][j] + x - 3. - a[L][i]) < cut ) value = 0;
    else value = 1;	
    
    switch(value) {
        case 0:
            result = etai * logeta;
        break;
        case 1:
            result = (etajx - etai)/(a[L][j] + x - 3. - a[L][i]);
        break;
    }
    return (result);
}


double EvolBsmm::R(unsigned int i, unsigned int j, int x, double mu, double M, double nf)

{
    int value = 0;
    int L = 6 - (int) nf;

    double etai = pow(eta, a[L][i]);
    double etajx = pow(eta, a[L][j] + x - 3.);
        
    double cut = 0.000000001;
    double result = 0.;

    if(fabs(a[L][j] + x - 3. - a[L][i]) < cut ) value = 0;
    else value = 1;	
    
    switch(value) {
        case 0:
            result = etai * logeta * logeta * 0.5;
        break;
        case 1:
            result = (etajx * logeta - ((etajx - etai)/(a[L][j] + x - 3. - a[L][i])))/(a[L][j] + x - 3. - a[L][i]);
        break;
    }
return (result);
}

double EvolBsmm::G(unsigned int i, unsigned int p, unsigned int j, int x, int y, double mu, double M, double nf)

{
    int value = 0;
    int L = 6 - (int) nf;
    
    double etai = pow(eta, a[L][i]);
    double etapx = pow(eta, a[L][p]+ x - 3.);
    double etajxy = pow(eta, a[L][j]+ x - 3. + y - 3.);
        
    double cut = 0.000000001;
    double result = 0.;

    if(fabs(a[L][j] + y - 3. - a[L][p]) < cut && (a[L][p] + x - 3. - a[L][i]) < cut ) value = 0;
    else if(fabs(a[L][j] + y - 3. - a[L][p]) < cut && (a[L][p] + x - 3. - a[L][i]) > cut) value = 1;	
    else if(fabs(a[L][j] + y - 3. - a[L][p]) > cut && fabs(a[L][j] + y - 3. + x - 3. - a[L][i]) < cut &&  fabs(a[L][p] + x - 3. - a[L][i]) < cut ) value = 2;	
    else if(fabs(a[L][j] + y - 3. - a[L][p]) > cut && fabs(a[L][j] + y - 3. + x - 3. - a[L][i]) > cut &&  fabs(a[L][p] + x - 3. - a[L][i]) > cut ) value = 3;	
    else if(fabs(a[L][j] + y - 3. - a[L][p]) > cut && fabs(a[L][j] + y - 3. + x - 3. - a[L][i]) < cut &&  fabs(a[L][p] + x - 3. - a[L][i]) > cut ) value = 4;	
    else if(fabs(a[L][j] + y - 3. - a[L][p]) > cut && fabs(a[L][j] + y - 3. + x - 3. - a[L][i]) > cut &&  fabs(a[L][p] + x - 3. - a[L][i]) < cut ) value = 5;	

    switch(value) {
        case 0:
            result = etai * logeta * log(eta) * 0.5;
        break;
        case 1:
            result = ((etapx * logeta - ((etapx - etai)/(a[L][p] + x - 3. - a[L][i])))/(a[L][p] + x - 3. - a[L][i]));
        break;
        case 2:
            result = (etai * logeta - etai * logeta)/(a[L][j] + y - 3. - a[L][p]);
        break;
        case 3:
            result = (((etajxy - etai)/(a[L][j] + x - 3. + y - 3. - a[L][i])) - ((etapx - etai)/(a[L][p] + x - 3. - a[L][i])))/(a[L][j] + y - 3. - a[L][p]);
        break;
        case 4:
            result = (etai * logeta - ((etapx - etai)/(a[L][p] + x - 3. - a[L][i])))/(a[L][j] + y - 3. - a[L][p]);
        break;
        case 5:
            result = (((etajxy - etai)/(a[L][j] + x - 3. + y - 3. - a[L][i])) - etai * log(eta))/(a[L][j] + y - 3. - a[L][p]);
        break;
        }
    return (result);
}

double EvolBsmm::H(unsigned int i, unsigned int p, unsigned int q, unsigned int j, int x, int y, int z, double mu, double M, double nf)

{
    int value = 0;
    int L = 6 - (int) nf;

    double etai = pow(eta, a[L][i]);
    double etapx = pow(eta, a[L][p] + x - 3.);
    double etaqx = pow(eta, a[L][q] + x - 3.);

    double cut = 0.000000001;
    double result = 0.;

    if(fabs(a[L][p] + x - 3. - a[L][i]) < cut && fabs(a[L][q] + y - 3. - a[L][p]) < cut && fabs(a[L][j] + z - 3. - a[L][q]) < cut) value = 0;
    else if(fabs(a[L][p] + x - 3. - a[L][i]) > cut  && fabs(a[L][q] + y - 3. - a[L][p]) < cut && fabs(a[L][j] + z - 3. - a[L][q]) < cut) value = 1;
    else if((a[L][q] + x - 3. + y + 3 - a[L][i] ) < cut && fabs(a[L][j] + z - 3. - a[L][q]) < cut){
        if(fabs(a[L][p] + x - 3. - a[L][i]) < cut) value = 2;
        else value = 3;
    }
    else if((a[L][q] + x - 3. + y + 3 - a[L][i]) > cut && fabs(a[L][j] + z - 3. - a[L][q]) < cut){
        if(fabs(a[L][p] + x - 3. - a[L][i]) > cut) value = 4;
        else value = 5;
    }  
    else if(fabs(a[L][j] + z - 3. - a[L][q]) > cut) value = 6;

    switch(value) {
        case 0:
            result = etai * logeta * logeta * logeta /6.;
        break;
        case 1:
            result = (0.5 * etapx * logeta * logeta - ((etapx * logeta 
                    - ((etapx - etai)/(a[L][p] + x - 3. - a[L][i])))/(a[L][p] + x - 3. - a[L][i])))/(a[L][p] + x - 3. - a[L][i]);
        break;
        case 2:
            result = (etai * logeta * logeta * 0.5 - ((etai * logeta - etai * logeta)/(a[L][q] + y - 3. - a[L][p])))/(a[L][q] + y - 3. - a[L][p]);
        break;
        case 3:
            result = (etai * logeta * logeta * 0.5 - ((etai * logeta- ((etapx - etai)/(a[L][p] + x - 3. - a[L][i])) )/(a[L][q] + y - 3. - a[L][p])))/(a[L][q] + y - 3. - a[L][p]);
        break;
        case 4:
            result = (((etaqx * logeta - ((etaqx - etai)/(a[L][q] + x - 3. - a[L][i])))/(a[L][q] + y + 3. + x - 3. - a[L][i])) - ((etai * logeta - etai * logeta)/(a[L][q] + y - 3. - a[L][p])))/(a[L][q] + y - 3. - a[L][p]);
        break;
        case 5:
            result = (((etaqx * logeta - ((etaqx - etai)/(a[L][q] + x - 3. - a[L][i])))/(a[L][q] + y + 3. + x - 3. - a[L][i])) - ((etai * logeta - ((etapx - etai)/(a[L][p] + x - 3. - a[L][i])) )/(a[L][q] + y - 3. - a[L][p])))/(a[L][q] + y - 3. - a[L][p]);
        break;
        case 6:
            result = (G(i, p, j, x , y + z - 3., mu, M, nf) - G(i, p, q, x , y, mu, M, nf))/(a[L][j] + z - 3. - a[L][q]);		
        break;
    }
return (result);
}

double EvolBsmm::alphatilde_e(double mu)

{     // also the running is only for nf = 5
    
    //double mu_0 = 91.1876;
    double mu_0 = model.getMz();
    //double alphatilde_e = 1./(127.751 * 4. * M_PI); // alpha_e at mu_0 = 91.1876 Gev
    double alphatilde_e = model.alphaMz()/4./M_PI;
    //double alphatilde_s = 0.1184/(4.* M_PI); // alpha_s at mu_0 = 91.1876 Gev
    double alphatilde_s = model.getAlsMz()/4./M_PI;
    unsigned int nf = 5;

    double B00S = model.Beta0(nf), B10S = model.Beta1(nf);//B20S = model.Beta2(nf), 
    double B01S = -22./9.; //B11S = -308./27.; 

    double B00E = 80./9., B01E = 176./9., B10E = 464./27.;

    //double b1 = B10S/(2. * B00S * B00S); //b2 = B20S/(4. * B00S * B00S * B00S) - b1 * b1 , 
    //double b3 = B01S/(2. * B00S * B00E);// b4 = B11S /(4. * B00S * B00S * B00E) - 2 * b1 * b3, 
            //b5 = B01E/(2. * B00S * B00E) - b1;

    double vs= 1. + 2. * B00S * alphatilde_s * log(mu/ mu_0);
    double ve= 1. - 2. * B00E * alphatilde_e * log(mu/ mu_0);
    //double ps= B00S * alphatilde_s /(B00S * alphatilde_s + B00E * alphatilde_e);
    double pe= B00E * alphatilde_e /(B00S * alphatilde_s + B00E * alphatilde_e);

    double logve = log(ve);
    double logvs = log(vs);
    double logeos = log(ve/vs);
    double asovs = alphatilde_s/vs;
    double aeove = alphatilde_e/ve;

    double result = 0;

    result = aeove - pow(aeove, 2) * (logve * B10E/ B00E - logvs * B01E/B00S) 
            + pow(aeove, 2) * (asovs) * ((logvs - vs + 1.) * B01E * B10S/(B00S * B00S)
            + logve * vs * pe * B01E * B10E/(B00E * B00E) +(logeos * vs * pe - logve) * B01E * B01S/(B00S * B00E));
    return (result);
}

double EvolBsmm::alphatilde_s(double mu)

{  // also the running is only for nf = 5
    
    //double mu_0 = 91.1876;
    double mu_0 = model.getMz();
    //double alphatilde_e = 1./(127.751 * 4. * M_PI); // alpha_e at mu_0 = 91.1876 Gev
    double alphatilde_e = model.alphaMz()/4./M_PI;
    //double alphatilde_s = 0.1184/(4.* M_PI); // alpha_s at mu_0 = 91.1876 Gev
    double alphatilde_s = model.getAlsMz()/4./M_PI;
    unsigned int nf = 5;

    double B00S = model.Beta0(nf), B10S = model.Beta1(nf), B20S = model.Beta2(nf), B30S = gsl_sf_zeta_int(3) * 352864./81. - 598391./1458,
            B01S = -22./9., B11S = -308./27., B02S = 4945./243.; 

    double B00E = 80./9., B01E = 176./9., B10E = 464./27.; 

    //double B00S2 = B00S * B00S;
    double B10soB00s = B10S / B00S;
    double B01soB00e = B01S/B00E;

    //double b1 = B10soB00s/(2. * B00S), b2 = B20S/(4. * B00S2 * B00S) - b1 * b1 , 
    //        b3 = B01soB00e/(2. * B00S ), b4 = B11S /(4. * B00S2 * B00E) - 2 * b1 * b3, 
    //        b5 = B01E/(2. * B00S * B00E) - b1;

    double vs= 1. + 2. * B00S * alphatilde_s * log(mu/ mu_0);
    double ve= 1. - 2. * B00E * alphatilde_e * log(mu/ mu_0);
    double ps= B00S * alphatilde_s /(B00S * alphatilde_s + B00E * alphatilde_e);
    //double pe= B00E * alphatilde_e /(B00S * alphatilde_s + B00E * alphatilde_e);

    double logve = log(ve);
    double logvs = log(vs);
    double logeos = log(ve/vs);
    double logsoe = log(vs/ve);
    double asovs = alphatilde_s/vs;
    double aeove = alphatilde_e/ve;

    double result = 0;

    result = asovs - pow(asovs, 2) * (logvs * B10soB00s - logve * B01soB00e) 
            +  pow(asovs, 3) * ((1. - vs) * B20S / B00S + B10soB00s * B10soB00s * (logvs * logvs - logvs
            + vs - 1.) + B01soB00e * B01soB00e * logve * logve + (-2. * logvs * logve 
            + ps * ve * logve) * B01S * B10S/(B00E * B00S)) 
            +  pow(asovs, 4) * (0.5 * B30S *(1. - vs * vs)/ B00S + ((2. * vs - 3.) * logvs + vs * vs 
            - vs) * B20S * B10soB00s /(B00S) + B10soB00s * B10soB00s * B10soB00s * (- pow(logvs,3) 
            + 5. * pow(logvs,2) / 2. + 2. * (1. - vs) * logvs - (vs - 1.) * (vs - 1.)* 0.5)) 
            + pow(asovs, 2) * (aeove) * ((ve - 1.) * B02S / B00E 
            + ps * ve * logeos * B11S /B00S +(logve - ve + 1.) * B01soB00e * B10E/(B00S) 
            + logvs * ve * ps * B01S * B10soB00s/(B00S) +(logsoe * ve * ps - logvs) * B01soB00e * B01E/( B00S));
    return (result);
}
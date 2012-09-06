/* 
 * File:   sigmamuLEP2.cpp
 * Author: giovannigrilli
 * 
 * Created on 30 agosto 2012, 11.40
 */

#include "sigmamuLEP2.h"
    
    //const double sqrt_s[] = {130.,136.,161.,172.,183.,189.,192.,196.,200.,202.,205.,207.};

    static int SigmaLEP2_Integrand(const int *ndim, const double x[],
                                      const int *ncomp, double dsigma[], void *userdata) {
        
        double W,X,Y;      
        double cos_theta = -1.+2.*x[0];//sign
       
        if(myEW.checkModelForSTU()){
            W = myEW.W();
            X = myEW.X();
            Y = myEW.Y();
        } else {
            W = 0.;
            X = 0.;
            Y = 0.;
        }
    
        double s = sqrt_s_global*sqrt_s_global;
        dsigma[0] = 2.*myEW.dsigma_lLEP2(SM.MU,s,W,X,Y,cos_theta);
    
        return 0;
    }



double sigmamuLEP2::getThValue() {   
    
    sqrt_s_global = sqrt_s_local;

    const int NDIM = 1;
    const int NCOMP = 1;
    const double EPSREL = 1e-3;
    const double EPSABS = 1e-12;
    const int SEED = 0;
    const int MINEVAL = 0;
    const int MAXEVAL = 50000;
    const int NSTART = 1000;
    const int NINCREASE = 500;
    const int NBATCH = 1000; 
    const int GRIDNO = 0;
    double integral[NCOMP], error[NCOMP], prob[NCOMP];
    
    int verbose = 0;// No intermediate output
    //int verbose = 3;
    // outputs
    int neval, fail;
    
        Vegas(NDIM, NCOMP, SigmaLEP2_Integrand, NULL,
          EPSREL, EPSABS, verbose, SEED,
          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          GRIDNO, NULL,
          &neval, &fail, integral, error, prob); 
        
    return (integral[0]);
}


    
//    double sigmamu[12];
//    
//    double sqrt_s[12] = {130.,136.,161.,172.,183.,189.,192.,196.,200.,202.,205.,207.};
//    double s[12] = {130.*130.,136.*136.,161.*161.,172.*172.,183.*183.,189.*189.,
//                    192.*192.,196.*196.,200.*200.,202.*202.,205.*205.,207.*207.}; 
//    int i;

/* 
 * File:   sigmamuLEP2.cpp
 * Author: giovannigrilli
 * 
 * Created on 30 agosto 2012, 11.40
 */

#include "sigmamuLEP2.h"


static int Integrand(const int *ndim, const double theta, const double sqrt_s,
                     const double W, const double X, const double Y,
                     const int *ncomp, double dsigma[], void *userdata) {
    
    double s = sqrt_s*sqrt_s;
    dsigma[] = myEW.dsigma_lLEP2(SM.MU,s,W,X,Y,theta);
    
    return 0;
}


double sigmamuLEP2::getThValue() {   
    
    double sigmamu[12];
    double W,X,Y;
    //double sqrt_s[12] = {130.,136.,161.,172.,183.,189.,192.,196.,200.,202.,205.,207.};
    double s[12] = {130.*130.,136.*136.,161.*161.,172.*172.,183.*183.,189.*189.,
                    192.*192.,196.*196.,200.*200.,202.*202.,205.*205.,207.*207.}; 
    int i;
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
    double integral[NCOMP], error[NCOMP], prob[NCOMP];
    
    for(i=0;i<12;i++){
         
        if(myEW.checkModelForSTU()){
            W = myEW.W();
            X = myEW.X();
            Y = myEW.Y();
        } else {
            W = 0.;
            X = 0.;
            Y = 0.;
        }
        
        
        Vegas(NDIM, NCOMP, Integrand, NULL,
          EPSREL, EPSABS, verbose, SEED,
          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          GRIDNO, NULL,
          &neval, &fail, integral, error, prob); 
        
        sigmamu[i] = integral;
        
    }
    
    return (sigmamu);
    
}




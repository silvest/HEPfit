/* 
 * File:   sigmamuLEP2.cpp
 * Author: giovannigrilli
 * 
 * Created on 30 agosto 2012, 11.40
 */

#include "sigmamuLEP2.h"
    
    //const double sqrt_s[] = {130.,136.,161.,172.,183.,189.,192.,196.,200.,202.,205.,207.};
    class sigmamuLEP2 *global_this;
    
    double sigmamuLEP2::sigmaLEP2_Integrand(double x) {
        
        double W,X,Y;      
        double cos_theta = -1.+2.*x;
        
        double s = sqrt_s_local*sqrt_s_local;
       
        if(myEW.checkModelForSTU()){
            W = myEW.W();
            X = myEW.X();
            Y = myEW.Y();
        } else {
            W = 0.;
            X = 0.;
            Y = 0.;
        }
        
        return (2.*myEW.dsigma_lLEP2(SM.MU,s,W,X,Y,cos_theta));
        
        //return (2*cos_theta*cos_theta);
    }
    
    int sigmamuLEP2::Integrand(const int *ndim, const double x[],
                                      const int *ncomp, double dsigma[], void *userdata) {
        
        dsigma[0] = global_this->sigmaLEP2_Integrand(x[0]);
        
        return 0;
    }


double sigmamuLEP2::getThValue() {   
    
    global_this = this;

    integrand_t Integrand_pointer = &sigmamuLEP2::Integrand;
    
    const int NDIM = 1;
    const int NCOMP = 1;
    const double EPSREL = 1e-5;
    const double EPSABS = 1e-10;
    const int SEED = 0;
    const int MINEVAL = 0;
    const int MAXEVAL = 5000;
    const int NSTART = 1000;
    const int NINCREASE = 500;
    const int NBATCH = 1000; 
    const int GRIDNO = 0;
    double integral[NCOMP], error[NCOMP], prob[NCOMP];
    
    //int verbose = 0;// No intermediate output
    int verbose = 3;
    // outputs
    int neval, fail;
    //int dsigma = 2.*myEW.dsigma_lLEP2(SM.MU,s,W,X,Y,cos_theta);
  
    Vegas(NDIM, NCOMP, Integrand_pointer, NULL,
          EPSREL, EPSABS, verbose, SEED,
          MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
          GRIDNO, NULL,
          &neval, &fail, integral, error, prob); 
     
//    ////////////////////////////////////////
//    double s = sqrt_s_local*sqrt_s_local;
//    double cos_theta = 0.5;
//    double W = 0.;
//    double X = 0.;
//    double Y = 0.;
//
//    
//    return 2.*myEW.dsigma_lLEP2(SM.MU,s,W,X,Y,cos_theta);
//    //////////////////////////////////////////////////////
    
    std::cout << "error = " << error[0]*myEW.GeVminus2_to_nb()*1000. << "\tprob = " << prob[0] <<std::endl;
    
     return (integral[0]*myEW.GeVminus2_to_nb()*1000.);
}


    
//    double sigmamu[12];
//    
//    double sqrt_s[12] = {130.,136.,161.,172.,183.,189.,192.,196.,200.,202.,205.,207.};
//    double s[12] = {130.*130.,136.*136.,161.*161.,172.*172.,183.*183.,189.*189.,
//                    192.*192.,196.*196.,200.*200.,202.*202.,205.*205.,207.*207.}; 
//    int i;

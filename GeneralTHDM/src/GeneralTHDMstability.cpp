/* 
 * Copyright (C) 2016 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMstability.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"

stability_GTHDM::stability_GTHDM(const StandardModel& SM_i)
: myGTHDM(static_cast<const GeneralTHDM&> (SM_i)), vecMinus1(2,-1.), vecStability(2,0.),
        LambmatE(4,4,0.), Lambeigvec(4,4,0.), Lambeigval(4,0.)
{}

stability_GTHDM::~stability_GTHDM() 
{}





bool stability_GTHDM::CalcStabeigen(gslpp::matrix<gslpp::complex>& Stabeigvec_i, gslpp::vector<gslpp::complex>& Stabeigval_i)
{
    double lambda1 = myGTHDM.getlambda1();
    double lambda2 = myGTHDM.getlambda2();
    double lambda3 = myGTHDM.getlambda3();
    double lambda4 = myGTHDM.getlambda4();
    double Relambda5 = myGTHDM.getRelambda5();
    double Imlambda5 = myGTHDM.getImlambda5();
    double Relambda6 = myGTHDM.getRelambda6();
    double Relambda7 = myGTHDM.getRelambda7();
    double Imlambda6 = myGTHDM.getImlambda6();
    double Imlambda7 = myGTHDM.getImlambda7();
   
    LambmatE.assign(0,0, ((lambda1+lambda2)/2.0 + lambda3)/2.0);
    LambmatE.assign(0,1, (Relambda6 + Relambda7)/2.0);
    LambmatE.assign(0,2, -(Imlambda6 + Imlambda7)/2.0);
    LambmatE.assign(0,3, (lambda1 - lambda2)/4.0);
    LambmatE.assign(1,0, -(Relambda6 + Relambda7)/2.0);
    LambmatE.assign(1,1, -(lambda4 + Relambda5)/2.0);
    LambmatE.assign(1,2,  Imlambda5/2.0);
    LambmatE.assign(1,3, -(Relambda6 - Relambda7)/2.0);
    LambmatE.assign(2,0, (Imlambda6 + Imlambda7)/2.0);
    LambmatE.assign(2,1, Imlambda5/2.0);
    LambmatE.assign(2,2, -(lambda4 - Relambda5)/2.0);
    LambmatE.assign(2,3, (Imlambda6 - Imlambda7)/2.0);
    LambmatE.assign(3,0, -(lambda1 - lambda2)/4.0);
    LambmatE.assign(3,1, -(Relambda6 - Relambda7)/2.0);
    LambmatE.assign(3,2, (Imlambda6 - Imlambda7)/2.0);
    LambmatE.assign(3,3, (-(lambda1+lambda2)/2.0 + lambda3)/2.0);
   
    LambmatE.eigensystem(Stabeigvec_i, Stabeigval_i);
   
    return true;
}





gslpp::vector<double> stability_GTHDM::getStability()
{
    
    //We will include the conditions from Ivanov (hep-ph/0609018) for the boundness from below
    //The vacuum stability constraints are obtained similarly as showned by Ivanov and Silva (1507.05100)
    //Look at 2210.00024 for a recent review.
    
    CalcStabeigen(Lambeigvec,Lambeigval);
    
    //The eigenvectors must be real
    if((Lambeigval(0)).imag()<10^(-15) && (Lambeigval(1)).imag()<10^(-15) && (Lambeigval(2)).imag()<10^(-15) && (Lambeigval(3)).imag()<10^(-15)){
        
        //Now we need to find the position of the timelike eigenvector.
        //First we need to compute the modulus.
        double mink_mod_sq_0=(Lambeigvec(0,0).real())*(Lambeigvec(0,0).real())-((Lambeigvec(1,0).real())*(Lambeigvec(1,0).real())+(Lambeigvec(2,0).real())*(Lambeigvec(2,0).real())+(Lambeigvec(3,0).real())*(Lambeigvec(3,0).real()));
        double mink_mod_sq_1=(Lambeigvec(0,1).real())*(Lambeigvec(0,1).real())-((Lambeigvec(1,1).real())*(Lambeigvec(1,1).real())+(Lambeigvec(2,1).real())*(Lambeigvec(2,1).real())+(Lambeigvec(3,1).real())*(Lambeigvec(3,1).real()));
        double mink_mod_sq_2=(Lambeigvec(0,2).real())*(Lambeigvec(0,2).real())-((Lambeigvec(1,2).real())*(Lambeigvec(1,2).real())+(Lambeigvec(2,2).real())*(Lambeigvec(2,2).real())+(Lambeigvec(3,2).real())*(Lambeigvec(3,2).real()));
        double mink_mod_sq_3=(Lambeigvec(0,3).real())*(Lambeigvec(0,3).real())-((Lambeigvec(1,3).real())*(Lambeigvec(1,3).real())+(Lambeigvec(2,3).real())*(Lambeigvec(2,3).real())+(Lambeigvec(3,3).real())*(Lambeigvec(3,3).real()));
        
        //Now we find here the position of the timelike. 
        //There are for sure ways of ordering a vector but since it's short let's do it as follows which will be more efficient. This can be modified later.
        //Note that we really don't need to order them but just to find the time-like
        int postemp=-1;
	int possp1=-1;
	int possp2=-1;
	int possp3=-1;
	if(mink_mod_sq_0>0){
		postemp=0;
		possp1=1;
		possp2=2;
		possp3=3;
	}
	else{
		if(mink_mod_sq_1>0){
			possp1=0;
			postemp=1;
			possp2=2;
			possp3=3;
		}
		else{
			if(mink_mod_sq_2>0){
				possp1=0;
				possp2=1;
				postemp=2;
				possp3=3;
			}
			else{
				if(mink_mod_sq_3>0){
					possp1=0;
					possp2=1;
					possp3=2;
					postemp=3;
				}
			}
		}
        }     
        if(postemp==-1){
            return vecMinus1;
        }
        else{
            double Lambda0=(Lambeigval(postemp)).real();
            double Lambda1=(Lambeigval(possp1)).real();
            double Lambda2=(Lambeigval(possp2)).real();
            double Lambda3=(Lambeigval(possp3)).real();
            
            if(Lambda0>0){
                if(Lambda0>Lambda1 && Lambda0>Lambda2 && Lambda0>Lambda3){
                
                    vecStability(0)=Lambda0;
                    vecStability(1)=-1;
                    
                    double mHp1 = myGTHDM.getmHp();
                    double vh = myGTHDM.v();
                    
                    double xi_n=mHp1*mHp1/(vh*vh);
                    double DetS=-((Lambda0-xi_n)*(Lambda1-xi_n)*(Lambda2-xi_n)*(Lambda3-xi_n));
                    
                    if(DetS>0){
                        vecStability(1)= DetS;
                    }  
                    else{
                        if(xi_n-Lambda0>0){
                            vecStability(1)= xi_n-Lambda0;
                        }
                    }
                    
                    return vecStability;
                }
                else{
                    return vecMinus1;
                }
            }
            else{
                return vecMinus1;
            }
        } 
    }
    else{
        return vecMinus1;
    }
    
}









bounded_from_below_GTHDM::bounded_from_below_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double bounded_from_below_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getStability())(0);
}





vacuum_stability_GTHDM::vacuum_stability_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double vacuum_stability_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getStability())(1);
}
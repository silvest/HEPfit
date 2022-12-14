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
: myGTHDM(static_cast<const GeneralTHDM&> (SM_i)), vecMinus1(7,-1.), vecStability(7,0.),
        Lambmat(4,4,0.), Lambeigvec(4,4,0.), Lambeigval(4,0.)
{}

stability_GTHDM::~stability_GTHDM() 
{}

gslpp::vector<double> stability_GTHDM::getStability()
{
    //double lambda1 = myGTHDM.getMyGTHDMCache()->lambda1;
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
    bool flag_CPconservation=myGTHDM.getCPconservationflag();
    
    //std::cout<<"\033[1;34m   lambda1 =\033[0m "<< lambda1 <<std::endl;
    //std::cout<<"\033[1;34m   lambda2 =\033[0m "<< lambda2 <<std::endl;
    //std::cout<<"\033[1;34m   lambda3 =\033[0m "<< lambda3 <<std::endl;
    //std::cout<<"\033[1;34m   lambda4 =\033[0m "<< lambda4 <<std::endl;
    //std::cout<<"\033[1;34m   Relambda5 =\033[0m "<< Relambda5 <<std::endl;
    //std::cout<<"\033[1;34m   Imlambda5 =\033[0m "<< Imlambda5 <<std::endl;
    //std::cout<<"\033[1;34m   Relambda6 =\033[0m "<< Relambda6 <<std::endl;
    //std::cout<<"\033[1;34m   Imlambda6 =\033[0m "<< Imlambda6 <<std::endl;
    //std::cout<<"\033[1;34m   Relambda7 =\033[0m "<< Relambda7 <<std::endl;
    //std::cout<<"\033[1;34m   Imlambda7 =\033[0m "<< Imlambda7 <<std::endl;
    //std::cout<<"\033[1;34m   flag_CPconservation =\033[0m "<< flag_CPconservation <<std::endl;
    
    //vecStability(0) = lambda3 + sqrt(lambda1*lambda2);
    //vecStability(1) = lambda3 + lambda4 - sqrt(Relambda5*Relambda5+Imlambda5*Imlambda5) + sqrt(lambda1*lambda2);
    //vecStability(2) = (lambda1+lambda2)/2. + lambda3 + lambda4 + Relambda5 - 2.*fabs(Relambda6 + Relambda7);
    //vecStability(3) = (lambda1+lambda2)/2. + lambda3 + lambda4 - Relambda5 - 2.*fabs(Imlambda6 + Imlambda7);
    
    //The CP conserving expressions are taken from 2203.11462. Note that the lambda1 and lambda2 in that paper are 2*our lambdas, we need to devide by two
    
    if(flag_CPconservation){
        if (lambda1>0 && lambda2>0){
            if(Relambda6==0 && Relambda7==0){
                vecStability(0)=lambda3+2*sqrt((lambda1*lambda2)/4);
                vecStability(1)=lambda3+lambda4-abs(Relambda5)+sqrt(lambda1*lambda2);//we remove the factor 2 since we have a factor 1/4 inside the sqrt
                vecStability(2)= 1.;//In this case we have less conditions, let's just repeat the only two conditions 
                vecStability(3)= 1.;
                vecStability(4)= 1.;
                vecStability(5)= 1.;
                vecStability(6)= 1.;
                return vecStability;
            }
            else{
                double Relam6tRelam7=Relambda6*Relambda7;//This appears inside a sqrt, it must be higher or equal to zero
                double lam1tlam2=lambda1*lambda2;
                double lam3plam4pRelam5=lambda3+lambda3+Relambda5;
                double radicandplus=(lam1tlam2*lam3plam4pRelam5+lam1tlam2*sqrt(lam1tlam2));//This radicand appears inside a sqrt so it must be positive//In our notation this should be devided by 4 but we'll take it into account later since the expressions include 2*sqrt(radicandplus)
                if(Relam6tRelam7>=0 && radicandplus>=0){
                    //In this expressions we'll try to avoid the pow() function since the product is much faster (although it makes our expressions longer)
                    vecStability(0)=Relam6tRelam7;
                    vecStability(1)=-(lambda4+Relambda5)+2*sqrt(Relam6tRelam7);
                    vecStability(2)=-Relambda5+sqrt(Relam6tRelam7);
                    vecStability(3)=4*(4*lam1tlam2-12*Relam6tRelam7+lam3plam4pRelam5*lam3plam4pRelam5)*(4*lam1tlam2-12*Relam6tRelam7+lam3plam4pRelam5*lam3plam4pRelam5)*(4*lam1tlam2-12*Relam6tRelam7+lam3plam4pRelam5*lam3plam4pRelam5)
                            -(18*lam1tlam2*lam3plam4pRelam5+36*Relam6tRelam7*lam3plam4pRelam5-2*lam3plam4pRelam5*lam3plam4pRelam5*lam3plam4pRelam5-54*lambda1*Relambda7*Relambda7-54*Relambda6*Relambda6*lambda2)*(18*lam1tlam2*lam3plam4pRelam5+36*Relam6tRelam7*lam3plam4pRelam5-2*lam3plam4pRelam5*lam3plam4pRelam5*lam3plam4pRelam5-54*lambda1*Relambda7*Relambda7-54*Relambda6*Relambda6*lambda2);
                    vecStability(4)=sqrt(radicandplus)-abs(Relambda6*sqrt(lambda2/2)-Relambda7*sqrt(lambda1/2));
                    if(lam3plam4pRelam5-3*sqrt(lam1tlam2)<=0){
                        vecStability(5)=3*sqrt(lam1tlam2)-lam3plam4pRelam5;
                        vecStability(6)=sqrt(lam1tlam2)+lam3plam4pRelam5;
                        return vecStability;
                    }
                    else{
                        double radicandminus=(lam1tlam2*lam3plam4pRelam5-lam1tlam2*sqrt(lam1tlam2));
                        if(radicandminus>=0){
                            vecStability(5)=lam3plam4pRelam5-3*sqrt(lam1tlam2);
                            vecStability(6)=radicandminus-abs(Relambda6*sqrt(lambda2/2)+Relambda7*sqrt(lambda1/2));
                            return vecStability;
                        }
                        else
                            return vecMinus1;
                    }
                }
                else
                    return vecMinus1;
            }
        }
        else{
            return vecMinus1;
        }
    }
    else{
        throw std::runtime_error("The stability constraints are only included for the CP-conserving case.");
        return vecMinus1;
    }
}



bool stability_GTHDM::CalcStabeigen(gslpp::matrix<gslpp::complex>& Stabeigvec_i, gslpp::vector<double>& Stabeigval_i)
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
   
    Lambmat.assign(0,0, ((lambda1+lambda2)/2.0 + lambda3)/2.0);
    Lambmat.assign(0,1, (Relambda6 + Relambda7)/2.0);
    Lambmat.assign(0,2, -(Imlambda6 + Imlambda7)/2.0);
    Lambmat.assign(0,3, (lambda1 - lambda2)/4.0);
    Lambmat.assign(1,0, -(Relambda6 + Relambda7)/2.0);
    Lambmat.assign(1,1, -(lambda4 + Relambda5)/2.0);
    Lambmat.assign(1,2,  Imlambda5/2.0);
    Lambmat.assign(1,3, -(Relambda6 - Relambda7)/2.0);
    Lambmat.assign(2,0, (Imlambda6 + Imlambda7)/2.0);
    Lambmat.assign(2,1, Imlambda5/2.0);
    Lambmat.assign(2,2, -(lambda4 - Relambda5)/2.0);
    Lambmat.assign(2,3, (Imlambda6 - Imlambda7)/2.0);
    Lambmat.assign(3,0, -(lambda1 - lambda2)/4.0);
    Lambmat.assign(3,1, -(Relambda6 - Relambda7)/2.0);
    Lambmat.assign(3,2, (Imlambda6 - Imlambda7)/2.0);
    Lambmat.assign(3,3, (-(lambda1+lambda2)/2.0 + lambda3)/2.0);
   
    Lambmat.eigensystem(Stabeigvec_i, Stabeigval_i);
   
    return true;
}


double stability_GTHDM::getVacuumStability()
{
    //Let's first check that the stability conditions (boundness from below) are satisfied
    if((getStability())(0)>=0 && (getStability())(1)>=0 && (getStability())(2)>=0 && (getStability())(3)>=0 && (getStability())(4)>=0 && (getStability())(5)>=0 && (getStability())(6)>=0){
        CalcStabeigen(Lambeigvec,Lambeigval);
        

        double mHp1 = myGTHDM.getmHp();
        double vh = myGTHDM.v();
        
        double Lambda_0=Lambeigval.max();
        double xi_n=mHp1*mHp1/(vh*vh);
        double DetS=(-(xi_n-Lambeigval(0))*(xi_n-Lambeigval(1))*(xi_n-Lambeigval(2))*(xi_n-Lambeigval(3)));
        
        if(DetS>0){
            return DetS;
        }  
        else{
            if(xi_n-Lambda_0>0){
                return xi_n-Lambda_0;
            }
            else
                return -1.;
        }
    }
    else
        return -1;

}




stability1_GTHDM::stability1_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double stability1_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getStability())(0);
}


stability2_GTHDM::stability2_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double stability2_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getStability())(1);
}


stability3_GTHDM::stability3_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double stability3_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getStability())(2);
}


stability4_GTHDM::stability4_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double stability4_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getStability())(3);
}


stability5_GTHDM::stability5_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double stability5_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getStability())(4);
}


stability6_GTHDM::stability6_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double stability6_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getStability())(5);
}



stability7_GTHDM::stability7_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double stability7_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getStability())(6);
}



vacuumstability_GTHDM::vacuumstability_GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i), mystability_GTHDM(SM_i)
{}

double vacuumstability_GTHDM::computeThValue()
{
    return (mystability_GTHDM.getVacuumStability());
}
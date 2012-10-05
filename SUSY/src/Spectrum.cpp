/* 
 * File:   Spectrum.cpp
 * Author: silvest
 * 
 * Created on October 5, 2012, 4:22 PM
 */

#include "Spectrum.h"
#include "SUSY.h"
#include <gslpp.h>

using namespace gslpp;

Spectrum::Spectrum() {
}

void Spectrum::CalcSpectrum(SUSY & SUSY){
    matrix<complex> CKM(SUSY.getVCKM());
    matrix<double> Mu(3,3,0.);
    
    Mu(1,1) = SUSY.Mrun(SUSY.GetQ(),SUSY.getQuarks(QCD::CHARM).getMass());
    Mu(2,2) = SUSY.Mrun(SUSY.GetQ(),SUSY.getQuarks(QCD::TOP).getMass());
    //Up-type squark masses
    matrix<complex> uLL(CKM*SUSY.MsQ2*CKM.hconjugate()+Mu*Mu+matrix<complex>::Id(3)*(.5-2./3.*SUSY.s02())*SUSY.getMz()*SUSY.getMz()*cos(2.*atan(SUSY.getTanb())));
    matrix<complex> uRR(SUSY.MsU2+Mu*Mu+matrix<complex>::Id(3)*(-1./3.*SUSY.s02())*SUSY.getMz()*SUSY.getMz()*cos(2.*atan(SUSY.getTanb())));
    matrix<complex> uLR(SUSY.v1()/sqrt(2.)*SUSY.GetTU().hconjugate()-SUSY.getMuH()*Mu/SUSY.getTanb());
    
    matrix<complex> M2(6,6,0.);
    
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++){
            M2.assign(i,j,uLL(i,j));
            M2.assign(i+3,j,uLR(i,j));
            M2.assign(i,j+3,uLR(j,i).conjugate());
            M2.assign(i+3,j+3,uRR(i,j));
        }
    M2.eigensystem(SUSY.Ru,SUSY.Msu2);
    
}

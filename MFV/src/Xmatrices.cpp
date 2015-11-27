/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Xmatrices.h"

Xmatrices::Xmatrices(): X1(3,3,0.), X2(3,3,0.), 
        X3(3,3,0.), X4(3,3,0.), X5(3,3,0.), X6(3,3,0.), X9(3,3,0.), X13(3,3,0.) 
{
    X1.assign(2,2,1.);
    X2.assign(1,1,1.);
    X3.assign(2,1,1.);
    X4.assign(1,2,1.);
}

void Xmatrices::Update(const CKM & CKM_in){
    if(CKM_in.getA() == myA && CKM_in.getRho() == myrho && CKM_in.getEta() == 
            myeta && CKM_in.getLambda() == mylambda)
        return;
    
    gslpp::matrix<gslpp::complex> ckm(3,3,0.);
    CKM_in.getCKM(ckm);
    
    for (int i = 0; i < 3; i++){
        X5.assign(2,i,ckm(2,i));
        X6.assign(1,i,ckm(1,i));
//        X7.assign(2,i,ckm(1,i));
//        X8.assign(1,i,ckm(2,i));
        X9.assign(i,2,ckm(2,i).conjugate());
//        X10.assign(i,1,ckm(1,i).conjugate());
//        X11.assign(i,1,ckm(2,i).conjugate());
//        X12.assign(i,2,ckm(1,i).conjugate());
        for (int j = 0; j < 3; j++ ){
            X13.assign(i,j,ckm(2,i).conjugate()*ckm(2,j));
//            X14.assign(i,j,ckm(1,i).conjugate()*ckm(1,j));
//            X15.assign(i,j,ckm(2,i).conjugate()*ckm(1,j));
//            X16.assign(i,j,ckm(1,i).conjugate()*ckm(2,j));
        }
    }
    myA = CKM_in.getA();
    myrho = CKM_in.getRho();
    myeta = CKM_in.getEta();
    mylambda = CKM_in.getLambda();
}


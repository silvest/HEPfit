/* 
 * Copyright (C) 2017 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "gslpp.h"
#include "GMcache.h"

GMcache::GMcache(const StandardModel& SM_i)
        : unitarityeigenvalues(4, 0.), myGM(static_cast<const GeorgiMachacek*> (&SM_i))
{
//    myRunnerGM=new RunnerGM(SM_i);
}

GMcache::~GMcache()
{
//  delete myRunnerGM;
}
//
///////////////////////////////////////////////////////////////////////////////////////////////////
//
int GMcache::CacheCheck(const gslpp::complex cache[][CacheSize], 
                          const int NumPar, const double params[]) const {
    bool bCache;
    for(int i=0; i<CacheSize; i++) {
        bCache = true;
        for(int j=0; j<NumPar; j++)
            bCache &= (params[j] == cache[j][i].real());
        if (bCache) return i;
    }
    return -1;
}
//
//int GMcache::CacheCheckReal(const double cache[][CacheSize], 
//                          const int NumPar, const double params[]) const {
//    bool bCache;
//    for(int i=0; i<CacheSize; i++) {
//        bCache = true;
//        for(int j=0; j<NumPar; j++)
//            bCache &= (params[j] == cache[j][i]);
//        if (bCache) return i;
//    }
//    return -1;
//}

void GMcache::CacheShift(gslpp::complex cache[][CacheSize], const int NumPar, 
                           const double params[], const gslpp::complex newResult) const {
    // shift old parameters and result
    for(int i=CacheSize-1; i>0; i--)
        for(int j=0; j<NumPar+1; j++)
            cache[j][i] = cache[j][i-1];

    // store new parameters and result
    for(int j=0; j<NumPar; j++) {
        cache[j][0] = gslpp::complex(params[j], 0.0, false);
        cache[NumPar][0] = newResult;
    }
}

gslpp::complex GMcache::I_h_U(const double mHl2, const double Mu, const double Mc, const double Mt) const {
    int NumPar = 4;
    double params[] = {mHl2, Mu, Mc, Mt};

    int i = CacheCheck(I_h_U_cache, NumPar, params);
    if (i>=0) {
        return ( I_h_U_cache[NumPar][i] );
    } else {
    	double TAUu=4.0*Mu*Mu/mHl2;
    	double TAUc=4.0*Mc*Mc/mHl2;
    	double TAUt=4.0*Mt*Mt/mHl2;
        gslpp::complex newResult = -(8./3.)*(TAUu*(1.0+(1.0-TAUu)*f_func(TAUu))
                         +TAUc*(1.0+(1.0-TAUc)*f_func(TAUc))+TAUt*(1.0+(1.0-TAUt)*f_func(TAUt)));
        CacheShift(I_h_U_cache, NumPar, params, newResult);
        return newResult;
    }
}
//
//gslpp::complex GMcache::I_HH_U(const double mHh2, const double Mc, const double Mt) const {
//    int NumPar = 3;
//    double params[] = {mHh2, Mc, Mt};
//
//    int i = CacheCheck(I_HH_U_cache, NumPar, params);
//    if (i>=0) {
//        return ( I_HH_U_cache[NumPar][i] );
//    } else {
//    	double TAUc=4.0*Mc*Mc/mHh2;
//    	double TAUt=4.0*Mt*Mt/mHh2;
//        gslpp::complex newResult = -(8./3.)*(TAUc*(1.0+(1.0-TAUc)*f_func(TAUc))
//                      +TAUt*(1.0+(1.0-TAUt)*f_func(TAUt)));
//        CacheShift(I_HH_U_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//gslpp::complex GMcache::I_A_U(const double mA2, const double Mc, const double Mt) const {
//    int NumPar = 3;
//    double params[] = {mA2, Mc, Mt};
//
//    int i = CacheCheck(I_A_U_cache, NumPar, params);
//    if (i>=0) {
//        return ( I_A_U_cache[NumPar][i] );
//    } else {
//    	double TAUc=4.0*Mc*Mc/mA2;
//    	double TAUt=4.0*Mt*Mt/mA2;
//        gslpp::complex newResult = -(8./3.)*(TAUc*f_func(TAUc)+TAUt*f_func(TAUt));
//        CacheShift(I_A_U_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}

gslpp::complex GMcache::I_h_D(const double mHl2, const double Md, const double Ms, const double Mb) const {
    int NumPar = 4;
    double params[] = {mHl2, Md, Ms, Mb};

    int i = CacheCheck(I_h_D_cache, NumPar, params);
    if (i>=0) {
        return ( I_h_D_cache[NumPar][i] );
    } else {
    	double TAUd=4.0*Md*Md/mHl2;
    	double TAUs=4.0*Ms*Ms/mHl2;
    	double TAUb=4.0*Mb*Mb/mHl2;
        gslpp::complex newResult = -(2./3.)*(TAUd*(1.0+(1.0-TAUd)*f_func(TAUd))
                         +TAUs*(1.0+(1.0-TAUs)*f_func(TAUs))+TAUb*(1.0+(1.0-TAUb)*f_func(TAUb)));
        CacheShift(I_h_D_cache, NumPar, params, newResult);
        return newResult;
    }
}
//
//gslpp::complex GMcache::I_HH_D(const double mHh2, const double Ms, const double Mb) const {
//    int NumPar = 3;
//    double params[] = {mHh2, Ms, Mb};
//
//    int i = CacheCheck(I_HH_D_cache, NumPar, params);
//    if (i>=0) {
//        return ( I_HH_D_cache[NumPar][i] );
//    } else {
//    	double TAUs=4.0*Ms*Ms/mHh2;
//    	double TAUb=4.0*Mb*Mb/mHh2;
//        gslpp::complex newResult = -(2./3.)*(TAUs*(1.0+(1.0-TAUs)*f_func(TAUs))
//                      +TAUb*(1.0+(1.0-TAUb)*f_func(TAUb)));
//        CacheShift(I_HH_D_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//gslpp::complex GMcache::I_A_D(const double mA2, const double Ms, const double Mb) const {
//    int NumPar = 3;
//    double params[] = {mA2, Ms, Mb};
//
//    int i = CacheCheck(I_A_D_cache, NumPar, params);
//    if (i>=0) {
//        return ( I_A_D_cache[NumPar][i] );
//    } else {
//    	double TAUs=4.0*Ms*Ms/mA2;
//    	double TAUb=4.0*Mb*Mb/mA2;
//        gslpp::complex newResult = -(2./3.)*(TAUs*f_func(TAUs)+TAUb*f_func(TAUb));
//        CacheShift(I_A_D_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}

gslpp::complex GMcache::I_h_L(const double mHl2, const double Me, const double Mmu, const double Mtau) const {
    int NumPar = 4;
    double params[] = {mHl2, Me, Mmu, Mtau};

    int i = CacheCheck(I_h_L_cache, NumPar, params);
    if (i>=0) {
        return ( I_h_L_cache[NumPar][i] );
    } else {
    	double TAUe=4.0*Me*Me/mHl2;
    	double TAUmu=4.0*Mmu*Mmu/mHl2;
    	double TAUtau=4.0*Mtau*Mtau/mHl2;
        gslpp::complex newResult = -2.0*(TAUe*(1.0+(1.0-TAUe)*f_func(TAUe))
                         +TAUmu*(1.0+(1.0-TAUmu)*f_func(TAUmu))
                         +TAUtau*(1.0+(1.0-TAUtau)*f_func(TAUtau)));
        CacheShift(I_h_L_cache, NumPar, params, newResult);
        return newResult;
    }
}
//
//gslpp::complex GMcache::I_HH_L(const double mHh2, const double Mmu, const double Mtau) const {
//    int NumPar = 3;
//    double params[] = {mHh2, Mmu, Mtau};
//
//    int i = CacheCheck(I_HH_L_cache, NumPar, params);
//    if (i>=0) {
//        return ( I_HH_L_cache[NumPar][i] );
//    } else {
//    	double TAUmu=4.0*Mmu*Mmu/mHh2;
//    	double TAUtau=4.0*Mtau*Mtau/mHh2;
//        gslpp::complex newResult = -2.0*(TAUmu*(1.0+(1.0-TAUmu)*f_func(TAUmu))+
//                           TAUtau*(1.0+(1.0-TAUtau)*f_func(TAUtau)));
//        CacheShift(I_HH_L_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//gslpp::complex GMcache::I_A_L(const double mA2, const double Mmu, const double Mtau) const {
//    int NumPar = 3;
//    double params[] = {mA2, Mmu, Mtau};
//
//    int i = CacheCheck(I_A_L_cache, NumPar, params);
//    if (i>=0) {
//        return ( I_A_L_cache[NumPar][i] );
//    } else {
//    	double TAUmu=4.0*Mmu*Mmu/mA2;
//    	double TAUtau=4.0*Mtau*Mtau/mA2;
//        gslpp::complex newResult = -2.0*(TAUmu*f_func(TAUmu)+TAUtau*f_func(TAUtau));
//        CacheShift(I_A_L_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}

gslpp::complex GMcache::I_H_W(const double mH, const double MW) const {
    int NumPar = 2;
    double params[] = {mH, MW};

    int i = CacheCheck(I_H_W_cache, NumPar, params);
    if (i>=0) {
        return ( I_H_W_cache[NumPar][i] );
    } else {
        double TAUw=4.0*MW*MW/(mH*mH);
        gslpp::complex newResult = 2.0 + 3.0*TAUw + 3.0*TAUw*(2.0-TAUw)*f_func(TAUw);
        CacheShift(I_H_W_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GMcache::I_H_Hp(const double mHp2, const double mH) const {
    int NumPar = 2;
    double params[] = {mHp2, mH};

    int i = CacheCheck(I_H_Hp_cache, NumPar, params);
    if (i>=0) {
        return ( I_H_Hp_cache[NumPar][i] );
    } else {
        double TAUhp=4.0*mHp2/(mH*mH);
        gslpp::complex newResult = -TAUhp*(1.0-TAUhp*f_func(TAUhp));
        CacheShift(I_H_Hp_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GMcache::A_h_U(const double mHl2, const double cW2, const double Mu, const double Mc, const double Mt, const double MZ) const {
    int NumPar = 6;
    double params[] = {mHl2, cW2, Mu, Mc, Mt, MZ};

    int i = CacheCheck(A_h_U_cache, NumPar, params);
    if (i>=0) {
        return ( A_h_U_cache[NumPar][i] );
    } else {
    	double TAUu=4.0*Mu*Mu/mHl2;
    	double TAUc=4.0*Mc*Mc/mHl2;
    	double TAUt=4.0*Mt*Mt/mHl2;
    	double LAMu=4.0*Mu*Mu/(MZ*MZ);
    	double LAMc=4.0*Mc*Mc/(MZ*MZ);
    	double LAMt=4.0*Mt*Mt/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = -4.0*(1.0/2.0-4.0/3.0*sW2)*(Int1(TAUu,LAMu)+Int1(TAUc,LAMc)
                           +Int1(TAUt,LAMt)-Int2(TAUu,LAMu)-Int2(TAUc,LAMc)-Int2(TAUt,LAMt));
        CacheShift(A_h_U_cache, NumPar, params, newResult);
        return newResult;
    }
}
//
//gslpp::complex GMcache::A_HH_U(const double mHh2, const double cW2, const double Mc, const double Mt, const double MZ) const {
//    int NumPar = 5;
//    double params[] = {mHh2, cW2, Mc, Mt, MZ};
//
//    int i = CacheCheck(A_HH_U_cache, NumPar, params);
//    if (i>=0) {
//        return ( A_HH_U_cache[NumPar][i] );
//    } else {
//    	double TAUc=4.0*Mc*Mc/mHh2;
//    	double TAUt=4.0*Mt*Mt/mHh2;
//    	double LAMc=4.0*Mc*Mc/(MZ*MZ);
//    	double LAMt=4.0*Mt*Mt/(MZ*MZ);
//	double sW2=1.0-cW2;
//        gslpp::complex newResult = -4.0*(1.0/2.0-4.0/3.0*sW2)*(Int1(TAUc,LAMc)-Int2(TAUc,LAMc)
//                                         +Int1(TAUt,LAMt)-Int2(TAUt,LAMt));
//        CacheShift(A_HH_U_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//gslpp::complex GMcache::A_A_U(const double mA2, const double cW2, const double Mc, const double Mt, const double MZ) const {
//    int NumPar = 5;
//    double params[] = {mA2, cW2, Mc, Mt, MZ};
//
//    int i = CacheCheck(A_A_U_cache, NumPar, params);
//    if (i>=0) {
//        return ( A_A_U_cache[NumPar][i] );
//    } else {
//    	double TAUc=4.0*Mc*Mc/mA2;
//    	double TAUt=4.0*Mt*Mt/mA2;
//    	double LAMc=4.0*Mc*Mc/(MZ*MZ);
//    	double LAMt=4.0*Mt*Mt/(MZ*MZ);
//	double sW2=1.0-cW2;
//        gslpp::complex newResult = -4.0*(1.0/2.0-4.0/3.0*sW2)*(-Int2(TAUc,LAMc)-Int2(TAUt,LAMt))/sqrt(sW2*cW2);
//        CacheShift(A_A_U_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}

gslpp::complex GMcache::A_h_D(const double mHl2, const double cW2, const double Md, const double Ms, const double Mb, const double MZ) const {
    int NumPar = 6;
    double params[] = {mHl2, cW2, Md, Ms, Mb, MZ};

    int i = CacheCheck(A_h_D_cache, NumPar, params);
    if (i>=0) {
        return ( A_h_D_cache[NumPar][i] );
    } else {
    	double TAUd=4.0*Md*Md/mHl2;
    	double TAUs=4.0*Ms*Ms/mHl2;
    	double TAUb=4.0*Mb*Mb/mHl2;
    	double LAMd=4.0*Md*Md/(MZ*MZ);
    	double LAMs=4.0*Ms*Ms/(MZ*MZ);
	double LAMb=4.0*Mb*Mb/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0/3.0*sW2)*(Int1(TAUd,LAMd)+Int1(TAUs,LAMs)
                           +Int1(TAUb,LAMb)-Int2(TAUd,LAMd)-Int2(TAUs,LAMs)-Int2(TAUb,LAMb));
        CacheShift(A_h_D_cache, NumPar, params, newResult);
        return newResult;
    }
}
//
//gslpp::complex GMcache::A_HH_D(const double mHh2, const double cW2, const double Ms, const double Mb, const double MZ) const {
//    int NumPar = 5;
//    double params[] = {mHh2, cW2, Ms, Mb, MZ};
//
//    int i = CacheCheck(A_HH_D_cache, NumPar, params);
//    if (i>=0) {
//        return ( A_HH_D_cache[NumPar][i] );
//    } else {
//    	double TAUs=4.0*Ms*Ms/mHh2;
//    	double TAUb=4.0*Mb*Mb/mHh2;
//    	double LAMs=4.0*Ms*Ms/(MZ*MZ);
//	double LAMb=4.0*Mb*Mb/(MZ*MZ);
//	double sW2=1.0-cW2;
//        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0/3.0*sW2)*(Int1(TAUs,LAMs)-Int2(TAUs,LAMs)
//                                          +Int1(TAUb,LAMb)-Int2(TAUb,LAMb));
//        CacheShift(A_HH_D_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//gslpp::complex GMcache::A_A_D(const double mA2, const double cW2, const double Ms, const double Mb, const double MZ) const {
//    int NumPar = 5;
//    double params[] = {mA2, cW2, Ms, Mb, MZ};
//
//    int i = CacheCheck(A_A_D_cache, NumPar, params);
//    if (i>=0) {
//        return ( A_A_D_cache[NumPar][i] );
//    } else {
//    	double TAUs=4.0*Ms*Ms/mA2;
//    	double TAUb=4.0*Mb*Mb/mA2;
//    	double LAMs=4.0*Ms*Ms/(MZ*MZ);
//	double LAMb=4.0*Mb*Mb/(MZ*MZ);
//	double sW2=1.0-cW2;
//        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0/3.0*sW2)*(-Int2(TAUs,LAMs)-Int2(TAUb,LAMb))/sqrt(sW2*cW2);
//        CacheShift(A_A_D_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}

gslpp::complex GMcache::A_h_L(const double mHl2, const double cW2, const double Me, const double Mmu, const double Mtau, const double MZ) const {
    int NumPar = 6;
    double params[] = {mHl2, cW2, Me, Mmu, Mtau, MZ};

    int i = CacheCheck(A_h_L_cache, NumPar, params);
    if (i>=0) {
        return ( A_h_L_cache[NumPar][i] );
    } else {
    	double TAUe=4.0*Me*Me/mHl2;
    	double TAUmu=4.0*Mmu*Mmu/mHl2;
    	double TAUtau=4.0*Mtau*Mtau/mHl2;
    	double LAMe=4.0*Me*Me/(MZ*MZ);
    	double LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
	double LAMtau=4.0*Mtau*Mtau/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0*sW2)*(Int1(TAUe,LAMe)+Int1(TAUmu,LAMmu)
                            +Int1(TAUtau,LAMtau)-Int2(TAUe,LAMe)-Int2(TAUmu,LAMmu)
                            -Int2(TAUtau,LAMtau));
        CacheShift(A_h_L_cache, NumPar, params, newResult);
        return newResult;
    }
}
//
//gslpp::complex GMcache::A_HH_L(const double mHh2, const double cW2, const double Mmu, const double Mtau, const double MZ) const {
//    int NumPar = 5;
//    double params[] = {mHh2, cW2, Mmu, Mtau, MZ};
//
//    int i = CacheCheck(A_HH_L_cache, NumPar, params);
//    if (i>=0) {
//        return ( A_HH_L_cache[NumPar][i] );
//    } else {
//    	double TAUmu=4.0*Mmu*Mmu/mHh2;
//    	double TAUtau=4.0*Mtau*Mtau/mHh2;
//    	double LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
//	double LAMtau=4.0*Mtau*Mtau/(MZ*MZ);
//	double sW2=1.0-cW2;
//        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0*sW2)*(Int1(TAUmu,LAMmu)-Int2(TAUmu,LAMmu)
//                                      +Int1(TAUtau,LAMtau)-Int2(TAUtau,LAMtau));
//        CacheShift(A_HH_L_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}
//
//gslpp::complex GMcache::A_A_L(const double mA2, const double cW2, const double Mmu, const double Mtau, const double MZ) const {
//    int NumPar = 5;
//    double params[] = {mA2, cW2, Mmu, Mtau, MZ};
//
//    int i = CacheCheck(A_A_L_cache, NumPar, params);
//    if (i>=0) {
//        return ( A_A_L_cache[NumPar][i] );
//    } else {
//    	double TAUmu=4.0*Mmu*Mmu/mA2;
//    	double TAUtau=4.0*Mtau*Mtau/mA2;
//    	double LAMmu=4.0*Mmu*Mmu/(MZ*MZ);
//	double LAMtau=4.0*Mtau*Mtau/(MZ*MZ);
//	double sW2=1.0-cW2;
//        gslpp::complex newResult = 2.0*(-1.0/2.0+2.0*sW2)*(-Int2(TAUmu,LAMmu)-Int2(TAUtau,LAMtau))/sqrt(sW2*cW2);
//        CacheShift(A_A_L_cache, NumPar, params, newResult);
//        return newResult;
//    }
//}

gslpp::complex GMcache::A_H_W(const double mH, const double cW2, const double MW, const double MZ) const {
    int NumPar = 4;
    double params[] = {mH, cW2, MW, MZ};

    int i = CacheCheck(A_H_W_cache, NumPar, params);
    if (i>=0) {
        return ( A_H_W_cache[NumPar][i] );
    } else {
        double TAUw=4.0*MW*MW/(mH*mH);
        double LAMw=4.0*MW*MW/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = -sqrt(cW2/sW2)*(4.0*(3.0-sW2/cW2)*Int2(TAUw,LAMw)
                            +((1.0+2.0/TAUw)*sW2/cW2-(5.0+2.0/TAUw))*Int1(TAUw,LAMw));
        CacheShift(A_H_W_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GMcache::A_H_Hp(const double mHp2, const double mH, const double cW2, const double MZ) const {
    int NumPar = 4;
    double params[] = {mHp2, mH, cW2, MZ};

    int i = CacheCheck(A_H_Hp_cache, NumPar, params);
    if (i>=0) {
        return ( A_H_Hp_cache[NumPar][i] );
    } else {
        double TAUhp=4.0*mHp2/(mH*mH);
        double LAMhp=4.0*mHp2/(MZ*MZ);
	double sW2=1.0-cW2;
        gslpp::complex newResult = (1.0-2.0*sW2)/sqrt(cW2*sW2)*Int1(TAUhp,LAMhp);
        CacheShift(A_H_Hp_cache, NumPar, params, newResult);
        return newResult;
    }
}

gslpp::complex GMcache::f_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    return -pow(log((1.0+sqrt(1.0-x))/(1.0-sqrt(1.0-x)))+z,2)/4.0;
    }
    else {
        return pow(asin(sqrt(1.0/x)),2);
    }
}

gslpp::complex GMcache::g_func(const double x) const{
    if(x<1) {
    gslpp::complex z = -gslpp::complex::i()*M_PI;
    gslpp::complex gs1 = sqrt(1.0-x)*(log((1.0+sqrt(1.0-x))/(1.0-sqrt(1.0-x)))+z)/2.0;
    return gs1;
    }
    else {
        gslpp::complex gg1 = sqrt(x-1.0)*asin(sqrt(1.0/x));
        return gg1;
    }
}

gslpp::complex GMcache::Int1(const double tau, const double lambda) const{
    return tau*lambda/(tau-lambda)/2.0+tau*tau*lambda*lambda/((tau-lambda)
           *(tau-lambda))/2.0*(f_func(tau)-f_func(lambda))+tau*tau*lambda/((tau-lambda)
           *(tau-lambda))*(g_func(tau)-g_func(lambda));
}

gslpp::complex GMcache::Int2(const double tau, const double lambda) const{
    return -tau*lambda/(tau-lambda)/2.0*(f_func(tau)-f_func(lambda));
}

void GMcache::computeSignalStrengthQuantities()
{
    double Mt = myGM->getQuarks(QCD::TOP).getMass();
    double Mb = myGM->getQuarks(QCD::BOTTOM).getMass();
    double Mc = myGM->getQuarks(QCD::CHARM).getMass();
    double Ms = myGM->getQuarks(QCD::STRANGE).getMass();
    double Mu = myGM->getQuarks(QCD::UP).getMass();
    double Md = myGM->getQuarks(QCD::DOWN).getMass();
    double Mtau = myGM->getLeptons(StandardModel::TAU).getMass();
    double Mmu = myGM->getLeptons(StandardModel::MU).getMass();
    double Me = myGM->getLeptons(StandardModel::ELECTRON).getMass();
    double MW = MWGM(myGM->Mw_tree());
    double cW2 = cW2GM(myGM->c02());

    double BrSM_htobb = myGM->computeBrHtobb();
    double BrSM_htotautau = myGM->computeBrHtotautau();
    double BrSM_htogaga = myGM->computeBrHtogaga();
    double BrSM_htoWW = myGM->computeBrHtoWW();
    double BrSM_htoZZ = myGM->computeBrHtoZZ();
    double BrSM_htogg = myGM->computeBrHtogg();
    double BrSM_htoZga = myGM->computeBrHtoZga();
    double BrSM_htocc = myGM->computeBrHtocc();

//    depending on which Higgs is heavier

    double gh_H3=0., gh_H5=0.;

    if(mH1sq>=mHl2)
    {
        rh_ff = cosa*cosa/(sinb*sinb);
        rh_VV = cosa*cosa*sinb*sinb - sqrt(32.0/3.0)*cosa*sina*cosb*sinb + 8.0/3.0*sina*sina*cosb*cosb;
        gh_H3 = ghH3pH3m;
        gh_H5 = ghH5pH5m;
    }
    else
    {
        rh_ff = sina*sina/(sinb*sinb);
        rh_VV = sina*sina*sinb*sinb + sqrt(32.0/3.0)*sina*cosa*cosb*sinb + 8.0/3.0*cosa*cosa*cosb*cosb;
        gh_H3 = gHH3pH3m;
        gh_H5 = gHH5pH5m;
    }

    //rh_gaga formula = abs(I_h_F+I_h_W+I_h_Hp)^2 / abs(I_hSM_F+I_hSM_W)^2

    gslpp::complex I_hSM_F = I_h_U(mHl2,Mu,Mc,Mt)+I_h_D(mHl2,Md,Ms,Mb)+I_h_L(mHl2,Me,Mmu,Mtau);
    gslpp::complex I_hSM_W = I_H_W(mHl2,MW);
    gslpp::complex I_h_F = sqrt(rh_ff)*I_hSM_F;
    gslpp::complex I_h_W = sqrt(rh_VV)*I_hSM_W;
    gslpp::complex I_h_Hp = 0.5*vev*( gh_H3*I_H_Hp(mAsq,mHl2)/mAsq + 5.0*gh_H5*I_H_Hp(mH5sq,mHl2)/mH5sq );

    double ABSgagaGM=(I_h_F+I_h_W+I_h_Hp).abs2();
    double ABSgagaSM=(I_hSM_F+I_hSM_W).abs2();
    rh_gaga=ABSgagaGM/ABSgagaSM;

    //rh_Zga formula = abs(A_h_F+A_h_W+A_h_Hp)^2 / abs(A_hSM_F+A_hSM_W)^2

    gslpp::complex A_hSM_F = A_h_U(mHl2,cW2,Mu,Mc,Mt,MZ)+A_h_D(mHl2,cW2,Md,Ms,Mb,MZ)+A_h_L(mHl2,cW2,Me,Mmu,Mtau,MZ);
    gslpp::complex A_hSM_W = A_H_W(mHl2,cW2,MW,MZ);
    gslpp::complex A_h_F = sqrt(rh_ff)*A_hSM_F;
    gslpp::complex A_h_W = sqrt(rh_VV)*A_hSM_W;
    gslpp::complex A_h_Hp = 0.5*vev*( gh_H3*A_H_Hp(mAsq,mHl2,cW2,MZ)/mAsq + 5.0*gh_H5*A_H_Hp(mH5sq,mHl2,cW2,MZ)/mH5sq );

    double ABSZgaGM=(A_h_F+A_h_W+A_h_Hp).abs2();
    double ABSZgaSM=(A_hSM_F+A_hSM_W).abs2();
    rh_Zga=ABSZgaGM/ABSZgaSM;

    rh_gg=rh_ff;

    sumModBRs = rh_ff*(BrSM_htobb+BrSM_htotautau+BrSM_htocc) + rh_VV*(BrSM_htoWW+BrSM_htoZZ) 
                + rh_gaga*BrSM_htogaga + rh_gg*BrSM_htogg + rh_Zga*BrSM_htoZga;

    Gamma_h = sumModBRs*myGM->computeGammaHTotal();

//    GM_BR_h_bb = rh_ff*BrSM_htobb/sumModBRs;
//    GM_BR_h_gaga = rh_gaga*BrSM_htogaga/sumModBRs;
//    GM_BR_h_tautau = rh_ff*BrSM_htotautau/sumModBRs;
//    GM_BR_h_WW = rh_VV*BrSM_htoWW/sumModBRs;
//    GM_BR_h_ZZ = rh_VV*BrSM_htoZZ/sumModBRs;
}

//void GMcache::runGMparameters()
//{
////    vev=myGM->v();
////    double cosb=myGM->getcosb();
//
//    std::string RGEorder=myGM->getRGEorderflag();
//    //flag will be used to transport information about model and RGEorder to the Runner:
//    //flag=0 for LO, 1 for approxNLO (and 2 for NLO - not implemented yet)
//    int flag;
//    if( RGEorder == "LO" ) flag=0;
//    else if( RGEorder == "approxNLO" ) flag=1;
////    else if( RGEorder == "NLO" ) flag=2;
//    else {
//        throw std::runtime_error("RGEorder can be only any of \"LO\", \"approxNLO\" or \"NLO\"");
//    }
//
//    double lambda1_at_MZ=lambda1;
//    double lambda2_at_MZ=lambda2;
//    double lambda3_at_MZ=lambda3;
//    double lambda4_at_MZ=lambda4;
//    double mu1_at_MZ=mu1;
//    double mu3_at_MZ=mu3;
//    double mu4_at_MZ=mu4;
//    double nu1_at_MZ=nu1;
//    double omega1_at_MZ=omega1;
//    double kappa1_at_MZ=kappa1;
//    double nu2_at_MZ=nu2;
//    double omega2_at_MZ=omega2;
//    double kappa2_at_MZ=kappa2;
//    double nu4_at_MZ=nu4;
//    double omega4_at_MZ=omega4;
//    double NLOuniscale=myGM->getNLOuniscaleGM();
//
//    if(fabs(Q_GM-log10(MZ))<0.005)   //at MZ scale
//    {
//        Q_cutoff=log10(MZ);
//
//        lambda1_at_Q = lambda1_at_MZ;
//        lambda2_at_Q = lambda2_at_MZ;
//        lambda3_at_Q = lambda3_at_MZ;
//        lambda4_at_Q = lambda4_at_MZ;
//        mu1_at_Q = mu1_at_MZ;
//        mu3_at_Q = mu3_at_MZ;
//        mu4_at_Q = mu4_at_MZ;
//        nu1_at_Q = nu1_at_MZ;
//        omega1_at_Q = omega1_at_MZ;
//        kappa1_at_Q = kappa1_at_MZ;
//        nu2_at_Q = nu2_at_MZ;
//        omega2_at_Q = omega2_at_MZ;
//        kappa2_at_Q = kappa2_at_MZ;
//        nu4_at_Q = nu4_at_MZ;
//        omega4_at_Q = omega4_at_MZ;
//    }
//    else   //at some other scale
//    {
//        double InitVals[15];
//        InitVals[0]=lambda1_at_MZ;
//        InitVals[1]=lambda2_at_MZ;
//        InitVals[2]=lambda3_at_MZ;
//        InitVals[3]=lambda4_at_MZ;
//        InitVals[4]=mu1_at_MZ;
//        InitVals[5]=mu3_at_MZ;
//        InitVals[6]=mu4_at_MZ;
//        InitVals[7]=nu1_at_MZ;
//        InitVals[8]=omega1_at_MZ;
//        InitVals[9]=kappa1_at_MZ;
//        InitVals[10]=nu2_at_MZ;
//        InitVals[11]=omega2_at_MZ;
//        InitVals[12]=kappa2_at_MZ;
//        InitVals[13]=nu4_at_MZ;
//        InitVals[14]=omega4_at_MZ;
//
//        Q_cutoff=myRunnerGM->RGERunnerGM(InitVals, 15, log10(MZ), Q_GM, flag, RpepsGM, NLOuniscale);  //Running up to Q_cutoff<=Q_GM
//
//        lambda1_at_Q = InitVals[0];
//        lambda2_at_Q = InitVals[1];
//        lambda3_at_Q = InitVals[2];
//        lambda4_at_Q = InitVals[3];
//        mu1_at_Q=InitVals[4];
//        mu3_at_Q=InitVals[5];
//        mu4_at_Q = InitVals[6];
//        nu1_at_Q = InitVals[7];
//        omega1_at_Q = InitVals[8];
//        kappa1_at_Q = InitVals[9];
//        nu2_at_Q = InitVals[10];
//        omega2_at_Q = InitVals[11];
//        kappa2_at_Q = InitVals[12];
//        nu4_at_Q = InitVals[13];
//        omega4_at_Q = InitVals[14];
//    }
//
//}

void GMcache::computeUnitarity()
{
//    std::string ModelType=myGM->getModelTypeGMflag();
//    if( ModelType != "custodial1" )
//    {
//        throw std::runtime_error("GM unitarity constraints are only implemented for the \"custodial1\" model.");
//    }
//
    double pi=M_PI;
//    /*
//    *******   LO part   *************
//    */

    unitarityeigenvalues.assign(0, (fabs(6.0*lambda1+11.0*lambda2+7.0*lambda3)+sqrt((6.0*lambda1-11.0*lambda2-7.0*lambda3)*(6.0*lambda1-11.0*lambda2-7.0*lambda3)+36.0*lambda4*lambda4))/(4.0*pi));
    unitarityeigenvalues.assign(1, fabs(lambda4-lambda5)/(2.0*pi));
    unitarityeigenvalues.assign(2, fabs(lambda2+2.0*lambda3)/pi);
    unitarityeigenvalues.assign(3, (fabs(2.0*lambda1+2.0*lambda2-lambda3)+sqrt((2.0*lambda1-2.0*lambda2+lambda3)*(2.0*lambda1-2.0*lambda2+lambda3)+lambda5*lambda5))/(4.0*pi));

//    
//    
//    /*
//    *******   NLO part   *************
//    */
//
//    double blambda1=(12.0*lambda1*lambda1 + 4.0*lambda3*lambda3 + 4.0*lambda3*lambda4 + 4.0*lambda4*lambda4 
//                     + 8.0*nu1*nu1 + 8.0*nu1*nu2 + 8.0*nu2*nu2)/(16.0*pi*pi);
//    double blambda2=(12.0*lambda2*lambda2 + 4.0*lambda3*lambda3 + 4.0*lambda3*lambda4 + 4.0*lambda4*lambda4
//                     + 8.0*omega1*omega1 + 8.0*omega1*omega2 + 8.0*omega2*omega2)/(16.0*pi*pi);
//    double blambda3=(4.0*lambda3*lambda3 + 4.0*lambda4*lambda4 + (lambda1+lambda2)*(6.0*lambda3+2.0*lambda4) 
//                     + 8.0*kappa2*kappa2 + 8.0*nu1*omega1 + 4.0*nu2*omega1 + 4.0*nu1*omega2)/(16.0*pi*pi);
//    double blambda4=(lambda1*lambda4 + lambda2*lambda4 + 4.0*lambda3*lambda4 + 6.0*lambda4*lambda4
//                     + 4.0*kappa1*kappa1 + 4.0*kappa1*kappa2 + 2.0*kappa2*kappa2 + 2.0*nu2*omega2)/(8.0*pi*pi);
//    double bmu1=(11.0*mu1*mu1 + 3.0*mu1*mu4 + mu1*(2.0*mu1+6.0*mu3+3.0*mu4)
//               + 3.0*nu4*nu4 + 3.0*omega4*omega4)/(16.0*pi*pi);
//    double bmu3=(18.0*kappa1*kappa1 + 18.0*kappa1*kappa2 + 134.0*mu1*mu1 + 6.0*mu1*(39.0*mu3 + 22.0*mu4)
//               + 3.0*(30.0*mu3*mu3 + 39.0*mu3*mu4 + 9.0*mu4*mu4 
//                      + 3.0*nu1*nu1 + 3.0*nu1*nu2 - 5.0*nu4*nu4
//                      + 3.0*omega1*omega1 + 3.0*omega1*omega2 - 5.0*omega4*omega4))/(72.0*pi*pi);
//    double bmu4=(18.0*kappa2*kappa2 + 4.0*mu1*mu1 + 156.0*mu1*mu4 + 54.0*mu3*mu4 + 144.0*mu4*mu4
//               + 9.0*nu2*nu2 + 6.0*nu4*nu4 + 9.0*omega2*omega2 + 6.0*omega4*omega4)/(144.0*pi*pi);
//    double bnu1=(6.0*kappa1*kappa1 + 6.0*kappa2*kappa2 + 18.0*lambda1*nu1
//               + 78.0*mu1*nu1 + 51.0*mu3*nu1 + 39.0*mu4*nu1 + 6.0*nu1*nu1
//               + 6.0*lambda1*nu2 + 32.0*mu1*nu2 + 24.0*mu3*nu2 + 6.0*mu4*nu2
//               + 6.0*nu2*nu2 + 10.0*nu4*nu4
//               + 12.0*lambda3*omega1 + 6.0*lambda4*omega1 + 6.0*lambda3*omega2)/(48.0*pi*pi);
//    double bomega1=(6.0*kappa1*kappa1 + 6.0*kappa2*kappa2 
//               + 12.0*lambda3*nu1 + 6.0*lambda4*nu1 + 6.0*lambda3*nu2
//               + 18.0*lambda2*omega1 + 78.0*mu1*omega1 + 51.0*mu3*omega1 + 39.0*mu4*omega1 + 6.0*omega1*omega1
//               + 6.0*lambda2*omega2 + 32.0*mu1*omega2 + 24.0*mu3*omega2 + 6.0*mu4*omega2 + 6.0*omega2*omega2
//               + 10.0*omega4*omega4)/(48.0*pi*pi);
//    double bkappa1=(6.0*kappa1*(2.0*lambda3 + 10.0*lambda4 + 18.0*mu1 + 17.0*mu3 + 13.0*mu4 + 2.0*nu1 + 2.0*omega1)
//               + kappa2*(24.0*lambda4 + 64.0*mu1 + 48.0*mu3 + 24.0*mu4 + 9.0*nu2 + 9.0*omega2)
//               + 20.0*nu4*omega4)/(96.0*pi*pi);
//    double bnu2=(4.0*kappa1*kappa2 + 6.0*kappa2*kappa2 + 2.0*lambda1*nu2 + ((14.0*mu1)/3.0 + mu3 + 9.0*mu4)*nu2 
//                + 4.0*nu1*nu2 + 6.0*nu2*nu2 + (25.0*nu4*nu4)/3.0 + 2.0*lambda4*omega2)/(16.0*pi*pi);
//    double bomega2=(4.0*kappa1*kappa2 + 6.0*kappa2*kappa2 + 2.0*lambda4*nu2 + 2.0*lambda2*omega2 
//                + ((14.0*mu1)/3.0 + mu3 + 9.0*mu4)*omega2 + 4.0*omega1*omega2 + 6.0*omega2*omega2 
//                + (25.0*omega4*omega4)/3.0)/(16.0*pi*pi);
//    double bkappa2=(kappa2*(6.0*lambda3 + 6.0*lambda4 + 14.0*mu1 + 3.0*mu3 + 27.0*mu4
//                     + 6.0*nu1 + 12.0*nu2 + 6.0*omega1 + 12.0*omega2)
//                + 6.0*kappa1*(nu2 + omega2) + 42.0*nu4*omega4)/(48.0*pi*pi);
//    double bnu4=(11.0*mu1*nu4 + 3.0*mu3*nu4 + 9.0*mu4*nu4 + 3.0*nu1*nu4 + 9.0*nu2*nu4 
//                + 3.0*kappa1*omega4 + 9.0*kappa2*omega4)/(16.0*pi*pi);
//    double bomega4=(3.0*kappa1*nu4 + 9.0*kappa2*nu4 
//                + (11.0*mu1 + 3.0*(mu3 + 3.0*mu4 + omega1 + 3.0*omega2))*omega4)/(16.0*pi*pi);
//
//    Sbmatrix1.assign(0,0, 3.0*blambda1/(16.0*pi));
//    Sbmatrix1.assign(0,1, (2.0*blambda3+blambda4)/(16.0*pi));
//    Sbmatrix1.assign(1,0, Sbmatrix1(0,1));
//    Sbmatrix1.assign(0,3, (2.0*bnu1+bnu2)/(8.0*sqrt(2.0)*pi));
//    Sbmatrix1.assign(3,0, Sbmatrix1(0,3));
//    Sbmatrix1.assign(1,1, 3.0*blambda2/(16.0*pi));
//    Sbmatrix1.assign(1,3, (2.0*bomega1+bomega2)/(8.0*sqrt(2.0)*pi));
//    Sbmatrix1.assign(3,1, Sbmatrix1(1,3));
//    Sbmatrix1.assign(2,2, (blambda3+5.0*blambda4)/(16.0*pi));
//    Sbmatrix1.assign(2,3, (4.0*bkappa1+2.0*bkappa2)/(16.0*pi));
//    Sbmatrix1.assign(3,2, Sbmatrix1(2,3));
//    Sbmatrix1.assign(3,3, (26.0*bmu1+17.0*bmu3+13.0*bmu4)/(32.0*pi));
//
//    Sbmatrix2.assign(0,0, blambda1/(16.0*pi));
//    Sbmatrix2.assign(0,1, blambda4/(16.0*pi));
//    Sbmatrix2.assign(1,0, Sbmatrix2(0,1));
//    Sbmatrix2.assign(0,3, bnu2/(8.0*sqrt(2.0)*pi));
//    Sbmatrix2.assign(3,0, Sbmatrix2(0,3));
//    Sbmatrix2.assign(1,1, blambda2/(16.0*pi));
//    Sbmatrix2.assign(1,3, bomega2/(8.0*sqrt(2.0)*pi));
//    Sbmatrix2.assign(3,1, Sbmatrix2(1,3));
//    Sbmatrix2.assign(2,2, (blambda3+blambda4)/(16.0*pi));
//    Sbmatrix2.assign(2,3, bkappa2/(8.0*pi));
//    Sbmatrix2.assign(3,2, Sbmatrix2(2,3));
//    Sbmatrix2.assign(3,3, (14.0*bmu1+3.0*bmu3+27.0*bmu4)/(96.0*pi));
//
//    Seigenvectors1T=Seigenvectors1.hconjugate();
//    Seigenvectors2T=Seigenvectors2.hconjugate();
//
//    for (int i=0; i < 4; i++) {
//        for (int k=0; k < 4; k++) {
//            for (int l=0; l < 4; l++) {
//                Sbeigenvalues1.assign(i, Sbeigenvalues1(i) + Seigenvectors1T(i,k) * Sbmatrix1(k,l) * Seigenvectors1(l,i) );
//                Sbeigenvalues2.assign(i, Sbeigenvalues2(i) + Seigenvectors2T(i,k) * Sbmatrix2(k,l) * Seigenvectors2(l,i) );
//            }                
//        }
//        betaeigenvalues.assign(i, -1.5 * Sbeigenvalues1(i));
//        betaeigenvalues.assign(i+4, -1.5 * Sbeigenvalues2(i));
//    }
//
//    betaeigenvalues.assign(8, -1.5 * (blambda3-blambda4)/(16.0*pi));
//    betaeigenvalues.assign(9, -1.5 * sqrt(15.0)*bnu4/(16.0*pi));
//    betaeigenvalues.assign(10, -1.5 * sqrt(15.0)*bomega4/(16.0*pi));
//
//    for (int i=0; i < 11; i++) {
//        NLOunitarityeigenvalues.assign(i, -(gslpp::complex::i()-1.0/pi)*unitarityeigenvalues(i)*unitarityeigenvalues(i) + betaeigenvalues(i) );
//    }
}

double GMcache::cW2GM(const double c02) const{
    return c02;
}



double GMcache::MWGM(const double MW) const{
    return MW;
}

double GMcache::updateCache()
{
//    GMmodel=myGM->getModelTypeGMflag();
    Q_GM=myGM->getQ_GM();
    vev=myGM->v();
    mHl=myGM->getMHl();
    mHl2=mHl*mHl;
    MZ=myGM->getMz();
    tanb=myGM->gettanb();
    sinb=myGM->getsinb();
    cosb=myGM->getcosb();
    sina=myGM->getsina();
    cosa=myGM->getcosa();
    mH1sq=myGM->getinputmHh2();
    mAsq=myGM->getmAsq();
    mH5sq=myGM->getmH5sq();
    Mu1=myGM->getMu1();
    Mu2=myGM->getMu2();
    M1sq=-vev*Mu1/(sqrt(2.0)*cosb);
    M2sq=-sqrt(18.0)*cosb*vev*Mu2;

    vPhi=vev*sinb;
    vDelta=vev*cosb/sqrt(8.0);
    double mHlsq, mHhsq;
    if(mH1sq>=mHl2)
    {
        mHlsq=mHl2;
        mHhsq=mH1sq;
    }
    else
    {
        mHlsq=mH1sq;
        mHhsq=mHl2;
    }
    double cos2a=cosa*cosa-sina*sina;
    double cos2b=cosb*cosb-sinb*sinb;

    lambda1 = (mHlsq+mHhsq+(mHlsq-mHhsq)*cos2a)/(16.0*vev*vev*sinb*sinb);
    lambda2 = (-M1sq+M2sq+mAsq-2.0/3.0*mH5sq+(M1sq-mAsq)*cos2b
               +(mHlsq+mHhsq+(mHhsq-mHlsq)*cos2a)/3.0)/(2.0*vev*vev*cosb*cosb);
    lambda3 = (mH5sq-M2sq+(2.0*M1sq-3.0*mAsq)*sinb*sinb)/(vev*vev*cosb*cosb);
    lambda4 = (mAsq-0.5*M1sq+(mHhsq-mHlsq)*sina*cosa/(sqrt(6.0)*sinb*cosb))/(vev*vev);
    lambda5 = 2.0*(M1sq-mAsq)/(vev*vev);

    //triple scalar couplings
    ghhh = -(2.0*sina*sina*sinb*(3.0*cosa+sqrt(6.0)*sina*tanb))/vev * M1sq
           +(sqrt(2.0/3.0)*sina*sina*sina)/(vev*cosb) * M2sq
           +(-3.0*cosa*cosa*cosa/sinb+2.0*sqrt(6.0)*sina*sina*sina/cosb)/vev * mHlsq;
    ghhH = sina*sinb*(4.0-6.0*sina*sina+2.0*sqrt(6.0)*sina*cosa*tanb)/vev * M1sq
            -(sqrt(2.0/3.0)*cosa*sina*sina)/(vev*cosb) * M2sq
            -sina*(cosa*cosa/sinb+2.0*sqrt(2.0/3.0)*sina*cosa/cosb)/vev * (2.0*mHlsq+mHhsq);
    ghHH = cosa*sinb*(4.0-6.0*cosa*cosa-2.0*sqrt(6.0)*sina*cosa*tanb)/vev * M1sq
            +(sqrt(2.0/3.0)*cosa*cosa*sina)/(vev*cosb) * M2sq
            +cosa*(-sina*sina/sinb+2.0*sqrt(2.0/3.0)*sina*cosa/cosb)/vev * (mHlsq+2.0*mHhsq);
    gHHH = (2.0*cosa*cosa*sinb*(-3.0*sina+sqrt(6.0)*cosa*tanb))/vev * M1sq
           -(sqrt(2.0/3.0)*cosa*cosa*cosa)/(vev*cosb) * M2sq
           -(3.0*sina*sina*sina/sinb+2.0*sqrt(6.0)*cosa*cosa*cosa/cosb)/vev * mHhsq;
    ghH3H3 = -(2.0*sqrt(2.0/3.0)*sina)/(vev*cosb) * M1sq
             +(4.0*sqrt(2.0/3.0)*sina*cosb-2.0*cosa*sinb)/vev * mAsq
             +(2.0*sqrt(2.0/3.0)*sina*sinb*tanb-cosa*cosb/tanb)/vev * mHlsq;
    gHH3H3 = (2.0*sqrt(2.0/3.0)*cosa)/(vev*cosb) * M1sq
             -(4.0*sqrt(2.0/3.0)*cosa*cosb+2.0*sina*sinb)/vev * mAsq
             -(2.0*sqrt(2.0/3.0)*cosa*sinb*tanb+sina*cosb/tanb)/vev * mHhsq;
    ghH3pH3m = ghH3H3;
    gHH3pH3m = gHH3H3;
    ghH5H5 = 2.0*sinb*(2.0*cosa+sqrt(6.0)*sina*tanb)/vev * (M1sq-2.0*mAsq)
             +(2.0*sqrt(2.0/3.0)*sina)/(vev*cosb) * (mHlsq + 2.0*mH5sq - M2sq)
             +(2.0*cosa*sinb)/vev * mAsq;
    gHH5H5 = 2.0*sinb*(2.0*sina-sqrt(6.0)*cosa*tanb)/vev * (M1sq-2.0*mAsq)
             -(2.0*sqrt(2.0/3.0)*cosa)/(vev*cosb) * (mHhsq + 2.0*mH5sq - M2sq)
             +(2.0*sina*sinb)/vev * mAsq;
    ghH5pH5m = ghH5H5;
    gHH5pH5m = gHH5H5;
    ghH5ppH5mm = ghH5H5;
    gHH5ppH5mm = gHH5H5;
    gH3H3H5 = 2.0*((mH5sq-2.0*mAsq)*cosb-(2.0*M1sq-3.0*mAsq+mH5sq)/cosb)/(sqrt(3.0)*vev);
    gH5H5H5 = 2.0*((3.0*mH5sq-4.0*M2sq)/cosb+(6.0*M1sq-9.0*mAsq)*sinb*tanb)/(sqrt(3.0)*vev);
    gH3H3pH5m = gslpp::complex::i()*(4.0*M1sq-4.0*mAsq+mH5sq+(2.0*mAsq-mH5sq)*(cosb*cosb-sinb*sinb))/(2.0*vev*cosb);
    gH5H3pH3m = (4.0*M1sq-4.0*mAsq+mH5sq+(2.0*mAsq-mH5sq)*(cosb*cosb-sinb*sinb))/(2.0*sqrt(3.0)*vev*cosb);
    gH5H5pH5m = 0.5*gH5H5H5;
    gH5H5ppH5mm = -gH5H5H5;
    gH5ppH3mH3m = sqrt(6.0)*gH5H3pH3m;
    gH5ppH5mH5m = gH5H5pH5m;

    //    runGMparameters();
    computeUnitarity();
    computeSignalStrengthQuantities();
    return mHl2;
}

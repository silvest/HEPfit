/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ScalarPotential.h"

SUSYScalarPotential::SUSYScalarPotential(const StandardModel& SM_i) : 
    mySUSY(static_cast<const SUSY&> (SM_i))

{}

gslpp::vector<double> SUSYScalarPotential::coefficients()
{
    gslpp::matrix<gslpp::complex> MsLhat2(3,3,0);
    gslpp::matrix<gslpp::complex> MsEhat2(3,3,0);
    MsLhat2 = mySUSY.getMsLhat2();
    MsEhat2 = mySUSY.getMsEhat2();
    double mstauLsq = MsLhat2(2,2).real();
    double msmuRsq = MsEhat2(1,1).real();
//	mhusq = ((mAsq/(1 + tb**2) + (self.mzsq/2)*((1 - tb**2)/(1 + tb**2)) - mu**2))
//	mhdsq = ((mAsq * tb**2/(1 + tb**2) - (self.mzsq/2)*((1 - tb**2)/(1 + tb**2)) - mu**2))
//	bmu = (mAsq) * tb/(1 + tb**2)
//    double mA = mySUSY.getMHa();
    double mA = 1000.;
    double tanb = mySUSY.getTanb();
    double MZ = mySUSY.getMz();
    gslpp::complex muH = mySUSY.getMuH();
    double musq = muH.abs2();
    double mHdsq = mA*mA*tanb*tanb/(1.0+tanb*tanb)-0.5*MZ*MZ*(1.0-tanb*tanb)/(1.0+tanb*tanb)-musq;
    double vd = mySUSY.v1();
    double Mw = mySUSY.Mw_tree();
    double sw2 = mySUSY.StandardModel::sW2(Mw);
    double ttw2 = sw2/(1.0 - sw2);
    double g2sq = 8. * mySUSY.getGF() / sqrt(2.0) * Mw * Mw;
    double g1sq = g2sq*ttw2;
    double mMU = mySUSY.getLeptons(StandardModel::MU).getMass();
    double mTAU = mySUSY.getLeptons(StandardModel::TAU).getMass();
    gslpp::matrix<gslpp::complex> TEhat = mySUSY.getTEhat();
    double Al23 = TEhat(1,2).abs();

    std::cout << "tanb = " << tanb << std::endl;
    std::cout << "MZ = " << MZ << std::endl;
    std::cout << "mHdsq = " << mHdsq << std::endl;
    std::cout << "g1 = " << sqrt(g1sq) << std::endl;
    std::cout << "g2 = " << sqrt(g2sq) << std::endl;
    std::cout << "Al23 = " << Al23 << std::endl;
    std::cout << "vd = " << vd << std::endl;
    std::cout << "ytau = " << mMU/vd*sqrt(2.) << std::endl;
    std::cout << "ymu = " << mTAU/vd*sqrt(2.) << std::endl;


    gslpp::vector<double> coefficients(35, 0.);

//          Let's call the scalar fields x1, x2 and x3.
//          The coefficients will we the ones corresponding to:
//          1,
//          x1^1, x2^1, x3^1,
//          x1^2, x1*x2, x1*x3, x2^2, x2*x3, x3^2,
//          x1^3, x1^2*x2, x1^2*x3, x1*x2^2, x1*x2*x3, x1*x3^2, x2^3, x2^2*x3, x2*x3^2, x3^3,
//          x1^4, x1^3*x2, x1^3*x3, x1^2*x2^2, x1^2*x2*x3, x1^2*x3^2, x1*x2^3, x1*x2^2*x3, x1*x2*x3^2, x1*x3^3, x2^4, x2^3*x3, x2^2*x3^2, x2*x3^3, x3^4

//          Let's take x1=Hd, x2=stauL, x3=smuR
    coefficients(0)=(mHdsq+musq+0.125*(g1sq+g2sq)*vd*vd)*vd*vd; //x1^0=x2^0=x3^0
    coefficients(1)=-vd*(2.0*(mHdsq+musq)+0.5*vd*vd*(g1sq+g2sq)); //x1^1
//    double x21 = 0.0;
//    double x31 = 0.0;
    coefficients(4)=mHdsq+musq+0.75*vd*vd*(g1sq+g2sq); //x1^2
//    double x11x21 = 0.0;
//    double x11x31 = 0.0;
    double x22 = mstauLsq+mTAU*mTAU+0.25*vd*vd*(g1sq-g2sq);
    coefficients(7)=x22;
    double x21x31 = 2.0*Al23*vd;
    coefficients(8)=x21x31;
    double x32 = msmuRsq+mMU*mMU-0.5*vd*vd*g1sq;
    coefficients(9)=x32;
    double x13 = -0.5*vd*(g1sq+g2sq);
    coefficients(10)=x13;
//    double x12x21 = 0.0;
//    double x12x31 = 0.0;
    double x11x22 = -2.0*mTAU*mTAU/vd-0.5*vd*(g1sq-g2sq);
    coefficients(13)=x11x22;
    double x11x21x31 = -2.0*Al23;
    coefficients(14)=x11x21x31;
    double x11x32 = -2.0*mMU*mMU/vd+vd*g1sq;
    coefficients(15)=x11x32;
//    double x23 = 0.0;
//    double x22x31 = 0.0;
//    double x21x32 = 0.0;
//    double x33 = 0.0;
    double x14 = 0.125*(g1sq+g2sq);
    coefficients(20)=x14;
//    double x13x21 = 0.0;
//    double x13x31 = 0.0;
    double x12x22 = mTAU*mTAU/(vd*vd)+0.25*(g1sq-g2sq);
    coefficients(23)=x12x22;
//    double x12x21x31 = 0.0;
    double x12x32 = mMU*mMU/(vd*vd)-0.5*g1sq;
    coefficients(25)=x12x32;
//    double x11x23 = 0.0;
//    double x11x22x31 = 0.0;
//    double x11x21x32 = 0.0;
//    double x11x33 = 0.0;
    double x24 = 0.125*(g1sq+g2sq);
    coefficients(30)=x24;
//    double x23x31 = 0.0;
    double x22x32 = -0.5*g1sq;
    coefficients(32)=x22x32;
//    double x21x33 = 0.0;
    double x34 = 0.5*g1sq;
    coefficients(34)=x34;

    return coefficients;
}

double SUSYScalarPotential::potential(gslpp::vector<double> coefficients, double field1, double field2, double field3)
{
    gslpp::vector<double> a=coefficients;
    double x1=field1;
    double x2=field2;
    double x3=field3;
    double pot=a(0)
               +a(1)*x1 +a(2)*x2 +a(3)*x3
               +a(4)*x1*x1 +a(5)*x1*x2 +a(6)*x1*x3 +a(7)*x2*x2 +a(8)*x2*x3 +a(9)*x3*x3
               +a(10)*x1*x1*x1 +a(11)*x1*x1*x2 +a(12)*x1*x1*x3 +a(13)*x1*x2*x2 +a(14)*x1*x2*x3 +a(15)*x1*x3*x3
                               +a(16)*x2*x2*x2 +a(17)*x2*x2*x3 +a(18)*x2*x3*x3 +a(19)*x3*x3*x3
               +a(20)*x1*x1*x1*x1 +a(21)*x1*x1*x1*x2 +a(22)*x1*x1*x1*x3 +a(23)*x1*x1*x2*x2 +a(24)*x1*x1*x2*x3
                                  +a(25)*x1*x1*x3*x3 +a(26)*x1*x2*x2*x2 +a(27)*x1*x2*x2*x3 +a(28)*x1*x2*x3*x3
                                  +a(29)*x1*x3*x3*x3 +a(30)*x2*x2*x2*x2 +a(31)*x2*x2*x2*x3 +a(32)*x2*x2*x3*x3
                                  +a(33)*x2*x3*x3*x3 +a(34)*x3*x3*x3*x3;
//    double pot=0.0017924;
    return pot;
}

gslpp::vector<double> SUSYScalarPotential::potentialderivative(gslpp::vector<double> coefficients, double field1, double field2, double field3)
{
    gslpp::vector<double> dV(3, 0.);
    gslpp::vector<double> a=coefficients;
    double x1=field1;
    double x2=field2;
    double x3=field3;
    std::cout << "coefficients = " << coefficients << std::endl;
    std::cout << "x1 = " << x1 << std::endl;
    std::cout << "x2 = " << x2 << std::endl;
    std::cout << "x3 = " << x3 << std::endl;

    //derivative wrt x1
    dV(0)=a(1) + 2.0*a(4)*x1 + 3.0*a(10)*x1*x1 + 4.0*a(20)*x1*x1*x1 
          + a(5)*x2 + 2.0*a(11)*x1*x2 + 3.0*a(21)*x1*x1*x2 + a(13)*x2*x2 + 2.0*a(23)*x1*x2*x2 + a(26)*x2*x2*x2
          + a(6)*x3 + 2.0*a(12)*x1*x3 + 3.0*a(22)*x1*x1*x3 + a(14)*x2*x3 + 2.0*a(24)*x1*x2*x3
          + a(27)*x2*x2*x3 + a(15)*x3*x3 + 2.0*a(25)*x1*x3*x3 + a(28)*x2*x3*x3 + a(29)*x3*x3*x3;
    //derivative wrt x2
    dV(1)=a(2) + a(5)*x1 + a(11)*x1*x1 + a(21)*x1*x1*x1
          + 2.0*a(7)*x2 + 2.0*a(13)*x1*x2 + 2.0*a(23)*x1*x1*x2 + 3.0*a(16)*x2*x2 + 3.0*a(26)*x1*x2*x2 + 4.0*a(30)*x2*x2*x2
          + a(8)*x3 + a(14)*x1*x3 + a(24)*x1*x1*x3 + 2.0*a(17)*x2*x3 + 2.0*a(27)*x1*x2*x3
          + 3.0*a(31)*x2*x2*x3 + a(18)*x3*x3 + a(28)*x1*x3*x3 + 2.0*a(32)*x2*x3*x3 + a(33)*x3*x3*x3;
    //derivative wrt x3
    dV(2)=a(3) + a(6)*x1 + a(12)*x1*x1 + a(22)*x1*x1*x1
          + a(8)*x2 + a(14)*x1*x2 + a(24)*x1*x1*x2 + a(17)*x2*x2 + a(27)*x1*x2*x2 + a(31)*x2*x2*x2
          + 2.0*a(9)*x3 + 2.0*a(15)*x1*x3 + 2.0*a(25)*x1*x1*x3 + 2.0*a(18)*x2*x3 + 2.0*a(28)*x1*x2*x3
          + 2.0*a(32)*x2*x2*x3 + 3.0*a(19)*x3*x3 + 3.0*a(29)*x1*x3*x3 + 3.0*a(33)*x2*x3*x3 + 4.0*a(34)*x3*x3*x3;
    return dV;
}

gslpp::vector<double> SUSYScalarPotential::secondpotentialderivative(gslpp::vector<double> coefficients, double field1, double field2, double field3)
{
    gslpp::vector<double> d2V(6, 0.);
    gslpp::vector<double> a=coefficients;
    double x1=field1;
    double x2=field2;
    double x3=field3;

    //derivative wrt x1*x1
    d2V(0)=2.0*a(4) + 6.0*a(10)*x1 + 12.0*a(20)*x1*x1 + 2.0*a(11)*x2 + 6.0*a(21)*x1*x2 
           + 2.0*a(23)*x2*x2 + 2.0*a(12)*x3 + 6.0*a(22)*x1*x3 + 2.0*a(24)*x2*x3 + 2.0*a(25)*x3*x3;
    //derivative wrt x1*x2
    d2V(1)=a(5) + 2.0*a(11)*x1 + 3.0*a(21)*x1*x1 + 2.0*a(13)*x2 + 4.0*a(23)*x1*x2 
           + 3.0*a(26)*x2*x2 + a(14)*x3 + 2.0*a(24)*x1*x3 + 2.0*a(27)*x2*x3 + a(28)*x3*x3;
    //derivative wrt x1*x3
    d2V(2)=a(6) + 2.0*a(12)*x1 + 3.0*a(22)*x1*x1 + a(14)*x2 + 2.0*a(24)*x1*x2 + a(27)*x2*x2 
           + 2.0*a(15)*x3 + 4.0*a(25)*x1*x3 + 2.0*a(28)*x2*x3 + 3.0*a(29)*x3*x3;
    //derivative wrt x2*x2
    d2V(3)=2.0*a(7) + 2.0*a(13)*x1 + 2.0*a(23)*x1*x1 + 6.0*a(16)*x2 + 6.0*a(26)*x1*x2 
           + 12.0*a(30)*x2*x2 + 2.0*a(17)*x3 + 2.0*a(27)*x1*x3 + 6.0*a(31)*x2*x3 + 2.0*a(32)*x3*x3;
    //derivative wrt x2*x3
    d2V(4)=a(8) + a(14)*x1 + a(24)*x1*x1 + 2.0*a(17)*x2 + 2.0*a(27)*x1*x2 + 3.0*a(31)*x2*x2 
           + 2.0*a(18)*x3 + 2.0*a(28)*x1*x3 + 4.0*a(32)*x2*x3 + 3.0*a(33)*x3*x3;
    //derivative wrt x3*x3
    d2V(5)=2.0*a(9) + 2.0*a(15)*x1 + 2.0*a(25)*x1*x1 + 2.0*a(18)*x2 + 2.0*a(28)*x1*x2 
           + 2.0*a(32)*x2*x2 + 6.0*a(19)*x3 + 6.0*a(29)*x1*x3 + 6.0*a(33)*x2*x3 + 12.0*a(34)*x3*x3;
    return d2V;
}

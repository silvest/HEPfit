/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdlib.h>
#include <iostream>
#include <gslpp_complex.h>
#include <gslpp_vector_double.h>
#include <gslpp_vector_complex.h>
#include <gslpp_matrix_double.h>
#include <gslpp_matrix_complex.h>

using namespace gslpp;

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "gslpptest test 1" << std::endl;
    // Tests of gslpp::complex

    std::cout << "Tests of gslpp::complex" << std::endl;
    complex z1(1.,1.,false);
    complex z2(2.,2.,false);
    complex z3(3., 3., true);
    complex z4(z2 / 3.);

    std::cout << "z1 = " << z1 << std::endl;
    std::cout << "z2 = " << z2 << std::endl;
    std::cout << "z3 = " << z3 << std::endl;
    z3.real() = 5.;
    z3.imag() = 7.;
    std::cout << "z3 = " << z3 << std::endl;
    std::cout << "z4 = " << z4 << std::endl;
    std::cout << "i^2 is real = " << (pow(complex::i(), 2.)).is_real() << std::endl;
    std::cout << "i^2 is imag = " << (pow(complex::i(), 2.)).is_imag() << std::endl;
    std::cout << "i is real = " << complex::i().is_real() << std::endl;
    std::cout << "i is imag = " << complex::i().is_imag() << std::endl;
    std::cout << "abs(z1) = " << z1.abs() << std::endl;
    std::cout << "arg(z1) = " << z1.arg() << std::endl;
    std::cout << "abs2(z1) = " << z1.abs2() << std::endl;
    std::cout << "log(abs(z1)) = " << z1.log_of_abs() << std::endl;
    complex z5 = z4;
    std::cout << "z5 = (z4) " << z5 << std::endl;
    std::cout << "z4 == z5 " << (z4 == z5) << std::endl;
    std::cout << "z4 != z5 " << (z4 != z5) << std::endl;
    std::cout << "z4 == z3 " << (z4 == z3) << std::endl;
    std::cout << "z4 != z3 " << (z4 != z3) << std::endl;
    z5 = 4.;
    std::cout << "z5 = (4) " << z5 << std::endl;
    z5 = -z4;
    std::cout << "z5 = (-z4) " << z5 << std::endl;
    std::cout << "1/z1 = " << z1.inverse() << std::endl;
    std::cout << "z1+z2 = " << z1 + z2 << std::endl;
    std::cout << "z1-z2 = " << z1 - z2 << std::endl;
    std::cout << "z1*z2 = " << z1 * z2 << std::endl;
    std::cout << "z1/z2 = " << z1 / z2 << std::endl;
    std::cout << "z1+2 = " << z1 + 2. << std::endl;
    std::cout << "z1-2 = " << z1 - 2. << std::endl;
    std::cout << "z1*2 = " << z1 * 2. << std::endl;
    std::cout << "z1/2 = " << z1 / 2. << std::endl;
    z5 = 2. * complex::i() + 3.;
    z5 += z1;
    std::cout << "z5(+=z1) = " << z5 << std::endl;
    z5 -= z1;
    std::cout << "z5(-=z1) = " << z5 << std::endl;
    z5 *= z1;
    std::cout << "z5(*=z1) = " << z5 << std::endl;
    z5 /= z1;
    std::cout << "z5(/=z1) = " << z5 << std::endl;
    z5 = 2. * complex::i() + 3.;
    z5 += 2.;
    std::cout << "z5(+=2.) = " << z5 << std::endl;
    z5 -= 2.;
    std::cout << "z5(-=2.) = " << z5 << std::endl;
    z5 *= 2.;
    std::cout << "z5(*=2.) = " << z5 << std::endl;
    z5 /= 2.;
    std::cout << "z5(/=2.) = " << z5 << std::endl;
    std::cout << "exp(z1) = " << exp(z1) << std::endl;
    std::cout << "log(z1) = " << log(z1) << std::endl;
    std::cout << "log10(z1) = " << log10(z1) << std::endl;
    std::cout << "log(z1,z2) = " << log(z1, z2) << std::endl;
    std::cout << "sqrt(z1) = " << sqrt(z1) << std::endl;
    std::cout << "pow(z1,z2) = " << pow(z1, z2) << std::endl;
    std::cout << "pow(z1,3.) = " << pow(z1, 3.) << std::endl;
    std::cout << "sin(z1) = " << sin(z1) << std::endl;
    std::cout << "cos(z1) = " << cos(z1) << std::endl;
    std::cout << "tan(z1) = " << tan(z1) << std::endl;
    std::cout << "sec(z1) = " << sec(z1) << std::endl;
    std::cout << "csc(z1) = " << csc(z1) << std::endl;
    std::cout << "cot(z1) = " << cot(z1) << std::endl;
    std::cout << "arcsin(z1) = " << arcsin(z1) << std::endl;
    std::cout << "arccos(z1) = " << arccos(z1) << std::endl;
    std::cout << "arctan(z1) = " << arctan(z1) << std::endl;
    std::cout << "arcsec(z1) = " << arcsec(z1) << std::endl;
    std::cout << "arccsc(z1) = " << arccsc(z1) << std::endl;
    std::cout << "arccot(z1) = " << arccot(z1) << std::endl;
    std::cout << "sinh(z1) = " << sinh(z1) << std::endl;
    std::cout << "cosh(z1) = " << cosh(z1) << std::endl;
    std::cout << "tanh(z1) = " << tanh(z1) << std::endl;
    std::cout << "sech(z1) = " << sech(z1) << std::endl;
    std::cout << "csch(z1) = " << csch(z1) << std::endl;
    std::cout << "coth(z1) = " << coth(z1) << std::endl;
    std::cout << "arcsinh(z1) = " << arcsinh(z1) << std::endl;
    std::cout << "arccosh(z1) = " << arccosh(z1) << std::endl;
    std::cout << "arctanh(z1) = " << arctanh(z1) << std::endl;
    std::cout << "arcsech(z1) = " << arcsech(z1) << std::endl;
    std::cout << "arccsch(z1) = " << arccsch(z1) << std::endl;
    std::cout << "arccoth(z1) = " << arccoth(z1) << std::endl;

    // Tests of gslpp::vector<double>
    std::cout << "\nTests of gslpp::vector<double>" << std::endl;
    vector<double> vd1(3, 3.);
    std::cout << "vd1(=(3,3,3)) = " << vd1 << std::endl;
    vector<double> vd2(vd1);
    std::cout << "vd2(=vd1) = " << vd2 << std::endl;
    std::cout << "size(vd2) = " << vd2.size() << std::endl;
    vector<double> vd3(3, 0.);
    std::cout << "vd3(=(0,0,0)) = " << vd3 << std::endl;
    vd3 = vd2;
    std::cout << "vd3(=vd2) = " << vd3 << std::endl;
    std::cout << "vd3(2) = " << vd3(2) << std::endl;
    std::cout << "|vd3| = " << vd3.mod() << std::endl;
    std::cout << "vd1+vd2 = " << vd1 + vd2 << std::endl;
    std::cout << "vd1-vd2 = " << vd1 - vd2 << std::endl;
    std::cout << "vd1*vd2 = " << vd1 * vd2 << std::endl;
    vd2 += vd1;
    std::cout << "vd2(+=vd1) = " << vd2 << std::endl;
    vd2 -= vd1;
    std::cout << "vd2(-=vd1) = " << vd2 << std::endl;
    std::cout << "vd1+2 = " << vd1 + 2. << std::endl;
    std::cout << "vd1-2 = " << vd1 - 2. << std::endl;
    std::cout << "vd1*2 = " << vd1 * 2. << std::endl;
    std::cout << "vd1/2 = " << vd1 / 2. << std::endl;
    std::cout << "2+vd2 = " << 2. + vd1 << std::endl;
    std::cout << "2-vd2 = " << 2. - vd1 << std::endl;
    std::cout << "2*vd2 = " << 2. * vd1 << std::endl;
    vd2 += 2.;
    std::cout << "vd2(+=2) = " << vd2 << std::endl;
    vd2 -= 2.;
    std::cout << "vd2(-=2) = " << vd2 << std::endl;
    vd2 *= 2.;
    std::cout << "vd2(*=2) = " << vd2 << std::endl;
    vd2 /= 2.;
    std::cout << "vd1+(2+i) = " << vd1 + (2. + complex::i()) << std::endl;
    std::cout << "vd1-(2+i) = " << vd1 - (2. + complex::i()) << std::endl;
    std::cout << "vd1*(2+i) = " << vd1 * (2. + complex::i()) << std::endl;
    std::cout << "vd1/(2+i) = " << vd1 / (2. + complex::i()) << std::endl;

    // Tests of gslpp::vector<complex>
    std::cout << "\nTests of gslpp::vector<complex>" << std::endl;
    vector<complex> vc1(3, 3. + complex::i());
    std::cout << "vc1(=(3+i,3+i,3+i)) = " << vc1 << std::endl;
    vector<complex> vc2(vc1);
    std::cout << "vc2(=vc1) = " << vc2 << std::endl;
    std::cout << "size(vc2) = " << vc2.size() << std::endl;
    vector<complex> vc3(3, 0.);
    std::cout << "vc3(=(0,0,0)) = " << vc3 << std::endl;
    vc3 = vc2;
    std::cout << "vc3(=vc2) = " << vc3 << std::endl;
    std::cout << "vc3(2) = " << vc3(2) << std::endl;
    vc3.assign(1, 6.);
    std::cout << "vc3 = " << vc3 << std::endl;
    std::cout << "|vc3| = " << vc3.mod() << std::endl;
    std::cout << "vc1+vc2 = " << vc1 + vc2 << std::endl;
    std::cout << "vc1-vc2 = " << vc1 - vc2 << std::endl;
    std::cout << "vc1*vc2 = " << vc1 * vc2 << std::endl;
    vc2 += vc1;
    std::cout << "vc2(+=vc1) = " << vc2 << std::endl;
    vc2 -= vc1;
    std::cout << "vc2(-=vc1) = " << vc2 << std::endl;
    std::cout << "vc1+(2-i) = " << vc1 + (2. - complex::i()) << std::endl;
    std::cout << "vc1-(2-i) = " << vc1 - (2. - complex::i()) << std::endl;
    std::cout << "vc1*(2-i) = " << vc1 * (2. - complex::i()) << std::endl;
    std::cout << "vc1/(2-i) = " << vc1 / (2. - complex::i()) << std::endl;
    std::cout << "(2-i)+vc2 = " << (2. - complex::i()) + vc1 << std::endl;
    std::cout << "(2-i)-vc2 = " << (2. - complex::i()) - vc1 << std::endl;
    std::cout << "(2-i)*vc2 = " << (2. - complex::i()) * vc1 << std::endl;
    vc2 += 2. - complex::i();
    std::cout << "vc2(+=2-i) = " << vc2 << std::endl;
    vc2 -= 2. - complex::i();
    std::cout << "vc2(-=2-i) = " << vc2 << std::endl;
    vc2 *= 2. - complex::i();
    std::cout << "vc2(*=2-i) = " << vc2 << std::endl;
    vc2 /= 2. - complex::i();
    std::cout << "vc2(/=2-i) = " << vc2 << std::endl;
    std::cout << "vc1+2 = " << vc1 + 2. << std::endl;
    std::cout << "vc1-2 = " << vc1 - 2. << std::endl;
    std::cout << "vc1*2 = " << vc1 * 2. << std::endl;
    std::cout << "vc1/2 = " << vc1 / 2. << std::endl;
    std::cout << "2+vc2 = " << 2. + vc1 << std::endl;
    std::cout << "2-vc2 = " << 2. - vc1 << std::endl;
    std::cout << "2*vc2 = " << 2. * vc1 << std::endl;
    vc2 += 2.;
    std::cout << "vc2(+=2) = " << vc2 << std::endl;
    vc2 -= 2.;
    std::cout << "vc2(-=2) = " << vc2 << std::endl;
    vc2 *= 2.;
    std::cout << "vc2(*=2) = " << vc2 << std::endl;
    vc2 /= 2.;
    std::cout << "vc2(/=2) = " << vc2 << std::endl;

// Tests of gslpp::matrix<double>
  std::cout << "Id(3) = " << matrix<double>::Id(3) << std::endl;
  std::cout << "\nTests of gslpp::matrix<double>" << std::endl;
  matrix<double> md1(3,3,3.);
  std::cout << "md1(=(3)) = " << md1 << std::endl;
  matrix<double> md2(md1);
  std::cout << "md2(=md1) = " << md2 << std::endl;
  matrix<double> mdsub(2,2,2.);
  std::cout << "md2sub = " << mdsub << std::endl;
  md2.assign(1,1,mdsub);
  std::cout << "md2.assign(1,1,mdsub) = " << md2 << std::endl;
  std::cout << "size_i(md2) = " << md2.size_i() << std::endl;
  std::cout << "size_j(md2) = " << md2.size_j() << std::endl;
  matrix<double> md3(3,3,0.);
  std::cout << "md3(=(0)) = " << md3 << std::endl;
  md3 = md2;
  std::cout << "md3(=md2) = " << md3 << std::endl;
  std::cout << "md3(2,1) = " << md3(2,1) << std::endl;
  md3(2,1)=5.;
  md3(0,2)=1.;
  std::cout << "md3^T = " << md3.transpose() << std::endl;
  std::cout << "md3^-1 = " << md3.inverse() << std::endl;
  std::cout << "md3^-1*md3 = " << md3.inverse()*md3 << std::endl;
  std::cout << "md3*md3^-1 = " << md3*md3.inverse() << std::endl;
  std::cout << "md1+md2 = " << md1+md2 << std::endl;
  std::cout << "md1-md2 = " << md1-md2 << std::endl;
  std::cout << "md1*md2 = " << md1*md2 << std::endl;
  md2+=md1;
  std::cout << "md2(+=md1) = " << md2 << std::endl;
  md2-=md1;
  std::cout << "md2(-=md1) = " << md2 << std::endl;
  md2*=md1;
  std::cout << "md2(*=md1) = " << md2 << std::endl;
  std::cout << "md1+2 = " << md1+2. << std::endl;
  std::cout << "md1-2 = " << md1-2. << std::endl;
  std::cout << "md1*2 = " << md1*2. << std::endl;
  std::cout << "md1/2 = " << md1/2. << std::endl;
  std::cout << "2+md2 = " << 2.+md1 << std::endl;
  std::cout << "2-md2 = " << 2.-md1 << std::endl;
  std::cout << "2*md2 = " << 2.*md1 << std::endl;
  md2+=2.;
  std::cout << "md2(+=2) = " << md2 << std::endl;
  md2-=2.;
  std::cout << "md2(-=2) = " << md2 << std::endl;
  md2*=2.;
  std::cout << "md2(*=2) = " << md2 << std::endl;
  md2/=2.;
  std::cout << "md2(/=2) = " << md2 << std::endl;
  std::cout << "md1*vd1 = " << md1*vd1 << std::endl;
  std::cout << "vd2*md2 = " << vd2*md2 << std::endl;
  std::cout << "md1*vc1 = " << md1*vc1 << std::endl;
  std::cout << "vc2*md2 = " << vc2*md2 << std::endl;

// Tests of gslpp::matrix<complex>
  std::cout << "\nTests of gslpp::matrix<complex>" << std::endl;
  std::cout << "Id(3) = " << matrix<complex>::Id(3) << std::endl;
  matrix<complex> mc1(3,3,3.+2.*complex::i());
  std::cout << "mc1(=(3+2i)) = " << mc1 << std::endl;
  matrix<complex> mc2(mc1);
  std::cout << "mc2(=mc1) = " << mc2 << std::endl;
  matrix<complex> mcsub(2,2,2.+complex::i());
  std::cout << "mcsub = " << mcsub << std::endl;
  mc2.assign(1,1,mcsub);
  std::cout << "mc2.assign(1,1,mcsub) = " << mc2 << std::endl;
  mc2.assign(1,1,mdsub);
  std::cout << "mc2.assign(1,1,mdsub) = " << mc2 << std::endl;
  std::cout << "size_i(mc2) = " << mc2.size_i() << std::endl;
  std::cout << "size_j(mc2) = " << mc2.size_j() << std::endl;
  matrix<complex> mc3(3,3,0.);
  std::cout << "mc3(=(0)) = " << mc3 << std::endl;
  mc3 = mc2;
  std::cout << "mc3(=mc2) = " << mc3 << std::endl;
  std::cout << "mc3(2,1) = " << mc3(2,1) << std::endl;
  mc3.assign(2,1,5.+4.*complex::i());
  mc3.assign(0,2,1.+4.*complex::i());
  std::cout << "mc3^T = " << mc3.transpose() << std::endl;
  std::cout << "mc3^* = " << mc3.hconjugate() << std::endl;
  std::cout << "mc3^-1 = " << mc3.inverse() << std::endl;
  std::cout << "mc3^-1 mc3 = " << mc3.inverse()*mc3 << std::endl;
  std::cout << "mc3 mc3^-1 = " << mc3*mc3.inverse() << std::endl;
  std::cout << "mc1+mc2 = " << mc1+mc2 << std::endl;
  std::cout << "mc1-mc2 = " << mc1-mc2 << std::endl;
  std::cout << "mc1*mc2 = " << mc1*mc2 << std::endl;
  mc2+=mc1;
  std::cout << "mc2(+=mc1) = " << mc2 << std::endl;
  mc2-=mc1;
  std::cout << "mc2(-=mc1) = " << mc2 << std::endl;
  mc2*=mc1;
  std::cout << "mc2(*=mc1) = " << mc2 << std::endl;
  std::cout << "mc1+2 = " << mc1+2. << std::endl;
  std::cout << "mc1-2 = " << mc1-2. << std::endl;
  std::cout << "mc1*2 = " << mc1*2. << std::endl;
  std::cout << "mc1/2 = " << mc1/2. << std::endl;
  std::cout << "2+mc2 = " << 2.+mc1 << std::endl;
  std::cout << "2-mc2 = " << 2.-mc1 << std::endl;
  std::cout << "2*mc2 = " << 2.*mc1 << std::endl;
  mc2+=2.;
  std::cout << "mc2(+=2) = " << mc2 << std::endl;
  mc2-=2.;
  std::cout << "mc2(-=2) = " << mc2 << std::endl;
  mc2*=2.;
  std::cout << "mc2(*=2) = " << mc2 << std::endl;
  mc2/=2.;
  std::cout << "mc2(/=2) = " << mc2 << std::endl;
  std::cout << "mc1*vc1 = " << mc1*vc1 << std::endl;
  std::cout << "vc2*mc2 = " << vc2*mc2 << std::endl;
  std::cout << "mc1*vd1 = " << mc1*vd1 << std::endl;
  std::cout << "vd2*mc2 = " << vd2*mc2 << std::endl;

  std::cout << "md1 + mc1 = " << md1 + mc1 << std::endl;
  std::cout << "md1 - mc1 = " << md1 - mc1 << std::endl;
}

void test2() {
    std::cout << "gslpptest test 2" << std::endl;
    matrix<complex> m(3,3,0.);
    m.assign(0,0,1.);
    m.assign(0,0,1.+0.*complex::i());
    m.assign(0,1,complex(1.435,.7548,true));
    m.assign(0,2,complex(0.2,.25,true));
    m.assign(1,0,complex(1.435,-.7548,true));
    m.assign(2,0,complex(0.2,-.5,true));
    m.assign(1,1,complex(-0.12,0.,true));
    m.assign(1,2,complex(7.12,0.07,true));
    m.assign(2,1,complex(2.12,-0.07,true));
    m.assign(2,2,complex(25.,0.,true));

    std::cout << m << std::endl;
    matrix<complex> U(m);
    matrix<complex> V(m);
    vector<double> S(3,0.);
    m.singularvalue(U,V,S);
    std::cout << S << std::endl;
 //   std::cout << U.hconjugate()*U << std::endl;
    std::cout << (U.hconjugate()*m)*V << std::endl;

//    std::cout << "%TEST_FAILED% time=0 testname=test2 (gslpptest) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% gslpptest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (gslpptest)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (gslpptest)" << std::endl;

    std::cout << "%TEST_STARTED% test2 (gslpptest)\n" << std::endl;
    test2();
    std::cout << "%TEST_FINISHED% time=0 test2 (gslpptest)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

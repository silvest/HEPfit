#include <iostream>
#include "gslpp_complex.h"
#include "gslpp_matrix_double.h"
#include "gslpp_matrix_complex.h"
#include "gslpp_vector_double.h"
#include "gslpp_vector_complex.h"
#include "Expanded.h"

int main(void) {

   gslpp::complex zi=gslpp::complex::i();
   std::vector<double> sd={10.,5.,1.};
   std::vector<gslpp::complex> sc={30.+zi,2.+3.*zi,1.+4.*zi};
   gslpp::matrix<double> md1(2,2);
   md1(0,0)=10.;md1(1,0)=20.;md1(0,1)=5.;md1(1,1)=1.;
   gslpp::matrix<double> md2(2,2);
   md2(0,0)=-3.;md2(1,0)=30.;md2(0,1)=-5.;md2(1,1)=4.;
   gslpp::matrix<gslpp::complex> mc1(2,2);
   mc1.assign(0,0,9.+2*zi);
   mc1.assign(1,0,19.-4*zi);
   mc1.assign(0,1,4.-6*zi);
   mc1.assign(1,1,3.+2*zi);
   gslpp::matrix<gslpp::complex> mc2(2,2);
   mc2.assign(0,0,-8.+3*zi);
   mc2.assign(1,0,11.-5*zi);
   mc2.assign(0,1,2.+5*zi);
   mc2.assign(1,1,-3.+4*zi);
   gslpp::vector<double> vd1(2);
   gslpp::vector<double> vd2(2);
   vd1(0)=14.;vd1(1)=-5;
   vd2(0)=-9.;vd2(1)=3;

   gslpp::vector<gslpp::complex> vc1(2);
   gslpp::vector<gslpp::complex> vc2(2);
   vc1.assign(0,4.-2.*zi);vc1.assign(1,6.-8.*zi);
   vc2.assign(0,5.+3.*zi);vc2.assign(1,4.-12.*zi);

   Expanded<double> esd(sd);
   Expanded<gslpp::complex> esc(sc);
   std::vector<gslpp::matrix<double> > mdv={md1,md2};
   Expanded<gslpp::matrix<double> > emd(mdv);
   std::vector<gslpp::matrix<gslpp::complex> > mcv={mc1,mc2};
   Expanded<gslpp::matrix<gslpp::complex> > emc(mcv);
   std::vector<gslpp::vector<double> > vdv={vd1,vd2};
   Expanded<gslpp::vector<double> > evd(vdv);
   std::vector<gslpp::vector<gslpp::complex> > vcv={vc1,vc2};
   Expanded<gslpp::vector<gslpp::complex> > evc(vcv);

// Print Input
   std::cout << std::endl;
   std::cout << "esd " << esd << std::endl;
   std::cout << "esc " << esc << std::endl;
   std::cout << "evd " << evd << std::endl;
   std::cout << "evc " << evc << std::endl;
   std::cout << "emd " << emd << std::endl;
   std::cout << "emc " << emc << std::endl;
   std::cout << "--------------"<< std::endl;



// Print Output
   

  // Zero
//  std::cout << "--------- zeri ----------" << std::endl;
//  std::cout << sdsc - scsd << std::endl;
//  std::cout << sdmd - mdsd << std::endl;
//  std::cout << scmd - mdsc << std::endl;
//  std::cout << scmc - mcsc << std::endl;
//  std::cout << sdmc - mcsd << std::endl;
//  std::cout << "-------------------" << std::endl;

//
std::cout << "sd*sd " << esd*esd << std::endl;
std::cout << "sd*sc " << esd*esc << std::endl;
std::cout << "sd*vd " << esd*evd << std::endl;
std::cout << "sd*vc " << esd*evc << std::endl;
std::cout << "sd*md " << esd*emd << std::endl;
std::cout << "sd*mc " << esd*emc << std::endl;

std::cout << "sc*sd " << esc*esd << std::endl;
std::cout << "sc*sc " << esc*esc << std::endl;
std::cout << "sc*vd " << esc*evd << std::endl;
std::cout << "sc*vc " << esc*evc << std::endl;
std::cout << "sc*md " << esc*emd << std::endl;
std::cout << "sc*mc " << esc*emc << std::endl;

std::cout << "vd*sd " << evd*esd << std::endl;
std::cout << "vd*sc " << evd*esc << std::endl;
std::cout << "vd*vd " << evd*evd << std::endl;
std::cout << "vd*vc " << evd*evc << std::endl;
std::cout << "vd*md " << evd*emd << std::endl;
std::cout << "vd*mc " << evd*emc << std::endl;

std::cout << "vc*sd " << evc*esd << std::endl;
std::cout << "vc*sc " << evc*esc << std::endl;
std::cout << "vc*vd " << evc*evd << std::endl;
std::cout << "vc*vc " << evc*evc << std::endl;
std::cout << "vc*md " << evc*emd << std::endl;
std::cout << "vc*mc " << evc*emc << std::endl;

std::cout << "md*sd " << emd*esd << std::endl;
std::cout << "md*sc " << emd*esc << std::endl;
std::cout << "md*vd " << emd*evd << std::endl;
std::cout << "md*vc " << emd*evc << std::endl;
std::cout << "md*md " << emd*emd << std::endl;
std::cout << "md*mc " << emd*emc << std::endl;

std::cout << "mc*sd " << emc*esd << std::endl;
std::cout << "mc*sc " << emc*esc << std::endl;
std::cout << "mc*vd " << emc*evd << std::endl;
std::cout << "mc*vc " << emc*evc << std::endl;
std::cout << "mc*md " << emc*emd << std::endl;
std::cout << "mc*mc " << emc*emc << std::endl;

   
//  std::cout << "sd + sd " << esd + esd << std::endl;
//  std::cout << "sd + sc " << esd + esc << std::endl;
//  std::cout << "sc + sd " << esc + esd << std::endl;
//  std::cout << "sc + sc " << esc + esc << std::endl;
//  std::cout << "vd + vd " << evd + evd << std::endl;
//  std::cout << "vd + vc " << evd + evc << std::endl;
//  std::cout << "vc + vd " << evc + evd << std::endl;
//  std::cout << "vc + vc " << evc + evc << std::endl;
//  std::cout << "md + md " << emd + emd << std::endl;
//  std::cout << "md + mc " << emd + emc << std::endl;
//  std::cout << "mc + md " << emc + emd << std::endl;
//  std::cout << "mc + mc " << emc + emc << std::endl;
//
//  std::cout << "sd - sd " << esd - esd << std::endl;
//  std::cout << "sd - sc " << esd - esc << std::endl;
//  std::cout << "sc - sd " << esc - esd << std::endl;
//  std::cout << "sc - sc " << esc - esc << std::endl;
//  std::cout << "vd - vd " << evd - evd << std::endl;
//  std::cout << "vd - vc " << evd - evc << std::endl;
//  std::cout << "vc - vd " << evc - evd << std::endl;
//  std::cout << "vc - vc " << evc - evc << std::endl;
//  std::cout << "md - md " << emd - emd << std::endl;
//  std::cout << "md - mc " << emd - emc << std::endl;
//  std::cout << "mc - md " << emc - emd << std::endl;
//  std::cout << "mc - mc " << emc - emc << std::endl;

  
//

//  std::cout << "esd*5 = " << esd*5. << std::endl;
//  std::cout << "esd*(4-2I) = " << esd*(4.-2.*zi) << std::endl;
//  std::cout << "esd*md1 = " << esd*md1 << std::endl;
//  std::cout << "esd*mc1 = " << esd*mc1 << std::endl;
//
//
//  std::cout << "esc*5 = " << esc*5. << std::endl;
//  std::cout << "esd*(4-2I) = " << esd*(4.-2.*zi) << std::endl;
//  std::cout << "esc*md1 = " << esc*md1 << std::endl;
//  std::cout << "esc*mc1 = " << esc*mc1 << std::endl;
//
//  std::cout << "emd*5 = " << emd*5. << std::endl;
//  std::cout << "emd*(4-2I) = " << emd*(4.-2.*zi) << std::endl;
//  std::cout << "emd*md1 = " << emd*md1 << std::endl;
//  std::cout << "emd*mc1 = " << emd*mc1 << std::endl;
//
//  std::cout << "emc*5 = " << emc*5. << std::endl;
//  std::cout << "emc*(4-2I) = " << emc*(4.-2.*zi) << std::endl;
//  std::cout << "emc*md1 = " << emc*md1 << std::endl;
//  std::cout << "emc*mc1 = " << emc*mc1 << std::endl;
//  
//
//  std::cout << "zero1 = " << 5.*esd - esd*5. << std::endl;
//  std::cout << "zero1 = " << 5.*esc - esc*5. << std::endl;
//  std::cout << "zero1 = " << 5.*emd - emd*5. << std::endl;
//  std::cout << "zero1 =" <<  5.*emc - emc*5. << std::endl;
//
//  std::cout << "zero2 = " << (4-2*zi)*esd - esd*(4-2*zi) << std::endl;
//  std::cout << "zero2 = " << (4-2*zi)*esc - esc*(4-2*zi)<< std::endl;
//  std::cout << "zero2 = " << (4-2*zi)*emd - emd*(4-2*zi)<< std::endl;
//  std::cout << "zero2 = " << (4-2*zi)*emc - emc*(4-2*zi)<< std::endl;
//
//
//  std::cout << "zero3 = " << md1*esd - esd*md1 << std::endl;
//  std::cout << "zero3 = " << md1*esc - esc*md1<< std::endl;
//  std::cout << "[md1,emd] = " << md1*emd - emd*md1 << std::endl;
//  std::cout << "[md1,emc] = " << md1*emc - emc*md1<< std::endl;
//
//  std::cout << "zero4 = " << mc1*esd - esd*mc1 << std::endl;
//  std::cout << "zero4 = " << mc1*esc - esc*mc1<< std::endl;
//  std::cout << "[mc1,emd] = " << mc1*emd - emd*mc1 << std::endl;
//  std::cout << "[mc1,emc] = " << mc1*emc - emc*mc1<< std::endl;
//

  return(0);
}

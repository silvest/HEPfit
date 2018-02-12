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

   double qq=5.;

    std::cout << isType<gslpp::complex>(zi) << std::endl;
    std::cout << isType<gslpp::complex>(qq) << std::endl;

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
   std::cout << "emd " << emd << std::endl;
   std::cout << "emc " << emc << std::endl;
   std::cout << "--------------"<< std::endl;

   Expanded<double> sdsd = esd*esd;
   Expanded<gslpp::complex> sdsc = esd*esc;
   Expanded<gslpp::complex> scsd = esc*esd;   
   Expanded<gslpp::complex> scsc = esc*esc;
   Expanded<gslpp::matrix<double> > sdmd = esd*emd;
   Expanded<gslpp::matrix<double> > mdsd = emd*esd;
   Expanded<gslpp::matrix<double> > mdmd = emd*emd;
   Expanded<gslpp::matrix<gslpp::complex> > scmd = esc*emd;
   Expanded<gslpp::matrix<gslpp::complex> > mdsc = emd*esc;
   Expanded<gslpp::matrix<gslpp::complex> > scmc = esc*emc;
   Expanded<gslpp::matrix<gslpp::complex> > mcsc = emc*esc;
   Expanded<gslpp::matrix<gslpp::complex> > sdmc = esd*emc;
   Expanded<gslpp::matrix<gslpp::complex> > mcsd = emc*esd;
   Expanded<gslpp::matrix<gslpp::complex> > mdmc = emd*emc;
   Expanded<gslpp::matrix<gslpp::complex> > mcmd = emc*emd;
   Expanded<gslpp::matrix<gslpp::complex> > mcmc = emc*emc;


// Print Output
   
   std::cout << "sdsd " << sdsd << std::endl;
   std::cout << "sdsc " << sdsc << std::endl;
   std::cout << "scsd " << scsd << std::endl;
   std::cout << "scsc " << scsc << std::endl;
   std::cout << "sdmd " << sdmd << std::endl;
   std::cout << "mdsd " << mdsd << std::endl;
   std::cout << "mdmd " << mdmd << std::endl;
   std::cout << "scmd " << scmd << std::endl;
   std::cout << "mdsc " << mdsc << std::endl;
   std::cout << "scmc " << scmc << std::endl;
   std::cout << "mcsc " << mcsc << std::endl;
   std::cout << "sdmc " << sdmc << std::endl;
   std::cout << "mcsd " << mcsd << std::endl;
   std::cout << "mdmc " << mdmc << std::endl;
   std::cout << "mcmd " << mcmd << std::endl;
   std::cout << "mcmc " << mcmc << std::endl;

  // Zero
  std::cout << "--------- zeri ----------" << std::endl;
  std::cout << sdsc - scsd << std::endl;
  std::cout << sdmd - mdsd << std::endl;
  std::cout << scmd - mdsc << std::endl;
  std::cout << scmc - mcsc << std::endl;
  std::cout << sdmc - mcsd << std::endl;
  std::cout << "-------------------" << std::endl;

//

  std::cout << "sc - md " << esc - emd << std::endl;
  std::cout << "md - sc " << emd - esc << std::endl;

  
//

  std::cout << "esd*5 = " << esd*5. << std::endl;
  std::cout << "esd*(4-2I) = " << esd*(4.-2.*zi) << std::endl;
  std::cout << "esd*md1 = " << esd*md1 << std::endl;
  std::cout << "esd*mc1 = " << esd*mc1 << std::endl;


  std::cout << "esc*5 = " << esc*5. << std::endl;
  std::cout << "esd*(4-2I) = " << esd*(4.-2.*zi) << std::endl;
  std::cout << "esc*md1 = " << esc*md1 << std::endl;
  std::cout << "esc*mc1 = " << esc*mc1 << std::endl;

  std::cout << "emd*5 = " << emd*5. << std::endl;
  std::cout << "emd*(4-2I) = " << emd*(4.-2.*zi) << std::endl;
  std::cout << "emd*md1 = " << emd*md1 << std::endl;
  std::cout << "emd*mc1 = " << emd*mc1 << std::endl;

  std::cout << "emc*5 = " << emc*5. << std::endl;
  std::cout << "emc*(4-2I) = " << emc*(4.-2.*zi) << std::endl;
  std::cout << "emc*md1 = " << emc*md1 << std::endl;
  std::cout << "emc*mc1 = " << emc*mc1 << std::endl;
  

  std::cout << "zero1 = " << 5.*esd - esd*5. << std::endl;
  std::cout << "zero1 = " << 5.*esc - esc*5. << std::endl;
  std::cout << "zero1 = " << 5.*emd - emd*5. << std::endl;
  std::cout << "zero1 =" <<  5.*emc - emc*5. << std::endl;

  std::cout << "zero2 = " << (4-2*zi)*esd - esd*(4-2*zi) << std::endl;
  std::cout << "zero2 = " << (4-2*zi)*esc - esc*(4-2*zi)<< std::endl;
  std::cout << "zero2 = " << (4-2*zi)*emd - emd*(4-2*zi)<< std::endl;
  std::cout << "zero2 = " << (4-2*zi)*emc - emc*(4-2*zi)<< std::endl;


  std::cout << "zero3 = " << md1*esd - esd*md1 << std::endl;
  std::cout << "zero3 = " << md1*esc - esc*md1<< std::endl;
  std::cout << "[md1,emd] = " << md1*emd - emd*md1 << std::endl;
  std::cout << "[md1,emc] = " << md1*emc - emc*md1<< std::endl;

  std::cout << "zero4 = " << mc1*esd - esd*mc1 << std::endl;
  std::cout << "zero4 = " << mc1*esc - esc*mc1<< std::endl;
  std::cout << "[mc1,emd] = " << mc1*emd - emd*mc1 << std::endl;
  std::cout << "[mc1,emc] = " << mc1*emc - emc*mc1<< std::endl;


  std::cout << "esd*evd = " << esd*evd <<std::endl;
  std::cout << "esd*evc = " << esd*evc <<std::endl;  
//  std::cout << "esc*evd = " << esc*evd <<std::endl;
  std::cout << "esc*evc = " << esc*evc <<std::endl;
  std::cout << "emd*evd = " << emd*evd <<std::endl;
  std::cout << "emc*evc = " << emc*evc <<std::endl;
  std::cout << "emd*evc = " << emd*evc <<std::endl;  
//  std::cout << "emc*evd = " << emc*evd <<std::endl;    
  return(0);
}

/* 
 * File:   CKM.h
 * Author: silvest
 *
 * Created on March 2, 2011, 2:07 PM
 */

#ifndef CKM_H
#define	CKM_H

#include <math.h>
#include <gslpp_complex.h>
#include <gslpp_matrix_complex.h>

class CKM {
public:
  CKM();
  CKM(const CKM&);
  ~CKM();

  void setWolfenstein(double, double, double, double);
  void setCKM(double, double, double, double);

  void getCKM(gslpp::matrix<gslpp::complex> &);

  // Wolfenstein parameters
  double getRho();
  double getEta();
  double getLambda();
  double getA();
  double getRhoNB();
  double getEtaNB();

  // Gilman parameterization
  double gets12();
  double gets13();
  double gets23();
  double getc12();
  double getc13();
  double getc23();
  double getdelta();

  // J_CP
  double getJcp();

  //Absolute values of CKM elements
  double getVud();
  double getVus();
  double getVub();
  double getVcd();
  double getVcs();
  double getVcb();
  double getVtd();
  double getVts();
  double getVtb();

  //Phases of CKM elements
  double getArgVud();
  double getArgVus();
  double getArgVub();
  double getArgVcd();
  double getArgVcs();
  double getArgVcb();
  double getArgVtd();
  double getArgVts();
  double getArgVtb();

  //Complex values of CKM elements
  gslpp::complex V_ud();
  gslpp::complex V_us();
  gslpp::complex V_ub();
  gslpp::complex V_cd();
  gslpp::complex V_cs();
  gslpp::complex V_cb();
  gslpp::complex V_td();
  gslpp::complex V_ts();
  gslpp::complex V_tb();

  // Angles
  double getBeta();
  double getGamma();
  double getAlpha();
  double getBetas();

// Lambda_q
  gslpp::complex getlamt();
  gslpp::complex getlamc();
  gslpp::complex getlamu();

  gslpp::complex getlamt_d();
  gslpp::complex getlamc_d();
  gslpp::complex getlamu_d();

  gslpp::complex getlamt_s();
  gslpp::complex getlamc_s();
  gslpp::complex getlamu_s();

// Sides
  double getRt();
  double getRts();
  double getRb();

private:
  double Rho, Eta, Lambda, A;
  double s12, s13, s23, delta;
  double c12, c23, c13;

  gslpp::complex Vud, Vcd, Vtd;
  gslpp::complex Vus, Vcs, Vts;
  gslpp::complex Vub, Vcb, Vtb;

};

#endif	/* CKM_H */


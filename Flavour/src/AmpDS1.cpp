/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AmpDS1.h"
#include "StandardModel.h"
#include "HeffDS1.h"
#include <sstream>


AmpDS1::AmpDS1(const StandardModel& SM_i)
: mySM(SM_i)
{
    mySM.initializeBParameter("BKd1");
    mySM.initializeBParameter("BKd3");
}

gslpp::complex AmpDS1::AmpDS1pp0(orders order)
{
    double GF = mySM.getGF();
    gslpp::complex Vud = mySM.getCKM().getV_ud();
    gslpp::complex Vus = mySM.getCKM().getV_us();
    //CKM parameter for the Wilson coeff.
    gslpp::complex tau = -((mySM.getCKM().getV_ts()).conjugate() * mySM.getCKM().getV_td())/((mySM.getCKM().getV_us()).conjugate() * mySM.getCKM().getV_ud());

    if (mySM.getFlavour().getHDS1().getCoeffDS1PP().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("AmpDK1::computeThValue(): requires cofficient of order"
                                 + out.str() + "not computed");
    }

    //Finding Wilson coefficients of the form C=z+tau*y where tau is a CKM parameter (scheme set to MSbar)
    gslpp::vector<gslpp::complex> ** allcoeffv = mySM.getFlavour().ComputeCoeffDS1PPv(
            mySM.getBKd1().getMu(), NDR);

    gslpp::vector<gslpp::complex> ** allcoeffz = mySM.getFlavour().ComputeCoeffDS1PPz(
            mySM.getBKd1().getMu(), NDR);

    gslpp::vector<gslpp::complex> allcoeffzLO = (*allcoeffz[LO]) + (*allcoeffz[LO_QED]);
    gslpp::vector<gslpp::complex> allcoeffzNLO = (*allcoeffz[NLO]) + (*allcoeffz[NLO_QED11]);
    gslpp::vector<gslpp::complex> allcoeffyLO = (*allcoeffv[LO]) + (*allcoeffv[LO_QED]);
    gslpp::vector<gslpp::complex> allcoeffyNLO = (*allcoeffv[NLO]) + (*allcoeffv[NLO_QED11]);
    for(int i = 0; i<2; i++){
        allcoeffzLO.assign(i,allcoeffyLO(i));
        allcoeffzNLO.assign(i,allcoeffyNLO(i));
        allcoeffyLO.assign(i,0.);
        allcoeffyNLO.assign(i,0.);
    }
    for(int i = 2; i<10; i++){
        allcoeffyLO.assign(i,allcoeffyLO(i)-allcoeffzLO(i));
        allcoeffyNLO.assign(i,allcoeffyNLO(i)-allcoeffzNLO(i));
    }

    gslpp::vector<double> meBKd1(mySM.getBKd1().getBpars());

    switch(order) {
        case NLO:
          if ( meBKd1(7) == 0 && meBKd1(8) == 0 && meBKd1(9) == 0 && mySM.getBKd1().getScheme() == LAT ){
            //If the matrix elements are given in the chiral basis on the lattice
            //Lattice Wilson coefficients
            gslpp::vector<gslpp::complex> ReW = ( (allcoeffzLO+allcoeffzNLO) + tau.real() * (allcoeffyLO+allcoeffyNLO) )
                                                * getChiralMatrixpp0() * getRISMOMTransMatrix(mySM.getBKd1().getMu(), FULLNLO) * getRIMatrixpp0();
            gslpp::vector<gslpp::complex> ImW = tau.imag() * (allcoeffyLO+allcoeffyNLO) * getChiralMatrixpp0()
                                                * getRISMOMTransMatrix(mySM.getBKd1().getMu(), FULLNLO) * getRIMatrixpp0();
            double CKMprod = (Vus.conjugate()*Vud).real();
            gslpp::complex ImA0 = M_SQRT1_2 * GF * CKMprod * ( (ImW * meBKd1 ) - ImW(2) * meBKd1(2) ) +
                          ( ImW(2)/ReW(2) ) * (getReA0() - M_SQRT1_2 * GF * CKMprod * ( ( ReW * meBKd1 ) - ReW(2) * meBKd1(2) ) );

            //The returned value for the amplitude is using the expt value of ReA0 to minimize the error on ImA0 (cfr. ArXiv:2004.09440)
            return gslpp::complex(getReA0() , ImA0.real() );

          } else if ( meBKd1(7) == 0 && meBKd1(8) == 0 && meBKd1(9) == 0 && mySM.getBKd1().getScheme() == LRI) {
            //If the me are given in the chiral basis renormalised in SMOM
            //Lattice Wilson coefficients
            gslpp::vector<gslpp::complex> ReW = ( (allcoeffzLO+allcoeffzNLO) + tau.real() * (allcoeffyLO+allcoeffyNLO) )
                                                * getChiralMatrixpp0() * getRISMOMTransMatrix(mySM.getBKd1().getMu(), FULLNLO);
            gslpp::vector<gslpp::complex> ImW = tau.imag() * (allcoeffyLO+allcoeffyNLO) * getChiralMatrixpp0()
                                                * getRISMOMTransMatrix(mySM.getBKd1().getMu(), FULLNLO);
            double CKMprod = (Vus.conjugate()*Vud).real();
            gslpp::complex ImA0 = M_SQRT1_2 * GF * CKMprod * ( (ImW * meBKd1 ) - ImW(2) * meBKd1(2) ) +
                          ( ImW(2)/ReW(2) ) * (getReA0() - M_SQRT1_2 * GF * CKMprod * ( ( ReW * meBKd1 ) - ReW(2) * meBKd1(2) ) );

            //The returned value for the amplitude is using the expt value of ReA0 to minimize the error on ImA0 (cfr. ArXiv:2004.09440)
            return gslpp::complex(getReA0() , ImA0.real() );

          } else {
            //If the me are given in the 10 basis renormalised in MSbar
            return M_SQRT1_2 * GF * (Vus.conjugate() * Vud) * ( (allcoeffzLO + allcoeffzNLO) + tau * (allcoeffyLO + allcoeffyNLO) ) * meBKd1;
          }
        case LO:
          if ( meBKd1(7) == 0 && meBKd1(8) == 0 && meBKd1(9) == 0 && mySM.getBKd1().getScheme() == LAT ){
            //Lattice Wilson coefficients
            gslpp::vector<gslpp::complex> ReW = ( allcoeffzLO + tau.real() * allcoeffyLO )
                                                * getChiralMatrixpp0() * getRISMOMTransMatrix(mySM.getBKd1().getMu(), FULLNLO) * getRIMatrixpp0();
            gslpp::vector<gslpp::complex> ImW = tau.imag() * allcoeffyLO * getChiralMatrixpp0()
                                                * getRISMOMTransMatrix(mySM.getBKd1().getMu(), FULLNLO) * getRIMatrixpp0();

            gslpp::complex ImA0 = M_SQRT1_2 * GF * ( Vus.conjugate() * Vud ).real() * ( (ImW * meBKd1 ) - ImW(2) * meBKd1(2) ) +
                          ( ImW(2)/ReW(2) ) * (getReA0() - M_SQRT1_2 * GF * ( Vus.conjugate() * Vud ).real() * ( ( ReW * meBKd1 ) - ReW(2) * meBKd1(2) ) );

            //The returned value for the amplitude is using the expt value of ReA0 to minimize the error on ImA0 (cfr. ArXiv:2004.09440)
            return gslpp::complex(getReA0() , ImA0.real() );

          } else if ( meBKd1(7) == 0 && meBKd1(8) == 0 && meBKd1(9) == 0 && mySM.getBKd1().getScheme() == LRI) {
            //If the me are given in the chiral basis renormalised in SMOM
            //Lattice Wilson coefficients
            gslpp::vector<gslpp::complex> ReW = ( allcoeffzLO + tau.real() * allcoeffyLO )
                                                * getChiralMatrixpp0() * getRISMOMTransMatrix(mySM.getBKd1().getMu(), FULLNLO);
            gslpp::vector<gslpp::complex> ImW = tau.imag() * allcoeffyLO * getChiralMatrixpp0()
                                                * getRISMOMTransMatrix(mySM.getBKd1().getMu(), FULLNLO);
            double CKMprod = (Vus.conjugate()*Vud).real();
            gslpp::complex ImA0 = M_SQRT1_2 * GF * CKMprod * ( (ImW * meBKd1 ) - ImW(2) * meBKd1(2) ) +
                          ( ImW(2)/ReW(2) ) * (getReA0() - M_SQRT1_2 * GF * CKMprod * ( ( ReW * meBKd1 ) - ReW(2) * meBKd1(2) ) );

            //The returned value for the amplitude is using the expt value of ReA0 to minimize the error on ImA0 (cfr. ArXiv:2004.09440)
            return gslpp::complex(getReA0() , ImA0.real() );

          } else {
            //If the me are given in the 10 basis renormalised in MSbar
            return M_SQRT1_2 * GF * (Vus.conjugate() * Vud) * ( allcoeffzLO + tau * allcoeffyLO ) * meBKd1;
          }
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("AmpDK1::AmpDK(): order " + out.str() + "not implemented");
    }
}

gslpp::complex AmpDS1::AmpDS1pp2(orders order)
{
    if (mySM.getFlavour().getHDS1().getCoeffDS1PP().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("AmpDK1::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }

    double GF = mySM.getGF();
    gslpp::complex Vud = mySM.getCKM().getV_ud();
    gslpp::complex Vus = mySM.getCKM().getV_us();
    //CKM parameter for the Wilson coeff.
    gslpp::complex tau = -((mySM.getCKM().getV_ts()).conjugate() * mySM.getCKM().getV_td())/((mySM.getCKM().getV_us()).conjugate() * mySM.getCKM().getV_ud());


    //Finding Wilson coefficients of the form C=z+tau*y where tau is a CKM parameter (scheme set to MSbar)
    gslpp::vector<gslpp::complex> ** allcoeffv = mySM.getFlavour().ComputeCoeffDS1PPv(
            mySM.getBKd3().getMu(), NDR);

    gslpp::vector<gslpp::complex> ** allcoeffz = mySM.getFlavour().ComputeCoeffDS1PPz(
            mySM.getBKd3().getMu(), NDR);

    gslpp::vector<gslpp::complex> allcoeffzLO = (*allcoeffz[LO]) + (*allcoeffz[LO_QED]);
    gslpp::vector<gslpp::complex> allcoeffzNLO = (*allcoeffz[NLO]) + (*allcoeffz[NLO_QED11]);
    gslpp::vector<gslpp::complex> allcoeffyLO = (*allcoeffv[LO]) + (*allcoeffv[LO_QED]);
    gslpp::vector<gslpp::complex> allcoeffyNLO = (*allcoeffv[NLO]) + (*allcoeffv[NLO_QED11]);
    for(int i = 0; i<2; i++){
        allcoeffzLO.assign(i,allcoeffyLO(i));
        allcoeffzNLO.assign(i,allcoeffyNLO(i));
        allcoeffyLO.assign(i,0.);
        allcoeffyNLO.assign(i,0.);
    }
    for(int i = 2; i<10; i++){
        allcoeffyLO.assign(i,allcoeffyLO(i)-allcoeffzLO(i));
        allcoeffyNLO.assign(i,allcoeffyNLO(i)-allcoeffzNLO(i));
    }


    gslpp::vector<double> meBKd3(mySM.getBKd3().getBpars());

    //If the me are given in the RI scheme, we convert it to MSbar
    if( mySM.getBKd3().getScheme() == LRI){
      meBKd3 =  getRISMOMTransMatrix(mySM.getBKd3().getMu(), FULLNLO) * meBKd3;
    }

    //We need to check if the B parameters are given in the chiral basis or not. If so, we convert it to the 10 op. basis
    if( meBKd3(1) == 0 && meBKd3(2) == 0 && meBKd3(3) == 0 && meBKd3(4) == 0 && meBKd3(8) == 0 && meBKd3(9) == 0 ){
      //The 1/sqrt(3) is the Clebsh-Gordan for the case where the me are evaluated for the virtual process K+->pi+ pi+
      meBKd3 = (1./sqrt(3.)) * getChiralMatrixpp2() * meBKd3;
    }

    switch(order) {
        case NLO:
           return M_SQRT1_2 * GF * (Vud.conjugate() * Vus) * ((allcoeffzLO + allcoeffzNLO) + tau * (allcoeffyLO + allcoeffyNLO)) * meBKd3;
        case LO:
            return M_SQRT1_2 * GF * (Vud.conjugate() * Vus) * (allcoeffzLO + tau * allcoeffyLO) * meBKd3;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("AmpDK1::AmpDK(): order " + out.str() + "not implemented");;
    }
}

double AmpDS1::getReA0(){
  double MP0 = mySM.getMesons(QCD::P_0).getMass();
  double MPp = mySM.getMesons(QCD::P_P).getMass();
  double MK0 = mySM.getMesons(QCD::K_0).getMass();

  //Evaluate ReA0 from expt input
  double GammaKstotal(mySM.getMesons(QCD::K_S).computeWidth());
  double GammaKSP0P0 = mySM.getOptionalParameter("Br_Ks_P0P0") * GammaKstotal;
  double GammaKSPpPm = mySM.getOptionalParameter("Br_Ks_PpPm") * GammaKstotal;

  double phasespacePpPm = 0.5 * sqrt(MK0*MK0 - 4.*MPp*MPp);
  double phasespaceP0P0 = 0.5 * sqrt(MK0*MK0 - 4.*MP0*MP0);

  double APpPm = sqrt(GammaKSPpPm * 8. * M_PI * MK0*MK0 / phasespacePpPm);
  double AP0P0 = sqrt(GammaKSP0P0 * 16. * M_PI * MK0*MK0 / phasespaceP0P0);

  double ReA0 = sqrt((APpPm * APpPm + 0.5 * AP0P0 * AP0P0)/2.);

  return ReA0;
}

double AmpDS1::getReA2(){
  double MP0 = mySM.getMesons(QCD::P_0).getMass();
  double MPp = mySM.getMesons(QCD::P_P).getMass();
  double MKp = mySM.getMesons(QCD::K_P).getMass();

  double GammaKptotal(mySM.getMesons(QCD::K_P).computeWidth());
  double GammaKp = mySM.getOptionalParameter("Br_Kp_P0Pp") * GammaKptotal;

  double phasespacePpP0 = sqrt((MKp*MKp / 4.) - (MPp*MPp + MP0*MP0)/2. + (MP0*MP0 - MPp*MPp)*(MP0*MP0 - MPp*MPp)/(4.*MKp*MKp));

  double ReA2 = sqrt((2./3.) * 8. * M_PI * GammaKp * MKp * MKp / phasespacePpP0);

  return ReA2;
}

gslpp::matrix<double> AmpDS1::getChiralMatrixpp0() const{
  //Converion from Chiral basis to 10 op. basis
  gslpp::matrix<double> chiral_to_10(10,10,0);

  chiral_to_10(0,0) = 1./5.;
  chiral_to_10(0,1) = 1.;
  chiral_to_10(1,0) = 1./5.;
  chiral_to_10(1,2) = 1.;
  chiral_to_10(2,1) = 3.;
  chiral_to_10(2,2) = 2.;
  chiral_to_10(3,1) = 2.;
  chiral_to_10(3,2) = 3.;
  chiral_to_10(4,3) = 1.;
  chiral_to_10(5,4) = 1.;
  chiral_to_10(6,5) = 1.;
  chiral_to_10(7,6) = 1.;
  chiral_to_10(8,0) = 3./10.;
  chiral_to_10(8,2) = -1.;
  chiral_to_10(9,0) = 3./10.;
  chiral_to_10(9,1) = -1.;

  return chiral_to_10;
}

gslpp::matrix<double> AmpDS1::getChiralMatrixpp2() const{
  //Conversion to the 10 op. basis from the chiral one
  gslpp::matrix<double> chiral_to_10(10,10,0);
  chiral_to_10(0,0) =1./3.;
  chiral_to_10(1,0) = 1./3.;
  chiral_to_10(6,5) = 1./2.;
  chiral_to_10(7,6) = 1./2.;
  chiral_to_10(8,0) = 1./2.;
  chiral_to_10(9,0) = 1./2.;

  return chiral_to_10;
}

gslpp::matrix<double> AmpDS1::getRIMatrixpp0() const{
  //Renormalization matrix for the bare lattice matrix elements of the chiral basis in the SMOMqq scheme
  gslpp::matrix<double> lat_SMOMqq(10,10,0);
  lat_SMOMqq(0,0) = mySM.getOptionalParameter("Zqq00");
  lat_SMOMqq(1,1) = mySM.getOptionalParameter("Zqq11");
  lat_SMOMqq(1,2) = mySM.getOptionalParameter("Zqq12");
  lat_SMOMqq(1,3) = mySM.getOptionalParameter("Zqq13");
  lat_SMOMqq(1,4) = mySM.getOptionalParameter("Zqq14");
  lat_SMOMqq(2,1) = mySM.getOptionalParameter("Zqq21");
  lat_SMOMqq(2,2) = mySM.getOptionalParameter("Zqq22");
  lat_SMOMqq(2,3) = mySM.getOptionalParameter("Zqq23");
  lat_SMOMqq(2,4) = mySM.getOptionalParameter("Zqq24");
  lat_SMOMqq(3,1) = mySM.getOptionalParameter("Zqq31");
  lat_SMOMqq(3,2) = mySM.getOptionalParameter("Zqq32");
  lat_SMOMqq(3,3) = mySM.getOptionalParameter("Zqq33");
  lat_SMOMqq(3,4) = mySM.getOptionalParameter("Zqq34");
  lat_SMOMqq(4,1) = mySM.getOptionalParameter("Zqq41");
  lat_SMOMqq(4,2) = mySM.getOptionalParameter("Zqq42");
  lat_SMOMqq(4,3) = mySM.getOptionalParameter("Zqq43");
  lat_SMOMqq(4,4) = mySM.getOptionalParameter("Zqq44");
  lat_SMOMqq(5,5) = mySM.getOptionalParameter("Zqq55");
  lat_SMOMqq(5,6) = mySM.getOptionalParameter("Zqq56");
  lat_SMOMqq(6,5) = mySM.getOptionalParameter("Zqq65");
  lat_SMOMqq(6,6) = mySM.getOptionalParameter("Zqq66");

  return lat_SMOMqq;
}

gslpp::matrix<double> AmpDS1::getRISMOMTransMatrix(double mu, orders order) const {
  gslpp::matrix<double> smom_ms_conversion(10,10,0);
  int Nc = 3; //Number of colors
  double C0 = 2.34391;
  smom_ms_conversion(0,0) = -12.*log(2.)/Nc+12.*log(2.)+9./Nc-9.;
  smom_ms_conversion(1,1) = -12.*log(2.)/Nc+8.*Nc/5.+9./Nc-12./5.;
  smom_ms_conversion(1,2) = 12.*log(2.)+12.*Nc/5.-53./5.;
  smom_ms_conversion(2,1) = 12.*log(2.)-12.*Nc/5.+2./(3.*Nc)-263./45.;
  smom_ms_conversion(2,2) = -12.*log(2.)/Nc-18.*Nc/5+85./(9.*Nc)+26./15.;
  smom_ms_conversion(2,3) = 2./(9.*Nc);
  smom_ms_conversion(2,4) = -2./9.;
  smom_ms_conversion(3,3) = 3.*C0/(2.*Nc)-2.*log(2.)/Nc-2./Nc;
  smom_ms_conversion(3,4) = -3.*C0/2.+2.*log(2.)+2.;
  smom_ms_conversion(4,1) = 5./Nc-10./3.;
  smom_ms_conversion(4,2) = 10./(3.*Nc)-5.;
  smom_ms_conversion(4,3) = 2.*log(2.)+5./(3.*Nc)-2.;
  smom_ms_conversion(4,4) = -3.*C0*Nc/2.+3.*C0/(2.*Nc)-2.*log(2.)/Nc+4.*Nc-2./Nc-5./3.;
  smom_ms_conversion(5,5) = 3.*C0/(2.*Nc)-2.*log(2.)/Nc-2./Nc;
  smom_ms_conversion(5,6) = -3.*C0/2.+2.*log(2.)+2.;
  smom_ms_conversion(6,5) = 2.*log(2.)-2.;
  smom_ms_conversion(6,6) = -3.*C0*Nc/2.+3.*C0/(2.*Nc)-2.*log(2.)/Nc+4.*Nc-2./Nc;

  smom_ms_conversion *= mySM.Als(mu, order)/(4.*M_PI);
  smom_ms_conversion += gslpp::matrix<double>::Id(10);

  return smom_ms_conversion;
}

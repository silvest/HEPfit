/*
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

/**
 * @example libmode_config.cpp
 * This is an example of how to compute observables from the input parameters
 * defined in a model configuration file.
 *
 */

#include <iostream>
#include <ComputeObservables.h>
#include "HeffDF1.h"
#include "HeffDB1.h"
#include "EvolDB1Mll.h"
#include "EvolDB1bsg.h"
#include "EvolBsmm.h"
#include "BXqll.h"
#include <gsl/gsl_sf_zeta.h>


int main(void) {

    /* Define the model configuration file.                        */
    /* Here it is passed as the first argument to the executable.  */
    /* The model configuration file provides the default values of */
    /* the mandatory model parameters.                             */
    std::string ModelConf = "StandardModel.conf";

    /* Define a map for the parameters to be varied. */
    std::map<std::string, double> DPars;

    /* Create objects of the classes ModelFactory and ThObsFactory */
    ModelFactory ModelF;
    ThObsFactory ThObsF;

    /* register user-defined model named ModelName defined in class ModelClass using the following syntax: */
    /* ModelF.addModelToFactory(ModelName, boost::factory<ModelClass*>() ) */

    /* register user-defined ThObservable named ThObsName defined in class ThObsClass using the following syntax: */
    /* ThObsF.addObsToFactory(ThObsName, boost::factory<ThObsClass*>() )*/

    /* Create an object of the class ComputeObservables. */
    ComputeObservables CO(ModelF, ThObsF, ModelConf);
    StandardModel& mySM = *CO.getModel();

    double mub = 5.;
    double myMW = mySM.Mw();
    std::cout << "MW: " << myMW <<  std::endl;
    std::cout << "sW2: " << mySM.sW2() <<  std::endl;
    std::cout << "1/alphaMz: " << 1./mySM.alphaMz() <<  std::endl;
    std::cout << "Als5: " << mySM.Als(5.,FULLNNNLO, true) <<  std::endl;
    std::cout << "Als5: " << mySM.Als(5.,FULLNNNLO) <<  std::endl;
    std::cout << "Alstilde5*4*pi: " <<mySM.Alstilde5(5.) * 4. * M_PI  <<  std::endl;
    std::cout << "Als120: " << mySM.Als(120.,FULLNNNLO, true) <<  std::endl;
    std::cout << "Alstilde120*4*pi: " << mySM.Alstilde5(120.) * 4. * M_PI  <<  std::endl;
    std::cout << "Ale5: " << mySM.Ale(5.,FULLNLO) <<  std::endl;
    std::cout << "Ale120: " << mySM.Ale(120.,FULLNLO) <<  std::endl;
    

    HeffDF1 Heff("CPL", mySM, NNLO, NLO_QED22);
    HeffDB1 HDB1(mySM);

//    CO.AddObservable("Rlow_BXsee");
//    /* Get the map of observables if necessary. */
//    std::map<std::string, double> DObs = CO.getObservables();
//
//    DObs = CO.compute(DPars);
//    
//    std::cout << DObs.at("Rlow_BXsee") << std::endl;

    std::cout << "Ale5_smm: " << HDB1.getUBsmm().alphatilde_e(5.) * 4. * M_PI <<  std::endl;    
    std::cout << "Ale120_smm: "  << HDB1.getUBsmm().alphatilde_e(120.) * 4. * M_PI <<  std::endl;    

    std::cout << "%SUITE_STARTING% Evolutor" << std::endl;
    std::cout << "%SUITE_STARTED% *****" << std::endl;
    gslpp::vector<gslpp::complex> ** allcoeff;
    gslpp::vector<gslpp::complex> ** allcoeff_smm;
    gslpp::vector<gslpp::complex> ** allcoeff_BMll;

    gslpp::matrix<gslpp::complex> myVCKM(mySM.getVCKM());
    double sw = sqrt( (M_PI * mySM.getAle() ) / ( sqrt(2.) * mySM.getGF() * mySM.Mw() * mySM.Mw() ) );
//    double sw = sqrt(mySM.sW2());
    double as5 =  mySM.Als(mub,FULLNNNLO, true) / 4. / M_PI;
    double ae5 = mySM.Ale(mub,FULLNLO) / 4. / M_PI;
    double k5 = ae5/as5;
//
//    std::cout << "%SUITE_STARTED% *****" << std::endl;
//
      allcoeff = Heff.ComputeCoeff(mub);
      allcoeff_smm = HDB1.ComputeCoeffsmumu(mub, NDR);
      allcoeff_BMll = HDB1.ComputeCoeffBMll(mub, StandardModel::MU, NDR);
      
//    allcoeff1 = HDB1.ComputeCoeffsgamma(5.);
//      std::cout << Heff.getEvol().DF1Evol(5., 120., LO, NDR)+Heff.getEvol().DF1Evol(5., 120., LO_QED, NDR) << std::endl; 

      std::cout << std::endl << "00:" << std::endl;
      std::cout << Heff.LowScaleCoeff(00) <<  std::endl;
      std::cout << (*(allcoeff_smm[LO])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;
      std::cout << Heff.LowScaleCoeff(00)-(*(allcoeff_smm[LO])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;

      std::cout << std::endl << "10:" << std::endl;
      std::cout << Heff.LowScaleCoeff(10) <<  std::endl;
      std::cout << (*(allcoeff_smm[NLO])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;
      std::cout << Heff.LowScaleCoeff(10)-(*(allcoeff_smm[NLO])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;

      std::cout << std::endl << "20:" << std::endl;
      std::cout << Heff.LowScaleCoeff(20) <<  std::endl;
      std::cout << (*(allcoeff_smm[NNLO])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;
      std::cout << Heff.LowScaleCoeff(20)-(*(allcoeff_smm[NNLO])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;
      
      std::cout << std::endl << "01:" << std::endl;
      std::cout << Heff.LowScaleCoeff(01) <<  std::endl;
      std::cout << (*(allcoeff_smm[LO_QED])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;
      std::cout << Heff.LowScaleCoeff(01)-(*(allcoeff_smm[LO_QED])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;

      std::cout << std::endl << "11:" << std::endl;
      std::cout << Heff.LowScaleCoeff(11) <<  std::endl;
      std::cout << (*(allcoeff_smm[NLO_QED11])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;
      std::cout << Heff.LowScaleCoeff(11)-(*(allcoeff_smm[NLO_QED11])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;

     
      std::cout << std::endl << "02:" << std::endl;
      std::cout << Heff.LowScaleCoeff(02) <<  std::endl;
      std::cout << (*(allcoeff_smm[NLO_QED02])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;
      std::cout <<  Heff.LowScaleCoeff(02) -(*(allcoeff_smm[NLO_QED02])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;

      std::cout << std::endl << "12:" << std::endl;
      std::cout << Heff.LowScaleCoeff(12) <<  std::endl;
      std::cout << (*(allcoeff_smm[NLO_QED12])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;
      std::cout << Heff.LowScaleCoeff(12)-(*(allcoeff_smm[NLO_QED12])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;

      std::cout << std::endl << "21:" << std::endl;
      std::cout << Heff.LowScaleCoeff(21) <<  std::endl;
      std::cout << (*(allcoeff_smm[NLO_QED21])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;
      std::cout << Heff.LowScaleCoeff(21)-(*(allcoeff_smm[NLO_QED21])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;

      std::cout << std::endl << "22:" << std::endl;
      std::cout << Heff.LowScaleCoeff(22) <<  std::endl;
      std::cout << (*(allcoeff_smm[NLO_QED22])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;
      std::cout << Heff.LowScaleCoeff(22)-(*(allcoeff_smm[NLO_QED22])) / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;

      gslpp::complex C10_low = (Heff.LowScaleCoeff(00))(7) + as5*(Heff.LowScaleCoeff(10))(7) +
                                as5*as5*(Heff.LowScaleCoeff(20))(7) + k5*(Heff.LowScaleCoeff(01))(7) +
                                as5*k5*(Heff.LowScaleCoeff(11))(7) + as5*as5*k5*(Heff.LowScaleCoeff(21))(7) +
                                k5*k5*(Heff.LowScaleCoeff(02))(7) + as5*k5*k5*(Heff.LowScaleCoeff(12))(7) +
                                as5*as5*k5*k5*(Heff.LowScaleCoeff(22))(7);
      
      gslpp::complex C10_student = (*(allcoeff_smm[LO]))(7)  +  as5*(*(allcoeff_smm[NLO]))(7) +
                                as5*as5*(*(allcoeff_smm[NNLO]))(7) + k5*(*(allcoeff_smm[LO_QED ]))(7) +
                                as5*k5*(*(allcoeff_smm[NLO_QED11]))(7) + k5*k5*(*(allcoeff_smm[NLO_QED02]))(7) +
                                as5*as5*k5*(*(allcoeff_smm[NLO_QED21]))(7) + as5*k5*k5*(*(allcoeff_smm[NLO_QED12]))(7) +
                                as5*as5*k5*k5*(*(allcoeff_smm[NLO_QED22]))(7);
      
      double C10_huber_11 = -4.222;
      double C10_huber_02 = 0.498e-2;
      double C10_huber_12 = -3.798;
      double C10_huber_21 = 6.380;
      double C10_huber_22 = -36.090;
      gslpp::complex C10_huber = as5*k5 * C10_huber_11 + k5*k5 * C10_huber_02 + as5*as5*k5 * C10_huber_21 +
                                as5*k5*k5 * C10_huber_12 + as5*as5*k5*k5 * C10_huber_22;
      
      std::cout << std::endl;
      std::cout << "C10: " << C10_low / ae5 << std::endl;
      std::cout << "C10_student: " << C10_student / ae5 / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;
      std::cout << "C10_BMll: " << (*(allcoeff_BMll[LO]))(9) + (*(allcoeff_BMll[NLO]))(9) << std::endl;
      std::cout << "C10_huber: " << C10_huber / ae5 << std::endl;
      
      gslpp::complex C9_low = (Heff.LowScaleCoeff(00))(6) + as5*(Heff.LowScaleCoeff(10))(6) +
                                as5*as5*(Heff.LowScaleCoeff(20))(6) + k5*(Heff.LowScaleCoeff(01))(6) +
                                as5*k5*(Heff.LowScaleCoeff(11))(6) + as5*as5*k5*(Heff.LowScaleCoeff(21))(6) +
                                k5*k5*(Heff.LowScaleCoeff(02))(6) + as5*k5*k5*(Heff.LowScaleCoeff(12))(6) +
                                as5*as5*k5*k5*(Heff.LowScaleCoeff(22))(6);
      gslpp::complex C9_student = (*(allcoeff_smm[LO]))(6)  +  as5*(*(allcoeff_smm[NLO]))(6) +
                                as5*as5*(*(allcoeff_smm[NNLO]))(6) + k5*(*(allcoeff_smm[LO_QED ]))(6) +
                                as5*k5*(*(allcoeff_smm[NLO_QED11]))(6) + k5*k5*(*(allcoeff_smm[NLO_QED02]))(6) +
                                as5*as5*k5*(*(allcoeff_smm[NLO_QED21]))(6) + as5*k5*k5*(*(allcoeff_smm[NLO_QED12]))(6) +
                                as5*as5*k5*k5*(*(allcoeff_smm[NLO_QED22]))(6);
      double C9_huber_01 = 3.722e-2;
      double C9_huber_11 = 1.934;
      double C9_huber_02 = 0.208e-2;
      double C9_huber_12 = -4.317;
      double C9_huber_21 = 3.538;
      double C9_huber_22 = 27.320;
      gslpp::complex C9_huber = k5 * C9_huber_01 + as5*k5 * C9_huber_11 + k5*k5 * C9_huber_02 +
                                as5*as5*k5 * C9_huber_21 + as5*k5*k5 * C9_huber_12 + as5*as5*k5*k5 * C9_huber_22;
      std::cout << std::endl;
      std::cout << "C9: " << C9_low / ae5 << std::endl;
      std::cout << "C9_student: " << C9_student / ae5 / (sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))) << std::endl;
      std::cout << "C9_BMll: " << (*(allcoeff_BMll[LO]))(8) + (*(allcoeff_BMll[NLO]))(8) << std::endl;
      std::cout << "C9_huber: " << C9_huber / ae5 << std::endl;
      
//    std::cout << (*(allcoeff[LO])) - *(allcoeff1[LO]) <<  std::endl;
//   
//    std::cout << "NLO:" << std::endl;
//    std::cout << (*(allcoeff[NLO])) <<  std::endl;
//    std::cout << (*(allcoeff[NLO])) - (*(allcoeff1[NLO])) <<  std::endl;
//
//    std::cout << "NNLO:" << std::endl;
//    std::cout << (*(allcoeff[NNLO]))  <<  std::endl;    
//    std::cout << (*(allcoeff[NNLO])) - (*(allcoeff1[NNLO])) <<  std::endl;

// test di CPL   
//    std::cout << "studente: " << as5*4.*M_PI << ", grande unificazione: " << mySM.Als(5., FULLNNNLO, true) << std::endl;
//    std::cout << "LO:" << std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[LO]));
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[LO])) - *(allcoeff_smm[LO]) <<  std::endl;
//    
//    std::cout << "NLO:" << std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[NLO])) - as5*(*(allcoeff_smm[NLO])) <<  std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[NLO])) <<  std::endl;
//
//    std::cout << "NNLO:" << std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[NNLO])) <<  std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[NNLO])) - as5*as5*(*(allcoeff_smm[NNLO])) <<  std::endl;
//
//    std::cout << "LO_QED:" << std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[LO_QED])) <<  std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[LO_QED])) - ae5/as5*(*(allcoeff_smm[LO_QED])) <<  std::endl;
//
//    std::cout << "NLO_QED11:" << std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[NLO_QED11])) <<  std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[NLO_QED11])) - ae5*(*(allcoeff_smm[NLO_QED11])) <<  std::endl;
//
//    std::cout << "NLO_QED02:" << std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[NLO_QED02])) <<  std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[NLO_QED02])) - ae5*ae5/as5/as5*(*(allcoeff_smm[NLO_QED02])) <<  std::endl;
//
//    std::cout << "NLO_QED12:" << std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[NLO_QED12])) <<  std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[NLO_QED12]))- ae5*ae5/as5*(*(allcoeff_smm[NLO_QED12])) <<  std::endl;
//
//    std::cout << "NLO_QED21:" << std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[NLO_QED21]))<<  std::endl;
//    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(allcoeff[NLO_QED21]))- ae5*as5*(*(allcoeff_smm[NLO_QED21])) <<  std::endl;
//
//    std::cout << "NLO_QED11+22:" << std::endl;
//    std::cout << myVCKM(2,2).conjugate() * myVCKM(2,1) * (*(allcoeff[NLO_QED11]) +
//            *(allcoeff[NLO_QED22])) - mySM.getGF() / sqrt(8.) * mySM.Mw() * mySM.Mw() / M_PI /M_PI *
//            (*(allcoeff_smm[NLO_QED11]) +  ae5 * (*(allcoeff_smm[NLO_QED22]))) <<  std::endl;
    
//   std::cout << *(allcoeff[NLO]) <<  std::endl;
//    std::cout << *(allcoeff[NNLO]) <<  std::endl;
//    std::cout << *(allcoeff[LO_QED]) <<  std::endl;    
//    std::cout << Heff.LowScaleCoeff(22) <<  std::endl;    
    
//    std::cout << *(mySM.getMatching().CMDF1("C",2)[0].getCoeff(LO)) <<  std::endl;    
//    std::cout << *(mySM.getMatching().CMDF1("C",2)[0].getCoeff(NLO)) <<  std::endl;    
//    std::cout << *(mySM.getMatching().CMDF1("C",2)[0].getCoeff(NNLO)) <<  std::endl;    

//    gslpp::matrix<gslpp::complex> myVCKM(mySM.getVCKM());
     sw = sqrt( (M_PI * mySM.getAle() ) / ( sqrt(2.) * mySM.getGF() * mySM.Mw() * mySM.Mw() ) );
//   sw = sqrt(mySM.sW2());
//    double as5 =  mySM.Alstilde5(mySM.getMuw());
//    double ae5 = mySM.getAle() / 4. / M_PI;
     as5 =  mySM.Als(120.,FULLNNNLO, true) / 4. / M_PI;
     ae5 = mySM.Ale(120.,FULLNLO) / 4. / M_PI;

    std::cout << std::endl;
    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(mySM.getMatching().CMDF1("CPL",8)[0].getCoeff(LO_QED))) <<  std::endl;    
    std::cout << *(mySM.getMatching().CMbsmm()[0].getCoeff(LO_QED)) << std::endl;
    std::cout << std::endl;
    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(mySM.getMatching().CMDF1("CPL",8)[0].getCoeff(NLO_QED11))) <<  std::endl;    
    std::cout << ae5*(*(mySM.getMatching().CMbsmm()[0].getCoeff(NLO_QED11))) << std::endl;
    std::cout << std::endl;
    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(mySM.getMatching().CMDF1("CPL",8)[0].getCoeff(NLO_QED02))) <<  std::endl;    
    std::cout << ae5*ae5/as5/as5*(*(mySM.getMatching().CMbsmm()[0].getCoeff(NLO_QED02))) << std::endl;
    std::cout << std::endl;
    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(mySM.getMatching().CMDF1("CPL",8)[0].getCoeff(NLO_QED21))) <<  std::endl;    
    std::cout << as5*ae5*(*(mySM.getMatching().CMbsmm()[0].getCoeff(NLO_QED21))) << std::endl;
    std::cout << std::endl;
    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(mySM.getMatching().CMDF1("CPL",8)[0].getCoeff(NLO_QED12))) <<  std::endl;    
    std::cout << ae5*ae5/as5*(*(mySM.getMatching().CMbsmm()[0].getCoeff(NLO_QED12))) << std::endl;
    std::cout << std::endl;
    std::cout << sw*sw*(myVCKM(2,2).conjugate() * myVCKM(2,1))*(*(mySM.getMatching().CMDF1("CPL",8)[0].getCoeff(NLO_QED22))) <<  std::endl;    
    std::cout << ae5*ae5*(*(mySM.getMatching().CMbsmm()[0].getCoeff(NLO_QED22))) << std::endl;    

//test OK  (up to 10-14)
//    std::cout << Heff.getEvol().AnomalousDimension(10, 2, 3) - HDB1.getUBsmm().AnomalousDimension(10, 2, 3) << std::endl;
//    std::cout << Heff.getEvol().AnomalousDimension(20, 2, 3) - HDB1.getUBsmm().AnomalousDimension(20, 2, 3) << std::endl;
//    std::cout << Heff.getEvol().AnomalousDimension(30, 2, 3) - HDB1.getUBsmm().AnomalousDimension(30, 2, 3) << std::endl;          
//    std::cout << Heff.getEvol().AnomalousDimension(01, 2, 3) - HDB1.getUBsmm().AnomalousDimension(01,2,3) <<  std::endl;
//    std::cout << Heff.getEvol().AnomalousDimension(02, 2, 3) - HDB1.getUBsmm().AnomalousDimension(02,2,3) <<  std::endl;
//    std::cout << Heff.getEvol().AnomalousDimension(11, 2, 3) - HDB1.getUBsmm().AnomalousDimension(11,2,3) <<  std::endl;
//    std::cout << Heff.getEvol().AnomalousDimension(21, 2, 3) - HDB1.getUBsmm().AnomalousDimension(21,2,3) <<  std::endl;

//
//      std::cout << Heff.getEvol().DF1Evol(5., 90., LO) << std::endl;
////    std::cout << HDB1.getUDF1BMll().Df1EvolMll(5., 90., LO) << std::endl;
//      std::cout << HDB1.getUDB1bsg().Df1Evolbsg(5., 90., LO) << std::endl;
//
//    std::cout << Heff.getEvol().DF1Evol(5., 90., NLO, NO_QED) << std::endl;
////    std::cout << HDB1.getUDF1BMll().Df1EvolMll(5., 90., NLO) << std::endl;
//    std::cout << HDB1.getUDB1bsg().Df1Evolbsg(5., 90., NLO) << std::endl;
//    
//    std::cout << Heff.getEvol().DF1Evol(5., 90., NNLO, NO_QED) << std::endl;
////    std::cout << HDB1.getUDF1BMll().Df1EvolMll(5., 90., NLO) << std::endl;
//    std::cout << HDB1.getUDB1bsg().Df1Evolbsg(5., 90., NNLO) << std::endl;
//
////    std::cout << Heff.getEvol().DF1Evol(2., 90., LO, NO_QED) << std::endl;
////    std::cout << HDB1.getUDF1BMll().Df1EvolMll(2., 90., LO) << std::endl;
////    std::cout << HDB1.getUDB1bsg().Df1Evolbsg(5., 90., LO) << std::endl;
////    std::cout << Heff.getEvol().DF1Evol(2., 90., NLO, NO_QED) << std::endl;
////    std::cout << HDB1.getUDF1BMll().Df1EvolMll(2., 90., NLO) << std::endl;
//    
////    std::cout << Heff.getEvol().DF1Evol(1., 90., LO, NO_QED) << std::endl;
////    std::cout << HDB1.getUDF1BMll().Df1EvolMll(1., 90., LO) << std::endl;
////    std::cout << HDB1.getUDB1bsg().Df1Evolbsg(5., 90., LO) << std::endl;
////    std::cout << Heff.getEvol().DF1Evol(1., 90., NLO, NO_QED) << std::endl;
////    std::cout << HDB1.getUDF1BMll().Df1EvolMll(1., 90., NLO) << std::endl;

    std::cout << std::endl << "LO:";
    std::cout << Heff.getEvol().DF1Evol(mub, 120., LO) << std::endl;
    std::cout << HDB1.getUBsmm().Df1Evol(mub, 120., LO, NO_QED) << std::endl;
    
    std::cout << std::endl << "NLO:";
    std::cout << Heff.getEvol().DF1Evol(mub, 120., NLO) / as5 << std::endl;
    std::cout << HDB1.getUBsmm().Df1Evol(mub, 120., NLO, NO_QED) << std::endl;
    
    std::cout << std::endl << "NNLO:";
    std::cout << Heff.getEvol().DF1Evol(mub, 120., NNLO) / (as5*as5) << std::endl;
    std::cout << HDB1.getUBsmm().Df1Evol(mub, 120., NNLO, NO_QED) << std::endl;
    
    std::cout << std::endl << "LO_QED:";
    std::cout << Heff.getEvol().DF1Evol(mub, 120., LO_QED) / (ae5/as5) << std::endl;
    std::cout << HDB1.getUBsmm().Df1Evol(mub, 120., NNLO, LO_QED) << std::endl;
    
    std::cout << std::endl << "NLO_QED11:";
    std::cout << Heff.getEvol().DF1Evol(mub, 120., NLO_QED11) / (ae5) << std::endl;
    std::cout << HDB1.getUBsmm().Df1Evol(mub, 120., NNLO, NLO_QED11) << std::endl;
    
    std::cout << std::endl << "NLO_QED21: <----- NO MATCH";
    std::cout << Heff.getEvol().DF1Evol(mub, 120., NLO_QED21) / (as5*ae5) << std::endl;
    std::cout << HDB1.getUBsmm().Df1Evol(mub, 120., NNLO, NLO_QED21) << std::endl;
    
    std::cout << std::endl << "NLO_QED02:";
    std::cout << Heff.getEvol().DF1Evol(mub, 120., NLO_QED02) / (ae5*ae5/as5/as5) << std::endl;
    std::cout << HDB1.getUBsmm().Df1Evol(mub, 120., NNLO, NLO_QED02) << std::endl;
    
    std::cout << std::endl << "NLO_QED12:";
    std::cout << Heff.getEvol().DF1Evol(mub, 120., NLO_QED12) / (ae5*ae5/as5) << std::endl;
    std::cout << HDB1.getUBsmm().Df1Evol(mub, 120., NNLO, NLO_QED12) << std::endl;
    
    std::cout << std::endl << "NLO_QED22:";
    std::cout << Heff.getEvol().DF1Evol(mub, 120., NLO_QED22) / (ae5*ae5) << std::endl;
    std::cout << HDB1.getUBsmm().Df1Evol(mub, 120., NNLO, NLO_QED22) << std::endl;
    
//    double eta =  as5 * 4. * M_PI / mySM.Als(mub,FULLNNNLO, true);
////    double eta = HDB1.getUBsmm().eta;
//    std::cout << std::endl << "eta: " << eta << std::endl;
//    std::cout << "eta_smm: " << HDB1.getUBsmm().eta << std::endl;
//    std::cout << "-2,1: " << Heff.getEvol().f_g(0,1,1,1,-2,1,eta) << std::endl;
//    std::cout << "-2,1: " << HDB1.getUBsmm().G(1,1,1,1,4,5.,120.,5.) << std::endl;
//    std::cout << "0,-1: " << Heff.getEvol().f_g(0,1,1,1,0,-1,eta) << std::endl;
//    std::cout << "0,-1: " << HDB1.getUBsmm().G(1,1,1,3,2,5.,120.,5.) << std::endl;
//    std::cout << "-1,0: " << Heff.getEvol().f_g(0,1,1,1,-1,0,eta) << std::endl;
//    std::cout << "-1,0: " << HDB1.getUBsmm().G(1,1,1,2,3,5.,120.,5.) << std::endl;
}
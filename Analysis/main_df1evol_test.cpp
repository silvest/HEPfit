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

    HeffDF1 Heff(8, "CPM", mySM, NNLO, NO_QED);
    HeffDB1 HDB1(mySM);

    std::cout << "%SUITE_STARTING% Evolutor" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

//    std::cout << *(mySM.getMatching().CMDF1("C",2)[0].getCoeff(LO)) <<  std::endl;    
//    std::cout << *(mySM.getMatching().CMDF1("C",2)[0].getCoeff(NLO)) <<  std::endl;    
//    std::cout << *(mySM.getMatching().CMDF1("C",2)[0].getCoeff(NNLO)) <<  std::endl;    

    std::cout << *(mySM.getMatching().CMDF1("CPM",8)[0].getCoeff(LO)) - *(mySM.getMatching().CMbsg()[0].getCoeff(LO))  <<  std::endl;    
    std::cout << *(mySM.getMatching().CMDF1("CPM",8)[0].getCoeff(NLO)) - *(mySM.getMatching().CMbsg()[0].getCoeff(NLO)) <<  std::endl;    
    std::cout << *(mySM.getMatching().CMDF1("CPM",8)[0].getCoeff(NNLO)) - *(mySM.getMatching().CMbsg()[0].getCoeff(NNLO)) <<  std::endl;    
    
//    std::cout << Heff.getEvol().AnomalousDimension(30, 2, 3) - HDB1.getUDB1bsg().AnomalousDimension_M(NNLO,2,3) <<  std::endl;
//
//    std::cout << Heff.getEvol().DF1Evol(5., 90., LO, NO_QED) << std::endl;
////    std::cout << HDB1.getUDF1BMll().Df1EvolMll(5., 90., LO) << std::endl;
//    std::cout << HDB1.getUDB1bsg().Df1Evolbsg(5., 90., LO) << std::endl;
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
}
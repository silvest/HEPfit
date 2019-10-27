/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ModelFactory.h"
#include "NPSTU.h"
#include "NPSTUZbbbarLR.h"
#include "NPEpsilons.h"
#include "NPEpsilons_pureNP.h"
#include "NPZbbbar.h"
#include "NPZbbbarLinearized.h"
#include "NPSMEFTd6.h"
#include "NPSMEFT6dtopquark.h"
#include "HiggsKvKf.h"
#include "HiggsKvKfgen.h" //
#include "HiggsKvgenKfgen.h"
#include "HiggsKigen.h"
#include "HiggsChiral.h"
#include "SUSY.h"
#include "THDM.h"
#include "FlavourWilsonCoefficient.h"
#include "FlavourWilsonCoefficient_DF2.h"
#include "RealWeakEFTLFV.h"
#include "RealWeakEFTCC.h"
#include "LoopMediators.h"
#include "NPDF2.h"
#include "CMFV.h"
#include "GeorgiMachacek.h"
#include "GeneralTHDM.h"
#include "THDMW.h"
#include <boost/bind.hpp>

ModelFactory::ModelFactory()
{
    modelFactory["StandardModel"] = boost::factory<StandardModel*>();
    modelFactory["NPSTU"] = boost::factory<NPSTU*>();
    modelFactory["NPSTUZbbbarLR"] = boost::factory<NPSTUZbbbarLR*>();
    modelFactory["NPEpsilons"] = boost::factory<NPEpsilons*>();
    modelFactory["NPEpsilons_pureNP"] = boost::factory<NPEpsilons_pureNP*>();
    modelFactory["NPZbbbar"] = boost::bind(boost::factory<NPZbbbar*>(), false);
    modelFactory["NPZbbbarLR"] = boost::bind(boost::factory<NPZbbbar*>(), true);
    modelFactory["NPZbbbarLinearized"] = boost::bind(boost::factory<NPZbbbarLinearized*>(), false);
    modelFactory["NPZbbbarLinearizedLR"] = boost::bind(boost::factory<NPZbbbarLinearized*>(), true);
    modelFactory["NPSMEFTd6"] = boost::bind(boost::factory<NPSMEFTd6*>(), false, false);
    modelFactory["NPSMEFT6dtopquark"] = boost::factory<NPSMEFT6dtopquark*>();
    modelFactory["NPSMEFTd6_LFU_QFU"] = boost::bind(boost::factory<NPSMEFTd6*>(), true, true);
    modelFactory["HiggsKvKf"] = boost::factory<HiggsKvKf*>();
    modelFactory["HiggsKvKfgen"] = boost::factory<HiggsKvKfgen*>();
    modelFactory["HiggsKvgenKfgen"] = boost::factory<HiggsKvgenKfgen*>();
    modelFactory["HiggsKigen"] = boost::factory<HiggsKigen*>();
    modelFactory["HiggsChiral"] = boost::factory<HiggsChiral*>();
    modelFactory["NPDF2"] = boost::factory<NPDF2*>();
    modelFactory["CMFV"] = boost::factory<CMFV*>();
    modelFactory["FlavourWilsonCoefficient"] = boost::factory<FlavourWilsonCoefficient*>();
    modelFactory["FlavourWilsonCoefficient_DF2"] = boost::factory<FlavourWilsonCoefficient_DF2*>();
    modelFactory["RealWeakEFTLFV"] = boost::factory<RealWeakEFTLFV*>();
    modelFactory["RealWeakEFTCC"] = boost::bind(boost::factory<RealWeakEFTCC*>(), 0);
    modelFactory["RealWeakEFTCCPM"] = boost::bind(boost::factory<RealWeakEFTCC*>(), 1);
    modelFactory["LoopMediators"] = boost::factory<LoopMediators*>();
    modelFactory["THDM"] = boost::factory<THDM*>();
    modelFactory["GeorgiMachacek"] = boost::factory<GeorgiMachacek*>();
    modelFactory["GeneralTHDM"] = boost::factory<GeneralTHDM*>();
    modelFactory["THDMW"] = boost::factory<THDMW*>();
}

void ModelFactory::addModelToFactory(const std::string name, boost::function<StandardModel*() > funct)
{
    modelFactory[name] = funct;
}

StandardModel* ModelFactory::CreateModel(const std::string& name)
{
    if (modelFactory.find(name) == modelFactory.end())
        throw std::runtime_error("ERROR: Wrong model " + name + " passed to ModelFactory.\n");
    return (modelFactory[name]());
}

/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ModelFactory.h"
#include "NPSTU.h"
#include "NPSTUZbbbarLR.h"
#include "NPEpsilons.h"
#include "NPEpsilons_pureNP.h"
#include "NPHiggs.h"
#include "NPZbbbar.h"
#include "NPZbbbarLinearized.h"
#include "NPEffectiveBS.h"
#include "NPEffectiveGIMR.h"
#include "HiggsKvKf.h"
#include "HiggsKvKfgen.h"
#include "HiggsKvgenKf.h"
/** BEGIN: REMOVE FROM THE PACKAGE **/
#include "pMSSM.h"
#include "MFV.h"
#include "NPSTUVWXY.h"
#include "GeneralSUSY.h"
#include "SUSYMassInsertion.h"
#include "SUSY.h"
#include "THDM.h"
#include "FlavourWilsonCoefficient.h"
/** END: REMOVE FROM THE PACKAGE **/
#include <boost/bind.hpp>

ModelFactory::ModelFactory()
{
    modelFactory["StandardModel"] = boost::factory<StandardModel*>();
    modelFactory["NPSTU"] = boost::factory<NPSTU*>();
    modelFactory["NPSTUZbbbarLR"] = boost::factory<NPSTUZbbbarLR*>();
    modelFactory["NPEpsilons"] = boost::factory<NPEpsilons*>();
    modelFactory["NPEpsilons_pureNP"] = boost::factory<NPEpsilons_pureNP*>();
    modelFactory["NPHiggs"] = boost::factory<NPHiggs*>();
    modelFactory["NPZbbbar"] = boost::bind(boost::factory<NPZbbbar*>(), false);
    modelFactory["NPZbbbarLR"] = boost::bind(boost::factory<NPZbbbar*>(), true);
    modelFactory["NPZbbbarLinearized"] = boost::bind(boost::factory<NPZbbbarLinearized*>(), false);
    modelFactory["NPZbbbarLinearizedLR"] = boost::bind(boost::factory<NPZbbbarLinearized*>(), true);
    modelFactory["NPEffectiveBS"] = boost::bind(boost::factory<NPEffectiveBS*>(), false, false);
    modelFactory["NPEffectiveBS_LFU"] = boost::bind(boost::factory<NPEffectiveBS*>(), true, false);
    modelFactory["NPEffectiveBS_QFU"] = boost::bind(boost::factory<NPEffectiveBS*>(), false, true);
    modelFactory["NPEffectiveBS_LFU_QFU"] = boost::bind(boost::factory<NPEffectiveBS*>(), true, true);
    modelFactory["NPEffectiveGIMR"] = boost::bind(boost::factory<NPEffectiveGIMR*>(), false, false);
    //modelFactory["NPEffectiveGIMR_LFU"] = boost::bind(boost::factory<NPEffectiveGIMR*>(), true, false);
    //modelFactory["NPEffectiveGIMR_QFU"] = boost::bind(boost::factory<NPEffectiveGIMR*>(), false, true);
    modelFactory["NPEffectiveGIMR_LFU_QFU"] = boost::bind(boost::factory<NPEffectiveGIMR*>(), true, true);
    modelFactory["HiggsKvKf"] = boost::factory<HiggsKvKf*>();
    modelFactory["HiggsKvKfgen"] = boost::factory<HiggsKvKfgen*>();
    modelFactory["HiggsKvgenKf"] = boost::factory<HiggsKvgenKf*>();
    /** BEGIN: REMOVE FROM THE PACKAGE **/
    modelFactory["MFV"] = boost::factory<MFV*>();
    modelFactory["pMSSM"] = boost::factory<pMSSM*>();
    modelFactory["NPSTUVWXY"] = boost::factory<NPSTUVWXY*>();
    modelFactory["GeneralSUSY"] = boost::factory<GeneralSUSY*>();
    modelFactory["SUSYMassInsertion"] = boost::factory<SUSYMassInsertion*>();
    modelFactory["THDM"] = boost::factory<THDM*>();
    modelFactory["FlavourWilsonCoefficient"] = boost::factory<FlavourWilsonCoefficient*>();
    /** END: REMOVE FROM THE PACKAGE **/
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

#include "ModelFactory.h"
#include <boost/bind.hpp>
#include <NPSTU.h>
#include <NPSTUVWXY.h>
#include <NPEpsilons.h>
#include <NPEpsilons_pureNP.h>
#include <NPHiggs.h>
#include <NPZbbbar.h>
#include <NPZbbbarLinearized.h>
#include <NPEffectiveBS.h>
#include <NPEffectiveGIMR.h>
#include <GeneralSUSY.h>
#include <pMSSM.h>
#include <SUSYMassInsertion.h>
#include <MFV.h>
#include <SUSY.h>
#include <THDM.h>
#include <HiggsKvKf.h>
#include <HiggsKvKfgen.h>
#include <HiggsKvgenKf.h>

ModelFactory::ModelFactory()
{
    modelFactory["StandardModel"] = boost::factory<StandardModel*>();
    modelFactory["NPSTU"] = boost::factory<NPSTU*>();
    modelFactory["NPSTUVWXY"] = boost::factory<NPSTUVWXY*>();
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
    modelFactory["MFV"] = boost::factory<MFV*>();
    modelFactory["GeneralSUSY"] = boost::factory<GeneralSUSY*>();
    modelFactory["pMSSM"] = boost::factory<pMSSM*>();
    modelFactory["SUSYMassInsertion"] = boost::factory<SUSYMassInsertion*>();
    modelFactory["THDM"] = boost::factory<THDM*>();
    modelFactory["HiggsKvKf"] = boost::factory<HiggsKvKf*>();
    modelFactory["HiggsKvKfgen"] = boost::factory<HiggsKvKfgen*>();
    modelFactory["HiggsKvgenKf"] = boost::factory<HiggsKvgenKf*>();
}

void ModelFactory::addModelToFactory(const std::string name, boost::function<StandardModel*() > funct)
{
    modelFactory[name] = funct;
}

StandardModel* ModelFactory::CreateModel(const std::string& name)
{
    if (modelFactory.find(name) == modelFactory.end())
        throw std::runtime_error("ERROR: Wrong model " + name + " passed to ModelFactory");
    return (modelFactory[name]());
}

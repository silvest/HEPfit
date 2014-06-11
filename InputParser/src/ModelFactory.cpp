#include "ModelFactory.h"
#include <boost/bind.hpp>
#include <NPSTU.h>
#include <NPSTUVWXY.h>
#include <NPEpsilons.h>
#include <NPEpsilons_pureNP.h>
#include <NPHiggs.h>
#include <NPZbbbar.h>
#include <NPZbbbarLinearized.h>
#include <NPEffective1.h>
#include <NPEffective2.h>
#include <GeneralSUSY.h>
#include <pMSSM.h>
#include <SUSYMassInsertion.h>
#include <MFV.h>
#include <SUSY.h>
#include <THDM.h>
#include <HiggsKvKf.h>

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
    modelFactory["NPEffective1"] = boost::factory<NPEffective1*>();
    modelFactory["NPEffective2"] = boost::factory<NPEffective2*>();
    modelFactory["MFV"] = boost::factory<MFV*>();
    modelFactory["GeneralSUSY"] = boost::factory<GeneralSUSY*>();
    modelFactory["pMSSM"] = boost::factory<pMSSM*>();
    modelFactory["SUSYMassInsertion"] = boost::factory<SUSYMassInsertion*>();
    modelFactory["THDM"] = boost::factory<THDM*>();
    modelFactory["HiggsKvKf"] = boost::factory<HiggsKvKf*>();
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

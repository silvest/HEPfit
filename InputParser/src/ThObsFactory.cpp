/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ThObsFactory.h"
#include "ThObservable.h"
#include "ParamObs.h"

ThObsFactory::ThObsFactory()
{
    registerEWObservables();
    registerDiBosonObservables();
    registerCouplingObservables();
    registerHiggsObservables();
    registerHiggsLeptonObservables();
    registerTopQuarkObservables();
    registerLEP2Observables();
    registerFlavourObservables();
    registerBSMObservables();
}

ThObsFactory::ThObsFactory(const ThObsFactory& orig)
{
    obsThFactory = orig.obsThFactory;
}

void ThObsFactory::addObsToFactory(const std::string name, std::function<ThObservable*(const StandardModel&) > funct)
{
    if (obsThFactory.find(name) == obsThFactory.end()) obsThFactory[name] = funct;
    else throw std::runtime_error("ERROR: Observable named: " + name + " already exists. Please give a different name to the user-defined observable " + name + ".");
}

ThObservable * ThObsFactory::CreateThMethod(const std::string& name, StandardModel& model) const
{
    if (model.isModelParam(name))
        return new ParamObs(model, name);
    if (obsThFactory.find(name) == obsThFactory.end())
        throw std::runtime_error("ERROR: Wrong observable " + name + " passed to ThObsFactory.\nIf " + name + " is a parameter that is specific to an observable, please list it after the observable in the configuration file.\n");
    ThObservable * myThObs = obsThFactory.at(name)(model);
    if (!myThObs->getParametersForObservable().empty()) model.addParameters(myThObs->getParametersForObservable());
    return (myThObs);
}

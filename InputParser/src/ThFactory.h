/* 
 * File:   ThFactory.h
 * Author: silvest
 *
 * Created on April 19, 2011, 12:23 PM
 */

#ifndef THFACTORY_H
#define	THFACTORY_H

#include <ThObservable.h>
#include <StandardModel.h>
#include <StandardModelMatching.h>
#include <FlavourObservables.h>
#include <EW.h>

class ThFactory {
public:
    ThFactory(const StandardModel& myModel, StandardModelMatching& myModelMatching);
    virtual ~ThFactory();
    ThObservable* getThMethod(const std::string& name);
private:
    std::map<std::string,ThObservable *> thobs;
    Flavour myFlavour;
    EW myEW;
};

#endif	/* THFACTORY_H */


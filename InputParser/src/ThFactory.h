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
#include <Flavour.h>
#include <EW.h>
#include <ZFitter.h>

class ThFactory {
public:
    ThFactory(const StandardModel& myModel);
    virtual ~ThFactory();
    ThObservable* getThMethod(const std::string& name);
private:
    std::map<std::string,ThObservable *> thobs;
    Flavour myFlavour;
    EW myEW;
    ZFitter myZFitter;
};

#endif	/* THFACTORY_H */

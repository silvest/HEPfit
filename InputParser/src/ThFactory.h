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

class ThFactory {
public:
    ThFactory(StandardModel * myModel);
    ThFactory(const ThFactory& orig);
    virtual ~ThFactory();
    ThObservable * getThMethod(const std::string& name);
private:
    std::map<std::string,ThObservable *> thobs;
};

#endif	/* THFACTORY_H */


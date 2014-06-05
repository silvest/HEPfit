/* 
 * Copyright (C) 2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef MODELFACTORY_H
#define	MODELFACTORY_H

#include <boost/functional/factory.hpp>
#include <boost/function.hpp>
#include <StandardModel.h>
#include <map>


class ModelFactory {
public:
    ModelFactory();
    ModelFactory(const ModelFactory& orig);
    virtual ~ModelFactory(){};
    
    void addModelToFactory (const std::string name, boost::function<StandardModel*() >);
    
    StandardModel* getModel(const std::string& ModelName);
private:
    std::map<std::string, boost::function<StandardModel* ()> > modelFactory;

};

#endif	/* MODELFACTORY_H */


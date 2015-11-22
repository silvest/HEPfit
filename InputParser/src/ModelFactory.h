/* 
 * Copyright (C) 2014 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */


#ifndef MODELFACTORY_H
#define	MODELFACTORY_H

#include "StandardModel.h"
#include <boost/functional/factory.hpp>
#include <boost/function.hpp>
#include <map>

/**
 * @class ModelFactory
 * @ingroup InputParser 
 * @brief A class for 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details
 */
class ModelFactory {
public:
    ModelFactory();
    ModelFactory(const ModelFactory& orig);
    virtual ~ModelFactory(){};
    
    void addModelToFactory (const std::string name, boost::function<StandardModel*() >);
    
    StandardModel* CreateModel(const std::string& ModelName);
private:
    std::map<std::string, boost::function<StandardModel* ()> > modelFactory;

};

#endif	/* MODELFACTORY_H */


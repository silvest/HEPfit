/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MODELPARAMETER_H
#define	MODELPARAMETER_H

#include <string>
#include <iostream>

/**
 * @class ModelParameter
 * @ingroup Observable
 * @brief A class for model parameters. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class ModelParameter {
public:
    ModelParameter(std::string, double, double, double);
    virtual ~ModelParameter();
    std::string name;
    double ave,errg,errf,min,max;
    friend std::ostream& operator<<(std::ostream& output, const ModelParameter& m);
private:
};

#endif	/* MODELPARAMETER_H */


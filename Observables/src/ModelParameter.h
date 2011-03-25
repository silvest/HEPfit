/* 
 * File:   ModelParameter.h
 * Author: silvest
 *
 * Created on March 15, 2011, 2:59 PM
 */

#ifndef MODELPARAMETER_H
#define	MODELPARAMETER_H

#include <string>
#include <iostream>

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


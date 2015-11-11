/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ModelParameter.h"

ModelParameter::ModelParameter(std::string name_in, double ave_in, double errg_in, double errf_in) 
{
    name = name_in;
    ave = ave_in;
    errg = errg_in;
    errf = errf_in;
    min = ave - errf - 5.*errg;
    max = ave + errf + 5.*errg;
    cgp_name.clear();
    if(errg == 0. && errf == 0.)
        isFixed = true;
    else
        isFixed = false;
}

ModelParameter::ModelParameter() 
{}

ModelParameter::~ModelParameter() 
{
}

std::ostream& operator<<(std::ostream& output, const ModelParameter& m)
{
    output << "ModelParameter name, average, gaussian error, flat error" << std::endl;
    output << m.name << " " << m.ave << " " << m.errg << " " << m.errf << std::endl;
    return output;
}

boost::tokenizer<boost::char_separator<char> >::iterator & ModelParameter::ParseModelParameter(boost::tokenizer<boost::char_separator<char> >::iterator & beg)
{

    name = *beg;
    ++beg;
    ave = atof((*beg).c_str());
    ++beg;
    errg = atof((*beg).c_str());
    ++beg;
    errf = atof((*beg).c_str());
    ++beg;
    min = ave - errf - 5.*errg;
    max = ave + errf + 5.*errg;
    
    cgp_name.clear();
    if(errg == 0. && errf == 0.)
        isFixed = true;
    else
        isFixed = false;
    
    return beg;
}
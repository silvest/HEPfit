/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CLAUSENFUNCTIONS_H
#define	CLAUSENFUNCTIONS_H

#include <gslpp.h>
#include "BernoulliNumbers.h"
using namespace gslpp;

/**
 * @class ClausenFunctions
 * @ingroup LoopFunctions 
 * @brief A class for Clausen functions. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class ClausenFunctions : public BernoulliNumbers {
public:
    
    /**
     * @brief  ClausenFunctions constructor. 
     */
    ClausenFunctions();
    
    ////////////////////////////////////////////////////////////////////////

    /**
     * @brief Clausen function of index 2.
     * @param[in] phi Real variable.
     * @return Clausen function of index 2, Cl_2(phi).
     */
    double Cl2(const double phi) const;    
    
    /**
     * @brief Clausen function of index 3.
     * @param[in] phi Real variable. 
     * @return Clausen function of index 3, Cl_3(phi). 
     * @attention applicable for 0<=phi<=Pi. 
     */
    double Cl3(const double phi) const;
        
    ////////////////////////////////////////////////////////////////////////    
    
private:
  
};

#endif	/* CLAUSENFUNCTIONS_H */


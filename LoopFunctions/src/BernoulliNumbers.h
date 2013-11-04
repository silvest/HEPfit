/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BERNOULLINUMBERS_H
#define	BERNOULLINUMBERS_H

/**
 * @addtogroup LoopFunctions
 * @brief A project for loop functions. 
 * @{
 */

/**
 * @class BernoulliNumbers
 * @brief A class for Bernoulli numbers. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class BernoulliNumbers {
public:

    /**
     * @brief BernoulliNumbers constructor. 
     */
    BernoulliNumbers();
    
    ////////////////////////////////////////////////////////////////////////

protected:
    double B[19]; /* Bernoulli numbers */
        
    ////////////////////////////////////////////////////////////////////////
    
private:

};

/** 
 * @}
 */

#endif	/* BERNOULLINUMBERS_H */


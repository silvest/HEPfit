/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef BERNOULLINUMBERS_H
#define	BERNOULLINUMBERS_H

/**
 * @addtogroup LoopFunctions
 * @brief A module for loop functions.
 * @{
 */

/**
 * @class BernoulliNumbers
 * @brief A class for the Bernoulli numbers.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details This class handles the Bernoulli numbers, which are used to
 * compute the trilogarithm function Polylogarithms::Li3() and
 * the Clausen function of index three ClausenFunctions::Cl3.
 * A list of the Bernoulli numbers can be found in @cite tHooft:1978xw.
 */
class BernoulliNumbers {
public:

    /**
     * @brief The default constructor.
     */
    BernoulliNumbers();
    
    ////////////////////////////////////////////////////////////////////////
protected:
    double B[19];///< the Bernoulli numbers

};

/** 
 * @}
 */

#endif	/* BERNOULLINUMBERS_H */


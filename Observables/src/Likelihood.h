/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LIKELIHOOD_H
#define	LIKELIHOOD_H

/**
 * @class Likelihood
 * @ingroup Observable
 * @brief A class for likelihood function (Not used). 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details  
 */
class Likelihood {
public:
    
    /**
     * @brief The default constructor.
     */
    Likelihood();
    
    /**
     * @brief The copy constructor.
     */
    Likelihood(const Likelihood& orig);
    
    /**
     * @brief The default destructor.
     */
    virtual ~Likelihood();
    
    /**
     * @brief
     * @param
     */
    double getLikelihood(const double) const;
    
    /**
     * @brief
     * @param
     */
    double getLikelihood(const double, const double) const;
private:

};

#endif	/* LIKELIHOOD_H */


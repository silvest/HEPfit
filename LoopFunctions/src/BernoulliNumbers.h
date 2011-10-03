/* 
 * File:   BernoulliNumbers.h
 * Author: mishima
 */

#ifndef BERNOULLINUMBERS_H
#define	BERNOULLINUMBERS_H


class BernoulliNumbers {
public:

    /**
     * @brief BernoulliNumbers constructor
     */
    BernoulliNumbers();

    //BernoulliNumbers(const BernoulliNumbers& orig);

    /**
     * @brief BernoulliNumbers destructor
     */
    virtual ~BernoulliNumbers();

    ////////////////////////////////////////////////////////////////////////

protected:
    double B[19]; /* Bernoulli numbers */
    
    
    ////////////////////////////////////////////////////////////////////////
    
private:

    
};

#endif	/* BERNOULLINUMBERS_H */


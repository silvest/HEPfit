/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RGEVOLUTOR_H
#define	RGEVOLUTOR_H

#include <gslpp.h>
#include "OrderScheme.h"
#include "WilsonTemplate.h"

/**
 * @class RGEvolutor
 * @ingroup StandardModel
 * @brief A class for the RG evolutor of the Wilson coefficients. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */

class RGEvolutor : public WilsonTemplate<gslpp::matrix<double> > {
public:
    
    /**
     *
     * @brief constructor
     * @param[in] dim dimension of the operator basis
     * @param[in] scheme renormalizations scheme
     * @param[in] order order of QCD coupling
     */
    RGEvolutor(unsigned int dim, schemes scheme, orders order);
    
    /**
     *
     * @brief constructor
     * @param[in] dim dimension of the operator basis
     * @param[in] scheme renormalizations scheme
     * @param[in] order order of QCD coupling
     * @param[in] order_ew order of Electroweak coupling
     */
    RGEvolutor(unsigned int dim, schemes scheme, orders order, orders_ew order_ew);
    
    /**
     *
     * @brief destructor
     */
    virtual ~RGEvolutor();
    
    /**
     *
     * @brief
     * @param[in] i
     * @param[in] j
     * @param[in] x
     * @param[in] order_i order of QCD coupling
     */
    void setEvol(unsigned int i, unsigned int j, double x, orders order_i);
    
    /**
     *
     * @brief
     * @param[in] i
     * @param[in] j
     * @param[in] x
     * @param[in] order_i order of QCD coupling
     * @param[in] order_ew order of Electroweak coupling
     */
    void setEvol(unsigned int i, unsigned int j, double x, orders order_i, orders_ew order_ew);
    
    /**
     *
     * @brief
     * @param[in] m Evolution matrix
     * @param[in] order_i order of QCD coupling
     */
    void setEvol(const gslpp::matrix<double>& m, orders order_i);
    
    /**
     *
     * @brief
     * @param[in] m
     * @param[in] order_ew_i order of Electroweak coupling
     */
    void setEvol(const gslpp::matrix<double>& m, orders_ew order_ew_i);
    
    /**
     *
     * @brief
     * @return
     */
    gslpp::matrix<double>** getEvol() const;

    /**
     *
     * @brief Retrieve the upper scale of the Wilson Coefficients
     * @return M The scale of the Wilson Coefficients set by the model
     */
    double getM() const;

    /**
     *
     * @brief Sets the upper and lower scale for the running of the Wilson Coefficients
     * @param[in] mu Lower RGE running scale
     * @param[in] M Upper RGE running scale
     */
    void setScales(double mu, double M);

    /**
     *
     * @brief Sets the upper scale for the running of the Wilson Coefficients
     * @param[in] M Upper RGE running scale
     */
    void setM(double M);
    
    /**
     *
     * @brief Sets the lower scale for the running of the Wilson Coefficients
     * @param[in] mu Lower RGE running scale
     */
    void setMu(double mu);

    /**
     *
     * @brief Evolution matrix set at a fixed order of QCD coupling
     * @param[in] order order of QCD coupling
     * @return The RGE evolution matrix at a fixed order of QCD coupling
     */
    gslpp::matrix<double>* Evol(orders order);
    
    /**
     *
     * @brief Evolution matrix set at a fixed order of Electroweak coupling
     * @param[in] order_ew order of Electroweak coupling
     * @return The RGE evolution matrix at a fixed order of Electroweak coupling
     */
    gslpp::matrix<double>* Evol(orders_ew order_ew);
    
protected:
    double M;
};

#endif	/* RGEVOLUTOR_H */


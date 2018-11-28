/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef RGEVOLUTORNEW_H
#define	RGEVOLUTORNEW_H

#include <gslpp.h>
#include "OrderScheme.h"
#include "WilsonTemplateNew.h"

/**
 * @class RGEvolutor
 * @ingroup StandardModel
 * @brief A class for the RG evolutor of the Wilson coefficients. 
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 * @details 
 */

class RGEvolutorNew : public WilsonTemplateNew<gslpp::matrix<double> > {
public:
    
    /**
     *
     * @brief constructor
     * @param[in] dim dimension of the operator basis
     * @param[in] scheme renormalizations scheme
     * @param[in] order order of QCD coupling
     * @param[in] order_qed order of Electroweak coupling
     */
    RGEvolutorNew(unsigned int dim, schemes scheme, qcd_orders order_qcd_i, qed_orders order_qed_i = QED0);
    
    /**
     *
     * @brief destructor
     */
    virtual ~RGEvolutorNew() {};
    
    /**
     *
     * @brief
     * @param[in] i
     * @param[in] j
     * @param[in] x
     * @param[in] order_i order of QCD coupling
     * @param[in] order_qed order of Electroweak coupling
     */
    void setEvol(unsigned int i, unsigned int j, double x, qcd_orders order_i, qed_orders order_qed = QED0);
    
    /**
     *
     * @brief
     * @param[in] m
     * @param[in] order_qed_i order of Electroweak coupling
     */
    void setEvol(const gslpp::matrix<double>& m, qcd_orders order_qcd_i, qed_orders order_qed_i = QED0);
    
    /**
     *
     * @brief
     * @return
     */
    const Expanded<gslpp::matrix<double> >& getEvol() const;

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
     * @brief Evolution matrix set at a fixed order of Electroweak coupling
     * @param[in] order_qed order of Electroweak coupling
     * @return The RGE evolution matrix at a fixed order of Electroweak coupling
     */
    const gslpp::matrix<double>& Evol(qcd_orders order_qcd, qed_orders order_qed = QED0) const;
    
protected:
    double M;
};

#endif	/* RGEVOLUTOR_H */


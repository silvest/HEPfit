/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef CPENGUINBOX_H
#define	CPENGUINBOX_H


#include <StandardModel.h>
#include <sstream>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_clausen.h>

/**
 * 
 * @class: CPenguinBox
 * @details: compute the charm contributions to the rare decays of K+ meson
 * including matching at top scale, RG evolution to the charm scale and 
 * matching at the charm scale. hep-ph/0603079
 */
class CPenguinBox {
    
public:
    
    /**
     * 
     * @brief constructor
     */
    CPenguinBox(const StandardModel& model_i);
    
    /**
     * 
     * @brief destructor
     */
    ~CPenguinBox();
    
    /**
     * 
     * @param order, perturbation theory order
     * @param nf, number of flavours
     * @return the LO RG evolution matrix for the Z-penguin contribution
     */
    gslpp::matrix<double> RGevolP(int nf);
    
    /**
     * 
     * @param order, perturbation theory order
     * @param nf, number of flavours
     * @return the LO RG evolution matrix for the Box contribution
     */
    gslpp::matrix<double> RGevolB(int nf);
    
    /**
     *  
     * @param x, x_t
     * @param order, perturbation theory order
     * @return matching Wilson coefficients at the top scale for the Z-penguin contribution
     */    
    gslpp::vector<double> WMatchP(orders order);
    
    /**
     * 
     * @param x, x_t
     * @param order, perturbation theory order
     * @return matching Wilson coefficients at the W scale for the Box contribution
     */    
    gslpp::vector<double> WMatchB(orders order);
    
    /**
     * 
     * @param order
     * @return return the scale-indipendent Wilson coefficient below the Chram scale
     */
    double Cmatch(orders order);
    
    /**
     * 
     * @brief NNLO QCD correction from  hep-ph/0603079, two loop EW corrections from hep-ph/1009.0947v2
     * @param order
     * @return chram contribuition times the model.getlam_c().real plus the 
     * isospin-breacking contribution.
     */
    double CT_tot(orders order, orders_ew order_ew);
    
protected:
    
private:
    const StandardModel& model;
    const StandardModelMatching& modelmatching; 
    
    gslpp::matrix<double> up, j5p1, j4p1, j5p2, j4p2;
    gslpp::matrix<double> ub, j5b1, j4b1, j5b2, j4b2;
    
    gslpp::matrix<double> up4, up5, ub4, ub5;
    
    gslpp::vector<double> p, b;
    
    gslpp::vector<double> p0, p1, p2;
    gslpp::vector<double> b0, b1, b2;
    
    gslpp::vector<double> r1, r2;
    
    gslpp::vector<double> CPM0, CPM1, CPM2;
    gslpp::vector<double> CBM0, CBM1, CBM2;
    
    double r1tau, r2tau, r1e, r2e;
    //double tbp, tbb;
    double eta;
    
    double cpm, cbmt, cbme;

    /**
     * 
     * @param order
     * @return matching at SM::getMub() threshold for the Z-penguin contribution
     */
    double BmatchP(orders order);
    
    /**
     * 
     * @param order
     * @return matching at SM::getMub() threshold for the Box contribution
     */
    double BmatchB(orders order);    
};

#endif	/* CPENGUINBOX_H */

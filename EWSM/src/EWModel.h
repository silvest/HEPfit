/* 
 * File:   EWModel.h
 * Author: mishima
 */

#ifndef EWMODEL_H
#define	EWMODEL_H

#include <gslpp.h>
using namespace gslpp;


class EWModel {
public:
    
    /**
     * @return the total radiative corrections to alpha at Mz
     */
    virtual double DeltaAlpha() {
        throw "EWModel::DeltaAlpha() is undefined.";
    };

    /**
     * @return Delta r
     */
    virtual double DeltaR() {
        throw "EWModel::DeltaR() is undefined.";
    };
    
    /**
     * @brief effective coupling g_V^l
     * @param[in] l name of a lepton 
     * @return g_V^l for lepton "l" 
     */
    virtual complex deltaGV_l(const StandardModel::lepton l) {
        throw "EWModel::deltaGV_l() is undefined.";
    };

    /**
     * @brief effective coupling g_V^q
     * @param[in] q name of a quark
     * @return g_V^q for quark "q" 
     */    
    virtual complex deltaGV_q(const StandardModel::quark q) {
        throw "EWModel::deltaGV_q() is undefined.";
    };
    
    /**
     * @brief effective coupling g_A^l
     * @param[in] l name of a lepton 
     * @return g_A^l for lepton "l" 
     */
    virtual complex deltaGA_l(const StandardModel::lepton l) {
        throw "EWModel::deltaGA_l() is undefined.";
    };
    
    /**
     * @brief effective coupling g_A^q
     * @param[in] q name of a quark
     * @return g_A^q for quark "q" 
     */
    virtual complex deltaGA_q(const StandardModel::quark q) {
        throw "EWModel::deltaGA_q() is undefined.";
    };

};

#endif	/* EWMODEL_H */


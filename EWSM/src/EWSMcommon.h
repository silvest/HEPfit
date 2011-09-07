/* 
 * File:   EWSMcommon.h
 * Author: mishima
 */

#ifndef EWSMCOMMON_H
#define	EWSMCOMMON_H

#include <StandardModel.h>


class EWSMcommon : public StandardModel {

public:

    /**
     * @brief EWSMcommon constructor
     * @param[in] SM_i reference to a StandardModel object
     */
    EWSMcommon(const StandardModel& SM_i);

    /**
     * @brief EWSMcommon copy constructor
     * @param[in] orig reference to an EWSMcommon object
     */
    //EWSMcommon(const EWSMcommon& orig);

    /**
     * @brief EWSMcommon destructor
     */
    virtual ~EWSMcommon();

    
    //////////////////////////////////////////////////////////////////////// 
    
    /**
     * @brief computes common variables
     */
    void SetConstants();

    /**
     * @ computes variables which depend on Mw
     */
    void Compute(const double Mw_i);

    
    //////////////////////////////////////////////////////////////////////// 
    
    /**
     * @return the W boson mass
     */
    double GetMW() const {
        return Mw;
    }

    /**
     * @return c_W^2
     */
    double GetCW2() const {
        return cW2;
    }

    /**
     * @return s_W^2
     */
    double GetSW2() const {
        return sW2;
    }
    
    /**
     * @return the conversion factor from alpha to GF
     */
    double GetF_AlphaToGF() const {
        return f_AlphaToGF;
    }
    
    /**
     * @return the zeta function zeta(2)
     */
    double GetZeta2() const {
        return zeta2;
    }

    /**
     * @return the zeta function zeta(3)
     */
    double GetZeta3() const {
        return zeta3;
    }

    /**
     * @return the zeta function zeta(5)
     */
    double GetZeta5() const {
        return zeta5;
    }

    /**
     * @return the logarithmic function log(2)
     */
    double GetLog2() const {
        return log2;
    }

    /**
     * @return the log of Mz/me
     */
    double GetLogMZtoME() const {
        return logMZtoME;
    }

    /**
     * @return the log of Mz/mmu
     */
    double GetLogMZtoMMU() const {
        return logMZtoMMU;
    }

    /**
     * @return the log of Mz/mtau
     */
    double GetLogMZtoMTAU() const {
        return logMZtoMTAU;
    }
    
    /**
     * @return the log of Mz/mtop
     */
    double GetLogMZtoMTOP() const {
        return logMZtoMTOP;
    }
    
    
    //////////////////////////////////////////////////////////////////////// 

protected:
    double Mw;
    double sW2, cW2;
    
    double f_AlphaToGF;
    
    double zeta2, zeta3, zeta5, log2;
    double logMZtoME, logMZtoMMU, logMZtoMTAU, logMZtoMTOP;    
    
};

#endif	/* EWSMCOMMON_H */


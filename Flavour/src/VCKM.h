/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef VCKM_H
#define	VCKM_H

#include "ThObservable.h"

class VCKM : public ThObservable {
public:
    VCKM(const StandardModel& SM_i, unsigned int obsFlag_1,  unsigned int obsFlag_2);

    double computeThValue();
    
    virtual ~VCKM();
    
private:
    unsigned int obs_1;
    unsigned int obs_2;
        
};

class CKM_Alpha : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    CKM_Alpha(const StandardModel& SM_i);
    
    /**
     * @return The CKM angle @f$\alpha @f$ in degrees
     */
    double computeThValue ();
      
};

class CKM_Beta : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_Beta(const StandardModel& SM_i);
    
    /**
     * @return The CKM angle @f$\beta @f$ in degrees
     */
    double computeThValue ();
      
};

class CKM_Betas : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_Betas(const StandardModel& SM_i);
    
    /**
     * @return The CKM angle @f$\beta_s @f$ in degrees
     */
    double computeThValue ();
      
};

class CKM_Gamma : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_Gamma(const StandardModel& SM_i);
    
    /**
     * @return The CKM angle @f$\gamma @f$ in degrees
     */
    double computeThValue ();
      
};

class CKM_2BpG : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_2BpG(const StandardModel& SM_i);
    
    /**
     * @return The sum of CKM angles @f$2 \beta + \gamma @f$ in degrees
     */
    double computeThValue ();
      
};

class CKM_S2Beta : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_S2Beta(const StandardModel& SM_i);
    
    /**
     * @return @f$\sin 2 \beta @f$
     */
    double computeThValue ();
      
};

class CKM_C2Beta : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_C2Beta(const StandardModel& SM_i);
    
    /**
     * @return @f$\cos 2 \beta @f$
     */
    double computeThValue ();
      
};

class CKM_SinTheta12 : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_SinTheta12(const StandardModel& SM_i);
    
    /**
     * @return @f$\sin \theta_{12} @f$ in the PDG parameterization of the CKM matrix
     */
    double computeThValue ();
      
};

class CKM_SinTheta13 : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_SinTheta13(const StandardModel& SM_i);
    
    /**
     * @return @f$\sin \theta_{13} @f$ in the PDG parameterization of the CKM matrix
     */
    double computeThValue ();
      
};

class CKM_SinTheta23 : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_SinTheta23(const StandardModel& SM_i);
    
    /**
     * @return @f$\sin \theta_{13} @f$ in the PDG parameterization of the CKM matrix
     */
    double computeThValue ();
      
};

class CKM_Delta : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_Delta(const StandardModel& SM_i);
    
    /**
     * @return the phase @f$\delta @f$ in the PDG parameterization of the CKM matrix
     */
    double computeThValue ();
      
};

class J_CP : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    J_CP(const StandardModel& SM_i);
    
    /**
     * @return the Jarlskog determinant
     */
    double computeThValue ();
      
};

class CKM_Rt : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_Rt(const StandardModel& SM_i);
    
    /**
     * @brief @f$R_t=|(V_{td} V_{tb}^*)/(V_{cd}V_{cb}^*)|@f$.
     * @return @f$R_t@f$
     */
    double computeThValue ();
      
};

class CKM_Rts : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_Rts(const StandardModel& SM_i);
    
    /**
     * @brief @f$R_{ts}=|(V_{ts} V_{tb}^*)/(V_{cs}V_{cb}^*)|@f$.
     * @return @f$R_{ts}@f$
     */
    double computeThValue ();
      
};

class CKM_Rb : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_Rb(const StandardModel& SM_i);
    
    /**
     * @brief @f$R_{b}=|(V_{ud} V_{ub}^*)/(V_{ud}V_{ub}^*)|@f$.
     * @return @f$R_{b}@f$
     */
    double computeThValue ();
      
};

class CKM_VtdoVts : public ThObservable {
public:
      
    /**
     * @brief Constructor.
     */
    CKM_VtdoVts(const StandardModel& SM_i);
    
    /**
     * @brief @f$|V_{td}/V_{ts}|@f$.
     * @return @f$|V_{td}/V_{ts}|@f$
     */
    double computeThValue ();
      
};


class Abslam_t : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Abslam_t(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{td} V_{ts}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

class Abslam_c : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Abslam_c(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{cd} V_{cs}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

class Abslam_u : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Abslam_u(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{ud} V_{us}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

class Abslam_td : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Abslam_td(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{td} V_{tb}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

class Abslam_cd : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Abslam_cd(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{cd} V_{cb}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

class Abslam_ud : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Abslam_ud(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{ud} V_{ub}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

class Abslam_ts : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Abslam_ts(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{ts} V_{tb}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

class Abslam_cs : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Abslam_cs(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{cs} V_{cb}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

class Abslam_us : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Abslam_us(const StandardModel& SM_i);
    
    /**
     * @return @f$\vert V_{us} V_{ub}^* \vert @f$
     */
    double computeThValue ();
    
private:
    
};

class Imlam_t : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Imlam_t(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{td} V_{ts}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Imlam_c : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Imlam_c(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{cd} V_{cs}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Imlam_u : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Imlam_u(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{ud} V_{us}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Imlam_td : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Imlam_td(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{td} V_{tb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Imlam_cd : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Imlam_cd(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{cd} V_{cb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Imlam_ud : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Imlam_ud(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{ud} V_{ub}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Imlam_ts : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Imlam_ts(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{ts} V_{tb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Imlam_cs : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Imlam_cs(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{cs} V_{cb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Imlam_us : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Imlam_us(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Im}( V_{us} V_{ub}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Relam_t : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Relam_t(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{td} V_{ts}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Relam_c : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Relam_c(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{cd} V_{cs}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Relam_u : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Relam_u(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{ud} V_{us}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Relam_td : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Relam_td(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{td} V_{tb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Relam_cd : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Relam_cd(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{cd} V_{cb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Relam_ud : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Relam_ud(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{ud} V_{ub}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Relam_ts : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Relam_ts(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{ts} V_{tb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Relam_cs : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Relam_cs(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{cs} V_{cb}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};

class Relam_us : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Relam_us(const StandardModel& SM_i);
    
    /**
     * @return @f$\mathrm{Re}( V_{us} V_{ub}^* ) @f$
     */
    double computeThValue ();
    
private:
    
};


#endif	/* VCKM_H */


/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef VCKM_H
#define	VCKM_H

#include <ThObservable.h>

class VCKM : public ThObservable {
public:
    VCKM(const StandardModel& SM_i, unsigned int obsFlag_1,  unsigned int obsFlag_2);

    double computeThValue();
    
    virtual ~VCKM();
    
private:
    unsigned int obs_1;
    unsigned int obs_2;
        
};

class Abslam_t : public ThObservable {
public:
    
    /**
     * @brief Constructor.
     */
    Abslam_t(const StandardModel& SM_i);
    
    /**
     * @return muVBFtata
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
     * @return muVBFtata
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
     * @return muVBFtata
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
     * @return muVBFtata
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
     * @return muVBFtata
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
     * @return muVBFtata
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
     * @return muVBFtata
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
     * @return muVBFtata
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
     * @return muVBFtata
     */
    double computeThValue ();
    
private:
    
};


#endif	/* VCKM_H */


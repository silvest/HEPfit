/*
 * Copyright (C) 2023 HEPfit Collaboration
 * 
 * 
 * For the licensing terms see doc/COPYING.
 */


#ifndef DGAMMA_H
#define DGAMMA_H
#include "ThObservable.h"
#include "AmpDB2.h"

class DGamma_d : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_d(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_d();
    
    double computeThValue();

    private:
};

class DGamma_s : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_s(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_s();
    
    double computeThValue();

    private:
};

class DGamma_d_NLO : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_d_NLO(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_d_NLO();
    
    double computeThValue();

    private:
};

class DGamma_s_NLO : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_s_NLO(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_s_NLO();
    
    double computeThValue();

    private:
};

class DGamma_d_NLO1 : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_d_NLO1(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_d_NLO1();
    
    double computeThValue();

    private:
};

class DGamma_s_NLO1 : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_s_NLO1(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_s_NLO1();
    
    double computeThValue();

    private:
};

#endif /* DGAMMA_H */


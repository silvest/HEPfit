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

class DGamma_d_pole : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_d_pole(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_d_pole();
    
    double computeThValue();

    private:
};

class DGamma_s_pole : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_s_pole(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_s_pole();
    
    double computeThValue();

    private:
};

class DGamma_d_MSbar : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_d_MSbar(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_d_MSbar();
    
    double computeThValue();

    private:
};

class DGamma_s_MSbar : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_s_MSbar(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_s_MSbar();
    
    double computeThValue();

    private:
};

class DGamma_d_PS : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_d_PS(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_d_PS();
    
    double computeThValue();

    private:
};

class DGamma_s_PS : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_s_PS(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_s_PS();
    
    double computeThValue();

    private:
};

class DGamma_d_MSbar_NLO : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_d_MSbar_NLO(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_d_MSbar_NLO();
    
    double computeThValue();

    private:
};

class DGamma_s_MSbar_NLO : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_s_MSbar_NLO(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_s_MSbar_NLO();
    
    double computeThValue();

    private:
};

class DGamma_d_PS_NLO : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_d_PS_NLO(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_d_PS_NLO();
    
    double computeThValue();

    private:
};

class DGamma_s_PS_NLO : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_s_PS_NLO(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_s_PS_NLO();
    
    double computeThValue();

    private:
};

class DGamma_d_MSbar_LO : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_d_MSbar_LO(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_d_MSbar_LO();
    
    double computeThValue();

    private:
};

class DGamma_s_MSbar_LO : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_s_MSbar_LO(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_s_MSbar_LO();
    
    double computeThValue();

    private:
};

class DGamma_d_PS_LO : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_d_PS_LO(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_d_PS_LO();
    
    double computeThValue();

    private:
};

class DGamma_s_PS_LO : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    DGamma_s_PS_LO(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~DGamma_s_PS_LO();
    
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


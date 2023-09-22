/* 
 * Copyright (C) 2023 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ASL_H
#define ASL_H
#include "ThObservable.h"
#include "AmpDB2.h"


class Asl_d_pole : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Asl_d_pole(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~Asl_d_pole();
    
    double computeThValue ();

    private:
};

class Asl_s_pole : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Asl_s_pole(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~Asl_s_pole();
    
    double computeThValue();

    private:
};

class Asl_d_MSbar : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Asl_d_MSbar(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~Asl_d_MSbar();
    
    double computeThValue ();

    private:
};

class Asl_s_MSbar : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Asl_s_MSbar(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~Asl_s_MSbar();
    
    double computeThValue();

    private:
};

class Asl_d_PS : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Asl_d_PS(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~Asl_d_PS();
    
    double computeThValue ();

    private:
};

class Asl_s_PS : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Asl_s_PS(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~Asl_s_PS();
    
    double computeThValue();

    private:
};

class Asl_d_MSbar_NLO : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Asl_d_MSbar_NLO(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~Asl_d_MSbar_NLO();
    
    double computeThValue ();

    private:
};

class Asl_s_MSbar_NLO : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Asl_s_MSbar_NLO(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~Asl_s_MSbar_NLO();
    
    double computeThValue();

    private:
};

class Asl_d_MSbar_NLO_tradBasis : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Asl_d_MSbar_NLO_tradBasis(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~Asl_d_MSbar_NLO_tradBasis();
    
    double computeThValue ();

    private:
};

class Asl_s_MSbar_NLO_tradBasis : public ThObservable, AmpDB2{
public:
    /**
     * @brief Constructor.
     * @param[in] SM_i a reference to an object of type StandardModel
     */
    Asl_s_MSbar_NLO_tradBasis(const StandardModel& SM_i);
    
    /**
     * @brief Destructor.
     */
    virtual ~Asl_s_MSbar_NLO_tradBasis();
    
    double computeThValue();

    private:
};
 
#endif /* ASL_H */


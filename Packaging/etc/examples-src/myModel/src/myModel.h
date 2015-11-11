/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MYMODEL_H
#define	MYMODEL_H

#include <HEPfit.h>

/**
 * @class myModel
 * @brief My own Model.
 */
class myModel: public StandardModel {
public:

    static const int NmyModelvars = 4; /* Define number of mandatory parameters in the model. */
    static const std::string myModelvars[NmyModelvars]; /* Vector of model variable names. */
    
    /**
     * @brief myModel constructor
     */
    myModel();
    
    /**
     * @brief myModel destructor
     */
    ~myModel();
    
    virtual bool InitializeModel();
    
    virtual bool Init(const std::map<std::string, double>& DPars);
    
    virtual bool PreUpdate();
    
    virtual bool Update(const std::map<std::string, double>& DPars);
    
    virtual bool PostUpdate();
    
    virtual bool CheckParameters(const std::map<std::string, double>& DPars);
    
    virtual bool setFlag(const std::string name, const bool value);
    
    
    /**
     * 
     * @return the coupling ct
     */
    double getc1() const
    {
        return c1;
    };

    /**
     *
     * @return the coupling cg
     */
    double getc2() const
    {
        return c2;
    };
    
    /**
     *
     * @return the coupling cV
     */
    double getc3() const
    {
        return c3;
    };
    
    /**
     *
     * @return the coupling cA
     */
    double getc4() const
    {
        return c4;
    };
    
    /**
     *
     * @return the coupling cA
     */
    bool get_condition_flag() const
    {
        return condition;
    };


protected:
    
    virtual void setParameter(const std::string, const double&);

private:
    
    double c1, c2, c3, c4; /* Model Parameters */
    bool condition;
    
};

#endif	/* MYMODEL_H */


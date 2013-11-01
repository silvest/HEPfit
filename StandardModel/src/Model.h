/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef MODEL_H
#define	MODEL_H

#include <map>

/**
 * @class Model
 * @ingroup StandardModel
 * @brief A class for the template of models. 
 * @author SusyFit Collaboration
 * @copyright GNU General Public License
 * @details 
 */
class Model {
public:
    Model(){
        ModelInitialized = false; 
    };
    virtual bool Update(const std::map<std::string, double>&) = 0;
    virtual bool PreUpdate() = 0;
    virtual bool PostUpdate() = 0;
    virtual bool CheckParameters(const std::map<std::string, double>&) = 0;
    virtual bool SetFlag(const std::string, const bool&) = 0;
    virtual void setParameters(const std::string, const double&) = 0;
    virtual bool Init(const std::map<std::string, double>&) = 0;
    
    virtual std::string ModelName() const = 0;
    
    bool IsModelInitialized() const 
    {
        return ModelInitialized;
    }

    void setModelInitialized(bool ModelInitialized) 
    {
        this->ModelInitialized = ModelInitialized;
    }
    
    
    bool IsUpdateError() const 
    {
        return UpdateError;
    }

    void setUpdateError(bool UpdateError)
    {
        this->UpdateError = UpdateError;
    }
    
    
protected:
    
    bool UpdateError;
    
private:

    bool ModelInitialized;
    
};

#endif	/* MODEL_H */


/* 
 * File:   Model.h
 * Author: marco
 *
 * Created on February 23, 2011, 3:48 PM
 */

#ifndef MODEL_H
#define	MODEL_H

#include <map>

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
    virtual void SetParameter(const std::string, const double&) = 0;
    virtual bool Init(const std::map<std::string, double>&) = 0;
    
    virtual std::string ModelName() const = 0;
    
    bool IsModelInitialized() const {
        return ModelInitialized;
    }

    void SetModelInitialized(bool ModelInitialized) {
        this->ModelInitialized = ModelInitialized;
    }
    
    
    bool IsUpdateError() const {
        return UpdateError;
    }

    void SetUpdateError(bool UpdateError) {
        this->UpdateError = UpdateError;
    }
    
    
protected:
    
    bool UpdateError;
    
private:

    bool ModelInitialized;
    
};

#endif	/* MODEL_H */


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
    
    /**
     * @brief The default constructor.
     */
    Model(){
        ModelInitialized = false; 
    };
    
    /**
     * @brief
     * @details
     * @param[in]
     * @return
     */
    virtual bool Update(const std::map<std::string, double>&) = 0;
    
    /**
     * @brief
     * @details
     * @param[in]
     * @return
     */
    virtual bool PreUpdate() = 0;
    
    /**
     * @brief
     * @details
     * @param[in]
     * @return
     */
    virtual bool PostUpdate() = 0;
    
    /**
     * @brief
     * @details
     * @param[in]
     * @return
     */
    virtual bool CheckParameters(const std::map<std::string, double>&) = 0;
    
    /**
     * @brief
     * @details
     * @param[in]
     * @return
     */
    virtual bool setFlag(const std::string, const bool&) = 0;
    
    /**
     * @brief
     * @details
     * @param[in]
     * @return
     */
    virtual bool CheckFlags() const = 0;
    
    /**
     * @brief
     * @details
     * @param[in]
     * @return
     */
    virtual void setParameter(const std::string, const double&) = 0;
    
    /**
     * @brief
     * @details
     * @param[in]
     * @return
     */
    virtual bool Init(const std::map<std::string, double>&) = 0;
    
    /**
     * @brief
     * @details
     * @param[in]
     * @return
     */
    virtual std::string ModelName() const = 0;
    
    /**
     * @brief
     * @details
     * @param[in]
     * @return
     */
    bool IsModelInitialized() const
    {
        return ModelInitialized;
    }
    
    /**
     * @brief
     * @details
     * @param[in]
     * @return
     */
    void setModelInitialized(bool ModelInitialized)
    {
        this->ModelInitialized = ModelInitialized;
    }
    
    /**
     * @brief
     * @details
     * @param[in]
     * @return
     */
    bool IsUpdateError() const
    {
        return UpdateError;
    }
    
    /**
     * @brief
     * @details
     * @param[in]
     * @return
     */
    void setUpdateError(bool UpdateError)
    {
        this->UpdateError = UpdateError;
    }
    
    
protected:
    
    bool UpdateError; /**< */
    
private:

    bool ModelInitialized; /**< */
    
};

#endif	/* MODEL_H */


/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/cppFiles/class.h to edit this template
 */

/* 
 * File:   TopQuarkObservables.h
 * Author: silvest
 *
 * Created on 21 settembre 2023, 16.45
 */

#ifndef TOPQUARKOBSERVABLES_H
#define TOPQUARKOBSERVABLES_H

#include "NPSMEFTd6General.h"
#include "ThObservable.h"

class TopQuarkObservables {
public:

    
    
    TopQuarkObservables(const NPSMEFTd6General& NP_i);
//    TopQuarkObservables(const TopQuarkObservables& orig);
    virtual ~TopQuarkObservables(){};
    

    const NPSMEFTd6General& GetNP() const {
        return NP;
    }

    
    inline double ewgc(const std::string name) const
    {
        return NP.ewgc(name);
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return NP.ewgc(name, i, j);
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return NP.ewgc(name, i, j, k, l);
    }
    
    
protected:
        


private:
    
    const NPSMEFTd6General& NP;
  
};
    
    
    /**
    * @class F0
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class F0_LO: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    F0_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    inline double ewgc(const std::string name) const
    {
        return mytopobs.ewgc(name);
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.ewgc(name, i, j);
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.ewgc(name, i, j, k, l);
    }
    
    
    private:

    const TopQuarkObservables mytopobs;
    
    

    };


    /**
    * @class FL
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class FL_LO: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    FL_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    inline double ewgc(const std::string name) const
    {
        return mytopobs.ewgc(name);
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.ewgc(name, i, j);
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.ewgc(name, i, j, k, l);
    }
    
    
    private:

    const TopQuarkObservables mytopobs;
    
    
    

    };





    /**
    * @class Test_direct
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class Test_direct: public ThObservable {
    public:

    /**
     * @brief Test_direct constructor.
     */
    Test_direct(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();

    
    inline double ewgc(const std::string name) const
    {
        return mytopobs.ewgc(name);
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.ewgc(name, i, j);
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.ewgc(name, i, j, k, l);
    }
    
    
    private:

    const TopQuarkObservables mytopobs;


    
    };
    
    



#endif /* TOPQUARKOBSERVABLES_H */


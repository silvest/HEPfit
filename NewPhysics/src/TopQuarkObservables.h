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

    
    inline double getSMEFTCoeffEW(const std::string name) const
    {
        return NP.getSMEFTCoeffEW(name);
    }
    
    inline double getSMEFTCoeffEW(const std::string name, int i, int j) const
    {
        return NP.getSMEFTCoeffEW(name, i, j);
    }
    
    inline double getSMEFTCoeffEW(const std::string name, int i, int j, int k, int l) const
    {
        return NP.getSMEFTCoeffEW(name, i, j, k, l);
    }
    
    
protected:
        


private:
    
    const NPSMEFTd6General& NP;
  
};

    



    /**
    * @class FB_asymmetry_Tevatron_tt_diff_mtt_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class FB_asymmetry_Tevatron_tt_diff_mtt_LO: public ThObservable {
    public:

    /**
     * @brief charge_asymmetry_tt_diff_mtt_LO constructor.
     */
    FB_asymmetry_Tevatron_tt_diff_mtt_LO(const StandardModel& SM_i);

    /**
     * @return The value of charge_asymmetry_tt_diff_mtt_LO
     */
    double computeThValue();
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const TopQuarkObservables mytopobs;

    };




    
    /**
    * @class sigma_ttbar_diff
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tt_diff_LO: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    sigma_tt_diff_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    
    
    
    
    /**
    * @class charge_asymmetry_tt_diff_mtt_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class charge_asymmetry_tt_diff_mtt_LO: public ThObservable {
    public:

    /**
     * @brief charge_asymmetry_tt_diff_mtt_LO constructor.
     */
    charge_asymmetry_tt_diff_mtt_LO(const StandardModel& SM_i);

    /**
     * @return The value of charge_asymmetry_tt_diff_mtt_LO
     */
    double computeThValue();
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    
    
    
    
    /**
    * @class sigma_ttz_diff
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_ttz_diff_LO_CMS_190711270: public ThObservable {
    public:

    /**
     * @brief sigma_ttz_diff_LO constructor.
     */
    sigma_ttz_diff_LO_CMS_190711270(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    
    
    
    
    
    
    
    
    /**
    * @class sigma_ttz_diff
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_ttz_diff_LO_ATLAS_210312603: public ThObservable {
    public:

    /**
     * @brief sigma_ttz_diff_LO constructor.
     */
    sigma_ttz_diff_LO_ATLAS_210312603(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    
    
    
    
    
    
    /**
    * @class sigma_tta_diff
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tta_diff_LO_ATLAS_emu_200706946: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    sigma_tta_diff_LO_ATLAS_emu_200706946(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    

    
    
    
    
    
    /**
    * @class sigma_tta_diff
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tta_diff_LO_CMS_semileptonic_210701508: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    sigma_tta_diff_LO_CMS_semileptonic_210701508(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    
    
    
    
    /**
    * @class sigma_tta_diff
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tta_diff_LO_CMS_dilepton_220107301: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    sigma_tta_diff_LO_CMS_dilepton_220107301(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
     
        double b_min = 0.;
        double b_max = 0.;
        
        
        const TopQuarkObservables mytopobs;

    };
    

    
    
    
    
    
    
        
    /**
    * @class s-channel production of top quarks, sigma_tb_13_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tb_8_LO: public ThObservable {
    public:

    /**
     * @brief s-channel production of top quarks at 8 TeV constructor.
     */
    sigma_tb_8_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    
    /**
    * @class s-channel production of top quarks, sigma_tb_13_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tb_13_LO: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    sigma_tb_13_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    
    
    
    
    /**
    * @class t-channel production of top quarks, sigma_tq_7_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tq_7_LO: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    sigma_tq_7_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    
    
    
    
    /**
    * @class t-channel production of top quarks, sigma_tq_13_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tq_8_LO: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    sigma_tq_8_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    
    
    /**
    * @class t-channel production of top quarks, sigma_tq_13_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tq_13_LO: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    sigma_tq_13_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
        
        const TopQuarkObservables mytopobs;

    };
    

    
    
    
    
    
    /**
    * @class sigma_taq_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_taq_LO_CMS: public ThObservable {
    public:

    /**
     * @brief sigma_taq_LO constructor.
     */
    sigma_taq_LO_CMS(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
        
        const TopQuarkObservables mytopobs;

    };
    
    
    /**
    * @class sigma_taq_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_taq_LO_ATLAS: public ThObservable {
    public:

    /**
     * @brief sigma_taq_LO constructor.
     */
    sigma_taq_LO_ATLAS(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    
    
    
    
    /**
    * @class sigma_tzq_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tzq_LO: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    sigma_tzq_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    
    
    
    
    /**
    * @class sigma_tw_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tw_7_LO: public ThObservable {
    public:

    /**
     * @brief sigma_tw_7_LO constructor.
     */
    sigma_tw_7_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    
    
    
    
    
    
    
    /**
    * @class sigma_tw_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tw_8_LO: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    sigma_tw_8_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    
    
    
    
    
    /**
    * @class sigma_tw_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_tw_13_LO: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    sigma_tw_13_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
        
        const TopQuarkObservables mytopobs;

    };
    

    
    
    
    
    /**
    * @class sigma_ttw_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class sigma_ttw_LO: public ThObservable {
    public:

    /**
     * @brief sigma_ttw_LO constructor.
     */
    sigma_ttw_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
        
        const TopQuarkObservables mytopobs;

    };
    
    
    
    

    /**
    * @class R_ttw_LO
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class R_ttw_LO: public ThObservable {
    public:

    /**
     * @brief R_ttw_LO constructor.
     */
    R_ttw_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }

    private:
        
        const TopQuarkObservables mytopobs;

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
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
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
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }
    
    
    private:

    const TopQuarkObservables mytopobs;
    
    
    

    };


    
    
    
    
    
    /**
    * @class FR
    * @ingroup TopQuarkObservables 
    * @brief Test Observable
    */
    class FR_LO: public ThObservable {
    public:

    /**
     * @brief FL constructor.
     */
    FR_LO(const StandardModel& SM_i);

    /**
     * @return The value of Test_direct
     */
    double computeThValue();
    
    
    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double ewgc(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double ewgc(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
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

    //We need to multiply by the square of the scale used to obtain the parametrisation (1000 GeV)
    inline double getSMEFTCoeffEW(const std::string name) const
    {
        return mytopobs.getSMEFTCoeffEW(name)*1000000;
    }
    
    inline double getSMEFTCoeffEW(const std::string name, int i, int j) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j)*1000000;
    }
    
    inline double getSMEFTCoeffEW(const std::string name, int i, int j, int k, int l) const
    {
        return mytopobs.getSMEFTCoeffEW(name, i, j, k, l)*1000000;
    }
    
    
    private:

    const TopQuarkObservables mytopobs;


    
    };
    
    



#endif /* TOPQUARKOBSERVABLES_H */

/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef GMDIRECTSEARCHES_H
#define	GMDIRECTSEARCHES_H

#include "ThObservable.h"

class GeorgiMachacek;


/**
 * @class GMDirectSearches
 * @ingroup GeorgiMachacek
 * @brief Contains several classes for direct searches constraining
 * the Gerogi Machachek model.
 * @author HEPfit Collaboration
 * @copyright GNU General Public License
 */

/**
 * @class BR_H1_tt_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_tt_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_tt_GM constructor.
     */
    BR_H1_tt_GM(const StandardModel& SM_i);
    
    /**
     * @return  @f$BR(H_1 \to tt)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_bb_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_bb_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_bb_GM constructor.
     */
    BR_H1_bb_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to b\bar b)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_tautau_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_tautau_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_tautau_GM constructor.
     */
    BR_H1_tautau_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to \tau \tau)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_WW_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_WW_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_WW_GM constructor.
     */
    BR_H1_WW_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to WW)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_ZZ_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_ZZ_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_ZZ_GM constructor.
     */
    BR_H1_ZZ_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to ZZ)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_gaga_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_gaga_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_gaga_GM constructor.
     */
    BR_H1_gaga_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to \gamma \gamma)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_Zga_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_Zga_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_Zga_GM constructor.
     */
    BR_H1_Zga_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to Z \ganna)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_H3Z_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_H3Z_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_H3Z_GM constructor.
     */
    BR_H1_H3Z_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to AZ)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_H3pW_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_H3pW_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_H3pW_GM constructor.
     */
    BR_H1_H3pW_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to H^{\pm} W^{\mp)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_hh_GM
 * @ingroup GeorgiMachacek
 * @brief Branching ratio of @f$H_1@f$ to two @f$h@f$ in the GeorgiMachacek model.
 */
class BR_H1_hh_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_hh_GM constructor.
     */
    BR_H1_hh_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to hh)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_H3H3_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_H3H3_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_H3H3_GM constructor.
     */
    BR_H1_H3H3_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to AA)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_H3pH3m_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_H3pH3m_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_H3pH3m_GM constructor.
     */
    BR_H1_H3pH3m_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to H^+ H^-)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_H5H5_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_H5H5_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_H5H5_GM constructor.
     */
    BR_H1_H5H5_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to H_5 H_5)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_H5pH5m_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_H5pH5m_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_H5pH5m_GM constructor.
     */
    BR_H1_H5pH5m_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to H_5^+ H_5^-)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H1_H5ppH5mm_GM
 * @ingroup GeorgiMachacek
 */
class BR_H1_H5ppH5mm_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H1_H5ppH5mm_GM constructor.
     */
    BR_H1_H5ppH5mm_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_1 \to H_5^{++} H_5^{--})@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3_tt_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3_tt_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3_tt_GM constructor.
     */
    BR_H3_tt_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_3 \to tt)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3_bb_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3_bb_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3_bb_GM constructor.
     */
    BR_H3_bb_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_3 \to b\bar b)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3_tautau_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3_tautau_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3_tautau_GM constructor.
     */
    BR_H3_tautau_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_3 \to \tau \tau)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3_gaga_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3_gaga_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3_gaga_GM constructor.
     */
    BR_H3_gaga_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_3 \to \gamma \gamma)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3_Zga_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3_Zga_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3_Zga_GM constructor.
     */
    BR_H3_Zga_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_3 \to Z \gamma)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3_hZ_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3_hZ_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3_hZ_GM constructor.
     */
    BR_H3_hZ_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_3 \to hZ)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3_H1Z_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3_H1Z_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3_H1Z_GM constructor.
     */
    BR_H3_H1Z_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_3 \to H_1 Z)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3_H5Z_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3_H5Z_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3_H5Z_GM constructor.
     */
    BR_H3_H5Z_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_3 \to H_5 Z)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3_H5pW_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3_H5pW_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3_H5pW_GM constructor.
     */
    BR_H3_H5pW_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_3 \to H_5^+ W^-)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3p_taunu_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3p_taunu_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3p_taunu_GM constructor.
     */
    BR_H3p_taunu_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_3^+ \to \tau^+ \nu)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3p_tb_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3p_tb_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3p_tb_GM constructor.
     */
    BR_H3p_tb_GM(const StandardModel& SM_i);
    
    /**
     * @return  @f$BR(H_3^+ \to tb)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3p_hW_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3p_hW_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3p_hW_GM constructor.
     */
    BR_H3p_hW_GM(const StandardModel& SM_i);
    
    /**
     * @return  @f$BR(H_3^+ \to h W)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3p_H1W_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3p_H1W_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3p_H1W_GM constructor.
     */
    BR_H3p_H1W_GM(const StandardModel& SM_i);
    
    /**
     * @return  @f$BR(H_3^+ \to H_1 W^+)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3p_H5pZ_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3p_H5pZ_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3p_H5pZ_GM constructor.
     */
    BR_H3p_H5pZ_GM(const StandardModel& SM_i);
    
    /**
     * @return  @f$BR(H_3^+ \to H_5^+ Z)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3p_H5W_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3p_H5W_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3p_H5W_GM constructor.
     */
    BR_H3p_H5W_GM(const StandardModel& SM_i);
    
    /**
     * @return  @f$BR(H_3^+ \to H_5 W^+)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H3p_H5ppW_GM
 * @ingroup GeorgiMachacek
 */
class BR_H3p_H5ppW_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H3p_H5ppW_GM constructor.
     */
    BR_H3p_H5ppW_GM(const StandardModel& SM_i);
    
    /**
     * @return  @f$BR(H_3^+ \to H_5^{++} W^-)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5_WW_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5_WW_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5_WW_GM constructor.
     */
    BR_H5_WW_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5 \to WW)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5_ZZ_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5_ZZ_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5_ZZ_GM constructor.
     */
    BR_H5_ZZ_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5 \to ZZ)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5_gaga_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5_gaga_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5_gaga_GM constructor.
     */
    BR_H5_gaga_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5 \to \gamma \gamma)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5_Zga_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5_Zga_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5_Zga_GM constructor.
     */
    BR_H5_Zga_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5 \to Z \gamma)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5_H3Z_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5_H3Z_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5_H3Z_GM constructor.
     */
    BR_H5_H3Z_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5 \to H_3 Z)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5_H3pW_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5_H3pW_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5_H3pW_GM constructor.
     */
    BR_H5_H3pW_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5 \to H_3^+ W^-)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5_H3pH3m_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5_H3pH3m_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5_H3pH3m_GM constructor.
     */
    BR_H5_H3pH3m_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5 \to H_3^+ H_3^-)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5_H3H3_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5_H3H3_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5_H3H3_GM constructor.
     */
    BR_H5_H3H3_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5 \to H_3 H_3)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5p_WZ_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5p_WZ_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5p_WZ_GM constructor.
     */
    BR_H5p_WZ_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5^+ \to WZ)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5p_H3W_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5p_H3W_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5p_H3W_GM constructor.
     */
    BR_H5p_H3W_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5^+ \to H_3^ W^-)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5p_H3pZ_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5p_H3pZ_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5p_H3pZ_GM constructor.
     */
    BR_H5p_H3pZ_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5^+ \to H_3^+ Z)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5p_H3pH3_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5p_H3pH3_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5p_H3pH3_GM constructor.
     */
    BR_H5p_H3pH3_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5^+ \to H_3^+ H_3)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5pp_WW_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5pp_WW_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5pp_WW_GM constructor.
     */
    BR_H5pp_WW_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5^{++} \to WW)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5pp_H3pW_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5pp_H3pW_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5pp_H3pW_GM constructor.
     */
    BR_H5pp_H3pW_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5^{++} \to H_3^+ W^+)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class BR_H5pp_H3pH3p_GM
 * @ingroup GeorgiMachacek
 */
class BR_H5pp_H3pH3p_GM : public ThObservable {
public:
    
    /**
     * @brief BR_H5pp_H3pH3p_GM constructor.
     */
    BR_H5pp_H3pH3p_GM(const StandardModel& SM_i);
    
    /**
     * @return @f$BR(H_5^{++} \to H_3^+ H_3^+)@f$ in the GeorgiMachacek model.
     */
    double computeThValue ();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_tt_H1_tt_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$tt\to H_1\to tt@f$.
 */

class Hobs_tt_H1_tt_ATLAS13: public ThObservable {
public:

     /**
     * @brief Hobs_tt_H1_tt_ATLAS13 constructor.
     */
    Hobs_tt_H1_tt_ATLAS13(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{tt\to H_1}\cdot BR^{\text{GM}}(H_1\to tt)]_{\text{theo}} / [\sigma_{tt\to H_1}\cdot BR(H_1\to tt)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H1_tt_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_1\to tt@f$.
 */
class Hobs_bb_H1_tt_ATLAS13: public ThObservable {
public:

     /**
     * @brief Hobs_bb_H1_tt_ATLAS13 constructor.
     */
    Hobs_bb_H1_tt_ATLAS13(const StandardModel& SM_i);
     /**
     * @return @f$[\sigma^{\text{GM}}_{bbt\to H_1}\cdot BR^{\text{GM}}(H_1\to tt)]_{\text{theo}} / [\sigma_{b\bar b\to H_1}\cdot BR(H_1\to tt)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_tt_H3_tt_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$tt\to H_3\to tt@f$.
 */
class Hobs_tt_H3_tt_ATLAS13: public ThObservable {
public:

     /**
     * @brief Hobs_tt_H3_tt_ATLAS13 constructor.
     */
    Hobs_tt_H3_tt_ATLAS13(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{tt\to H_3}\cdot BR^{\text{GM}}(H_3\to tt)]_{\text{theo}} / [\sigma_{tt\to H_3}\cdot BR(H_3\to tt)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H3_tt_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_3\to tt@f$.
 */
class Hobs_bb_H3_tt_ATLAS13: public ThObservable {
public:

     /**
     * @brief Hobs_bb_H3_tt_ATLAS13 constructor.
     */
    Hobs_bb_H3_tt_ATLAS13(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_3}\cdot BR^{\text{GM}}(H_3\to tt)]_{\text{theo}} / [\sigma_{b\bar b\to H_3}\cdot BR(H_3\to tt)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H1_bb_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_1\to b\bar b@f$.
 */
class Hobs_bb_H1_bb_CMS8: public ThObservable {
public:

     /**
     * @brief Hobs_bb_H1_bb_CMS8 constructor.
     */
    Hobs_bb_H1_bb_CMS8(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to b\bar b)]_{\text{theo}} / [\sigma_{b\bar b\to H_1}\cdot BR(H_1\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_bb_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_1\to b\bar b@f$.
 */
class Hobs_gg_H1_bb_CMS8: public ThObservable {
public:

     /**
     * @brief Hobs_gg_H1_bb_CMS8 constructor.
     */
    Hobs_gg_H1_bb_CMS8(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to b\bar b)]_{\text{theo}} / [\sigma_{gg\to H_1}\cdot BR(H_1\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_bb_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_1\to b\bar b@f$.
 */
class Hobs_pp_H1_bb_CMS13: public ThObservable {
public:

     /**
     * @brief Hobs_pp_H1_bb_CMS13 constructor.
     */
    Hobs_pp_H1_bb_CMS13(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to b\bar b)]_{\text{theo}} / [\sigma_{pp\to H_1}\cdot BR(H_1\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H1_bb_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_1\to b\bar b@f$.
 */
class Hobs_bb_H1_bb_CMS13: public ThObservable {
public:

     /**
     * @brief Hobs_bb_H1_bb_CMS13 constructor.
     */
    Hobs_bb_H1_bb_CMS13(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to b\bar b)]_{\text{theo}} / [\sigma_{b\bar b\to H_1}\cdot BR(H_1\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H3_bb_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_3\to b\bar b@f$.
 */
class Hobs_bb_H3_bb_CMS8: public ThObservable {
public:

     /**
     * @brief Hobs_bb_H3_bb_CMS8 constructor.
     */
    Hobs_bb_H3_bb_CMS8(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_3}\cdot BR^{\text{GM}}(H_3\to b\bar b)]_{\text{theo}} / [\sigma_{b\bar b\to H_3}\cdot BR(H_3\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_bb_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_3\to b\bar b@f$.
 */
class Hobs_gg_H3_bb_CMS8: public ThObservable {
public:

     /**
     * @brief Hobs_gg_H3_bb_CMS8 constructor.
     */
    Hobs_gg_H3_bb_CMS8(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to b\bar b)]_{\text{theo}} / [\sigma_{gg\to H_3}\cdot BR(H_3\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H3_bb_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_3\to b\bar b@f$.
 */
class Hobs_pp_H3_bb_CMS13: public ThObservable {
public:

     /**
     * @brief Hobs_pp_H3_bb_CMS13 constructor.
     */
    Hobs_pp_H3_bb_CMS13(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_3}\cdot BR^{\text{GM}}(H_3\to b\bar b)]_{\text{theo}} / [\sigma_{pp\to H_3}\cdot BR(H_3\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H3_bb_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_3\to b\bar b@f$.
 */
class Hobs_bb_H3_bb_CMS13: public ThObservable {
public:

     /**
     * @brief Hobs_bb_H3_bb_CMS13 constructor.
     */
    Hobs_bb_H3_bb_CMS13(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_3}\cdot BR^{\text{GM}}(H_3\to b\bar b)]_{\text{theo}} / [\sigma_{b\bar b\to H_3}\cdot BR(H_3\to b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_tautau_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_1\to \tau \tau@f$.
 */
class Hobs_gg_H1_tautau_CMS8: public ThObservable {
public:

     /**
     * @brief Hobs_gg_H1_tautau_CMS8 constructor.
     */
    Hobs_gg_H1_tautau_CMS8(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau \tau)]_{\text{theo}} / [\sigma_{gg\to H_1}\cdot BR(H_1\to \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H1_tautau_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_1\to \tau \tau@f$.
 */
class Hobs_bb_H1_tautau_CMS8: public ThObservable {
public:

     /**
     * @brief Hobs_bb_H1_tautau_CMS8 constructor.
     */
    Hobs_bb_H1_tautau_CMS8(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau \tau)]_{\text{theo}} / [\sigma_{b\bar b\to H_1}\cdot BR(H_1\to \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_tautau_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_1\to \tau \tau@f$.
 */
class Hobs_gg_H1_tautau_ATLAS13: public ThObservable {
public:

     /**
     * @brief Hobs_gg_H1_tautau_ATLAS13 constructor.
     */
    Hobs_gg_H1_tautau_ATLAS13(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau \tau)]_{\text{theo}} / [\sigma_{gg\to H_1}\cdot BR(H_1\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_tautau_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_1\to \tau \tau@f$.
 */
class Hobs_gg_H1_tautau_CMS13: public ThObservable {
public:

     /**
     * @brief Hobs_gg_H1_tautau_CMS13 constructor.
     */
    Hobs_gg_H1_tautau_CMS13(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau \tau)]_{\text{theo}} / [\sigma_{gg\to H_1}\cdot BR(H_1\to \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H1_tautau_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_1\to \tau \tau@f$.
 */
class Hobs_bb_H1_tautau_ATLAS13: public ThObservable {
public:

     /**
     * @brief Hobs_bb_H1_tautau_ATLAS13 constructor.
     */
    Hobs_bb_H1_tautau_ATLAS13(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau \tau)]_{\text{theo}} / [\sigma_{b\bar b\to H_1}\cdot BR(H_1\to \tau \tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H1_tautau_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_1\to \tau \tau@f$.
 */
class Hobs_bb_H1_tautau_CMS13: public ThObservable {
public:

     /**
     * @brief Hobs_bb_H1_tautau_CMS13 constructor.
     */
    Hobs_bb_H1_tautau_CMS13(const StandardModel& SM_i);
    /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau \tau)]_{\text{theo}} / [\sigma_{b\bar b\to H_1}\cdot BR(H_1\to \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_tautau_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section
 *  times branching ratio of the process @f$gg\to H_1\to \tau\tau@f$.
 */
class Hobs_gg_H1_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_tautau_ATLAS8 constructor.
     */
    Hobs_gg_H1_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]_{\text{theo}} / [\sigma_{gg\to H_1}\cdot BR(H_1\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H1_tautau_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_1\to \tau\tau@f$.
 */
class Hobs_bb_H1_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_bb_H1_tautau_ATLAS8 constructor.
     */
    Hobs_bb_H1_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]_{\text{theo}} / [\sigma_{b\bar b\to H_1}\cdot BR(H_1\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_tautau_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_3\to \tau\tau@f$.
 */
class Hobs_gg_H3_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_tautau_ATLAS8 constructor.
     */
    Hobs_gg_H3_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_3}\cdot BR^{\text{GM}}(H_3\to \tau\tau)]_{\text{theo}} / [\sigma_{gg \to H_3}\cdot BR(H_3\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_tautau_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_3\to \tau\tau@f$.
 */
class Hobs_gg_H3_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_tautau_CMS8 constructor.
     */
    Hobs_gg_H3_tautau_CMS8(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_3}\cdot BR^{\text{GM}}(H_3\to \tau\tau)]_{\text{theo}} / [\sigma_{gg \to H_3}\cdot BR(H_3\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H3_tautau_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$bb \to H_3\to \tau\tau@f$.
 */
class Hobs_bb_H3_tautau_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_bb_H3_tautau_ATLAS8 constructor.
     */
    Hobs_bb_H3_tautau_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{bb \to H_3}\cdot BR^{\text{GM}}(H_3\to \tau\tau)]_{\text{theo}} / [\sigma_{bb \to H_3}\cdot BR(H_3\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};


/**
 * @class Hobs_bb_H3_tautau_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$bb \to H_3\to \tau\tau@f$.
 */
class Hobs_bb_H3_tautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_bb_H3_tautau_CMS8 constructor.
     */
    Hobs_bb_H3_tautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{bb \to H_3}\cdot BR^{\text{GM}}(H_3\to \tau\tau)]_{\text{theo}} / [\sigma_{bb \to H_3}\cdot BR(H_3\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_tautau_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_3\to \tau\tau@f$.
 */
class Hobs_gg_H3_tautau_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_tautau_ATLAS13 constructor.
     */
    Hobs_gg_H3_tautau_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_3}\cdot BR^{\text{GM}}(H_3\to \tau\tau)]_{\text{theo}} / [\sigma_{gg \to H_3}\cdot BR(H_3\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_tautau_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_3\to \tau\tau@f$.
 */
class Hobs_gg_H3_tautau_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_tautau_CMS13 constructor.
     */
    Hobs_gg_H3_tautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_3}\cdot BR^{\text{GM}}(H_3\to \tau\tau)]_{\text{theo}} / [\sigma_{gg \to H_3}\cdot BR(H_3\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H3_tautau_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$bb \to H_3\to \tau\tau@f$.
 */
class Hobs_bb_H3_tautau_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_H3_tautau_ATLAS13 constructor.
     */
    Hobs_bb_H3_tautau_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{bb \to H_3}\cdot BR^{\text{GM}}(H_3\to \tau\tau)]_{\text{theo}} / [\sigma_{bb \to H_3}\cdot BR(H_3\to \tau\tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H3_tautau_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$bb \to H_3\to \tau\tau@f$.
 */
class Hobs_bb_H3_tautau_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_H3_tautau_CMS13 constructor.
     */
    Hobs_bb_H3_tautau_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{bb \to H_3}\cdot BR^{\text{GM}}(H_3\to \tau\tau)]_{\text{theo}} / [\sigma_{bb \to H_3}\cdot BR(H_3\to \tau\tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_gaga_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$bb \to H_1 \to \gamma \gamma@f$.
 */
class Hobs_gg_H1_gaga_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_gaga_ATLAS8 constructor.
     */
    Hobs_gg_H1_gaga_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_1}\cdot BR^{\text{GM}}(H_1\to \gamma \gamma)]_{\text{theo}} / [\sigma_{gg \to H_1}\cdot BR(H_1\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_gaga_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_1 \to \gamma \gamma@f$.
 */
class Hobs_pp_H1_gaga_ATLAS13: public ThObservable {
public:

     /**
     * @brief Hobs_pp_H1_gaga_ATLAS13 constructor.
     */
    Hobs_pp_H1_gaga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_1}\cdot BR^{\text{GM}}(H_1\to \gamma \gamma)]_{\text{theo}} / [\sigma_{pp \to H_1}\cdot BR(H_1\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_gaga_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to \gamma \gamma@f$.
 */
class Hobs_gg_H1_gaga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_gaga_CMS13 constructor.
     */
    Hobs_gg_H1_gaga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_1}\cdot BR^{\text{GM}}(H_1\to \gamma \gamma)]_{\text{theo}} / [\sigma_{gg \to H_1}\cdot BR(H_1\to \gamma \gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_gaga_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_3 \to \gamma \gamma@f$.
 */
class Hobs_gg_H3_gaga_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_gaga_ATLAS8 constructor.
     */
    Hobs_gg_H3_gaga_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_3}\cdot BR^{\text{GM}}(H_3\to \gamma \gamma)]_{\text{theo}} / [\sigma_{gg \to H_3}\cdot BR(H_3\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H3_gaga_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_3 \to \gamma \gamma@f$.
 */
class Hobs_pp_H3_gaga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H3_gaga_ATLAS13 constructor.
     */
    Hobs_pp_H3_gaga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_3}\cdot BR^{\text{GM}}(H_3\to \gamma \gamma)]_{\text{theo}} / [\sigma_{pp \to H_3}\cdot BR(H_3\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_gaga_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_3 \to \gamma \gamma@f$.
 */
class Hobs_gg_H3_gaga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_gaga_CMS13 constructor.
     */
    Hobs_gg_H3_gaga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_3}\cdot BR^{\text{GM}}(H_3\to \gamma \gamma)]_{\text{theo}} / [\sigma_{gg \to H_3}\cdot BR(H_3\to \gamma \gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5_gaga_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_5 \to \gamma \gamma@f$.
 */
class Hobs_pp_H5_gaga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H5_gaga_ATLAS13 constructor.
     */
    Hobs_pp_H5_gaga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_5}\cdot BR^{\text{GM}}(H_5\to \gamma \gamma)]_{\text{theo}} / [\sigma_{pp \to H_5}\cdot BR(H_5\to \gamma \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_Zga_llga_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_1 \to Z \gamma \to \ell \ell \gamma @f$.
 */
class Hobs_pp_H1_Zga_llga_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_Zga_llga_ATLAS8 constructor.
     */
    Hobs_pp_H1_Zga_llga_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_1}\cdot BR^{\text{GM}}(H_1\to Z\gamma \to \ell \ell \gamma)]_{\text{theo}} / [\sigma_{pp \to H_1}\cdot BR(H_1\to Z\gamma \to \ell \ell \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_Zga_llga_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_1 \to Z \gamma \to \ell \ell \gamma@f$.
 */
class Hobs_pp_H1_Zga_llga_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_Zga_llga_CMS8 constructor.
     */
    Hobs_pp_H1_Zga_llga_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_1}\cdot BR^{\text{GM}}(H_1\to Z\gamma \to \ell \ell \gamma)]_{\text{theo}} / [\sigma_{pp \to H_1}\cdot BR(H_1\to Z\gamma \to \ell \ell \gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_Zga_llga_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to Z \gamma@f$.
 */
class Hobs_gg_H1_Zga_llga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_Zga_llga_ATLAS13 constructor.
     */
    Hobs_gg_H1_Zga_llga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_1}\cdot BR^{\text{GM}}(H_1\to Z\gamma \to \ell \ell \gamma)]_{\text{theo}} / [\sigma_{gg \to H_1}\cdot BR(H_1\to Z\gamma \to \ell \ell \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_Zga_qqga_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to Z \gamma \to q  q\bar \gamma@f$.
 */
class Hobs_gg_H1_Zga_qqga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_Zga_qqga_ATLAS13 constructor.
     */
    Hobs_gg_H1_Zga_qqga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_1}\cdot BR^{\text{GM}}(H_1\to Z\gamma \to q  q\bar \gamma)]_{\text{theo}} / [\sigma_{gg \to H_1}\cdot BR(H_1\to Z\gamma \to q  q\bar \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_Zga_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to Z \gamma@f$.
 */
class Hobs_gg_H1_Zga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_Zga_CMS13 constructor.
     */
    Hobs_gg_H1_Zga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$\left[\sigma^{\text{GM}}_{gg \to H_1}\cdot BR^{\text{GM}}(H_1\to Z\gamma)\right]_{\text{theo}} / \left[\sigma_{gg \to H_1}\cdot BR(H_1\to Z\gamma)\right]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H3_Zga_llga_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_3 \to Z \gamma  \to \ell \ell \gamma@f$.
 */
class Hobs_pp_H3_Zga_llga_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H3_Zga_llga_ATLAS8 constructor.
     */
    Hobs_pp_H3_Zga_llga_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$\left[\sigma^{\text{GM}}_{pp \to H_3}\cdot BR^{\text{GM}}(H_3\to Z\gamma \to \ell \ell \gamma)\right]_{\text{theo}} / \left[\sigma_{pp \to H_3}\cdot BR(H_3\to Z\gamma \to \ell \ell \gamma)\right]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H3_Zga_llga_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_3 \to Z \gamma \to \ell \ell \gamma@f$.
 */
class Hobs_pp_H3_Zga_llga_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H3_Zga_llga_CMS8 constructor.
     */
    Hobs_pp_H3_Zga_llga_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_3}\cdot BR^{\text{GM}}(H_3\to Z\gamma \to \ell \ell \gamma)]_{\text{theo}} / [\sigma_{pp \to H_3}\cdot BR(H_3\to Z\gamma \to \ell \ell \gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_Zga_llga_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_3 \to Z \gamma \to \ell \ell \gamma@f$.
 */
class Hobs_gg_H3_Zga_llga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_Zga_llga_ATLAS13 constructor.
     */
    Hobs_gg_H3_Zga_llga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_3}\cdot BR^{\text{GM}}(H_3\to Z\gamma \to \ell \ell \gamma)]_{\text{theo}} / [\sigma_{gg \to H_3}\cdot BR(H_3\to Z\gamma \to \ell \ell \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_Zga_qqga_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_3 \to Z \gamma \to q q\bar \gamma@f$.
 */
class Hobs_gg_H3_Zga_qqga_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_Zga_qqga_ATLAS13 constructor.
     */
    Hobs_gg_H3_Zga_qqga_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_3}\cdot BR^{\text{GM}}(H_3\to Z\gamma \to \ell \ell \gamma \to q q\bar \gamma)]_{\text{theo}} / [\sigma_{gg \to H_3}\cdot BR(H_3\to Z\gamma \to \ell \ell \gamma \to q q\bar \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_Zga_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_3 \to Z \gamma@f$.
 */
class Hobs_gg_H3_Zga_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_Zga_CMS13 constructor.
     */
    Hobs_gg_H3_Zga_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_3}\cdot BR^{\text{GM}}(H_3\to Z\gamma)]_{\text{theo}} / [\sigma_{gg \to H_3}\cdot BR(H_3\to Z\gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5_Zga_llga_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_5 \to Z \gamma \to \ell \ell \gamma@f$.
 */
class Hobs_pp_H5_Zga_llga_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H5_Zga_llga_ATLAS8 constructor.
     */
    Hobs_pp_H5_Zga_llga_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_5}\cdot BR^{\text{GM}}(H_5\to Z\gamma \to \ell \ell \gamma)]_{\text{theo}} / [\sigma_{pp \to H_5}\cdot BR(H_5\to Z\gamma \to \ell \ell \gamma)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5_Zga_llga_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_5 \to Z \gamma \to \ell \ell \gamma@f$.
 */
class Hobs_pp_H5_Zga_llga_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H5_Zga_llga_CMS8 constructor.
     */
    Hobs_pp_H5_Zga_llga_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_5}\cdot BR^{\text{GM}}(H_5\to Z\gamma \to \ell \ell \gamma)]_{\text{theo}} / [\sigma_{pp \to H_5}\cdot BR(H_5\to Z\gamma \to \ell \ell \gamma)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_ZZ_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to ZZ@f$.
 */
class Hobs_gg_H1_ZZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_ZZ_ATLAS8 constructor.
     */
    Hobs_gg_H1_ZZ_ATLAS8(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_1}\cdot BR^{\text{GM}}(H_1\to ZZ)]_{\text{theo}} / [\sigma_{gg \to H_1}\cdot BR(H_1\to ZZ )]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H1_ZZ_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$VV \to H_1 \to ZZ@f$.
 */
class Hobs_VV_H1_ZZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H1_ZZ_ATLAS8 constructor.
     */
    Hobs_VV_H1_ZZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV \to H_1}\cdot BR^{\text{GM}}(H_1\to ZZ)]_{\text{theo}} / [\sigma_{VV \to H_1}\cdot BR(H_1\to ZZ )]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_ZZ_llllnunu_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to ZZ \to \ell \ell \ell \ell \nu \nu@f$.
 */
class Hobs_gg_H1_ZZ_llllnunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_ZZ_llllnunu_ATLAS13 constructor.
     */
    Hobs_gg_H1_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_1}\cdot BR^{\text{GM}}(H_1\to ZZ)]_{\text{theo}} / [\sigma_{gg \to H_1}\cdot BR(H_1\to ZZ )]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H1_ZZ_llllnunu_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$VV \to H_1 \to ZZ@f$.
 */
class Hobs_VV_H1_ZZ_llllnunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H1_ZZ_llllnunu_ATLAS13 constructor.
     */
    Hobs_VV_H1_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV \to H_1}\cdot BR^{\text{GM}}(H_1\to ZZ)]_{\text{theo}} / [\sigma_{VV \to H_1}\cdot BR(H_1\to ZZ )]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_ZZ_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to ZZ@f$.
 */
class Hobs_gg_H1_ZZ_qqllnunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_ZZ_qqllnunu_ATLAS13 constructor.
     */
    Hobs_gg_H1_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_1}\cdot BR^{\text{GM}}(H_1\to ZZ)]_{\text{theo}} / [\sigma_{gg \to H_1}\cdot BR(H_1\to ZZ )]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H1_ZZ_qqllnunu_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to ZZ \to q  q\bar \ell \ell \nu \nu@f$.
 */
class Hobs_VV_H1_ZZ_qqllnunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H1_ZZ_qqllnunu_ATLAS13 constructor.
     */
    Hobs_VV_H1_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV \to H_1}\cdot BR^{\text{GM}}(H_1\to ZZ \to q  q\bar \ell \ell \nu \n)]_{\text{theo}} / [\sigma_{VV \to H_1}\cdot BR(H_1\to ZZ \to q  q\bar \ell \ell \nu \n)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_ZZ_llqqnunull_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to ZZ@f$.
 */
class Hobs_pp_H1_ZZ_llqqnunull_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_ZZ_llqqnunull_CMS13 constructor.
     */
    Hobs_pp_H1_ZZ_llqqnunull_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_1}\cdot BR^{\text{GM}}(H_1\to ZZ)]_{\text{theo}} / [\sigma_{pp \to H_1}\cdot BR(H_1\to ZZ )]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_ZZ_qqnunu_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_1 \to ZZ \to q  q\bar \nu \nu@f$.
 */
class Hobs_pp_H1_ZZ_qqnunu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_ZZ_qqnunu_CMS13 constructor.
     */
    Hobs_pp_H1_ZZ_qqnunu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_1}\cdot BR^{\text{GM}}(H_1\to ZZ \to q  q\bar \nu \nu)]_{\text{theo}} / [\sigma_{pp \to H_1}\cdot BR(H_1\to ZZ \to q  q\bar \nu \nu )]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H5_ZZ_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$VV \to H_1 \to ZZ@f$.
 */
class Hobs_VV_H5_ZZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H5_ZZ_ATLAS8 constructor.
     */
    Hobs_VV_H5_ZZ_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV \to H_5}\cdot BR^{\text{GM}}(H_5\to ZZ)]_{\text{theo}} / [\sigma_{VV \to H_5}\cdot BR(H_5\to ZZ )]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H5_ZZ_llllnunu_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$VV \to H_5 \to ZZ@f$.
 */
class Hobs_VV_H5_ZZ_llllnunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H5_ZZ_llllnunu_ATLAS13 constructor.
     */
    Hobs_VV_H5_ZZ_llllnunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV \to H_5}\cdot BR^{\text{GM}}(H_5\to ZZ)]_{\text{theo}} / [\sigma_{VV \to H_5}\cdot BR(H_5\to ZZ )]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H5_ZZ_qqllnunu_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$VV \to H_5 \to ZZ@f$.
 */
class Hobs_VV_H5_ZZ_qqllnunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H5_ZZ_qqllnunu_ATLAS13 constructor.
     */
    Hobs_VV_H5_ZZ_qqllnunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV \to H_5}\cdot BR^{\text{GM}}(H_5\to ZZ)]_{\text{theo}} / [\sigma_{VV \to H_5}\cdot BR(H_5\to ZZ )]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5_ZZ_llqqnunull_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_5 \to ZZ@f$.
 */
class Hobs_pp_H5_ZZ_llqqnunull_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H5_ZZ_llqqnunull_CMS13 constructor.
     */
    Hobs_pp_H5_ZZ_llqqnunull_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_5}\cdot BR^{\text{GM}}(H_5\to ZZ)]_{\text{theo}} / [\sigma_{pp \to H_5}\cdot BR(H_5\to ZZ )]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5_ZZ_qqnunu_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_5 \to ZZ \to q  q\bar \nu \nu@f$.
 */
class Hobs_pp_H5_ZZ_qqnunu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H5_ZZ_qqnunu_CMS13 constructor.
     */
    Hobs_pp_H5_ZZ_qqnunu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_5}\cdot BR^{\text{GM}}(H_5\to ZZ \to q  q\bar \nu \nu)]_{\text{theo}} / [\sigma_{pp \to H_5}\cdot BR(H_5\to ZZ \to q  q\bar \nu \nu)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_WW_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to WW@f$.
 */
class Hobs_gg_H1_WW_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_WW_ATLAS8 constructor.
     */
    Hobs_gg_H1_WW_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_1}\cdot BR^{\text{GM}}(H_1\to WW )]_{\text{theo}} / [\sigma_{gg \to H_1}\cdot BR(H_1\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H1_WW_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$VV \to H_1 \to WW@f$.
 */
class Hobs_VV_H1_WW_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H1_WW_ATLAS8 constructor.
     */
    Hobs_VV_H1_WW_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV \to H_1}\cdot BR^{\text{GM}}(H_1\to WW )]_{\text{theo}} / [\sigma_{VV \to H_1}\cdot BR(H_1\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_WW_enumunu_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to WW@f$.
 */
class Hobs_gg_H1_WW_enumunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_WW_enumunu_ATLAS13 constructor.
     */
    Hobs_gg_H1_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_1}\cdot BR^{\text{GM}}(H_1\to WW )]_{\text{theo}} / [\sigma_{gg \to H_1}\cdot BR(H_1\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H1_WW_enumunu_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$VV \to H_1 \to WW@f$.
 */
class Hobs_VV_H1_WW_enumunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H1_WW_enumunu_ATLAS13 constructor.
     */
    Hobs_VV_H1_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV \to H_1}\cdot BR^{\text{GM}}(H_1\to WW )]_{\text{theo}} / [\sigma_{VV \to H_1}\cdot BR(H_1\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_WW_lnuqq_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to WW@f$.
 */
class Hobs_gg_H1_WW_lnuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_WW_lnuqq_ATLAS13 constructor.
     */
    Hobs_gg_H1_WW_lnuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_1}\cdot BR^{\text{GM}}(H_1\to WW )]_{\text{theo}} / [\sigma_{gg \to H_1}\cdot BR(H_1\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H1_WW_lnuqq_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$VV \to H_1 \to WW@f$.
 */
class Hobs_VV_H1_WW_lnuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H1_WW_lnuqq_ATLAS13 constructor.
     */
    Hobs_VV_H1_WW_lnuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV \to H_1}\cdot BR^{\text{GM}}(H_1\to WW )]_{\text{theo}} / [\sigma_{VV \to H_1}\cdot BR(H_1\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_ggVV_H1_WW_lnulnu_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$(gg+VV) \to H_1 \to WW@f$.
 */
class Hobs_ggVV_H1_WW_lnulnu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggVV_H1_WW_lnulnu_CMS13 constructor.
     */
    Hobs_ggVV_H1_WW_lnulnu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{(gg +VV) \to H_1}\cdot BR^{\text{GM}}(H_1\to WW )]_{\text{theo}} / [\sigma_{(gg+VV) \to H_1}\cdot BR(H_1\to WW)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_WW_lnuqq_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_1 \to WW@f$.
 */
class Hobs_pp_H1_WW_lnuqq_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_WW_lnuqq_CMS13 constructor.
     */
    Hobs_pp_H1_WW_lnuqq_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_1}\cdot BR^{\text{GM}}(H_1\to WW )]_{\text{theo}} / [\sigma_{pp \to H_1}\cdot BR(H_1\to WW)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H5_WW_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$VV \to H_5 \to WW@f$.
 */
class Hobs_VV_H5_WW_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H5_WW_ATLAS8 constructor.
     */
    Hobs_VV_H5_WW_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV \to H_5}\cdot BR^{\text{GM}}(H_5\to WW )]_{\text{theo}} / [\sigma_{VV \to H_5}\cdot BR(H_5\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H5_WW_enumunu_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$VV \to H_5 \to WW@f$.
 */
class Hobs_VV_H5_WW_enumunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H5_WW_enumunu_ATLAS13 constructor.
     */
    Hobs_VV_H5_WW_enumunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV \to H_5}\cdot BR^{\text{GM}}(H_5\to WW )]_{\text{theo}} / [\sigma_{VV \to H_5}\cdot BR(H_5\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H5_WW_lnuqq_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$VV \to H_5 \to WW@f$.
 */
class Hobs_VV_H5_WW_lnuqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H5_WW_lnuqq_ATLAS13 constructor.
     */
    Hobs_VV_H5_WW_lnuqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV \to H_5}\cdot BR^{\text{GM}}(H_5\to WW )]_{\text{theo}} / [\sigma_{VV \to H_5}\cdot BR(H_5\to WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_ggVV_H5_WW_lnulnu_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$(gg+VV) \to H_5 \to WW@f$.
 */
class Hobs_ggVV_H5_WW_lnulnu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_ggVV_H5_WW_lnulnu_CMS13 constructor.
     */
    Hobs_ggVV_H5_WW_lnulnu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{(gg+VV) \to H_5}\cdot BR^{\text{GM}}(H_5\to WW )]_{\text{theo}} / [\sigma_{(gg+VV) \to H_5}\cdot BR(H_5\to WW)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5_WW_lnuqq_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_5 \to WW@f$.
 */
class Hobs_pp_H5_WW_lnuqq_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H5_WW_lnuqq_CMS13 constructor.
     */
    Hobs_pp_H5_WW_lnuqq_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_5}\cdot BR^{\text{GM}}(H_5\to WW )]_{\text{theo}} / [\sigma_{pp \to H_5}\cdot BR(H_5\to WW)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_mu_pp_H1_VV_CMS8
 * @ingroup GeorgiMachacek
  * @brief Ratio of the prediction and CMS upper limit for the signal strength 
 * of the process @f$pp \to H_1\to VV@f$.
 */
class Hobs_mu_pp_H1_VV_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_mu_pp_H1_VV_CMS8 constructor.
     */
    Hobs_mu_pp_H1_VV_CMS8(const StandardModel& SM_i);

     /**
     * @return @f$[\mu_H^{\text{GM}}(H_1\to VV)]_{\text{theo}} / [\mu_H_1(H\to VV)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_VV_qqqq_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_1 \to VV@f$.
 */
class Hobs_pp_H1_VV_qqqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_VV_qqqq_ATLAS13 constructor.
     */
    Hobs_pp_H1_VV_qqqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_1}\cdot BR^{\text{GM}}(H_1\to VV )]_{\text{theo}} / [\sigma_{pp \to H_1}\cdot BR(H_1\to VV)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_mu_pp_H5_VV_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the signal strength 
 * of the process @f$pp \to H_5\to VV@f$.
 */
class Hobs_mu_pp_H5_VV_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_mu_pp_H5_VV_CMS8 constructor.
     */
    Hobs_mu_pp_H5_VV_CMS8(const StandardModel& SM_i);

     /**
     * @return @f$[\mu_H^{\text{THDM}}(H_5\to VV)]_{\text{GM}} / [\mu_H(H_5\to VV)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5_VV_qqqq_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_5 \to WW@f$.
 */
class Hobs_pp_H5_VV_qqqq_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H5_VV_qqqq_ATLAS13 constructor.
     */
    Hobs_pp_H5_VV_qqqq_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_5}\cdot BR^{\text{GM}}(H_5\to VV )]_{\text{theo}} / [\sigma_{pp \to H_5}\cdot BR(H_5\to VV)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_hh_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to hh@f$.
 */
class Hobs_gg_H1_hh_ATLAS8: public ThObservable {
public:

     /**
     * @brief Hobs_gg_H1_hh_ATLAS8 constructor.
     */
    Hobs_gg_H1_hh_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_1}\cdot BR^{\text{GM}}(H_1\to hh )]_{\text{theo}} / [\sigma_{gg \to H_1}\cdot BR(H_1\to hh)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_hh_bbbb_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_1 \to hh \to b\bar b b\bar b@f$.
 */
class Hobs_pp_H1_hh_bbbb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_hh_bbbb_CMS8 constructor.
     */
    Hobs_pp_H1_hh_bbbb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_1}\cdot BR^{\text{GM}}(H_1\to hh \to b\bar b b\bar b )]_{\text{theo}} / [\sigma_{pp \to H_1}\cdot BR(H_1\to hh \to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_hh_gagabb_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_1 \to hh \to \gamma\gamma b\bar b@f$.
 */
class Hobs_pp_H1_hh_gagabb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_hh_gagabb_CMS8 constructor.
     */
    Hobs_pp_H1_hh_gagabb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_1}\cdot BR^{\text{GM}}(H_1\to hh  \to \gamma \gamma b\bar b)]_{\text{theo}} / [\sigma_{pp \to H_1}\cdot BR(H_1\to hh \to \gamma \gamma b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_hh_bbtautau_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg \to H_1 \to hh \to b\bar b \tau \tau@f$.
 */
class Hobs_gg_H1_hh_bbtautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_hh_bbtautau_CMS8 constructor.
     */
    Hobs_gg_H1_hh_bbtautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg \to H_1}\cdot BR^{\text{GM}}(H_1\to hh \to b\bar b \tau \tau )]_{\text{theo}} / [\sigma_{gg \to H_1}\cdot BR(H_1\to hh \to b\bar b \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_hh_bbtautau_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_1 \to hh \to b\bar b \tau \tau@f$.
 */
class Hobs_pp_H1_hh_bbtautau_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_hh_bbtautau_CMS8 constructor.
     */
    Hobs_pp_H1_hh_bbtautau_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_1}\cdot BR^{\text{GM}}(H_1\to hh \to b\bar b \tau \tau)]_{\text{theo}} / [\sigma_{pp \to H_1}\cdot BR(H_1\to hh \to b\bar b \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_hh_bbbb_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp \to H_1 \to hh \to  b\bar b b\bar b @f$.
 */
class Hobs_pp_H1_hh_bbbb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_hh_bbbb_ATLAS13 constructor.
     */
    Hobs_pp_H1_hh_bbbb_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp \to H_1}\cdot BR^{\text{GM}}(H_1\to hh \to  b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp \to H_1}\cdot BR(H_1\to hh \to  b\bar b b\bar b)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_hh_bbbb_1_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_1\to hh\to b\bar b b\bar b@f$.
 */
class Hobs_pp_H1_hh_bbbb_1_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_hh_bbbb_1_CMS13 constructor.
     */
    Hobs_pp_H1_hh_bbbb_1_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{THDM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to H_1}\cdot BR(H_1\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_hh_bbbb_2_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_1\to hh\to b\bar b b\bar b@f$.
 */
class Hobs_pp_H1_hh_bbbb_2_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_hh_bbbb_2_CMS13 constructor.
     */
    Hobs_pp_H1_hh_bbbb_2_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{pp\to H_1}\cdot BR(H_1\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_hh_bbbb_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_1\to hh\to b\bar b b\bar b@f$.
 */
class Hobs_gg_H1_hh_bbbb_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_hh_bbbb_CMS13 constructor.
     */
    Hobs_gg_H1_hh_bbbb_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to b\bar b b\bar b)]_{\text{theo}} / [\sigma_{gg\to H_1}\cdot BR(H_1\to hh\to b\bar b b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_hh_gagabb_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_1\to hh\to \gamma\gamma  b\bar b@f$.
 */
class Hobs_pp_H1_hh_gagabb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_hh_gagabb_ATLAS13 constructor.
     */
    Hobs_pp_H1_hh_gagabb_ATLAS13(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to \gamma \gamma b b\bar b)]_{\text{theo}} / [\sigma_{pp\to H_1}\cdot BR(H_1\to hh\to  \gamma \gamma b\bar b)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_hh_gagabb_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_1\to hh\to \gamma\gamma  b\bar b@f$.
 */
class Hobs_pp_H1_hh_gagabb_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_hh_gagabb_CMS13 constructor.
     */
    Hobs_pp_H1_hh_gagabb_CMS13(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to \gamma \gamma b\bar b)]_{\text{theo}} / [\sigma_{pp\to H_1}\cdot BR(H_1\to hh\to \gamma \gamma b\bar b)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_hh_bbtautau_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_1\to hh\to b\bar b \tau \tau b@f$.
 */
class Hobs_pp_H1_hh_bbtautau_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_hh_bbtautau_ATLAS13 constructor.
     */
    Hobs_pp_H1_hh_bbtautau_ATLAS13(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to b\bar b \tau \tau)]_{\text{theo}} / [\sigma_{pp\to H_1}\cdot BR(H_1\to hh\to b\bar b \tau \tau)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_hh_bbtautau_1_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_1\to hh\to b\bar \tau \tau @f$.
 */
class Hobs_pp_H1_hh_bbtautau_1_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_hh_bbtautau_1_CMS13 constructor.
     */
    Hobs_pp_H1_hh_bbtautau_1_CMS13(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to b\bar \tau \tau)]_{\text{theo}} / [\sigma_{pp\to H_1}\cdot BR(H_1\to hh\to b\bar \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_hh_bbtautau_2_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_1\to hh\to b\bar b \tau \tau b@f$.
 */
class Hobs_pp_H1_hh_bbtautau_2_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_hh_bbtautau_2_CMS13 constructor.
     */
    Hobs_pp_H1_hh_bbtautau_2_CMS13(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to b\bar \tau \tau)]_{\text{theo}} / [\sigma_{pp\to H_1}\cdot BR(H_1\to hh\to b\bar \tau \tau)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_hh_bblnulnu_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_1\to hh\to b\bar b \ell \nu \ell \nu@f$.
 */
class Hobs_pp_H1_hh_bblnulnu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_hh_bblnulnu_CMS13 constructor.
     */
    Hobs_pp_H1_hh_bblnulnu_CMS13(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to b\bar b \ell \nu \ell \nu)]_{\text{theo}} / [\sigma_{pp\to H_1}\cdot BR(H_1\to hh\to b\bar b \ell \nu \ell \nu)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H1_hh_gagaWW_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_1\to hh\to \gamma \gamma WW @f$.
 */
class Hobs_gg_H1_hh_gagaWW_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H1_hh_gagaWW_ATLAS13 constructor.
     */
    Hobs_gg_H1_hh_gagaWW_ATLAS13(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to \gamma \gamma WW)]_{\text{theo}} / [\sigma_{gg\to H_1}\cdot BR(H_1\to hh\to \gamma \gamma WW)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_hZ_bbZ_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_3\to hZ\to b\bar b Z@f$.
 */
class Hobs_gg_H3_hZ_bbZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_hZ_bbZ_ATLAS8 constructor.
     */
    Hobs_gg_H3_hZ_bbZ_ATLAS8(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ\to b\bar b Z)]_{\text{theo}} / [\sigma_{gg\to H_3}\cdot BR(H_3\to hZ \to b\bar b Z)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_hZ_bbll_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_3\to hZ\to b\bar b \ell \ell@f$.
 */
class Hobs_gg_H3_hZ_bbll_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_hZ_bbll_CMS8 constructor.
     */
    Hobs_gg_H3_hZ_bbll_CMS8(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ\to b\bar b \ell \ell)]_{\text{theo}} / [\sigma_{gg\to H_3}\cdot BR(H_3\to hZ\to b\bar b \ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_hZ_tautauZ_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_3\to hZ\to \tau \tau Z@f$.
 */
class Hobs_gg_H3_hZ_tautauZ_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_hZ_tautauZ_ATLAS8 constructor.
     */
    Hobs_gg_H3_hZ_tautauZ_ATLAS8(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ\to \tau \tau Z)]_{\text{theo}} / [\sigma_{gg\to H_3}\cdot BR(H_3\to hZ\to \tau \tau Z)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_hZ_tautaull_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_3\to hZ\to \tau\tau\ell\ell@f$.
 */
class Hobs_gg_H3_hZ_tautaull_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_hZ_tautaull_CMS8 constructor.
     */
    Hobs_gg_H3_hZ_tautaull_CMS8(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ\to hZ\to \tau\tau\ell\ell)]_{\text{theo}} / [\sigma_{gg\to H_3}\cdot BR(H_3\to hZ\to hZ\to \tau\tau\ell\ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_hZ_bbZ_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_3\to hZ\to  b\bar b Z@f$.
 */
class Hobs_gg_H3_hZ_bbZ_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_hZ_bbZ_ATLAS13 constructor.
     */
    Hobs_gg_H3_hZ_bbZ_ATLAS13(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ\to  b\bar b Z)]_{\text{theo}} / [\sigma_{gg\to H_3}\cdot BR(H_3\to hZ\to  b\bar b Z)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H3_hZ_bbZ_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_3\to hZ\to  b\bar b Z@f$.
 */
class Hobs_bb_H3_hZ_bbZ_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_H3_hZ_bbZ_ATLAS13 constructor.
     */
    Hobs_bb_H3_hZ_bbZ_ATLAS13(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ\to  b\bar b Z)]_{\text{theo}} / [\sigma_{gg\to H_3}\cdot BR(H_3\to hZ\to  b\bar b Z)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_hZ_bbZ_1_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_3\to hZ\to  b\bar b Z@f$.
 */
class Hobs_gg_H3_hZ_bbZ_1_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_hZ_bbZ_1_CMS13 constructor.
     */
    Hobs_gg_H3_hZ_bbZ_1_CMS13(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ\to  b\bar b Z)]_{\text{theo}} / [\sigma_{gg\to H_3}\cdot BR(H_3\to hZ\to  b\bar b Z)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H3_hZ_bbZ_1_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_3\to hZ\to b\bar b Z@f$.
 */
class Hobs_bb_H3_hZ_bbZ_1_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_H3_hZ_bbZ_1_CMS13 constructor.
     */
    Hobs_bb_H3_hZ_bbZ_1_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ\to  b\bar b Z)]_{\text{theo}} / [\sigma_{b\bar b\to H_3}\cdot BR(H_3\to hZ\to  b\bar b Z)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_hZ_bbZ_2_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_3\to hZ\to b\bar b Z@f$.
 */
class Hobs_gg_H3_hZ_bbZ_2_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_hZ_bbZ_2_CMS13 constructor.
     */
    Hobs_gg_H3_hZ_bbZ_2_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ\to  b\bar b Z)]_{\text{theo}} / [\sigma_{gg\to H_3}\cdot BR(H_3\to hZ\to  b\bar b Z)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H3_hZ_bbZ_2_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_3\to hZ\to b\bar b Z@f$.
 */
class Hobs_bb_H3_hZ_bbZ_2_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_H3_hZ_bbZ_2_CMS13 constructor.
     */
    Hobs_bb_H3_hZ_bbZ_2_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ\to  b\bar b Z)]_{\text{theo}} / [\sigma_{b\bar b\to H_3}\cdot BR(H_3\to hZ\to  b\bar b Z)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H3_H1Z_bbll_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_3\to hZ\to b\bar b \ell \ell@f$.
 */
class Hobs_pp_H3_H1Z_bbll_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H3_H1Z_bbll_CMS8 constructor.
     */
    Hobs_pp_H3_H1Z_bbll_CMS8(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_3}\cdot BR^{\text{GM}}(H_3\to H_1 Z\to  b\bar b \ell \ell)]_{\text{theo}} / [\sigma_{b\bar b\to H_3}\cdot BR(H_3\to H_1 Z\to  b\bar b \ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_gg_H3_H1Z_bbll_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$gg\to H_3\to H_1 Z\to b\bar b \ell \ell@f$.
 */
class Hobs_gg_H3_H1Z_bbll_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_gg_H3_H1Z_bbll_ATLAS13 constructor.
     */
    Hobs_gg_H3_H1Z_bbll_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to H_1 Z\to  b\bar b \ell \ell)]_{\text{theo}} / [\sigma_{gg\to H_3}\cdot BR(H_3\to H_1 Z\to  b\bar b \ell \ell)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_bb_H3_H1Z_bbll_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$b\bar b\to H_13to H_1 Z \\to b\bar b \ell \ell@f$.
 */
class Hobs_bb_H3_H1Z_bbll_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_bb_H3_H1Z_bbll_ATLAS13 constructor.
     */
    Hobs_bb_H3_H1Z_bbll_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{b\bar b\to H_3}\cdot BR^{\text{GM}}(H_3\to H_1 Z\to  b\bar b \ell \ell)]_{\text{theo}} / [\sigma_{b\bar b\to H_3}\cdot BR(H_3\to H_1 Z\to  b\bar b \ell \ell)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H1_H3Z_bbll_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_1\to H_3 Z\to b\bar b \ell \ell@f$.
 */
class Hobs_pp_H1_H3Z_bbll_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H1_H3Z_bbll_CMS8 constructor.
     */
    Hobs_pp_H1_H3Z_bbll_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_3}\cdot BR^{\text{GM}}(H_3\to H_1 Z\to  b\bar b \ell \ell)]_{\text{theo}} / [\sigma_{pp\to H_3}\cdot BR(H_3\to H_1 Z\to  b\bar b \ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5_H3Z_bbll_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_1\to H_3 Z\to b\bar b \ell \ell@f$.
 */
class Hobs_pp_H5_H3Z_bbll_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H5_H3Z_bbll_CMS8 constructor.
     */
    Hobs_pp_H5_H3Z_bbll_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_5}\cdot BR^{\text{GM}}(H_5\to H_3 Z\to  b\bar b \ell \ell)]_{\text{theo}} / [\sigma_{pp\to H_5}\cdot BR(H_5\to H_3 Z\to  b\bar b \ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H3pm_taunu_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_3^\pm\to \tau \nu@f$.
 */
class Hobs_pp_H3pm_taunu_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H3pm_taunu_ATLAS8 constructor.
     */
    Hobs_pp_H3pm_taunu_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_3^\pm}\cdot BR^{\text{GM}}( H_3^\pm \to \tau \nu)]_{\text{theo}} / [\sigma_{pp\to  H_3^\pm}\cdot BR( H_3^\pm \to \tau \nu)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H3p_taunu_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_3^\pm\to \tau \nu@f$.
 */
class Hobs_pp_H3p_taunu_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H3p_taunu_CMS8 constructor.
     */
    Hobs_pp_H3p_taunu_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_3^\pm}\cdot BR^{\text{GM}}( H_3^\pm \to \tau \nu)]_{\text{theo}} / [\sigma_{pp\to  H_3^\pm}\cdot BR( H_3^\pm \to \tau \nu)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H3pm_taunu_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_3^\pm\to \tau \nu@f$.
 */
class Hobs_pp_H3pm_taunu_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H3pm_taunu_ATLAS13 constructor.
     */
    Hobs_pp_H3pm_taunu_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_3^\pm}\cdot BR^{\text{GM}}( H_3^\pm \to \tau \nu)]_{\text{theo}} / [\sigma_{pp\to  H_3^\pm}\cdot BR( H_3^\pm \to \tau \nu)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H3pm_taunu_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_3^\pm \to \tau \nu@f$.
 */
class Hobs_pp_H3pm_taunu_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H3pm_taunu_CMS13 constructor.
     */
    Hobs_pp_H3pm_taunu_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_3^\pm}\cdot BR^{\text{GM}}( H_3^\pm \to \tau \nu)]_{\text{theo}} / [\sigma_{pp\to  H_3^\pm}\cdot BR( H_3^\pm \to \tau \nu)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H3pm_tb_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_3^\pm \to tb@f$.
 */
class Hobs_pp_H3pm_tb_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H3pm_tb_ATLAS8 constructor.
     */
    Hobs_pp_H3pm_tb_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_3^\pm}\cdot BR^{\text{GM}}( H_3^\pm \to tb)]_{\text{theo}} / [\sigma_{pp\to  H_3^\pm}\cdot BR( H_3^\pm \to \tau \nu)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H3p_tb_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_3^\pm \to tb@f$.
 */
class Hobs_pp_H3p_tb_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H3p_tb_CMS8 constructor.
     */
    Hobs_pp_H3p_tb_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_3^\pm}\cdot BR^{\text{GM}}( H_3^\pm \to tb)]_{\text{theo}} / [\sigma_{pp\to  H_3^\pm}\cdot BR( H_3^\pm \to \tau \nu)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H3pm_tb_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$pp\to H_3^\pm \to tb @f$.
 */
class Hobs_pp_H3pm_tb_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H3pm_tb_ATLAS13 constructor.
     */
    Hobs_pp_H3pm_tb_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_3^\pm}\cdot BR^{\text{GM}}( H_3^\pm \to tb)]_{\text{theo}} / [\sigma_{pp\to  H_3^\pm}\cdot BR( H_3^\pm \to \tau \nu)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_WZ_H5pm_WZ_qqll_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$WZ\to H_5^\pm \to WZ \to qq \ell \ell @f$.
 */
class Hobs_WZ_H5pm_WZ_qqll_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_WZ_H5pm_WZ_qqll_ATLAS8 constructor.
     */
    Hobs_WZ_H5pm_WZ_qqll_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{WZ\to H_5^\pm}\cdot BR^{\text{GM}}( H_5^\pm \to WZ \to qq \ell \ell)]_{\text{theo}} / [\sigma_{WZ\to  H_5^\pm}\cdot BR( H_5^\pm \to qq \ell \ell)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_WZ_H5pm_WZ_lnull_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * times branching ratio of the process @f$WZ\to H_5^\pm \to WZ \to \ell \nu \ell \ell @f$.
 */
class Hobs_WZ_H5pm_WZ_lnull_ATLAS13: public ThObservable {
public:

     /**
     * @brief Hobs_WZ_H5pm_WZ_lnull_ATLAS13 constructor.
     */
    Hobs_WZ_H5pm_WZ_lnull_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{WZ\to H_5^\pm}\cdot BR^{\text{GM}}( H_5^\pm \to WZ \to \ell \nu \ell \ell)]_{\text{theo}} / [\sigma_{WZ\to  H_5^\pm}\cdot BR( H_5^\pm \to \ell \nu \ell \ell)]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Robs_WZ_H5pm_WZ_lnull_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Observable for the implementation of the ATLAS upper limit on signal 
 * strength of the process @f$WZ\to H_5^\pm\to WZ \to \ell \nu \ell \ell@f$ assuming a Gaussian likelihood.
 */
class Robs_WZ_H5pm_WZ_lnull_ATLAS13: public ThObservable {
public:

    /**
     * @brief Robs_WZ_H5pm_WZ_lnull_ATLAS13 constructor.
     */
    Robs_WZ_H5pm_WZ_lnull_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$1.96 + 1.96\left( [\sigma^{\text{GM}}_{WZ\to H_5^\pm}\cdot BR^{\text{GM}}(H_5^\pm\to WZ \to \ell \nu \ell \ell)]_{\text{theo}} - [\sigma^{\text{GM}}_{WZ\to H_5^\pm}\cdot BR^{\text{GM}}(H_5^\pm \to WZ \to \ell \nu \ell \ell)]_{\text{ATLAS,95\% observed}} \right) / [\sigma^{\text{THDM}}_{gg\to H}\cdot BR^{\text{THDM}}(H\to ZZ)]_{\text{ATLAS,95\% expected}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_WZ_H5pm_WZ_lnull_1_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$WZ\to H_5^\pm \to WZ \to \ell \nu \ell \ell @f$.
 */
class Hobs_WZ_H5pm_WZ_lnull_1_CMS13: public ThObservable {
public:

     /**
     * @brief Hobs_WZ_H5pm_WZ_lnull_1_CMS13 constructor.
     */
    Hobs_WZ_H5pm_WZ_lnull_1_CMS13(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{WZ\to H_5^\pm}\cdot BR^{\text{GM}}( H_5^\pm \to WZ \to \ell \nu \ell \ell)]_{\text{theo}} / [\sigma_{WZ\to  H_5^\pm}\cdot BR( H_5^\pm \to \ell \nu \ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_WZ_H5pm_WZ_lnull_2_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$WZ\to H_5^\pm \to WZ \to \ell \nu \ell \ell @f$.
 */
class Hobs_WZ_H5pm_WZ_lnull_2_CMS13: public ThObservable {
public:

     /**
     * @brief Hobs_WZ_H5pm_WZ_lnull_2_CMS13 constructor.
     */
    Hobs_WZ_H5pm_WZ_lnull_2_CMS13(const StandardModel& SM_i);

     /**
     * @return @f$[\sigma^{\text{GM}}_{WZ\to H_5^\pm}\cdot BR^{\text{GM}}( H_5^\pm \to WZ \to \ell \nu \ell \ell)]_{\text{theo}} / [\sigma_{WZ\to  H_5^\pm}\cdot BR( H_5^\pm \to \ell \nu \ell \ell)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5ppmmH5mmpp_eeee_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * of the process @f$pp\to H_5^\pm\pm H_5^\mp\mp @f$ measured by the decays @f$ H_5^\pm\pm \to ee @f$ .
 */
class Hobs_pp_H5ppmmH5mmpp_eeee_ATLAS8: public ThObservable {
public:

     /**
     * @brief Hobs_pp_H5ppmmH5mmpp_eeee_ATLAS8 constructor.
     */
    Hobs_pp_H5ppmmH5mmpp_eeee_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_5^\pm\pm H_5^\mp\mp}\cdot /[\sigma_{pp\to H_5^\pm\pm H_5^\mp\mp}]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5ppmmH5mmpp_emuemu_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * of the process @f$pp\to H_5^\pm\pm H_5^\mp\mp @f$ measured by the decays @f$ H_5^\pm\pm \to e \nu @f$ .
 */
class Hobs_pp_H5ppmmH5mmpp_emuemu_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H5ppmmH5mmpp_emuemu_ATLAS8 constructor.
     */
    Hobs_pp_H5ppmmH5mmpp_emuemu_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_5^\pm\pm H_5^\mp\mp}\cdot /[\sigma_{pp\to H_5^\pm\pm H_5^\mp\mp}]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5ppmmH5mmpp_mumumumu_ATLAS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * of the process @f$pp\to H_5^\pm\pm H_5^\mp\mp @f$ measured by the decays @f$ H_5^\pm\pm \to \mu \mu  @f$ .
 */
class Hobs_pp_H5ppmmH5mmpp_mumumumu_ATLAS8: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H5ppmmH5mmpp_mumumumu_ATLAS8 constructor.
     */
    Hobs_pp_H5ppmmH5mmpp_mumumumu_ATLAS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_5^\pm\pm H_5^\mp\mp}\cdot /[\sigma_{pp\to H_5^\pm\pm H_5^\mp\mp}]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5ppmmH5mmpp_llll_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * of the process @f$pp\to H_5^\pm\pm H_5^\mp\mp @f$ measured by the decays @f$ H_5^\pm\pm \to \ell \ell @f$ .
 */
class Hobs_pp_H5ppmmH5mmpp_llll_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H5ppmmH5mmpp_llll_ATLAS13 constructor.
     */
    Hobs_pp_H5ppmmH5mmpp_llll_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_5^\pm\pm H_5^\mp\mp}\cdot /[\sigma_{pp\to H_5^\pm\pm H_5^\mp\mp}]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_pp_H5ppmmH5mmpp_WWWW_ATLAS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and ATLAS upper limit for the cross section 
 * of the process @f$pp\to H_5^\pm\pm H_5^\mp\mp @f$ measured by the decays @f$ H_5^\pm\pm \to WW @f$ .
 */
class Hobs_pp_H5ppmmH5mmpp_WWWW_ATLAS13: public ThObservable {
public:

    /**
     * @brief Hobs_pp_H5ppmmH5mmpp_WWWW_ATLAS13 constructor.
     */
    Hobs_pp_H5ppmmH5mmpp_WWWW_ATLAS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{pp\to H_5^\pm\pm H_5^\mp\mp}\cdot /[\sigma_{pp\to H_5^\pm\pm H_5^\mp\mp}]_{\text{ATLAS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H5ppmm_WW_jjll_CMS8
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$VV\to H_5^\pm\pm \to WW @f$.
 */
class Hobs_VV_H5ppmm_WW_jjll_CMS8: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H5ppmm_WW_jjll_CMS8 constructor.
     */
    Hobs_VV_H5ppmm_WW_jjll_CMS8(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV\to H_5^\pm\pm}\cdot BR^{\text{GM}}( H_5^\pm\pm \to WW )]_{\text{theo}} / [\sigma_{VV\to  H_5^\pm\pm}\cdot BR( H_5^\pm\pm \to WW)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class Hobs_VV_H5ppmm_WW_jjll_CMS13
 * @ingroup GeorgiMachacek
 * @brief Ratio of the prediction and CMS upper limit for the cross section 
 * times branching ratio of the process @f$VV\to H_5^\pm\pm \to WW @f$.
 */
class Hobs_VV_H5ppmm_WW_jjll_CMS13: public ThObservable {
public:

    /**
     * @brief Hobs_VV_H5ppmm_WW_jjll_CMS8 constructor.
     */
    Hobs_VV_H5ppmm_WW_jjll_CMS13(const StandardModel& SM_i);

    /**
     * @return @f$[\sigma^{\text{GM}}_{VV\to H_5^\pm\pm}\cdot BR^{\text{GM}}( H_5^\pm\pm \to WW )]_{\text{theo}} / [\sigma_{VV\to  H_5^\pm\pm}\cdot BR( H_5^\pm\pm \to WW)]_{\text{CMS,95\%}}@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

 /**
 * @class log10_tt_H1_tt_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process 
  * @f$tt\to H_1\to tt @f$ at 13 TeV.
 */  
class log10_tt_H1_tt_TH13: public ThObservable {
public:

    /**
     * @brief log10_tt_H1_tt_TH13 constructor.
     */
    log10_tt_H1_tt_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{tt\to H_1}\cdot BR^{\text{GM}}(H_1\to tt)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

 /**
 * @class log10_bb_H1_tt_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process 
  * @f$b\bar b\to H_1\to tt @f$ at 13 TeV.
 */  
class log10_bb_H1_tt_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_H1_tt_TH13 constructor.
     */
    log10_bb_H1_tt_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to tt)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

 /**
 * @class log10_tt_H3_tt_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process 
  * @f$tt\to H_3\to tt @f$ at 13 TeV.
 */  
class log10_tt_H3_tt_TH13: public ThObservable {
public:

    /**
     * @brief log10_tt_H3_tt_TH13 constructor.
     */
    log10_tt_H3_tt_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{tt\to H_3}\cdot BR^{\text{GM}}(H_3\to tt)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

 /**
 * @class log10_bb_H3_tt_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process 
  * @f$b\bar b\to H_3\to tt @f$ at 13 TeV.
 */  
class log10_bb_H3_tt_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_H3_tt_TH13 constructor.
     */
    log10_bb_H3_tt_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{b\bar b\to H_3}\cdot BR^{\text{GM}}(H_3\to tt)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

 /**
 * @class log10_bb_H1_bb_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process 
  * @f$b\bar b\to H_1\to b\bar b @f$ at 8 TeV.
 */  
class log10_bb_H1_bb_TH8: public ThObservable {
public:

    /**
     * @brief log10_bb_H1_bb_TH8 constructor.
     */
    log10_bb_H1_bb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

 /**
 * @class log10_gg_H1_bb_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process 
  * @f$gg\to H_1\to b\bar b @f$ at 8 TeV.
 */  
class log10_gg_H1_bb_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_bb_TH8 constructor.
     */
    log10_gg_H1_bb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

 /**
 * @class log10_pp_H1_bb_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process 
  * @f$pp\to H_1\to b\bar b @f$ at 13 TeV.
 */  
class log10_pp_H1_bb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_bb_TH13 constructor.
     */
    log10_pp_H1_bb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

 /**
 * @class log10_bb_H1_bb_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process 
  * @f$b\bar b\to H_1\to b\bar b @f$ at 13 TeV.
 */  
class log10_bb_H1_bb_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_H1_bb_TH13 constructor.
     */
    log10_bb_H1_bb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

 /**
 * @class log10_bb_H3_bb_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process 
  * @f$b\bar b\to H_3\to b\bar b @f$ at 8 TeV.
 */  
class log10_bb_H3_bb_TH8: public ThObservable {
public:

    /**
     * @brief log10_bb_H3_bb_TH8 constructor.
     */
    log10_bb_H3_bb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{b\bar b\to H_3}\cdot BR^{\text{GM}}(H_3\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

 /**
 * @class log10_gg_H3_bb_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process 
  * @f$gg\to H_3\to b\bar b @f$ at 8 TeV.
 */  
class log10_gg_H3_bb_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H3_bb_TH8 constructor.
     */
    log10_gg_H3_bb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

 /**
 * @class log10_pp_H3_bb_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process 
  * @f$pp\to H_3\to b\bar b @f$ at 13 TeV.
 */  
class log10_pp_H3_bb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H3_bb_TH13 constructor.
     */
    log10_pp_H3_bb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_3}\cdot BR^{\text{GM}}(H_3\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

 /**
 * @class log10_bb_H3_bb_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process 
  * @f$b\bar b\to H_3\to b\bar b @f$ at 13 TeV.
 */  
class log10_bb_H3_bb_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_H3_bb_TH13 constructor.
     */
    log10_bb_H3_bb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{b\bar b\to H_3}\cdot BR^{\text{GM}}(H_3\to b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H1_tautau_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the process @f$gg\to H_1\to \tau\tau@f$ at 8 TeV.
 */
class log10_gg_H1_tautau_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_tautau_TH8 constructor.
     */
    log10_gg_H1_tautau_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};


/**
 * @class log10_bb_H1_tautau_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the 
 * process @f$b\bar b\to H_1\to \tau\tau@f$ at 8 TeV.
 */class log10_bb_H1_tautau_TH8: public ThObservable {
public:

     /**
     * @brief log10_bb_H1_tautau_TH8 constructor.
     */
    log10_bb_H1_tautau_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};


/**
 * @class log10_gg_H1_tautau_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the 
 * process @f$gg\to H_1\to \tau\tau@f$ at 13 TeV.
 */
class log10_gg_H1_tautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_tautau_TH13 constructor.
     */
    log10_gg_H1_tautau_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_bb_H1_tautau_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the 
 * process @f$b\bar b\to H_1\to \tau\tau@f$ at 13 TeV.
 */
class log10_bb_H1_tautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_H1_tautau_TH13 constructor.
     */
    log10_bb_H1_tautau_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{b\bar b\to H_1}\cdot BR^{\text{GM}}(H_1\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H3_tautau_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of 
 * the process @f$gg\to H_3\to \tau\tau@f$ at 8 TeV.
 */
class log10_gg_H3_tautau_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H3_tautau_TH8 constructor.
     */
    log10_gg_H3_tautau_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_bb_H3_tautau_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$b\bar b\to H_3\to \tau\tau@f$ at 8 TeV.
 */
class log10_bb_H3_tautau_TH8: public ThObservable {
public:

    /**
     * @brief log10_bb_H3_tautau_TH8 constructor.
     */
    log10_bb_H3_tautau_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{b\bar b\to H_3}\cdot BR^{\text{GM}}(H_3\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H3_tautau_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the 
 * process @f$gg\to H_3\to \tau\tau@f$ at 13 TeV.
 */
class log10_gg_H3_tautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_H3_tautau_TH13 constructor.
     */
    log10_gg_H3_tautau_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_bb_H3_tautau_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$b\bar b\to H_3\to \tau\tau@f$ at 13 TeV.
 */
class log10_bb_H3_tautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_H3_tautau_TH13 constructor.
     */
    log10_bb_H3_tautau_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{b\bar b\to H_3}\cdot BR^{\text{GM}}(H_3\to \tau\tau)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};


/**
 * @class log10_gg_H1_gaga_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_1\to \gamma\gamma@f$ at 8 TeV.
 */
class log10_gg_H1_gaga_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_gaga_TH8 constructor.
     */
    log10_gg_H1_gaga_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_gaga_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to \gamma\gamma@f$ at 13 TeV.
 */
class log10_pp_H1_gaga_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_gaga_TH13 constructor.
     */
    log10_pp_H1_gaga_TH13(const StandardModel& SM_i);

      /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H1_gaga_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_1\to \gamma\gamma@f$ at 13 TeV.
 */
class log10_gg_H1_gaga_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_gaga_TH13 constructor.
     */
    log10_gg_H1_gaga_TH13(const StandardModel& SM_i);

      /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H3_gaga_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_3\to \gamma\gamma@f$ at 8 TeV.
 */
class log10_gg_H3_gaga_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H3_gaga_TH8 constructor.
     */
    log10_gg_H3_gaga_TH8(const StandardModel& SM_i);

      /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H3_gaga_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_3\to \gamma\gamma@f$ at 13 TeV.
 */
class log10_pp_H3_gaga_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H3_gaga_TH13 constructor.
     */
    log10_pp_H3_gaga_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_3}\cdot BR^{\text{GM}}(H_3\to \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H3_gaga_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_3\to \gamma\gamma@f$ at 13 TeV.
 */
class log10_gg_H3_gaga_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_H3_gaga_TH13 constructor.
     */
    log10_gg_H3_gaga_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H5_gaga_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_5\to \gamma\gamma@f$ at 13 TeV.
 */
class log10_pp_H5_gaga_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H5_gaga_TH13 constructor.
     */
    log10_pp_H5_gaga_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_5}\cdot BR^{\text{GM}}(H_5\to \gamma \gamma)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};


/**
 * @class log10_pp_H1_Zga_llga_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to Z\gamma \to \ell \ell \gamma@f$ at 8 TeV.
 */
class log10_pp_H1_Zga_llga_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_Zga_llga_TH8 constructor.
     */
    log10_pp_H1_Zga_llga_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to Z \gamma  \to \ell \ell \gamma)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H1_Zga_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_1\to Z\gamma @f$ at 13 TeV.
 */
class log10_gg_H1_Zga_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_Zga_TH13 constructor.
     */
    log10_gg_H1_Zga_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to Z \gamma )]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H3_Zga_llga_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_3\to Z\gamma \to \ell \ell \gamma@f$ at 8 TeV.
 */
class log10_pp_H3_Zga_llga_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H3_Zga_llga_TH8 constructor.
     */
    log10_pp_H3_Zga_llga_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_3}\cdot BR^{\text{GM}}(H_3\to Z \gamma  \to \ell \ell \gamma)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H3_Zga_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_3\to Z\gamma@f$ at 13 TeV.
 */
class log10_gg_H3_Zga_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_H3_Zga_TH13 constructor.
     */
    log10_gg_H3_Zga_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to Z \gamma  )]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H5_Zga_llga_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_5\to Z\gamma \to \ell \ell \gamma@f$ at 8 TeV.
 */
class log10_pp_H5_Zga_llga_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H5_Zga_llga_TH8 constructor.
     */
    log10_pp_H5_Zga_llga_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_5}\cdot BR^{\text{GM}}(H_5\to Z \gamma  \to \ell \ell \gamma)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};


/**
 * @class log10_gg_H1_ZZ_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_1\to ZZ@f$ at 8 TeV.
 */
class log10_gg_H1_ZZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_ZZ_TH8 constructor.
     */
    log10_gg_H1_ZZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_VV_H1_ZZ_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$VV\to H_1\to ZZ@f$ at 8 TeV.
 */
class log10_VV_H1_ZZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_VV_H1_ZZ_TH8 constructor.
     */
    log10_VV_H1_ZZ_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{VV\to H_1}\cdot BR^{\text{GM}}(H_1\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H1_ZZ_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_1\to ZZ@f$ at 13 TeV.
 */
class log10_gg_H1_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_ZZ_TH13 constructor.
     */
    log10_gg_H1_ZZ_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_VV_H1_ZZ_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$VV\to H_1\to ZZ@f$ at 13 TeV.
 */
class log10_VV_H1_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_VV_H1_ZZ_TH13 constructor.
     */
    log10_VV_H1_ZZ_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{VV\to H_1}\cdot BR^{\text{GM}}(H_1\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_ZZ_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to ZZ@f$ at 13 TeV.
 */
class log10_pp_H1_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_ZZ_TH13 constructor.
     */
    log10_pp_H1_ZZ_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_VV_H5_ZZ_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$VV\to H_5\to ZZ@f$ at 8 TeV.
 */
class log10_VV_H5_ZZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_VV_H5_ZZ_TH8 constructor.
     */
    log10_VV_H5_ZZ_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{VV\to H_5}\cdot BR^{\text{GM}}(H_5\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_VV_H5_ZZ_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$VV\to H_5\to ZZ@f$ at 13 TeV.
 */
class log10_VV_H5_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_VV_H5_ZZ_TH13 constructor.
     */
    log10_VV_H5_ZZ_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{VV\to H_5}\cdot BR^{\text{GM}}(H_5\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H5_ZZ_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_5\to ZZ@f$ at 13 TeV.
 */
class log10_pp_H5_ZZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H5_ZZ_TH13 constructor.
     */
    log10_pp_H5_ZZ_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_5}\cdot BR^{\text{GM}}(H_5\to ZZ)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H1_WW_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_1\to WW@f$ at 8 TeV.
 */
class log10_gg_H1_WW_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_WW_TH8 constructor.
     */
    log10_gg_H1_WW_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_VV_H1_WW_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$VV\to H_1\to WW@f$ at 8 TeV.
 */
class log10_VV_H1_WW_TH8: public ThObservable {
public:

    /**
     * @brief log10_VV_H1_WW_TH8 constructor.
     */
    log10_VV_H1_WW_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{VV\to H_1}\cdot BR^{\text{GM}}(H_1\to WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H1_WW_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_1\to WW@f$ at 13 TeV.
 */
class log10_gg_H1_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_WW_TH13 constructor.
     */
    log10_gg_H1_WW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_VV_H1_WW_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$VV\to H_1\to WW@f$ at 13 TeV.
 */
class log10_VV_H1_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_VV_H1_WW_TH13 constructor.
     */
    log10_VV_H1_WW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{VV\to H_1}\cdot BR^{\text{GM}}(H_1\to WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_ggVV_H1_WW_lnulnu_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$(gg+vV)\to H_1\to WW  \to \el \nu \el \nu@f$ at 13 TeV.
 */
class log10_ggVV_H1_WW_lnulnu_TH13: public ThObservable {
public:

    /**
     * @brief log10_ggVV_H1_WW_lnulnu_TH13 constructor.
     */
    log10_ggVV_H1_WW_lnulnu_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{(gg+VV)\to H_1}\cdot BR^{\text{GM}}(H_1\to WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_WW_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to WW@f$ at 13 TeV.
 */
class log10_pp_H1_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_WW_TH13 constructor.
     */
    log10_pp_H1_WW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_VV_H5_WW_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$VV\to H_5\to WW@f$ at 8 TeV.
 */
class log10_VV_H5_WW_TH8: public ThObservable {
public:

    /**
     * @brief log10_VV_H5_WW_TH8 constructor.
     */
    log10_VV_H5_WW_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{VV\to H_5}\cdot BR^{\text{GM}}(H_5\to WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_VV_H5_WW_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$VV\to H_5\to WW@f$ at 13 TeV.
 */
class log10_VV_H5_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_VV_H5_WW_TH13 constructor.
     */
    log10_VV_H5_WW_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{VV\to H_5}\cdot BR^{\text{GM}}(H_5\to WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_ggVV_H5_WW_lnulnu_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$(gg+VV)\to H_5\to WW \to \ell \nu \ell \nu@f$ at 13 TeV.
 */
class log10_ggVV_H5_WW_lnulnu_TH13: public ThObservable {
public:

    /**
     * @brief log10_ggVV_H5_WW_lnulnu_TH13 constructor.
     */
    log10_ggVV_H5_WW_lnulnu_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{(gg+VV)\to H_5}\cdot BR^{\text{GM}}(H_5\to WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H5_WW_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_5\to WW@f$ at 13 TeV.
 */
class log10_pp_H5_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H5_WW_TH13 constructor.
     */
    log10_pp_H5_WW_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_5}\cdot BR^{\text{GM}}(H_5\to WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_VV_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to VV@f$ at 8 TeV.
 */
class log10_pp_H1_VV_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_VV_TH8 constructor.
     */
    log10_pp_H1_VV_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to VV)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_VV_H5_WW_TH8
 * @ingroup GeorgiMachacek
  * @brief Decadic logarithm of the signal strength of 
 * the process @f$pp\to H_1\to VV@f$ at 8 TeV.
 */
class log10_mu_pp_H1_VV_TH8: public ThObservable {
public:

    /**
     * @brief log10_mu_pp_H1_VV_TH8 constructor.
     */
    log10_mu_pp_H1_VV_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\mu_H_1^{\text{THDM}}(H_1\to VV)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_VV_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to VV@f$ at 13 TeV.
 */
class log10_pp_H1_VV_TH13: public ThObservable {
public:

     /**
     * @brief log10_pp_H1_VV_TH13 constructor.
     */
    log10_pp_H1_VV_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to VV)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H5_VV_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_5\to VV@f$ at 8 TeV.
 */
class log10_pp_H5_VV_TH8: public ThObservable {
public:

     /**
     * @brief log10_pp_H5_VV_TH8 constructor.
     */
    log10_pp_H5_VV_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_5}\cdot BR^{\text{GM}}(H_5\to VV)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_mu_pp_H5_VV_TH8
 * @ingroup GeorgiMachacek
  * @brief Decadic logarithm of the signal strength of the 
 * process @f$pp\to H_5\to VV@f$ at 8 TeV.
 */
class log10_mu_pp_H5_VV_TH8: public ThObservable {
public:

     /**
     * @brief log10_mu_pp_H5_VV_TH8 constructor.
     */
    log10_mu_pp_H5_VV_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\mu_H_5^{\text{THDM}}(H_5\to VV)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H5_VV_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_5\to VV@f$ at 13 TeV.
 */
class log10_pp_H5_VV_TH13: public ThObservable {
public:

     /**
     * @brief log10_pp_H5_VV_TH13 constructor.
     */
    log10_pp_H5_VV_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_5}\cdot BR^{\text{GM}}(H_5\to VV)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H1_hh_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$hh\to H_1\to hh@f$ at 8 TeV.
 */
class log10_gg_H1_hh_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_hh_TH8 constructor.
     */
    log10_gg_H1_hh_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to hh)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_hh_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to hh@f$ at 8 TeV.
 */
class log10_pp_H1_hh_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_hh_TH8 constructor.
     */
    log10_pp_H1_hh_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_hh_bbbb_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to hh \to  b\bar b  b\bar b@f$ at 8 TeV.
 */
class log10_pp_H1_hh_bbbb_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_hh_bbbb_TH8 constructor.
     */
    log10_pp_H1_hh_bbbb_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh to  b\bar b  b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_hh_gagabb_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to hh \to \gamma \gamma  b\bar b@f$ at 8 TeV.
 */
class log10_pp_H1_hh_gagabb_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_hh_gagabb_TH8 constructor.
     */
    log10_pp_H1_hh_gagabb_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh \to \gamma \gamma  b\bar b\to \gamma \gamma  b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H1_hh_bbtautau_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg \to H_1\to b\bar b \tau \tau @f$ at 8 TeV.
 */
class log10_gg_H1_hh_bbtautau_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_hh_bbtautau_TH8 constructor.
     */
    log10_gg_H1_hh_bbtautau_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to hh to b\bar b \tau \tau)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_hh_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to hh@f$ at 13 TeV.
 */
class log10_pp_H1_hh_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_hh_TH13 constructor.
     */
    log10_pp_H1_hh_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H1_hh_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_1\to hh@f$ at 13 TeV.
 */
class log10_gg_H1_hh_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_hh_TH13 constructor.
     */
    log10_gg_H1_hh_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to hh)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_hh_bbbb_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the 
 * process @f$pp\to H_1\to hh\to b\bar b b\bar b@f$ at 13 TeV.
 */
class log10_pp_H1_hh_bbbb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_hh_bbbb_TH13 constructor.
     */
    log10_pp_H1_hh_bbbb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh\to b\bar b b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H1_hh_bbbb_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_1\to hh \to  b\bar b b\bar b @f$ at 13 TeV.
 */
class log10_gg_H1_hh_bbbb_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_hh_bbbb_TH13 constructor.
     */
    log10_gg_H1_hh_bbbb_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to hh \to  b\bar b b\bar b )]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_hh_gagabb_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to hh \to \gamma \gamma b\bar b @f$ at 13 TeV.
 */
class log10_pp_H1_hh_gagabb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_hh_gagabb_TH13 constructor.
     */
    log10_pp_H1_hh_gagabb_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh \to \gamma \gamma b\bar b)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_hh_bbtautau_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to hh \to b\bar b \tau \tau@f$ at 13 TeV.
 */
class log10_pp_H1_hh_bbtautau_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_hh_bbtautau_TH13 constructor.
     */
    log10_pp_H1_hh_bbtautau_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh \to b\bar b \tau \tau)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_hh_bblnulnu_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to hh \to b\bar b \ell \nu \ell \nu \@f$ at 13 TeV.
 */
class log10_pp_H1_hh_bblnulnu_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_hh_bblnulnu_TH13 constructor.
     */
    log10_pp_H1_hh_bblnulnu_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to hh  \to b\bar b \ell \nu \ell \nu)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H1_hh_gagaWW_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_1\to hh \to \gamma \gamma WW@f$ at 13 TeV.
 */
class log10_gg_H1_hh_gagaWW_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_H1_hh_gagaWW_TH13 constructor.
     */
    log10_gg_H1_hh_gagaWW_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_1}\cdot BR^{\text{GM}}(H_1\to hh \to \gamma \gamma WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H3_hZ_bbZ_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_3  \to h Z \to b\bar b Z@f$ at 8 TeV.
 */
class log10_gg_H3_hZ_bbZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H3_hZ_bbZ_TH8 constructor.
     */
    log10_gg_H3_hZ_bbZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ \to b\bar bZ)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H3_hZ_bbll_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_3 \to h Z \to b\bar b \ell \ell@f$ at 13 TeV.
 */
class log10_gg_H3_hZ_bbll_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H3_hZ_bbll_TH8 constructor.
     */
    log10_gg_H3_hZ_bbll_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ \to b\bar b\ell \ell)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H3_hZ_tautauZ_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_3\to h Z \to \tau \tau b Z@f$ at 8 TeV.
 */
class log10_gg_H3_hZ_tautauZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H3_hZ_tautauZ_TH8 constructor.
     */
    log10_gg_H3_hZ_tautauZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ \to \tau \tau Z)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H3_hZ_tautaull_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_3 \to hZ \to \tau \tau \ell@f$ at 8 TeV.
 */
class log10_gg_H3_hZ_tautaull_TH8: public ThObservable {
public:

    /**
     * @brief log10_gg_H3_hZ_tautaull_TH8 constructor.
     */
    log10_gg_H3_hZ_tautaull_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ \to \tau \tau \ell \ell)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H3_hZ_bbZ_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$gg\to H_3\to h Z \to b\bar b Z@f$ at 13 TeV.
 */
class log10_gg_H3_hZ_bbZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_gg_H3_hZ_bbZ_TH13 constructor.
     */
    log10_gg_H3_hZ_bbZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{gg\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ \to  b\bar b Z)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_gg_H3_hZ_bbZ_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$b\bar b\to H_3\to h Z \to b\bar b Z@f$ at 13 TeV.
 */
class log10_bb_H3_hZ_bbZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_bb_H3_hZ_bbZ_TH13 constructor.
     */
    log10_bb_H3_hZ_bbZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{b\bar b\to H_3}\cdot BR^{\text{GM}}(H_3\to hZ \to  b\bar b Z)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H3_H1Z_bbll_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_3\to H_1 Z \to b\bar b \ell \ell@f$ at 8 TeV.
 */
class log10_pp_H3_H1Z_bbll_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H3_H1Z_bbll_TH8 constructor.
     */
    log10_pp_H3_H1Z_bbll_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_3}\cdot BR^{\text{GM}}(H_3\to H_1 Z \to b\bar b \ell \ell)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H3_H5Z_bbll_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_3\to H_5 Z \to b\bar b \ell \ell@f$ at 8 TeV.
 */
class log10_pp_H3_H5Z_bbll_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H3_H5Z_bbll_TH8 constructor.
     */
    log10_pp_H3_H5Z_bbll_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_3}\cdot BR^{\text{GM}}(H_3\to H_5 Z \to b\bar b \ell \ell)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H1_H3Z_bbll_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_1\to H_3 Z \to b\bar b \ell \ell@f$ at 8 TeV.
 */
class log10_pp_H1_H3Z_bbll_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H1_H3Z_bbll_TH8 constructor.
     */
    log10_pp_H1_H3Z_bbll_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_1}\cdot BR^{\text{GM}}(H_1\to H_3 Z \to b\bar b \ell \ell)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H5_H3Z_bbll_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_5\to H_3 Z \to b\bar b \ell \ell@f$ at 8 TeV.
 */
class log10_pp_H5_H3Z_bbll_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H5_H3Z_bbll_TH8 constructor.
     */
    log10_pp_H5_H3Z_bbll_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_5}\cdot BR^{\text{GM}}(H_5\to H_3 Z \to b\bar b \ell \ell)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H3pm_taunu_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_3^\pm\to \tau \nu@f$ at 8 TeV.
 */
class log10_pp_H3pm_taunu_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H3pm_taunu_TH8 constructor.
     */
    log10_pp_H3pm_taunu_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_3^\pm}\cdot BR^{\text{GM}}(H_3^\pm \to \tau \nu)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H3p_taunu_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_3^+\to \tau \nu@f$ at 8 TeV.
 */
class log10_pp_H3p_taunu_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H3p_taunu_TH8 constructor.
     */
    log10_pp_H3p_taunu_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_3^+}\cdot BR^{\text{GM}}(H_3^+\to \tau \nu)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H3pm_taunu_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_3^\pm\to \tau \nu@f$ at 13 TeV.
 */
class log10_pp_H3pm_taunu_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H3pm_taunu_TH13 constructor.
     */
    log10_pp_H3pm_taunu_TH13(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_3^\pm}\cdot BR^{\text{GM}}(H_3^\pm \to \tau \nu)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H3pm_tb_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_3^\pm\to tb@f$ at 8 TeV.
 */
class log10_pp_H3pm_tb_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H3pm_tb_TH8 constructor.
     */
    log10_pp_H3pm_tb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_3^\pm}\cdot BR^{\text{GM}}(H_3^\pm \to tb)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};


/**
 * @class log10_pp_H3p_tb_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_3^+\to tb@f$ at 8 TeV.
 */
class log10_pp_H3p_tb_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H3p_tb_TH8 constructor.
     */
    log10_pp_H3p_tb_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_3^+\cdot BR^{\text{GM}}(H_3^+ \to tb)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H3pm_tb_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_3^\pm\to tb@f$ at 13 TeV.
 */
class log10_pp_H3pm_tb_TH13: public ThObservable {
public:

    /**
     * @brief log10_pp_H3pm_tb_TH13 constructor.
     */
    log10_pp_H3pm_tb_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_3^\pm}\cdot BR^{\text{GM}}(H_3^\pm \to tb)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_WZ_H5pm_WZ_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$WZ\to H_5^\pm \to WZ@f$ at 8 TeV.
 */
class log10_WZ_H5pm_WZ_TH8: public ThObservable {
public:

    /**
     * @brief log10_WZ_H5pm_WZ_TH8 constructor.
     */
    log10_WZ_H5pm_WZ_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{WZ\to H_5^\pm}\cdot BR^{\text{GM}}(H_5^\pm \to WZ)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_WZ_H5pm_WZ_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$WZ\to H_5^\pm \to WZ@f$ at 13 TeV.
 */
class log10_WZ_H5pm_WZ_TH13: public ThObservable {
public:

    /**
     * @brief log10_WZ_H5pm_WZ_TH13 constructor.
     */
    log10_WZ_H5pm_WZ_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{WZ\to H_5^\pm}\cdot BR^{\text{GM}}(H_5^\pm \to WZ)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H5ppmmH5mmpp_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section  of the
 *  process @f$pp\to H_5^\pm\pm  H_5^\mp\mp@f$ at 8 TeV.
 */
class log10_pp_H5ppmmH5mmpp_TH8: public ThObservable {
public:

    /**
     * @brief log10_pp_H5ppmmH5mmpp_TH8 constructor.
     */
    log10_pp_H5ppmmH5mmpp_TH8(const StandardModel& SM_i);

     /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_5^\pm\pm H_5^\mp\mp}]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H5ppmmH5mmpp_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section  of the
 *  process @f$pp\to H_5^\pm\pm  H_5^\mp\mp @f$ at 13 TeV.
 */
class log10_pp_H5ppmmH5mmpp_TH13: public ThObservable {
public:

     /**
     * @brief log10_pp_H5ppmmH5mmpp_TH13 constructor.
     */
    log10_pp_H5ppmmH5mmpp_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to H_5^\pm\pm H_5^\mp\mp}]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_pp_H5ppmmH5mmpp_WWWW_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$pp\to H_5^\pm\pm  H_5^\mp\mp \to WWWW@f$ at 13 TeV.
 */
class log10_pp_H5ppmmH5mmpp_WWWW_TH13: public ThObservable {
public:

     /**
     * @brief log10_pp_H5ppmmH5mmpp_WWWW_TH13 constructor.
     */
    log10_pp_H5ppmmH5mmpp_WWWW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{pp\to  H_5^\pm\pm  H_5^\mp\mp}\cdot 2 BR^{\text{GM}}(H_5^\pm \to WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_VV_H5ppmm_WW_TH8
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$VV\to H_5^\pm\pm \to WW@f$ at 8 TeV.
 */
class log10_VV_H5ppmm_WW_TH8: public ThObservable {
public:

    /**
     * @brief log10_VV_H5ppmm_WW_TH8 constructor.
     */
    log10_VV_H5ppmm_WW_TH8(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{VV\to  H_5^\pm\pm }\cdot  BR^{\text{GM}}(H_5^\pm\pm  \to WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

/**
 * @class log10_VV_H5ppmm_WW_TH13
 * @ingroup GeorgiMachacek
 * @brief Decadic logarithm of the cross section times branching ratio of the
 *  process @f$VV\to H_5^\pm\pm \to WW@f$ at 13 TeV.
 */
class log10_VV_H5ppmm_WW_TH13: public ThObservable {
public:

    /**
     * @brief log10_VV_H5ppmm_WW_TH13 constructor.
     */
    log10_VV_H5ppmm_WW_TH13(const StandardModel& SM_i);

    /**
     * @return @f$\log_{10}[\sigma^{\text{GM}}_{VV\to  H_5^\pm\pm }\cdot  BR^{\text{GM}}(H_5^\pm\pm  \to WW)]@f$
     */
    double computeThValue();
private:
    const GeorgiMachacek& myGM;
};

#endif	/* GMDIRECTSEARCHES_H */

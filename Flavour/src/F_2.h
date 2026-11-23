/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

/* Code derived from arXiv:0810.4077
 * authors: Christoph Greub, Volker Pilipp, Christof Schüpbach
 */

#ifndef F_2_H
#define F_2_H

/**
 * @class F_2
 * @brief Two-loop matrix elements F_27, F_29 (and DeltaF) of b -> s l+ l- as functions of s^ = q^2/m_b^2 and z = m_c^2/m_b^2.
 * @details The expansion of hep-ph/0109140 is used for s^ < min(0.4, 4z) and that of arXiv:0810.4077 for s^ >= 4.5z.
 * Near the charm threshold, for min(0.4, 4z) <= s^ < 4.5z, neither expansion converges and all functions return NaN
 * (before, the two branches were switched at s^ = 0.4, where they differ by O(100%)). See the note in F_2.cpp.
 */
class F_2 {
public:
    /**
     * @brief The default constructor.
     */
    F_2();
    
    /**
     * @brief The default destructor.
     */
    ~F_2();

    double F_27re(double muh, double z, double sh, int maxpow = 20);

    double F_27im(double muh, double z, double sh, int maxpow = 20);

    double F_29re(double muh, double z, double sh, int maxpow = 20);

    double F_29im(double muh, double z, double sh, int maxpow = 20);

    double DeltaF_29re(double muh, double z, double sh, int maxpow = 20);

    double DeltaF_29im(double muh, double z, double sh, int maxpow = 20);

};

#endif /* F_2_H */
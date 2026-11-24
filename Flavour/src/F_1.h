/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

/* Code derived from arXiv:0810.4077
 * authors: Christoph Greub, Volker Pilipp, Christof Schüpbach
 */

#ifndef F_1_H
#define F_1_H

/**
 * @class F_1
 * @brief Two-loop matrix elements F_17, F_19 (and DeltaF) of b -> s l+ l- as functions of s^ = q^2/m_b^2 and z = m_c^2/m_b^2.
 * @details The expansion of hep-ph/0109140 is used for s^ < min(0.4, 4z) and that of arXiv:0810.4077 for s^ >= 4.5z.
 * Near the charm threshold, for min(0.4, 4z) <= s^ < 4.5z, neither expansion converges and all functions return NaN
 * (before, the two branches were switched at s^ = 0.4, where they differ by O(100%)). See the note in F_1.cpp.
 */
class F_1 {
public:

    /**
     * @brief The default constructor.
     */
    F_1();

    /**
     * @brief The default destructor.
     */
    ~F_1();

    double F_17re(double muh, double z, double sh, int maxpow = 20);

    double F_17im(double muh, double z, double sh, int maxpow = 20);

    double F_19re(double muh, double z, double sh, int maxpow = 20);

    double F_19im(double muh, double z, double sh, int maxpow = 20);

    double DeltaF_19re(double muh, double z, double sh, int maxpow = 20);

    double DeltaF_19im(double muh, double z, double sh, int maxpow = 20);

};
#endif /* F_1_H */

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

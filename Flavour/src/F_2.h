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
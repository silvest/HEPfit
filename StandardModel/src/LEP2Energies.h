/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2ENERGIES_H
#define LEP2ENERGIES_H

namespace LEP2Energies {

struct Energy {
    int label;
    double sqrtS; // GeV
};

constexpr Energy inclusive[] = {
    {130, 130.2}, {136, 136.2}, {161, 161.3}, {172, 172.1},
    {183, 182.7}, {189, 188.6}, {192, 191.6}, {196, 195.5},
    {200, 199.5}, {202, 201.6}, {205, 204.9}, {207, 206.7}
};

constexpr Energy heavyFlavour[] = {
    {133, 133.2}, {167, 166.6}, {183, 182.7}, {189, 188.6},
    {192, 191.6}, {196, 195.5}, {200, 199.5}, {202, 201.6},
    {205, 204.9}, {207, 206.7}
};

} // namespace LEP2Energies

#endif // LEP2ENERGIES_H

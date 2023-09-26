/*
 * Copyright (C) 2023 HEPfit Collaboration
 * 
 * 
 * For the licensing terms see doc/COPYING.
 */

#include "OrderScheme.h"

orders getHighest(orders order){
    if (order < FULLNLO) return order;
    return orders(order - FULLNLO + NLO);
}


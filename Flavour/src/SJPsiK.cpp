/*
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "SJPsiK.h"

double SJPsiK::getThValue() {
    return sin(-AmpBd(FULLNLO).arg());
}

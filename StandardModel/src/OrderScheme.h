/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef ORDERSCHEME_H
#define	ORDERSCHEME_H

#define MAXORDER NNLO
#define MAXORDER_EW NLO_ew

enum schemes {NDR=0, HV, LRI};
enum orders {LO=0, NLO, FULLNLO, NNLO, FULLNNLO};
enum orders_ew {NULL_ew = orders(MAXORDER+1), LO_ew, NLO_ew};

#endif	/* ORDERSCHEME_H */


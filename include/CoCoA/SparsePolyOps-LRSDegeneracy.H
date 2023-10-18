#ifndef CoCoA_SparsePolyOps_LRSDegeneracy_H
#define CoCoA_SparsePolyOps_LRSDegeneracy_H

//   Copyright (c)  2022  John Abbott,  Anna M. Bigatti
//   Original author: Nico Mexis (2022)

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.


#include "CoCoA/ring.H"
#include "CoCoA/VerificationLevel.H"

#include <vector>
// using std::vector;


namespace CoCoA
{

  /**
   * Algorithm 1 (naive factor-approach)
   * OBSOLESCENT: NOT UPDATED AND NOT SUPPORTED!
   */
//  unsigned long IsLRSDegenerateNaive(RingElem f);

  /**
   * Algorithm 2 (Graeffe iteration)
   * OBSOLESCENT: NOT UPDATED AND NOT SUPPORTED!
   */
//  unsigned long IsLRSDegenerateIter(RingElem f);

  /**
   * Checks for LRS-degeneracy of the given order using the given verification level
   */
  bool IsLRSDegenerateOrderMod(RingElem f, const unsigned long n, VerificationLevel VerLev);

  /**
   * Checks for LRS-degeneracy of the given order
   */
  bool IsLRSDegenerateOrder(ConstRefRingElem f, const unsigned long n);

  /**
   * Modular approach using the given verification level
   */
  unsigned long IsLRSDegenerateMod(RingElem f, VerificationLevel VerLev);

  /**
   * Chooses the right approach for the given polynomial
   */
  unsigned long LRSDegeneracyOrder(ConstRefRingElem f);
  std::vector<unsigned long> LRSDegeneracyOrders(ConstRefRingElem f);

  bool IsLRSDegenerate(RingElem f);

} // end of namespace CoCoA

// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-LRSDegeneracy.H,v 1.6 2023/03/06 22:06:56 abbott Exp $
// $Log: SparsePolyOps-LRSDegeneracy.H,v $
// Revision 1.6  2023/03/06 22:06:56  abbott
// Summary: Added new fn LRSDegeneracyOrders (returns a vector)
//
// Revision 1.5  2023/02/23 20:49:10  abbott
// Summary: Added new fn LRSDegeneracyOrder (was IsLRSDegenerate); major code revision
//
// Revision 1.4  2023/01/31 12:33:34  abbott
// Summary: Added IsLRSDegenerateOrderMod; some tidying
//
// Revision 1.3  2022/11/30 15:19:43  abbott
// Summary: Removed unused flag arg
//
// Revision 1.2  2022/11/24 21:11:29  abbott
// Summary: Added copyright notice
//


#endif

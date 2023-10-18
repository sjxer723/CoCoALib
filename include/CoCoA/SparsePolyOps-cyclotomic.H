#ifndef CoCoA_SparsePolyOps_cyclotomic_H
#define CoCoA_SparsePolyOps_cyclotomic_H

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


#include "CoCoA/BigInt.H"
#include "CoCoA/factorization.H"
#include "CoCoA/ring.H"

#include <vector>
// using std::vector;

namespace CoCoA
{

  RingElem cyclotomic(long n, ConstRefRingElem x);  ///< need n >= 1 and x must be an indet

  //----------------------
  // Below are fns to detect if a poly is cyclo

  // If f is (monic) cyclo, returns its index; o/w returns 0
  unsigned long CyclotomicTest(ConstRefRingElem f);  // 

  unsigned long CyclotomicIndex(ConstRefRingElem f);  // if f cyclo, return its index; o/w return some integer


  // Returns a factorization containing all cyclotomic factors of f
  // Factorization is empty iff f has no cyclotomic factors.
  // Internally uses variant of Smyth+Beukers.
  factorization<RingElem> CyclotomicFactors(ConstRefRingElem f);
  std::vector<long> CyclotomicFactorIndexes(ConstRefRingElem f);


} // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/include/CoCoA/SparsePolyOps-cyclotomic.H,v 1.6 2023/06/06 19:20:59 abbott Exp $
// $Log: SparsePolyOps-cyclotomic.H,v $
// Revision 1.6  2023/06/06 19:20:59  abbott
// Summary: Now have CyclotomicTest and CyclotomicIndex
//
// Revision 1.5  2023/05/11 19:27:13  abbott
// Summary: Added CycloFactorIndexes; some other temp changes
//
// Revision 1.4  2023/03/15 21:32:36  abbott
// Summary: Renamed CycloFactors to CyclotomicFactors; also refined semantics
//
// Revision 1.3  2023/02/23 20:50:47  abbott
// Summary: Added new fn CycloFactors (poor name)
//
// Revision 1.2  2023/01/31 12:32:56  abbott
// Summary: Moved CyclotomicIndex (& unchecked version) here from SparsePolyOps-Cyclotomicity
//
// Revision 1.1  2023/01/01 11:33:43  abbott
// Summary: Improved impl of cyclotomic; also cyclotomic moved to new file
//
//


#endif


//   Copyright (c)  2022  John Abbott and Anna M. Bigatti
//   Code developed from original CoCoA-3 code by Fabio Rossi (1999)
//   Major contributions from Nico Mexis (2022-2023)

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


#include "CoCoA/SparsePolyOps-cyclotomic.H"
#include "CoCoA/BigInt.H"
#include "CoCoA/BigIntOps.H"
#include "CoCoA/DUPFp.H"
#include "CoCoA/error.H"
#include "CoCoA/factor.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/NumTheory-factor.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/NumTheory-misc.H"
#include "CoCoA/NumTheory-modular.H"
#include "CoCoA/PPMonoid.H"
#include "CoCoA/ring.H"
#include "CoCoA/RingElemOps-CoprimeFactorBasis.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SmallFpImpl.H"
#include "CoCoA/SparsePolyIter.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/verbose.H"

#include "CoCoA/time.H"

#include <functional>
//#include <unordered_map>
// using std::unordered_map;
// #include <vector>
using std::vector;

namespace CoCoA
{

  namespace // anonymous
  {
    // Returns const ref to a table of height bounds for first few coeffs of cyclo polys:
    // coeff of x^k has height at most tbl[k]
    const std::vector<int>& CyclotomicCoeffHeightTable()
    {
      // Table was computed using program ..... ???? How to specify this?
      // Currently entries up to deg 120
      static vector<int> tbl{1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 4, 4, 5, 4, 4, 4, 5, 5, 6, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 9, 9, 7, 8, 8, 10, 13, 12, 10, 12, 9, 11, 15, 13, 13, 14, 15, 13, 16, 15, 15, 14, 16, 24, 17, 21, 21, 16, 22, 28, 26, 23, 28, 26, 25, 35, 34, 33, 28, 34, 36, 37, 49, 43, 33, 44, 48, 49, 55, 53, 53, 48, 60, 70, 66, 65, 70, 65, 68, 91, 86, 78, 87, 86, 86, 109, 110, 98, 104, 108, 116, 124, 136, 136};
      return tbl;
    }
    
  } // end of namespace anonymous

  namespace // anonymous
  {

  // SUPER DODGY CODE -- fiddles with GMP internal repr!!!
  void AddShift(BigInt& n, long c, unsigned long shift)
  {
    if (c==0) return;
    mpz_t& repr = mpzref(n);
    constexpr unsigned int BitsPerLimb = 64;  // FIX ME
    const long LimbSkip = shift/BitsPerLimb;
    const unsigned int ExtraBits = shift - BitsPerLimb*LimbSkip;
//    std::clog<<"LimbSkip="<<LimbSkip<<endl;
//    std::clog<<"ExtraBits="<<ExtraBits<<endl;
    BigInt C(c);
    mpz_mul_2exp(mpzref(C), mpzref(C), ExtraBits);
    repr[0]._mp_d += LimbSkip;
    repr[0]._mp_size -= LimbSkip; // better be positive!
    n += C;//     mpz_add(mpzref(n), mpzref(n), mpzref(C));
    repr[0]._mp_size += LimbSkip;
    repr[0]._mp_d -= LimbSkip;
  }


  BigInt EvalAt2_dodgy(const vector<long>& c)
  {
    const long n = len(c);
    BigInt N = power(2,n-1); // assume LC=1
    for (long j=n-2; j >= 0; --j)
    {
      if (c[j] != 0)
        AddShift(N, c[j], j);
    }
    return N;
  }

    // "binary" evaluation of poly in ZZ[x] at integer pt x=a
    BigInt EvalAt2_bin(vector<BigInt> C, int exp)
    {
      while (len(C) > 1)
      {
        if (IsOdd(len(C))) C.push_back(BigInt(0));
        BigInt pwr2(2);
//        BigInt tmp;
        long n = len(C);
        for (long i=0; i < n/2; ++i)
        {
///          C[i] = C[2*i]+a*C[2*i+1];  // maybe use mpz_add_mul  ??? worth it ???
          mpz_addmul(mpzref(C[2*i]), mpzref(C[2*i+1]), mpzref(pwr2));
          mpz_swap(mpzref(C[i]), mpzref(C[2*i]));
///          mpz_mul_2exp(mpzref(C[2*i+1]), mpzref(C[2*i+1]), exp);
///          mpz_add(mpzref(C[2*i+1]), mpzref(C[2*i+1]), mpzref(C[2*i]));
///          mpz_swap(mpzref(C[i]), mpzref(C[2*i+1])); // really assignment
        }
        C.resize(n/2);
        pwr2 *= pwr2;
//        exp *= 2;
      }
      return C[0];
    }


    // "binary" evaluation of poly in ZZ[x] at integer pt x=a
    BigInt EvalAt_bin(vector<BigInt> C, BigInt a)
    {
      while (len(C) > 1)
      {
        if (IsOdd(len(C))) C.push_back(BigInt(0));
        long n = len(C);
        for (long i=0; i < n/2; ++i)
          C[i] = C[2*i]+a*C[2*i+1];  // maybe use mpz_add_mul  ??? worth it ???
        C.resize(n/2);
        a *= a;
      }
      return C[0];
    }


    // "binary" evaluation of poly in ZZ[x] at rational pt x=a/b
    // RETURNS NUMERATOR!!
    BigInt EvalAt_bin(vector<BigInt> C, BigInt a, BigInt b)
    {
      if (IsOne(b)) return EvalAt_bin(C, a);
      BigInt ExcessFactor(1);
      while (len(C) > 1)
      {
        if (IsOdd(len(C)))  { C.push_back(BigInt(0));  ExcessFactor *= b; }
        long n = len(C);
        for (long i=0; i < n/2; ++i)
          C[i] = b*C[2*i]+a*C[2*i+1];
        C.resize(n/2);
        a *= a;
        b *= b;
      }
      return C[0]/ExcessFactor;
    }

    
    // Eval f(n)   ASSUMES f is univariate with integer coeffs
    BigInt EvalAt2_dodgy(ConstRefRingElem f)
    {
      vector<long> C(1+deg(f));
      for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
      {
        C[deg(PP(it))] = ConvertTo<long>(coeff(it));
      }
      return EvalAt2_dodgy(C); // 1 is exponent: i.e. 2 = 2^1
    }

    // Eval f(n)   ASSUMES f is univariate with integer coeffs
    BigInt EvalAt(ConstRefRingElem f, long n)
    {
      vector<BigInt> C(1+deg(f));
      for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
      {
        C[deg(PP(it))] = ConvertTo<BigInt>(coeff(it));
      }
//      if (n == 2) return EvalAt2_bin(C,1); // 1 is exponent: i.e. 2 = 2^1
      return EvalAt_bin(C, BigInt(n));
    }


    // evaluation of poly in ZZ[x] at rational pt x=a/b
    // RETURNS NUMERATOR!!
    BigInt EvalAt(ConstRefRingElem f, long n, long d)
    {
      vector<BigInt> C(1+deg(f));
      for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
      {
        C[deg(PP(it))] = ConvertTo<BigInt>(coeff(it));
      }
      return EvalAt_bin(C, BigInt(n), BigInt(d));
    }



    // The bounds below were taken from a poster by Arnold+Monagan (CECM 2009)
    // Result is a prime modulus greater than 2*CoeffHeight(cyclotomic(n)).
    // Here index is the "index" of the cyclotomic poly.
    long long cyclotomic_modulus(long index)
    {
      if (index < 40755L) return 127LL;
      if (index < 327845L) return 2371LL;
      if (index < 707455L) return 62039LL;
      if (index < 1181895L) return 119633LL;
      if (index < 10163195L) return 150000001LL;
      if (index < 40324935L) return 5398416817471LL;
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "cyclotomic_modulus");
    }


    /**
     * Returns the highest coefficient a cyclotomic polynomial
     * of squarefree index and given degree could possibly have.
     * The comment after the line corresponds to the cyclotomic
     * polynomial with the smallest degree that breaks the
     * respective case.
     */
    long CycloCoeffHeightBound(long d)
    {
      CoCoA_ASSERT(d > 1);
      if (d == 1) return 1;
      if (IsOdd(d)) return 0;
      if (d%8 != 0) return 1; // at most 2 distinct primes [!!CHECK!!  well-known]
//      if (d%16 == 8) return (2.0/3.0)*cbrt(d);  // at most 3 distinct primes "corrected sister Beiter conjecture"  SEE  arxiv:0910.2770.  Needs fixing?

      // Values below here have been generated using a
      // brute-force approach and are definitely correct.
      if (d < 48) return 1; // Phi(105)
      if (d < 240) return 2; // Phi(385)
      if (d < 576) return 3; // Phi(1365)
      if (d < 768) return 4; // Phi(1785)
      if (d < 1280) return 5; // Phi(2805)
      if (d < 1440) return 6; // Phi(3135)
      if (d < 3840) return 7; // Phi(6545)
      if (d < 5760) return 9; // Phi(15015)
      if (d < 8640) return 23; // Phi(21945)
      if (d < 10368) return 25; // Phi(25935)
      if (d < 10560) return 27; // Phi(26565)
      if (d < 17280) return 59; // Phi(40755)
      if (d < 50688) return 359; // Phi(106743)
      if (d < 82944) return 397; // Phi(171717)
      if (d < 92160) return 434; // Phi(255255)
      // Values below here have been generated using data
      // by Arnold and Monagan, see here:
      // http://wayback.cecm.sfu.ca/~ada26/cyclotomic/
      // TODO: Data below is correct if and only if all heights are correct
      if (d < 103680) return 532; // Phi(285285)
      if (d < 126720) return 1182; // Phi(345345)
      if (d < 138240) return 1311; // Phi(373065)
      if (d < 193536) return 5477; // Phi(327845)
      if (d < 387072) return 31010; // Phi(983535)
      if (d < 483840) return 59518; // Phi(1181895)
      if (d < 725760) return 14102773; // Phi(1752465)
      if (d < 1824768) return 14703509; // Phi(3949491)
      if (d < 3732480) return 56938657; // Phi(10555545)
      if (d < 4147200) return 88835350; // Phi(11565015)
      if (d < 4354560) return 197756850; // Phi(12267255)
      if (d < 5806080) return 310102051; // Phi(10163195)
      if (d < 7741440) return 1376877780831; // Phi(13441645)
      if (d < 8709120) return 1475674234751; // Phi(15069565)
      if (d < 11612160) return 1666495909761; // Phi(30489585)
      // // TODO: Data below is most definitely incomplete
      //   if (d < 15482880) return 2201904353336; // Phi(40324935)
      //   if (d < 17418240) return 2699208408726; // Phi(43730115)
      //   if (d < 104509440) return 862550638890874931; // Phi(306110805)
      //   if (d < 231469056) return 4722828832054556497; // Phi(497111433)
      //   if (d < 237828096) return 8171111062118177960; // Phi(516742863)
      //   if (d < 240869376) return 8768227953282038629; // Phi(522080013)
      //   if (d < 268240896) return 9038988754691465073; // Phi(582923523)
      //   if (d < 321159168) return 9118090447189969651; // Phi(693722757)
      //   if (d < 1072341504) return 9164566312887510757; // Phi(2583303555)
      CoCoA_THROW_ERROR(ERR::ArgTooBig,"CycloCoeffHeightBound");
    }


    /**
     * Returns the highest coefficient a cyclotomic polynomial
     * of squarefree index <= k.
     * The comment after the line corresponds to the cyclotomic
     * polynomial with the smallest degree that breaks the
     * respective case.
     */
    long CycloCoeffHeightBound_index(long k)
    {
      CoCoA_ASSERT(k > 0);
      if (IsEven(k)) k /= 2;
      if (k == 1) return 1;
      // Values below here have been generated using a
      // brute-force approach and are definitely correct.
      if (k < 105) return 1;
      if (k < 385) return 2;
      if (k < 1365) return 3;
      if (k < 1785) return 4;
      if (k < 2805) return 5;
      if (k < 3135) return 6;
      if (k < 6545) return 7;
      if (k < 10465) return 9;
      if (k < 11305) return 14;
      if (k < 17255) return 23;
      if (k < 20615) return 25;
      if (k < 26565) return 27;
      if (k < 40755) return 59;
      if (k < 106743) return 359;
      if (k < 171717) return 397;
      if (k < 255255) return 434;
      // Values below here have been generated using data by Arnold and Monagan
      // see here:   http://wayback.cecm.sfu.ca/~ada26/cyclotomic/
      // TODO: Data below is correct if and only if all heights are correct
      if (k < 279565) return 532;
      if (k < 327845) return 1182;
      if (k < 707455) return 31010;
      if (k < 886445) return 35111;
      if (k < 983535) return 44125;
      if (k < 1181895) return 59518;
      if (k < 1752465) return 14102773;
      if (k < 3949491) return 14703509;
      if (k < 8070699) return 56938657;
      if (k < 10163195) return 74989473;
      if (k < 13441645) return 1376877780831;
      if (k < 15069565) return 1475674234751;
      if (k < 30489585) return 1666495909761;
      if (k < 40324935) return 2201904353336;
      if (k < 43730115) return 2699208408726;
      if (k < 306110805) return 862550638890874931;
      if (k < 497111433) return 4722828832054556497;
      if (k < 516742863) return 8171111062118177960;
      if (k < 522080013) return 8768227953282038629;
      if (k < 582923523) return 9038988754691465073;
      if (k < 693722757) return 9118090447189969651;
      if (k < 2583303555) return 9164566312887510757; // 64-bit limit!
      CoCoA_THROW_ERROR(ERR::ArgTooBig,"CycloCoeffHeightBound");
//???      return 0; // means index too big for current table
    }



    namespace // anonymous
    {

      // !!! SUGGESTION !!!
      // Make this into a template fn: replace "unsigned long" by T;
      // change onlt type decl of cyclo, and return type.
      // Then we can instantiate with unsigned int types or BigInt!
      
      // Compute the first Dmax+1 coeffs of the cyclotomic poly (in a std::vector)
      // Dmax == 0 implies use half the degree; also if Dmax < 0 then
      // we compute the prefix of reciprocal of the requested cyclo
      std::vector<unsigned long> CycloPrefix(const std::vector<long>& plist,  long Dmax=0)
      {
        // assume Dmax == 0 || std::abs(Dmax) > 1  (maybe higher LWB???)
        // assume plist not empty, entries are DISTINCT PRIMES -- NOT CHECKED!!!
        vector<long> even;
        vector<long> odd;
        even.push_back(1);
        long phi = 1;
        bool Prime2Present = false;
        // This loop computes all products of even/odd subsets of plist; prime 2 is "ignored"
        for (int i=0; i < len(plist); ++i)
        {
          const long p = plist[i];
          if (Dmax != 0 && p > std::abs(Dmax)) continue;
          if (p == 2) { Prime2Present = true; continue; }
          phi *= (p-1);
          const long OVERFLOWz = std::abs(Dmax)/p; // "OVERFLOW" used in MacOS
          const int len_even = len(even);
          const int len_odd = len(odd);
          for (int j=0; j < len_even; ++j)
            if (Dmax==0 || even[j] <= OVERFLOWz)
              odd.push_back(p*even[j]);
          for (int j=0; j < len_odd; ++j)
            if (Dmax==0 || odd[j] <= OVERFLOWz)
              even.push_back(p*odd[j]);
        }
        const bool SwapOver = (Dmax < 0) ^ Prime2Present ^ IsOdd(len(plist));

        if ( SwapOver )
          swap(odd, even);
        // Now "even" corr to numerator; "odd" corr to denominator
        sort(even.begin(), even.end());
        const long UPB = (Dmax==0) ? phi/2 : std::abs(Dmax); //??? or maybe std::min(std::abs(Dmax),phi/2)
        vector<unsigned long> cyclo(UPB+1); // unsigned important to avoid UB upon overflow
        cyclo[0] = 1;
        long d = 0; // current deg, but at most UPB
        // This loop computes numerator product:
        for (const long i: even)
        {
          if (i > UPB) continue; // ??? probably no longer necessary
          // mult by (1+x^i)
          long top = d;
          if (top+i > UPB)
            top = UPB-i; // >= 0 by check above
          for (long j=top; j >= 0; --j)
            cyclo[j+i] += cyclo[j];
          d += i; if (d > UPB) d = UPB; // EQUIV:  d = min(d+i, UPB);
        }
        // Map x |-> -x
        for (long i=1; i <= UPB; i += 2)
          cyclo[i] = -cyclo[i];

        // This loop "divides" by denominator, one factor at a time.
        // Division is multiplication by truncasted inverse power series:
        // e.g. 1/(1-x) = 1+x+x^2+x^3+...
        for (const long i: odd)
        {
          if (i > UPB) continue; // ??? probably no longer necessary
          // div by (1-x^i)
          long top = UPB - i; // >= 0 by check above
          for (long j=0; j <= top; ++j)
            cyclo[j+i] += cyclo[j];
        }
        if (Prime2Present)
        {
          // Map x |-> -x
          for (long i=1; i <= UPB; i += 2)
            cyclo[i] = -cyclo[i];
        }
        return cyclo;    
      }
  

    } // end of namespace anonymous





    // Compute coeff vec of cyclotomic(n) mod modulus, where n = product(plist)
    // plist is a list of odd primes in incr order.
    // Let v be the result then v[k] is coeff of x^k in cyclotomic(n) mod modulus.
    // Entries in v are symmetric remainders.
    std::vector<long> cyclotomic_modp(long modulus, const vector<long>& plist)
    {
      // assume plist contains ODD primes in incr order
      if (plist.empty()) return vector<long>(2,1);
      long n=1; for (long p: plist) n *= p;
      SmallFpImpl ModP(modulus);
      DUPFp OnePoly(0,ModP);  AssignOne(OnePoly);
      DUPFp tmp(n,ModP);

      long p = plist[0];
//    for (long k=0; k < p; ++k) ShiftAdd(tmp, OnePoly,one(SmallFp),k*(n/p));
      for (long k=p-1; k >= 0; --k) ShiftAdd(tmp, OnePoly,one(SmallFp),k*(n/p));

      DUPFp phi = tmp;
      long i=1;
      while (i < len(plist))
      {
        CheckForInterrupt("cyclotomic_modp: gcd loop");
        const long p = plist[i];
        AssignZero(tmp);
//      for (long k=0; k < p; ++k) ShiftAdd(tmp, OnePoly,one(SmallFp),k*(n/p));
        for (long k=p-1; k >= 0; --k) ShiftAdd(tmp, OnePoly,one(SmallFp),k*(n/p));
        phi = gcd(phi, tmp); // may take a long time
        ++i;
      }
      const long d = deg(phi);
      vector<long> ans(1+d);
      for (long i=0; i <= d; ++i)
        ans[i] = ModP.myExportSymm(phi.myCoeffs[i]);
      return ans;
    }



    // Naive impl of CRT reconstruction -- probably good enough here
    long CRT2(long r1, long m1, long r2, long m2)
    {
      long m1inv = InvMod(m1,m2);
      long m2inv = InvMod(m2,m1);
      long R = m2*SymmRemainder(r1*m2inv, m1) + m1*SymmRemainder(r2*m1inv, m2);
      return SymmRemainder(R, m1*m2);
    }



  // Impl below is "curious": it computes the coeffs as machine "long int"
  // then maps the result into the given ring.

  // Developed from CoCoA code by Fabio Rossi (date 1999?)
  RingElem cyclotomic_FabioRossi(long n, const RingElem& x)
  {
    if (n < 1)
      CoCoA_THROW_ERROR("arg 1: n must be >= 1","cyclotomic");
    if (!IsIndet(x))
      CoCoA_THROW_ERROR("arg 2: x must be an indet", "cyclotomic");
    if (n == 1) return x-1;

    int power2 = 0;
    while (IsEven(n)) { ++power2; n /= 2; }
    const factorization<long> facs = factor(n);
    long radn = 1; for (long j: facs.myFactors()) radn *= j;
    const long deg = EulerTotient(radn);
    const vector<long> plist = facs.myFactors();
    const long long M = cyclotomic_modulus(radn);
    long p =  (M < 66000) ? M : NextPrime(1024);
    long long CRT_modulus = 1;
    vector<long> CoeffVec; // no need to reserve space
    while (M > CRT_modulus)
    {
      auto CoeffVec_modp = cyclotomic_modp(p, plist);
      if (CRT_modulus == 1)
      { // here only for the 1st iter
        CRT_modulus = p;
        swap(CoeffVec, CoeffVec_modp); /*swap is assignment*/;
        p = NextPrime(p);
        continue;
      }
      // Here every iter after the first one
      for (int i=0; i <= deg; ++i)
      {
        CoeffVec[i] = CRT2(CoeffVec[i], CRT_modulus, CoeffVec_modp[i],p);
      }
      CRT_modulus *= p;
      p = NextPrime(p);
    }

    if (power2 > 0)
    {
      for (int i=1; i < deg; i += 2)
        CoeffVec[i] = -CoeffVec[i];
    }

    // Convert result into a poly in ring of x
    const long xpower = (power2 <= 1) ? (n/radn) : ((n/radn) << (power2-1));
    const ring& P = owner(x);
    const ring& R = CoeffRing(P);
    const PPMonoidElem X = power(LPP(x), xpower);
    RingElem ans(P);
    for (int i=0; i <= deg; ++i)
    {
      if (CoeffVec[i] == 0) continue;
      PushFront(ans, RingElem(R,CoeffVec[i]), power(X,i));
    }
    return ans;
  }

  } // end of namespace anonymous


  // namespace // anonymous
  // {
  //   // 2022-12-30
  //   // This is an unused impl: it was slower than the one above;
  //   // keeping it here in case it should ever become useful.

  //   // plist is a list of distinct odd primes.
  //   RingElem cyclo_OddSqFr(const std::vector<long>& plist, const RingElem& x)
  //   {
  //     const long n = len(plist);
  //     const ring& R = owner(x);
  //     RingElem EvenProd = one(R);
  //     RingElem OddProd = one(R);
  //     SubsetIter it(n);

  //     while (!IsEnded(it))
  //     {
  //       const vector<long>& facs = *it;
  //       long prod=1; for (long j: facs) prod *= p[j]; // cannot overflow (because ...)
  //       if (IsEven(len(facs)))
  //         EvenProd *= power(x,prod)-1;
  //       else
  //         OddProd *=  power(x,prod)-1;
  //       ++it;
  //     }
  //     if (IsEven(n)) return EvenProd/OddProd;
  //     return OddProd/EvenProd;
  //   }

  // } // end of namespace anonymous


  // Compute cyclotomic polys via "truncated power series"

//   namespace // anonymous
//   {
    
//     // compute the first half of the coeffs of the cyclotomic poly (in a std::vector)
//     std::vector<unsigned long> HalfCyclo(const std::vector<long>& plist)
//     {
//       // assume plist not empty, entries are DISTINCT ODD PRIMES
//       vector<long> even;
//       vector<long> odd;
//       even.push_back(1);
//       long phi = 1;
//       for (int i=0; i < len(plist); ++i)
//       {
//         const long p = plist[i];
//         phi *= (p-1);
//         const int len_even = len(even);
//         const int len_odd = len(odd);
//         for (int j=0; j < len_even; ++j)
//           odd.push_back(p*even[j]);
//         for (int j=0; j < len_odd; ++j)
//           even.push_back(p*odd[j]);
//       }
//       if (IsOdd(len(plist)))
//         swap(odd, even);
//       // "even" corr to numerator; "odd" corr to denominator
//       sort(even.begin(), even.end());
// //NOT NEEDED    sort(odd.begin(), odd.end(), [](long a, long b) { return (a>b); });
//       const long UPB = phi/2;
//       vector<unsigned long> cyclo(UPB+1);
//       cyclo[0] = 1;
//       long d = 0; // current deg, but at most UPB
//       for (const long i: even)
//       {
//         // mult by (1+x^i)
//         long top = d;
//         if (top+i > UPB)
//           top = UPB-i; // might be negative!
//         for (long j=top; j >= 0; --j)
//           cyclo[j+i] += cyclo[j];
//         d += i; if (d > UPB) d = UPB; // d = min(d+i, UPB);
//       }
//       // Map x |-> -x
//       for (long i=1; i <= UPB; i += 2)
//         cyclo[i] = -cyclo[i];

//       for (const long i: odd)
//       {
//         // div by (1-x^i)  maybe?
//         long top = UPB - i;
//         for (long j=0; j <= top; ++j)
//           cyclo[j+i] += cyclo[j];
//       }
//       return cyclo;    
//     }
  
// }//end anon

// // Is there a smarter way?  Duplicate loop below is annoying!
//   long RecognizeCyclo_BruteForce(ConstRefRingElem f)
//   {
//     // assume f is univariate, palindromic, monic, with integer coeffs, and height < MaxLong; also deg(f) > 2
//     const long degf = deg(f);
//     vector<unsigned long> C(degf/2+1);
//     for (SparsePolyIter it=BeginIter(f); ; ++it)
//     {
//       const long d = deg(PP(it));
//       if (d < degf/2) break;
//       C[degf-d] = ConvertTo<long>(coeff(it));  // deliberate implicit conversion to "unsigned long"  !!
//     }
//     if (C[1] != 1 && C[1] != -1UL)  return 0;
//     const int mu = (C[1] == 1) ? -1 : 1;
//     const vector<long> CandidateIndexes = InvTotient(degf, InvTotientMode::SqFreePreimages);
// ///std::clog<<"Cand: " << CandidateIndexes<<std::endl;
// ///std::clog<<"mu: " << mu<<std::endl;
//     if (CandidateIndexes.empty()) return 0;  // definitely not cyclo
//     for (long k:  CandidateIndexes)
//     {
//       if (IsEven(k)) continue;
//       const factorization<long> facs = factor(k);
//       const int mu_k = (IsEven(len(facs.myFactors()))) ? 1 : -1;
// ///      std::clog<<"facs: "<<facs<<std::endl;
// ///      std::clog<<"mu_k: "<<mu_k<<std::endl;
//       if (mu_k != mu)  continue;
//       const vector<unsigned long> phi_k = HalfCyclo(facs.myFactors());
//       if (C == phi_k)  return k;
//     }
//     for (long i=1; i <= degf/2; i+=2)  C[i] = - C[i]; // YES, even though type is "unsigned long"
// ///std::clog<<"mu: " << mu<<std::endl;
//     for (long k:  CandidateIndexes)
//     {
//       if (IsEven(k)) continue;
//       const factorization<long> facs = factor(k);
//       const int mu_k = (IsEven(len(facs.myFactors()))) ? -1 : 1; // flip sign to account for "hidden" factor of 2
// ///      std::clog<<"facs: "<<facs<<std::endl;
// ///      std::clog<<"mu_k: "<<mu_k<<std::endl;
//       if (mu_k != mu)  continue;
//       const vector<unsigned long> phi_k = HalfCyclo(facs.myFactors());
//       if (C == phi_k)  return 2*k;
//     }
//     return 0;
//   }


//  } // end of namespace anonymous

  // Calls CycloPrefix then builds the actual poly.
  // TODO make a version which uses CRT to allow larger indexes?
  RingElem cyclotomic(long n, ConstRefRingElem x)
  {
    if (n < 1)
      CoCoA_THROW_ERROR("arg 1: n must be >= 1","cyclotomic");
    if (!IsIndet(x))
      CoCoA_THROW_ERROR("arg 2: x must be an indet", "cyclotomic");
    if (n == 1) return x-1;
    const long OddPart = CoprimeFactor(n,2);
    if (OddPart == 1)  return power(x,n/2) + 1;
    const bool EvenIndex = IsEven(n);
    const long EvenPart = n/OddPart; // power of 2
    factorization<long> facs = factor(OddPart);
    long PhiOddPrimes = 1;
    n = OddPart;
    long OddRadical = 1;
    for (long p: facs.myFactors())
    { OddRadical *= p; PhiOddPrimes *= (p-1); } // neither prod can overflow!
    n /= OddRadical;
    // Safeguard against coeff overflow (since we use long int internally)
    if (sizeof(unsigned long) == 4 && OddRadical >= 10163195L)
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "cyclotomic");
    if (sizeof(unsigned long) == 8 && OddRadical >= 169828113L)
      CoCoA_THROW_ERROR(ERR::ArgTooBig, "cyclotomic");
    
    const auto coeff = CycloPrefix(facs.myFactors());
    const long UPB = PhiOddPrimes/2; // exact division
    const long exponent = EvenIndex ? (n*EvenPart/2) : n;
    const PPMonoidElem X = power(LPP(x), exponent);
    { /*not used*/ const PPMonoidElem Xpwr = power(X, UPB); /* might ExpTooBig trigger error*/ if (deg(Xpwr)/UPB != deg(X)) CoCoA_THROW_ERROR(ERR::ExpTooBig, "cyclotomic"); }
    const vector<long> Xexpv = exponents(X);
    const int nvars = len(Xexpv);
    const SparsePolyRing P = owner(x);
    const ring& k = CoeffRing(P);
    RingElem coef(k);
    vector<long> expv(nvars);
    PPMonoidElem Xpwr = one(owner(X));
    RingElem ans(P);
    for (long i=0; i <= UPB; ++i)
    {
      if (coeff[i] == 0) continue;
///      const PPMonoidElem t = power(X,i);
      long C = ULong2Long(coeff[i]);
      if (EvenIndex && IsOdd(i))
        C = -C;
//      PushFront(ans, RingElem(k,C), t);
      coef = C; for(int j=0;j<nvars;++j)expv[j]=i*Xexpv[j];
      P->myPushFront(raw(ans), raw(coef), expv);
    }

    for (long i=UPB-1; i >= 0; --i)
    {
      if (coeff[i] == 0) continue;
      const PPMonoidElem t = power(X, 2*UPB-i);
      long C = ULong2Long(coeff[i]);
      if (EvenIndex && IsOdd(i))
        C = -C;
      coef = C; for(int j=0;j<nvars;++j)expv[j]=(2*UPB-i)*Xexpv[j];
///      PushFront(ans, RingElem(k,C), t);
      P->myPushFront(raw(ans), raw(coef), expv);
    }

    return ans;
  }




  //--------------------------------------------
  // Fns to check whether a poly is cyclo


  // This is a "naughty/ugly" function: it does several things at once
  //  (a) checks cyclo by matching against "prefixes" (& only this is not DoFullCheck)
  //  (b) checks that f is univariate, palindromic, below height bound H
  //  (c) checks for special cases...
  // Input poly f
  // Output:
  //  (a) if input is "obviously not cyclo" return  0 
  //  (b) special case: index is p^k or 2*p^k
  //  (c) special case: f(x) = g(x^r) for some r > 1
  unsigned long CyclotomicTest(ConstRefRingElem f, bool DoFullCheck)
  {
    const char* const FnName = "CyclotomicTest";
    constexpr unsigned long DefinitelyNotCyclo = 0;  // "impossible" index
    if (IsConstant(f))  return DefinitelyNotCyclo;
    if (!IsOne(LC(f)))  return DefinitelyNotCyclo;
    const long degf = deg(LPP(f)); // assume even, monic, not of form g(x^r), and coeff of x^(degf-1) is +-1
    SparsePolyIter it = BeginIter(f);
    ++it;  // skip LPP -- already handled it above
    if (IsEnded(it))  return DefinitelyNotCyclo; // cannot be cyclo with just a single summand!
    // Do special case deg(f) == 1:
    if (degf == 1)
    {
      if (!IsOne(PP(it)))  return DefinitelyNotCyclo;
      if (IsMinusOne(coeff(it)))  return 1;
      if (IsOne(coeff(it)))  return 2;
      return DefinitelyNotCyclo;
    }
    if (IsOdd(degf))  return DefinitelyNotCyclo;
    const PPMonoidElem x = radical(LPP(f));
    if (deg(x) != 1)
      CoCoA_THROW_ERROR(ERR::NotUnivariate, FnName);
    if (IsOne(PP(it)))
    {
      // Poly is of form x^k+nzconst
      if (IsOne(coeff(it)) && CoprimeFactor(degf,2) == 1)
        return 2*degf; // ???BUG??? overflow???
      else
        return DefinitelyNotCyclo;
    }
    
    const long r = degf - deg(PP(it));
    if (r <= 0) // cannot be univariate
      CoCoA_THROW_ERROR(ERR::NotUnivariate, FnName);
    if (degf%r != 0)  return DefinitelyNotCyclo; // not of form g(x^r)
    if (!IsOne(coeff(it)) && !IsMinusOne(coeff(it)))  return DefinitelyNotCyclo; 
    const long mu = (IsOne(coeff(it)))? -1 : 1; // minus 2nd coeff
    const long radr = radical(r);
    const long degr = degf/r; // "reduced degree"
///    std::clog<<"degr = " << degr << endl;
///    std::clog<<"r = " << r << endl;
///    std::clog<<"mu = " << mu << endl;
    // Obtain list of candidate indexes to try: in 3 stages
    // InvTotient (sqfr preimages), with correct MoebiusFn value, & must be mult of r
    vector<long> cand = InvTotient(degr, InvTotientMode::SqFreePreimages);
///    std::clog<<"len cand(init) " << len(cand) << endl;
    const auto WrongMu = [mu](long n){ return (MoebiusFn(n) != mu); };
    { auto it = std::remove_if(cand.begin(), cand.end(), WrongMu);  cand.erase(it, cand.end()); } // C++20 erase_if(...)
///    std::clog<<"len cand(mu check) " << len(cand) << endl;
    const auto NotMultOfRadr = [radr](long n){return (n%radr != 0);};
//C++20    if (radr != 1)  erase_if(cand, [radr](long n){return (n%radr != 0);});
    if (radr != 1)
    { auto it = std::remove_if(cand.begin(), cand.end(), NotMultOfRadr);  cand.erase(it, cand.end()); }
///    std::clog<<"len cand(r check) " << len(cand) << std::endl;
    if (cand.empty())  return DefinitelyNotCyclo;
    if (!DoFullCheck && len(cand) == 1)  return r*cand[0];
    const vector<int>& CycloCoeffHeightTbl = CyclotomicCoeffHeightTable();
    vector<long> CandHeightBound;
///    const long H = CycloCoeffHeightBound(degr);
    long H=0;  // overall height bound -- 0 means "to be computed"
    vector<long> C(1+degr/2);  C[0] = 1;
    long prev = 2;
    long thresh = degr/2;  while (thresh >= 64)  { thresh = 1+thresh/4; }
    ///std::clog<<"init thresh="<<thresh<<std::endl;
    // This loop handles the "upper half"
    while (!IsEnded(it))
    {
      const long d = deg(PP(it));
      if (d < degf/2)  break;
      if (radical(PP(it)) != x)
        CoCoA_THROW_ERROR(ERR::NotUnivariate, FnName);
      const long dr = (degf-d)/r;
      if (degf-d != r*dr)  return DefinitelyNotCyclo;
      long c;
      if (!IsConvertible(c, coeff(it)))  {/*std::clog<<"Not conv dr="<<dr<<std::endl;*/return DefinitelyNotCyclo;}
      if (dr < len(CycloCoeffHeightTbl))
      { if (std::abs(c) > CycloCoeffHeightTbl[dr])  {/*std::clog<<"Global coeff ht "<<c<<"  "<<CycloCoeffHeightTbl[dr]<<std::endl;*/return DefinitelyNotCyclo; } }
      else { if (H==0) { /*std::clog<<"Recompute H"<<std::endl;*/long kmax=0;for (long k: cand) kmax = std::max(kmax,k); /*std::clog<<"kmax="<<kmax<<std::endl;*/H = CycloCoeffHeightBound_index(kmax); }
        if (std::abs(c) > H)  { /*        std::clog<<" Ht "<<c<<"  "<<H<<std::endl;*/ return DefinitelyNotCyclo; } }
      C[dr] = c;
      ++it;
      // OR MAYBE:   if (dr < thresh) continue;
      if (dr >= thresh)
      {
        /*std::clog<<"Prefix check at dr="<<dr<<std::endl;*/
        vector<long> NewCand;
        for (long k: cand)
        {
          const factorization<long> facs = factor(k);
          const vector<unsigned long> prefix = CycloPrefix(facs.myFactors(), dr);
//          std::clog<<"k="<<k<<"  prefix="<<prefix<<endl;
          bool eq = true;
          for (int i=prev; i <= dr; ++i)
            if (prefix[i] != static_cast<unsigned long>(C[i])) { /*std::clog<<"MISMATCH i="<<i<<"  "<<prefix[i]<< "  "<< C[i]<<std::endl;*/eq=false; break; }
          if (eq)  NewCand.push_back(k);
        }
        if (NewCand.empty())  return DefinitelyNotCyclo;
        if (len(cand) != len(NewCand)) { swap(cand, NewCand); H = 0/*force recomputation*/; }
        if (!DoFullCheck && len(cand) == 1)  return r*cand[0];
        /*std::clog<<"Updated cand="<<cand<<std::endl;*/
        prev = dr+1;
        do {thresh *= 4;} while (thresh <= 2*dr);
        thresh = std::min(thresh, degr/2);
        //std::clog<<"new thresh="<<thresh<<std::endl;
      }
    }
    if (IsEnded(it)) return DefinitelyNotCyclo;
//    std::clog<<"PP="<<PP(it)<<"  coeff="<<coeff(it)<<endl;
    // Next loop checks that f is palindromic
    // std::clog<<"palin check   with  C[0]="<<C[0]<<std::endl;
    // Loop below checks that f is palindromic (actually a bit slow)
    long PredictedDegr = degr/2;
    while (!IsEnded(it))
    {
      while (C[--PredictedDegr] == 0); // SAFE!  PredictedDegr can never become negative since C[0] != 0
//      std::clog<<"PredictedDegr="<<PredictedDegr<<endl;
     const long d = deg(PP(it));
//     std::clog<<"obs deg "<<d<<endl;
      if (d != r*PredictedDegr)  return DefinitelyNotCyclo;
      if (d > 0 && radical(PP(it)) != x)  // because not univariate
        CoCoA_THROW_ERROR(ERR::NotUnivariate, FnName);
//      if (coeff(it) != C[PredictedDegr])  return DefinitelyNotCyclo;  // SLOWER THAN CODE BELOW
      long c;
      if (!IsConvertible(c,coeff(it)) || c != C[PredictedDegr])  return DefinitelyNotCyclo;
      ++it;
    }
//    std::clog<<"HERE RTN "<<r*cand[0]<<endl;
    return r*cand[0];
  }




  //------------------------------------------------------------------
  // Impl of CyclotomicFactors (developed from Smyth+Beukers)

  namespace // anonymous
  {
    // Input is f(x^2); output is f(x)
    RingElem sqrtx(ConstRefRingElem f)
    {
      // Assume f is univariate;  non-zero???
      if (IsConstant(f)) return f;
      const ring& P = owner(f); // assume sparse poly ring
      PPMonoidElem x = radical(LPP(f));
      RingElem ans(P);
      for (SparsePolyIter it=BeginIter(f); !IsEnded(it); ++it)
      {
        PushBack(ans, coeff(it), power(x,deg(PP(it))/2));  // assume exact division!
      }
      return ans;
    }



    void CONCAT_move(vector<RingElem>& v, vector<RingElem> v2)
    {
      // Taken from https://stackoverflow.com/questions/201718/concatenating-two-stdvectors
      v.insert(
        v.end(),
        std::make_move_iterator(v2.begin()),
        std::make_move_iterator(v2.end())
               );
    }


  
    std::vector<RingElem> SmythBeukersOpG(ConstRefRingElem f, long Xindex)
    {
      VerboseLog VERBOSE("SmythBeukersOpG");
      // assume f is univariate (non-const?)
      const ring& P = owner(f);
      const long nvars = NumIndets(P);
      vector<RingElem> ImagesNegate(nvars, zero(P));
      ImagesNegate[Xindex] = -indet(P,Xindex);
      RingHom NegateX = PolyAlgebraHom(P,P,ImagesNegate);
      RingElem fneg = NegateX(f);
      VERBOSE(85) << "Computing even factor" << std::endl;
      RingElem f2 = gcd(f, fneg);
      vector<RingElem> ans;
      if (IsConstant(f)) { VERBOSE(85) << "Trivial even factor" << std::endl; return ans; } // empty list
      vector<RingElem> ImagesSquare(nvars, zero(P));
      ImagesSquare[Xindex] = power(indet(P,Xindex),2);
      RingHom SquareX = PolyAlgebraHom(P,P, ImagesSquare);
      RingElem g1 = SquareX(f/f2);
      RingElem g2 = SquareX(f2);
      RingElem g3 = SquareX(fneg/f2);
      VERBOSE(85) << "Computing 3 gcds..." << std::endl;
      VERBOSE(89) << "GCD 1" << std::endl;
      VERBOSE(89) << "f = " << f << std::endl;
      VERBOSE(89) << "g1 = " << g1 << std::endl;
      VERBOSE(89) << "g2 = " << g2 << std::endl;
      VERBOSE(89) << "g3 = " << g3 << std::endl;
      RingElem h = gcd(f,g1);
      VERBOSE(89) << "gcd(f,g1) = " << h << std::endl;
      if (!IsConstant(h)) ans.push_back(h);
      VERBOSE(89) << "GCD 2" << std::endl;
      h = gcd(f,g2);
      VERBOSE(89) << "gcd(f,g2) = " << h << std::endl;
      if (!IsConstant(h)) ans.push_back(h);
      VERBOSE(89) << "GCD 3" << std::endl;
      h = gcd(f,g3);
      VERBOSE(89) << "gcd(f,g3) = " << h << std::endl;
      if (!IsConstant(h)) ans.push_back(h);
      VERBOSE(85) << "Finished" << std::endl;
      return ans;
    }

    std::vector<RingElem> SmythBeukersOpC(ConstRefRingElem f, long Xindex)
    {
      VerboseLog VERBOSE("SmythBeukersOpC");
      VERBOSE(85) << "Applying OpG" << std::endl;
      vector<RingElem> L = SmythBeukersOpG(f, Xindex);
      if (!L.empty())
      {
        long d=0; for (const RingElem& g: L) d += deg(g);
        if (d == deg(f)) { VERBOSE(85) << "Deg check short cut" << std:: endl; return L; }
      }
      // Copy-paste from above (redundant!!)
      const ring& P = owner(f);
      const long nvars = NumIndets(P);
      vector<RingElem> ImagesNegate(nvars, zero(P));
      ImagesNegate[Xindex] = -indet(P,Xindex);
      RingHom NegateX = PolyAlgebraHom(P,P,ImagesNegate);
      vector<RingElem> ImagesSquare(nvars, zero(P));
      ImagesSquare[Xindex] = power(indet(P,Xindex),2);
      RingHom SquareX = PolyAlgebraHom(P,P, ImagesSquare);

      VERBOSE(85) << "Checking even factor" << std::endl;
      RingElem f2 = gcd(f, NegateX(f));
      if (IsConstant(f2)) { VERBOSE(85) << "Trivial even factor" << std::endl; return L; }
      RingElem g = sqrtx(f2);
      VERBOSE(85) << "Applying OpC to even factor" << std::endl;
      vector<RingElem> Lg = SmythBeukersOpC(g,Xindex);
      for (int i=0; i < len(Lg); ++i) Lg[i] = SquareX(Lg[i]);

//    vector<RingElem> V;
      for (const RingElem& ff: L) { CONCAT_move(Lg, SmythBeukersOpC(ff, Xindex)); }
///    Lg = CoprimeFactorBasis(Lg);
      VERBOSE(85) << "Computing CoprimeFactorBasis" << std::endl;
      CoprimeFactorBasis_RingElem CFB;  CFB.myAddInfo(Lg);
      Lg = FactorBase(CFB);
      return Lg;
    }

  } // end of namespace anonymous



  // BUG???  Put rest of f into RemainingFactor???
  factorization<RingElem> CycloFactors_SmythBeukers(ConstRefRingElem f)
  {
    VerboseLog VERBOSE("SmythBeukers");
    if (IsZero(f)) CoCoA_THROW_ERROR(ERR::NotNonZero, "CycloFactors");
    if (IsConstant(f))  return factorization<RingElem>(f);
    const long x = UnivariateIndetIndex(f);
    if (x < 0) CoCoA_THROW_ERROR("Arg must be univariate", "CycloFactors");
    factorization<RingElem> ans(one(owner(f)));
    if (IsConstant(f)) return ans;// vector<RingElem>();
    VERBOSE(85) << "Extract palindromic factor" << std::endl;
    const RingElem pal = PalindromicFactor(f);
    if (IsConstant(pal)) return ans; //vector<RingElem>(); // no cyclo factors
    VERBOSE(85) << "Compute SquareFree factors" << std::endl;
    const auto facs = SqFreeFactor(pal);
    VERBOSE(85) << "Num factors: " << len(facs.myFactors()) << std::endl;

    for (int i=0; i < len(facs.myFactors()); ++i)
    {
      const long m = facs.myMultiplicities()[i];
      VERBOSE(85) << "Doing factor " << i << " (mult=" << m << ", deg="<< deg(facs.myFactors()[i]) << ")" << std::endl;
      for (const RingElem& g: SmythBeukersOpC(facs.myFactors()[i],x))
        ans.myAppend(g, m);
    }
    return ans;
  }



  namespace // anonymous
  {


// class EvalPoly
// {
// public:
//   EvalPoly(ConstRefRingElem f): myF(f) {}
//   BigInt operator()(long n) const;
// private:  // data members
//   RingElem myF;
// };

// BigInt EvalPoly::operator()(long n) const
// {
//   if (IsZero(myF)) return BigInt();
//   if (IsConstant(myF)) return ConvertTo<BigInt>(ConstantCoeff(myF));
//   return EvalAt(myF,n);  
// }

class EvalPolyAtRAT
{
public:
  EvalPolyAtRAT(ConstRefRingElem f): myF(f) {}
  BigInt operator()(long n, long d) const;
private:  // data members
  RingElem myF;
};

    BigInt EvalPolyAtRAT::operator()(long n, long d) const
{
  if (IsZero(myF)) return BigInt();
  if (IsConstant(myF)) return ConvertTo<BigInt>(ConstantCoeff(myF));
  return EvalAt(myF, n,d);
}


    // "next biggest rational"
    void  GotoNextEvalPt(long& n, long& d)
    {
      // assume n >= 2 and 1 <= d < n
//      constexpr long UPB = 3;
      if (d == n-1) { ++n; d = 1; return; }
      do { ++d; } while (gcd(n,d) != 1);
    }



    
//  CandidateCycloIndexes is a list of integers >= 2 being indexes of
//  possible cyclotomic factors -- see also CycloIndexesUptoDeg (below).
  std::vector<long> FindCycloFactor(vector<long> CandidateCycloIndexes, std::function<BigInt /**/(long n, long d)> EvalF)
  {
    // Assume f is non-constant with integer coeffs (ideally squarefree, content=1, and f palindromic)
    VerboseLog VERBOSE("FindCycloFactor");

    const ring ZZx = NewPolyRing(RingZZ(), symbols("x")); // used only to evaluate a cyclo poly
    const RingElem& x = indet(ZZx,0);                    //


    VERBOSE(80) <<"NUM INIT CAndidates: " << len(CandidateCycloIndexes) << std::endl;
    int TargetNumIters = 2;
    long EvalPtNumer = 1;  long EvalPtDenom = 0;
    for (int i=0; i < TargetNumIters; ++i)
    {
      BigInt Valf; // will be set in loop below
///      RingHom evaluate = IdentityHom(P); // will be set in loop below
      int NumConsectiveZeroes = 0;
///      clog<<"DO LOOP START"<<std::endl;
      do
      {
        if (NumConsectiveZeroes < 1) { GotoNextEvalPt(EvalPtNumer, EvalPtDenom); }
        else { long skip = EvalPtNumer + RandomLong(1, EvalPtNumer); for (long i=0; i < skip; ++i)  GotoNextEvalPt(EvalPtNumer, EvalPtDenom); NumConsectiveZeroes = 0; }

///        clog<<" Try EvalPt = " << EvalPtNumer << "/" << EvalPtDenom << endl;
///        evaluate = EvalHom(P, vector<RingElem>(nvars, RingElem(CoeffRing(P), EvalPt)));
///        Valf = ConvertTo<BigInt>(evaluate(f));
        Valf = EvalF(EvalPtNumer, EvalPtDenom);
///        std::clog<<"Valf = " << Valf << endl;
        if (IsZero(Valf)) std::clog<<"Z";
        if(IsZero(Valf)) ++NumConsectiveZeroes; else NumConsectiveZeroes = 0;
      } while (IsZero(Valf));
///      clog<<"DO LOOP END"<<std::endl;
      VERBOSE(80) <<"Chosen EvalPt="<<EvalPtNumer << "/" << EvalPtDenom <<std::endl;
      BigInt Valf_reduced = CoprimeFactor(Valf, EvalPtNumer*(EvalPtNumer*EvalPtNumer-EvalPtDenom*EvalPtDenom)); // BUG?????  to avoid overflow need EvalPt^3 < MaxLong (or MaxULong)

      vector<long> ReducedCandidateList;
      for (long k:  CandidateCycloIndexes)
      {
        if (EvalPtNumer == 2 && EvalPtDenom == 1 && k == 6) { ReducedCandidateList.push_back(k); continue; } // exception in Zsygmondi's Thm
//???        /*const*/ BigInt g = PowerMod(EvalPtNumer, k, Valf_reduced) - PowerMod(EvalPtDenom,k,Valf_reduced);
        /*const*/ BigInt g = power(EvalPtNumer, k) - power(EvalPtDenom,k);
        g = gcd(g, Valf_reduced);
        if (IsOne(g)) { /*clog << '*';*/ continue; }
        Valf_reduced = CoprimeFactor(Valf_reduced, g);
        const BigInt ValCyclo = EvalAt(cyclotomic(k,x),EvalPtNumer,EvalPtDenom);///ConvertTo<BigInt>(evaluate(cyclotomic(k,x)));
        if (IsDivisible(Valf,ValCyclo))
          ReducedCandidateList.push_back(k);
        if (IsOne(Valf_reduced)) { if(EvalPtNumer==2&&EvalPtDenom==1&&k<6)ReducedCandidateList.push_back(6);break; }
      }
      if (len(ReducedCandidateList) < len(CandidateCycloIndexes))
      {
        TargetNumIters = i+3; // if no indexes removed in 3 consecutive iters, we accept the list as "probably correct" (i.e. only very few false positives)
///        clog << "iter is " << i << "   New iter target is " << TargetNumIters << endl;
        swap(CandidateCycloIndexes, ReducedCandidateList); // really assignment
        VERBOSE(80)<<"CandidateCycloIndexes="<<CandidateCycloIndexes<<std::endl;
if (CandidateCycloIndexes.empty()) break;
      }
    }
    return CandidateCycloIndexes;
  }



// //  std::vector<long> FindCycloFactor(long d, BigInt EvalF(long))
//   std::vector<long> FindCycloFactor(vector<long> CandidateCycloIndexes, std::function<BigInt /**/(long)> EvalF)
//   {
//     VerboseLog VERBOSE("FindCycloFactor");
//     // Assume f is non-constant with integer coeffs (ideally squarefree, content=1, and f palindromic)
//     const ring ZZx = NewPolyRing(RingZZ(), symbols("x")); // used only to evaluate a cyclo poly
//     const RingElem& x = indet(ZZx,0);                    //

// //    const long degf = deg(f);
// //    long UPB = InvTotientBoundUpto_ulong(degf);
//     // vector<long> PossibleCycloIndexes;  PossibleCycloIndexes.reserve(2*degf); //??
//     // // ignore cyclos with deg < 10
//     // for (long k=3; k <= UPB; ++k)
//     // {
//     //   if (EulerTotient(k) > degf) continue;
//     //   PossibleCycloIndexes.push_back(k);
//     // }

// ///    const ring& P = owner(f);
// ///    const long nvars = NumIndets(P);
// ///    RingElem x = indet(P, UnivariateIndetIndex(f));

// ///std::clog<<"INIT CAndidates: " << CandidateCycloIndexes << std::endl;
//     VERBOSE(80) << "Num candidate indexes: " << len(CandidateCycloIndexes) << std::endl;
//     int TargetNumIters = 2;
//     long EvalPt = 1;
//     for (int i=0; i < TargetNumIters; ++i)
//     {
//       BigInt Valf; // will be set in loop below
// ///      RingHom evaluate = IdentityHom(P); // will be set in loop below
//       int NumConsectiveZeroes = 0;
//       do
//       {
//         if (NumConsectiveZeroes < 1) {
//           if (true)        EvalPt=NextPrime(EvalPt); else   ++EvalPt;}
//         else { EvalPt += EvalPt/2+RandomLong(1, EvalPt/2); // random jump to betw (3/2)* and 2*
//           NumConsectiveZeroes = 0; }
//         VERBOSE(85) <<" Testing EvalPt = " << EvalPt << std::endl;
// ///        evaluate = EvalHom(P, vector<RingElem>(nvars, RingElem(CoeffRing(P), EvalPt)));
// ///        Valf = ConvertTo<BigInt>(evaluate(f));
//         Valf = abs(EvalF(EvalPt));
// //        std::clog<<"Valf = " << Valf << endl;
//         if (IsZero(Valf)) std::clog<<"Z";
//         if(IsZero(Valf)) ++NumConsectiveZeroes; else NumConsectiveZeroes = 0;
//       } while (IsZero(Valf));
//       VERBOSE(80) <<"Chosen EvalPt="<<EvalPt<<std::endl;
//       BigInt Valf_reduced = CoprimeFactor(Valf, EvalPt*(EvalPt*EvalPt-1)); // BUG?????  to avoid overflow need EvalPt^3 < MaxLong (or MaxULong)

//       vector<long> ReducedCandidateList;
//       for (long k : CandidateCycloIndexes)
//       {
//         if (EvalPt == 2 && k == 6) { ReducedCandidateList.push_back(k); continue; } // exception in Zsygmondi's Thm
//         if (IsOne(Valf_reduced)) break; // cannot be any more cyclo factors
//         /*const*/ BigInt g = PowerMod(EvalPt, k, Valf_reduced) - 1;
//         g = gcd(g, Valf_reduced);
//         if (IsOne(g)) { /*clog << '*';*/ continue; }
//         Valf_reduced = CoprimeFactor(Valf_reduced, g);
//         const BigInt ValCyclo = EvalAt(cyclotomic(k,x),EvalPt);///ConvertTo<BigInt>(evaluate(cyclotomic(k,x)));
//         if (IsDivisible(Valf,ValCyclo))
//           ReducedCandidateList.push_back(k);
//         if (IsOne(Valf_reduced)) { if(EvalPt==2&&k<6)ReducedCandidateList.push_back(6);break; }
//       }
//       if (len(ReducedCandidateList) < len(CandidateCycloIndexes))
//       {
//         TargetNumIters = i+3; // if no indexes removed in 3 consecutive iters, we accept the list as "probably correct" (i.e. only very few false positives)
// ///        clog << "iter is " << i << "   New iter target is " << TargetNumIters << endl;
//         swap(CandidateCycloIndexes, ReducedCandidateList); // really assignment
//         VERBOSE(80) << "new CandidateCycloIndexes="<<CandidateCycloIndexes<<std::endl;
//       }
//     }
//     return CandidateCycloIndexes;
//   }



std::vector<long> CycloIndexesUptoDeg(long d)
{
    vector<long> CandidateIndexes;
    const long UPB = InvTotientBoundUpto_ulong(d);
///    clog << "CycloIndexesUptoDeg:  UPB = " << UPB << endl;
    for (long k=3; k <= UPB; ++k)
    {
      if (IsOdd(k) && k > UPB/2) continue;
      const long phi = EulerTotient(k);
      if (/*phi > 3 &&*/ phi <= d)
        CandidateIndexes.push_back(k);
    }
    return CandidateIndexes;
}






  } // end of namespace anonymous


std::vector<long> CyclotomicFactorIndexes(ConstRefRingElem f)
{
  return FindCycloFactor(CycloIndexesUptoDeg(deg(f)), EvalPolyAtRAT(f));
}

  factorization<RingElem> CyclotomicFactors(ConstRefRingElem f)
  {
    VerboseLog VERBOSE("CyclotomicFactors");
    // assume f is univariate over ZZ/QQ
    VERBOSE(80) << "PALINDROMIC: " << deg(PalindromicFactor(f)) << std::endl;
    VERBOSE(80) << "Running Smyth+Beukers" << std::endl;
    const factorization<RingElem> facs = CycloFactors_SmythBeukers(PalindromicFactor(f));
//    const factorization<RingElem> facs = SqFreeFactor(PalindromicFactor(f));
    VERBOSE(80) << "Finished Smyth+Beukers" << std::endl;
///    cout << "facs = " << facs << endl;
    const vector<RingElem>& F = facs.myFactors();
    const vector<long>& m = facs.myMultiplicities();
    const ring& P = owner(f);
    const RingElem& x = indet(P, UnivariateIndetIndex(f));
    factorization<RingElem> ans(one(P)); // remaning factor
    for (int i=0; i < len(F); ++i)
    {
      VERBOSE(85) << "LOOP("<<i<<"): doing " << LPP(F[i]) << std::endl;
      if (CyclotomicIndex(F[i]) != 0)
      { VERBOSE(85) << "RECOGNIZED: "<<CyclotomicIndex(F[i])<<std::endl; ans.myAppend(F[i], m[i]); continue; }
      const vector<long> index = CyclotomicFactorIndexes(F[i]);
      VERBOSE(85) << "LOOP("<<i<<")  CandidateIndexes=" << index << std::endl;
      long D = 0; for (long k: index) D += EulerTotient(k);
      bool MustCheck = (D != deg(F[i]));
      for (long k: index)
      {
        RingElem phi = cyclotomic(k,x);
        if (MustCheck && !IsDivisible(f,phi)) continue;
        ans.myAppend(phi, m[i]);
      }
    }
    return ans;
  }










  namespace ////OBSOLETE_CODE
  {

    unsigned long IsCyclotomicGraeffe(ConstRefRingElem f);

    /**
     * Checks if f has a square root and then if that
     * square root is cyclotomic.  Returns true iff
     * both conditions are met.
     * For private use.  Hence, no arg checking.
     */
    bool HasCyclotomicSqrt(ConstRefRingElem f)
    {
      VerboseLog VERBOSE("HasCyclotomicSqrt");
      factorization<RingElem> fac = SqFreeFactor(f);
      VERBOSE(80) << "fac = " << fac << std::endl;
      if (!IsOne(fac.myRemainingFactor()))
        return false;
      RingElem sqrt = one(owner(f));

      const std::vector<RingElem> &factors = fac.myFactors();
      const std::vector<long> &mults = fac.myMultiplicities();
      for (size_t i = 0; i < factors.size(); ++i)
      {
        CheckForInterrupt("HasCyclotomicSqrt");

        const long mult = mults[i];

        if (IsOdd(mult))
          return false;

        sqrt *= power(factors[i], mult / 2);
      }
      VERBOSE(50) << "sqrt = " << sqrt << std::endl;

      return IsCyclotomicGraeffe(sqrt);
    }


    /**
     * Graeffe approach (see Bradford+Davenport, Beukers+Smyth)
     * OBSOLESCENT: NOT UPDATED AND NOT SUPPORTED!
     */
    unsigned long IsCyclotomicGraeffe(ConstRefRingElem f)
    {
      const char* const FnName = "IsCyclotomicGraeffe";
      VerboseLog VERBOSE(FnName);
      const ring &Px = owner(f);
      if (!IsSparsePolyRing(Px))
        CoCoA_THROW_ERROR(ERR::NotSparsePolyRing, FnName);
      const ring &P = CoeffRing(Px);
      const long IndetIndex = UnivariateIndetIndex(f);
      if (IndetIndex < 0)
        CoCoA_THROW_ERROR(ERR::NotUnivariate, FnName);
      if (!IsZZ(P) && !IsQQ(P))
        CoCoA_THROW_ERROR(ERR::BadRing, FnName);
      const RingElem &x = indet(Px, IndetIndex);

      // Cyclotomic polynomials have LC = 1 and CC = +-1
      if (!IsOne(LC(f)))
        return 0;
      if (!IsOne(abs(ConstantCoeff(f))))
        return 0;

      const long degf = deg(f);
      // deg = 1 -> x +- 1 is cyclotomic
      if (degf == 1)
        return ConstantCoeff(f) > 0 ? 2 : 1;
      // Constant polynomials and ones with odd degree
      // cannot be cyclotomic (except deg = 1; see above!)
      if (degf < 1 || IsOdd(degf))
        return 0;
      // Cyclotomic polynomials are always irreducible
      if (!IsIrred(f))
        return 0;

      const RingElem f1 = graeffe(f);
      VERBOSE(80) << "graeffe = " << f1 << std::endl;

      if (f1 == f)
        return 1;

      const RingElem f2 = PolyAlgebraHom(Px, Px, std::vector<RingElem>(NumIndets(Px), -x))(f);
      VERBOSE(80) << "f2 = " << f2 << std::endl;
      if (f1 == f2 && IsCyclotomicGraeffe(f2))
        return 1;
      VERBOSE(40) << "SqrtTest" << std::endl;
      if (HasCyclotomicSqrt(f1))
        return 1;

      return 0;
    }

  } // end of namespace OBSOLETE_CODE

  
  // -------------------------------------------------------
  // CyclotomicIndex

  namespace  // anonymous
  {

    // Gets 2nd highest term in f -- caller wants coeff and deg
    RingElem SecondMonomial(ConstRefRingElem f)
    {
      // Assumes NumTerms(f) >= 2
      CoCoA_ASSERT(NumTerms(f) >= 2);

      SparsePolyIter it = BeginIter(f);
      ++it;
      return monomial(owner(f), coeff(it), PP(it));
    }


    unsigned long CheckDegAndSecondCoeff(const long snd, const long degf,
                                         const unsigned long n)
    {
      if (snd == -MoebiusFn(n) && EulerTotient(n) == degf) return n;
      return 0;
    }


    RingElem RootX(RingElem f, long n) /*noexcept*/
    {
      // ASSUMES input is non constant, univariate.
      // ASSUMES input is of form g(x^n)
      CoCoA_ASSERT(n > 1);
      const ring& P = owner(f);
      CoCoA_ASSERT(IsPolyRing(P));
      CoCoA_ASSERT(!IsConstant(f));
      CoCoA_ASSERT(UnivariateIndetIndex(f) >= 0);

      PPMonoidElem x = radical(LPP(f));
      RingElem ret = zero(P);
      for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
      {
        const long d = deg(PP(it));
        ret += monomial(P, coeff(it), power(x,d/n)); // division ASSUMED exact!
      }
      return ret;
    }
  }


//     /**
//      * Inverse Phi method
//      */
//     unsigned long CyclotomicIndex_unchecked(ConstRefRingElem f)
//     {
//       const char* const FnName = "CyclotomicIndex_unchecked";
//       VerboseLog VERBOSE(FnName);
//       const ring &Px = owner(f);
//       if (!IsSparsePolyRing(Px))
//         CoCoA_THROW_ERROR(ERR::NotSparsePolyRing, FnName);
//       const ring &P = CoeffRing(Px);

//       // Dispose of trivial case
//       if (IsConstant(f)) return 0;

//       // Cyclotomic polynomials have LC = 1
//       if (!IsOne(LC(f))) return 0;

//       const long IndetIndex = UnivariateIndetIndex(f);
//       if (IndetIndex < 0)
//         CoCoA_THROW_ERROR(ERR::NotUnivariate, FnName);
//       if (!IsZZ(P) && !IsQQ(P)) // HINT: could allow other coeff rings if all coeffs are integer (???)
//         CoCoA_THROW_ERROR(ERR::BadRing, FnName);


//       // Handle the 2 special cases of deg=1:
//       const long degf = deg(f);
//       if (!IsOne(ConstantCoeff(f)))
//         return (degf == 1 && IsMinusOne(ConstantCoeff(f))); // !! auto convert true->1, false->0

//       // deg = 1 -> Only possibility is now Phi_2
//       if (degf == 1) return 2;

//       // Polynomials with odd degree cannot be cyclotomic
//       // (except deg = 1; see above!)
//       if (IsOdd(degf)) return 0;

//       // Second coefficient equal to -Moebius function (of radical of index)
//       long snd;
//       RingElem sndMonomial = SecondMonomial(f);
//       if (!IsConvertible(snd, LC(sndMonomial))) return 0; // coeff obviously not 1 or -1 or 0
//       if (uabs(snd) > 1) return 0;

//       // deg > 1 -> Always palindromic
//       if (!IsPalindromic(f)) return 0;

//       // Check coefficient height (for "small" degrees)
//       const RingElem H = CoeffHeight(f);
//       if (degf < 11612160  &&  H > CycloCoeffHeightBound(degf)) return 0;

//       // Cyclotomic polynomials are always irreducible
//       // But we do not check this as it is too costly.

//       // Now check if poly is of form g(x^d); and d must be deg difference between first 2 terms.
//       const long degDiff = degf - deg(sndMonomial);
//       if (degDiff != 1) // Non-square-free index
//       {
//         if (degf%degDiff != 0) return 0;
//         long lastDeg = degf + degDiff;
//         for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
//         {
//           const long d = deg(PP(it));
//           if ((lastDeg - d) % degDiff != 0) return 0;
//           lastDeg = d;
//         }

//         // OK Poly is of form g(x^d): get cyclo index of g
//         const unsigned long index = CyclotomicIndex_unchecked(RootX(f, degDiff));

//         if (index == 0) return 0;
//         // Next lines check that every prime dividing degDiff also divides index; if not, return 0 [not cyclo]
//         long r = degDiff;
//         for (long g=gcd(r,index); g > 1; g = gcd(r,index))
//         { r/=g; }
//         if (r != 1) return 0;
//         return degDiff * index;
//       } // end of case f(x) = g(x^d)

//       // From here onwards we know that if is cyclo then its index is sqfree.
//       const vector<long> InvTot = InvTotient(degf,  InvTotientMode::SqFreePreimages);
//       if (InvTot.empty())  return 0;  // e.g. 54 has no sqfr preimages

//       // See https://arxiv.org/pdf/1611.06783.pdf
//       BigInt evalAtPos1;
//       BigInt evalAtNeg1;
//       for (SparsePolyIter it = BeginIter(f); !IsEnded(it); ++it)
//       {
//         const long d = deg(PP(it));
//         const BigInt c = ConvertTo<BigInt>(coeff(it));
//         evalAtPos1 += c;
//         if (IsEven(d)) evalAtNeg1 += c; else evalAtNeg1 -= c;
//       }

//       if (!IsOne(evalAtPos1))
//       {
//         // We already know the cyc index to be square-free at this point,
//         // so it can only be the prime index corresponding to evalAtPos1
//         if (evalAtPos1 != 1+degf || !IsOne(H)) return 0;
//         if (!IsPrime(1+degf)) return 0; // ??? BUG: OVERFLOW ???
//         return 1+degf; // we are certain this is correct
// //        return CheckDegAndSecondCoeff(snd, degf, evalAtPos1);
//       }

//       if (!IsOne(evalAtNeg1))
//       {
//         // We already know the cyc index to be square-free at this point,
//         // so it can only be 2*evalAtNeg1 (2*p^e for some e in general)
//         if (evalAtNeg1 != 1+degf || !IsOne(H)) return 0;
//         if (!IsPrime(1+degf)) return 0; // ??? BUG: OVERFLOW ???
//         return 2*(1+degf); // we are certain this is correct:  ??? BUG:  OVERFLOW ???
//       }

// //??? if InvTot contains only few values, maybe use brute force???
      
//       // Now we evaluate at x|->2 and check that the result is plausible
// //      const std::vector<RingElem> images(NumIndets(Px), 2*one(P));
// //      const BigInt ValAt2 = ConvertTo<BigInt>(EvalHom(Px, images)(f));
// const double t0=CpuTime();
// const BigInt ValAt2 = EvalAt2_dodgy(f);//EvalAt(f,2);
// std::clog<<"EvalAt2 time: "<<CpuTime()-t0<<std::endl;
//       if (sign(ValAt2) != 1 /*|| IsEven(ValAt2)*/) return 0;
//       const long Log2Eval = FloorLog2(ValAt2);
//       if (Log2Eval > degf+1 || Log2Eval < degf-1) return 0;

//       const unsigned long MaxExp = InvTotientBound_ulong(degf); // InvTotientBound_long ???
//       if (MaxExp == 0) return 0; // not cyclo
//       if (MaxExp == 1)
//         CoCoA_THROW_ERROR(ERR::ArgTooBig, "IsCyclotomicIter");
//       // const BigInt bound = InvTotientBound(degf); // InvTotientBound_long ???
//       // unsigned long MaxExp;
//       // if (!IsConvertible(MaxExp, bound))
//       //   CoCoA_THROW_ERROR(ERR::ArgTooBig, "IsCyclotomicIter");

// //???      VERBOSE(20) << "Exp bound = " << MaxExp << std::endl;

        
//       // IMPORTANT: must handle phi_3 = x^2+x+1 and phi_6 = x^2-x+1 as special cases!
//       if (degf==2) { /*assume palindromic*/; if (snd==1) return 3; else return 6;}

//       unsigned long exp = degf; /// NOT REALLY NECESSARY:  if (snd==-1) ++exp;
//       BigInt pwr2mod = power(2,2*degf);
//       while (exp < MaxExp)
//       {
//         exp += degf;
//         mod(pwr2mod, pwr2mod, ValAt2);
//         if (IsPowerOf2(pwr2mod))
//         {
//           const long log_pwr2mod = FloorLog2(pwr2mod);
//           const long index = exp-log_pwr2mod; // candidate index!
// std::clog<<"TIME (EvalAt2 + checking PowerMod values): "<<CpuTime()-t0<<std::endl;
//           return CheckDegAndSecondCoeff(snd, degf, index);
//         }
//         pwr2mod *= power(2,degf);
// //        mpz_mul_2exp(mpzref(pwr2mod), mpzref(pwr2mod), degf);   // mpz_mul_2exp is not appreciably faster ?!?
//       }
//       // OLD IMPL
//       // while (exp < MaxExp)
//       // {
//       //   exp += degf; if (exp > MaxExp) exp = MaxExp;
//       //   const BigInt rem = PowerMod(2,exp,ValAt2);
//       //   if (IsPowerOf2(rem))
//       //   {
//       //     const long log_rem = FloorLog2(rem);
//       //     const long index = exp-log_rem; // candidate index!
//       //     return CheckDegAndSecondCoeff(snd, degf, index);
//       //   }
//       // }
//       return 0;
//     }


  // If f is monic cyclo with index n, return n; o/w return 0.
  // Uses CyclotomicIndex_unchecked (above) to do most of the work
    unsigned long CyclotomicIndex(ConstRefRingElem f)
    {
      return CyclotomicTest(f, false/*no full check*/);
      // const unsigned long CandidateIndex = CyclotomicIndex_unchecked(f);
      // if (CandidateIndex == 0) return 0;
      // // Nothing more to do for cases where _unchecked actually gives a guaranteed answer:
      // if (CandidateIndex <= 2 || IsPrime(CandidateIndex)) return CandidateIndex;
      // if (CandidateIndex%4 == 2 && IsPrime(CandidateIndex/2)) return CandidateIndex;
      // const ring& P = owner(f);
      // const RingElem& x = indet(P, UnivariateIndetIndex(f)); // we are sure that f is univariate
      // if (f == cyclotomic(CandidateIndex, x))
      //   return CandidateIndex;
      // return 0;
    }

    unsigned long CyclotomicTest(ConstRefRingElem f)
    {
      return CyclotomicTest(f, true/*do full check*/);
    }






}  // end of namespace CoCoA


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-cyclotomic.C,v 1.20 2023/06/28 12:10:06 bigatti Exp $
// $Log: SparsePolyOps-cyclotomic.C,v $
// Revision 1.20  2023/06/28 12:10:06  bigatti
// Summary: changed OVERFLOW into OVERFLOWz for MacOS problem
//
// Revision 1.19  2023/06/23 13:54:52  abbott
// Summary: Added check for exp overflow in cyclotomic
//
// Revision 1.18  2023/06/16 07:56:38  abbott
// Summary: Updated CyclotomicTest(2args); still needs work
//
// Revision 1.17  2023/06/08 17:11:13  abbott
// Summary: Fixed silly bug (not on redmine)
//
// Revision 1.16  2023/06/06 20:18:45  abbott
// Summary: Now throws if input is not univariate
//
// Revision 1.15  2023/06/06 19:20:59  abbott
// Summary: Now have CyclotomicTest and CyclotomicIndex
//
// Revision 1.14  2023/05/16 07:19:06  abbott
// Summary: Some tidying
//
// Revision 1.13  2023/05/11 19:27:33  abbott
// Summary: Added CycloFactorIndexes; some other temp changes
//
// Revision 1.12  2023/04/12 20:13:55  abbott
// Summary: Added lots of verbosity stuff (mostly for debugging) Remove later?
//
// Revision 1.11  2023/03/28 14:44:54  abbott
// Summary: Fixed minor (embarrassing) bug
//
// Revision 1.10  2023/03/15 21:32:36  abbott
// Summary: Renamed CycloFactors to CyclotomicFactors; also refined semantics
//
// Revision 1.9  2023/03/13 19:51:13  abbott
// Summary: Corrected PushFront into PushBack (oops)
//
// Revision 1.8  2023/03/09 22:36:31  abbott
// Summary: Changed order initial easy tests -- hope it is faster
//
// Revision 1.7  2023/02/23 20:50:47  abbott
// Summary: Added new fn CycloFactors (poor name)
//
// Revision 1.6  2023/02/14 20:01:40  abbott
// Summary: Removed table-lookup; minor adjustment after changing InvTotientBound
//
// Revision 1.5  2023/02/02 19:00:55  abbott
// Summary: Tidying up (& some improvements)
//
// Revision 1.4  2023/01/31 21:32:33  abbott
// Summary: Tidying up
//
// Revision 1.3  2023/01/31 12:31:34  abbott
// Summary: Moved CyclotomicIndex (& unchecked version) here from SparsePolOps-Cyclotomicity
//
// Revision 1.2  2023/01/01 11:58:02  abbott
// Summary: Commented out old code which required a dependency not needed by the main code
//
// Revision 1.1  2023/01/01 11:33:44  abbott
// Summary: Improved impl of cyclotomic; also cyclotomic moved to new file
//
//

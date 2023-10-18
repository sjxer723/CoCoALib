//   Copyright (c)  2022  John Abbott and Anna M. Bigatti
//   Original author:  Nico Mexis  (2022)

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


#include "CoCoA/SparsePolyOps-LRSDegeneracy.H"
#include "CoCoA/CanonicalHom.H"
#include "CoCoA/factor.H"
#include "CoCoA/factorization.H"
#include "CoCoA/interrupt.H"
#include "CoCoA/NumTheory-gcd.H"
#include "CoCoA/NumTheory-misc.H"
#include "CoCoA/NumTheory-modular.H"
#include "CoCoA/NumTheory-PrimeModSeq.H"
#include "CoCoA/QuotientRing.H"
#include "CoCoA/ring.H"
#include "CoCoA/RingQQ.H"
#include "CoCoA/RingZZ.H"
#include "CoCoA/SparsePolyOps-cyclotomic.H"
#include "CoCoA/SparsePolyOps-Graeffe.H"
#include "CoCoA/SparsePolyOps-resultant.H"
#include "CoCoA/SparsePolyOps-RingElem.H"
#include "CoCoA/SparsePolyRing.H"
#include "CoCoA/verbose.H"
#include "CoCoA/VerificationLevel.H"

#include <vector>
using std::vector;
#include <algorithm>
using std::min;

namespace CoCoA
{

  namespace // anonymous
  {

    // Check input poly: throw if not univariate over char 0
    // Returns the radical (divided by x if possible) -- if this has deg < 2 then obv not LRSdegen
    // Error if f not univariate: return radical(f), and divide by x if possible.
    RingElem LRSDegenerate_check_and_preprocess(ConstRefRingElem f, const ErrorContext& OrigContext)
    {
      const ring &Px = owner(f);
      if (!IsSparsePolyRing(Px))
        CoCoA_THROW_ERROR_WITH_CONTEXT(ERR::NotSparsePolyRing, OrigContext);
      if (!IsZero(characteristic(Px)))
        CoCoA_THROW_ERROR_WITH_CONTEXT("Characteristic must be 0", OrigContext);
      if (IsConstant(f))
        CoCoA_THROW_ERROR_WITH_CONTEXT("Polynomial must be non-constant", OrigContext);
      const long IndetIndex = UnivariateIndetIndex(f);
      if (IndetIndex < 0)
        CoCoA_THROW_ERROR_WITH_CONTEXT(ERR::NotUnivariate, OrigContext);

      // Divide out factors of x:
      RingElem f_reduced = prim(f);
      if (IsZero(ConstantCoeff(f)))
        f_reduced /= gcd(f, IndetPower(Px, IndetIndex, deg(f)));

      // Extract square-free part
      f_reduced /= gcd(f_reduced, deriv(f_reduced, IndetIndex));

      return f_reduced;
    }

  } // end of namespace anonymous


  namespace OBSOLETE
  {

    /**
     * Algorithm 1 (naive factor-approach)
     * OBSOLESCENT: NOT UPDATED AND NOT SUPPORTED!
     */
    unsigned long IsLRSDegenerateNaive(RingElem f)
    {
      const char* const FnName = "IsLRSDegenerateNaive";
      VerboseLog VERBOSE(FnName);
      const ErrorContext context = CoCoA_ERROR_CONTEXT;
      f = LRSDegenerate_check_and_preprocess(f, context);
      if (deg(f) < 2) return 0;
//        CoCoA_THROW_ERROR_WITH_CONTEXT("Polynomial does not have two distinct non-zero roots.", OrigContext);
      const ring& Px = owner(f);
      const ring& ZZ = RingZZ(); // just an abbr
      const SparsePolyRing ZZxy = NewPolyRing(ZZ, symbols("x,y"));
      vector<RingElem> images(NumIndets(Px), indet(ZZxy,0));
      const RingHom toZZxy = PolyRingHom(Px, ZZxy, CanonicalHom(CoeffRing(Px),ZZ), images);
      f = toZZxy(prim(f));

      const RingElem& x = indet(ZZxy,0);
      vector<RingElem> NegImages{-x, zero(ZZxy)};
      RingHom NegX = PolyAlgebraHom(ZZxy, ZZxy, NegImages);
      if (!IsConstant(gcd(f, NegX(f))))
        return 2;

      const RingElem& y = indet(ZZxy,1);

      // g is a polynomial whose roots are the quotients of the roots of f
      RingHom x2xy = PolyAlgebraHom(ZZxy, ZZxy, std::vector<RingElem>{x*y, zero(ZZxy)});
      RingElem g = resultant(f, x2xy(f), 0); // 0 is index of "x"
      g /= power(y-1, deg(f));

      const auto facs = CyclotomicFactors(f);
      VERBOSE(50) << "facs = " << facs << std::endl;
      return !(facs.myFactors().empty());
      // f is LRS-degenerate if g has a cyclotomic factor
      // other than x-1, which would indicate the quotient
      // of a root and itself
    }


    /**
     * Algorithm 2 (Graeffe iteration)
     * OBSOLESCENT: NOT UPDATED AND NOT SUPPORTED!
     */
    unsigned long LRSDegeneracyOrder_Iter(RingElem f)
    {
      const char* const FnName = "LRSDegeneracyOrder_Iter";
      VerboseLog VERBOSE(FnName);
      const ErrorContext context = CoCoA_ERROR_CONTEXT;
      f = LRSDegenerate_check_and_preprocess(f, context);
      if (deg(f) < 2) return 0;  // or ERROR????
      // if (!IsSparsePolyRing(Px))
      //   CoCoA_THROW_ERROR(ERR::NotSparsePolyRing, FnName);
      // const ring &P = CoeffRing(Px);
      // const long IndetIndex = UnivariateIndetIndex(f);
      // if (IndetIndex < 0)
      //   CoCoA_THROW_ERROR(ERR::NotUnivariate, FnName);

      // if (IsQQ(P))
      // {
      //   const ring &ZZ = RingZZ();
      //   const SparsePolyRing ZZx = NewPolyRing(ZZ, NewSymbols(NumIndets(Px)));
      //   const RingHom toZZx = PolyRingHom(Px, ZZx, QQEmbeddingHom(ZZ), indets(ZZx));
      //   return IsLRSDegenerateIter(toZZx(f / content(f)));
      // }

      // if (!IsZZ(P))
      //   CoCoA_THROW_ERROR(ERR::BadRing, FnName);
      // const RingElem &x = indet(Px, IndetIndex);

      // // Must have non-zero constant coefficient
      // // TODO: Better or faster way?
      // f /= gcd(f, IndetPower(Px, IndetIndex, deg(f)));

      // // Extract square-free part
      // f /= gcd(f, deriv(f, IndetIndex));

      // const long degf = deg(f);
      // if (degf < 2)
      //   CoCoA_THROW_ERROR("Polynomial does not have two distinct non-zero roots.", FnName);

      const ring& Px = owner(f);
      const RingElem& x = indet(Px,0);
      if (!IsConstant(gcd(f, PolyAlgebraHom(Px, Px, std::vector<RingElem>(NumIndets(Px), -x))(f))))
        return 2;

      BigInt bound = InvTotientBoundUpto(deg(f)); // * (degf-1) ??? TODO
      unsigned long UPB;
      if (!IsConvertible(UPB, bound))
        CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);

      // Skip first two Graeffe iterations
      GraeffeSeq graeffeSeq(f, UPB);
      ++graeffeSeq;
      ++graeffeSeq;

      VERBOSE(20) << "n bound = " << UPB << std::endl;
      for (unsigned long n = 3; n <= UPB; ++n)
      {
        VERBOSE(100) << "n = " << n << std::endl;

        ++graeffeSeq;
        if (!IsSqFree(*graeffeSeq))
          return n;
      }
      return 0;
    }

  } // end of namespace OBSOLETE
  

  namespace // anonymous
  {
    /**
     * Used in modular approach to decide whether
     * the n-th power Graeffe transform of f needs
     * to be checked for proving its LRS-degeneracy.
     * For private use.  Hence, no arg checking.
     */
    bool IsLRSDegenerateOrder_ModularCheck(ConstRefRingElem f, const unsigned long n, const long NumPrimesToTry)
    {
      // ASSUME: f univariate, sqfree, f(0) non-zero, deg(f) >= 2, coeffs integer
      const char* const FnName = "IsLRSDegenerateOrder_ModularCheck";
      VerboseLog VERBOSE(FnName);

      const ring &Px = owner(f);
      const ring &P = CoeffRing(Px);
      const long nvars = NumIndets(Px);
      const BigInt LCf = ConvertTo<BigInt>(LC(f));

      const unsigned long gcdBound = n/2;

      static const auto X = NewSymbols(1); // do this once to avoid "burning through" lots of anon symbs

      PrimeSeq1ModN primeSeq(n); while (*primeSeq < 16*n) ++primeSeq; // skip past "very small primes"
      for (long i = 0; i < NumPrimesToTry; ++i)
      {
        SmallPrime p = *primeSeq;
        while (p != 0 && LCf%p == 0) { ++primeSeq; p = *primeSeq; }

        if (IsZero(p)) // upper prime bound has been reached
          break;

        // zeta will be a primitive root of unity in F_p
        const QuotientRing Fp = NewZZmod(p);
        const SparsePolyRing Fpx = NewPolyRing(Fp, X);
        const RingElem& x = indet(Fpx,0);
        const std::vector<RingElem> images(nvars, x);
        const RingElem f_p = PolyRingHom(Px, Fpx, CanonicalHom(P, Fp), images)(f);
        const unsigned long zeta = PowerMod(PrimitiveRoot(p), (p-1)/n, p);
        VERBOSE(120) << "zeta: " << zeta << std::endl;

        for (unsigned long k = 1; k <= gcdBound; ++k)
        {
          CheckForInterrupt(FnName);
          if (IsCoprime(n, k) &&
              IsConstant(gcd(f_p, PolyAlgebraHom(Fpx, Fpx, {PowerMod(zeta, k, p) * x})(f_p))))
          {
            return false;
          }
        }
        ++primeSeq;
      }

      return true;
    }



    /**
     * Used in modular approach to decide whether
     * the n-th power Graeffe transform of f needs
     * to be checked for proving its LRS-degeneracy.
     * For private use.  Hence, no arg checking.
     */

    // Input: D is the degree (assumed >= 2)
    // Return list of all indexes k s.t. cyclo(k) could divide LRS test resultant
    std::vector<unsigned long> GetKToTryList(const unsigned long D/*, const unsigned long Kbound*/)
    {
      CoCoA_ASSERT(D > 1);
      CoCoA_ASSERT(std::numeric_limits<unsigned long>::max()/D > D); // D*D does not overflow

      const unsigned long Kbound = InvTotientBoundUpto_ulong(D*(D-1));
//???    VERBOSE(20) << "K bound = " << Kbound << std::endl;

      // Create list of values which phi(k) must divide
      // List is repr as a "boolean mask": mask[d] == 1 iff d is in list
      std::vector<unsigned char> mask(D*D);  // vector<bool>  ???
      for (unsigned long d=2; d <= D; ++d)
        mask[d*(d-1)] = 1;

      for (unsigned long a = 1; a < D; ++a)
      {
        unsigned long UPB = min(a, D-a);
        for (unsigned long b = 1; b <= UPB; ++b)
        {
          if (IsEven(a * b))
            mask[a * b] = 1;
        }
      }

      // AllEvenFactorsMask is really a boolean array:
      // entry d is true iff d is even and divides one of the values in "mask"
      std::vector<unsigned char> AllEvenFactorsMask(D*D);  // vector<bool>  ???
      for (unsigned long d=D*(D-1); d >= 2; d -= 2) // start from biggest, only even values.
      {
        if (mask[d] == 0) continue;
        if (AllEvenFactorsMask[d] != 0) continue;
        // Now mark all even factors of d:
        const unsigned long half = d/2;
        for (unsigned long fac=1; fac*fac <= half; ++fac)
          if (half%fac == 0)
          { AllEvenFactorsMask[2*fac] = 1; AllEvenFactorsMask[d/fac] = 1; }
      }

      // Find all k values such that phi(k) divides some elem of NewL
      std::vector<unsigned long> KToTry;
      for (unsigned long k = 3; k <= Kbound; ++k)
      {
        const unsigned long phi = EulerTotient(k);
        if (phi < D*D && AllEvenFactorsMask[phi] != 0)
          KToTry.push_back(k);
      }
      return KToTry;
    }

  }  // end of namespace anonymous



  /**
   * Checks for LRS-degeneracy of the given order, without arg checks
   */
  bool IsLRSDegenerateOrderMod_NoArgChecks(ConstRefRingElem f, const unsigned long n, const long IndetIndex, VerificationLevel VerifLev)
  {
    // ASSUME: f is elem of ZZ[x], sqfr, f(0) != 0, deg(f) >= 2
    const char* const FnName = "IsLRSDegenerateOrderMod_NoArgChecks";
    VerboseLog VERBOSE(FnName);
    CoCoA_ASSERT(n >= 2);
///    if (n < 2)
///      CoCoA_THROW_ERROR(ERR::OutOfRange, FnName);

    if (n == 2)
    {
      const ring &Px = owner(f);
      const RingElem &x = indet(Px, IndetIndex);
      // Equiv to !IsCoprime
      return (!IsConstant(gcd(f, PolyAlgebraHom(Px, Px, std::vector<RingElem>(NumIndets(Px), -x))(f))));
    }

    const long NumPrimesToTry = IsGuaranteed(VerifLev) ? 3 : level(VerifLev); // magic number
    if (!IsLRSDegenerateOrder_ModularCheck(f, n, NumPrimesToTry)) return false;

      VERBOSE(20) << "Exponent " << n << " could not be excluded" << std::endl;

      if (!IsGuaranteed(VerifLev)) return true;

      // Also possible: Res_y(f(y), cyclotomic(n,y))   BUT SEEMS TO BE MUCH SLOWER!!!
      const RingElem Pn = GraeffeN(f, n); // equal to resultant(eval(f, [y]), y^n-x, y);
      return (!IsSqFree(Pn));
  }



  /**
   * Checks for LRS-degeneracy of the given order using the given verification level
   */
  bool IsLRSDegenerateOrderMod(RingElem f, const unsigned long n, VerificationLevel VerifLev)
  {
    const char* const FnName = "IsLRSDegenerateOrderMod";
    VerboseLog VERBOSE(FnName);
    const ErrorContext context = CoCoA_ERROR_CONTEXT;
    f = LRSDegenerate_check_and_preprocess(f, context);
    if (deg(f) < 2) return false;
//        CoCoA_THROW_ERROR_WITH_CONTEXT("Polynomial does not have two distinct non-zero roots.", OrigContext);

    const ring& Px = owner(f);
    const ring& ZZ = RingZZ(); // just an abbr.
    const SparsePolyRing ZZx = NewPolyRing(ZZ, symbols("x"));
    vector<RingElem> images(NumIndets(Px), indet(ZZx,0));
    const RingHom toZZx = PolyRingHom(Px, ZZx, CanonicalHom(CoeffRing(Px),ZZ), images);
    f = toZZx(prim(f));

    return IsLRSDegenerateOrderMod_NoArgChecks(f, n, 0/*IndetIndex*/, VerifLev);
  }


  /**
   * Checks for LRS-degeneracy of the given order
   */
  bool IsLRSDegenerateOrder(ConstRefRingElem f, const unsigned long n)
  {
    return IsLRSDegenerateOrderMod(f, n, guaranteed());
  }


  enum class FirstOrAll { JustFirst, FindAll };


  /**
   * Modular approach using the given verification level
   */
  std::vector<unsigned long> LRSDegeneracyOrder_modular(RingElem f, VerificationLevel VerifLev, FirstOrAll StopCriterion)
  {
    const char* const FnName = "LRSDegeneracyOrder_modular";
    VerboseLog VERBOSE(FnName);

    const ErrorContext context = CoCoA_ERROR_CONTEXT;
    f = LRSDegenerate_check_and_preprocess(f, context);
    vector<unsigned long> ListOfOrders;
    if (deg(f) < 2)  return ListOfOrders;
//        CoCoA_THROW_ERROR_WITH_CONTEXT("Polynomial does not have two distinct non-zero roots.", OrigContext);
    const ring& Px = owner(f);
    const ring& ZZ = RingZZ(); // just an abbr
    const SparsePolyRing ZZx = NewPolyRing(ZZ, symbols("x"));
    vector<RingElem> images(NumIndets(Px), indet(ZZx,0));
    const RingHom toZZx = PolyRingHom(Px, ZZx, CanonicalHom(CoeffRing(Px),ZZ), images);
    f = toZZx(prim(f));

    const RingElem& x = indet(ZZx,0);
    vector<RingElem> NegImages{-x};
    RingHom NegX = PolyAlgebraHom(ZZx, ZZx, NegImages);
    if (!IsConstant(gcd(f, NegX(f))))
    {
      ListOfOrders.push_back(2);
      if (StopCriterion == FirstOrAll::JustFirst) return ListOfOrders;
    }

//???    if (IsLRSDegenerateOrderMod_NoArgChecks(f, 2, 0/*IndetIndex*/, VerifLev))
//???      return 2;

    const long degf = deg(f);
    if (degf > std::numeric_limits<long>::max()/degf)  // I'd be amazed if this ever happens!
      CoCoA_THROW_ERROR("Deg too big", "LRSDegeneracyOrder_modular");
//     const unsigned long Kbound = InvTotientBoundUpto_ulong(degf*(degf-1));
// //     const BigInt KBOUND = InvTotientBoundUpto(degf*(degf-1));
// // //    const BigInt KBOUND = InvTotientBoundUpto(degf); //  BUG!!! depends on conjecture
// //     unsigned long Kbound;
// //     if (!IsConvertible(Kbound, KBOUND))
// //       CoCoA_THROW_ERROR(ERR::ArgTooBig, FnName);

//    VERBOSE(20) << "K bound = " << Kbound << std::endl;

    const std::vector<unsigned long> KToTry = GetKToTryList(degf/*, Kbound*/);
    VERBOSE(80) << "How many k to try? " << len(KToTry) << std::endl;
    for (unsigned long k : KToTry)
    {
      CheckForInterrupt("LRSDegeneracyOrder_modular");
      VERBOSE(100) << "Checking k = " << k << std::endl;
      if (IsLRSDegenerateOrderMod_NoArgChecks(f, k, 0/*IndetIndex*/, VerifLev))
      {
        ListOfOrders.push_back(k);
        if (StopCriterion == FirstOrAll::JustFirst) return ListOfOrders;
      }
    }
    return ListOfOrders;
  }


  /**
   * Chooses the right approach for the given polynomial
   */
  unsigned long LRSDegeneracyOrder(ConstRefRingElem f)
  {
    vector<unsigned long> orders = LRSDegeneracyOrder_modular(f, guaranteed(), FirstOrAll::JustFirst);
    if (!orders.empty()) return orders[0];
    return 0;
  }

  std::vector<unsigned long> LRSDegeneracyOrders(ConstRefRingElem f)
  {
//    return LRSDegeneracyOrder_modular(f, VerificationLevel(3), FirstOrAll::FindAll);
    return LRSDegeneracyOrder_modular(f, guaranteed(), FirstOrAll::FindAll);
  }


  bool IsLRSDegenerate(RingElem f)
  {
      f = LRSDegenerate_check_and_preprocess(f, CoCoA_ERROR_CONTEXT);
      // We do a quick check for cyclo factors before trying general method.
      const auto facs = CyclotomicFactors(f);
      if (!facs.myFactors().empty())
      {
        // Definitely LRS-degenerate unless the only cyclo factor is x-1
        // That's what we're checking for here... terribly long-winded!!
        if (len(facs.myFactors()) > 1 ||
            deg(facs.myFactors()[0]) > 1 ||
            IsOne(ConstantCoeff(facs.myFactors()[0])))
          return true;
      }
    return (LRSDegeneracyOrder(f) != 0);
  }

}


// RCS header/log in the next few lines
// $Header: /Volumes/Home_1/cocoa/cvs-repository/CoCoALib-0.99/src/AlgebraicCore/SparsePolyOps-LRSDegeneracy.C,v 1.14 2023/05/11 19:27:53 abbott Exp $
// $Log: SparsePolyOps-LRSDegeneracy.C,v $
// Revision 1.14  2023/05/11 19:27:53  abbott
// Summary: Changes spacing
//
// Revision 1.13  2023/04/12 20:11:56  abbott
// Summary: GetKToTryList now computes the Kbound (instead of taking it as input)
//
// Revision 1.12  2023/03/28 14:45:35  abbott
// Summary: Made interruptible
//
// Revision 1.11  2023/03/15 21:32:36  abbott
// Summary: Renamed CycloFactors to CyclotomicFactors; also refined semantics
//
// Revision 1.10  2023/03/06 22:06:56  abbott
// Summary: Added new fn LRSDegeneracyOrders (returns a vector)
//
// Revision 1.9  2023/02/27 21:10:01  abbott
// Summary: Several minor improvements & bug fixes (sigh!)
//
// Revision 1.8  2023/02/23 20:49:10  abbott
// Summary: Added new fn LRSDegeneracyOrder (was IsLRSDegenerate); major code revision
//
// Revision 1.7  2023/02/16 20:48:28  abbott
// Summary: Replaced calls to InvTotientBound (wrong) by InvTotientBoundUpto (correct)
//
// Revision 1.6  2023/01/31 12:33:59  abbott
// Summary: Major tidying
//
// Revision 1.5  2022/11/30 20:38:08  abbott
// Summary: Fixed subtle bug
//
// Revision 1.4  2022/11/30 15:19:12  abbott
// Summary: Many changes; revised indentation
//
// Revision 1.3  2022/11/24 21:11:29  abbott
// Summary: Added copyright notice
//
//

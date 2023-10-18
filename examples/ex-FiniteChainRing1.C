#include "CoCoA/library.H"

#include <algorithm>
using namespace std;

namespace CoCoA
{
    void program()
    {
        GlobalManager CoCoAFoundations;

        ring Fp = NewZZmod(31);         // coefficient ring
        ring Fpwyz = NewPolyRing_DMPII(Fp, symbols("w,y,z"), lex);
        RingElem w = RingElem(Fpwyz, symbol("w"));
        RingElem y = RingElem(Fpwyz, symbol("y"));
        RingElem z = RingElem(Fpwyz, symbol("z"));
        ideal J = ideal(power(w,2), power(y,3), power(z,3), w*y*z*z);
        cout << "J  = " << J << endl;
        cout << "GBasis(J)  = " << GBasis(J) << endl;  
    }
}

int main()
{
  try
  {
    CoCoA::program();
    return 0;
  }
  catch (const CoCoA::ErrorInfo& err)
  {
    cerr << "***ERROR***  UNCAUGHT CoCoA error";
    ANNOUNCE(cerr, err);
  }
  catch (const std::exception& exc)
  {
    cerr << "***ERROR***  UNCAUGHT std::exception: " << exc.what() << endl;
  }
  catch(...)
  {
    cerr << "***ERROR***  UNCAUGHT UNKNOWN EXCEPTION" << endl;
  }

  CoCoA::BuildInfo::PrintAll(cerr);
  return 1;
}
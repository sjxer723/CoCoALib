#include "CoCoA/library.H"

#include <algorithm>
using namespace std;

namespace CoCoA
{
    void program()
    {
        GlobalManager CoCoAFoundations;

        ring Fp = NewZZmod(31);         // coefficient ring
        ring Fp1 = NewZZmod(32);         // coefficient ring
        ring Fpwyz = NewPolyRing_DMPII(Fp, symbols("w,y,z"), lex);
        ring Fgalois = NewGaloisRing(2, 10);   
        cout << "Is Ring Fp" << IsRingFp(Fgalois) << endl;
        cout << "Is Ring Fp" << IsRingFp(Fp) << endl;
        cout << "Is Ring Fp" << IsRingFp(Fp1) << endl;
        
        ring Fgaloispwyz = NewPolyRing(Fgalois, symbols("w1,y1,z1"), lex);
        cout << IsPowerOf2(Fgalois-> myCharacteristic()) << endl;
        cout << "Is Finite Field" << IsFiniteField(Fgalois) << endl;
        cout << "Is Field" << IsField(Fgalois) << endl;
        RingElem w = RingElem(Fpwyz, symbol("w"));
        RingElem y = RingElem(Fpwyz, symbol("y"));
        RingElem z = RingElem(Fpwyz, symbol("z"));
        cout << "21" << endl;
        RingElem w1 = RingElem(Fgaloispwyz, symbol("w1"));
        cout << "23" << endl;
        RingElem y1 = RingElem(Fgaloispwyz, symbol("y1"));
        cout << "25" << endl;
        RingElem z1 = RingElem(Fgaloispwyz, symbol("z1"));
        
        ideal J = ideal(w - 2*y, w + y);
        cout << "J" << endl;
        ideal J1 = ideal(512*w1*w1 + 2*w1, w1 - 2*y1, w1 + y1);
        // cout << "J  = " << J << endl;
        // cout << "GBasis(J)  = " << GBasis(J) << endl;  
        cout << "J1  = " << J1 << endl;
        cout << "GBasis(J1)  = " << GBasis(J1) << endl;  
      
    }
}

int main()
{
  CoCoA::SetVerbosityLevel(200); // change arg to get more/less verbose logging.
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
#include "CoCoA/library.H"

#include <algorithm>
using namespace std;

namespace CoCoA
{  
  // Test GetReverseOverPowerOf2 Function
  void program1()
  {
    BigInt res;
    GetReverseOverPowerOf2(res, BigInt(MachineInt(5)), 8);
    std::cout << res << std::endl;    
  }

  void program2()
  {
    GlobalManager CoCoAFoundations;

    ring Fgalois = NewGaloisRing(2, 8);   
    ring Fgaloispwyz = NewPolyRing(Fgalois, symbols("w,y,z"), lex);
    cout << IsPowerOf2(Fgalois-> myCharacteristic()) << endl;
    // cout << "Is Finite Field" << IsFiniteField(Fgalois) << endl;
    // cout << "Is Field" << IsField(Fgalois) << endl;
    RingElem w = RingElem(Fgaloispwyz, symbol("w"));
    RingElem y = RingElem(Fgaloispwyz, symbol("y"));
    RingElem z = RingElem(Fgaloispwyz, symbol("z"));
    
    ideal J = ideal(16*w + 8*y +4*z, 6*y);
    cout << "GBasis(J)  = " << GBasis(J) << endl; 

    ideal J1 = ideal(16*w*w + 8*y +4*z, 8*w + 4*z*z, 6*y);
    cout << "GBasis(J1)  = " << GBasis(J1) << endl;  
  
    ideal J2 = ideal(w*w + 8*y +4*z, 8*w + 4*z*z, 6*y);
    cout << "GBasis(J2)  = " << GBasis(J2) << endl;  

    ideal J3 = ideal(8*w + 2, 6*w);
    cout << "GBasis(J3)  = " << GBasis(J3) << endl;  

    std::vector<RingElem> v = {16*w*y + 4*y +2*z, 8*y + 4*z*z, 32*z};
    ideal J4 = ideal(v);
    cout << "GBasis(J4)  = " << GBasis(J4) << endl;  

    ideal J5 = ideal(3*w - 2*y*y, 3*w + y);
    cout << "GBasis(J5)  = " << GBasis(J5) << endl;

    ideal J6 = ideal(3*w*y - 2*y*y*z, 3*w + y, 5*z);
    cout << "GBasis(J6)  = " << GBasis(J6) << endl;   

    ideal J7 = ideal(85*w*y, 64*w + 4*z, 12*z);
    cout << "GBasis(J7)  = " << GBasis(J7) << endl;    
  }

  void program3()
  {
    BigInt ans1, ans2;
    BigInt x(36);
    BigInt y(48);
    LcmDivisor(ans1, ans2, x, y, 8);    
    cout << "ans 1: " << ans1 << " " << "ans2: " << ans2 << endl;
  }
}

int main()
{
  // CoCoA::SetVerbosityLevel(200); // change arg to get more/less verbose logging.
  try
  {
    CoCoA::program2();
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
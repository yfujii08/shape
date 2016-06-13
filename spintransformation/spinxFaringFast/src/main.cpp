// =============================================================================
// SpinYForm is for spintransformation-based fairing
// that is based on combination of DDG and spinxform by Keenan Crane
// Ryo Yamada
// Feburuary 18, 2017
// main.cpp
//

#include <iostream>
#include <sstream>
#include <fstream>

#include "MeshY.h"
#include "Application.h"



using namespace std;
using namespace DDG;

// argv[1] : obj file. hogehoge.obj containing vertices v\sx\sy\sz and edges f\x\sy\sz
// argv[2] : the step parameter during fairing process. 0.9 means make it plane by changing curvature 10%.
// argv[3] : the number of steps.
// argv[4] : output file name. hoge outputs hoge0.obj, hoge1.obj, ...
// ./fairing hogei.obj can be viewed.
int main( int argc, char **argv )
{

   int n = atoi( argv[3] ); //‘æ3ˆø”‚ğ®”Œ^‚É•ÏŠ·‚µ‚Än‚É‘ã“ü
   double step = atof( argv[2] ); //‘æ2ˆø”‚ğ­”Œ^‚É•ÏŠ·‚µ‚Äx‚É‘ã“ü
   // load mesh
   spiny::Mesh meshY;
   meshY.read( argv[1] ); //‘æ1ˆø”‚Íobjƒtƒ@ƒCƒ‹–¼

   // processing with Application.h
   int counter = 0;
   Application app;
   for(int i=0;i<n;i++){
     app.process(step,meshY);
     std::ostringstream os;
     os << argv[4] << i << ".obj";
     std::string name=os.str();
     meshY.write(name);
     counter += 1;
   }


   return 0;
}


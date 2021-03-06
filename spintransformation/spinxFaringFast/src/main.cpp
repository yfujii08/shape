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

   int n = atoi( argv[3] ); //第3引数を整数型に変換してnに代入
   double step = atof( argv[2] ); //第2引数を少数型に変換してxに代入
   // load mesh
   spiny::Mesh meshY;
   meshY.read( argv[1] ); //第1引数はobjファイル名

   // processing with Application.h
   int counter = 0;
   Application app;
   for(int i=0;i<n;i++){
     app.process(step,meshY);
     std::ostringstream os;
     os << argv[4] << i+1 << ".obj";
     std::string name=os.str();
     meshY.write(name);
     counter += 1;
   }

   std::vector<double> data(meshY.vertices.size(), 0.28209479177387814);  // 0.28209479177387814 // 1/2/sqrt(pi)
   std::vector<double> coef = meshY.SH_innerProduct(30);//, &data=rho);
   std::ofstream ofs("coef.txt");
   for (auto c : coef)
   {
      ofs << c << " ";
   }

   return 0;
}


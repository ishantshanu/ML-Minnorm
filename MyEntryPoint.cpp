//#include "stdafx.h"
//#include "MyEntryPoint.h"
//#include <conio.h>
//
////#include "test/HandcraftedTest.h"
//#include "test/DeblurringTest.h"
//#include "test/DenoisingTest.h"
//#include "test/CliquePotential.h"
//#include "test/StereoTest.h"
//#include <fstream>
//using namespace std;
//
////void KSTest()
////{
////    Petter::PseudoBoolean<double> f;
////    f.add_clique(0,20,20);
////    f.add_clique(1,20,10);
////    f.add_clique(0,1,0,10,10,100);
////    
////    vector<Petter::label> x(2); // Holds the solution x
////    int labeled; // Will hold the number of labeled variables
////    f.minimize(x, labeled, Petter::GRD);
////    cout << "Optimal relaxation labeled " << labeled << " variables\n" << x[0] << "," << x[1] << endl;
////    f.minimize(x, labeled, Petter::GRD_heur);
////    cout << "Heuristic relaxation labeled " << labeled << " variables\n" << x[0] << "," << x[1] << endl;
////    f.minimize(x, labeled, Petter::HOCR);
////    cout << "HOCR labeled " << labeled << " variables\n" << x[0] << "," << x[1] << endl;
////}
//
//void TestSubmodularity()
//{
//    CliquePotential cp(3);
//    memset(cp._labelingCosts, 0, 8*sizeof(float));
//    cp._labelingCosts[6] = 10;
//    if(cp.IsSubmodular())
//        printf("potential is submodular\n");
//    else
//        printf("potential is NOT submodular\n");
//}
//
//void MyEntryPoint::MyMain(vector< vector< int > >cliqueinfo, int total_nodes, std::vector< double >x_neg_V)
//{
//    //STATE_FILE = fopen("./out-images/log.txt", "w");change made
//
//    //HandcraftedTest::Test();
//
//    //TestSubmodularity();
//
//    //StereoTest::Test();
//    DenoisingTest::Test(cliqueinfo,total_nodes,x_neg_V);        //change made
//    
//    //DeblurringTest::Test();
//
//    //fclose(STATE_FILE);Change made
//
////#ifdef _DEBUG
////    printf("press any key to continue....\n");
////    _getch();
////#endif
//}

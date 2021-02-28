// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include <iostream>
#include <unordered_map>
#include <stdio.h>
#include <unordered_set>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <vector>
#include <list>
#include "opencv2/opencv.hpp"
#include "Eigen/Dense"
#include <chrono>
#include <fstream>
#include <stdio.h>
#include <iomanip>

using namespace std;
using namespace cv;
using namespace Eigen;


typedef int NodeIndex;
typedef vector<int> Clique;

const int MAX_NUM_CLIQUES = 25000;
const double EPSILON = 0.000001;
const double INF = 10000;
extern double EPSILON2;
const double UNARY_FACTOR2 = .5;
const double UNARY_FACTOR1 = .5;
//
//class Oracle;
//
//extern Oracle *ORACLE;


// TODO: reference additional headers your program requires here

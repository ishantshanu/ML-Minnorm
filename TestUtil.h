#pragma once
#include "Oracle.h"
class TestUtil {
public:
    cv::Mat GetSolutionAsImage(std::unordered_set<int> &minimizer_set, double &dual, double &primal,
                               std::vector <Clique> numclique, map<int, int> labelmapback);
	TestUtil(Oracle* ORACLE);
	~TestUtil();
	Oracle* ORACLE;
	uchar* RED;// = new uchar[21]{ 0, 128, 0, 128, 0, 128, 0, 128, 64, 192, 64, 192, 64, 192, 64, 192, 0, 128, 0, 128, 0 };
	uchar* GREEN;// = new uchar[21]{ 0, 0, 128, 128, 0, 0, 128, 128, 0, 0, 128, 128, 0, 0, 128, 128, 64, 64, 192, 192, 64 };
	uchar* BLUE;// = new uchar[21]{ 0, 0, 0, 0, 128, 128, 128, 128, 0, 0, 0, 0, 128, 128, 128, 128, 0, 0, 0, 0, 128 };
};

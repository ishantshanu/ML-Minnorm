#include "stdafx.h"
#include "GenerateTable.h"
#include "TestSegmentation.h"
#include <fstream>
int main(int argc, char* argv[])
{

	std::ifstream file("imagenames.txt");
	std::string line;
	while (getline(file, line))
	{
		std::istringstream iss(line);
		string x;
		iss >> x;
		cout << x << endl;
		TestSegmentation *test = new TestSegmentation;
		std::cout << "time = " << test->Test(x);
		delete test;
	}
	//string x = "2007_000799";
	//TestSegmentation test;
	//test.Test(x,classname);
	
	return 0;
}

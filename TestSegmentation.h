#pragma once

class TestSegmentation
{
public:
	uchar* RED;// = new uchar[21]{ 0, 128, 0, 128, 0, 128, 0, 128, 64, 192, 64, 192, 64, 192, 64, 192, 0, 128, 0, 128, 0 };
	uchar* GREEN;// = new uchar[21]{ 0, 0, 128, 128, 0, 0, 128, 128, 0, 0, 128, 128, 0, 0, 128, 128, 64, 64, 192, 192, 64 };
	uchar* BLUE;// = new uchar[21]{ 0, 0, 0, 0, 128, 128, 128, 128, 0, 0, 0, 0, 128, 128, 128, 128, 0, 0, 0, 0, 128 };
	~TestSegmentation();
	TestSegmentation();

	double Test(string name);
	double TestSegmentation::GetIou(cv::Mat min_pred, cv::Mat segnet_pred, cv::Mat gt);
	int TestSegmentation::computenumlabels(string ss);
};

//
//TestSegmentation::TestSegmentation()
//{
//	RED = new uchar[21];// { 0, 128, 0, 128, 0, 128, 0, 128, 64, 192, 64, 192, 64, 192, 64, 192, 0, 128, 0, 128, 0 };
//	GREEN = new uchar[21];// { 0, 0, 128, 128, 0, 0, 128, 128, 0, 0, 128, 128, 0, 0, 128, 128, 64, 64, 192, 192, 64 };
//	BLUE = new uchar[21];// { 0, 0, 0, 0, 128, 128, 128, 128, 0, 0, 0, 0, 128, 128, 128, 128, 0, 0, 0, 0, 128 };
//}

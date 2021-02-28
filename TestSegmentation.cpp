#include "stdafx.h"
#include "TestSegmentation.h"
#include "TestUtil.h"
#include "Minnorm.h"
#include <time.h>
#include <fstream>
double EPSILON2;


TestSegmentation::~TestSegmentation()
{
	delete[] RED;
	delete[] GREEN;
	delete[] BLUE;
}


TestSegmentation::TestSegmentation()
{
	RED = new uchar[21]{ 0, 128, 0, 128, 0, 128, 0, 128, 64, 192, 64, 192, 64, 192, 64, 192, 0, 128, 0, 128, 0 };
	GREEN = new uchar[21]{ 0, 0, 128, 128, 0, 0, 128, 128, 0, 0, 128, 128, 0, 0, 128, 128, 64, 64, 192, 192, 64 };
	BLUE = new uchar[21]{ 0, 0, 0, 0, 128, 128, 128, 128, 0, 0, 0, 0, 128, 128, 128, 128, 0, 0, 0, 0, 128 };
}


int TestSegmentation::computenumlabels(string ss)
{
	cv::Mat unary = imread("pred_deep/" + ss + "_pred.png");
	int count = 1;
	bool *visitlabel = new bool[21];
	for (int label = 0; label < 21; label++)
		visitlabel[label] = 0;
	for(int i=0; i< unary.rows ; i++)
		for (int j = 0; j < unary.cols; j++) {
			for (int label = 0; label < 21; label++) {
				if (unary.at<Vec3b>(i, j) == Vec3b(BLUE[label], GREEN[label], RED[label])) {
					if (!visitlabel[label]) {
						visitlabel[label] = 1; 
						count++;
					}
				}
			}
		}
	delete[] visitlabel;
	unary.release();
	return count;
}

double TestSegmentation::Test(string name)
{
	string ss = name;
	//cout <<"\n \n name = "<< ss << endl;
	int img_cols = 200, img_rows = 200;	
	int numLabel = computenumlabels(name);
	cout <<"numLabel  "<< numLabel<<endl;
	int MAX_NUM_CLIQUES = 100000;
	EPSILON2 = 10000;
	 
	double **UF;
	UF = new double*[numLabel + 2];
	for (int i = 0; i < numLabel + 2; i++) {
		UF[i] = new double[numLabel + 2];
	}
	// unary factor
	for (int i = 0; i < numLabel; i++)
		for (int j = 0; j < numLabel; j++)
			UF[i][j] = 5;

	cv::Mat gt;
	gt = cv::imread("gt/" + name + "_gt.png");

	// make images for color reference
	for (int i = 0; i < 21; i++)
	{
		cv::Mat imgtmp = cv::Mat(img_rows, img_cols, CV_8UC3, cv::Scalar(0));
		uchar blue, green, red;
		blue = BLUE[i];
		green = GREEN[i];
		red = RED[i];
		for (int x = 0; x < img_rows; x++)
		{
			for (int y = 0; y < img_cols; y++)
			{
				imgtmp.at<Vec3b>(x, y) = Vec3b(blue, green, red);
			}
		}
		cv::imwrite("label/label_" + std::to_string(i) + ".png", imgtmp);
		imgtmp.release();
	}


	cv::Mat image; 

	std::string line;
	
	/*
		------ READ CLIQUES -------
	*/

	// --------- EPSILON2 -------
	std::vector<Clique> temp_cliques(MAX_NUM_CLIQUES);
	std::ifstream infile("segments/segment_" + ss + ".txt");
	int x, y, clique_index;


	while (getline(infile, line))
	{
		std::istringstream iss(line);
		iss >> y;
		iss >> x;
		iss >> clique_index;
		temp_cliques[clique_index].push_back(y * img_cols + x);
	}

	std::vector<Clique> cliques;
	for (int i = 0; i < MAX_NUM_CLIQUES; ++i)
	{
		if (!temp_cliques[i].empty())
			cliques.push_back(temp_cliques[i]);
	}
	temp_cliques.clear();

	int maxz = 0;
	for (int i = 0; i < (int)cliques.size(); i++) {
		if (maxz < (int)cliques[i].size())
			maxz = (int)cliques[i].size();
	}
	std::cout << "max clique size = " << maxz << std::endl;
	cout << "number of cliques= " << cliques.size() << endl;
	vector<int> clsize;

	// expand the cliques for multilable
	std::vector<Clique> tmp(int(cliques.size()));
	for (int i = 0; i<int(cliques.size()); i++)
	{
		for (int j = 0; j<int(cliques[i].size()); j++)
		{
			int idx = cliques[i][j];
			for (int k = 1; k < (numLabel); k++)
				tmp[i].push_back(idx*(numLabel - 1) + (k - 1));
		}
	}
	cliques = tmp;
	tmp.clear();
	/*
		------- READ UNARY ------
	*/


	vector<vector<vector<double> > > unary(img_rows, vector<vector<double> >(img_cols, vector<double>(numLabel)));
	for (int i = 0; i < img_rows; i++)
		for (int j = 0; j < img_cols; j++)
			for (int k = 0; k < numLabel; k++)
				unary[i][j][k] = 0;

	std::ifstream file("unary/unary" + ss + ".txt");
	std::string line;
	cv::Mat imagetmp = cv::Mat(img_rows, img_cols, CV_8UC3, cv::Scalar(0));
	map<int, int> labelmap;
	map<int, int> labelmapback;
	// starts with 2 one for bacground one for class
	int labelcount = 1;
	// assign background label 0
	labelmap[0] = 0;
	labelmapback[0] = 0;
	while (getline(file, line))
	{
		std::istringstream iss(line);
		int x, y,idx;
		double value;
		iss >> x;
		iss >> y;
		iss >> value;
		iss >> idx;
		if (labelmap.find(idx) == labelmap.end()) {
			labelmap[idx] = labelcount;
			labelmapback[labelcount] = idx;
			labelcount++;
		}

		value = 255.0;

		for (int i = 0; i < numLabel; i++)
		{
			if (i == labelmap[idx])
				unary[x][y][i] = (255 - value);   // cost for label idx
			else
				unary[x][y][i] = value - 0;       // cost for lebel i
		}
		for (int i = 0; i < numLabel; i++)
			unary[x][y][i] = unary[x][y][i] * UF[labelmap[idx]][i];   // cost for label idx


		uchar blue, green, red;
		blue = BLUE[idx];
		green = GREEN[idx];
		red = RED[idx];
		imagetmp.at<Vec3b>(x, y) = Vec3b(blue, green, red);
		

		imwrite("label/res_unary_" + ss + ".png", imagetmp);
		Oracle *ORACLE = new Oracle(3, clsize, cliques, unary);
		int64 t1 = cvGetTickCount();
		double dual = 0;
		cout << "Entering minnorm" << endl;
		Minnorm minnorm(img_cols*img_rows*(numLabel - 1), numLabel, cliques,ORACLE);
		std::unordered_set<NodeIndex> minimizerset = minnorm.DoInference(dual, ss);
		cout << "Exited minnorm" << endl;
		int64 t2 = cvGetTickCount();
		double timet = (t2 - t1) / (cvGetTickFrequency()*1000.0) / 1000.0;
		TestUtil* testutil = new TestUtil(ORACLE);
		double primal;

		//// function to calculate the primal and form the final output image
		cv::Mat pred_image = testutil->GetSolutionAsImage(minimizerset, dual, primal, cliques, labelmapback);
		
		double iouvalue = GetIou(pred_image, imagetmp, gt);
		cout << "iouvalue = " << iouvalue<<endl; 
		

		delete ORACLE;
		pred_image.release();
		imagetmp.release();
		delete testutil;
		labelmap.clear();
		labelmapback.clear();
		cout << "time = " << timet << endl; 
	}


	gt.release();
	for (int i = 0; i < numLabel + 2; i++) {
		delete[] UF[i];
	}
	delete[] UF;
	cliques.clear();
	clsize.clear();

	return 0;
}


double TestSegmentation::GetIou(cv::Mat min_pred, cv::Mat segnet_pred, cv::Mat gt)
{

	double *intersectioncount  = new double[21];
	double *unioncount = new double[21];
	double *ioulabel = new double[21];

	for (int label = 0; label < 21; label++) {
		intersectioncount[label] = 0;
		unioncount[label] = 0;
		ioulabel[label] = 0;
	}
	double avgiouseg = 0,avgioumin = 0;
	for (int label = 0; label < 21; label++) {
		for(int i=0; i < min_pred.rows; i++)
			for (int j = 0; j < min_pred.cols; j++) {
				Vec3b gtv = gt.at<Vec3b>(i, j);
				Vec3b minpredv = min_pred.at<Vec3b>(i, j);
				if(std::tie(gtv[0], gtv[1], gtv[2]) == std::tie(minpredv[0], minpredv[1], minpredv[2]) && std::tie(minpredv[0], minpredv[1], minpredv[2]) == std::tie(BLUE[label], GREEN[label], RED[label]))
					intersectioncount[label]++;
				if (std::tie(minpredv[0], minpredv[1], minpredv[2]) == std::tie(BLUE[label], GREEN[label], RED[label]) || std::tie(gtv[0], gtv[1], gtv[2]) == std::tie(BLUE[label], GREEN[label], RED[label]))
					unioncount[label]++;
			}
		if (unioncount[label] == 0)
			ioulabel[label] = 0;
		else 
			ioulabel[label] = intersectioncount[label] *1.0 / unioncount[label];
		avgioumin += ioulabel[label];
	}
	avgioumin /= 21;

	for (int label = 0; label < 21; label++) {
		intersectioncount[label]=0;
		unioncount[label]=0;
		ioulabel[label]=0;
	}
	for (int label = 0; label < 21; label++) {
		for (int i = 0; i < segnet_pred.rows; i++)
			for (int j = 0; j < segnet_pred.cols; j++) {
				Vec3b gtv = gt.at<Vec3b>(i, j);
				Vec3b minpredv = segnet_pred.at<Vec3b>(i, j);
				if (std::tie(gtv[0], gtv[1], gtv[2]) == std::tie(minpredv[0], minpredv[1], minpredv[2]) && std::tie(minpredv[0], minpredv[1], minpredv[2]) == std::tie(BLUE[label], GREEN[label], RED[label]))
					intersectioncount[label]++; //variable name abuse
				if (std::tie(minpredv[0], minpredv[1], minpredv[2]) == std::tie(BLUE[label], GREEN[label], RED[label]) || std::tie(gtv[0], gtv[1], gtv[2]) == std::tie(BLUE[label], GREEN[label], RED[label]))
					unioncount[label]++;
			}
		if (unioncount[label] == 0)
			ioulabel[label] = 0;
		else
			ioulabel[label] = intersectioncount[label] * 1.0 / unioncount[label];
		avgiouseg += ioulabel[label];
	}
	avgiouseg /= 21;
	delete[] intersectioncount;
	delete[] unioncount;
	delete[] ioulabel;

	return avgioumin - avgiouseg;
}
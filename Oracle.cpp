#include "stdafx.h"
#include "Oracle.h"

  

Oracle::Oracle(int type, std::vector<int> clsize,vector<vector<int>> cliques , vector<vector<vector<double> > > unary)
	:_unary(unary), _clsize(clsize), _potentialtype(type), _cliques(cliques)
{
	UNARY_WIEGHT_FACTOR = 0.5;
	_imgHieght = unary.size();
	_imgWidth = unary[0].size();
	_imgLabel = unary[0][0].size(); 
}

Oracle::~Oracle()
{
	_unary.clear();
	_clsize.clear();
	_cliques.clear();
	_submodfunc.clear();
}

double Oracle::GetCost(std::vector<int> indicator,int clindex)
{
	if (_potentialtype == 0)
	{
		return (double) edgebasedpotential(indicator);
	}
	else if(_potentialtype == 1)
	{
		// ORACLE HAS EXTENDED CLIQUES -- TAKE CARE WHILE CODING !!!

		// when running with Graphcut or DP... oracle for illegal state is never called so 'cnt' should be zero then

		// this is Implementation of |x^bar| - |x|
		int num_pixel_clique = (int(_cliques[clindex].size()))/(_imgLabel-1);
		std::vector<int> mp1(num_pixel_clique,0);
		std::vector<int> mp2(num_pixel_clique,0);
		for(int i=0;i<int(num_pixel_clique);i++)
		{
			mp1[i]=0;
			mp2[i]=0;
		}

		for(auto w:indicator)
		{
			int p = w/(_imgLabel-1);   // pixel p corresponding to w
			int lab = ( w%(_imgLabel-1) )+1;	// these two lines essentially takes x^bar
			mp1[p]=max(mp1[p],lab);
			mp2[p]++;
		}

		double cnt = 0;
		for(int i=0;i<int(num_pixel_clique);i++)
		{
			cnt = cnt + (mp1[i]-mp2[i]);
		}

		return ( (concavepotential(indicator.size()+cnt, clindex)) + cnt*INF ); // no of 1's is indicator.size() + cnt since we are using f(x^bar)
	}
	else if(_potentialtype ==2)
	{
		int num_pixel_clique = (int(_cliques[clindex].size()))/(_imgLabel-1);
		std::vector<int> mp1(num_pixel_clique,0);
		std::vector<int> mp2(num_pixel_clique,0);
		for(int i=0;i<int(num_pixel_clique);i++)
		{
			mp1[i]=0;
			mp2[i]=0;
		}

		for(auto w:indicator)
		{
			int p = w/(_imgLabel-1);
			int lab = ( w%(_imgLabel-1) ) +1;
			mp1[p]=max(mp1[p],lab);
			mp2[p]++;
		}


		double ans=0;
		double cnt = 0;

		map<int,int> mp;		// mp stores frequency of each label
		for(int i=0;i<_imgLabel;i++)
		mp[i]=0;

		for(int i=0;i<int(num_pixel_clique);i++)
		{
			mp[mp1[i]]++;
		}

		for(int i=0;i<_imgLabel;i++)
		{
			for(int j=0;j<i;j++)
			{
				ans = ans + (mp[i]*mp[j]*(i-j));		// calculates sigma |a-b|
			}
		}
		for(int i=0;i<int(num_pixel_clique);i++)
		{
			cnt = cnt + (mp1[i]-mp2[i]);
		}
		if (cnt != 0) {
			//cout << "cnt!=0";
			//getchar();
		}
		//cout <<"ans   "<< ans << endl;
		return ans + cnt* INF;

	}
	if (_potentialtype == 3) {

		float alpha = .5;
		int num_pixel_clique = (int(_cliques[clindex].size())) / (_imgLabel - 1);
		std::vector<int> mp1(num_pixel_clique, 0);
		std::vector<int> mp2(num_pixel_clique, 0);
		for (int i = 0; i<int(num_pixel_clique); i++)
		{
			mp1[i] = 0;
			mp2[i] = 0;
		}
		for (auto w : indicator)
		{
			int p = w / (_imgLabel - 1);
			int lab = (w % (_imgLabel - 1)) + 1;
			mp1[p] = max(mp1[p], lab);
			mp2[p]++;
		}


		double ans = 0;
		double cnt = 0;

		map<int, int> mp;		// mp stores frequency of each label
		for (int i = 0; i < _imgLabel; i++)
			mp[i] = 0;

		for (int i = 0; i<int(num_pixel_clique); i++)
		{
			mp[mp1[i]]++;
		}

		for (int i = 0; i < _imgLabel; i++)
		{
			ans = ans + pow((num_pixel_clique - mp[i]), alpha);		// calculates sigma |a-b|
		}
		return ans*100;
	}

}

double Oracle::edgebasedpotential(std::vector<int> indicator)
{
	int index = 0;
	for (int i = 0; i < indicator.size(); i++)
		index+= 1 << indicator[i];
	return _submodfunc[index];
}

double Oracle::concavepotential(int size, int clindex)
{
	return (_cliques[clindex].size() - size)*size;
}
/*
bool Oracle::checkifsubmodular(int potentialtype)
{
	int clindex=1;
	vector<double> a_c = {};
	for (int i = 0; i < (1 << (_clsize[0] * _clsize[1])); i++){
		for (int j = 0; j < (1 << (_clsize[0] * _clsize[1])); j++){
			vector<int> indicator1, indicator2,indicator3,indicator4;
			for (int k = 0; k < _clsize[0] * _clsize[1]; k++)
			if (i & (1 << k))
				indicator1.push_back(k);
			for (int k = 0; k < _clsize[0] * _clsize[1]; k++)
			if (j & (1 << k))
				indicator2.push_back(k);
			for (int k = 0; k < _clsize[0] * _clsize[1]; k++)
			if ((i&j) & (1 << k))
				indicator3.push_back(k);
			for (int k = 0; k < _clsize[0] * _clsize[1]; k++)
			if ((i|j) & (1 << k))
				indicator4.push_back(k);
			if (GetCost(indicator1, clindex) + GetCost(indicator2, clindex) < GetCost(indicator3, clindex) + GetCost(indicator4, clindex)){
				cout << "NOT SUBMODULAR";
				getchar();
			}
		}
	}
	return 1;
}
*/
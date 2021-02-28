#pragma once

class Oracle
{
public:
	Oracle(int potential_type, std::vector<int> cliqusize, std::vector<vector<int>> cliques, vector<vector<vector<double> > >unary);
	~Oracle();
	std::vector<vector<vector<double> > > _unary;
	int _imgHieght, _imgWidth, _imgLabel;
	double _C(int clindex);
	//for real experiments
	//double Unarycost(int i, int j, int k) { return _unary[i][j][k] -((k == 0) ? 0 : _unary[i][j][k - 1]); }    //change
	//for artificial
	double Unarycost(int i, int j, int k) {
		if (k == 0)
			return _unary[i][j][k];
		else
			return _unary[i][j][k] - _unary[i][j][k - 1];	// telescopic 
	}
    //change
	double edgebasedpotential(std::vector<int> indicator);
	double concavepotential(int ind_size, int clindex);
	//double pairAbsolutePotential(std::vector<int> indicator); //TODO
	double GetCost(std::vector<int> indicator, int clindex);
	int _potentialtype;
	std::vector<int> _clsize;
	std::vector<vector<int>> _cliques;
	std::vector<int> _submodfunc;
	double UNARY_WIEGHT_FACTOR;
};

#pragma once
#include "Oracle.h"


class Ordering
{
public:
	~Ordering();
public:
	std::vector<int> _globalnodeindex;
	std::vector<int> _cliquenodeindex;
};

class ExtremeBase
{
public:
	ExtremeBase(int clique_index, Ordering ordering, Oracle* ORACLE);
	~ExtremeBase();

public:
	Ordering _ordering;
	Oracle *ORACLE;
	std::vector<double> _coordinateValues;
	double _lambda;
	int _cliqueIndex;

private:
	std::vector<double> GenerateCoordinatesUsingEdmondsAlgorithm();
};

class CliqueBase
{
public:
	CliqueBase(int clique_index, Clique& clique, Oracle* ORACLE);
	~CliqueBase();
	void ReduceOrdering(int clique_index);
	double DotProd(std::vector<double>, std::vector<double>);
	std::vector<double> Scalevector(std::vector<double>, double);
	std::vector<double> Differencofvectors(std::vector<double>, std::vector<double>);
	std::vector<double> GetConvexCombination();

public:
	std::vector<ExtremeBase*> _extremeBases;
	std::vector<double> _invalidEbCombination;


private:
};

class SoSBase
{
public:
	SoSBase(int num_nodes,int numLabel, std::vector<Clique>& cliques, Oracle* ORACLE);
	~SoSBase();

public:
	double DotProd(std::vector<double> A, std::vector<double> B);
	std::vector<double> Scalevector(std::vector <double> A, double scalar);
	std::vector<double> Differencofvectors(std::vector<double> A, std::vector<double> B);
	std::vector<double> Addofvectors(std::vector<double> A, std::vector<double> B);
	void AddUnary(void);
	int Addordering(int cliqueIndex, bool minnorm_only_flag, int &num_of_legal_bases, int &num_of_illegal_bases);
	std::vector<double> X_clique(int cliqueindex);
	std::vector<double> Computexminus(int cliqueindex);
	std::vector<CliqueBase*> _bases;
	Oracle* ORACLE;
	int _numPixel, _numLabel;
	vector<Clique> _cliques;
	double* _coordinates;
};

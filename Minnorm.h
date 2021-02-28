#include "stdafx.h"
#include "Base.h"

class Minnorm{

public:
	double _initialEnergy;
	int _numCliques;
	Oracle *ORACLE;
	int _numPixel, _imgHeight, _imgWidth, _numLabel;
	SoSBase* _sosbase;
	vector< set<int> > clique_pixel;	// store corresponding pixels of extended cliques

public:
	Minnorm(int numPixel,int numLabel, std::vector<Clique> cliques, Oracle *ORACLE);
	~Minnorm();
	MatrixXd computeminnormpoint(MatrixXd Bases);
	std::unordered_set<NodeIndex> DoInference(double& dual, string name);
};

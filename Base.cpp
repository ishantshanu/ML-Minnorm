#include "stdafx.h"
#include "Base.h"
#include "Oracle.h"

ExtremeBase::ExtremeBase(int clique_index, Ordering ordering,Oracle* ORACLE)
        : _cliqueIndex(clique_index), _ordering(ordering), ORACLE(ORACLE) {
    _coordinateValues = GenerateCoordinatesUsingEdmondsAlgorithm();
}

ExtremeBase::~ExtremeBase()
{
	_coordinateValues.clear();
}

Ordering::~Ordering()
{
	_cliquenodeindex.clear();
	_globalnodeindex.clear();
}

std::vector<double> ExtremeBase::GenerateCoordinatesUsingEdmondsAlgorithm()   
{
    double prev_value = 0, curr_value;
    int orderingsize = _ordering._cliquenodeindex.size();
    std::vector<double> res(orderingsize);
    std::vector<int> indicator; 
    for (int j = 0; j < orderingsize; j++) {
        indicator.push_back(_ordering._cliquenodeindex[j]);
        curr_value = ORACLE->GetCost(indicator, _cliqueIndex);
        res[_ordering._cliquenodeindex[j]] = curr_value - prev_value; 
        prev_value = curr_value;
    }
    return res;
}

CliqueBase::CliqueBase(int clique_index, Clique &clique, Oracle* ORACLE) {
    Ordering ordering;
    for (int i = 0; i < clique.size(); i++) {
        ordering._cliquenodeindex.push_back(i);
        ordering._globalnodeindex.push_back(clique[i]);
        _invalidEbCombination.push_back(0);
    }

    ExtremeBase *eb = new ExtremeBase(clique_index, ordering, ORACLE);
    eb->_lambda = 1;
    _extremeBases.push_back(eb);
}

CliqueBase::~CliqueBase() {
	for (auto it = _extremeBases.begin(); it != _extremeBases.end(); ++it) {
		delete *it;
	}
	_extremeBases.clear();
}

std::vector<double> CliqueBase::GetConvexCombination() {
    int clique_size = (*_extremeBases.begin())->_coordinateValues.size();

    std::vector<double> out(clique_size);
    for (int i = 0; i < clique_size; ++i)
        out[i] = 0;

    int num_extreme_bases = (int) _extremeBases.size();
    for (auto it = _extremeBases.begin(); it != _extremeBases.end(); ++it) {
        ExtremeBase *eb = *it;
        for (int j = 0; j < clique_size; ++j) {
            out[j] += eb->_lambda * eb->_coordinateValues[j];
        }
    }
    return out;
}


SoSBase::SoSBase(int num_nodes, int numLabel, std::vector<Clique> &cliques, Oracle* oracle)
        : _numPixel(num_nodes), _cliques(cliques), _numLabel(numLabel), ORACLE(oracle) {
    _coordinates = new double[num_nodes];
    memset(_coordinates, 0, num_nodes * sizeof(double));

    int size = (int) cliques.size();
    for (int i = 0; i < size; ++i) {
        Clique &clique = cliques[i];
        CliqueBase *b = new CliqueBase(i, clique, ORACLE);
        _bases.push_back(b);
        std::vector<double> clique_coordinates = b->GetConvexCombination();
        int clique_size = (int) clique.size();
        for (int j = 0; j < clique_size; ++j) {
            int node_index = clique[j];
            _coordinates[node_index] += clique_coordinates[j];
        }

    }
}

SoSBase::~SoSBase() {
	for (std::vector<CliqueBase*>::iterator it = _bases.begin(); it != _bases.end(); ++it)
	{
		delete (*it);
	}
	_bases.clear();
	_cliques.clear();
	delete[] _coordinates;
}

void SoSBase::AddUnary(void) {
    int _imgHeight = ORACLE->_imgHieght;
    int _imgWidth = ORACLE->_imgWidth;
    int _imgLabel = ORACLE->_imgLabel; 
    for (int j = 0; j < _imgHeight; j++) {
        for (int i = 0; i < _imgWidth; i++) {
            for (int k = 1; k < (_imgLabel); k++) {

                int node_index = (j * _imgWidth + i) * (_imgLabel - 1) + (k - 1);
                double unaryCost = ORACLE->Unarycost(j, i, k);
                _coordinates[node_index] += unaryCost; 
            }
        }
    }

}


double SoSBase::DotProd(std::vector<double> A, std::vector<double> B) {
    double res = 0;
    for (int i = 0; i < A.size(); i++) {
        res += A[i] * B[i];
    }
    return res;
}

std::vector<double> SoSBase::Scalevector(std::vector<double> A, double scalar) {
    for (int i = 0; i < A.size(); i++) {
        A[i] *= scalar;
    }

    return A;
}

std::vector<double> SoSBase::Differencofvectors(std::vector<double> A, std::vector<double> B) {
    for (int i = 0; i < A.size(); i++) {
        A[i] -= B[i];
    }
    return A;
}


std::vector<double> SoSBase::Addofvectors(std::vector<double> A, std::vector<double> B) {
    for (int i = 0; i < A.size(); i++) {
        A[i] += B[i];
    }
    return A;
}

std::vector<double> SoSBase::Computexminus(int cliqueindex) {
    std::vector<double> xminus;
    std::vector<double> x_clique = _bases[cliqueindex]->GetConvexCombination();
	for (int i = 0; i < x_clique.size(); i++)
		xminus.push_back(
			_coordinates[_cliques[cliqueindex][i]] - x_clique[i] -_bases[cliqueindex]->_invalidEbCombination[i]);
    return xminus;
}

std::vector<double> SoSBase::X_clique(int cliqueindex) {
    std::vector<double> x_c_total;
    for (int i = 0; i < _cliques[cliqueindex].size(); i++)
        x_c_total.push_back(_coordinates[_cliques[cliqueindex][i]]);
    return x_c_total;
} 


void printvector(vector<double> v) { 
    int flag = 0;
    for (int i = 0; i < v.size(); i++) {
        cout << v[i] << " ";
    }
    cout << endl; 
}

// this comparator ensures that for 2 equal elements, element with lower index comes first...to ensure that it is legal.
bool comp(const pair<double, int> &a, const pair<double, int> &b) {
    if (fabs(a.first - b.first) < 1e-10)
        return a.second < b.second;
    return a.first < b.first;

}

int SoSBase::Addordering(int cliqueIndex, bool minnorm_only_flag, int &num_of_legal_bases, int &num_of_illegal_bases)  
{ 
    vector <pair<double, int>> greedysort;
    Ordering ordering;
    int temp_cliquesize = _cliques[cliqueIndex].size();
    std::vector<double> x_c_total = X_clique(cliqueIndex); 
    for (int i = 0; i < temp_cliquesize; i++) {
        greedysort.push_back(make_pair(x_c_total[i], i));
    }

    sort(greedysort.begin(), greedysort.end(), comp);

    for (int i = 0; i < _cliques[cliqueIndex].size(); i++) {
        ordering._cliquenodeindex.push_back(greedysort[i].second);
        ordering._globalnodeindex.push_back(_cliques[cliqueIndex][greedysort[i].second]);
    }

    // following for loop checks if the ordering is illegal

    bool illegal_flag = 0;
    map<int, int> mp;
    for (int i = 0; i < greedysort.size(); i++) {
        int p = greedysort[i].second;
        int lb = p % (_numLabel - 1);
        p = p / (_numLabel - 1);
        if (mp[p] > lb) {
            illegal_flag = 1;
        }
        mp[p] = max(mp[p], lb);
    }

    if(illegal_flag){
        num_of_illegal_bases++;
    }
    else{
        num_of_legal_bases++;
    } 

    if(illegal_flag && !minnorm_only_flag){
        return 3;
    }


    ExtremeBase *newebase = new ExtremeBase(cliqueIndex, ordering, ORACLE);
    newebase->_lambda = 0;
    std::vector<double> newebasetranslated = Addofvectors(newebase->_coordinateValues, Computexminus(cliqueIndex)); 
    if (DotProd(x_c_total, newebasetranslated) + EPSILON2 >= DotProd(x_c_total, x_c_total))
        return 0;

    for (int i = 0; i < _bases[cliqueIndex]->_extremeBases.size(); i++) {
        if (_bases[cliqueIndex]->_extremeBases[i]->_coordinateValues == newebase->_coordinateValues)
            return 0;  // changed
    }
    _bases[cliqueIndex]->_extremeBases.push_back(newebase);
    return 1;

}

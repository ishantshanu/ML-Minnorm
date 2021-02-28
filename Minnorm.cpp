#include "stdafx.h"
#include "Minnorm.h"
#include <fstream>
#include "TestUtil.h"
#include <ctime>
Minnorm::Minnorm(int numPixels, int numLabel, std::vector<Clique> cliques, Oracle* ORACLE)
        : _numPixel(numPixels), _numLabel(numLabel) {

    _numCliques = cliques.size();
    _sosbase = new SoSBase(numPixels, numLabel, cliques, ORACLE);
    _sosbase->AddUnary();

    clique_pixel.resize(cliques.size());    // storing corresponding pixels of extended cliques in clique_pixel
    for (int i = 0; i < int(cliques.size()); i++) {
        for (auto w:cliques[i]) {
            int p = w / (numLabel - 1);
            clique_pixel[i].insert(p);
        }
    }
    //cout << "RUN>>>>>>>>>>";
}

Minnorm::~Minnorm()
{
	//cout << "Inside minnorm Destructor" << endl;
	clique_pixel.clear();
	delete _sosbase;
	//cout << "Going Out of minnorm Destructor" << endl;
}

/*

S0 = S(:,2:end)-S(:,ones(1,size(S,2)-1)); % subspace after translating by S(:,1)

y = S(:,1)- S0*((S0'*S0)\(S0'*S(:,1))); %now y is min norm

*/

MatrixXd Minnorm::computeminnormpoint(MatrixXd Bases) {
    MatrixXd translatedbases(Bases.rows(), Bases.cols() - 1);
    for (int i = 0; i < Bases.rows(); i++) {
        for (int j = 0; j < Bases.cols() - 1; j++) {
            translatedbases(i, j) = Bases(i, j + 1) - Bases(i, 0);
        }
    }
    MatrixXd firstcolumn(Bases.rows(), 1);
    for (int i = 0; i < Bases.rows(); i++) {
        firstcolumn(i, 0) = Bases(i, 0);
    }
    HouseholderQR<MatrixXd> qr(translatedbases.transpose() * translatedbases);
    MatrixXd temp_matrix = (translatedbases.transpose() * translatedbases).fullPivLu().solve(
            translatedbases.transpose() * firstcolumn);
    MatrixXd temp_matrix2 = translatedbases * temp_matrix;
    for (int i = 0; i < Bases.rows(); i++) {
        firstcolumn(i, 0) -= temp_matrix2(i, 0);
    }
    return firstcolumn;
}

std::unordered_set<NodeIndex> Minnorm::DoInference(double &dual, string name) {
    _initialEnergy = 0;
    for (int i = 0; i < _numPixel; i++) {
		//cout << _sosbase->_coordinates[i] << " ";
        if (_sosbase->_coordinates[i] < 0)
            _initialEnergy += _sosbase->_coordinates[i];
    }
    std::cout << "init E = " << _initialEnergy << endl;
    //getchar();
    int terminatecount = 0;
    int cliqueIndex = -1;
    double tempdotprod = DBL_MAX;
    vector<double> prev_normclique(_numCliques);
    for (int i = 0; i < _numCliques; i++)
        prev_normclique[i] = 100000000000000;
    int num_of_iteration = 0;
    int num_of_legal_bases = 0;
    int num_of_illegal_bases = 0;
	bool only_minnorm_flag = 0; 
	clock_t time_req;
	time_req = clock();
    while (1) {
        if (terminatecount == _numCliques)
            break;
        cliqueIndex++;
        if (cliqueIndex == _numCliques) {
            cliqueIndex = 0;
            terminatecount = 0;
        }
        double prevDotprod = DBL_MAX;
        double tempdotprodm = DBL_MAX;
        int countchangeincliques = 0;
        int mazorcount = 0;
        double prev_dual = 0;
        for (int i = 0; i < _sosbase->_numPixel; i++) {
            if (_sosbase->_coordinates[i] < -EPSILON) {
                prev_dual += _sosbase->_coordinates[i];
            }
        }
		while (1) //MAJOR cycle
		{
			mazorcount++;
			std::vector<double> xCorrespondingToClique = _sosbase->X_clique(cliqueIndex);
			double prev_norm = _sosbase->DotProd(xCorrespondingToClique, xCorrespondingToClique);
			std::vector<double> x_minus = _sosbase->Computexminus(cliqueIndex);
            int flag = _sosbase->Addordering(cliqueIndex, only_minnorm_flag, num_of_legal_bases, num_of_illegal_bases);
            num_of_iteration++;
            double current_primal = 0;
            int cliquesize = _sosbase->_cliques[cliqueIndex].size();
            if (flag == 3) {
                // DP Implementation
				vector<double> temp_xc;
				for (int i = 0; i < _sosbase->_cliques[cliqueIndex].size(); i++) {
					_sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] -= _sosbase->_bases[cliqueIndex]->_invalidEbCombination[i];
					temp_xc.push_back(_sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]]);
					_sosbase->_bases[cliqueIndex]->_invalidEbCombination[i] = 0;
				}
                for (auto w:clique_pixel[cliqueIndex]) {
                    for (int i = _numLabel - 2; i >= 0; i--) {
                        double sum = _sosbase->_coordinates[w * (_numLabel - 1) + i];
                        double denom = 1;
                        int idx = i;
                        for (int j = i + 1; j < (_numLabel - 1); j++) {
                            if (sum / denom < _sosbase->_coordinates[w * (_numLabel - 1) + j])    // termination
                            {
                                break;
                            }
                            idx = j;
                            sum += _sosbase->_coordinates[w * (_numLabel - 1) + j];
                            denom++;
                        }
                        for (int k = i; k <= idx; k++)
                            _sosbase->_coordinates[w * (_numLabel - 1) + k] = sum / denom;
                    }
                }
				for (int i = 0; i < _sosbase->_cliques[cliqueIndex].size(); i++) {
					_sosbase->_bases[cliqueIndex]->_invalidEbCombination[i] = _sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] - temp_xc[i];
				}
				temp_xc.clear();
				countchangeincliques++;
                continue;
            }// end of flag 3

            if (!flag) {
				if (!countchangeincliques ) {
					terminatecount++;
				}
				prev_normclique[cliqueIndex] = prev_norm;
                break;
            }

            int minorcount = 0;
            while (1) //MINOR CYCLE // *************************************
            {
                minorcount++;
                if (minorcount > 50) { 
                    break;
                }
                // "minor"
                countchangeincliques++;
                int numofbasesinclique = _sosbase->_bases[cliqueIndex]->_extremeBases.size();
				if (numofbasesinclique < 2)
                    break;
                MatrixXd B(cliquesize, numofbasesinclique);
                for (int i = 0; i < cliquesize; i++)
                    for (int j = 0; j < numofbasesinclique; j++) {
                        B(i, j) = _sosbase->_bases[cliqueIndex]->_extremeBases[j]->_coordinateValues[i] + x_minus[i]+ _sosbase->_bases[cliqueIndex]->_invalidEbCombination[i];
                    }
                MatrixXd affineMinimizer = computeminnormpoint(B);
                affineMinimizer.conservativeResize(affineMinimizer.rows() + 1, affineMinimizer.cols());
                affineMinimizer(cliquesize, 0) = 1;
                B.conservativeResize(B.rows() + 1, B.cols());
                for (int i = 0; i < numofbasesinclique; i++) {
                    B(B.rows() - 1, i) = 1;
                }
                MatrixXd alpha = B.fullPivLu().solve(affineMinimizer);
                MatrixXd xRegTempVec(cliquesize + 1, 1);
                for (int i = 0; i < cliquesize; i++)
                    xRegTempVec(i, 0) = xCorrespondingToClique[i];
                xRegTempVec(cliquesize, 0) = 1;
                int flag_alpha = 0;
                for (int i = 0; i < numofbasesinclique; i++) {
                    if (alpha(i, 0) < 0.0)
                        flag_alpha = 1;
                }
                if (flag_alpha == 0) {
                    vector<int> erase_list1;
                    std::vector<double> lambda;
                    for (int j = 0; j < numofbasesinclique; j++) { 
                        if (alpha(j, 0) < EPSILON)
                            erase_list1.push_back(j);
                        else
                            lambda.push_back(alpha(j, 0));
                    }
                    int index_shift1 = 0;
                    for (int i = 0; i < erase_list1.size(); i++) {
						ExtremeBase *del = _sosbase->_bases[cliqueIndex]->_extremeBases.at(erase_list1[i] - index_shift1);
						delete del;
                        _sosbase->_bases[cliqueIndex]->_extremeBases.erase(
                                _sosbase->_bases[cliqueIndex]->_extremeBases.begin() + (erase_list1[i] - index_shift1));
                        index_shift1++;
                    }
                    numofbasesinclique = _sosbase->_bases[cliqueIndex]->_extremeBases.size();
                    for (int i = 0; i < cliquesize; i++) {
                        _sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] = affineMinimizer(i, 0); 
                    } 
                    for (int i = 0; i < numofbasesinclique; i++) {
                        _sosbase->_bases[cliqueIndex]->_extremeBases[i]->_lambda = lambda[i];
                    }
					erase_list1.clear();
					lambda.clear();
                    break;
                }
                double beta = 2.0;
                double lambda;
                for (int i = 0; i < numofbasesinclique; i++) {
                    lambda = _sosbase->_bases[cliqueIndex]->_extremeBases[i]->_lambda;
                    if (alpha(i, 0) < 0 && lambda != alpha(i, 0)) {
                        if (lambda / (lambda - alpha(i, 0)) >= 0.0) {
                            beta = min(beta, lambda / (lambda - alpha(i, 0)));
                        }
                    }
                } 
                if (beta == 2)
                    beta = 0; 
				MatrixXd projected_x; 
				projected_x = beta * affineMinimizer + (1 - beta) * xRegTempVec; 
				for (int i = 0; i < cliquesize; i++) { 
                    xCorrespondingToClique[i] = projected_x(i, 0);
                    _sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] = projected_x(i, 0);
                }
                vector<int> erase_list;
                vector<double> new_lambda;
                for (int i = 0; i < numofbasesinclique; i++) {
                    lambda = _sosbase->_bases[cliqueIndex]->_extremeBases[i]->_lambda;
                    if ((lambda * (1 - beta) + alpha(i, 0) * beta) <= EPSILON) {
                        erase_list.push_back(i);
                    } else {
                        new_lambda.push_back((lambda * (1 - beta) + alpha(i, 0) * beta));
                    }
                }
                //cout<<endl;
                int index_shift = 0;
                for (int i = 0; i < erase_list.size(); i++) {
					ExtremeBase *del = _sosbase->_bases[cliqueIndex]->_extremeBases.at(erase_list[i] - index_shift);
					delete del;
                    _sosbase->_bases[cliqueIndex]->_extremeBases.erase(
                            _sosbase->_bases[cliqueIndex]->_extremeBases.begin() + (erase_list[i] - index_shift));
                    index_shift++;
                }
				
				for (int i = 0; i < new_lambda.size(); i++) {
					_sosbase->_bases[cliqueIndex]->_extremeBases[i]->_lambda = new_lambda[i];
				}
				//deallocate memory
				erase_list.clear();
				new_lambda.clear();
				//delete projected_x;

            }// end of minor cycle

			//TestUtil testutil;
			std::unordered_set<NodeIndex> minimizer_set_tmp;
			prev_dual = 0;
			for (int i = 0; i < _sosbase->_numPixel; i++) {
				if (_sosbase->_coordinates[i] < -EPSILON) {
					minimizer_set_tmp.insert(i);
					prev_dual += _sosbase->_coordinates[i];
				}
			}
		
			if (flag == 2) {
                if (prev_norm == _sosbase->DotProd(xCorrespondingToClique, xCorrespondingToClique)) {
					if (!countchangeincliques) {
						terminatecount++;
					}
                    break;
                }
            }
            // put the value for current iteration in a file
            // break if dual did not change
            double current_dual = 0;
            for (int i = 0; i < _sosbase->_numPixel; i++) {
                if (_sosbase->_coordinates[i] < -EPSILON) {
                    current_dual += _sosbase->_coordinates[i];
                }
            }
            prev_dual = current_dual;
			if (!only_minnorm_flag) {
				// DP Implementation
				vector<double> temp_xc;
				for (int i = 0; i < _sosbase->_cliques[cliqueIndex].size(); i++) {
					_sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] -= _sosbase->_bases[cliqueIndex]->_invalidEbCombination[i];
					temp_xc.push_back(_sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]]);
					_sosbase->_bases[cliqueIndex]->_invalidEbCombination[i] = 0;
				}

				for (auto w : clique_pixel[cliqueIndex]) {
					for (int i = _numLabel - 2; i >= 0; i--) {
						double sum = _sosbase->_coordinates[w * (_numLabel - 1) + i];
						double denom = 1;
						int idx = i;
						for (int j = i + 1; j < (_numLabel - 1); j++) {
							if (sum / denom < _sosbase->_coordinates[w * (_numLabel - 1) + j])    // termination
							{
								break;
							}
							idx = j;
							sum += _sosbase->_coordinates[w * (_numLabel - 1) + j];
							denom++;
						}
						for (int k = i; k <= idx; k++)
							_sosbase->_coordinates[w * (_numLabel - 1) + k] = sum / denom;
					}
				}
				for (int i = 0; i < _sosbase->_cliques[cliqueIndex].size(); i++) {
					_sosbase->_bases[cliqueIndex]->_invalidEbCombination[i] = _sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] - temp_xc[i];
				}
				temp_xc.clear();
			}
		xCorrespondingToClique.clear();
		x_minus.clear();
        }// end of major cycle
    }// end of outter loop
	
	std::ofstream outfilex;
	outfilex.open("xvec/xvec"+name+".txt", std::ios_base::trunc);
	for (int i = 0; i < _sosbase->_numPixel; i++) {
		outfilex << _sosbase->_coordinates[i] << " \n";
		//cout << _sosbase->_coordinates[i] << " ";
	}
    cout << "FINAL X   ";
    std::unordered_set<NodeIndex> minimizer_set;
    //cout << endl;
    for (int i = 0; i < _sosbase->_numPixel; i++) {
       // cout << _sosbase->_coordinates[i] << "  ";
        if (_sosbase->_coordinates[i] < -EPSILON) {
           // cout << " " << i;
            dual += _sosbase->_coordinates[i];
            minimizer_set.insert(i);
        }
    }
	prev_normclique.clear();

	return minimizer_set;
}

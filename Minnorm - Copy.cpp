#include "stdafx.h"
#include "Minnorm.h"
#include <fstream>
#include "TestUtil.h"
#include <ctime>
Minnorm::Minnorm(int numPixels, int numLabel, std::vector<Clique> cliques)
        : _numPixel(numPixels), _numLabel(numLabel) {

    _numCliques = cliques.size();
    _sosbase = new SoSBase(numPixels, numLabel, cliques);
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
   // cout << "num of cliques " << _numCliques << endl;
	bool only_minnorm_flag = 0;

	std::ofstream outfile;
	if(only_minnorm_flag)
		outfile.open("stats_minnorm.txt", std::ios_base::trunc);
	else
	outfile.open("stats_ours.txt", std::ios_base::trunc);

	clock_t time_req;

	// Using pow function
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
        //cout <<"clique index  "<< cliqueIndex << endl;

		//cout << _sosbase->_cliques[cliqueIndex].size() << endl;
        // coordinate value update
        //std::vector<double> y_clique = _sosbase->_bases[cliqueIndex]->GetConvexCombination();
        //for (int i = 0; i<_sosbase->_cliques[cliqueIndex].size(); i++)
        //	_sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] -= y_clique[i];
        ////for(int i=0;i< _sosbase->_bases[cliqueIndex]->_extremeBases.size();i++){
        ////    ExtremeBase *del=_sosbase->_bases[cliqueIndex]->_extremeBases.at(i);
        ////    delete del;

        //for (auto it = _sosbase->_bases[cliqueIndex]->_extremeBases.begin(); it != _sosbase->_bases[cliqueIndex]->_extremeBases.end(); ++it){
        //	delete *it;
        //}
        //_sosbase->_bases[cliqueIndex]->_extremeBases.clear();
        //Ordering ordering;
        //for (int i = 0; i <_sosbase->_cliques[cliqueIndex].size(); i++){
        //	ordering._cliquenodeindex.push_back(i);
        //	ordering._globalnodeindex.push_back(_sosbase->_cliques[cliqueIndex][i]);
        //}
        //ExtremeBase* newebase = new ExtremeBase(cliqueIndex, ordering);
        //newebase->_lambda = 1;
        //_sosbase->_bases[cliqueIndex]->_extremeBases.push_back(newebase);
        //y_clique = _sosbase->_bases[cliqueIndex]->GetConvexCombination();
        //for (int i = 0; i<_sosbase->_cliques[cliqueIndex].size(); i++)
        //	_sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] += y_clique[i];
        // coordinate value update
        //cout<<"Ebs removed \n";
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
			/*	if (mazorcount > 50) {
				  terminatecount++;
				  cout <<cliqueIndex<<"  "<< "yes1 \n";
				  break;
			  }*/

			  //cout << "terminate count: " << terminatecount << endl;
			  //cout << "clique index: " << cliqueIndex << endl;
			  //	cout << "Mazor \n";
			  //	cout <<cliqueIndex << "   MAZOR \n";
			  /*if (countchangeincliques >2)
				  break;*/
				  //cout << "before Addordering \n";
				  //for (int j = 0; j < _sosbase->_bases[cliqueIndex]->_extremeBases.size(); j++) {
				  //	cout << _sosbase->_bases[cliqueIndex]->_extremeBases[j]->_lambda << "       ";
				  //	for (int i = 0; i < _sosbase->_cliques[cliqueIndex].size(); i++) {
				  //		cout << /*_sosbase->_bases[cliqueIndex]->_extremeBases[j]->_lambda * */_sosbase->_bases[cliqueIndex]->_extremeBases[j]->_coordinateValues[i] << " ";
				  //	}
				  //	cout << endl;
				  //}
				  //cout << "X values \n";
				  //for (int i = 0; i < _sosbase->_numPixel; i++) {
				  //	cout << _sosbase->_coordinates[i] << "  ";
				  //}
				  //cout << endl;

			std::vector<double> xCorrespondingToClique = _sosbase->X_clique(cliqueIndex);
			double prev_norm = _sosbase->DotProd(xCorrespondingToClique, xCorrespondingToClique);
			std::vector<double> x_minus = _sosbase->Computexminus(cliqueIndex);
			// parameter #2 == 1 for minnorm only
			// parameter #2 == 0 for our algo

			////cout << "before Addordering \n";
			//for (int j = 0; j < _sosbase->_bases[cliqueIndex]->_extremeBases.size(); j++) {
			//	//cout << _sosbase->_bases[cliqueIndex]->_extremeBases[j]->_lambda << "       ";
			//	for (int i = 0; i < _sosbase->_cliques[cliqueIndex].size(); i++) {
			//		cout << _sosbase->_bases[cliqueIndex]->_extremeBases[j]->_lambda * _sosbase->_bases[cliqueIndex]->_extremeBases[j]->_coordinateValues[i] << " ";
			//	}
			//	cout << endl;
			//}
			/*for (int i = 0; i < _sosbase->_bases[cliqueIndex]->_invalidEbCombination.size(); i++) {
				cout << _sosbase->_bases[cliqueIndex]->_invalidEbCombination[i] << " ";
			}
			cout << endl;*/

            int flag = _sosbase->Addordering(cliqueIndex, only_minnorm_flag, num_of_legal_bases, num_of_illegal_bases);
           // cout << flag << " = = flag \n";

            // write values in the stats file for plot
            num_of_iteration++;
           
            double current_primal = 0;

            
//            cout << "X values \n";
           /* for (int i = 0; i < _sosbase->_numPixel; i++) {
            	cout << _sosbase->_coordinates[i] << "  ";
            }
            cout << endl;*/
            int cliquesize = _sosbase->_cliques[cliqueIndex].size();
			TestUtil testutil;
            // ordering was illegal... call GC / DP
            if (flag == 3) {
                // DP Implementation
				vector<double> temp_xc;
				//cout << "calling Dp\n";
			/*	for (int i = 0; i < _sosbase->_cliques[cliqueIndex].size(); i++) {
					cout << _sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] << " ";
				}
				cout << endl;*/
				for (int i = 0; i < _sosbase->_cliques[cliqueIndex].size(); i++) {
					_sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] -= _sosbase->_bases[cliqueIndex]->_invalidEbCombination[i];
					temp_xc.push_back(_sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]]);
					_sosbase->_bases[cliqueIndex]->_invalidEbCombination[i] = 0;
					//cout << _sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] << " ";
				}
				//cout << endl;
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
					//cout << _sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] << " ";
				}
				//cout << endl;
				/*for (int i = 0; i < _sosbase->_cliques[cliqueIndex].size(); i++) {
					cout << _sosbase->_bases[cliqueIndex]->_invalidEbCombination[i] << " ";
				}
				cout << endl;*/

				//std::unordered_set<NodeIndex> minimizer_set_tmp;
				//for (int i = 0; i < _sosbase->_numPixel; i++) {
				//	if (_sosbase->_coordinates[i] < -EPSILON) {
				//		minimizer_set_tmp.insert(i);
				//		//cout << i << "  ";
				//	}
				//}
				////cout << endl;
				//prev_dual = 0;
				//for (int i = 0; i < _sosbase->_numPixel; i++) {
				//	if (_sosbase->_coordinates[i] < -EPSILON) {
				//		prev_dual += _sosbase->_coordinates[i];
				//	}
				//}
				//testutil.GetSolutionAsImage(minimizer_set_tmp, prev_dual, current_primal, _sosbase->_cliques);
				////cout << current_primal << "   cp  " << prev_dual << endl;
				countchangeincliques++;
                continue;
            }// end of flag 3

            //cout << " cliqueIndex   flag  " << cliqueIndex <<"  "<<flag<< "\n";
            if (!flag) {
                //cout << prev_normclique[cliqueIndex] << "       prev_normclique[cliqueIndex]     prev_norm        " << prev_norm<<endl;
				if (!countchangeincliques /*&& abs(prev_normclique[cliqueIndex] - prev_norm)<EPSILON2*/) {
					//cout << "terminatecount++1 \n";
					terminatecount++;
				}
                //cout << "terminatecount...............................  " << terminatecount<<endl;
                prev_normclique[cliqueIndex] = prev_norm;
                break;
            }
            //cout << " prev_norm  " << prev_norm << endl;
            //for (int i = 0; i < _sosbase->_bases[cliqueIndex]->_extremeBases.size(); i++) {
            //cout << "  " << _sosbase->_bases[cliqueIndex]->_extremeBases[i]->_lambda << endl;
            //}
			
            int minorcount = 0;
            while (1) //MINOR CYCLE // *************************************
            {
                minorcount++;
                if (minorcount > 50) {
                   // cout << "yes2 \n";
                    break;
                }
                //cout << "minor\n";
                countchangeincliques++;
                //cout << cliqueIndex <<"  "<< _sosbase->_bases[cliqueIndex]->_extremeBases.size() << "   minor \n";
                int numofbasesinclique = _sosbase->_bases[cliqueIndex]->_extremeBases.size();
                //cout << "hmmm";
                if (numofbasesinclique < 2)
                    break;
                MatrixXd B(cliquesize, numofbasesinclique);
                for (int i = 0; i < cliquesize; i++)
                    for (int j = 0; j < numofbasesinclique; j++) {
                        B(i, j) = _sosbase->_bases[cliqueIndex]->_extremeBases[j]->_coordinateValues[i] + x_minus[i]+ _sosbase->_bases[cliqueIndex]->_invalidEbCombination[i];
                    }
                //	cout << "hmmm1";
                MatrixXd affineMinimizer = computeminnormpoint(B);
                affineMinimizer.conservativeResize(affineMinimizer.rows() + 1, affineMinimizer.cols());
                affineMinimizer(cliquesize, 0) = 1;
                B.conservativeResize(B.rows() + 1, B.cols());
                for (int i = 0; i < numofbasesinclique; i++) {
                    B(B.rows() - 1, i) = 1;
                }
                MatrixXd alpha = B.fullPivLu().solve(affineMinimizer);

                //	cout << "hmmm2";
                MatrixXd xRegTempVec(cliquesize + 1, 1);
                for (int i = 0; i < cliquesize; i++)
                    xRegTempVec(i, 0) = xCorrespondingToClique[i];
                xRegTempVec(cliquesize, 0) = 1;
                int flag_alpha = 0;
              //  cout << "alpha   lambda \n";
                for (int i = 0; i < numofbasesinclique; i++) {
                   // cout << "  " << alpha(i, 0) << "   " << _sosbase->_bases[cliqueIndex]->_extremeBases[i]->_lambda << endl;
                    if (alpha(i, 0) < 0.0)
                        flag_alpha = 1;
                }
                //DEBUG STATEMENT
                /*	if (alpha(numofbasesinclique - 1, 0) < 0){
                        cout << "FOUND BUG.............................................."; getchar();
                    }*/


                if (flag_alpha == 0) {
                    vector<int> erase_list1;
                    std::vector<double> lambda;
                    for (int j = 0; j < numofbasesinclique; j++) {
						//cout << "  " << alpha(j, 0);
                        if (alpha(j, 0) < EPSILON)
                            erase_list1.push_back(j);
                        else
                            lambda.push_back(alpha(j, 0));
                    }
                    int index_shift1 = 0;
                    for (int i = 0; i < erase_list1.size(); i++) {
                        _sosbase->_bases[cliqueIndex]->_extremeBases.erase(
                                _sosbase->_bases[cliqueIndex]->_extremeBases.begin() + (erase_list1[i] - index_shift1));
                        index_shift1++;
                    }
                    numofbasesinclique = _sosbase->_bases[cliqueIndex]->_extremeBases.size();
                    for (int i = 0; i < cliquesize; i++) {
                        _sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] = affineMinimizer(i, 0);
						//cout << _sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] << " \n";
                    }
					//cout << endl;
                    for (int i = 0; i < numofbasesinclique; i++) {
                        _sosbase->_bases[cliqueIndex]->_extremeBases[i]->_lambda = lambda[i];
                    }

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
                	//cout << "beta = " << beta<<endl;
                if (beta == 2)
                    beta = 0;
				//cout << "projected_x ";
                MatrixXd projected_x = beta * affineMinimizer + (1 - beta) * xRegTempVec;
                for (int i = 0; i < cliquesize; i++) {
					//cout << projected_x(i, 0)<<"  ";
                    xCorrespondingToClique[i] = projected_x(i, 0);
                    _sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] = projected_x(i, 0);
                }
                vector<int> erase_list;
                vector<double> new_lambda;
                //cout<< "Erase index    ";
				//cout << "new lambda \n";
				//cout << "\n alpha   lambda \n";
                for (int i = 0; i < numofbasesinclique; i++) {
                    lambda = _sosbase->_bases[cliqueIndex]->_extremeBases[i]->_lambda;
					//cout << alpha(i, 0) <<"  "<< lambda * (1 - beta) + alpha(i, 0) * beta << "  \n";
                    if ((lambda * (1 - beta) + alpha(i, 0) * beta) <= EPSILON) {
                        erase_list.push_back(i);
                        //cout<<i <<"  ";
                    } else {
                        new_lambda.push_back((lambda * (1 - beta) + alpha(i, 0) * beta));
                    }
                }
                //cout<<endl;
                int index_shift = 0;
                for (int i = 0; i < erase_list.size(); i++) {
                    _sosbase->_bases[cliqueIndex]->_extremeBases.erase(
                            _sosbase->_bases[cliqueIndex]->_extremeBases.begin() + (erase_list[i] - index_shift));
                    index_shift++;
                }
				
				for (int i = 0; i < new_lambda.size(); i++) {
					_sosbase->_bases[cliqueIndex]->_extremeBases[i]->_lambda = new_lambda[i];
					//cout << new_lambda[i]<<"    ";
				}
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
			//testutil.GetSolutionAsImage(minimizer_set_tmp, prev_dual, current_primal, _sosbase->_cliques);
			//outfile << num_of_iteration << " " << num_of_legal_bases << " " << num_of_illegal_bases << " "
			//	<< current_primal << " " << prev_dual << endl;
			//outfile << (clock() - time_req)<< " " << num_of_legal_bases << " " << num_of_illegal_bases << " "	<< current_primal << " " << prev_dual << endl;

			

            if (flag == 2) {
                if (prev_norm == _sosbase->DotProd(xCorrespondingToClique, xCorrespondingToClique)) {
					if (!countchangeincliques) {
						//cout<<"terminatecount++2 \n";
						terminatecount++;
					}
                    break;
                    //cout << "flag=2 \n";
                }
                //cout.precision(17);
                //cout << prev_norm << "  " << _sosbase->DotProd(xCorrespondingToClique, xCorrespondingToClique) << "  " << terminatecount << " flag=2 \n";

            }
			//cout << prev_norm << "  " << _sosbase->DotProd(xCorrespondingToClique, xCorrespondingToClique) << "  " << terminatecount << " prev_norm    now_norm \n";
            //cout << "Num of EBs after minor cycle " << _sosbase->_bases[cliqueIndex]->_extremeBases.size() << endl<<endl;

            // put the value for current iteration in a file
            // break if dual did not change
            double current_dual = 0;
            for (int i = 0; i < _sosbase->_numPixel; i++) {
                if (_sosbase->_coordinates[i] < -EPSILON) {
                    current_dual += _sosbase->_coordinates[i];
                }
            }
          //  cout << setprecision(10) << prev_dual << " " << current_dual << endl;
           /* if (abs(prev_dual - current_dual) < EPSILON) {
                terminatecount++;
				cout << "terminatecount++3 \n";
                prev_dual = current_dual;
                break;
            }*/
            prev_dual = current_dual;
			if (!only_minnorm_flag) {
				// DP Implementation
				vector<double> temp_xc;
			//	cout << "calling Dp_up\n";
				/*for (int i = 0; i < _sosbase->_cliques[cliqueIndex].size(); i++) {
					cout << _sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] << " ";
				}
				cout << endl;*/
				for (int i = 0; i < _sosbase->_cliques[cliqueIndex].size(); i++) {
					_sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] -= _sosbase->_bases[cliqueIndex]->_invalidEbCombination[i];
					temp_xc.push_back(_sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]]);
					_sosbase->_bases[cliqueIndex]->_invalidEbCombination[i] = 0;
				//	cout << _sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] << " ";
				}
				//cout << endl;

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
				//	cout << _sosbase->_coordinates[_sosbase->_cliques[cliqueIndex][i]] << " ";
				}
				//cout << endl;
				/*for (int i = 0; i < _sosbase->_cliques[cliqueIndex].size(); i++) {
					cout<< _sosbase->_bases[cliqueIndex]->_invalidEbCombination[i]<<" ";
				}
				cout << endl;*/

			}



        }// end of major cycle
		//cout << "cout << terminatecount  << "  << "    " << terminatecount << endl;
    }// end of outter loop
	  
	 //for (int j = 0; j < _sosbase->_bases[0]->_extremeBases.size(); j++) {
		//cout << _sosbase->_bases[cliqueIndex]->_extremeBases[j]->_lambda << "       ";
		//for (int i = 0; i < _sosbase->_cliques[0].size(); i++) {
		//	cout << /*_sosbase->_bases[cliqueIndex]->_extremeBases[j]->_lambda * */_sosbase->_bases[0]->_extremeBases[j]->_coordinateValues[i] << " ";
		//}
		//cout << endl;
		//}



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
	//cout <<"dual inside  " << dual << endl;
    return minimizer_set;
}

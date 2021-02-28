// #include "stdafx.h"
// #include "GenerateTable.h"
// #include "TestBinarySegmentation.h"
//
// void GenerateTable::generateTable() {
//     TestBinarySegmentation test;
//     int choice = 1;
//     while (choice) {
//         cout << "Enter table number \n";
//         cin >> choice;
//         if (choice > 3) {
//             cout << "invalid choice please choose between 1 and 3 \n";
//             continue;
//         }
//         vector<double> time_result;
//         if (choice == 1) {
//             int cliquetype = 0;
//             vector<int> clsize(2), ImgSize(2);
//             clsize[0] = 2;
//             clsize[1] = 2;
//             for (int i = 4; i <= 12; i += 2) {
//                 ImgSize[0] = i;
//                 ImgSize[1] = i;
//                 time_result.push_back(
//                         test.Test_artificial(cliquetype, i * i / (clsize[0] * clsize[1]), clsize, ImgSize));
//             }
//         } else if (choice == 2) {
//             int cliquetype = 0;
//             vector<int> clsize(2), ImgSize(2);
//             ImgSize[0] = 20;
//             ImgSize[1] = 20;
//             clsize[0] = 2;
//             clsize[1] = 2;
//             time_result.push_back(
//                     test.Test_artificial(cliquetype, (ImgSize[0]) * (ImgSize[1]) / (clsize[0] * clsize[1]), clsize,
//                                          ImgSize));
//             for (int i = 2; i <= 6; i += 2) {
//                 clsize[0] = 4;
//                 clsize[1] = i;
//                 time_result.push_back(
//                         test.Test_artificial(cliquetype, (ImgSize[0]) * (ImgSize[1]) / (clsize[0] * clsize[1]), clsize,
//                                              ImgSize));
//             }
//         } else if (choice == 3) {
//             int cliquetype = 0;
//             vector<int> clsize(2), ImgSize(2);
//             clsize[0] = 4;
//             clsize[1] = 4;
//             for (int i = 10; i <= 40; i += 10) {
//                 ImgSize[0] = i;
//                 ImgSize[1] = i;
//                 time_result.push_back(
//                         test.Test_artificial(cliquetype, (ImgSize[0]) * (ImgSize[1]) / (clsize[0] * clsize[1]), clsize,
//                                              ImgSize));
//             }
//         }
//
//         cout << "Running time of table number " << choice << " is as following \n";
//
//         for (int i = 0; i < time_result.size(); i++) {
//             cout << time_result[i] << " ";
//         }
//         cout << endl;
//     }
// }
//
// void GenerateTable::debuggerfunction() {
//     int cliquetype = 1;
//     TestBinarySegmentation test;
//     vector<int> clsize(2), ImgSize(2);
//     ImgSize[0] = 2;
//     ImgSize[1] = 2;
//     clsize[0] = 2;
//     clsize[1] = 1;
//     int numclique = 1;
// 	for (int i = 0; i < 100000; i++) {
// 		cout << " for (int i = 0; i < 100000; i++)  " << i << endl;
// 		cout << test.Test_artificial(cliquetype, numclique, clsize, ImgSize);
// 	}
// }

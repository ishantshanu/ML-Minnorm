#include "stdafx.h"
#include "TestUtilm.h"
#include "Oracle.h"

cv::Mat TestUtil::GetSolutionAsImage(int64 key, std::unordered_set<int> &minimizer_set, double &dual,
                                     std::vector<Clique> cliques, int num_largeCliques,
                                     std::vector<int> pairwise_cliques) {
    int imgHeight = ORACLE->_imgHieght;
    int imgWidth = ORACLE->_imgWidth;
    int numLabel = ORACLE->_imgLabel;

    cv::Mat image = cv::Mat(imgHeight, imgWidth, CV_8UC3, cv::Scalar(0));
    double min_sum = dual;
    double function_value = 0;

    uchar BLUE[22] = {0, 255, 0, 0, 255, 255, 0, 128, 128, 0, 128, 0, 0, 128, 54, 0, 0, 54, 54, 0, 54, 255};
    uchar GREEN[22] = {0, 0, 255, 0, 255, 0, 255, 128, 0, 128, 0, 128, 0, 128, 0, 54, 0, 54, 0, 0, 54, 255};
    uchar RED[22] = {0, 0, 0, 255, 0, 255, 255, 0, 128, 128, 0, 0, 128, 128, 0, 0, 54, 0, 54, 54, 54, 255};

    // calculating primal
    std::vector<int> support(imgHeight * imgWidth * (numLabel), 0);
    for (auto w : minimizer_set) {
        int la = (w % (numLabel - 1)) + 1;
        w = w / (numLabel - 1);
        support[w] = max(support[w], la);
    }

    // TODO
    // add all the unaries till the max label
    for (int y = 0; y < imgHeight; y++) {
        for (int x = 0; x < imgWidth; x++) {
            for (int k = 1; k < (numLabel); k++) {
                if (minimizer_set.find((y * imgWidth + x) * (numLabel - 1) + k - 1) != minimizer_set.end()) {
                    double unaryCost = ORACLE->Unarycost(y, x, support[y * imgWidth + x]);
                    function_value += unaryCost;
                }
            }
        }
    }
    //std::cout << "unary cost Function value: " << function_value << std::endl;

    // FORMING THE IMAGE
    // for (auto it = minimizer_set.begin(); it != minimizer_set.end(); ++it) {
    //     int y_index = (int)(*it) / ORACLE->_imgHieght;
    //     int x_index = (*it) % ORACLE->_imgWidth;
    //     image.at<uchar>(y_index, x_index) = 255;
    // }
    for (int y = 0; y < imgHeight; y++) {
        for (int x = 0; x < imgWidth; x++) {
            //cout << support[y * imgWidth + x] << endl;
            //getchar();
            uchar blue, green, red;
            blue = BLUE[support[y * imgWidth + x]];
            green = GREEN[support[y * imgWidth + x]];
            red = RED[support[y * imgWidth + x]];
            blue = (255.0 / numLabel) * support[y * imgWidth + x];
            green = blue;
            red = blue;
            image.at<Vec3b>(y, x) = Vec3b(blue, green, red);
        }
    }

    for (int y = 0; y < cliques.size(); y++) {
        std::vector<int> indicator;
        for (unsigned int j = 0; j < cliques[y].size(); j++) {
            if (minimizer_set.find(cliques[y][j]) != minimizer_set.end()) {
                indicator.push_back(j);
            }
        }
        function_value += ORACLE->GetCost(indicator, y);
    }

    std::cout << "" << function_value;
    std::cout << " " << min_sum << std::endl;
    //if (function_value != min_sum){
    //    std::ofstream outfile;
    //    outfile.open("eval.txt", std::ios_base::app);
    //    outfile  << " final_x= " << min_sum << " energy= " << function_value  << std::endl << std::endl;
    //    //cout << "ASSERTION FAIL:\n Key == "<<key;
    //    //getchar();
    //}
    return image;
}

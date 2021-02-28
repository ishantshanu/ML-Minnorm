#pragma once

class TestUtil {
public:
    cv::Mat
    GetSolutionAsImage(int64 key, std::unordered_set<int> &minimizer_set, double &dual, std::vector <Clique> numclique,
                       int num_largeCliques, std::vector<int> pairwise_cliques);
};

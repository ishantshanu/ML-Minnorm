

## Paper 
This is project page for the paper title "An Inference Algorithm for Multi-Label MRF-MAP Problems with Clique Size 100" published in ECCV 2020. The paper is available at link: https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123650256.pdf. The proposed method is an inference algorithm for solving multi-label MRF-MAP problem. Please cite our paper if you find this code useful. 

'''
@inproceedings{shanu2020inference,
  title={An Inference Algorithm for Multi-Label MRF-MAP Problems with Clique Size 100},
  author={Shanu, Ishant and Bharti, Siddhant and Arora, Chetan and Maheshwari, SN},
  booktitle={European Conference on Computer Vision},
  pages={257--274},
  year={2020},
  organization={Springer}
}
'''

### Getting Started

1. Download the code from this repo.
2. Download PASCAL VOC12 dataset from http://host.robots.ox.ac.uk/pascal/VOC/voc2012/#data.

### Data Preprocess

Please run data_preprocess/runslic.py and data_preprocess/get_unary.py for generating and storing the clique and unary information in the required format. 
 
### Prerequisites

This project can be run on Visual Studio 2017 with the following packges installed. Please download VS2017 from https://visualstudio.microsoft.com/vs/older-downloads/. 

### Packages

1. Please download and setup Eigen3 as per information given on the link: https://phylogeny.uconn.edu/tutorial-v2/part-1-ide-project-v2/setting-up-the-eigen-library-v2/. 
2. Similarly download and setup openCV package as per instructions given at https://www.opencv-srf.com/2017/11/install-opencv-with-visual-studio.html. 

import math
import matplotlib.pyplot as plt
import numpy as np
from skimage.segmentation import slic, quickshift, felzenszwalb
from skimage.segmentation import mark_boundaries
from skimage.transform import resize
from skimage.util import img_as_float
from skimage import io
from PIL import Image
import cv2
import os

img_rows = 200
img_cols = 200



def write_in_file(f_in, h, w, segments_slic, cnt_supp):
    x, y = (img_rows, img_cols)
    segments_slic_1 = resize(segments_slic, (x, y), preserve_range='true')

    for x in range(0, h):
        for y in range(0, w):
            segtmp = segments_slic_1[x, y]
            values = str(x) + ' ' + str(y) + ' ' + str(int(math.floor(segtmp + cnt_supp + 0.5))) + '\n'
            # values = str(values)
            f_in.write(values)


def merge_super_pixel(img, segments_slic, X_O, Y_O):
    for outter_loopy in range(len(X_O)):

        X = X_O[outter_loopy]
        Y = Y_O[outter_loopy]

        mnn = 1000000

        for loopy in range(len(X)):
            x = X[loopy]
            y = Y[loopy]
            if segments_slic[x, y] < mnn:
                mnn = segments_slic[x, y]

        for loopy in range(len(X)):
            x = X[loopy]
            y = Y[loopy]
            lab = segments_slic[x, y]

            for a in range(500):
                for b in range(500):
                    if segments_slic[a, b] == lab:
                        segments_slic[a, b] = mnn

    plt.imshow(mark_boundaries(img, segments_slic, ))
    plt.show()


def solve():

	path = r"path to image folder"

	#fname = open('filenames.txt', 'w')
	for a,b,c in os.walk(os.path.join(path)):
		for filename in c:
			readfile = filename
			writefile = filename
			print(filename)
			print(filename[:-4])
			writefile = 'segment_' + str(filename[:-4]) + '.txt'
			
			#fname.write(str(filename[:-4]) + '\n')
			
			img = io.imread(os.path.join(path,readfile))
			# img = cv2.resize(img, (500, 500))
			f = open(os.path.join(path,'segments',writefile), 'w')

			#print writefile

			# layer one

			cnt_supp = 0

			# segments_slic = felzenszwalb(img, scale=50, sigma=3, min_size=10)
			# segments_slic = quickshift(img, kernel_size=10, max_dist=6, ratio=0.2)

			segments_slic = slic(img, n_segments=300, compactness=10, sigma=2, enforce_connectivity=True)

			# following 4 lines is for visualization, comment them if we are segmenting 1000s of images
			print('SLIC number of segments: {}'.format(len(np.unique(segments_slic))))

			# note that vertical is your x co-ordinate
			# note that horizontal is your y co-ordinate

			X = [[]]
			Y = [[]]
			# merge_super_pixel(img, segments_slic, X, Y)

			write_in_file(f, img_rows, img_cols, segments_slic, cnt_supp)

			# # layer two (uncomment if needed)..... add more layers if needed
			cnt_supp = cnt_supp + len(np.unique(segments_slic))
			#
			#
			segments_slic = slic(img, n_segments= 350 , compactness=10, sigma=2, enforce_connectivity=True)
			print('SLIC number of segments: {}'.format(len(np.unique(segments_slic))))
			#plt.imshow(mark_boundaries(img,segments_slic,))
			#plt.show()
			# X = [[]]
			# Y = [[]]
			# merge_super_pixel(img, segments_slic, X, Y)
			#
			write_in_file(f, img_rows, img_cols, segments_slic, cnt_supp)
			#
			#
			
			# # layer two (uncomment if needed)..... add more layers if needed
			cnt_supp = cnt_supp + len(np.unique(segments_slic))
			#
			#
			segments_slic = slic(img, n_segments= 250 , compactness=10, sigma=2, enforce_connectivity=True)
			print('SLIC number of segments: {}'.format(len(np.unique(segments_slic))))
			#plt.imshow(mark_boundaries(img,segments_slic,))
			#plt.show()
			# X = [[]]
			# Y = [[]]
			# merge_super_pixel(img, segments_slic, X, Y)
			#
			write_in_file(f, img_rows, img_cols, segments_slic, cnt_supp)


			
solve()


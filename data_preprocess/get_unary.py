import os

from skimage.io import imread
from numpy import array_equal

def get_label(pixel):
    label = 0
    if array_equal(pixel, [0,0,0]):
        label = 0
    elif array_equal(pixel, [128, 0,0]):
        label = 1
    elif array_equal(pixel, [0, 128,0]):
        label = 2
    elif array_equal(pixel, [128, 128,0]):
        label = 3
    elif array_equal(pixel, [0, 0,128]):
        label = 4
    elif array_equal(pixel, [128, 0,128]):
        label = 5
    elif array_equal(pixel, [0, 128,128]):
        label = 6
    elif array_equal(pixel, [128, 128,128]):
        label = 7
    elif array_equal(pixel, [64, 0,0]):
        label = 8
    elif array_equal(pixel, [192, 0,0]):
        label = 9
    elif array_equal(pixel, [64, 128,0]):
        label = 10
    elif array_equal(pixel, [192, 128,0]):
        label = 11
    elif array_equal(pixel, [64, 0,128]):
        label = 12
    elif array_equal(pixel, [192, 0,128]):
        label = 13
    elif array_equal(pixel, [64, 128,128]):
        label = 14
    elif array_equal(pixel, [192, 128,128]):
        label = 15
    elif array_equal(pixel, [0, 64,0]):
        label = 16
    elif array_equal(pixel, [128, 64,0]):
        label = 17
    elif array_equal(pixel, [0, 192,0]):
        label = 18
    elif array_equal(pixel, [128, 192,0]):
        label = 19
    elif array_equal(pixel, [0, 64,128]):
        label = 20
    return label


def write_unary(root, file):
    image_longname = os.path.join(root,'out', file)
    print(root)
    print(file[:-9])
    unary_longname = os.path.join(os.path.join(root,'unary'), 'unary'+str(file[:-9]) + '.txt')
    img = imread(image_longname)
    height, weight, _ = img.shape
    with open(unary_longname,'w+') as f:
        for i in range(height):
            for j in range(weight):
                f.write(f"{i} {j} {255} {get_label(img[i][j])}\n")


def batch_process(folder_name):
    for root, directory, files in os.walk(os.path.join(folder_name,'out'), topdown=True):
        for file in files:
            if '.png' in file or '.jpg' in file:
                write_unary(folder_name, file)


if __name__ == '__main__':
    folder_name = r'path to deeplabv3 prediction images'
    batch_process(folder_name)

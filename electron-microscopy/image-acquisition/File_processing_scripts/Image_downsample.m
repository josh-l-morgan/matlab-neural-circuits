image=imread('e:/SEM_Users/Richard/wafer2_Sec1_Montage/FijiStitched2_1.tif');
dimage = imresize(image, [2048 2048]);
imwrite(dimage, 'e:/SEM_Users/Richard/wafer2_Sec1_Montage/FijiStitched2_1d.tif');
% test_patchmatch.m
clear; % make;
clc; clear;
disp('------------------------------------');
img = double(imresize(imread('lena.png'),0.25)); 
patchsize = 5; nn = 10;
tic;
cnn = patchmatch(img, [], patchsize, nn, 20, [], 1, [], [], [], 2);
toc;

visualizer(img, cnn, patchsize);



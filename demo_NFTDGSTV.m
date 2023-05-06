%% 基于WSNM的RNIPT算法
tic
clc;
clear;
close all;
%% setup parameters
lambdaL =4;%
C=5;%5
L1=3;%L4 L
saveDir= '..\NFTDGSTV公开\results\1\';     % save patch
imgpath='..\NFTDGSTV公开\data\1\';      % Data input path
%%
imgDir = dir([imgpath '*.bmp']);
len = length(imgDir);
for i=1:len
    picname=[imgpath  num2str(i),'.bmp'];
    I=imread(picname);%
    [m,n]=size(I);
    [~, ~, ch]=size(I);
    if ch==3
        I=rgb2gray(I); 
    end
    D(:,:,i)=I;
end
tenD=double(D);
[n1,n2,n3]=size(tenD);
n_1=max(n1,n2);%n(1)
n_2=min(n1,n2);%n(2)
patch_frames=L1;% temporal slide parameter
patch_num=n3/patch_frames;
%% constrcut image tensor
for l=1:patch_num
    l
    for i=1:patch_frames
        temp(:,:,i)=tenD(:,:,patch_frames*(l-1)+i);
    end           
        T=C*sqrt(n1*n2);
        lambda4 =lambdaL / sqrt(min(n_1*patch_frames));
%% The proposed NFTDGSTV model
[tenB, tenT,change] = NFTDGSTV(temp, lambda4); 
%% recover the target and background image       
 for i=1:patch_frames
     tarImg=tenT(:,:,i);
     backImg=tenB(:,:,i);
     maxv = max(max(double(I)));
     tarImg=double(tarImg);
     backImg=double(backImg);
     E = uint8( mat2gray(tarImg)*maxv );
     A = uint8( mat2gray(backImg)*maxv );
%% save the results
imwrite(E, [saveDir 'target/' imgDir(i+patch_frames*(l-1)).name]);     % Save target image 
imwrite(A, [saveDir 'background/' imgDir(i+patch_frames*(l-1)).name]); % Save background image
end 
end
toc  
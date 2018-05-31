%The core implementation of "Single Image Defogging Based on Illumination
%Decomposition for Visual Maritime Surveillance",

clear all;
close all;
img = 'E:\QXH\picture\海上雾天\Sea Fog\fog\海雾图像数据库\测试集\4.jpg';
I = im2double(imread(img));
[H, W, D] = size(I);
figure,imshow(I);
% 
%%%%%%%%%%the illumination decomposition
% I: input
% lambda: corresponds to the alpha in Eq. (7) of the paper, control the sharpness of fog layer.
% lambda2: corresponds to the beta in Eq. (7) of the paper, control the illumination of G.
% lb: lower bound of the Layer 1,need to be same dimention with input I 
% I0: initialization of Layer 1, default as the input I
% fast: whether use the fast impletment, 1:fast 2:normal
%%%%%%%%%%the third parameter 
[LB, LR] = layer_decom(I, 50, 0.01, zeros(H,W,D)+0.01, I, 2);  
figure,imshow(LB);
figure,imshow(adjust(LR));
%%%%%%%%%%the fog layer defogging  process. We utilized the improved Berman's algorithm in the process.  
gamma = 1.3;
[out_Im, trans_refined,ind] = non_local_dehazing(uint8((LB)*255),LR, gamma);
out_Im = im2double(out_Im);
figure,imshow(trans_refined);
figure,imshow(out_Im);

%%%%%%%%%luminance conpensate process
% out_Im: input
% LR: the glow-shaped illumination layer.
% I: the original fog image.
% ga: the index in Eq.(17), which control the gamma transformation
[LR2,out_Im2] = luminance_com(out_Im,LR,I,2);
out_Im3 = out_Im2  + LR2;    
adj_percent = [0.0, 0.95];   
out_Im3 = adjust(out_Im3,adj_percent);
figure,imshow(out_Im3);

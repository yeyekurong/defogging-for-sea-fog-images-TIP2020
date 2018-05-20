
clear all;
% close all;
img = 'E:\QXH\Guoqiang\论文\实验图片\fig12\147.jpg';
I = im2double(imread(img));
[H, W, D] = size(I);
figure,imshow(I);

%%%%%%%%%%the illumination decomposition
%%%%%%%%%%process最重要的是迭代次数和亮度控制参数,
[LB, LR] = sepglow5(I, 50, zeros(H,W,D)+0.01, I);  
figure,imshow(LB);
figure,imshow(adjust(LR));
%%%%%%%%%%the fog layer defogging  process. We utilized the improved Berman's algorithm in the process.  
%%%%%%%%%alpha control the fog-removal and artifacts in sky regions
gamma = 1.3;
% [out_Im, trans_refined,ind] = non_local_dehazing_52(uint8((LB)*255),16,LR,A);    
[out_Im, trans_refined,ind] = non_local_dehazing(uint8((LB)*255),LR, gamma);
out_Im = im2double(out_Im);
% figure,imshow(ind);
figure,imshow(trans_refined);
figure,imshow(out_Im*2);

%%%%%%%%%luminance conpensate process
% [LR2,out_Im2] = relectance_coff(out_Im,LR,I);
% [LR2,out_Im2] = relectance_coff_gamma(out_Im,LR,I);
% [LR2,out_Im2] = relectance_coff_o(out_Im,LR,I,2);
[LR2,out_Im2] = relectance_coff2(out_Im,LR,I,1);

out_Im3 = out_Im2  + LR2;    
adj_percent = [0.01, 0.99];   
out_Im3 = adjust(out_Im3,adj_percent);

figure,imshow(out_Im3);
otg = [img,'_decom_50_0.001_1.3.jpg'];
otg_t = [img,'_decom_t.jpg'];
otg_LR = [img,'_decom_LR.jpg'];
otg_LB = [img,'_decom_LB.jpg'];
imwrite(out_Im3,otg);
imwrite(trans_refined,otg_t);
imwrite(LR,otg_LR);
imwrite(LB,otg_LB);
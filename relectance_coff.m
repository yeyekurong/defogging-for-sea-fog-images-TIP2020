function [LR3,out_q] = relectance_coff(out_Im,LR,I)
[H, W, D] = size(I);
eps = 10^-4;
L1_max = (max(out_Im,[],3));
pp = minmaxfilt(double(L1_max),round(H*0.02),'max','same');
localmax = fastguidedfilter(max(I,[],3), pp, uint16(H*0.1), eps,4);%
out_q = out_Im;
R = (out_q(:,:,1))./(localmax+~localmax);
G = (out_q(:,:,2))./(localmax+~localmax);
B = (out_q(:,:,3))./(localmax+~localmax);
LR2(:,:,1) = LR(:,:,1).*(R);
LR2(:,:,2) = LR(:,:,2).*(G);
LR2(:,:,3) = LR(:,:,3).*(B);
LR2(LR2>1) = 1;
LR2(LR2<0) = 0;
LR3 = double(power(LR2,(1/1.5)));%%%%%%%%%luminance parameters
figure,imshow(LR3*2);
end
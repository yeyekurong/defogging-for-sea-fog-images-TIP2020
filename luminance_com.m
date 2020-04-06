function [LR3,out_q] = luminance_com(out_Im,LR,I,alpha)

[H, W, D] = size(I);
eps = 10^-4;
L1_max = (max(out_Im,[],3));
L1_min = (min(out_Im,[],3));
L1_max(L1_max>1) = 1;
pp = minmaxfilt(double(L1_max),round(H*0.04),'max','same');
localmax = fastguidedfilter(max(I,[],3), pp, uint16(H*0.2), eps,4);%
ratio1=(L1_min)./(localmax+~localmax);    
ratio1(ratio1<0)= 0;
ratio1(ratio1>1)= 1;
lam = adjust(LR(:,:,1));  
c = max(min((1-(ratio1.*lam)),1),0.2);

out_hsv = rgb2hsv(out_Im);
out_hsv(:,:,3) = fastguidedfilter( out_hsv(:,:,3), out_hsv(:,:,3), uint16(H*0.01), 10^-4, 4);  %%%%%%%滤波器窗口大小设置 默认0.04      %%%%%%%%将去噪，光照补偿，去雾融合到一个过程中来。因为本来就是一个模型过程的降质。默认0.01
out_hsv(:,:,2) = out_hsv(:,:,2).*c;
out_hsv(out_hsv<0)=0;
out_hsv(out_hsv>1) = 1;
out_q = hsv2rgb(out_hsv);
R = (out_q(:,:,1))./(localmax+~localmax);
G = (out_q(:,:,2))./(localmax+~localmax);
B = (out_q(:,:,3))./(localmax+~localmax);
LR2(:,:,1) = LR(:,:,1).*(R);
LR2(:,:,2) = LR(:,:,2).*(G);
LR2(:,:,3) = LR(:,:,3).*(B);
LR2(LR2>1) = 1;
LR2(LR2<0) = 0;
LR3 = double(power(LR2,(1/alpha)));%%%%%%%%%luminance parameters
%figure,imshow(LR3);
end
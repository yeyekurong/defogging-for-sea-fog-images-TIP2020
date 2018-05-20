function [LR3,out_q] = relectance_coff2(out_Im,LR,I,alpha)
[H, W, D] = size(I);
eps = 10^-4;
L1_max = (max(out_Im,[],3));
L1_mean = rgb2gray(out_Im);
L1_min = (min(out_Im,[],3));

L1_mean_g = fastguidedfilter(max(I,[],3), L1_min, uint16(H*0.01), eps,4);  %%%%%%%滤波器窗口大小设置 默认0.04
% L1_mean_g = fastguidedfilter(max(I,[],3), L1_mean, uint16(H*0.04), eps,4);  %%%%%%%滤波器窗口大小设置 默认0.04
pp = minmaxfilt(double(L1_max),round(H*0.04),'max','same');
localmax = fastguidedfilter(max(I,[],3), pp, uint16(H*0.2), eps,4);%
ratio1=(L1_min)./(localmax+~localmax);    
ratio1(ratio1<0)= 0;
ratio1(ratio1>1)= 1;
% c = min(alpha*adjust(1-ratio1),1);
lam = adjust(LR(:,:,1));  
c = min((1-(ratio1.*lam).^alpha),1);


% lam = adjust(LR(:,:,1));  
% c = 1./power(exp(max(lam - 0.5, 0)),2);%%%%%%%%%这个参数如何选择，来控制天空区域的颜色问题， 处理后天空区域颜色变浅，不处理颜色比较深
out_hsv = rgb2hsv(out_Im);
% out_hsv(:,:,3) = qxfilter(out_hsv(:,:,3)*255,0.01,0.05)/255;       %%%%%%%%将去噪，光照补偿，去雾融合到一个过程中来。因为本来就是一个模型过程的降质。默认0.01
% out_hsv(:,:,3) = fastguidedfilter( max(I,[],3), out_hsv(:,:,3), uint16(H*0.01), eps, 4);  %%%%%%%滤波器窗口大小设置 默认0.04      %%%%%%%%将去噪，光照补偿，去雾融合到一个过程中来。因为本来就是一个模型过程的降质。默认0.01
out_hsv(:,:,3) = fastguidedfilter( out_hsv(:,:,3), out_hsv(:,:,3), uint16(H*0.01), 10^-4, 4);  %%%%%%%滤波器窗口大小设置 默认0.04      %%%%%%%%将去噪，光照补偿，去雾融合到一个过程中来。因为本来就是一个模型过程的降质。默认0.01
% out_hsv(:,:,3) = fastguidedfilter( out_hsv(:,:,3), out_hsv(:,:,3), uint16(H*0.01), eps, 4);  %%%%%%%滤波器窗口大小设置 默认0.04      %%%%%%%%将去噪，光照补偿，去雾融合到一个过程中来。因为本来就是一个模型过程的降质。默认0.01
% out_hsv(:,:,3) = wlsFilter(out_hsv(:,:,3),0.01, 2, out_hsv(:,:,3));          %%%%%%%%%%%%默认为1，1.2，高原原中用0.5，1.5更不平滑,这样边缘更锐利，L2应该是平滑的，所以取1.5 ， 1.2
out_s = out_hsv(:,:,2);
out_s = out_s.*c;
% out_s = fastguidedfilter(max(I,[],3), out_s, uint16(H*0.04), 10^-4, 4);  %%%%%%%滤波器窗口大小设置 默认0.01      %%%%%%%%将去噪，光照补偿，去雾融合到一个过程中来。因为本来就是一个模型过程的降质。默认0.01
% out_s = qxfilter(out_s*255,0.01,0.1)/255;       %%%%%%%%将去噪，光照补偿，去雾融合到一个过程中来。因为本来就是一个模型过程的降质。默认0.01
out_hsv(:,:,2) = out_s;
out_q = hsv2rgb(out_hsv);
% out_q = qxfilter(out_q*255,0.01,0.02)/255;       %%%%%%%%将去噪，光照补偿，去雾融合到一个过程中来。因为本来就是一个模型过程的降质。默认0.01
R = (out_q(:,:,1))./(localmax+~localmax);
G = (out_q(:,:,2))./(localmax+~localmax);
B = (out_q(:,:,3))./(localmax+~localmax);
LR2(:,:,1) = LR(:,:,1).*(R);
LR2(:,:,2) = LR(:,:,2).*(G);
LR2(:,:,3) = LR(:,:,3).*(B);
% out_q = (out_q*2+out_Im)/3;
LR2(LR2>1) = 1;
LR2(LR2<0) = 0;
LR3 = double(power(LR2,(1.3)));%%%%%%%%%luminance parameters
figure,imshow(LR3);
end
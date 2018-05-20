clear all;
% close all;
img = 'E:\QXH\picture\海上雾天\Sea Fog\49.jpg';
otg = [img,'_decom.jpg'];
I = im2double(imread(img));
[H, W, D] = size(I);
figure,imshow(I);

%%%%%%%%%%the illumination decomposition
%%%%%%%%%%process最重要的是迭代次数和亮度控制参数,
[LB, LR] = sepglow5(I, 500, zeros(H,W,D)+0.01, I);  
figure,imshow(LB);
figure,imshow(adjust(LR));
%%%%%%%%%%the fog layer defogging  process. We utilized the improved Berman's algorithm in the process.  
%%%%%%%%%alpha control the fog-removal and artifacts in sky regions
[out_Im, trans_refined,ind] = non_local_dehazing_52(uint8((LB)*255),3,LR);    
% figure,imshow(ind);
% figure,imshow(trans_refined);
% figure,imshow(out_Im);

    %out_Im = out_Im*2;       %%%%%%%%%如何利用光源信息叠加，来得到清晰无雾的自然亮度图片
    out_q = qxfilter(out_Im*255,0.01)/255;       %%%%%%%%将去噪，光照补偿，去雾融合到一个过程中来。因为本来就是一个模型过程的降质。默认0.01
    %out_q = out_Im;
    L1_max = (max(out_Im,[],3));
    %     L1_max_q = qxfilter(L1_max*255)/255;       %%%%%%%%将去噪，光照补偿，去雾融合到一个过程中来。因为本来就是一个模型过程的降质。
    L1_max_g = fastguidedfilter(max(I,[],3), L1_max, uint16(H*0.04), eps,4);  %%%%%%%滤波器窗口大小设置 默认0.04
    pp = minmaxfilt(double(L1_max),10,'max','same');
    localmax = fastguidedfilter(max(I,[],3), pp, 60, eps,4);%
    ratio1=L1_max_g./(localmax+~localmax);    
    figure,imshow(ratio1);


    lam = adjust(LR(:,:,1));  
    c = exp(max(lam - 0.5, 0));%%%%%%%%%这个参数如何选择，来控制天空区域的颜色问题， 处理后天空区域颜色变浅，不处理颜色比较深
    %%%%或者是白平衡，一开始进行白平衡，还是最后在进行白平衡。
    R = (out_q(:,:,1))./L1_max_g;
    G = (out_q(:,:,2))./L1_max_g;
    B = (out_q(:,:,3))./L1_max_g;
    R(R>1) = 1;
    R(R<0) = 0;
    G(G>1) = 1;
    G(G<0) = 0;
    B(B>1) = 1;
    B(B<0) = 0;
    R = R.^(1./(c.^3));
    G = G.^(1./(c.^3));
    B = B.^(1./(c.^3));

    ratio1_r = ratio1.*R;
    ratio1_g = ratio1.*G;
    ratio1_b = ratio1.*B;


    LR2(:,:,1) = LR(:,:,1).*(ratio1_r);
    LR2(:,:,2) = LR(:,:,2).*(ratio1_g);
    LR2(:,:,3) = LR(:,:,3).*(ratio1_b);
    LR2(LR2>1) = 1;
    LR2(LR2<0) = 0;
    LR3 = double(LR2.^(1.5));%%%%%%%%%参数
    figure,imshow(LR3);
%     imwrite(LR3,'compensate_light.bmp');
    out_Im3 = out_Im*2  + LR3;             

    adj_percent = [0, 0.995];
    out_Im3 = adjust(out_Im3,adj_percent);
    figure,imshow(out_Im3);
    imwrite(out_Im3,otg);

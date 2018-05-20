
function [L1, L2] = sepglow5(I,lambda,lb,hb,L1_0)
% Layer Separation using Relative Smoothness (general)
% [L1 L2] = septRelSmo(I,lambda,lb,hb,I0)
% I: input
% lambda: smoothness regularization parameter on the smoother layer Layer 2 
% lb: lower bound of the Layer 1,need to be same dimention with input I 
% hb: upper bound of the Layer 1,need to be same dimention with input I
% I0: initialization of Layer 1, default as the input I
pad_size = 64;
I = padarray(I,[pad_size,pad_size],'replicate','both');
lb = padarray(lb,[pad_size,pad_size],'replicate','both');
hb = padarray(hb,[pad_size,pad_size],'replicate','both');

lambda2 = 0.0001;   %%%%%%0.001 浓雾要大一些, 控制L2亮度大小的,越大，L2月亮 L2大，则天空部分透射率高，则LR差距大，则图像去雾更亮，但是过大会导致边缘artifacts
if ~exist('I0','var')
    L1_0 = I;
end
[N,M,D] = size(I);

% filters
f1 = [1, -1];
f2 = [1; -1];
f3 = [0, -1, 0; 
      -1, 4, -1;
      0, -1, 0];
f3 = [-1 -1 -1;
      -1 8 -1;
      -1 -1 -1];
f4 = [1 1 1;
      1 1 1;
      1 1 1];
f5 = fspecial('gaussian',10,3);
f4 = f4/sum(f4(:));
sizeI2D = [N,M];
otfFx = psf2otf(f1,sizeI2D);
otfFy = psf2otf(f2,sizeI2D);
otfL = psf2otf(f3,sizeI2D);
otfI = psf2otf(f5,sizeI2D);
% Weight = CalWeightFun(I);
% Weight = Weight.*Weight;
Weight = ones(sizeI2D);
Normin1 = repmat(abs(otfL),[1,1,D]).^2.*fft2(I).*repmat(Weight,[1,1,D]);
Denormin1 = abs(otfL).^2.*Weight;
Denormin2 = abs(otfFx).^2 + abs(otfFy ).^2;
intensity = lambda2*abs(otfI).^2;

if D>1
    Denormin1 = repmat(Denormin1,[1,1,D]);
    Denormin2 = repmat(Denormin2,[1,1,D]);
    intensity = repmat(intensity,[1,1,D]);
end

eps = 1e-16;
L1 = L1_0;

thr = 0.05;   %%%%%默认0.05
beta = 20;
for i = 1:8  %%%%%%%迭代次数的加大会导致L1高频部分更亮一些 8
    Denormin   = lambda*Denormin1 + beta*Denormin1;
    %% update g
    gFx = -imfilter(L1,f3,'circular');  %%相当imfilter(L1,[-1,1],'circular'); 
    gFx = max(abs(gFx) - 1/ beta, 0) .* sign(gFx); 
    Normin2 = imfilter(gFx,f3,'circular');  %%相当imfilter(L1,[-1,1],'circular'); 

    FL1 = (lambda*Normin1 + beta*fft2(Normin2))./(Denormin+eps+intensity);
    L1 = real(ifft2(FL1));
    %% normalize L1
    for c = 1:D%%%%%%%整体可以调整为更快的方法
        L1t = L1(:,:,c);
        channel = I(:,:,c);
        options = ['Algorithm','guasi-newton']; %%guasi-newton,,,trust-region
        [t] = fminunc(@fun1,1,options,L1t,channel);
        L1(:,:,c) = L1t+t;
    end
    t = L1<lb;
    L1(t) = lb(t);
    t = L1>hb;
    L1(t) = hb(t);
    figure,imshow(L1);
    L2 = I - L1;
    n1 = abs(imfilter(L2,f3,'same'));
    n2 = abs(imfilter(L1,f3,'same'));
    n3 = imfilter(L1,f5,'same');
    res(i) = lambda*power(norm(n1(:),2),2) + norm(n2(:),1) + lambda2*norm(n3(:),2);

end
L1 = L1(pad_size+1:N-pad_size,pad_size+1:M-pad_size,:);
L2 = L2(pad_size+1:N-pad_size,pad_size+1:M-pad_size,:);
I = I(pad_size+1:N-pad_size,pad_size+1:M-pad_size,:);

%%%%%用滤波来代替迭代，边缘效果会更自然，但是天空区域也会产生光圈。迭代效果一般却不会产生光圈效应
L2_coarse = min(L2,[],3);
% L2_g2 = fastguidedfilter(max(L1,[],3), L2_coarse, 10, 10^-4,4);
L2_g2 = wlsFilter(L2_coarse,1.5, 1.2, min(L1,[],3));          %%%%%%%%%%%%对边缘的收敛较好，默认为1，1.2，高原原中用0.5，1.5更不平滑,这样边缘更锐利，L2应该是平滑的，所以取1.5 ， 1.2
% L2_g2 = L2_coarse;
% L2_g2 = qxfilter(L2_coarse,0.01,0.05);

L2 = repmat(L2_g2,[1,1,D]);
L1 = I - L2;




% gF = imfilter(L2,f3,'circular');
% gF = gF.^2*lambda;
% gFr = gF(:,:,1);
% gFx2 = -imfilter(L1,f1,'circular');  %%相当imfilter(L1,[-1,1],'circular'); 
% gFy2 = -imfilter(L1,f2,'circular');
% gF2 = beta*(gFx2 + gFy2 - gFx - gFy).^2;
% %     gF2 = min(gFxxx.^2*10000,1);
% gF2r = gF2(:,:,1);
% por = gF2r./(gFr+eps);


%%
function WFun = CalWeightFun(HazeImg)
% parameters setting
L = rgb2gray(HazeImg);
sigma = 0.5; 

% calculating the weighting function
% weighting function
% method = 'circular';
% d_r = imfilter(HazeImg(:, :, 1), D, method);
WFun = 1./(L+eps);
% WFun = exp(-(L) / sigma / 2);

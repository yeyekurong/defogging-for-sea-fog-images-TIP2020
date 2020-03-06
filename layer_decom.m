%%%Part of code refers the code of the paper "Nighttime haze removal with glow and
%%%multiple light colors", ICCV 2015
function [L1, L2] = layer_decom(I,lambda,lambda2,lb,hb,fast)
% Layer Separation using Relative Smoothness (general)
% I: input
% lambda: smoothness regularization parameter on the smoother layer Layer 2 
% lb: lower bound of the Layer 1,need to be same dimention with input I 
% hb: upper bound of the Layer 1,need to be same dimention with input I
pad_size = 64;
I = padarray(I,[pad_size,pad_size],'replicate','both');
lb = padarray(lb,[pad_size,pad_size],'replicate','both');
hb = padarray(hb,[pad_size,pad_size],'replicate','both');
if ~exist('lambda2','var')
lambda2 = 0.01;   %%%%%%0.001 
end
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

thr = 0.05;   %%%%%0.05
if fast==1
    iter = 1;
else
    iter = 4;
end
for i = 1:iter  
    beta = 2^(i-1)/thr;
    Denormin   = lambda*Denormin1 + beta*Denormin2;
    %% update g
    gFx = -imfilter(L1,f1,'circular');  %%imfilter(L1,[-1,1],'circular'); 
    gFy = -imfilter(L1,f2,'circular');
    gFx = max(abs(gFx) - 1/ beta, 0) .* sign(gFx); 
    gFy = max(abs(gFy) - 1/ beta, 0) .* sign(gFy); 
    
    %% compute L1
    Normin2 = [gFx(:,end,:) - gFx(:, 1,:), -diff(gFx,1,2)];   
    Normin2 = Normin2 + [gFy(end,:,:) - gFy(1, :,:); -diff(gFy,1,1)];
    FL1 = (lambda*Normin1 + beta*fft2(Normin2))./(Denormin+eps+intensity);
    L1 = real(ifft2(FL1));

    %% normalize L1
    for c = 1:D
        L1t = L1(:,:,c);
        channel = I(:,:,c);
        options = optimoptions('fminunc','Algorithm','quasi-newton');
        [t] = fminunc(@fun1,1,options,L1t,channel);
        L1(:,:,c) = L1t+t;
    end
    t = L1<lb;
    L1(t) = lb(t);
    t = L1>hb;
    L1(t) = hb(t);
    L2 = I - L1;
end
L1 = L1(pad_size+1:N-pad_size,pad_size+1:M-pad_size,:);
L2 = L2(pad_size+1:N-pad_size,pad_size+1:M-pad_size,:);
I = I(pad_size+1:N-pad_size,pad_size+1:M-pad_size,:);

%%%%%The fast implement
L2_coarse = min(L2,[],3);
% L2_g2 = fastguidedfilter(max(L1,[],3), L2_coarse, 10, 10^-4,4);
if fast == 1
L2_g2 = wlsFilter(L2_coarse,0.5, 1.2, min(L1,[],3));    
else
L2_g2 = im2double(qxfilter(L2_coarse*255,0.01,0.02));
end

L2 = repmat(L2_g2,[1,1,D]);
L1 = I - double(L2);


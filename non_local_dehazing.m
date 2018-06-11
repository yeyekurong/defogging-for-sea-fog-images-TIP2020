<<<<<<< HEAD
function [img_dehazed, transmission, ind_dis] = non_local_dehazing(img_hazy, LR, gamma, air_light)
%The core implementation of "Non-Local Image Dehazing", CVPR 2016
% 
% The details of the algorithm are described in our paper: 
% Non-Local Image Dehazing. Berman, D. and Treibitz, T. and Avidan S., CVPR2016,
% which can be found at:
% www.eng.tau.ac.il/~berman/NonLocalDehazing/NonLocalDehazing_CVPR2016.pdf
% If you use this code, please cite the paper.
%
%   Input arguments:
%   ----------------
%	img_hazy     - A hazy image in the range [0,255], type: uint8
%	air_light    - As estimated by prior methods, normalized to the range [0,1]
%	gamma        - Radiometric correction. If empty, 1 is assumed
%
%   Output arguments:
%   ----------------
%   img_dehazed  - The restored radiance of the scene (uint8)
%   transmission - Transmission map of the scene, in the range [0,1]
%
% Author: Dana Berman, 2016. 
%
% The software code of the Non-Local Image Dehazing algorithm is provided
% under the attached LICENSE.md

show = 0
%% Validate input
[h,w,n_colors] = size(img_hazy);
if (n_colors ~= 3) % input verification
    error(['Non-Local Dehazing reuires an RGB image, while input ',...
        'has only ',num2str(n_colors),' dimensions']);
end

% if ~exist('air_light','var') || isempty(air_light) || (numel(air_light)~=3)
%     error('Dehazing on sphere requires an RGB airlight');
% end

if ~exist('gamma','var') || isempty(gamma), gamma = 1; end

img_hazy = im2double(img_hazy);
img_hazy_corrected = img_hazy.^gamma; % radiometric correction
if (~exist('air_light','var'))
    windowSize=round(min(h,w)*0.01);
    Imax = min(img_hazy_corrected,[],3);
    % Imax = (img_hazy_corrected(:,:,1) + img_hazy_corrected(:,:,2) + img_hazy_corrected(:,:,3))/3;
    pp = minmaxfilt(double(Imax),windowSize,'min','same');
    [count,scale]=imhist(pp);
    sumpoint=0;
    i=256;
    while (sumpoint<h*w*0.0001)
        sumpoint=sumpoint+count(i);
        i=i-1;
    end
    temp=im2bw(pp,double(i-1)/255);
    % temp = ones(h,w);
    AR=double(temp).*double(img_hazy_corrected(:,:,1));
    AG=double(temp).*double(img_hazy_corrected(:,:,2));
    AB=double(temp).*double(img_hazy_corrected(:,:,3));
    A1=double(max(max(AR))); 
    A2=double(max(max(AG))); 
    A3=double(max(max(AB))); 
    Amax = max(max(A1,A2),A3);
    air_light(:,:,1) = A1;
    air_light(:,:,2) = A2;
    air_light(:,:,3) = A3;
end

%% Find Haze-lines
% Translate the coordinate system to be air_light-centric (Eq. (3))
dist_from_airlight = double(zeros(h,w,n_colors));
for color_idx=1:n_colors
    dist_from_airlight(:,:,color_idx) = img_hazy_corrected(:,:,color_idx) - air_light(:,:,color_idx);
end

% Calculate radius (Eq. (5))
radius = sqrt( dist_from_airlight(:,:,1).^2 + dist_from_airlight(:,:,2).^2 +dist_from_airlight(:,:,3).^2 );

% Cluster the pixels to haze-lines
% Use a KD-tree impementation for fast clustering according to their angles
dist_unit_radius = reshape(dist_from_airlight,[h*w,n_colors]);
dist_norm = sqrt(sum(dist_unit_radius.^2,2));
dist_unit_radius = bsxfun(@rdivide, dist_unit_radius, dist_norm);
n_points = 1000;
% load pre-calculated uniform tesselation of the unit-sphere
fid = fopen(['TR',num2str(n_points),'.txt']);
points = cell2mat(textscan(fid,'%f %f %f')) ;
fclose(fid);
mdl = KDTreeSearcher(points);
ind = knnsearch(mdl, dist_unit_radius);
ind_dis = reshape(ind,h,w);
ind_h = ind_dis(1:round(h/3),:);   %%%%天空部分的数字
K = accumarray(ind_h(:),1,[n_points,1]);
[num,ind_ix] = max(K);
ind_dis = reshape(ind,h,w);
ind_dis = im2double(ind_dis)/1000;
% figure,imshow(ind_dis);
if(show)
    xlabel = 0:1:999;
    len = accumarray(ind,radius(:),[n_points,1],@length);
    ST = accumarray(ind,radius(:),[n_points,1],@std);
    figure,plot(xlabel,len);

    len = accumarray(ind,radius(:),[n_points,1],@length);
    xlabel=1:1:1000;
    figure,plot(xlabel,len);
    figure,scatter3(points(:,1),points(:,2),points(:,3),30,len,'filled');
    colormap(copper);
    view(40,35)
    hold on;%画结果
    [xx,yy,zz]=sphere(50); 
    h2=surf(xx,yy,zz); %画一个单位球做参考
    set(h2,'edgecolor','none','facecolor','k','facealpha',0.3);
    axis equal;
    axis([-1 1 -1 1 -1 0]);
    hold off;
%% Estimating Initial Transmission
end
% Estimate radius as the maximal radius in each haze-line (Eq. (11))
K = accumarray(ind,radius(:),[n_points,1],@max);
if(show)
    figure,plot(xlabel,K,'r');
    hold on;
    plot(xlabel,ST);
    legend('max','std');
    hold off;
end
radius_new = reshape( K(ind), h, w);
    
% Estimate transmission as radii ratio (Eq. (12))
transmission_estimation = radius./radius_new;
% figure,imshow(transmission_estimation);
% Limit the transmission to the range [trans_min, 1] for numerical stability
trans_min = 0.1;
transmission_estimation = min(max(transmission_estimation, trans_min),1);


%% Regularization

% Apply lower bound from the image (Eqs. (13-14))
trans_lower_bound = 1 - min(bsxfun(@rdivide,img_hazy_corrected,reshape(air_light,1,1,3)) ,[],3);
transmission_estimation = max(transmission_estimation, trans_lower_bound);

% Solve optimization problem (Eq. (15))
% find bin counts for reliability - small bins (#pixels<50) do not comply with 
% the model assumptions and should be disregarded
bin_count       = accumarray(ind,1,[n_points,1]);
bin_count_map   = reshape(bin_count(ind),h,w);   %图片每个像素点对应灰度的像素数量
bin_eval_fun    = @(x) min(1, x/50);

% Calculate std - this is the data-term weight of Eq. (15)
K_std = accumarray(ind,radius(:),[n_points,1],@std);   %每个line上的半径的标准差
radius_std = reshape( K_std(ind), h, w);              %图片每个像素点对应半径的标准差
radius_eval_fun = @(r) min(1, 3*max(0.001, r-0.1));
%%%%%
lam = adjust(LR(:,:,1));  
tr = transmission_estimation(ind==ind_ix);
tr = min(power(max(transmission_estimation/min(tr),1),1.5),100);  %****/mean(tr)
c2 = 2./(2+exp((lam-1).*tr));
transmission_estimation = 1-c2+c2.*transmission_estimation;
%%%%%%%
radius_reliability = radius_eval_fun(radius_std./max(radius_std(:)));%每个像素通过对应的半径的标准差来衡量可信度
data_term_weight   = bin_eval_fun(bin_count_map).*radius_reliability; %权重
lambda = 0.02;
transmission = wls_optimization(transmission_estimation, data_term_weight, img_hazy.^0.7, lambda);


%% Dehazing
% (Eq. (16))
img_dehazed = zeros(h,w,n_colors);
leave_haze = 1.06; % leave a bit of haze for a natural look (set to 1 to reduce all haze)
for color_idx = 1:3
    img_dehazed(:,:,color_idx) = ( img_hazy_corrected(:,:,color_idx) - ...
        (1-leave_haze.*transmission).*air_light(color_idx) )./ max(transmission,trans_min);
end

% Limit each pixel value to the range [0, 1] (avoid numerical problems)
img_dehazed(img_dehazed>1) = 1;
img_dehazed(img_dehazed<0) = 0;
img_dehazed = img_dehazed.^(1/gamma); % radiometric correction


end % function non_local_dehazing
=======
function [img_dehazed, transmission, ind_dis] = non_local_dehazing(img_hazy, LR, gamma, air_light)
%The core implementation of "Non-Local Image Dehazing", CVPR 2016
% 
% The details of the algorithm are described in our paper: 
% Non-Local Image Dehazing. Berman, D. and Treibitz, T. and Avidan S., CVPR2016,
% which can be found at:
% www.eng.tau.ac.il/~berman/NonLocalDehazing/NonLocalDehazing_CVPR2016.pdf
% If you use this code, please cite the paper.
%
%   Input arguments:
%   ----------------
%	img_hazy     - A hazy image in the range [0,255], type: uint8
%	air_light    - As estimated by prior methods, normalized to the range [0,1]
%	gamma        - Radiometric correction. If empty, 1 is assumed
%
%   Output arguments:
%   ----------------
%   img_dehazed  - The restored radiance of the scene (uint8)
%   transmission - Transmission map of the scene, in the range [0,1]
%
% Author: Dana Berman, 2016. 
%
% The software code of the Non-Local Image Dehazing algorithm is provided
% under the attached LICENSE.md

show = 0
%% Validate input
[h,w,n_colors] = size(img_hazy);
if (n_colors ~= 3) % input verification
    error(['Non-Local Dehazing reuires an RGB image, while input ',...
        'has only ',num2str(n_colors),' dimensions']);
end

% if ~exist('air_light','var') || isempty(air_light) || (numel(air_light)~=3)
%     error('Dehazing on sphere requires an RGB airlight');
% end

if ~exist('gamma','var') || isempty(gamma), gamma = 1; end

img_hazy = im2double(img_hazy);
img_hazy_corrected = img_hazy.^gamma; % radiometric correction
if (~exist('air_light','var'))
    windowSize=round(min(h,w)*0.01);
    Imax = min(img_hazy_corrected,[],3);
    % Imax = (img_hazy_corrected(:,:,1) + img_hazy_corrected(:,:,2) + img_hazy_corrected(:,:,3))/3;
    pp = minmaxfilt(double(Imax),windowSize,'min','same');
    [count,scale]=imhist(pp);
    sumpoint=0;
    i=256;
    while (sumpoint<h*w*0.0001)
        sumpoint=sumpoint+count(i);
        i=i-1;
    end
    temp=im2bw(pp,double(i-1)/255);
    % temp = ones(h,w);
    AR=double(temp).*double(img_hazy_corrected(:,:,1));
    AG=double(temp).*double(img_hazy_corrected(:,:,2));
    AB=double(temp).*double(img_hazy_corrected(:,:,3));
    A1=double(max(max(AR))); 
    A2=double(max(max(AG))); 
    A3=double(max(max(AB))); 
    Amax = max(max(A1,A2),A3);
    air_light(:,:,1) = A1;
    air_light(:,:,2) = A2;
    air_light(:,:,3) = A3;
end

%% Find Haze-lines
% Translate the coordinate system to be air_light-centric (Eq. (3))
dist_from_airlight = double(zeros(h,w,n_colors));
for color_idx=1:n_colors
    dist_from_airlight(:,:,color_idx) = img_hazy_corrected(:,:,color_idx) - air_light(:,:,color_idx);
end

% Calculate radius (Eq. (5))
radius = sqrt( dist_from_airlight(:,:,1).^2 + dist_from_airlight(:,:,2).^2 +dist_from_airlight(:,:,3).^2 );

% Cluster the pixels to haze-lines
% Use a KD-tree impementation for fast clustering according to their angles
dist_unit_radius = reshape(dist_from_airlight,[h*w,n_colors]);
dist_norm = sqrt(sum(dist_unit_radius.^2,2));
dist_unit_radius = bsxfun(@rdivide, dist_unit_radius, dist_norm);
n_points = 1000;
% load pre-calculated uniform tesselation of the unit-sphere
fid = fopen(['TR',num2str(n_points),'.txt']);
points = cell2mat(textscan(fid,'%f %f %f')) ;
fclose(fid);
mdl = KDTreeSearcher(points);
ind = knnsearch(mdl, dist_unit_radius);
ind_dis = reshape(ind,h,w);
ind_h = ind_dis(1:round(h/3),:);   %%%%天空部分的数字
K = accumarray(ind_h(:),1,[n_points,1]);
[num,ind_ix] = max(K);
ind_dis = reshape(ind,h,w);
ind_dis = im2double(ind_dis)/1000;
% figure,imshow(ind_dis);
if(show)
    xlabel = 0:1:999;
    len = accumarray(ind,radius(:),[n_points,1],@length);
    ST = accumarray(ind,radius(:),[n_points,1],@std);
    figure,plot(xlabel,len);

    len = accumarray(ind,radius(:),[n_points,1],@length);
    xlabel=1:1:1000;
    figure,plot(xlabel,len);
    figure,scatter3(points(:,1),points(:,2),points(:,3),30,len,'filled');
    colormap(copper);
    view(40,35)
    hold on;%画结果
    [xx,yy,zz]=sphere(50); 
    h2=surf(xx,yy,zz); %画一个单位球做参考
    set(h2,'edgecolor','none','facecolor','k','facealpha',0.3);
    axis equal;
    axis([-1 1 -1 1 -1 0]);
    hold off;
%% Estimating Initial Transmission
end
% Estimate radius as the maximal radius in each haze-line (Eq. (11))
K = accumarray(ind,radius(:),[n_points,1],@max);
if(show)
    figure,plot(xlabel,K,'r');
    hold on;
    plot(xlabel,ST);
    legend('max','std');
    hold off;
end
radius_new = reshape( K(ind), h, w);
    
% Estimate transmission as radii ratio (Eq. (12))
transmission_estimation = radius./radius_new;
% figure,imshow(transmission_estimation);
% Limit the transmission to the range [trans_min, 1] for numerical stability
trans_min = 0.1;
transmission_estimation = min(max(transmission_estimation, trans_min),1);


%% Regularization

% Apply lower bound from the image (Eqs. (13-14))
trans_lower_bound = 1 - min(bsxfun(@rdivide,img_hazy_corrected,reshape(air_light,1,1,3)) ,[],3);
transmission_estimation = max(transmission_estimation, trans_lower_bound);

% Solve optimization problem (Eq. (15))
% find bin counts for reliability - small bins (#pixels<50) do not comply with 
% the model assumptions and should be disregarded
bin_count       = accumarray(ind,1,[n_points,1]);
bin_count_map   = reshape(bin_count(ind),h,w);   %图片每个像素点对应灰度的像素数量
bin_eval_fun    = @(x) min(1, x/50);

% Calculate std - this is the data-term weight of Eq. (15)
K_std = accumarray(ind,radius(:),[n_points,1],@std);   %每个line上的半径的标准差
radius_std = reshape( K_std(ind), h, w);              %图片每个像素点对应半径的标准差
radius_eval_fun = @(r) min(1, 3*max(0.001, r-0.1));
%%%%%
lam = adjust(LR(:,:,1));  
tr = transmission_estimation(ind==ind_ix);
tr = min(power(max(transmission_estimation/min(tr),1),1.5),100);  %****/mean(tr)
c2 = 2./(2+exp((lam-1).*tr));
transmission_estimation = 1-c2+c2.*transmission_estimation;
%%%%%%%
radius_reliability = radius_eval_fun(radius_std./max(radius_std(:)));%每个像素通过对应的半径的标准差来衡量可信度
data_term_weight   = bin_eval_fun(bin_count_map).*radius_reliability; %权重
lambda = 0.02;
transmission = wls_optimization(transmission_estimation, data_term_weight, img_hazy.^0.7, lambda);


%% Dehazing
% (Eq. (16))
img_dehazed = zeros(h,w,n_colors);
leave_haze = 1.06; % leave a bit of haze for a natural look (set to 1 to reduce all haze)
for color_idx = 1:3
    img_dehazed(:,:,color_idx) = ( img_hazy_corrected(:,:,color_idx) - ...
        (1-leave_haze.*transmission).*air_light(color_idx) )./ max(transmission,trans_min);
end

% Limit each pixel value to the range [0, 1] (avoid numerical problems)
img_dehazed(img_dehazed>1) = 1;
img_dehazed(img_dehazed<0) = 0;
img_dehazed = img_dehazed.^(1/gamma); % radiometric correction


end % function non_local_dehazing
>>>>>>> abf47cccc469294dc9ac70fd4e5258faff5a1a04

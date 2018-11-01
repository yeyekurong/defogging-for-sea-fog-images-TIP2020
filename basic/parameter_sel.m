function [alpha, beta, pro] = parameter_sel(img_hazy)
I = img_hazy;
gray = max(I,[],3);
gray = adjust(gray);
f = fspecial('gaussian',10,10);
f3 = [0, -1, 0; 
      -1, 4, -1;
      0, -1, 0];
shap = abs(imfilter(gray,f3));
shap = shap(10:end-10,10:end-10);
shap = minmaxfilt(shap,20,'max','same');
gray = imfilter(gray,f);
gray = imfilter(gray,f);
gray = imfilter(gray,f);
gray = imfilter(gray,f);
gray = gray(10:end-10,10:end-10);
% figure,imshow(gray(10:end-10,10:end-10));
ratio = sum(gray(:))/sum(shap(:))
alpha = 500;
beta = 0.001;
pro = 1.6;
if(ratio>2.3)
    beta = 0.01;
    pro = 3.3;
end
end
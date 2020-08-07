function [imgOut,scale]=getMaskNormalizedImg(img,mask)
index=find(mask);
scale=median(img(index));
%scale=mean(img(index));
imgOut=img/scale;
function [imgOut,scale]=getOccNormalizedImg(img,occMask)
index=find(occMask);
scale=median(img(index));
%scale=mean(img(index));
imgOut=img/scale;
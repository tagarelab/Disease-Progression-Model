function [imgs_rescaled, scales, scale2s] = rescale_by_occ_lobe(imgs, occMask)
%RESCALE_BY_OCC_LOBE Summary of this function goes here
%   Detailed explanation goes here
nImgs=size(imgs,4);
occIndex=find(occMask>0.5);

imgs_rescaled = zeros(size(imgs));
scales = zeros(1, nImgs);
scale2s = zeros(1, nImgs);

for i=1:nImgs,
    tmp=squeeze(imgs(:,:,:,i));
    occ_values = tmp(occIndex);
    occMed=median(occ_values);
    occMean = mean(occ_values);
    tmp=tmp/occMed;
    
    imgs_rescaled(:,:,:,i) = tmp;
    scales(i) = occMed;
    scale2s(i) = occMean;
end

end
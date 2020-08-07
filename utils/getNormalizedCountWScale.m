function [counts,scales,imgs_rescaled, scale2s]=getNormalizedCountWScale(imgs,occMask,sMask)
nImgs=size(imgs,4);

[imgs_rescaled, scales, scale2s] = rescale_by_occ_lobe(imgs, occMask);

for i=1:nImgs,
    tmp=squeeze(imgs_rescaled(:,:,:,i));
    tmp=tmp.*sMask;
    counts(i)=sum(tmp(:));
end




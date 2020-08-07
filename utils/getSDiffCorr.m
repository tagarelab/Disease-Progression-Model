function counts=getSDiffCorr(imgs,occMask,sMask)
nImgs=size(imgs,4);
occIndex=find(occMask>0.5);

for i=1:nImgs,
    tmp=squeeze(imgs(:,:,:,i));
    occMed=median(tmp(occIndex));
    tmp=tmp/occMed;
    tmp=tmp.*sMask;
    counts(i)=sum(tmp(:));
end




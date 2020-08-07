function cc=getMaskNormalizedImg(img1,img2,mask)
% imgSize=size(img1);
% img1=img1.*mask;
% img2=img2.*mask;
% tmp=[reshape(img1,[prod(imgSize) 1]) reshape(img2,[prod(imgSize) 1])];
% corr=corrcoef(tmp);
% cc=(corr(1,2));

 index=find(mask>0.5);
% tmp1=img1(index);
% tmp1=tmp1/norm(tmp1,2);
% tmp2=img2(index);
% tmp2=tmp2/norm(tmp2,2);
% cc=tmp1'*tmp2;

tmp=[img1(index) img2(index)];
cc=corrcoef(tmp);
cc=cc(1,2);




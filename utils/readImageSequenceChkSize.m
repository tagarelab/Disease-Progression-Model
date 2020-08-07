function [imgs,elapsedTime,imgs_raw,flags,firstDate]=readImageSequenceChkSize(dirName,chkSize)

scanDir=strcat(dirName,'\Reconstructed_DaTSCAN');
dataDirs=dir(scanDir);
nImgs= numel(dataDirs)-2;
imgs=zeros([chkSize nImgs]);
elapsedTime=zeros(nImgs,1);
valid=zeros(nImgs,1);

elapsedTime=[];
imgs=[];
imgs_raw = [];
flags = struct('isNegative', {}, 'sex', {});

for i=3:numel(dataDirs);
    curDir=strcat(scanDir,'\',dataDirs(i).name);
    SDir=dir(curDir);
    curDir=strcat(curDir,'\',SDir(3).name);
    dicomFile=dir(strcat(curDir,'\*.dcm'));
    dicomFileName=strcat(curDir,'\',dicomFile(1).name);
    
    %class(dicomread(dicomFileName))
    
    info = dicominfo(dicomFileName);
    sex = info.PatientSex;
    tmpImg=double(squeeze((dicomread(dicomFileName))));
    img_raw = tmpImg;
    % Flip left and right since increasing x points to the left of the
    % patient. Since the center of the MNI atlas is at the 45th column, we
    % flip the image along the 45th column.
    tmpImg = fliplr(tmpImg); tmpImg = tmpImg(:,[3:91,91,91],:);
    img_raw = fliplr(img_raw); img_raw = img_raw(:,[3:91,91,91],:);
    
    test_sample = tmpImg(1:5,1:5,5);
    
    % Checked all the images. Once there are negative values, the minimum
    % will always be -32,768.
    if mean(test_sample(:)) < 0 % range is in -32,768 to 32,767, need to rescale
%         tmpImg=double(squeeze(im2uint16(dicomread(dicomFileName))));
        tmpImg = (tmpImg + 32768) / 2;
        isNeg = 1;
    else
        isNeg = 0;
    end

    if (numel(tmpImg)==prod(chkSize))
        imgs(:,:,:,i-2)=tmpImg;
        imgs_raw(:,:,:,i-2) = img_raw;
        elapsedTime(i-2)= getDate(dataDirs(i).name(1:10));
        valid(i-2)=1;
        flags(i-2).isNegative = isNeg;
        flags(i-2).sex = sex;
    end
end
%Get the valid images
index=find(valid);
if isempty(index)==0
    imgs=imgs(:,:,:,index);
    imgs_raw = imgs_raw(:,:,:,index);
    flags = flags(index);
    elapsedTime=elapsedTime(index);
    firstDate = elapsedTime(1);
    elapsedTime=elapsedTime-elapsedTime(1);
end


function date=getDate(dateStr)
yr=str2num(dateStr(1:4));
month=str2num(dateStr(6:7));
day=str2num(dateStr(9:10));
date=datenum(yr,month,day);
function [imgs,elapsedTime]=readImageSequence(dirName)

scanDir=strcat(dirName,'\Reconstructed_DaTSCAN');
dataDirs=dir(scanDir);

elapsedTime=[];
imgs=[];
for i=3:numel(dataDirs);
    curDir=strcat(scanDir,'\',dataDirs(i).name);
    SDir=dir(curDir);
    curDir=strcat(curDir,'\',SDir(3).name);
    dicomFile=dir(strcat(curDir,'\*.dcm'));
    dicomFileName=strcat(curDir,'\',dicomFile(1).name);
    
    tmpImg=double(squeeze(dicomread(dicomFileName)));
   % i
    if i==3
        imgs=zeros([size(tmpImg) numel(dataDirs)-2]);
    end
    
    
        imgs(:,:,:,i-2)=tmpImg;
    elapsedTime=[elapsedTime; getDate(dataDirs(i).name(1:10))];
end
elapsedTime=elapsedTime-elapsedTime(1);


function date=getDate(dateStr)
yr=str2num(dateStr(1:4));
month=str2num(dateStr(6:7));
day=str2num(dateStr(9:10));
date=datenum(yr,month,day);
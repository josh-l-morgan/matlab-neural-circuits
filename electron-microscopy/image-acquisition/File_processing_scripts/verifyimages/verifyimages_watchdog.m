function verifyimages_watchdog(sourcedir,sourcefiletemplate)
%verifyimages.m
%A script for online EM image checking
%Give a directory that should be watched and a template for files that should be investigated
%This assumes that the source files will be images which matlab can open
%Runs continuously and checks the directory for new files and scans them
%Press ESC to exit
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, March 2010

%close all;
%clear all;

%sourcedir='z:/Thousand+sectionRunCortexRoto/Thousand_high_res/Thousand_w04_high_res_st1-8/';
%sourcefiletemplate='*.tif';

targetdir_downscaled='c:/Daniel/verifyimages/downscaled/';
downscaledtargetsize=2048;
downscaledfileextension='_downscaled.png';
downscaledfiletype='png';
targetdir_cropped='c:/Daniel/verifyimages/center/';
croppedfileextension='_center.png';
croppedfiletype='png';
croppedtargetsize=1024;
croppedxoffset=1024;
croppedyoffset=1024;

  
outlierspos=[]; %for scanfaults
overallfilelist={};
overallptr=0;
exitloop=0;

while ~exitloop
  %Get list of filenames in watched directory
  txt=sprintf('Scanning for files matching "%s" in directory %s ...',sourcefiletemplate,sourcedir); disp(txt);
  filelist=getfilelist(sourcedir,sourcefiletemplate);
  filelist=getnewnames(filelist,overallfilelist); %only scan new files  
  
  nroffiles=size(filelist,2);
  if (nroffiles>0)
    txt=sprintf('%d new files matching "%s" found in directory %s; security wait 60 seconds',nroffiles,sourcefiletemplate,sourcedir);
    disp(txt);
  end;
%   txt=sprintf('%d new files matching "%s" found in directory %s',nroffiles,sourcefiletemplate,sourcedir);
%   disp(txt);
  
  %After scanning for files we wait 60 seconds, which should make every new files being written completely
  k=getkeywait(2); %Wait 60 seconds for keypress
  if k==27 %if ESC was pressed
    exitloop=1;
    disp('ESC PRESSED - Program will exit after any remaining images have been processed...');
  end;
  
  %Then we process these new files
  for f=1:1:nroffiles
    filenamewithpath=filelist{f};
    filename=filenameonly(filenamewithpath);
    
    if exist(filenamewithpath,'file')==0
      txt=sprintf('Warning: File %s not found.',filenamewithpath); disp(txt);
    else
      txt=sprintf('Loading image %s ...',filenamewithpath); disp(txt);
      image=imread(filenamewithpath);
      
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       % check for scanfaults
%       
%       scanerrorthreshold=10;
%       
%       %%%% Compute differences of lines in the image
%       %disp('Computing ...');
%       diffimg=single(image(2:size(image,1),:))-single(image(1:size(image,1)-1,:));
%       vdiffmean=mean(abs(diffimg),2);
%       vdiffstd=std(diffimg,0,2);
%       
%       %%%% Outlier detection of STD measures
%       st=std(vdiffstd);
%       stdoutliers=sum(vdiffstd-mean(vdiffstd)>st*scanerrorthreshold);
%       if stdoutliers>0
%         outlierspos{overallptr+f}=find(vdiffstd-mean(vdiffstd)>st*scanerrorthreshold);
%         txt=sprintf('  %s:  %d scan errors found; at:',filename,stdoutliers);
%         disp(txt);
%         disp([outlierspos{f}]);
%       else
%         outlierspos{overallptr+f}=[];
%         txt=sprintf('  %s:  No scan errors found.',filename); disp(txt);
%       end;
%       
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % downsample and save
      %txt=sprintf('  Downscaling original to %d x %d ...',downscaledtargetsize, downscaledtargetsize); disp(txt);
      dimage=imresize(image,[downscaledtargetsize downscaledtargetsize]);
      targetname=[targetdir_downscaled filename(1:end-4) downscaledfileextension];
      txt=sprintf('  Writing downscaled image to %s',targetname); disp(txt);
      imwrite(dimage,targetname,downscaledfiletype);
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % crop center and save
      %txt=sprintf('  Cropping original to %d x %d ...',croppedtargetsize, croppedtargetsize); disp(txt);
      midx=size(image,2)/2;
      midy=size(image,1)/2;
      %xstart=floor(midx-(croppedtargetsize/2)); xend=xstart+croppedtargetsize-1;
      xstart=floor(midx-(croppedtargetsize/2)); 
      xstart=xstart+croppedxoffset;
      xend=xstart+croppedtargetsize-1;
      if (xstart<1)
        xstart=1;
      end;
      if (xend>size(image,2))
        xend=size(image,2);
      end;
      %ystart=floor(midy-(croppedtargetsize/2)); yend=ystart+croppedtargetsize-1;
      ystart=floor(midy-(croppedtargetsize/2));
      ystart=ystart+croppedyoffset;
      yend=ystart+croppedtargetsize-1;
      if (ystart<1)
        ystart=1;
      end;
      if (yend>size(image,1))
        yend=size(image,1);
      end;
      if (yend>=ystart)&&(xend>=xstart)
        dimage=image(ystart:yend,xstart:xend);
        targetname=[targetdir_cropped filename(1:end-4) croppedfileextension];
        txt=sprintf('  Writing cropped image to %s',targetname); disp(txt);
        imwrite(dimage,targetname,croppedfiletype);
      else
        disp('  WARNING: Cannot write cropped image because the selected FOV is outside of the image.');
      end;
    end;
  end;
  overallfilelist=[filelist overallfilelist];
  overallptr=size(overallfilelist,2);
end;

disp('---------------------------------------------------------------------');
disp('Writing out a text file that reports scan faults...');
disp('---------------------------------------------------------------------');
scanfaultfilename='scanfaults.txt';


% isfound=0;
% fid = fopen(scanfaultfilename, 'w+');
% if max(size(outlierspos)>0)
%   for f=1:1:size(outlierspos,2) %nroffiles
%     nrofscanfaults=max(size(outlierspos{f}));
%     if nrofscanfaults>0
%       name=overallfilelist{f}; %sprintf(param.baserawname,slice+param.firstslice-1);
%       if (nrofscanfaults==1)
%         fprintf(fid, '%d scan error  found in image %s, at:',nrofscanfaults,name);
%       else
%         fprintf(fid, '%d scan errors found in image %s, at:',nrofscanfaults,name);
%       end;
%       scanfaultline=outlierspos{f};
%       for i=1:1:nrofscanfaults
%         fprintf(fid,' %d',scanfaultline(i));
%       end;
%       fprintf(fid,'\r\n');
%       isfound=1;
%     end;
%   end;
%   if isfound==0
%     fprintf(fid, 'No scan errors have been found in directory %s and subdirectories.',sourcedir);
%   end;
% else
%   fprintf(fid, 'No scan errors have been found in directory %s and subdirectories.',sourcedir);
% end;
fclose(fid);

function gather_downsample_images(sourcedir,sourcefiletemplate)
%gather_images.m
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

%targetdir_gathered='E:/SEM_Users/Richard/Zhen_c-elegan/C-elegan_stack3_processed/cropped/'; 
targetdir_downscaled='E:/SEM_Users/Richard/Zhen_c-elegan/C-elegan_w03/Mosaics_downsampled/';
downscaledtargetsize=4096;
%gatheredfileextension='.tif';
downscaledfileextension='_downscaled.tif';
%gatheredfiletype='tif'; 
downscaledfiletype='tif';
% targetdir_cropped='./center/';
% croppedfileextension='_center.png';
% croppedfiletype='png';
% croppedtargetsize=2048;
% croppedxoffset=2048;
% croppedyoffset=2048;

  
% outlierspos=[]; %for scanfaults
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
  k=getkeywait(20); %Wait 60 seconds for keypress
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
      %targetname=[targetdir_gathered filename(1:end-4) gatheredfileextension];
      %imwrite(image,targetname,gatheredfiletype);
      
       txt=sprintf('  Downscaling original to %d x %d ...',downscaledtargetsize, downscaledtargetsize); disp(txt);
      dimage=imresize(image,[downscaledtargetsize downscaledtargetsize]);
      targetname=[targetdir_downscaled filename(1:end-4) downscaledfileextension];
      txt=sprintf('  Writing downscaled image to %s',targetname); disp(txt);
      imwrite(dimage,targetname,downscaledfiletype);
      
    end;
  end;
  overallfilelist=[filelist overallfilelist];
  overallptr=size(overallfilelist,2);
end;    
    
      
     
      
%fclose(fid);

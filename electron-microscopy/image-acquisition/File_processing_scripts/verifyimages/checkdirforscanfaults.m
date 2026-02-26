function checkdirforscanfaults(directory)
%Checks for scan faults in all images in the given directory and writes the result to 'scanfaults.txt'
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, March 2010

tiflist=gettiflist(directory);
txt=sprintf('%d TIFF files found in directory %s (incl. subdirectories).',size(tiflist,2),directory);
disp(txt);

disp('---------------------------------------------------------------------');
disp('Checking images for scan faults ...');
disp('---------------------------------------------------------------------');

for f=1:1:size(tiflist,2)
  filename=tiflist{f};
  
  lparam.scanerrorthreshold=10;
  %param.checkimages.scanerrorpos=pipeline_checkforscanerrors(param);
  outlierspos=[];
  
  %threshold=10; %lower this to detect more outliers
  %lparam=param.checkimages;
  
  %Use only slice number in name if there is only one row and column
  %outlierspos=zeros(param.nrofslices,1); %if we initialize this as an array, cell storage fails..
  %   for slice=1:param.nrofslices
  %     name=sprintf(param.baserawname,slice+param.firstslice-1); filename=sprintf('%s%s',param.rawdir,name);
  if exist(filename,'file')
    txt=sprintf('Loading image %s ...',filename);
    disp(txt);
    image=imread(filename);
    
    %%%% Compute differences of lines in the image
    disp('Computing ...');
    diffimg=single(image(2:size(image,1),:))-single(image(1:size(image,1)-1,:));
    vdiffmean=mean(abs(diffimg),2);
    vdiffstd=std(diffimg,0,2);
    
    %%%% Outlier detection of STD measures
    st=std(vdiffstd);
    stdoutliers=sum(vdiffstd-mean(vdiffstd)>st*lparam.scanerrorthreshold);
    if stdoutliers>0
      outlierspos{f}=find(vdiffstd-mean(vdiffstd)>st*lparam.scanerrorthreshold);
      txt=sprintf('  %d scan errors found; at:',stdoutliers);
      disp(txt);
      disp([outlierspos{f}]);
    else
      disp('  No scan errors found.');
    end;
  else
    txt=sprintf('Warning: File %s not found.',imagename);
    disp(txt);
  end;
end;


disp('---------------------------------------------------------------------');
disp('Writing out a text file that reports scan faults...');
disp('---------------------------------------------------------------------');
lparam.filename='scanfaults.txt';
isfound=0;


fid = fopen(lparam.filename, 'w+');
if max(size(outlierspos)>0)
  for f=1:1:size(tiflist,2) %slice=1:1:size(param.checkimages.scanerrorpos,2) %param.nrofslices
    nrofscanfaults=max(size(outlierspos{f}));
    if nrofscanfaults>0
      name=tiflist{f}; %sprintf(param.baserawname,slice+param.firstslice-1);
      if (nrofscanfaults==1)
        fprintf(fid, '%d scan error  found in image %s, at:',nrofscanfaults,name);
      else
        fprintf(fid, '%d scan errors found in image %s, at:',nrofscanfaults,name);
      end;
      scanfaultline=outlierspos{f};
      for i=1:1:nrofscanfaults
        fprintf(fid,' %d',scanfaultline(i));
      end;
      fprintf(fid,'\r\n');
      isfound=1;
    end;
  end;
  if isfound==0
    fprintf(fid, 'No scan errors have been found in directory %s and subdirectories.',directory);
  end;
else
  fprintf(fid, 'No scan errors have been found in directory %s and subdirectories.',directory);
end;
fclose(fid);
end

function tiflist=gettiflist(directory)
%returns a cell array of filenames
  d=dir(directory);
  nrtifs=0;
  tiflist={};
  nrf=size(d,1);
  for f=1:1:nrf
    if d(f).isdir
      if (strcmp(d(f).name,'.'))||(strcmp(d(f).name,'..'))
      else
        %go into subdirectory
        tiflist2=gettiflist([directory '/' d(f).name]);
        if size(tiflist2,1)>0
          tiflist=[tiflist tiflist2];
          nrtifs=size(tiflist,2);
        end;
      end;
    else
      ending=d(f).name(end-3:end);
      if strcmp(ending,'.tif')||strcmp(ending,'.TIF')
        nrtifs=nrtifs+1;
        tiflist{nrtifs}=[directory '/' d(f).name];
      end;
    end;
  end;
end

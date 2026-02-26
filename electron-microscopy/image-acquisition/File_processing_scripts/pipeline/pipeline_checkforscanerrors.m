function outlierspos=pipeline_checkforscanerrors(param)
%This function finds scan errors in EM images.
%
%Usage: outlierspos = checkforscanerrors(imagename)
%
%  - imagename is the name (including path) of the image to be analyzed.
%  - The results are displayed in a new figure.
%  - The Y-positions of detected faults (if any) are returned.
%
%This version is made for use in the processing pipeline
%by Daniel Berger for MIT-BCS (Seung) and Harvard (Lichtman), May 2009

outlierspos=[];

%threshold=10; %lower this to detect more outliers
lparam=param.checkimages;

if (param.nrofrows==1)&&(param.nrofcolumns==1)
  %Use only slice number in name if there is only one row and column
  %outlierspos=zeros(param.nrofslices,1); %if we initialize this as an array, cell storage fails..
  for slice=1:param.nrofslices
    name=sprintf(param.baserawname,slice+param.firstslice-1); filename=sprintf('%s%s',param.rawdir,name);
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
        outlierspos{slice}=find(vdiffstd-mean(vdiffstd)>st*lparam.scanerrorthreshold);
        txt=sprintf('  %d scan errors found; at:',stdoutliers);
        disp(txt);
        disp([outlierspos{slice}]);
      else
        disp('  No scan errors found.');
      end;
      
      %%%% Display results
%       figure;
%       title(imagename);
%       plot(vdiffstd,'r');
%       hold on;
%       plot(vdiffmean);
%       if max(size(stdoutliers))>0
%         plot(outlierspos,vdiffstd(outlierspos),'r*');
%       end;
%       hold off;
%       grid on;
%       title(imagename);
%       xlabel('Y position');
%       ylabel('Scanline difference (blue: absmean, red: std)');
    else
      txt=sprintf('Warning: File %s not found.',imagename);
      disp(txt);
    end;
  end;
else
  %Use slice number, row and column in name for several rows/columns
  %outlierspos=zeros(nrofslices,nrofrows,nrofcolumns);
  for slice=1:param.nrofslices
    for row=1:param.nrofrows
      for column=1:param.nrofcolumns
        name=sprintf(param.baserawname,slice,row,column); filename=sprintf('%s%s',param.rawdir,name);
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
            outlierspos{slice,row,column}=find(vdiffstd-mean(vdiffstd)>st*lparam.scanerrorthreshold);
            txt=sprintf('  %d scan errors found; at:',stdoutliers);
            disp(txt);
            disp([outlierspos{slice,row,column}]);
          else
            disp('  No scan errors found.');
          end;
        else
          txt=sprintf('Warning: File %s not found.',filename);
          disp(txt);
        end;
      end;
    end;
  end;
end;
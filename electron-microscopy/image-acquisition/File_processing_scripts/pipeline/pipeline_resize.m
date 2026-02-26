function [rawsize,scaledsize]=pipeline_resize(param)
%This function is part of the image processing pipeline, called by runpipeline.m 
%based on resize.m
%By Daniel Berger for MIT-BCS Seung, June 1 2009

lparam=param.downscale;

for row=1:1:param.nrofrows
  for column=1:1:param.nrofcolumns
    for slice=1:param.nrofslices
      if (param.nrofrows==1)&&(param.nrofcolumns==1)
        name=sprintf(param.baserawname,slice);
      else
        name=sprintf(param.baserawname,slice,row,column);
      end;
      filename=sprintf('%s%s',param.rawdir,name);
      %disp(filename);
      %name=sprintf(param.baserawname,slice); filename=sprintf('%s%s',param.rawdir,name);
      if exist(filename,'file')
        txt=sprintf('Loading image %s ...',filename);
        disp(txt);
        image=imread(filename);
        disp('  Scaling ...');
        if lparam.invert==1
          image=double(255-image)/255;
        else
          image=double(image)/255;
        end;
        rimage=imresize(image,1/lparam.scale);
        
        if (param.nrofrows==1)&&(param.nrofcolumns==1)
          name=sprintf(param.basescaledname,slice);
        else
          name=sprintf(param.basescaledname,slice,row,column);
        end;
        filename=sprintf('%s%s',param.scaleddir,name);
        %wname=sprintf(param.basescaledname,slice); filename=sprintf('%s%s',param.scaleddir,wname);
        txt=sprintf('  Writing image %s ...',filename);
        disp(txt);
        imwrite(rimage,filename,lparam.targetfileformat);
        %imwrite(rimage,wname,'jpg','Compression',50);
      else
        txt=sprintf('WARNING: File %s not found.', filename);
        disp(txt);
      end;
    end;
  end;
end;

rawsize=size(image);
scaledsize=size(rimage);

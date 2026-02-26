function param=pipeline_unwafer(param)
%A function to copy files from a wafer-slice naming, and possibly missing
%slices, to a single stack numbered by slices (copies from original to raw)
%for use with the pipeline
%By Daniel Berger for MIT-BCS Seung, June 7th 2009

lparam=param.unwafer;

nrofslices=0; targetslice=0;
for wafer=1:1:lparam.nrofwafers
  for slice=1:1:lparam.nrofslices(wafer)
    for row=1:1:param.nrofrows
      for column=1:1:param.nrofcolumns
        name=sprintf(param.baseorigname,wafer,slice,row,column);
        filename=sprintf('%s%s',param.origdir,name);
        if exist(filename,'file')
          if (row==1)&&(column==1)  %bad hack.. works if always all images of a slice are missing if one is missing
            targetslice=targetslice+1;
            nrofslices=nrofslices+1;
          end;
          
          txt=sprintf('Loading %s ...',filename);
          disp(txt);
          image=imread(filename);
          disp('  Processing ...');
          if lparam.invert==1
            image=double(255-image)/255;
          else
            image=double(image)/255;
          end;
          
          name=sprintf(param.baserawname,targetslice,row,column);
          filename=sprintf('%s%s',param.rawdir,name);
          txt=sprintf('  Writing image %s ...',filename);
          disp(txt);
          imwrite(image,filename,lparam.targetfileformat);
        else
          txt=sprintf('Warning: File %s not found.',filename);
          disp(txt);
        end;
      end;
    end;
  end;
end;
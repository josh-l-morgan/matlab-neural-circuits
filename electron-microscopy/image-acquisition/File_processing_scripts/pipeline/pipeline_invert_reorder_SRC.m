function pipeline_invert_reorder_SRC(param)
%This script is made to reorder and rename raw images
%For usage with the pipeline
%By Daniel Berger for MIT-BCS Seung, June 2 2009
lparam=param.reorder;

if lparam.dopart
  startslice=lparam.firstslice;
  endslice=lparam.lastslice;
else
  startslice=1;
  endslice=param.nrofslices;
end;

if (param.nrofrows==1)&&(param.nrofcolumns==1)
  for slice=startslice:1:endslice
    name=sprintf(param.baseorigname,lparam.sourceslice(slice)); filename=sprintf('%s%s',param.origdir,name);
    if exist(filename,'file')
      txt=sprintf('Loading image %s ...',filename);
      disp(txt);
      image=imread(filename);
      if lparam.invert(slice)==1
        disp('Inverting...');
        image=double(255-image)/255;
      else
        image=double(image)/255;
      end;
      name=sprintf(param.baserawname,slice-startslice+lparam.targetfirstslice); 
      filename=sprintf('%s%s',param.rawdir,name);
      txt=sprintf('Writing image %s ...',filename);
      disp(txt);
      imwrite(image,filename,lparam.targetfileformat);
    else
      txt=sprintf('Warning: image %s not found!',filename);
      disp(txt);
    end;
  end;
else
  for slice=startslice:1:endslice
    for row=1:1:param.nrofrows
      for column=1:1:param.nrofcolumns
        name=sprintf(param.baseorigname,lparam.sourceslice(slice),row,column); filename=sprintf('%s%s',param.origdir,name);
        if exist(filename,'file')
          txt=sprintf('Loading image %s ...',filename);
          disp(txt);
          image=imread(filename);
          if lparam.invert(slice)==1
            disp('Inverting...');
            image=double(255-image)/255;
          else
            image=double(image)/255;
          end;
          name=sprintf(param.baserawname,slice-startslice+lparam.targetfirstslice,row,column);
          filename=sprintf('%s%s',param.rawdir,name);
          txt=sprintf('Writing image %s ...',filename);
          disp(txt);
          imwrite(image,filename,lparam.targetfileformat);
        else
          txt=sprintf('Warning: image %s not found!',filename);
          disp(txt);
        end;
      end;
    end;
  end;
end;
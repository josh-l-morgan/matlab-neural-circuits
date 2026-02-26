function pipeline_writescanfaultfile(param)

lparam=param.writescanfaultfile;

if (param.nrofrows==1)&&(param.nrofcolumns==1)
  fid = fopen(lparam.filename, 'w+');
  for slice=1:1:size(param.checkimages.scanerrorpos,2) %param.nrofslices
    nrofscanfaults=max(size(param.checkimages.scanerrorpos{slice}));
    if nrofscanfaults>0
      name=sprintf(param.baserawname,slice+param.firstslice-1);
      if (nrofscanfaults==1)
        fprintf(fid, '%d scan error  found in image %s, at:',nrofscanfaults,name);
      else
        fprintf(fid, '%d scan errors found in image %s, at:',nrofscanfaults,name);
      end;
      scanfaultline=param.checkimages.scanerrorpos{slice};
      for i=1:1:nrofscanfaults
        fprintf(fid,' %d',scanfaultline(i));
      end;
      fprintf(fid,'\r\n');
    end;
  end;
  fclose(fid);
else
  fid = fopen(lparam.filename, 'w+');
  for slice=1:1:size(param.checkimages.scanerrorpos,1) %param.nrofslices
    for row=1:1:size(param.checkimages.scanerrorpos,2) %param.nrofrows
      for column=1:1:size(param.checkimages.scanerrorpos,3) %param.nrofcolumns
        nrofscanfaults=max(size(param.checkimages.scanerrorpos{slice,row,column}));
        if nrofscanfaults>0
          name=sprintf(param.baserawname,slice+param.firstslice-1,row,column);
          if (nrofscanfaults==1)
            fprintf(fid, '%d scan error  found in image %s, at:',nrofscanfaults,name);
          else
            fprintf(fid, '%d scan errors found in image %s, at:',nrofscanfaults,name);
          end;
          scanfaultline=param.checkimages.scanerrorpos{slice,row,column};
          for i=1:1:nrofscanfaults
            fprintf(fid,' %d',scanfaultline(i));
          end;
          fprintf(fid,'\r\n');
        end;
      end;
    end;
  end;
fclose(fid);
end;
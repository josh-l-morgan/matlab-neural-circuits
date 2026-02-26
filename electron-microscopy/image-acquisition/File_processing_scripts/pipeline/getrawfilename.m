function filename=getrawfilename(param,slice,row,column)

if (param.nrofrows==1)&&(param.nrofcolumns==1)
  name=sprintf(param.baserawname,slice);
else
  name=sprintf(param.baserawname,slice,row,column);
end;
filename=sprintf('%s%s',param.rawdir,name);
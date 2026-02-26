function filename=filenameonly(filenamewithpath)
%This utility function returns only a filename if it is given a filename with path (using /)
%By Daniel Berger for MIT-BCE Seung / Harverd Lichtman, March 2010

pos=size(filenamewithpath,2);

%scan from end of string backwards until first occurence of '/'
while (pos>0)&&(filenamewithpath(pos)~='/')
  pos=pos-1;
end;

filename=filenamewithpath(pos+1:end);
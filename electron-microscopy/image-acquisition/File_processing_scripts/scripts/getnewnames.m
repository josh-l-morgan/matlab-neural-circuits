function filelist=getnewnames(newfilelist,oldfilelist)
%This function checks for each filename in newfilelist whether it is 
%in oldfilelist and returns only the filenames which were not found.
%By Daniel Berger for MIT-BCS Seung / Harvard Lichtman, March 2010

nrofnewfiles=size(newfilelist,2);
nrofoldfiles=size(oldfilelist,2);

filelist={};

for nf=1:1:nrofnewfiles
  found=0;
  nffn=newfilelist{nf};
  for of=1:1:nrofoldfiles
    offn=oldfilelist{of};
    if strcmp(nffn,offn)==1
      found=1; 
      break;
    end;
  end;
  if found==0
    filelist=[filelist nffn];
  end;
end;
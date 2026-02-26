function image=generatesineweightimage(nrofrows,nrofcolumns)
%This function generates an image for sine weighting for image blending
%Please note that the inputs are (Y,X) and not (X,Y)!
%By Daniel Berger for MIT-BCS Seung, June 16th 2009

%generate a sine row and a sine column
sinerow=zeros(1,nrofcolumns);
for column=1:1:nrofcolumns
  sinerow(1,column)=sin((column-1)*pi/(nrofcolumns-1));
end;
sinecolumn=zeros(nrofrows,1);
for row=1:1:nrofrows
  sinecolumn(row,1)=sin((row-1)*pi/(nrofrows-1));
end;

%Multiply sine rows and columns
image=zeros(nrofrows,nrofcolumns);
for row=1:1:nrofrows
  for column=1:1:nrofcolumns
    image(row,column)=sinerow(1,column)*sinecolumn(row,1);
  end;
end;
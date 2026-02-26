function [cpsx,cpsy,cptx,cpty]=computeaffine_iterative_getcplist_nooutliers(slicecount,nrofrows,nrofcolumns,gridwidth,gridheight, as,arow,acolumn)
%computeaffine_iterative_getcplist.m
%A function for computeaffine_iterative 
%which extracts the list of corresponding points 
%from the 
%By Daniel Berger for MIT-BCS Seung, May 2 2009

global cp_gsx cp_gsy cp_gtx cp_gty;

%count number of corresponding point pairs connecting this tile to its
%neighbors
pcount=0;

%with previous slice
if as>1
  for row=arow-1:arow+1
    for column=acolumn-1:acolumn+1
      if (row>0) && (row<=nrofrows) && (column>0) && (column<=nrofcolumns)
        pcount=pcount+sum(sum((cp_gtx(as,arow,acolumn,1,row,column,:,:)~=0)&(cp_gty(as,arow,acolumn,1,row,column,:,:)~=0)));
        pcount=pcount+sum(sum((cp_gtx(as-1,row,column,3,arow,acolumn,:,:)~=0)&(cp_gty(as-1,row,column,3,arow,acolumn,:,:)~=0)));
      end;
    end;
  end;
end;
% %with same slice
for row=arow-1:arow+1
  for column=acolumn-1:acolumn+1
    if (row>0) && (row<=nrofrows) && (column>0) && (column<=nrofcolumns) && ((row~=arow) || (column~=acolumn))
      pcount=pcount+sum(sum((cp_gtx(as,arow,acolumn,2,row,acolumn,:,:)~=0)&(cp_gty(as,arow,acolumn,2,row,column,:,:)~=0)));
      pcount=pcount+sum(sum((cp_gtx(as,row,column,2,arow,acolumn,:,:)~=0)&(cp_gty(as,row,column,2,arow,acolumn,:,:)~=0)));
    end;
  end;
end;
%with next slice
if as<slicecount
  for row=arow-1:arow+1
    for column=acolumn-1:acolumn+1
      if (row>0) && (row<=nrofrows) && (column>0) && (column<=nrofcolumns)
        pcount=pcount+sum(sum((cp_gtx(as,arow,acolumn,3,row,column,:,:)~=0)&(cp_gty(as,arow,acolumn,3,row,column,:,:)~=0)));
        pcount=pcount+sum(sum((cp_gtx(as+1,row,column,1,arow,acolumn,:,:)~=0)&(cp_gty(as+1,row,column,1,arow,acolumn,:,:)~=0)));
      end;
    end;
  end;
end;

%Later, the Affine update is computed so that A*(cpsx,cpsy,1)' = (cptx,cpty,1);
%this means cpsx,cpsy are at (as,arow,acolumn) and cptx,cpty around it

cpsx=zeros(1,pcount);
cpsy=zeros(1,pcount);
cptx=zeros(1,pcount);
cpty=zeros(1,pcount);
pc=1;

%with previous slice
if as>1
  for row=arow-1:arow+1
    for column=acolumn-1:acolumn+1
      if (row>0) && (row<=nrofrows) && (column>0) && (column<=nrofcolumns)
        for gy=1:1:gridheight
          for gx=1:1:gridwidth
            if (cp_gtx(as,arow,acolumn,1,row,column,gy,gx)~=0)&&(cp_gty(as,arow,acolumn,1,row,column,gy,gx)~=0)
              cpsx(pc)=cp_gsx(as,arow,acolumn,1,row,column,gy,gx);
              cpsy(pc)=cp_gsy(as,arow,acolumn,1,row,column,gy,gx);
              cptx(pc)=cp_gtx(as,arow,acolumn,1,row,column,gy,gx);
              cpty(pc)=cp_gty(as,arow,acolumn,1,row,column,gy,gx);
              pc=pc+1;
            end;
          end;
        end;
        for gy=1:1:gridheight
          for gx=1:1:gridwidth
            if (cp_gtx(as-1,row,column,3,arow,acolumn,gy,gx)~=0)&&(cp_gty(as-1,row,column,3,arow,acolumn,gy,gx)~=0)
              cptx(pc)=cp_gsx(as-1,row,column,3,arow,acolumn,gy,gx);
              cpty(pc)=cp_gsy(as-1,row,column,3,arow,acolumn,gy,gx);
              cpsx(pc)=cp_gtx(as-1,row,column,3,arow,acolumn,gy,gx);
              cpsy(pc)=cp_gty(as-1,row,column,3,arow,acolumn,gy,gx);
              pc=pc+1;
            end;
          end;
        end;
      end;
    end;
  end;
end;

%with same slice
for row=arow-1:arow+1
  for column=acolumn-1:acolumn+1
    if (row>0) && (row<=nrofrows) && (column>0) && (column<=nrofcolumns) && ((row~=arow) || (column~=acolumn))
      for gy=1:1:gridheight
        for gx=1:1:gridwidth
          if (cp_gtx(as,arow,acolumn,2,row,column,gy,gx)~=0)&&(cp_gty(as,arow,acolumn,2,row,column,gy,gx)~=0)
            cpsx(pc)=cp_gsx(as,arow,acolumn,2,row,column,gy,gx);
            cpsy(pc)=cp_gsy(as,arow,acolumn,2,row,column,gy,gx);
            cptx(pc)=cp_gtx(as,arow,acolumn,2,row,column,gy,gx);
            cpty(pc)=cp_gty(as,arow,acolumn,2,row,column,gy,gx);
            pc=pc+1;
          end;
        end;
      end;
      for gy=1:1:gridheight
        for gx=1:1:gridwidth
          if (cp_gtx(as,row,column,2,arow,acolumn,gy,gx)~=0)&&(cp_gty(as,row,column,2,arow,acolumn,gy,gx)~=0)
            cptx(pc)=cp_gsx(as,row,column,2,arow,acolumn,gy,gx);
            cpty(pc)=cp_gsy(as,row,column,2,arow,acolumn,gy,gx);
            cpsx(pc)=cp_gtx(as,row,column,2,arow,acolumn,gy,gx);
            cpsy(pc)=cp_gty(as,row,column,2,arow,acolumn,gy,gx);
            pc=pc+1;
          end;
        end;
      end;
    end;
  end;
end;

%with next slice
if as<slicecount
  for row=arow-1:arow+1
    for column=acolumn-1:acolumn+1
      if (row>0) && (row<=nrofrows) && (column>0) && (column<=nrofcolumns)
        for gy=1:1:gridheight
          for gx=1:1:gridwidth
            if (cp_gtx(as,arow,acolumn,3,row,column,gy,gx)~=0)&&(cp_gty(as,arow,acolumn,3,row,column,gy,gx)~=0)
              cpsx(pc)=cp_gsx(as,arow,acolumn,3,row,column,gy,gx);
              cpsy(pc)=cp_gsy(as,arow,acolumn,3,row,column,gy,gx);
              cptx(pc)=cp_gtx(as,arow,acolumn,3,row,column,gy,gx);
              cpty(pc)=cp_gty(as,arow,acolumn,3,row,column,gy,gx);
              pc=pc+1;
            end;
          end;
        end;
        for gy=1:1:gridheight
          for gx=1:1:gridwidth
            if (cp_gtx(as+1,row,column,1,arow,acolumn,gy,gx)~=0)&&(cp_gty(as+1,row,column,1,arow,acolumn,gy,gx)~=0)
              cptx(pc)=cp_gsx(as+1,row,column,1,arow,acolumn,gy,gx);
              cpty(pc)=cp_gsy(as+1,row,column,1,arow,acolumn,gy,gx);
              cpsx(pc)=cp_gtx(as+1,row,column,1,arow,acolumn,gy,gx);
              cpsy(pc)=cp_gty(as+1,row,column,1,arow,acolumn,gy,gx);
              pc=pc+1;
            end;
          end;
        end;
      end;
    end;
  end;
end;

%outlier detection and removal.
dist=(cptx-cpsx).*(cptx-cpsx)+(cpty-cpsy).*(cpty-cpsy);
vd=(dist<3*std(dist));

cpsx=cpsx(vd);
cpsy=cpsy(vd);
cptx=cptx(vd);
cpty=cpty(vd);
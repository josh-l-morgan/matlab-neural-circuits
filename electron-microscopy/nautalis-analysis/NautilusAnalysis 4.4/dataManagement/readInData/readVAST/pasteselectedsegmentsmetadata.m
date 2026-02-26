%pasteselectedsegmentsmetadata.m
%Pastes metadata stored in variable 'copydata' to selected segmentation layer in VAST
%By Daniel Berger for VAST, January 2020

global vdata;
if (max(size(vdata))==0)
  disp('WARNING: Not connected to VAST. Please first run VAST, enable the API, then VastTools.m and connect.');
  %return;
end;

if (~exist('copydata','var'))
  disp('ERROR: Variable "copydata" not found!');
  return;
end;
if (~strcmp(class(copydata),'struct'))
  disp('ERROR: "segnames" is not a struct!');
  return;
end;
[selectedlayernr, selectedemlayernr, selectedsegmentlayernr, res] = vdata.vast.getselectedlayernr();
if (selectedsegmentlayernr==-1)
  disp('ERROR: No segmentation layer loaded in VAST!');
  return;
end;

%Append root segment of copied data to selected segmentation layer
rootrow=find(copydata.segdata(:,1)==copydata.rootnr);

%Find selected segment in selected segment layer and add new root after it
[segdata, res] = vdata.vast.getallsegmentdatamatrix();
isselected=bitand(bitshift(segdata(:,2),-16),3);
torefid=find(isselected==1);

tpid=addonesegment(torefid,0,copydata.segnames{rootrow},copydata.segdata(rootrow,:));
tparentidlist=[0];
sparentidlist=[0];

csid=copydata.segdata(rootrow,15);
nextorchild=1; trefid=tpid; srefid=copydata.rootnr;
while (csid~=0)
  csrow=find(copydata.segdata(:,1)==csid);
  if (min(size(csrow))==0)
    %There were segments following in the source but they are not stored in copydata. end paste
    csid=0;
  else
    %tpid=tparentidlist(end);
    tid=addonesegment(trefid,nextorchild,copydata.segnames{csrow},copydata.segdata(csrow,:));
    if (nextorchild==1)
      tparentidlist=[tparentidlist trefid];
      sparentidlist=[sparentidlist srefid];
    end;
    if (copydata.segdata(csrow,15)~=0)
      %This segment has a child. Go to child.
      csid=copydata.segdata(csrow,15);
      nextorchild=1;
      trefid=tparentidlist(end);
    else
      if (copydata.segdata(csrow,17)~=0)
        %This segment has no child but a next. Go to next.
        trefid=tid;
        csid=copydata.segdata(csrow,17);
        nextorchild=0;
      else
        %This segment has no child and no next. Recurse up
        found=0;
        while (found==0)
          srefid=sparentidlist(end);
          sparentidlist=sparentidlist(1:end-1);
          ctid=tparentidlist(end);
          tparentidlist=tparentidlist(1:end-1);
          if (srefid==0)
            found=1;
          else
            %check if csid has next
            csrow=find(copydata.segdata(:,1)==srefid);
            if (copydata.segdata(csrow,17)~=0)
              trefid=ctid;
              csid=copydata.segdata(csrow,17);
              found=1;
            end;
          end;
        end;
      end;
    end;
  end;
end;


function id=addonesegment(refid,nextorchild,name,data)
global vdata;
  id=vdata.vast.addsegment(refid,nextorchild,name);
  vdata.vast.setanchorpoint(id,data(11),data(12),data(13));
  vdata.vast.setsegmentcolor8(id,data(3),data(4),data(5),data(6),data(7),data(8),data(9),data(10));
end
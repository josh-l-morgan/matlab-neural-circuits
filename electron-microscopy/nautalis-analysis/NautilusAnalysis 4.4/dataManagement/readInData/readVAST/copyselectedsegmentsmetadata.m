%copyselectedsegmentsmetadata.m
%Copies metadata of selected segment and children to Matlab, stored in variable 'copydata'
%By Daniel Berger for VAST, January 2020

global vdata;
if (max(size(vdata))==0)
  disp('WARNING: Not connected to VAST. Please first run VAST, enable the API, then VastTools.m and connect.');
  %return;
end;

[selectedlayernr, selectedemlayernr, selectedsegmentlayernr, res] = vdata.vast.getselectedlayernr();
[linfo, res] = vdata.vast.getlayerinfo(selectedsegmentlayernr);
[segdata, res] = vdata.vast.getallsegmentdatamatrix();
[segnames, res] = vdata.vast.getallsegmentnames();
segnames(1)=[];
isselected=bitand(bitshift(segdata(:,2),-16),3);
copydata.segdata=segdata(isselected>0,:);
copydata.segnames=segnames(isselected>0);
copydata.rootnr=find(isselected==1);
disp(sprintf('Copied metadata of %d selected segments from selected segmentation layer',size(copydata.segdata,1)));
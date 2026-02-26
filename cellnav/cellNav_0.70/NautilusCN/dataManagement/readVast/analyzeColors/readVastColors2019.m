function[objInfo] = readVastColor2019(vdata);

global vdata

objInfo.info = vdata.vast.getinfo();
segData = vdata.vast.getallsegmentdata

for i  = 2:length(segData)
    
    objInfo.ids(i-1) = segData{i}.id;
    objInfo.names{i-1} = vdata.vast.getsegmentname(objInfo.ids(i-1));
    objInfo.anchors(i-1,:) = segData{i}.anchorpoint;%([2 1 3]);
    objInfo.flags(i-1,:) = segData{i}.flags;
    objInfo.boundBox(i-1,:) = segData{i}.boundingbox;
    objInfo.col1(i-1,:) = segData{i}.col1;
    objInfo.col2(i-1,:) = segData{i}.col2;
    objInfo.parent(i-1) = segData{i}.hierarchy(1);
    objInfo.children(i-1) = segData{i}.hierarchy(2);
    objInfo.hierarchy{i-1} = segData{i}.hierarchy;
    %objInfo.parentName(i-1) = segData{i}.hierarchy(1);
    
end





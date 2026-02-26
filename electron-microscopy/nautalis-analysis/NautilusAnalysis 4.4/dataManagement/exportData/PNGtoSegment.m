% function[rleDat] = getRLEtoSubs(TPN,miplevel);
%
% %Please connect to VAST with vasttools first!
%
% global vdata;
%
% tic

vast=VASTControlClass();
res = vast.connect('127.0.0.1',22081,1000)

%%
SPN = GetMyDir;
segName = [SPN 'import.vsv'];

%% makeImportFile
if 0 %below is script necessary to fetch color file info
    %if ~exist([SPN 'colData.mat'], 'file')
    [nr, res] = vast.getnumberofsegments();
    colData = vast.getallsegmentdata;
    names = vast.getallsegmentnames;
    for s = 1 : length(names)
        colData{s}.name = names{s};
    end
    save([SPN 'colData.mat'],'colData');
end

load([SPN 'colData.mat']);



if ~exist([SPN 'prog.mat'],'file')
    
    dSPN = dir([SPN '*.png'])
    cSPN = dir([SPN '*.txt'])
    inams = {dSPN.name};
    cnams = {dSPN.name};
    
    prog.SPN = SPN;
    prog.inams = inams;
    prog.inum = length(inams);
    prog.read = zeros(prog.inum,1);
    prog.imported = zeros(prog.inum,1);
    
    imageInfo = imfinfo([SPN inams{1}]);
    prog.bitDepth = imageInfo.BitDepth;
    prog.width = imageInfo.Width;
    prog.height = imageInfo.Height;
   
    save([SPN 'prog.mat'],'prog')
    
else
    load([SPN 'prog.mat'])
end

%% Create color

%%Make all segment names
for i = 1: length(colData)
    [nr, res] = vast.getnumberofsegments();
    dat = colData{i};
    id  = dat.id;
    
    if (nr -1) >= id
        res = vast.setsegmentname(id,dat.name);
    else       
        res = vast.addsegment(nr-1, 0, dat.name);
    end    
    if ~res
        disp('failed to make segment name')
    end
end
   

%%Set segment properties
for i = 1:length(colData)
    dat = colData{i};
    id  = dat.id;
      
    res = vast.setsegmentcolor32(id,dat.col1,dat.col2)
    res = vast.setanchorpoint(id, dat.anchorpoint(1),...
        dat.anchorpoint(2), dat.anchorpoint(3));
    res = vast.setsegmentbbox(id, dat.boundingbox(1),...
        dat.boundingbox(2), dat.boundingbox(3), dat.boundingbox(4),...
        dat.boundingbox(5), dat.boundingbox(6));
end

%%Move parent, first child, previous, next 0 if there is none.
%%res = movesegment(id, refid, nextorchild)
H = zeros(length(colData),5);
ids = zeros(length(colData),1);
for i = 1:length(colData)
    H(i,1:4) = colData{i}.hierarchy;
    H(i,5) = colData{i}.id;
end
H(1,1) = inf;

parents = unique(H(2:end,1));
for i = 1:length(parents)
    p = parents(i);
    moveS = find(H(:,1)==p);
    Hs = H(moveS,:);
    lastS = find(Hs(:,3)==0);
        
    if p>0
       res = vast.movesegment(Hs(lastS,5), p, 1);
    end
    
    for m = 1:length(moveS)-1
        newMove = find(Hs(:,3) == Hs(lastS,5));
        res = vast.movesegment(Hs(newMove,5), Hs(lastS,5), 0);
        lastS = newMove;
    end
    
end

if 0 % skripts to use
    res = setselectedlayernr(obj);
    res = setsegmentcolor8(id,r1,g1,b1,p1,r2,g2,b2,p2);
    res = setsegimageraw(miplevel,minx, maxx, miny, maxy, minz, maxz, segimage);
    res = setsegimageRLE(miplevel, minx, maxx, miny, maxy, minz, maxz, segimage);
    res = setviewcoordinates(x,y,z);
    res = setviewzoom(zoom);
end


%% Read in Data

toImport = find(prog.imported == 0);
Itimes = [];
for i = 1:length(toImport)
    Itic = tic;
    isize = prog.width;
    nam = prog.inams{toImport(i)};
    
    I = imread([prog.SPN nam]);
    prog.read(toImport(i)) = 1;
    
    %%parse name
    und = regexp(nam,'_');
    dots = regexp(nam,'.png');
    r = str2num(nam(und(1)+2:und(2)-1));
    c = str2num(nam(und(2)+2:und(3)-1));
    s = str2num(nam(und(3)+2:dots(1)-1));
    x = [isize * (c -1)  isize * (c -1) + isize - 1];
    y = [isize * (r -1)  isize * (r -1) + isize - 1];
    
    %I3 = cat(3,I,I,I);
    
    res = vast.setsegimageRLE(0,x(1), x(2), y(1), y(2), s, s, uint16(I)')
    %res = vast.setsegimageraw(0,x(1), x(2), y(1), y(2), s, s, uint16(I)');
    if res;  prog.imported(toImport(i)) = 1;  end
   
    %%Change view
    if 1
        res = vast.setviewcoordinates(round(mean(x)),round(mean(y)),round(mean(s)));
        res = vast.setviewzoom(-32);
    end
    
    Itoc = toc(Itic)/60;
    Itimes = [Itimes Itoc];
    meanTime = mean(Itimes);
    timeLeft = meanTime * (length(toImport)-i) / 60 ;
    disp(sprintf('%d of %d, Imin = %0.2f, hours left = %0.2f', i,length(toImport),Itoc,timeLeft))
    save([SPN 'prog.mat'],'prog')

end


%res = vast.savelayer(1,segName, 1)
















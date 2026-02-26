function[wif] = GetMyWafer

WPN = GetMyDir;  %waffer path name

%% Get Waffer folder information
clear wif
dWPN = dir(WPN); dWPN = dWPN(3:end);
SPN = {}; xmlNam = {}; logNam = {}; secDir = {};%create cell array of section names
for i = 1:length(dWPN)
    nam = dWPN(i).name;
    if isdir([WPN nam])
        SPN{length(SPN)+1,1} = nam;
        secDir{length(secDir)+1,1} = [WPN nam];
        q = strfind(nam,'w');
        und = strfind(nam,'_s');
        if ~isempty(q) & ~isempty(und)
        secPos(1) = str2num(nam(q(1)+2:und(1)-1));
        secPos(2) = str2num(nam(und(1)+3:und(2)-1));
        secPos(3) = str2num(nam(und(2)+4:end));
        else
            nums = [];
            for n = 1:length(nam)
                nums = [nums str2num(nam(n))]
            end
            secPos = nums;
        end
        sid(length(secDir),:) = secPos;
    elseif regexp(nam,'.log')
        logNam{length(logNam)+1,1} = [WPN nam];
    elseif regexp(nam, '.xml');
        xmlNam{length(xmlNam)+1,1} = [WPN nam];
    end
end

wif.dir = WPN;
wif.secNam = SPN;
wif.secDir = secDir;
wif.log = logNam;
wif.xml = xmlNam;  %waffer xml
wif.secID = sid; %numeric position of section (waffer,strip,section)
%% Get section xml name

for i = 1:length(secDir)
    dSec = dir(secDir{i}); dSec = dSec(3:end);
    tile = {}; 
    for s = 1:length(dSec);
        nam = dSec(s).name;
        if regexp(nam,'Tile_')
            tile{length(tile)+1,1} = [secDir{i} '\' nam];
            und = regexp(nam,'_');
            dash = regexp(nam,'-');
            rc(length(tile),:) = [str2num(nam(und(1)+2:dash(1)-1)) ...
                str2num(nam(dash(1)+2:und(2)-1))];
        elseif regexp(nam,'MosaicInfo')
           sec.xml =  [secDir{i} '\' nam]; 
        elseif regexp(nam,'StageStitchedOverview')
           sec.ov =  [secDir{i} '\' nam];
        elseif regexp(nam,'SiteLocation')
           sec.site =  [secDir{i} '\' nam];       
        end
    end

%     %%XML too slow    
%     [tree, rootname, dom]=xml_read(sec.xml); 
%     tiles = tree.Tiles.Tile;
%     for s = 1:length(tiles)
%         tile{s} = tiles(s).Filename;
%         rc(s,:) = [tiles(s).ATTRIBUTE.row tiles(s).ATTRIBUTE.col];
%     end
 
    %%Sort tiles
    maxC = max(rc(:,2));
    sortRC = rc(:,1) * maxC + rc(:,2);
    [sorted order] = sort(sortRC,'ascend');
    tile = tile(order);
    rc = rc(order,:);
    
    %%Record
    sec.tile = tile;   
    sec.rc = rc;
    wif.sec(i)= sec;
end
    

%% Get image info

inf = imfinfo(wif.sec(1).tile{1});
wif.imfinfo = inf;






    
    
%function[userEval] = manFocus(WPN,tile);

if ~(exist('WPN','var'))
    WPN = GetMyDir;
end

if exist([WPN 'mif.mat'],'file')
    load([WPN 'mif.mat'])
else
    mif = getMif(WPN);
end

w = 1;
fig = gcf;
sid = 1:length(mif.w(w).sections);
colormap gray(256)
xs = mif.Info.Width;
ys = mif.Info.Height;
q = 0;
readI = 1; readS = 1;
maxRes = 400;

brightness = -90;
contrast = 1.5;
button = 1;
t = 1; tNew = 1;
isize = 400;
qualScale = 10;

maxSize = mif.Info.Height;
subReg = [round(xs/2 - isize/2) round(xs/2 - isize/2) + isize-1];
subRegX = subReg; subRegY = subReg;
s = 1;
secNam = mif.w.sec(s);
maxSec = length(mif.w.sec);

numTile = length(mif.w.sec(s).tileNams);
good = zeros(1,numTile); rate = good;
I = 255-imread([WPN mif.w.sec(s).secFold '\' mif.w.sec(s).tileNams{t}],'PixelRegion',{subRegY,subRegX});

%% Make dummy mos

cNum = max(mif.w(w).sec(s).colIds);
rNum = max(mif.w(w).sec(s).rowIds);
qualSamp = zeros(rNum,cNum);



%imread([WPN mif.sec(s).name mif.sec(s).ov]);  %Get overview image
%downSamp = imread([mif.dir(1:end-1) 'shaped\downSamp\' mif.secNam{s} '.tif']);
tileFold = mif.w(w).sec(s).tileFolders{1};
slashes = regexp(tileFold,'\');
secNam = tileFold(slashes(end)+1:end-8)
if isfield(mif.w(w).sec(s),'ov')
    ov = [WPN mif.w(w).sec(s).secFold '\' mif.w(w).sec(s).ov];
    downSamp = imread(ov);
else
    downSamp = mif.w(w).sec(s).idMos;
end



%% Get images of mosic quality

if exist([WPN 'qual.mat'])
    load([WPN 'qual.mat'])
else
    qual = q2qual(TPN);
end
qualSamp = qual.w(w).sec(s).tile * qualScale;

% qual = wQual.sec(s);
% for i = 1:length(wQual.mos)
%    mosA(:,:,i) = wQual.mos{i};
% end
% showMos = uint8((mosA-2)*20);
% datDims = size(showMos);
% rNum = datDims(1);
% cNum = datDims(2);
% if length(datDims)>2
%
%     sNum = datDims(3);
% else
%     sNum = 1;
% end
% tileID = (sNum-1)*cNum*rNum + ((rNum-1)*cNum)+cNum;
% rL = sNum * cNum * rNum;
% for c = 1:3; %scale colors
%     vals = showMos(:,:,:,c);
%     showMos(:,:,:,c) = 255 * (showMos(:,:,:,c)-min(vals(:)))/(max(vals(:))/min(vals(:)));
% end
% showMos = uint8(showMos);
% qualSamp = showMos(:,:,s,:);
% % qualSamp = reshape(qualSamp,rNum,cNum,3);
% qualSamp = showMos(:,:,s);
% qualSamp = uint8(cat(3,qualSamp,qualSamp,qualSamp));
%
% [blankY blankX] = find(ones(rNum,cNum));

%%Image
downSampAx = subplot(2,3,1);
image(downSamp), hold on
mosAx = subplot(2,3,4);
image(qualSamp)
imageAx = subplot(2,3,[2 3 5 6]);
image((I + brightness) * contrast),hold off


%% get previous quality values
% if exist([mif.dir 'userEval.mat'],'file')
%     load([mif.dir 'userEval.mat'])
%     if isfield(userEval,'qualVals')
%         qualVals = userEval.qualVals;
%     else
%         qualVals = zeros(rNum,cNum,sNum);
%     end
% else
%     qualVals = zeros(rNum,cNum,sNum);
% end
% [rates fails thresh] = findThresh(mosA,qualVals);
%
% %qualVals = zeros(cNum,rNum,sNum);
% %Display
% downSampAx = subplot(2,3,1);
% image(downSamp)
% mosAx = subplot(2,3,4);
% image(qualSamp), hold on
% scatter(blankX,blankY,'.','SizeData',500,'MarkerEdgeColor',[.1 .1 .1])
% scatter(mif.sec(s).tile(t).col,mif.sec(s).tile(t).row);
% failed = fails{s};
% rated = rates{s};
% if ~isempty(rated),scatter(rated(:,2),rated(:,1),'*','g'),end
% if ~isempty(failed),scatter(failed(:,2),failed(:,1),'x','r'),end
% %scatter(rL(rL == 1,1)+.4,rL(rL==2,2),'.','g')
%
%
% rc = [[mif.sec(s).tile.row]' [mif.sec(s).tile.col]'];
rc = [mif.w(w).sec(s).rowIds' mif.w(w).sec(s).colIds'];
idMos = mif.w(w).sec(s).idMos;

%% Start



while ~q  % While not quitting
    figure(fig)
    mouseOK = 1;
    
    %pause(.03)
    [pressed Ax] = getMyIn;
    
    figure(fig)
    
    
    if isstruct(pressed) %was there a keypress?
        key = pressed.Key
        changeC = 0;
        if strcmp(key,'w'),      subRegY = subRegY - size(I,1)/2; readI = 1;
        elseif strcmp(key,'s'),  subRegY = subRegY + size(I,1)/2; readI = 1;
        elseif strcmp(key,'a'),  subRegX = subRegX - size(I,2)/2; readI = 1;
        elseif strcmp(key,'d'),  subRegX = subRegX + size(I,2)/2; readI = 1;
        elseif strcmp(key,'q') | strcmp(key,'leftbracket'), %zoom out
            mY = mean(subRegY); dif =diff(subRegY);
            subRegY = round([mY-dif mY+dif]);
            mX = mean(subRegX); dif =diff(subRegX);
            subRegX = round([mX-dif mX+dif]);
            readI = 1;
        elseif strcmp(key,'e') | strcmp(key,'rightbracket'), %zoom out
            mY = mean(subRegY); dif =diff(subRegY);
            subRegY = round([mY-dif/2.5 mY+dif/4]);
            mX = mean(subRegX); dif =diff(subRegX);
            subRegX = round([mX-dif/2.5 mX+dif/4]);
            readI = 1;
        elseif strcmp(key,'c'),subRegX = subReg; subRegY = subReg; readI = 1; %reset zoom
        elseif strcmp(key,'f'),readI = 1; %point to quality mosaic
            [c r] = ginput(1);
            t = find((rc(:,1) ==round(r))&(rc(:,2)==round(c)));
        elseif strcmp(key,'downarrow'), brightness = brightness - 10;
        elseif strcmp(key,'uparrow'), brightness = brightness + 10;
        elseif strcmp(key,'rightarrow'), contrast = contrast + .25;
        elseif strcmp(key,'leftarrow'), contrast = contrast - .25;
        elseif strcmp(key,'z'),s = s-1; %move back one section
            if s<1;s = 1; end;   readI = 1; readS = 1;
            rc = [mif.w(w).sec(s).rowIds' mif.w(w).sec(s).colIds'];
            idMos = mif.w(w).sec(s).idMos;
        elseif strcmp(key,'x'),s = s+1; %move forward one section
            if s>maxSec, s = maxSec; end; readI = 1; readS = 1;
            rc = [mif.w(w).sec(s).rowIds' mif.w(w).sec(s).colIds'];
            idMos = mif.w(w).sec(s).idMos;
        elseif strcmp(key,'space'), good(t) = 1;'good',readI = 1;
            userEval.sec(s).tile(t).retake = 0;
            rL(sub2ind([sNum rNum cNum],s,rc(t,1),rc(t,2))) = 0;
        elseif strcmp(key,'shift'), good(t) = 2;'bad', readI = 1;
            userEval.sec(s).tile(t).retake = 1;
            rL(sub2ind([sNum rNum cNum],s,rc(t,1),rc(t,2))) = 1;
        elseif ~isempty(str2num(key)),  readS = 1;
            rating = str2double(key);
            sprintf('you rated the image quality a %d out of 4',rating)
            qualVals(rc(t,1),rc(t,2),s) = rating
            userEval.qualVals = qualVals;
            userEval.sec(s).tile(s).rating = rating;
            [rates fails thresh] = findThresh(mosA,qualVals);
            userEval.thresh = thresh;
            userEval.rates = rates;
            userEval.fails = fails;
            save([mif.dir 'userEval.mat'],'userEval');
        elseif strcmp(key,'tab'),
            escape = verify;
            if escape
                break
            end
            
        elseif strcmp(key,'backslash'), readS
            qualVals = zeros(cNum,rNum,sNum)
        elseif strcmp(key,'capslock') | strcmp(key,'control') ,'saving'
            
        else
            
        end
        
        if readS
            rc = [mif.w(w).sec(s).rowIds' mif.w(w).sec(s).colIds'];
            if isfield(mif.w(w).sec(s),'ov')
                ov = [WPN mif.w(w).sec(s).secFold '\' mif.w(w).sec(s).ov];
                downSamp = imread(ov);
            else
                downSamp = mif.w(w).sec(s).idMos;
            end
            subplot(2,3,1),image(downSamp),pause(.01)
            
            %             [rates fails thresh] = findThresh(mosA,qualVals);
            qualSamp = qual.w(w).sec(s).tile * qualScale;
            
            %qualSamp = reshape(qualSamp,rNum,cNum);
            subplot(2,3,4),image(qualSamp),hold on
            
            %             scatter(blankX,blankY,'.','SizeData',500,'MarkerEdgeColor',[.1 .1 .1])
            %             scatter(mif.sec(s).tile(t).col,mif.sec(s).tile(t).row,'k');
            %             scatter(mif.sec(s).tile(t).col,mif.sec(s).tile(t).row,'.','w');
            %             failed = fails{s};
            %             rated = rates{s};
            % if ~isempty(rated),scatter(rated(:,2),rated(:,1),'*','g'),end
            % if ~isempty(failed),scatter(failed(:,2),failed(:,1),'x','r'),end
            %             pause(.01)
            %             hold off
            
        end
        
        
        if  readI
            
            subplot(2,3,[2 3 5 6])
            if t<1,t=1;end
            if t> numTile, t = numTile;end
            subRegY = round(subRegY);subRegX = round(subRegX);
            if subRegY(1)<0; subRegY = subRegY-subRegY(1)+1;   end
            if subRegX(1)<0; subRegX = subRegX-subRegX(1)+1;   end
            if subRegY(2)>maxSize; subRegY = subRegY+(maxSize - subRegY(2));   end
            if subRegX(2)>maxSize; subRegX = subRegX+(maxSize - subRegX(2));   end
            if diff(subRegY)>maxSize; subRegY(1) = 1; subRegX(1)=1; end
            
            
            rsize = diff(subRegY);
            if rsize>maxRes
                useRegY = [subRegY(1) fix(rsize/maxRes)+1 subRegY(2)];
                useRegX = [subRegX(1) fix(rsize/maxRes)+1 subRegX(2)];
            else
                useRegY = subRegY; useRegX = subRegX;
            end
            I = 255-imread([WPN mif.w(w).sec(s).secFold '\' mif.w(w).sec(s).tileNams{t}],...
                'PixelRegion',{useRegY,useRegX});
            readI = 0;
            
        end
        
        subplot(2,3,[2 3 5 6])
        image((I + brightness) * contrast)
        pause(.01)
        
    else
        
        if length(pressed)>1
            inP = round(pressed(1,:));
            if Ax == mosAx %if clicked on mosaic
                newt = find((rc(:,1) ==inP(2))&(rc(:,2)==inP(1)));
            elseif Ax == downSampAx
                c = fix(cNum * inP(1)/size(downSamp,2))+1;
                r = fix(rNum * inP(2)/size(downSamp,1))+1;
                newt = find((rc(:,1) ==r)&(rc(:,2)==c));
                
            else
                button = get(gcf, 'SelectionType');
                if ~strcmp(button,'open')
                    butt = button;
                end
                if strcmp(butt,'normal') %change colors
                    newt = t-1;
                elseif strcmp(butt,'alt') %edit labels
                    newt = t+1;
                end
            end%which axis clicked
        end
        if (newt>0) & (newt<=numTile)
            t = newt;
        end
        
        rsize = diff(subRegY);
        if rsize>maxRes
            useRegY = [subRegY(1) fix(rsize/maxRes)+1 subRegY(2)];
            useRegX = [subRegX(1) fix(rsize/maxRes)+1 subRegX(2)];
        else
            useRegY = subRegY; useRegX = subRegX;
        end
        
        rc = [mif.w(w).sec(s).rowIds' mif.w(w).sec(s).colIds'];
        if isfield(mif.w(w).sec(s),'ov')
            ov = [WPN mif.w(w).sec(s).secFold '\' mif.w(w).sec(s).ov];
            downSamp = imread(ov);
        else
            downSamp = mif.w(w).sec(s).idMos;
        end
        
        
        qualSamp = qual.w(w).sec(s).tile * qualScale;
        %         qualSamp = reshape(qualSamp,rNum,cNum);
        subplot(2,3,4),image(qualSamp),hold on
        %         scatter(blankX,blankY,'.','SizeData',500,'MarkerEdgeColor',[.1 .1 .1])
        %         scatter(mif.sec(s).tile(t).col,mif.sec(s).tile(t).row,'k');
        %         scatter(mif.sec(s).tile(t).col,mif.sec(s).tile(t).row,'.','w');
        %         failed = fails{s};
        %         rated = rates{s};
        % if ~isempty(rated),scatter(rated(:,2),rated(:,1),'*','g'),end
        % if ~isempty(failed),scatter(failed(:,2),failed(:,1),'x','r'),end
        %         pause(.01)
        %         hold off
        
        subplot(2,3,[2 3 5 6])
        I = 255-imread([WPN mif.w(w).sec(s).secFold '\' mif.w(w).sec(s).tileNams{t}],...
            'PixelRegion',{useRegY,useRegX});
        image((I + brightness) * contrast ) %.3 sec
        %        use = wQual.sec(s).tile(t).quality
        %use = focSec(s).tile(t).use
        pause(.01)
        
        %}
        
    end
end






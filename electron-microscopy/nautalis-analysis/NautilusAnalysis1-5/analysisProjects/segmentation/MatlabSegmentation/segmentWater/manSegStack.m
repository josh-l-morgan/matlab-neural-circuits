
%% Matlab 2d segmentor
%%-Josh Morgan

%%Notes
%{
    fix undo
    check save
%}
clear all
[TFN TPN] = GetMyFile;  %Get Tif
I = imread([TPN TFN]);
%I = I(600:900,600:900);

%Defaults
colorRange = 30;
brightness = 0;
contrast = .7;

memSize = 3; %number of undos available
if exist([TPN 'fused\fused' TFN])
    lI = double(imread([TPN 'fused\fused' TFN]));
    'fused'       
elseif exist([TPN 'edited\ed' TFN])
    lI = double(imread([TPN 'edited\ed' TFN]));
    'ed'
elseif exist([TPN 'labeled\lab' TFN])
    'lab'
    lI =double(imread([TPN 'labeled\lab' TFN]));
    if ~exist([TPN 'edited']),mkdir([TPN 'edited']); end
else
    'No labels found'
end



[ys, xs] = size(I);
I3 = cat(3,I,I,I); %Color image
subplot(1,2,1)
image(I)
subplot(1,2,2)

%% Remember changes
memLab = zeros(xs,ys,memSize,'uint16');
memLab(:,:,memSize) = lI;
currentMem = memSize;
lIb = lI; %back up copy of lI for single undo
%% Label
subplot(1,1,1)

for i = 1:1
    myCol = hsv(max(lI(:))+1)*colorRange;

    [r rix] = sort(rand(size(myCol,1)),1);
    myCol= myCol(rix,:);
    myCol = cat(1,[0 0 0],myCol);
    red = myCol(:,1); green = myCol(:,2); blue = myCol(:,3);
    skipCol = 5 + round(rand*10);

    colO = uint8(cat(3,red(lI+1),green(lI+1),blue(lI+1)));
    image((I3 + brightness) * contrast + colO),pause(.01)
    xax = get(gca,'xlim'); yax = get(gca,'ylim');
end



%% Get mouse inputs
%%Left click fuse, right click cut

%%MakeCuttingdisk
sRad = 3;
diskE = gaus3d([sRad * 2 + 1,sRad * 2 + 1,1],2);
diskE = diskE>=diskE(sRad,1);
pDisk = bwperim(diskE);
diskE(pDisk) = 0;
[dy dx] = find(diskE);
[py px] = find(pDisk);
dy = dy -sRad -1; dx = dx - sRad - 1;
py = py - sRad - 1; px = px - sRad -1;


%% Draw on figure
q = 0;
while ~q
    directions = 'e = edit, s = save, q = quit, rl = reload starting labels : '
    inp = input(directions,'s');
    inp = lower(inp);
    
    if strcmp(inp,'s')
        imwrite(uint16(lI),[TPN 'edited\ed' TFN],'Compression','none');
    elseif strcmp(inp,'q')
        q = 1;
    elseif strcmp(inp,'rl')
        lI =double(imread([TPN 'labeled\lab' TFN]));
    end
    if inp == 'e' | inp == 'rl'
        'editing...'
        
        image((I3 + brightness) * contrast + colO),pause(.01)
        % w = waitforbuttonpress
        % inp = CurrentCharacter


        flipCol = uint8(1);
        escapeRequest = 0;

        while ~escapeRequest
            xax = get(gca,'xlim'); yax = get(gca,'ylim'); %update axis vars
            %[mx my button] = ginput;
            pts = getLine;pause(.01)

            if isstruct(pts) %was there a keypress?
                key = lower(pts.Key)
                changeC = 0;
                if strcmp(key,'w'),     yax = yax - diff(yax)/5; ylim(yax);
                elseif strcmp(key,'s'), yax = yax + diff(yax)/5; ylim(yax);
                elseif strcmp(key,'a'), xax = xax - diff(xax)/5; xlim(xax);
                elseif strcmp(key,'d'), xax = xax + diff(xax)/5; xlim(xax);
                elseif strcmp(key,'q'), %zoom out
                    xdiff = diff(xax); ydiff = diff(yax);
                    xax(1) = xax(1) - xdiff/5; yax(1) = yax(1) - ydiff/5;
                    xax(2) = xax(2) + xdiff/5;  yax(2) = yax(2) + ydiff/5; 
                    xlim(xax); ylim(yax);              
                elseif strcmp(key,'e'), %zoom out
                    xdiff = diff(xax);  ydiff = diff(yax);
                    xax(1) = xax(1) + xdiff/5; yax(1) = yax(1) + ydiff/5;
                    xax(2) = xax(2) - xdiff/5; yax(2) = yax(2) - ydiff/5;
                    xlim(xax);  ylim(yax);
                elseif strcmp(key,'space'), flipCol = uint8(~flipCol)
%                     xax = get(gca,'xlim'); yax = get(gca,'ylim');
%                     image((I3 + brightness) * contrast + colO * flipCol) %.3 sec
%                     xlim([xax]);  ylim([yax]);,pause(.001)
                elseif strcmp(key,'q'), colorRange = colorRange - 5; changeC = 1;
                elseif strcmp(key,'e'), colorRange = colorRange + 5; changeC = 1;
                elseif strcmp(key,'c'), changeC = 1;
                elseif strcmp(key,'tab'),escapeRequest = 1;
                elseif strcmp(key,'downarrow'), brightness = brightness - 10;
                elseif strcmp(key,'uparrow'), brightness = brightness + 10;
                elseif strcmp(key,'rightarrow'), contrast = contrast + .1;
                elseif strcmp(key,'leftarrow'), contrast = contrast - .1;
                elseif strcmp(key,'space'), flipCol = ~flipCol; changeC = 1;
                elseif strcmp(key,'home'), colorRange = colorRange - 5; changeC = 1;
                elseif strcmp(key,'pageup'), colorRange = colorRange + 5; changeC = 1;
                elseif strcmp(key,'c'), changeC = 1;
                elseif strcmp(key,'tab'),escapeRequest = 1;
                elseif strcmp(key,'capslock'),'saving'
                   imwrite(uint16(lI),[TPN 'edited\ed' TFN],'Compression','none');
                    'saved'
                elseif strcmp(key,'r'), currentMem = max(1,currentMem -1);
                    %lI = memLab(:,:,currentMem);
                    lIc = lI; lI = lIb; lIb = lIc; %ugly juggle
                    colO = uint8(cat(3,red(lI+1),green(lI+1),blue(lI+1)));
                elseif strcmp(key,'f'), currentMem = min(memSize,currentMem+1);
                    lI = memLab(:,:,currentMem);
                else
                    
                end
                if changeC
                    myCol = hsv(max(lI(:))+1) * colorRange;
                    [r rix] = sort(rand(size(myCol,1)),1);
                    myCol= myCol(rix,:);
                    myCol = cat(1,[0 0 0],myCol);
                    red = myCol(:,1); green = myCol(:,2); blue = myCol(:,3);
                    skipCol = 5 + round(rand*10);
                    colO = uint8(cat(3,red(lI+1),green(lI+1),blue(lI+1)));
                end
                xax = get(gca,'xlim'); yax = get(gca,'ylim');
                image((I3 + brightness) * contrast + colO * flipCol)
                xlim([xax]);  ylim([yax]);
                pause(.01)

            else
                button = get(gcf, 'SelectionType')
                if strcmp(button,'open')
                    button = 2
                elseif strcmp(button,'normal')
                    button = 1;
                elseif strcmp(button,'extend')
                    button = 0;
                elseif strcmp(button,'alt')
                    button = 3;
                else
                    error('MATLAB:ginput:InvalidSelection', 'Invalid mouse selection.')
                end

                %% Respond to buttons
                if button == 2 %return to keyboard
                    escapeRequest = 1;
                elseif button == 0 %change colors
                    flipCol = ~flipCol
                    myCol = hsv(max(lI(:))+1) * 30 * flipCol;
                    [r rix] = sort(rand(size(myCol,1)),1);
                    myCol= myCol(rix,:);
                    myCol = cat(1,[0 0 0],myCol);
                    red = myCol(:,1); green = myCol(:,2); blue = myCol(:,3);
                    skipCol = 5 + round(rand*10);
                elseif size(pts,2)==2 %edit labels
                    mx = pts(:,1);my = pts(:,2); 
                    lIb = lI;
                    lI = changeIDs(lI,mx,my,button,dy,dx,py,px);
                    memLab(:,:,1:end-1)=memLab(:,:,2:end); %store changes
                    memLab(:,:,end) = lI; 
                    if ~isempty(button)
                        escapeRequest = button(1)==2;
                    end
                end

               
                tic
                colO = uint8(cat(3,red(lI+1),green(lI+1),blue(lI+1))); %.4 sec
                toc
                tic
                toc
                xax = get(gca,'xlim'); yax = get(gca,'ylim');
                toc
                image((I3 + brightness) * contrast + colO*flipCol) %.3 sec
                toc
                xlim([xax]);  ylim([yax]);
                pause(.01)
                toc
                %}
                
            end
        end
    end
end
%% Save file

imwrite(uint16(lI),[TPN 'edited\ed' TFN],'Compression','none');



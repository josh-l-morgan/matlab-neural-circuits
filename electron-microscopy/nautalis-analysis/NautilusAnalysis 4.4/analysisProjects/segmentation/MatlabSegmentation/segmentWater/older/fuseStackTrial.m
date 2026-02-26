clear all
TPN = GetMyDir;

Inams = getPics([TPN]);
if exist([TPN 'edited'])
    labelDir = 'edited\';
elseif exist([TPN 'labeled'])
    labelDir = 'labeled\';
else
    sprintf('No labels found.')
end

inams = getPics([TPN labelDir]);  %find all tifs

if ~exist([TPN 'fused'])
    mkdir([TPN 'fused']);
end



I1 = imread([TPN labelDir inams{1}]);
Iraw1 = imread([TPN Inams{1}]);
I1 = I1(1:400,1:400);
Iraw1 = Iraw1(1:400,1:400);
rProps =regionprops(I1,'Area');
Areas(:,1) = [rProps.Area]';
Areas(:,2) = Areas(:,1) * 0 + 1;

myCol = hsv(double(max(I1(:)))+1) * 220;
[r rix] = sort(rand(size(myCol,1)),1);
myCol= myCol(rix,:);
myCol = cat(1,[0 0 0],myCol);
red = myCol(:,1); green = myCol(:,2); blue = myCol(:,3);
colO = uint8(cat(3,red(I1+1),green(I1+1),blue(I1+1)));
image(colO)

id16 = 1:65535;
id16 = setdiff(id16,I1(:));

mergeThresh = .4;
for i = 2:length(inams)
    imwrite(uint16(I1),[TPN 'fused\fused' inams{i}],'Compression','none');
    storeI(:,:,i) = I1;
    memIds = I1;
    sprintf('running plane %d of %d',i,length(inams))
    maxI1 = max(I1(:));
    sprintf('Ids are at %d',maxI1)
    I2 = imread([TPN labelDir inams{i}]);
    I2 = I2(1:400,1:400);
    Iraw2 = imread([TPN Inams{i}]);
    Iraw2 = Iraw2(1:400,1:400);
    newI1 = I1;
    newI2 = I2;
    newI2(newI2>0) = id16(newI2(newI2>0));

    for r = 1:1
        I2props = regionprops(newI2,'Area');
        iAreas = [I2props.Area];
        allIds = find(iAreas>0);
        iAreas = iAreas(iAreas>0);
        %    [iAreas idx] = sort([I2props.Area],'ascend');

        for a = 1: length(allIds)
            %sprintf('running area %d of %d',a,length(allIds))
            overlapped = newI1(newI2 == allIds(a));
            overlapped = overlapped(overlapped>0);
            if ~isempty(overlapped)
                [mO F] = mode(double(overlapped)); %apply threshold
                %             F,iAreas(a)
                %             F >= (iAreas(a) * mergeThresh)
                %             pause
                if F>= (iAreas(a) * mergeThresh)
                    newI2(newI2==allIds(a)) = mO;
                end
            end
%             
%                     subplot(2,2,1)
%                     colO = uint8(cat(3,red(newI1+1),green(newI1+1),blue(newI1+1)));
%                     image(colO)
%                     subplot(2,2,2)
%                     colO = uint8(cat(3,red(newI2+1),green(newI2+1),blue(newI2+1)));
%                     bright = uint8((newI2 == allIds(a)) * 1000);
%                     image(colO+ cat(3,bright,bright,bright))
%                     subplot(2,2,3)
%                     colO = uint8(cat(3,red(newI2+1),green(newI2+1),blue(newI2+1)));
%                     image(colO)
%                     pause

        end  % end run merge forward

        %% run murge backward
        I1props = regionprops(newI1,'Area');
        iAreas = [I1props.Area];
        allIds = find(iAreas>0);
        iAreas = iAreas(iAreas>0);
        %    [iAreas idx] = sort([I2props.Area],'ascend');

        for a = 1: length(allIds)
            %sprintf('running area %d of %d',a,length(allIds))
            overlapped = newI2(newI1 == allIds(a));
            overlapped = overlapped(overlapped>0);
            if ~isempty(overlapped)
                [mO F] = mode(double(overlapped)); %apply threshold
                %             F,iAreas(a)
                %             F >= (iAreas(a) * mergeThresh)
                %             pause
                if F>= (iAreas(a) * mergeThresh)
                    newI1(newI1==allIds(a)) = mO;
                    newI2(newI2 == allIds(a)) = mO;
                end
            end
        end

        %I2(I2>0) = I2(I2>0) + max(I2(:));
        %newI2(newI2 == 0) = I2(newI2 == 0);

        %% Reintegrate with branch points.
        %% Conserve IDs 


        %% Display
        subplot(2,2,1)
        colO = uint8(cat(3,red(newI1+1),green(newI1+1),blue(newI1+1)));
        image(colO)
        subplot(2,2,2)
        colO = uint8(cat(3,red(newI2+1),green(newI2+1),blue(newI2+1)));
        image(colO)

        subplot(2,2,3)
        image(Iraw1)
        subplot(2,2,4)
        image(Iraw2)
        colormap gray(256)
        pause(.01)
    end

    id16 = setdiff(id16,[newI2(:); newI1(:)]);

    I1 = newI2;
    Iraw1 = Iraw2;


end
imwrite(uint16(I1),[TPN 'fused\fused' inams{i}],'Compression','none');


%%

inams = getPics([TPN 'fused']);  %find all tifs
subplot(1,1,1)
cmap = hsv(256);
cmap(1,:) = 0;
colormap(cmap)




for i = 1: 4%length(inams)
    I = imread([TPN 'fused\' inams{i}]);
    Imod = mod(I,256);
    Imod(~I) = 0;
    subplot(2,2,i)
    image(Imod)
    pause(1)
end

% %%  PlotSkel
% %Make all Seg
% 'Display results'
% subplot(1,1,1)
% hold off
% plot(1,1)
% whitebg('k')
% 'render cell'
% uCol = hsv(100);
% skipdist = 11 + round(rand  * 8);
% for w = 2:max(storeI(:))
%     wInd = storeI == w;
%     if length(wInd) >10
%         subplot(1,1,1)
%         w
%         'rendering'
%     p = patch(isosurface(storeI==w,.1));
%     set(p,'FaceColor',uCol(mod(w*skipdist,100)+1,:),'EdgeColor','none');
%     daspect([1/4 1/4 1/30])
%     alpha(.8)
%     view(3); axis tight
%     camlight
%     lighting gouraud
%     'rendered'
%     pause
%     end
% end

%%
cmap = hsv(256);
cmap(1,:) = 0;
colormap(cmap)
while 1
    for i = 1:length(storeI)
        image(mod(storeI(:,:,i),256)),pause(.01)
    end
end

%% Possible methods
%{
Thresholds
%}

function[] = tt_oif2Stack_MultiChan(app,CPN);

global globTT


try
    if strcmp(globTT.dir.TFN(end-3:end),'.oif')
        CPN = [globTT.dir.TPN globTT.dir.TFN '.files\'];
    else
        CPN = globTT.dir.TPN;
    end
catch
    CPN = globTT.dir.TPN;

end

zOne = app.z1Button.Value;
shouldWrite = app.toggle_shouldWrite.Value;
shouldMedFilt = 1;
shouldShow = 0;

TPN = [CPN(1:end-1) '_Processed\'];

%stkChan = {[3 2 1] [5 4  2]};

stkChan = {[app.edit_chanA1.Value app.edit_chanA2.Value app.edit_chanA3.Value];
    [app.edit_chanB1.Value app.edit_chanB2.Value app.edit_chanB3.Value]};

if zOne
    useSlices = [1];%[1:34];
else
    useSlices = [];
end

if shouldWrite
    mkdir(TPN)
end


%% Read and parse names
dCPN = dir([CPN '*.tif'])
inams = {dCPN.name};
clear color slice
if length(inams)>0
for i  = 1:length(inams)
    nam = inams{i};
    c = regexp(nam,'C');
    d = regexp(nam,'.tif');
    z = regexp(nam,'Z');
    
    if isempty(z)
        color(i) = str2num(nam(c+1:d-1));
        slice(i) = 1;
    else
        color(i) = str2num(nam(c+1:z-1));
        slice(i) = str2num(nam(z+1:d-1));
    end
end



%% Choose slices
if isempty(useSlices)
    slices = unique(slice);
else
    slices = useSlices;
end


%% Read as ..
if 0
    for s = 1:length(stkChan)
        TPN_I = sprintf('%sI_%d\\',TPN,s);
        if exist(TPN_I,'dir'), rmdir(TPN_I,'s'); end
        mkdir(TPN_I)
        for i = 1:length(slices)
            
            clear Ic
            targ = find((slice == slices(i)) & ( color == stkChan{s}(1)));
            if ~isempty(targ)
                I = imread([CPN inams{targ}]);
                I = shouldMedFilt2(I,[3 3]);
                I = I/256;
                Ic(:,:,1) = uint8(I);
            end
            
            targ = find((slice == slices(i)) & ( color ==  stkChan{s}(2)));
            if ~isempty(targ)
                I = imread([CPN inams{targ}]);
                I = I/256;
                I = shouldMedFilt2(I,[3 3]);
                Ic(:,:,2) = I;
            end
            
            targ = find((slice == slices(i)) & ( color ==  stkChan{s}(3)));
            if ~isempty(targ)
                I = imread([CPN inams{targ}]);
                I = shouldMedFilt2(I,[3 3]);
                I = I/256;
                Ic(:,:,3) = I;
            end
            if size(Ic,3)<3
                Ic(1,1,3) = 0;
            end
            
            image(uint8(Ic*100))
            i
            pause
            pause(.01)
            imwrite(uint8(Ic),sprintf('%s%05.0f.tif',TPN_I,slices(i)))
            
        end
    end
    
end


%% Read as ..
if 1
    testI = imread([CPN inams{1}]);
    [ys xs zs] = size(testI);
    for s = 1:length(stkChan)
        TPN_I = sprintf('%sI_%d\\',TPN,s);
        if exist(TPN_I,'dir'), rmdir(TPN_I,'s'); end
        mkdir(TPN_I)
        
        %% Select Display Tab
        if s == 1
            rawAx = app.rawAAxes;
            filtAx = app.filtAAxes;
            maxAx = app.maxAAxes;
            app.mainTab.SelectedTab = app.RawATab;
            
        elseif s == 2
            rawAx = app.rawBAxes;
            filtAx = app.filtBAxes;
            maxAx = app.maxBAxes;
            app.mainTab.SelectedTab = app.RawATab;
        end
        
        %% Make Stacks
        Ic1 = zeros(ys,xs,length(slices),'double');
        Ic2 = Ic1; Ic3 = Ic1;
        for i = 1:length(slices)
            
            txt = sprintf('reading mix %d of %d, slice %d of %d',s,length(stkChan),i,length(slices));
            app.textOut.Value = txt;
            pause(.001)
            
            targ = find((slice == slices(i)) & ( color == stkChan{s}(1)));
            if ~isempty(targ)
                I = double(imread([CPN inams{targ}]));
                if shouldMedFilt, I = medfilt2(I,[3 3]); end
                Ic1(:,:,i) = I;
            end
            
            targ = find((slice == slices(i)) & ( color ==  stkChan{s}(2)));
            if ~isempty(targ)
                I = double(imread([CPN inams{targ}]));
                if shouldMedFilt, I = medfilt2(I,[3 3]); end
                Ic2(:,:,i) = I;
            end
            
            targ = find((slice == slices(i)) & ( color ==  stkChan{s}(3)));
            if ~isempty(targ)
                I = double(imread([CPN inams{targ}]));
                if shouldMedFilt, I = medfilt2(I,[3 3]); end
                Ic3(:,:,i) = I;
            end
        end
        
        Ic1 = Ic1 * 256/max(max(Ic1(:)),1);
        Ic2 = Ic2 * 256/max(max(Ic2(:)),1);
        Ic3 = Ic3 * 256/max(max(Ic3(:)),1);
        
        Imax1 = max(Ic1,[],3);
        Imax2 = max(Ic2,[],3);
        Imax3 = max(Ic3,[],3);
        
        Imax = (cat(3,Imax1,Imax2,Imax3));
        
        Imean1 = mean(Ic1,3);
        Imean2 = mean(Ic2,3);
        Imean3 = mean(Ic3,3);
        
        Imean = cat(3,Imean1,Imean2,Imean3);
        
        if shouldWrite
            imwrite(uint8(Imax),sprintf('%smax_stk%d.tif',TPN,s))
        end
        image(maxAx,uint8(Imax*256/max(Imax(:))));
        
        %% Combine, write Display stacks
        TPN_I = sprintf('%sI_%d\\',TPN,s);
        if exist(TPN_I,'dir'), rmdir(TPN_I,'s'); end
        mkdir(TPN_I)
        for i = 1:length(slices)
            txt = sprintf('building display stacks, slice %d of %d',i,length(slices));
            app.textOut.Value = txt;
            drawnow
            IcSlice = cat(3,Ic1(:,:,i),Ic2(:,:,i),Ic3(:,:,i));
            Ic(:,:,:,i) = IcSlice;
            if shouldShow | (i == 1)
                image(rawAx,uint8(IcSlice))
                image(filtAx,uint8(IcSlice*(300/max(IcSlice(:)))))
                pause(.01)
            end
            if shouldWrite
                imwrite(uint8(Ic(:,:,:,i)),sprintf('%s%05.0f.tif',TPN_I,slices(i)))
            end
        end
        
        if s == 1
             globTT.I.tab{1} = Ic;
            globTT.I.tab{2} = Ic;
            globTT.I.tab{3} = Imax;
        elseif s == 2
            globTT.I.tab{4} = Ic;
            globTT.I.tab{5} = Ic;
            globTT.I.tab{6} = Imax;
        end
        
        
    end
    
    for i = 1:length(globTT.I.tab);
        
        globTT.I.max(i) = max(globTT.I.tab{i}(:));
        hRange = [0:globTT.I.max(i)];
        clear cHist
        for c = 1:3
            vals = globTT.I.tab{i}(:,:,c,:);
            cHist(:,c) = hist(vals(:),hRange);
        end        
        globTT.I.hist{i} = cHist;
        globTT.I.hRange{i} = hRange;
    end
    
    
    
end
else
    app.textOut.Value = 'No images found';
    pause(3)
end












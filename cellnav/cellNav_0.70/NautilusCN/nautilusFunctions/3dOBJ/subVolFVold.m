function[fv] = subVolFV(obSubRaw,saveName,renderProps)

if ~exist('renderProps','var')
    renderProps = [];
end

%% make vol
buf = 5;
obSub = round(obSubRaw);
% obSub(:,1) = obSub(:,1)- min(obSub(:,1))+1;
% obSub(:,2) = obSub(:,2)- min(obSub(:,2))+1;
% obSub(:,3) = obSub(:,3)- min(obSub(:,3))+1;

minSub = floor(min(obSub,[],1));


obSub = obSub - repmat(minSub,[size(obSub,1) 1])+buf;
obSiz = max(obSub,[],1)+buf;

minVec = minSub;
bufVec = buf;

obInd = sub2ind(obSiz,obSub(:,1),obSub(:,2),obSub(:,3));
%% make vol
obVol = zeros(obSiz,'logical');
obVol(obInd) = 1;



%% dilate
if isfield(renderProps,'dilate')
    if renderProps.dilate>0
        SE = strel('sphere',renderProps.dilate);
        obVol = imdilate(obVol,SE);
    end
end

%% smooth
if isfield(renderProps,'smooth')
    totVol = sum(obVol(:)>0);
    oldInd = find(obVol>0);
    if renderProps.smooth>0
        %kern = 2*ceil(renderProps.smooth)+1;
        for r = 1:renderProps.smooth
            obVol = smooth3(obVol,'gaussian',7,0.7);
        end
    end
    obVol = obVol>.08;
    
end

%% scale
scale = 1;
if isfield(renderProps,'resize')
    if renderProps.resize~=1
        scale = renderProps.resize;
        obVol = imresize3(double(obVol),scale,'method','cubic');
        minVec = minSub * scale;
        bufVec = buf * scale;
    end
end

%% smooth after scale
if isfield(renderProps,'smooth')
    totVol = sum(obVol(:)>0);
    oldInd = find(obVol>0);
    if renderProps.smooth>0
        %kern = 2*ceil(renderProps.smooth)+1;
        for r = 1:renderProps.smooth
            obVol = smooth3(obVol,'gaussian',7,.7);
        end
    end
    obVol = obVol>.8;
end

%%

fv=isosurface(obVol,0.5);
if ~exist('renderOb','var')
    renderOb = 0;
end

if ~isempty(fv.vertices)
    fv.vertices = fv.vertices + repmat(minVec([2 1 3]) ,[size(fv.vertices,1) 1])-bufVec;
    
    
    
    if isfield(renderProps,'smoothPatch')
        if renderProps.smoothPatch>0
            sMode = 1;
            lambda = .8;
            sigma = 1;
            itt = ceil(renderProps.smoothPatch);
            fv=smoothpatch(fv,sMode,itt,lambda,sigma);
        end
    end
    
    if isfield(renderProps,'renderOb')
        
        if renderProps.renderOb
            
            p = patch(fv);
            view(30,-15);
            axis vis3d;
            colormap copper
            set(p,'FaceColor','red','EdgeColor','none');
            daspect([1,1,1])
            view(3); axis tight
            camlight
            lighting gouraud
        end
    end
end
%%


%%
if exist('fileName','var')
    objectname='hi';
    disp(['Saving ' filename ' as Wavefront OBJ.....']);
    vertface2obj(v,f,filename,objectname);
end





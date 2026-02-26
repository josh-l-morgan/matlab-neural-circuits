function[] = showSub3D(subs,renderProps);

if ~isempty(subs)
    if size(subs,2) == 3
        
        downSamp = 4;
        renderProps.smooth = 0;
        renderProps.resize = 1;
        renderProps.smoothPatch = 0;
        
        smallSub = shrinkSub(subs,downSamp);
        %smallSub = smallSub(:,flipDim);
        
        fv = subVolFV(smallSub,[],renderProps);
        %fv.vertices = fv.vertices * downSamp;
        fv.vertices = fv.vertices(:,[2 1 3]) * downSamp;
        
        [p] = renderFV(fv,renderProps.col,renderProps.alph);
        viewAngle = [0 0];
        view(viewAngle)
        axis off
        lightAngle = [150 150];
        lightangle(lightAngle(1),lightAngle(2)) ;
        
    end
end
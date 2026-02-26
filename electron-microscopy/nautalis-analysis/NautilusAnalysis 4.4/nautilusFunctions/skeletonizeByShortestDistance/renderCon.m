function[fv] = renderCon(subs,con2,col,alph)

%subs = dropSubs(subs);
%%

if ~isempty(subs)
    
    if ~exist('col','var')
        col = [1 1 1];
    end
    
    if ~exist('alph','var')
        alph = .2;
    end
    
    if ~exist('con2','var') | isempty(con2)
        [subs ia ic] = uniqueSubs(subs);
        %col = col(ia,:);
        %alph = alph(ia);
        con2 = obj2con(subs);
    end
    con = con2(:,1:6);
    
    %% make cube of verts and faces organized by con position
    
    vertShift = [.5 .5 .5; .5 .5 -.5; .5 -.5 .5; .5 -.5 -.5];
    vertNeg = vertShift .* repmat([-1 1 1],[size(vertShift,1) 1]);
    oneFace = [1 2 3; 2 3 4];
    
    clear vert6
    vert6(:,:,1) = vertNeg(:,[2 3 1]);
    vert6(:,:,2) = vertNeg(:,[3 1 2]);
    vert6(:,:,3) = vertNeg;
    vert6(:,:,4) = vertShift;
    vert6(:,:,5) = vertShift(:,[3 1 2]);
    vert6(:,:,6) = vertShift(:,[2 3 1]);
    vert6 = permute(vert6,[3 2 1]); %shape like subs + corner
    
    %% surfaces
    surfInd = find(con==0);
    [surfY surfX] = ind2sub([size(con,1),size(con,2)],surfInd);
    sizeSurf = length(surfY);
    
    
    vert = vert6(surfX,:,:) + repmat(subs(surfY,:),[1 1 4]);
    vertSubs = cat(1,vert(:,:,1),vert(:,:,2),vert(:,:,3),vert(:,:,4));
    
    sL = [1:sizeSurf]';
    oL = ones(sizeSurf,1);
    
    face1 = zeros(sizeSurf,3);
    face2 = zeros(sizeSurf,3);
    for f = 1:3
        face1(:,f) = sub2ind([sizeSurf 4], sL, oL * oneFace(1,f));
        face2(:,f) = sub2ind([sizeSurf 4], sL, oL * oneFace(2,f));
    end
    faces = cat(1,face1,face2);
    
    clear vertSubs
    
    vertSubs = [];
    for i = 1:4
        vertSubs = cat(1,vertSubs,subs(surfY,:) + vert6(surfX,:,i));
    end
    
    axis off
    
    fv.vertices = vertSubs;
    fv.faces = faces;
    [p] = renderFV(fv,col,alph);
    
    
    
    p.FaceLighting = 'gouraud';
    p.AmbientStrength = 1;%.4; shadowColor
    p.DiffuseStrength = .8;%.1;
    p.SpecularStrength = .6;
    p.SpecularExponent = 3;
    p.BackFaceLighting = 'lit';
    p.Clipping = 'off';
    
    
    
    
end





































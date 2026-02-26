function[fv] = face4to3(fv)

%% input fv and change faces from a n by 4 patch to a n by 3 patch
if size(fv.faces,2)==4
    faces4 = fv.faces;
    faces3 = cat(1,faces4(:,[1 2 3]),faces4(:,[3 4 1]));
    fv.faces = faces3;
else
    'not n by 4 face'
end





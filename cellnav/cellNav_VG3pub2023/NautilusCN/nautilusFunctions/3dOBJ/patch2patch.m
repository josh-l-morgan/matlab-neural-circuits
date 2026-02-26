function[p2] = patch2patch(p)


p2 = patch('vertices',p.Vertices,'faces',p.Faces)

f = fields(p);

for i = 1:length(f)
    try
        setfield(p2,f{i},getfield(p,f{i}))
    catch err
    end
end
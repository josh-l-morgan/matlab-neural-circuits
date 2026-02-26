[y x vals] = find(connections)

preIn = [];
postIn = [];
for i = 1:length(y)
    for r = 1:vals(i)
        preIn(length(preIn)+1) = y(i);
        postIn(length(postIn)+1) = x(i);
    end
end
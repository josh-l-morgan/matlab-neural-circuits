function[allShared] = minSharedMat(con)

%%Find shared synapses in connectivity matrix (along first dimension)


for i = 1:size(con,1)
    for p = 1:size(con,1)
        if i~=p
            shared1(i,p) = sum(min(con(i,:),con(p,:)));
        end
    end
end


allShared = sum(shared1(:));




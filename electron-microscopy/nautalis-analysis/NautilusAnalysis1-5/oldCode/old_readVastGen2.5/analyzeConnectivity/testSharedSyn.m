function[] = testSharedSyn(con)

%%Find deviation from average number of shared synapses


for i = 1:size(con,1)
    for p = 1:size(con,1)
        if i~=p
            shared1(i,p) = sum(min(con(i,:),con(p,:)));
        end
    end
end

for i = 1:size(con,2)
    for p = 1:size(con,2)
        if i~=p
            shared2(i,p) = sum(min(con(:,i),con(:,p)));
        end
    end
end






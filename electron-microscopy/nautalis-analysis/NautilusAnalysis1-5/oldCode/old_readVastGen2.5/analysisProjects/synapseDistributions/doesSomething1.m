




for i = 1: 20

[postGraph3 postIDX] = reshuffleSimilarity(postGraph2);
collectPost{i} = postGraph3 * 15;

end


%%
for i = 1:length(collectPost)
    image(uint8(collectPost{i}))
    pause(.1)
end


GPNmany = [GPN 'many2\'];mkdir(GPNmany);
for i = 1:length(collectPost)

imwrite(uint8(collectPost{i}*10),[GPNmany 'sGraphPreType' num2str(i) '.png']);
end
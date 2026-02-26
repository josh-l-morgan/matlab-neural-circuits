%%Test for redundant synapses and shared synapses using simple shuffle of
%%edges


reps = 100;

for i = 1:reps

    newEdges = shuffleEdges(edges);
    con = edge2con(newEdges);
    sharedPre(i) = minSharedMat(con);
    sharedPost(i) = minSharedMat(con');
    redundant(i) = countRedundant(con);
    image(con),pause(.1);
    
end
    
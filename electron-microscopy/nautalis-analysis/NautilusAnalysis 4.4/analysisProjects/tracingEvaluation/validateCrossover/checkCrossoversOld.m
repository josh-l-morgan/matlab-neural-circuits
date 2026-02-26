
edges = obI.nameProps.edges(:,[2,1]);

pre108 = edges(:,2) == 108;
pre108 = unique(edges(pre108,1));


pre201 = edges(:,2) == 201;
pre201 = unique(edges(pre201,1));


pre170 = edges(:,2) == 170;
pre170 = unique(edges(pre170,1));


pre127 = edges(:,2) == 127;
pre127 = unique(edges(pre127,1));



cross201cross108 = intersect(pre108,pre201)

cross201cross170 = intersect(pre201,pre170)

cross108cross170 = intersect(pre108,pre170)

post2035 = edges(:,1) == 2035;
post2035 = unique(edges(post2035,2));


%%
% check 127 to 201 through 2035 and 127 to 108 through 1025
%   followed 2035 to relink 127 to 268 ((13563, 16045, 5016)) and 201
% 127 has glomeruli
% following 1025  split at (16672, 11071, 5940), traced to 107
% split (13888, 17848, 2715), traced to 108




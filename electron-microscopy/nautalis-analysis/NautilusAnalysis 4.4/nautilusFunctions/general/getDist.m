function[dist] = getDist(sub1,sub2)
%%creates distance matrix with first input along the first dimension and 
%%the second input along the second dimension

if size(sub1,2) == 3
dist = sqrt((sub1(:,1)-sub2(:,1)').^2 + (sub1(:,2)-sub2(:,2)').^2 + ...
    (sub1(:,3)-sub2(:,3)').^2);
else
    dist = sqrt((sub1(:,1)-sub2(:,1)').^2 + (sub1(:,2)-sub2(:,2)').^2);
end
    
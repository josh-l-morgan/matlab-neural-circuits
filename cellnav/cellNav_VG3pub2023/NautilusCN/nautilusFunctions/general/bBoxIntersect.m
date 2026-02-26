function[hit hits] = bBoxIntersect(box1,box2);


hits = zeros(size(box1,1), size(box2,1));
for b1 = 1:size(box1,1)
    for b2 = 1:size(box2,1)
        
        bx1 = box1(b1,:);
        bx2 = box2(b2,:);
        
        check1 = bx1(1:3)<bx2(4:6); % 1 too high
        check2 = bx2(1:3)<bx1(4:6); % 2 too high
        if (sum(check1) + sum(check2)) == 6
           hits(b1,b2) = 1; 
        else 
           hits(b1,b2) = 0;
        end
        
    end
end

hit = sum(hits(:));

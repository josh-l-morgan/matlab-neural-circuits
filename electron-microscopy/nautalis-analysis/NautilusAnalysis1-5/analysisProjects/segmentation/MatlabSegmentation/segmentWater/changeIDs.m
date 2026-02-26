function[lI] = chaneIDs(lI,mx, my, button,dy,dx,py,px)
%%Change Labelefield according to mouse output

escapeRequest = 0;

%[mx my button] = ginput;
myx = wall([my mx], size(lI));
mx = myx(:,2); my = myx(:,1);
lid = sub2ind(size(lI),round(my),round(mx));
labs = sort(lI(lid));
labs = labs(labs>0);
if (numel(button) >0 )
    if (button(1) ==3 ) ||(button(1) ==1 )
        %toOpen = imopen(toOpen,strel('disk',1));
        clearInd = []; grabInd = [];
        for s = 1:length(mx)-1;
            L = sqrt((mx(s+1)-mx(s)).^2 + (my(s+1)-my(s)).^2);
            if diff([mx(s+1) mx(s)]) % if there is some length
                step = (mx(s+1)-mx(s))/L * 3;
                lSlope = (my(s+1)-my(s))/(mx(s+1)-mx(s));
                lx = mx(s) :step: mx(s+1);
                ly = (lx-mx(s)) * lSlope + my(s);
            elseif diff([my(s+1) my(s)])
                step = (my(s+1)-my(s))/L * 3;
                lSlope = (mx(s+1)-mx(s))/(my(s+1)-my(s));
                ly = my(s) :step: my(s+1);
                lx = (ly-my(s)) * lSlope + mx(s);
            else %if just a point
                lx = mx(s);
                ly = my(s);
            end
            for l = 1: length(ly)
                lyx = wall(round([dy+ly(l) dx + lx(l)]),size(lI));
                clearInd = cat(1,clearInd,sub2ind(size(lI),lyx(:,1),lyx(:,2)));
                gyx = wall(round([py+ly(l) px + lx(l)]),size(lI));
                grabInd = cat(1,grabInd,sub2ind(size(lI),gyx(:,1),gyx(:,2)));

            end


        end
    end
    if ~isempty(grabInd)
        grabInd = grabInd(~isnan(grabInd));
        effected = unique(lI(grabInd));
        effected = effected(effected>0);

        if button(1) == 1
            for i = 2:length(effected) % color objects
                lI(lI == effected(i)) = effected(1);
            end
            %lI(toOpen >0) = toOpen(toOpen>0) + max(lI(:));
            lI(clearInd) = effected(1);
        else
            %%Break and label
            toOpen = lI * 0;

            toOpen(ismember(lI,effected)) = 1;
            %image(fitH(toOpen))
            toOpen(clearInd) = 0;
            toOpen = bwlabel(toOpen,4);
            Imax = max(lI(:));
            lI(toOpen >0) = uint16(toOpen(toOpen>0)) + Imax;
            lI(clearInd) = 0;
            if Imax > 60000
                'rearranging IDs';
                lI = uint16(bwlabel(lI,4));
            end
            

            
            
            
        end
    end

end

end %if button

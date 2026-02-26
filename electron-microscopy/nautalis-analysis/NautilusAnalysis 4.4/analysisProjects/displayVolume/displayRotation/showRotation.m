function[I_rot] = showRotation(viewProps,getDeg,movieDir)

    %%  Make rotation
    
    %getDeg = [0:60:360];
   % movieDir = 'D:\LGNs1\Analysis\movies\mixedClade\rot90_background_b\';
    
    if exist('movieDir','var')
        saveFile = 1;
            if ~exist(movieDir,'dir'),mkdir(movieDir),end

    else 
        saveFile = 0;
    end
    
    c = 0;
    for d = 1:length(getDeg)
        disp(sprintf('Rotating view %d of %d',d,length(getDeg)));
        %viewProps.degRot = getDeg(d);
        I = stereoCellsAndMoreFull(viewProps,getDeg(d));
        %I = uint8(I.^.7 * 4.5);
        image(uint8(I)),pause(.01)
        if saveFile
            fileName = sprintf('%srot_%05.0f_deg%d.png',movieDir,d, getDeg(d));
            imwrite(uint8(I),fileName)
                I_rot = I;

        else
            I_rot{d} = I;
        end
    end
%     
%     while 1
%         for d = 1:length(I_rot)
%             image(uint8(I_rot{d})),pause(.3)
%         end
%         for d = length(I_rot)-1:-1:2
%             image(uint8(I_rot{d})),pause(.3)
%         end
%         for d = 1:length(I_rot)
%             image(uint8(I_rot{d})),pause(.3)
%         end
%         for d = length(I_rot)-1:-1:2
%             image(uint8(I_rot{d})),pause(.3)
%         end
%         
%     end
    
    
function[] = playRot(I_rot)
    
    while 1
        for d = 1:length(I_rot)
            image(uint8(I_rot{d})),pause(.3)
        end
        for d = length(I_rot)-1:-1:2
            image(uint8(I_rot{d})),pause(.3)
        end
        for d = 1:length(I_rot)
            image(uint8(I_rot{d})),pause(.3)
        end
        for d = length(I_rot)-1:-1:2
            image(uint8(I_rot{d})),pause(.3)
        end
        
    end
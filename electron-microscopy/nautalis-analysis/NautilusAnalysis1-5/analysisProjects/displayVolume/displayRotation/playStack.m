function[] = playStack(I_stack)



 while 1
        for d = 1:size(I_stack,3)
            image(uint8(I_stack(:,:,d))),pause(.3)
        end
        for d = size(I_stack,3)-1:-1:2
            image(uint8(I_rI_stackot(:,:,d))),pause(.3)
        end
        for d = 1:size(I_stack,3)
            image(uint8(I_stack(:,:,d))),pause(.3)
        end
        for d = size(I_stack,3)-1:-1:2
            image(uint8(I_stack(:,:,d))),pause(.3)
        end
        
    end
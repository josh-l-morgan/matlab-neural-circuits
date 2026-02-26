

%% Set Block and knife

%%Block dimensions in nanometers
bY = 3000000;
bX = 1000000;
bZ = 10000000;
bF = ones(bY,1)*bZ; % blockFace

%%Block Position
pY = 0;
pZ = 0;

%%Knife Position
kY = 0;
kZset = bZ+100;

zstep = 35;

%% Set knife forces
flex = 1:100;
kSpring = flex.^3;
plot(kSpring)

%% Cut
for z = 1:1000
    bF = bF + zstep; % Step block forward
    kZ = kZset; % Return knife to starting position
    for y = 1 : bY ; % Start Cutting
        cutT = bF(1) - kZ; 
        secT(y) = cutT;
        flex = cutT/10;
        kSpring = flex^3;
        
    end % end y
    
end % end z

    
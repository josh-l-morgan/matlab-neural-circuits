clear
clc
uiopen
[m nExperiments] = size(OFFBSBIPHASIC);
g = waitbar(0, 'Experiments done...')
for k=1:nExperiments
    [m nCells] = size(OFFBSBIPHASIC(k).spikes);
    for i=1:nCells
        [z beta0 xy] = gausspreparation (OFFBSBIPHASIC(k).spikes(i).STA);
        gaussOptions = statset ('MaxIter', 10000);
        beta = nlinfit(xy, z, @gaussian2d, beta0, gaussOptions);
        %     beta = [amplitude xOffset yOffset theta sigmaX sigmaY
        %     background];
        %     ellipsedraw(beta(5), beta(6),beta(2), beta(3) , beta(4),'m')
        %     hold on
        OFFBSBIPHASIC(k).spikes(i).space.beta = beta;
    end

    SR = zeros(nCells);
    for i=1:nCells
        xI = OFFBSBIPHASIC(k).spikes(i).space.beta(2);
        yI = OFFBSBIPHASIC(k).spikes(i).space.beta(3);
        for j=1:nCells
            if j > i
                xJ = OFFBSBIPHASIC(k).spikes(j).space.beta(2);
                yJ = OFFBSBIPHASIC(k).spikes(j).space.beta(3);
                D = sqrt( (yJ-yI)^2 + (xJ-xI)^2 );
                slope = (yJ-yI)/(xJ-xI);
                X = linspace(0,80,1000);
                yzero = yJ - slope*xJ;
                yeighty = yJ - slope*(xJ-80);
                Y = linspace(yzero, yeighty, 1000);
                beta = OFFBSBIPHASIC(k).spikes(i).space.beta;
                [hEllipse XE YE]=...
                    ellipsedraw2(beta(5), beta(6),beta(2), beta(3) , beta(4),'m');
                nX = length(X);
                nXE = length(XE);
                E = zeros(nX, nXE);
                indVect = 1:nX*nXE;
                for d=1:nX
                    for e=1:nXE
                        E(d,e) = (X(d)-XE(e))^2 + (Y(d)-YE(e))^2;
                    end
                end
                IND = indVect(E==min(E(:)));
                [dc ec] = ind2sub([nX, nXE], IND);
                SI = sqrt( (yI-Y(dc))^2 + (xI-X(dc))^2 );

                beta = OFFBSBIPHASIC(k).spikes(j).space.beta;
                [hEllipse XE YE]=...
                    ellipsedraw2(beta(5), beta(6),beta(2), beta(3) , beta(4),'m');
                nXE = length(XE);
                E = zeros(nX, nXE);
                indVect = 1:nX*nXE;
                for d=1:nX
                    for e=1:nXE
                        E(d,e) = (X(d)-XE(e))^2 + (Y(d)-YE(e))^2;
                    end
                end
                IND = indVect(E==min(E(:)));
                [dc ec] = ind2sub([nX, nXE], IND);
                SJ = sqrt( (yJ-Y(dc))^2 + (xJ-X(dc))^2 );
                SR(i,j) = D/(SI+SJ);
            else
            end
        end
    end
    if exist('spaceRatio','var')
        spaceRatio = [spaceRatio; SR(SR>0)];
    else
        spaceRatio = SR(SR>0);
    end
    waitbar(k/nExperiments, g)
end
close(g)





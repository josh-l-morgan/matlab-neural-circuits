%% LOAD FILE
clear spikes
clc
uiopen

%% INPUTDLG
prompt = {'TYPE1','TYPE2'};
title = 'cell selection';
answer = inputdlg(prompt,title);

type1 = str2num(answer{1}); %#ok<ST2NM>
nType1 = length(type1);

type2 = str2num(answer{2}); %#ok<ST2NM>
nType2 = length(type2);

type1_spikes = spikes(type1);
type2_spikes = spikes(type2);

%% DETERMINE RF SPATIAL PARAMETERS

for i=1:nType1
    [z beta0 xy] = gausspreparation (type1_spikes(i).STA);
    gaussOptions = statset ('MaxIter', 10000);
    beta = nlinfit(xy, z, @gaussian2d, beta0, gaussOptions);
    type1_spikes(i).space.beta = beta;
end

for i=1:nType2
    [z beta0 xy] = gausspreparation (type2_spikes(i).STA);
    gaussOptions = statset ('MaxIter', 10000);
    beta = nlinfit(xy, z, @gaussian2d, beta0, gaussOptions);
    type2_spikes(i).space.beta = beta;
end

SR = zeros(nType1, nType2);
for i=1:nType1
    xI = type1_spikes(i).space.beta(2);
    yI = type1_spikes(i).space.beta(3);
    for j=1:nType2
        xJ = type2_spikes(j).space.beta(2);
        yJ = type2_spikes(j).space.beta(3);
        D = sqrt( (yJ-yI)^2 + (xJ-xI)^2 );
        slope = (yJ-yI)/(xJ-xI);
        if ~isnan(slope)
            X = linspace(0,80,1000);
            yzero = yJ - slope*xJ;
            yeighty = yJ - slope*(xJ-80);
            Y = linspace(yzero, yeighty, 1000);
            beta = type1_spikes(i).space.beta;
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

            beta = type2_spikes(j).space.beta;
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
    spaceRatio = [spaceRatio; SR(:)];
else
    spaceRatio = SR(:);
end






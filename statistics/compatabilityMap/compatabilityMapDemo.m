

m = 10;
s = 10;
n = 20;
tnum = 1;
clear ex ps cis
ax = gca;
cla(ax)


if 0

    for t = 1:tnum

        g = randn(n,1)*s+m;

        shiftH = [-10:.1:30];

        ps = zeros(length(shiftH),1);
        cis = zeros(length(shiftH),2);

        [a b c] = ttest(g);
        ex(t).cis  = c';

        for i = 1:length(shiftH)
            [a b c] = ttest(g-shiftH(i));
            ps(i) = b;
        end

        ex(t).ps = ps;

    end

    cla(ax)
    tcol = [0 0 1 ;0 0 0 ; 1 0 0 ;0 0 1 ;0 0 0 ;1 0 0 ;0 0 0 ] *0;
    tMarker = {'-','--',':','-.','-','--',':','-.'};
    tFilled = [0 0 0; 1 1 1; 1 1 1; 1 1 1]
    ax.NextPlot = 'add';
    for t = 1:tnum
        plot(ax,shiftH,ex(t).ps,'color',tcol(t,:),'LineStyle',tMarker{t})
        plot(ax,ex(t).cis, [.05 .05],'color',tcol(t,:),'LineStyle',tMarker{t})
        scatter(ax,ex(t).cis,[.05 .05],40,[0 0 0],'MarkerFaceColor',tFilled(t,:))
    end
    plot(ax,[m m],[0 1],'r')

    xlim(ax,[0 20])
    ylim(ax,[0 1])

else %Model 2

    pop = rand(1000000,1).^3;
    pop = pop-mean(pop);
    popE = sort(abs(pop),"ascend");
    popS = popE(round(length(pop)*.67));
    pop = pop * 1/popS;
    sortPop = sort(pop);
    lowBound = sortPop(floor(length(pop)*.025));
    upBound = sortPop(ceil(length(pop)*.975));

    

    for t = 1:tnum

        g = randsample(pop,n)+m;
        sampMean = mean(g);

        shiftH = [-10:.1:30];

        ps = zeros(length(shiftH),1);
       % ex(t).cis  = sampMean - [lowBound upBound];
        
        [a b c] = ttest(g);
        ex(t).cis  = c';

        for i = 1:length(shiftH)
            realDif = shiftH(i)-sampMean;
            if realDif>0
            ps(i) = mean(popDif >= realDif);
            else
            ps(i) = mean(popDif <= realDif);
            end
        end

        ex(t).ps = ps;

    end


    cla(ax)
    tcol = [0 0 1 ;0 0 0 ; 1 0 0 ;0 0 1 ;0 0 0 ;1 0 0 ;0 0 0 ] *0;
    tMarker = {'-','--',':','-.','-','--',':','-.'};
    tFilled = [0 0 0; 1 1 1; 1 1 1; 1 1 1]
    ax.NextPlot = 'add';
    for t = 1:tnum
        plot(ax,shiftH,ex(t).ps,'color',tcol(t,:),'LineStyle',tMarker{t})
        plot(ax,ex(t).cis, [.05 .05],'color',tcol(t,:),'LineStyle',tMarker{t})
        scatter(ax,ex(t).cis,[.05 .05],40,[0 0 0],'MarkerFaceColor',tFilled(t,:))
    end
    plot(ax,[m m],[0 1],'r')

    xlim(ax,[0 20])
    ylim(ax,[0 1])

end






rand(20,1).^4
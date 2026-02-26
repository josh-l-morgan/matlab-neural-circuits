
clear all
frac = 0.2%[0:.1:.6]
reps = 1000000;
ops = 3;
checkDifs = 0:1:ops;

for f = 1:length(frac)
        hits = rand(reps,ops)<=frac(f);
        sumHits = sum(hits,2);
        histHits = hist(sumHits,checkDifs);
        hitProp = histHits/reps;
        hist(sumHits)
        pause(.1)
        %{
        difs = abs(sumHits-frac(f));
        histDifs = hist(difs,checkDifs)/reps;
        hist(sumHits,checkDifs)
        
        for h = 1:length(histDifs)
            cumHistDifs(h) = sum(histDifs(h:end));            
        end
        %%
        plot(1./cumHistDifs)
        hold on
        plot((1/frac(f)).^(checkDifs.^2/ops),'r')
        hold off
        xlim([0 6]);
        ylim([0 50]);
        pause(.1)
        %}
        
end
    
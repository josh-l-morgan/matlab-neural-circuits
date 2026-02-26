
%Copy all or part of a vsvi database from one folder to another


SPN = 'Y:\Active\morganLab\DATA\pratyush\pratyush_mipmaps_jxQ\';
TPN = 'Y:\Active\sean\JxQ_forPMC\jxQ_mipMaps\';

Xrange = [3930 54977];
Yrange = [9896 27952];
Zrange = [21 250];
%Zrange = [0 20];

Mrange = [0 9];


%% Get image info
done = 0;
for m =  Mrange(1):Mrange(2)
    for z = Zrange(1):Zrange(2)
        fold = sprintf('%s%d\\%d\\',SPN,z,m);
        nams = dir([fold '*.png']);
        for n = 1:length(nams)
            Iinfo = imfinfo([fold nams(n).name]);
            Isize = Iinfo.Width;
            done = 1;
            break
        end
        if done,break,end
    end
    if done, break, end
end


%% 






for z = Zrange(1):Zrange(2)

    for m = Mrange(1):Mrange(2)
        Xds = Xrange/(2^m);
        Yds = Yrange/(2^m);
        cols = [floor(Xds(1)/Isize) ceil(Xds(2)/Isize)];
        rows = [floor(Yds(1)/Isize) ceil(Yds(2)/Isize)];
        nams = {};
        num = 0;
        for c = cols(1):cols(2)
            for r = rows(1):rows(2)
                num = num+1;
                nams{num,1} = sprintf('%d_%d.png',c,r);
            end
        end

        fold = sprintf('%d\\%d\\',z,m);
        if ~exist([TPN fold],'dir'),mkdir([TPN fold]); end
        for n = 1:length(nams)
            fromFile = sprintf('%s%s%s',SPN,fold,nams{n});
            if exist(fromFile,'file')

                toFile = sprintf('%s%s%s',TPN,fold,nams{n});
                disp(sprintf('Copying %s',fromFile))

                try
                    copyfile(fromFile,toFile);
                catch err
                    disp(sprintf('Failed to copy %s',fromFile))

                end
            end
        end
    end
end















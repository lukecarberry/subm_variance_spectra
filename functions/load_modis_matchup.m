function [str] = load_modis_matchup(dt2)


try
    files2 = dir(string(dt2) + '*.mat');
    for s = 1:length(files2)
    S = whos("-file",files2(s).name);
    siz(s) = prod(S(2).size);
    end
    [~,I] = max(siz);
    str = files2(I).name;
    %load(files2(I).name,'chlr','latr','lonr');
catch
    try
        files3 = dir(string(dt2-1) + '*.mat');
        for s = 1:length(files3)
        S = whos("-file",files3(s).name);
        siz(s) = prod(S(2).size);
        end
        [~,I] = max(siz);
        str = files3(I).name;
        %load(files2(I).name,'chlr','latr','lonr');
    catch
        try
            files4 = dir(string(dt2+1) + '*.mat');
            for s = 1:length(files4)
            S = whos("-file",files4(s).name);
            siz(s) = prod(S(2).size);
            end
            [~,I] = max(siz);
            str = files4(I).name;
            %load(files2(I).name,'chlr','latr','lonr');
        catch
            disp(" -> No matchup files found")%continue
        end
    end
end



function [lrs,lrl] = findgoodruns(chlr,dX,max)

% give the chlorophyll region chlr
% get row vectors of;
% lrs = long run starts = indices of first cell left to right of long run
% lrl = long run lengths = length of each long run

nanregion = isnan(chlr);

lrs = NaN(size(nanregion,1),5);
lrl = NaN(size(nanregion,1),5);
minrunlength = round(max/dX);        %1600 * 30m = 50km threshold run length

for i = 1:size(nanregion,1)
I = find(diff(nanregion(i,:))==-1); %going from isnan=1 to 0, start of run
J = find(diff(nanregion(i,:))==1); %going from isnan=0 to 1, end of run
    
    if isempty(I) || isempty(J)
        goodruns = find(nanregion(i,:)==0,1);
        goodrunlength = find(nanregion(i,:)==0,1,'last')-goodruns;
        
        if goodrunlength > minrunlength
            lrs(i,1:length(goodruns)) = goodruns; % long run starts
            lrl(i,1:length(goodruns)) = goodrunlength; % long run lengths
        else
            continue
        end
    else
        [allruns,id] = sort([I,J]);
        runlength = diff(allruns);
        f = find(id==1);
        goodruns = allruns(f:2:end)+1; %this is the index for the start of each good run
        goodrunlength = runlength(f:2:end);
    
        lr = find(goodrunlength>minrunlength); % long runs
    
        if lr > 0
            lrs(i,1:length(goodruns(lr))) = goodruns(lr); % long run starts
            lrl(i,1:length(goodruns(lr))) = goodrunlength(lr); % long run lengths
        else
            continue
        end
    end
end




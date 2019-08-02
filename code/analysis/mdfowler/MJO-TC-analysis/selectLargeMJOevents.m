% This function takes input of a 12x8 array of MJO indices and time (such
% as that stored in 'MJOindex&time_monthPhase.mat'), and with the user's
% choice of an integer month/phase, returns the appropriate number of indices 
% and dates to be used in extractDataForKerry.m. 
%
% For example, if the user wanted a subset of data for Phase 3 in January,
% there are 118 possible choices of dates to select. In this case, the MJO
% indices are sorted and the largest 31 of those would be selected and
% returned to the user. 
%
% On the other hand, if there are less days than are necessary for the
% month (e.g., Phase 3 in September has only 26 days), the appropriate
% number of days are randomly selected for repetition. 
%
% INPUT: MJOindex and MJOtime - from MJOindex&time_monthPhase.mat,
%           or otherwise 12x8 cell arrays that match the format of those. 
%        month - user specified month to get indices/dates for 
%        phase - user specified MJO phase to get indices/dates for
% OUTPUT: getTimes - array of times for which obs should be retrieved
%
% Meg D. Fowler, 2017-06-14


function [getTimes]=selectLargeMJOevents(MJOindex,MJOtime,month,phase)
    
    if nargin==2        %Default month/phase settings
       month = 1;
       phase = 1;
    end

    monDays   = [31,28,31,30,31,30,31,31,30,31,30,31];    %Number of days/month
    needDays  = monDays(month);                           %Number of days in selected month

    selectTime = MJOtime{month,phase};                    %Number of observations in month/phase   
    
    %Determine how to choose indices 
    if numel(selectTime)>needDays       %More observations than needed 
       
        MJOmag = MJOindex{month,phase};             %Magnitude of MJO index in month/phase
        [sortMJO,isort] = sort(MJOmag,'descend');   %Sort from high to low
        sortTime = selectTime(isort);               %Match time order to MJO sort
        
        selection = sortMJO(1:needDays);
        getTimes  = sortTime(1:needDays); 
    elseif numel(selectTime)<needDays   %Less observations than needed
        diff = needDays - numel(selectTime); %Number of times to resample
        r = randi([1 numel(selectTime)],1,diff);
        
        getTimes = [selectTime;selectTime(r)];        
    elseif numel(selectTime)==needDays
        getTimes = selectTime;
        
    end


end




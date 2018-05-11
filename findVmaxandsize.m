function [ Vmaxbyweek, Vmaxweekavg, ninweek, wknum, Vmaxall] = findVmaxandsize(vddata)
% Function tells you the average Max Viability per week (Vmaxbyweek), the
% average over all of the data (Vmaxweekavg), the number of dose-response
% curves in each week (so number of different cohorts in that week) the
% number of actual data points in each week of data (wknum(:,2), generally
% 12 x ninweek but some weeks have less data) and Vmaxall which puts the
% Vmax for each week (Vmax by week) alongside each dose and viability point
% to be used to calculate the fsens parameter for that week

 dose0ind = vddata(:,2) == 0; % finds all points where dose = 0
 Vmaxdata = vddata(dose0ind,:); % extracts only rows where dose = 0

wkidentifier = Vmaxdata(:,1); % identifies wk number of row where dose = 0
weeks = unique(Vmaxdata(:,1)); % identifies the unique values of the week identifier
Vmaxweek = zeros([length(weeks), 1]); % sets up the Vmax in each week (so 12 x1)

for i = 1:length(weeks)
    wklycounter(i,1) = nnz(wkidentifier==i); % iterates through Vmax data and counts
    % the number of times the wk = i which spans the length of the weeks
end

for t = 1:length(weeks)  % iterates through and finds when Vmaxdata corresponds to that week
   wkind = Vmaxdata(:,1) == weeks(t); 
   Vmaxweekadd(t) = sum(Vmaxdata(wkind,3)); % computes sum of the Vmaxes at each week
   Vmaxbyweek(t,1) = Vmaxweekadd(t)./wklycounter(t); % divides by the number of weeks and adds 
   % to Vmaxbyweek vector
end

Vmaxweekavg = sum(Vmaxweekadd)./ sum(wklycounter);
ninweek = wklycounter;

% Now find 
weekidentifierall = vddata(:,1);

wknumber = zeros([length(weeks) 1]);
wkcount = zeros([length(weeks) 1]);

for i = 1:length(weeks)
wkcount(i) = nnz(weekidentifierall==i);
wknumber(i,1) = i;
end

wknum = horzcat(wknumber, wkcount);


for i = 1:length(vddata)
    if vddata(i,1) == 1
    vddata(i,5) = Vmaxbyweek(1,1);
    end
    if vddata(i,1) == 2
    vddata(i,5) = Vmaxbyweek(2,1);
    end
    if vddata(i,1) == 3
    vddata(i,5) = Vmaxbyweek(3,1);
    end
    if vddata(i,1) == 4
    vddata(i,5) = Vmaxbyweek(4,1);
    end
    if vddata(i,1) == 5
    vddata(i,5) = Vmaxbyweek(5,1);
    end
    if vddata(i,1) == 6
    vddata(i,5) = Vmaxbyweek(6,1);
    end
    if vddata(i,1) == 7
    vddata(i,5) = Vmaxbyweek(7,1);
    end
    if vddata(i,1) == 8
    vddata(i,5) = Vmaxbyweek(8,1);
    end
    if vddata(i,1) == 9
    vddata(i,5) = Vmaxbyweek(9,1);
    end
    if vddata(i,1) == 10
    vddata(i,5) = Vmaxbyweek(10,1);
    end
    if vddata(i,1) == 11
    vddata(i,5) = Vmaxbyweek(11,1);
    end
    if vddata(i,1) == 12
    vddata(i,5) = Vmaxbyweek(12,1);
    end
    if vddata(i,1) == 13
    vddata(i,5) = Vmaxbyweek(13,1); 
    end
    if vddata(i,1) == 14
    vddata(i,5) = Vmaxbyweek(14,1);
    end
    if vddata(i,1) == 15
    vddata(i,5) = Vmaxbyweek(15,1);
end
Vmaxall = vddata(:,5);

% Next need to make a vector containing the Vmax by week alongside dose and
% var

  

end



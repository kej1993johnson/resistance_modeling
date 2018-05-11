function [ Vmaxbyweek, Vmaxweekavg, ninweek, wknum, Vmaxall] = findVmaxandsizetest15(vddata)

% This function performs the same task as in findVmaxandsizetest, but
% allows for weeks to be missing in computing Vmaxbyweek ( with Nans for
% missing weeks, average Vmax for all weeks in the testing set (cohort),
% tells the weeks with values (ninweek, should be 1 or 0), then computes
% the number of data points in each week (usually should be 12 but some
% missing data), and then puts Vmaxs next to all data points (Vmaxall)

 dose0ind = vddata(:,2) == 0;
 Vmaxdata = vddata(dose0ind,:);

wkidentifier = Vmaxdata(:,1);
weeks = unique(Vmaxdata(:,1));
Vmaxweek = zeros([15, 1]);
ninweek = zeros([15,1]);
time = 1:15;


for i = 1:length(time)
    wklycounter(i,1) = nnz(wkidentifier==i);
end

ninweek = wklycounter;

for t = 1:length(time)
   wkind = Vmaxdata(:,1) == time(t);
   Vmaxweekadd(t) = sum(Vmaxdata(wkind,3));
end
Vmaxweekavg = sum(Vmaxweekadd)./sum(ninweek);
Vmaxbyweek = Vmaxweekadd'./ninweek;
Vmaxbyweekdbl(:,1) = ninweek;
Vmaxbyweekdbl(:,2) = Vmaxbyweek;

vddata(:,5) = zeros([length(vddata) 1]);

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
end

Vmaxall = vddata(:,5);

% Now find 
weekidentifierall = vddata(:,1);

wknumber = zeros([length(time) 1]);
wkcount = zeros([length(time) 1]);

for i = 1:length(time)
wkcount(i) = nnz(weekidentifierall==i);
wknumber(i,1) = i;
end

wknum = horzcat(wknumber, wkcount);



  

end
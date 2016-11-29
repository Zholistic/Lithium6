xtofit832 = [80 36 20];
ytofit832 = [43.9 44.7 46.9];
ytofiterror832 = [0.1 0.5 0.5];
ytofit880 = [43];
ytofiterror880 = [1];
xtofit880 = [80];
ytofit860 = [43.82 45.62 50.5];
ytofiterror860 = [0.13 0 0.15];
xtofit860 = [80 20 10];



figure(555);
errorbar(xtofit832,ytofit832,ytofiterror832,'o');
hold on;
%errorbar(xtofit880,ytofit880,ytofiterror880,'o');
errorbar(xtofit860,ytofit860,ytofiterror860,'x');
line([10 90],[24.1*2 24.1*2]);
line([10 90],[24.1*((3)^(1/2)) 24.1*((3)^(1/2))]);
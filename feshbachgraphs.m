close all;

a0 = 5.29e-11;


if(1)
abg = -1405*a0./(1000*a0);
deltaB = 30;
B0 = 83.22;
alpha = 0.004;    
end

abg13 = -1727*a0./(1000*a0);
deltaB13 = 12.23;
B013 = 69.043;
alpha13 = 0.002;

xs = 1:0.1:2000;
a = abg .* (1+ (deltaB./(xs - B0))).*(1+alpha.*(xs - B0));
%a = abg .* (1+ (deltaB./(xs - B0)));
xs13 = 1:0.1:2000;
a13 = abg .* (1+ (deltaB13./(xs13 - B013))).*(1+alpha13.*(xs13 - B013));


plot(xs*10,a,'-');
grid on;
axis([550 1000 -20 20]);
hold on;
plot(xs13*10,a13,'-');


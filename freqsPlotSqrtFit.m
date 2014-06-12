freqs = [];

freqs(1,1) = 5;
freqs(1,2) = 4;
freqs(1,3) = 1.75;
freqs(1,4) = 0;
freqs(2,1) = 5711;
freqs(2,2) = 4866;
freqs(2,3) = 2927;
freqs(2,4) = 0;
freqs(3,1) = 375;
freqs(3,2) = 254;
freqs(3,3) = 21;
freqs(3,4) = 0;

freqsT(1,1) = 5;
freqsT(1,2) = 4;
freqsT(1,3) = 3;
freqsT(1,4) = 2;
freqsT(1,5) = 1;
freqsT(1,6) = 0;
freqsT(2,1) = 6030;
freqsT(2,2) = 5394;
freqsT(2,3) = 4671;
freqsT(2,4) = 3814;
freqsT(2,5) = 2696;
freqsT(2,6) = 0;

fg = @(p,x)(p(1)*x.^(1/2));
p0 = [1000];
lb = [100];
ub = [7000];


errorbarxy(freqs(1,:),freqs(2,:),zeros(4),freqs(3,:),zeros(4),freqs(3,:),'.');

[coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,freqsT(1,:),freqsT(2,:),lb,ub);
hold on;
largexs = 0.001:0.001:8;
plot(largexs,fg(coefs,largexs),'k','LineWidth',0.5);

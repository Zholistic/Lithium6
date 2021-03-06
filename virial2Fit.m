function [ outputFunc, outputGof, outputFoutput ] = virial2Fit( profileToFit, xs, vcoefs, omegaz)
%Fit Virial2 Function with T and m0 as parameters to fit too; takes virial
%co-efficients as constants.

b2 = vcoefs(1);
b3 = vcoefs(2);

massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
kB = 1.38e-23; %m^2 kg s^-2 K^-1
%omegazR = omegaz * 2 * pi;


%CONVERSION to units that won't break because e^1000 is too much for matlab
%to handle:

%profileToFit = profileToFit.*10e-11;
%xs = xs.*1;

    %Make sure profile is row vector:
    if(iscolumn(profileToFit))
        profileToFit = profileToFit';
    end
    
    if(iscolumn(xs))
        xs = xs';
    end

    %Function to fit with:
    %fg = @(p,x)(p(1).*exp(-p(2).*x).*sin(p(3).*x+p(4))-p(5)*x+p(6));
    fg =  @(p,x)((2/((2*pi*hbar^2)/(massL6*kB*p(1)))) * ...
        (log(1 + exp((1/(kB*p(1))) * (p(2) - x))) + ...
        (2*b2*exp(2*(1/(kB*p(1)))*(p(2) - x))) + ...
        (3*b3*exp(3*(1/(kB*p(1)))*(p(2) - x)))));
    
    fg =  @(p,x)((1/(massL6*omegaz/hbar))*(2/((2*pi*hbar^2)/(massL6*kB*(p(1)*1e-9)))) * ...
        (log(1 + exp((1/((hbar*omegaz)/2))*(1/(kB*(p(1)*1e-9))) * (p(2) - x))) + ...
        (2*b2*exp((1/((hbar*omegaz)/2))*2*(1/(kB*(p(1)*1e-9)))*(p(2) - x))) + ...
        (3*b3*exp((1/((hbar*omegaz)/2))*3*(1/(kB*(p(1)*1e-9)))*(p(2) - x)))));
 
    fg =  @(p,x)((1/(massL6*omegaz/hbar))*(2/((2*pi*hbar^2)/(massL6*kB*(p(1))))) * ...
        (log(1 + exp((1/((hbar*omegaz)/2))*(1/(kB*p(1))) * (p(2) - x))) + ...
        (2*b2*exp((1/((hbar*omegaz)/2))*2*(1/(kB*p(1)))*(p(2) - x))) + ...
        (3*b3*exp((1/((hbar*omegaz)/2))*3*(1/(kB*p(1)))*(p(2) - x)))));
    
    fg2 =  @(p)((1/(massL6*omegaz/hbar))*(2/((2*pi*hbar^2)/(massL6*kB*(p(1))))) * ...
        (log(1 + exp((1/((hbar*omegaz)/2))*(1/(kB*p(1))) * (p(2) - xs))) + ...
        (2*b2*exp((1/((hbar*omegaz)/2))*2*(1/(kB*p(1)))*(p(2) - xs))) + ...
        (3*b3*exp((1/((hbar*omegaz)/2))*3*(1/(kB*p(1)))*(p(2) - xs))))) - profileToFit;
    
    
    
            if(0)
                   (2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*a))) * ...
                       (log(1 + exp((1/((1.38e-23)*a)) * (b - x))) + ...
                       (2*0.5179*exp(2*(1/((1.38e-23)*a))*(b - x))) + ...
                       (3*(-0.4215)*exp(3*(1/((1.38e-23)*a))*(b - x))))
                   
                   (2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*a))) * (log(1 + exp((1/((1.38e-23)*a)) * (b - x))) + (2*0.5179*exp(2*(1/((1.38e-23)*a))*(b - x))) + (3*(-0.4215)*exp(3*(1/((1.38e-23)*a))*(b - x))))
                   
                   (2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*40e-9))) * (log(1 + exp((1/((1.38e-23)*40e-9)) * (1.4e-30 - x))) + (2*0.5179*exp(2*(1/((1.38e-23)*40e-9))*(1.4e-30 - x))) + (3*(-0.4215)*exp(3*(1/((1.38e-23)*40e-9))*(1.4e-30 - x))))
                   (2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*40e-9))) * (log(1 + exp((1/((1.38e-23)*40e-9)) * (1.4e-30 - xs))) + (2*0.5179*exp(2*(1/((1.38e-23)*40e-9))*(1.4e-30 - xs))) + (3*(-0.4215)*exp(3*(1/((1.38e-23)*40e-9))*(1.4e-30 - xs))))
                   ((1/(massL6*omegaz / hbar))*2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*40e-9))) * (log(1 + exp((1/((hbar*omegaz)/2))*(1/((1.38e-23)*40e-9)) * (1.4e-30 - xs))) + (2*0.5179*exp((1/((hbar*omegaz)/2))*2*(1/((1.38e-23)*40e-9))*(1.4e-30 - xs))) + (3*(-0.4215)*exp((1/((hbar*omegaz)/2))*3*(1/((1.38e-23)*40e-9))*(1.4e-30 - xs))))
                   plot(xs2,(2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*40e-9))) * (log(1 + exp((1/((1.38e-23)*40e-9)) * (1.4e-30 - xs2))) + (2*0.5179*exp(2*(1/((1.38e-23)*40e-9))*(1.4e-30 - xs2))) + (3*(-0.4215)*exp(3*(1/((1.38e-23)*40e-9))*(1.4e-30 - xs2)))) )
                   plot(xs(16:end),(2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*40e-9))) * (log(1 + exp((1/((1.38e-23)*40e-9)) * (1.4e-30 - xs(16:end)))) + (2*0.5179*exp(2*(1/((1.38e-23)*40e-9))*(1.4e-30 - xs(16:end)))) + (3*(-0.4215)*exp(3*(1/((1.38e-23)*40e-9))*(1.4e-30 - xs(16:end))))) )
                   ((1/(9.988e-27*3.6373e+04 / 1.05457e-34))*(2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*a)))) * (log(1 + exp((((1.05457e-34*3.6373e+04)/2))*(1/((1.38e-23)*a)) * (b - x))) + (2*1.4259*exp((((1.05457e-34*3.6373e+04)/2))*2*(1/((1.38e-23)*a))*(b - x))) + (3*(-1.0561)*exp(((1.05457e-34*3.6373e+04)/2)*3*(1/((1.38e-23)*a))*(b - x))))
                   
                   
                   close all;
                   xs2 = 1.9e-30:0.1e-31:2.7e-30;
                   xs2 = xs2.*(1/((hbar*omegaz)/2)); %conv to HO units
                   tempGuess2 = 28e-9;
                   mu0Guess2 = 1.75e-30;
                   mu0Guess2 = mu0Guess2.*(1/((hbar*omegaz)/2));
                   %plot(xs2,(2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*tempGuess2))) * (log(1 + exp((1/((1.38e-23)*tempGuess2)) * (mu0Guess2 - xs2))) + (2*b2*exp(2*(1/((1.38e-23)*tempGuess2))*(mu0Guess2 - xs2))) + (3*b3*exp(3*(1/((1.38e-23)*tempGuess2))*(mu0Guess2 - xs2)))),'r')
                   plot(xs2,((1/(massL6*omegaz / hbar))*(2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*tempGuess2)))) * (log(1 + exp((((hbar*omegaz)/2))*(1/((1.38e-23)*tempGuess2)) * (mu0Guess2 - xs2))) + (2*b2*exp((((hbar*omegaz)/2))*2*(1/((1.38e-23)*tempGuess2))*(mu0Guess2 - xs2))) + (3*b3*exp(((hbar*omegaz)/2)*3*(1/((1.38e-23)*tempGuess2))*(mu0Guess2 - xs2)))),'r')
                   %plot(xs2,((1/(massL6*omegaz / hbar))*2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*tempGuess2))) * (log(1 + exp((1/((1.38e-23)*tempGuess2)) * (mu0Guess2 - xs2))) + (2*b2*exp(2*(1/((1.38e-23)*tempGuess2))*(mu0Guess2 - xs2))) + (3*b3*exp(3*(1/((1.38e-23)*tempGuess2))*(mu0Guess2 - xs2)))),'r')
                   hold on; plot(xs,profileToFit,'.');
                    
                    
                    
                    %plot(xs2,((1/(massL6*omegaz / hbar))*2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*tempGuess2))) * (log(1 + exp((1/((hbar*omegaz)/2))*(1/((1.38e-23)*tempGuess2)) * (mu0Guess2 - xs2))) + (2*b2*exp((1/((hbar*omegaz)/2))*2*(1/((1.38e-23)*tempGuess2))*(mu0Guess2 - xs2))) + ((1/((hbar*omegaz)/2))*3*b3*exp(3*(1/((1.38e-23)*tempGuess2))*(mu0Guess2 - xs2)))),'r')
              end
    
    %Initial guesses:
    %[T mu0]
    p0 = [30e-9 0.3e-30];
    lb = [10e-9 0.01e-30];
    ub = [110e-9 2.5e-30]; 
    
    %CONV to H.O. units
    p0(2) = p0(2)*(1/((hbar*omegaz)/2));
    lb(2) = lb(2)*(1/((hbar*omegaz)/2));
    ub(2) = ub(2)*(1/((hbar*omegaz)/2));
    
    %p0(1) = p0(1)*1e9;
    %lb(1) = lb(1)*1e9;
    %ub(1) = ub(1)*1e9;
    
    x = xs;
    y = profileToFit;
    
    %virialEqn = '((1/(9.988e-27*3.6373e+04 / 1.05457e-34))*(2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*a)))) * (log(1 + exp((((1.05457e-34*3.6373e+04)/2))*(1/((1.38e-23)*a)) * (b - x))) + (2*1.4259*exp((((1.05457e-34*3.6373e+04)/2))*2*(1/((1.38e-23)*a))*(b - x))) + (3*(-1.0561)*exp(((1.05457e-34*3.6373e+04)/2)*3*(1/((1.38e-23)*a))*(b - x))))';
    virialEqn = ['((1/(9.988e-27*' num2str(omegaz) '/ 1.05457e-34))*(2/((2*pi*(1.05457e-34)^2)/((9.988e-27)*(1.38e-23)*a)))) * (log(1 + exp((((1.05457e-34*' num2str(omegaz) ')/2))*(1/((1.38e-23)*a)) * (b - x))) + (2*' num2str(vcoefs(1)) '*exp((((1.05457e-34*' num2str(omegaz) ')/2))*2*(1/((1.38e-23)*a))*(b - x))) + (3*(' num2str(vcoefs(2)) ')*exp(((1.05457e-34*' num2str(omegaz) ')/2)*3*(1/((1.38e-23)*a))*(b - x))))'];
    [f1,gof,foutput] = fit(x',y',virialEqn,'Start', p0,'Lower',lb,'Upper',ub,'MaxFunEvals',100000,'MaxIter',100000,'TolFun',1e-11,'TolX',1e-11);
    
    
    %Fitting:
    %xs = [];
    %curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',100000,'Display','off','TolFun',1e-11,'TolX',1e-11);
    %xs = 1:length(profileToFit);
    %[coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs,profileToFit,lb,ub,curvefitoptions);
    %[coefs,resnorm,r,exitflag,output,lambda,J] = lsqnonlin(fg,p0,lb,ub);
    
    outputFunc = f1;
    outputGof = gof;
    outputFoutput = foutput;
    %outputCoefError = nlparci(coefs,r,'jacobian',J);

end
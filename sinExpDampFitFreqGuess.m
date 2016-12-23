function [ outputCoefs, outputCoefError ] = sinExpDampFitFreqGuess( profileToFit, xs, freqguess )
%Takes a 1D array half-profile with gaussian shape and fits using lsqcurvefit.
%Returns the co-efficients of the fit. 

    %Make sure profile is row vector:
    if(iscolumn(profileToFit))
        profileToFit = profileToFit';
    end
    
    if(iscolumn(xs))
        xs = xs';
    end

    %Function to fit with:
    fg = @(p,x)(p(1).*exp(-p(2).*x).*sin(p(3).*x+p(4))+p(5));
    
    %Initial guesses:
    peak = max(profileToFit); %Sensitive to highest value in array (problem if noise)
    centerIndex = find(profileToFit == peak);
    magicWidth = ceil((1/3)*length(profileToFit));
    
    %For ~5kHz trap freq fits:
    if(0)
    p0 = [3.7 0.006 30000 -.18 mean(profileToFit)];
    lb = [0.1 0.0001 25000 -12 mean(profileToFit)-4];
    ub = [7 0.1 35000 12 mean(profileToFit)+4]; %0.5 is ~ 80Hz
    end
    
    %For PCA/50Hz fits    
    p0 = [0.28 0.05 freqguess 0 mean(profileToFit)];
    lb = [-3 -0.01 -1 -4 mean(profileToFit)-10];
    ub = [3 0.2 1 4 mean(profileToFit)+10]; %0.5 is ~ 80Hz 0.32 ~ 50Hz
    
    %Fitting:
    %xs = [];
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    %xs = 1:length(profileToFit);
    [coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs,profileToFit,lb,ub,curvefitoptions);
    
    outputCoefs = coefs;
    outputCoefError = nlparci(coefs,r,'jacobian',J);

end
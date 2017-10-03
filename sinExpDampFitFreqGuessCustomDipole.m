function [ outputCoefs, outputCoefError ] = sinExpDampFitFreqGuessCustomDipole( profileToFit, xs, freqguess, ampguess )
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
    
    %For PCA/50Hz fits    
    p0 = [ampguess 0 freqguess 0 mean(profileToFit)];
    lb = [(ampguess-(0.2*ampguess)) -0.001 (freqguess-0.04) -4 mean(profileToFit)-5];
    ub = [(ampguess+(0.2*ampguess)) 0.001 (freqguess+0.04) 1 mean(profileToFit)+5]; %0.5 is ~ 80Hz 0.32 ~ 50Hz
    
    %Fitting:
    %xs = [];
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    %xs = 1:length(profileToFit);
    [coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs,profileToFit,lb,ub,curvefitoptions);
    
    outputCoefs = coefs;
    outputCoefError = nlparci(coefs,r,'jacobian',J);

end
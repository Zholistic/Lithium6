function [ outputCoefs, outputCoefError ] = gausFit1DLockZero( profileToFit )
%Takes a 1D array profile with gaussian shape and fits using lsqcurvefit.
%Returns the co-efficients of the fit.

    %Make sure profile is row vector:
    if(iscolumn(profileToFit))
        profileToFit = profileToFit';
    end

    %Function to fit with:
    fg1d = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2))); 
    
    %Initial guesses:
    magicWidth = ceil((1/8)*length(profileToFit));
    peak = max(profileToFit(magicWidth*2:end-(magicWidth*2))); %Sensitive to highest value in array (problem if noise)
    centerIndex = find(profileToFit == peak,1);
    
    
    p0 = [peak centerIndex magicWidth];
    lb = [peak/2 centerIndex-magicWidth magicWidth/50];
    ub = [peak*2 centerIndex+magicWidth magicWidth*2];
    
    %Fitting:
    xs = [];
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    xs = 1:length(profileToFit);
    [coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg1d,p0,xs,profileToFit,lb,ub,curvefitoptions);
    
    if(length(coefs) ~= 3)
        display('Error on gaus fit, coefs length incorrect');
        coefs = [0 0 0];
    end
    
    outputCoefs = coefs;
    outputCoefError = nlparci(coefs,r,'jacobian',J);

end


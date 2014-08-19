function [ outputCoefs ] = gausFit1D( profileToFit )
%Takes a 1D array profile with gaussian shape and fits using lsqcurvefit.
%Returns the co-efficients of the fit.

    %Make sure profile is row vector:
    if(iscolumn(profileToFit))
        profileToFit = profileToFit';
    end

    %Function to fit with:
    fg = @(p,x)(p(1).*exp((-1).*((x-p(2)).^2) ./ (2.*p(3).^2)) + p(4)); 
    
    %Initial guesses:
    peak = max(profileToFit); %Sensitive to highest value in array (problem if noise)
    centerIndex = find(profileToFit == peak);
    magicWidth = ceil((1/8)*length(profileToFit));
    
    p0 = [peak centerIndex magicWidth 1];
    lb = [peak/2 centerIndex-magicWidth magicWidth/50 -magicWidth*2];
    ub = [peak*2 centerIndex+magicWidth magicWidth*2 magicWidth*2];
    
    %Fitting:
    xs = [];
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    xs = 1:length(profileToFit);
    [coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs,profileToFit,lb,ub,curvefitoptions);
    
    if(length(coefs) ~= 4)
        display('Error on gaus fit, coefs length incorrect');
        coefs = [0 0 0 0];
    end
    
    outputCoefs = coefs;

end


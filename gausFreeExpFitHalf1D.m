function [ outputCoefs ] = gausFreeExpFitHalf1D( profileToFit, xs )
%Takes a 1D array half-profile with gaussian shape and fits using lsqcurvefit.
%Returns the co-efficients of the fit. 

    %Make sure profile is row vector:
    if(iscolumn(profileToFit))
        profileToFit = profileToFit';
    end

    %Function to fit with:
    fg = @(p,x)(p(1).*exp((-1).*((x).^(p(2))) ./ (2.*p(3).^(2.5)))); 
    
    %Initial guesses:
    peak = max(profileToFit); %Sensitive to highest value in array (problem if noise)
    centerIndex = find(profileToFit == peak);
    magicWidth = ceil((1/3)*length(profileToFit));
    
    p0 = [peak 2.5 magicWidth];
    lb = [peak/2 2 magicWidth/10];
    ub = [peak*2 4 magicWidth*4];
    
    %Fitting:
    %xs = [];
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    %xs = 1:length(profileToFit);
    [coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs,profileToFit,lb,ub,curvefitoptions);
    
    outputCoefs = coefs;

end
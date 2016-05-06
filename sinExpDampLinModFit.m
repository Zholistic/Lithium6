function [ outputCoefs, outputCoefError ] = sinExpDampLinModFit( profileToFit, xs )
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
    fg = @(p,x)(p(1).*exp(-p(2).*x).*sin(p(3).*x+p(4))-p(5)*x+p(6));
    
    %Initial guesses:
    peak = max(profileToFit); %Sensitive to highest value in array (problem if noise)
    centerIndex = find(profileToFit == peak);
    magicWidth = ceil((1/3)*length(profileToFit));
    
    p0 = [3.7 216 35000 2 450 mean(profileToFit)];
    lb = [0.1 150 25000 0 420 mean(profileToFit)-4];
    ub = [7 230 40000 10 480 mean(profileToFit)+4]; %0.5 is ~ 80Hz
    
    %Fitting:
    %xs = [];
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    %xs = 1:length(profileToFit);
    [coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs,profileToFit,lb,ub,curvefitoptions);
    
    outputCoefs = coefs;
    outputCoefError = nlparci(coefs,r,'jacobian',J);

end
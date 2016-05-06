function [ outputCoefs ] = sinFit( profileToFit, xs )
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
    fg = @(p,x)(p(1).*sin(p(2).*x+p(3))+p(4));
    
    %Initial guesses:
    peak = max(profileToFit); %Sensitive to highest value in array (problem if noise)
    centerIndex = find(profileToFit == peak);
    magicWidth = ceil((1/3)*length(profileToFit));
    
    p0 = [3.7 0.24 -.18 mean(profileToFit)];
    lb = [0.1 0.01 -12 mean(profileToFit)-4];
    ub = [4 1 12 mean(profileToFit)+4];
    
    %Fitting:
    %xs = [];
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    %xs = 1:length(profileToFit);
    [coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs,profileToFit,lb,ub,curvefitoptions);
    
    outputCoefs = coefs;

end
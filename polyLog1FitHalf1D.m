function [ outputCoefs ] = polyLog1FitHalf1D( profileToFit, xs, camera )
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
    %fg = @(p,x)(p(1).*exp((-1).*((x).^(2.5)) ./ (2.*p(2).^(2.5))));
    fg = @(p,x)(p(1).*log(1+exp((p(2)+(-1).*x.^2)./(p(3).^2))));
    
    %Initial guesses:
    peak = max(profileToFit); %Sensitive to highest value in array (problem if noise)
    centerIndex = find(profileToFit == peak);
    magicWidth = ceil((1/3)*length(profileToFit));
    
    
    
    p0 = [0.3 2856 12];
    lb = [0 500 1];
    ub = [1.8 5000 35];
    if(strcmp(camera,'sidecam'))
        p0 = [100 6000 55];
        lb = [0 500 1];
        ub = [500 20000 100];
    end

    
    %Fitting:
    %xs = [];
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    %xs = 1:length(profileToFit);
    [coefs,resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fg,p0,xs,profileToFit,lb,ub,curvefitoptions);
   
    if(coefs(1) == lb(1) || coefs(2) == lb(2) || coefs(3) == lb(3))
        disp('Polylog fit hit lower bound!');
    end
    if(coefs(1) == ub(1) || coefs(2) == ub(2) || coefs(3) == ub(3))
        disp('Polylog fit hit upper bound!');
    end
    
    outputCoefs = coefs;

end
function [ outputCoefs ] = gausFit2D( imageToFit )
%Takes a 2D Image with gaussian shape and fits using lsqcurvefit.
%Returns the co-efficients of the fit.

    %Find center of the hump using 1-d fits:
    gcoefsX = gausFit1D(mean(imageToFit(:,:),1)); %mean averages over y
    gcoefsY = gausFit1D(mean(imageToFit(:,:),2)); %mean averages over x
    center = [ceil(gcoefsX(2)), ceil(gcoefsY(2))]; %center = [x y]

    MdataSizeX = length(imageToFit(1,:)); 
    MdataSizeY = length(imageToFit(:,1));
    %[X,Y] = meshgrid(-center(1):(MdataSizeX-center(1)-1),-(center(2)):(MdataSizeY-center(2)-1));
    [X,Y] = meshgrid(1:(MdataSizeX),1:(MdataSizeY));
    xdata = zeros(MdataSizeY,MdataSizeX,2);
    xdata(:,:,1) = X;
    xdata(:,:,2) = Y;
    
    %Initial guesses:
    peak = max(max(imageToFit)); %Sensitive to highest value in array (problem if noise)
    %TODO centerIndex = find(imageToFit == peak);
    %TODO magicWidth = ceil((1/8)*length(imageToFit));
    
    %[Amp,xo,wx,yo,wy,fi]
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    x0 = [peak,110,50,80,50,0]; %Inital guess parameters
    %lb = [peak/2,-40,1,-40,1,-pi/4];
    %ub = [peak*2,40,25,40,25,pi/4];
    lb = [peak/2,50,1,50,1,-pi/2];
    ub = [peak*2,150,50,150,50,pi/2];
    %lb = [0,-MdataSizeX/2,0,-MdataSizeY/2,0,-pi/4];
    %ub = [realmax('double'),MdataSizeX/2,(MdataSizeX/2)^2,MdataSizeY/2,(MdataSizeY/2)^2,pi/4];
    [coefs,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunctionRot,x0,xdata,imageToFit,lb,ub,curvefitoptions);
    
    
    %Plotting code:
    if(0)
    figure(1)
    %C = del2(imageToFit);
    mesh(X,Y,imageToFit) %plot data
    hold on;
    surface(X,Y,D2GaussFunctionRot(coefs,xdata)) %plot fit
    %axis([-MdataSizeX/2-0.5 MdataSizeX/2+0.5 -MdataSizeY/2-0.5 MdataSizeY/2+0.5])
    %axis([])
    alpha(0.8);
    hold off;
    end
   
    
    outputCoefs = coefs;

end


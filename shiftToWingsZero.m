function [ shiftedImage ] = shiftToWingsZero( imageToShift )
%shiftToWingsZero: fit a 1D gaussian to the profile and then move the wings
%to the zero point. 

gcoefsX = gausFit1D(mean(imageToShift,1)); %mean averages over y
gcoefsY = gausFit1D(mean(imageToShift,2)); %mean averages over x

xshift = gcoefsX(4);
yshift = gcoefsY(4);
avgshift = (xshift+yshift)/2;

if(avgshift > 0)
shiftedImage = imageToShift - abs(avgshift);
else 
shiftedImage = imageToShift + abs(avgshift);
end


end


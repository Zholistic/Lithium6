function [profileOutput] = radAverageBigSquare(inputImage)
%A function which takes a 2D image and radially averages. 
%Returns the 1D profile of that image from center to edge. 

%First find the center of the cloud:
gcoefsX = gausFit1D(mean(inputImage,1)); %mean averages over y
gcoefsY = gausFit1D(mean(inputImage,2)); %mean averages over x

cCenter = [round(gcoefsX(2)), round(gcoefsY(2))]; %[x y]


for i=1:length(inputImage(:,1)) %y
    for j=1:length(inputImage(1,:)) %x       
        rVectorsL(i,j) = sqrt((i-cCenter(2))*(i-cCenter(2)) + (j-cCenter(1))*(j-cCenter(1)));      
    end
end


%Binning:
reImageAverage = reshape(inputImage, 1, length(inputImage(:,1))*length(inputImage(1,:)));
reRVectors = reshape(rVectorsL, 1, length(rVectorsL(:,1))*length(rVectorsL(1,:)));

nbins = 100; binMean = [];
binEdges = linspace(min(min(rVectorsL)),max(max(rVectorsL)),nbins+1);
binLength = (max(max(rVectorsL))-min(min(rVectorsL)))/(nbins+1);

[h,whichBin] = histc(reRVectors, binEdges);

for i = 1:nbins
    flagBinMembers = (whichBin == i);
    binMembers     = reImageAverage(flagBinMembers);
    binMean(i)     = mean(binMembers);
end

radialProfile = [];

for i = 1:length(binMean)
    radialProfile(1,i) = binMean(i);
    radialProfile(2,i) = i*binLength; %pixel vector
end

profileOutput = radialProfile;


end
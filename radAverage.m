function [profileOutput] = radAverage(inputImage)
%A function which takes a 2D image and radially averages. 
%Returns the 1D profile of that image from center to edge. 

%First find the center of the cloud:
gcoefsX = gausFit1D(mean(inputImage,1)); %mean averages over y
gcoefsY = gausFit1D(mean(inputImage,2)); %mean averages over x

cCenter = [round(gcoefsX(2)), round(gcoefsY(2))]; %[x y]

%Take square region about center.
%First find closestEdge to center:
closestEdge = cCenter(1);
if  closestEdge > cCenter(2)
    closestEdge = cCenter(2);
end
if closestEdge > length(inputImage(:,1))-cCenter(2)
    closestEdge = length(inputImage(:,1))-cCenter(2);
end
if closestEdge > length(inputImage(1,:))-cCenter(1)
    closestEdge = length(inputImage(1,:))-cCenter(1);
end

closestEdge = floor(closestEdge);

if ~mod(closestEdge,2) %if even
    %To make sure even side length, closestEdge must be odd:
    closestEdge = closestEdge-1; 
end

squareCutImage = inputImage(cCenter(2)-closestEdge:cCenter(2)+closestEdge-1,...
    cCenter(1)-closestEdge:cCenter(1)+closestEdge-1);


%Fold square region into a quarter (average values):
qlength = size(squareCutImage,1);
qheight = size(squareCutImage,2);

if qlength/2 ~= qheight/2
    disp('Error sidelengths not equal');
end

squareLength = qlength/2;

quarteredImageSum = []; quarteredImageAvg = [];
for i=1:(squareLength) %y
    for j=1:(squareLength) %x
        quarteredImageSum(i,j) = squareCutImage(squareLength+i,squareLength+j)...
                            + squareCutImage(squareLength+i,squareLength-j+1)...
                            + squareCutImage(squareLength-i+1,squareLength-j+1)...
                            + squareCutImage(squareLength-i+1,squareLength+j);        
    end
end

quarteredImageAvg = quarteredImageSum./4; %4 quarters...


%Now bin on vector from top-left region:
rVectors = []; %Contains doubles of position

for i=1:length(quarteredImageAvg(:,1))
    for j=1:length(quarteredImageAvg(1,:))       
        rVectors(i,j) = sqrt(i*i + j*j);      
    end
end


%Binning:
reQuarteredImageAverage = reshape(quarteredImageAvg, 1, length(quarteredImageAvg(:,1))*length(quarteredImageAvg(1,:)));
reRVectors = reshape(rVectors, 1, length(rVectors(:,1))*length(rVectors(1,:)));

nbins = 100; binMean = [];
binEdges = linspace(min(min(rVectors)),max(max(rVectors)),nbins+1);
binLength = (max(max(rVectors))-min(min(rVectors)))/(nbins+1);

[h,whichBin] = histc(reRVectors, binEdges);

for i = 1:nbins
    flagBinMembers = (whichBin == i);
    binMembers     = reQuarteredImageAverage(flagBinMembers);
    binMean(i)     = mean(binMembers);
end

radialProfile = [];

for i = 1:length(binMean)
    radialProfile(1,i) = binMean(i);
    radialProfile(2,i) = i*binLength; %pixel vector
end

profileOutput = radialProfile;


end
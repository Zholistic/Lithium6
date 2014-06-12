function [profileOutput] = radAverageBigSquare(inputImage)
%A function which takes a 2D image and radially averages. 
%Returns the 1D profile of that image from center to edge. 

%First find the center of the cloud:
gcoefsX = gausFit1D(mean(inputImage,1)); %mean averages over y
gcoefsY = gausFit1D(mean(inputImage,2)); %mean averages over x

cCenter = [round(gcoefsX(2)), round(gcoefsY(2))]; %[x y]

%Take square region about center.
%First find furthestEdge to center:
if(0)
furthestEdge = cCenter(1);
if  furthestEdge < cCenter(2)
    furthestEdge = cCenter(2);
end
if furthestEdge < length(inputImage(:,1))-cCenter(2)
    furthestEdge = length(inputImage(:,1))-cCenter(2);
end
if furthestEdge < length(inputImage(1,:))-cCenter(1)
    furthestEdge = length(inputImage(1,:))-cCenter(1);
end

furthestEdge = floor(furthestEdge);

if ~mod(furthestEdge,2) %if even
    %To make sure even side length, furthestEdge must be odd:
    furthestEdge = furthestEdge-1;
end

%squareCutImage = inputImage(cCenter(2)-closestEdge:cCenter(2)+closestEdge-1,...
%    cCenter(1)-closestEdge:cCenter(1)+closestEdge-1);

%Inflate image to large square:
%Add buffer around each image:
originalSize = [length(inputImage(1,:)), length(inputImage(:,1))]; %[x y]

bigSquareImg = []; buffedImgArrayPre = [];
    paddingpre = [farthestEdge-cCenter(2) farthestEdge-cCenter(1)];%[y x] Padding depends on relative center
    paddingpost = [originalSize(2)-cCenter(2) ceil(originalSize(2)./2)-shift(1)];
    buffedImgArrayPre(:,:) = padarray(inputImage(:,:),paddingpre,0,'pre'); %Adds a border to each image of zeros
    bigSquareImg(:,:) = padarray(buffedImgArrayPre(:,:),paddingpost,0,'post');


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
rVectorsL = [];


for i=1:length(quarteredImageAvg(:,1))
    for j=1:length(quarteredImageAvg(1,:))       
        rVectors(i,j) = sqrt(i*i + j*j);      
    end
end

end

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
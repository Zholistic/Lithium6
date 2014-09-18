function [ imageOutput, domainOutput ] = binMeCenterAndAverageIncNaN( imgArray, inputX, nbins )
%Bins the image array inputY on the domain array inputX. Puts into nbins.
%Returns outputBinned of which outputBinned(1,:) is the range.
%scans over NaN values.

%inputY = imgArray(1,1,:);

binMean = []; binStd = []; binEdges=[]; imageOutput = []; finalBinned = [];
binEdgesR = [];
binLengths = []; binMembers = [];
binsStart = 2000; binsEnd = 60000;
binEdges = linspace(binsStart,binsEnd,nbins+1);
binLength = (binsEnd-binsStart)/(nbins);
[h,whichBin] = histc(inputX, binEdges);

j=1; binCount = 1;
for i = 1:nbins+1
    flagBinMembers = (whichBin == i);
    binMembers     = imgArray(:,:,flagBinMembers);
    %what bins to exclude? Nothing
    if (1)
        if ~isempty(binMembers)
            binMean(:,:,j) = centerAndAverage(binMembers);
        else
            binMean(:,:,j) = ones([length(imgArray(:,1,1)) length(imgArray(1,:,1))]).*(40000);
        end
        %binStd(j)      = std(binMembers);
        binLengths(j) = binCount*binLength;
        binEdgesR(j) = binEdges(i);
        j=j+1;
        binCount = 1;
        binMembers = [];
    else
        binCount = binCount +1;
    end
end

images = []; xOutput = [];
for i=1:length(binMean(1,1,:))
    images(:,:,i) = binMean(:,:,i);
    %Sits of left hand side of bin:
    %finalBinned(2,i) = binsStart - binLengths(1) + sum(binLengths(1:i)); 
    %offset half a binlength
    thisLength = binLengths(i);
    xOutput(i) = binsStart + (sum(binLengths(1:i))) - thisLength/2;    
    %finalBinned(3,i) = binStd(i); %error on the bin sums
    %finalBinned(4,i) = binEdgesR(i)+binLength;
end

imageOutput = images;
domainOutput = xOutput;

end
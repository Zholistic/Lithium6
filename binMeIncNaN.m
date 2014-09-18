function [ outputBinned ] = binMeIncNaN( inputY, inputX, nbins )
%Bins the range array inputY on the domain array inputX. Puts into nbins.
%Returns outputBinned of which outputBinned(1,:) is the range.
%scans over NaN values.

binMean = []; binStd = []; binEdges=[]; outputBinned = []; finalBinned = [];
binEdgesR = [];
binLengths = []; binMembers = [];
binsStart = 2000; binsEnd = 60000;
binEdges = linspace(binsStart,binsEnd,nbins+1);
binLength = (binsEnd-binsStart)/(nbins);
[h,whichBin] = histc(inputX, binEdges);

j=1; binCount = 1;
for i = 1:nbins+1
    flagBinMembers = (whichBin == i);
    binMembers     = inputY(flagBinMembers);
    %what bins to exclude? None
    if(1)
        binMean(j)     = mean(binMembers);
        binStd(j)      = std(binMembers);
        binLengths(j) = binCount*binLength;
        binEdgesR(j) = binEdges(i);
        j=j+1;
        binCount = 1;
        binMembers = [];
    else
        binCount = binCount +1;
    end
end


for i=1:length(binMean)
    finalBinned(1,i) = binMean(i);
    %Sits of left hand side of bin:
    %finalBinned(2,i) = min(inputX) - binLengths(1) + sum(binLengths(1:i)); 
    %offset half a binlength
    thisLength = binLengths(i);
    finalBinned(2,i) = binsStart + (sum(binLengths(1:i))) - thisLength/2;    
    finalBinned(3,i) = binStd(i); %error on the bin sums
    finalBinned(4,i) = binEdgesR(i)+binLength;
end

outputBinned = finalBinned;

end
function [ outputBinned ] = binMe( inputY, inputX, nbins )
%Bins the range array inputY on the domain array inputX. Puts into nbins.
%Returns outputBinned of which outputBinned(1,:) is the range.
%scans over NaN values.

binMean = []; binStd = []; binEdges=[]; outputBinned = []; finalBinned = [];
binEdgesR = []; binCounts = []; binSizes = [];
binLengths = []; binMembers = [];
binEdges = linspace(min(inputX),max(inputX),nbins+1);
binLength = (max(inputX)-min(inputX))/(nbins);
[h,whichBin] = histc(inputX, binEdges);

j=1; binCount = 1;
for i = 1:nbins+1
    flagBinMembers = (whichBin == i);
    binMembers     = inputY(flagBinMembers);
    binSizes(i) = length(binMembers);
    %what bins to exclude? NaN
    if ~isnan(binMembers)

        binMean(j)     = mean(binMembers);
        binStd(j)      = std(binMembers);
        binLengths(j) = binCount*binLength;
        binEdgesR(j) = binEdges(i);
        j=j+1;
        binCounts(i) = binCount;
        binCount = 1;
        binMembers = [];
    else
        binCount = binCount +1;
    end
    
    flagBinMembers = [];
    binMembers = [];
end


for i=1:length(binMean)
    finalBinned(1,i) = binMean(i);
    %Sits of left hand side of bin:
    %finalBinned(2,i) = min(inputX) - binLengths(1) + sum(binLengths(1:i)); 
    %offset half a binlength
    thisLength = binLengths(i);
    finalBinned(2,i) = min(inputX) + (sum(binLengths(1:i))) - thisLength/2;    
    finalBinned(3,i) = binStd(i); %error on the bin sums
    finalBinned(4,i) = binEdgesR(i)+binLength;
end

outputBinned = finalBinned;

end
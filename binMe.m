function [ outputBinned ] = binMe( inputY, inputX, nbins )
%Bins the range array inputY on the domain array inputX. Puts into nbins.
%Returns outputBinned of which outputBinned(1,:) is the range.
%scans over NaN values.

binMean = []; binStd = []; binEdges=[]; outputBinned = []; finalBinned = [];
binEdges = linspace(min(inputX),max(inputX),nbins+1);
binLength = (max(inputX)-min(inputX))/(nbins+1);
[h,whichBin] = histc(inputX, binEdges);

j=1; binCount = 1;
for i = 1:nbins
    flagBinMembers = (whichBin == i);
    binMembers     = inputY(flagBinMembers);
    %what bins to exclude? NaN
    if ~isnan(mean(binMembers))
        binMean(j)     = mean(binMembers);
        binStd(j)      = std(binMembers);
        binLengths(j) = binCount*binLength;
        j=j+1;
        binCount = 1;
    else
        binCount = binCount +1;
    end
end


for i=1:length(binMean)
    finalBinned(1,i) = binMean(i);
    finalBinned(2,i) = sum(binLengths(1:i));
    finalBinned(3,i) = binStd(i); %error on the bin sums
end

outputBinned = finalBinned;

end
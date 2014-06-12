function [ scaleFactor ] = makeScalefactorHL( highIntImg, lowIntImg )
%From a high intensity image representing the 'real' atom number and a low
%intensity image both expressed in optical density (OD) create a comparison
%spectrum that is OD(H)/OD(L) vs OD(L). 

%Reshape to 1D and sort the images:
reshapedHighIntImg = reshape(highIntImg, 1, length(highIntImg(:,1))*length(highIntImg(1,:)));
reshapedLowIntImg = reshape(lowIntImg, 1, length(lowIntImg(:,1))*length(lowIntImg(1,:)));

sortedHighInt = sort(reshapedHighIntImg);
sortedLowInt = sort(reshapedLowIntImg);

sortedHighInt = sortedHighInt+abs(sortedHighInt(1)); %Shifts so everything > 0
sortedLowInt = sortedLowInt+abs(sortedLowInt(1));

%minOD = min(min(sortedHighInt),min(sortedLowInt));

%Now need to bin on optical density:
nbins = 100; binIdxHigh = []; binIdxLow = [];
binArrayHigh = linspace(min(sortedHighInt),max(sortedHighInt),nbins+1);
[nHigh,binIdxHigh] = histc(sortedHighInt, [binArrayHigh(1:end-1) Inf]);
%nj = accumarray(binIdx, 1, [nbins 1], @sum);

binArrayLow = linspace(min(sortedLowInt),max(sortedLowInt),nbins+1);
[nLow,binIdxLow] = histc(sortedLowInt, [binArrayLow(1:end-1) Inf]);
%nj = accumarray(binIdx, 1, [nbins 1], @sum);

dividedArray = binArrayHigh./binArrayLow;
scaleFactor = mean(dividedArray(10:end-10)); %move 10 each way into array to avoic NaN's and Infs etc



end


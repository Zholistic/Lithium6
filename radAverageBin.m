function [profileOutput,center] = radAverageBin(inputImage)
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

%[sortedImageData,indexs] = sort(reImageAverage,'descend');
%sortedreRVectors= reRVectors(indexs);

[sortedreRVectors,indexs] = sort(reRVectors,'ascend');
sortedImageData = reImageAverage(indexs);

%plot(sortedreRVectors,sortedImageData)
%plot(meanNelements(sortedreRVectors,5),meanNelements(sortedImageData,5),'.')
%n = 5; plot(meanNelements(sortedreRVectors,n),meanNelements(sortedImageData,n),'.');
%figure(2); plot(sortedImageData)
%figure(3); plot(sortedreRVectors)
%figure(4); plot(sortedreRVectors,sortedImageData); hold on; 

%binResult = binMe(sortedImageData,sortedreRVectors,500);
%plot(binResult(2,:),binResult(1,:));

if(0)
figure(1000);
subplot(1,5,1);
n = 1; plot(meanNelements(sortedreRVectors,n),meanNelements(sortedImageData,n));
axis([0 90 0 2.2]);
title([ num2str(n) ' averaged']);

subplot(1,5,2);
n = 10; plot(meanNelements(sortedreRVectors,n),meanNelements(sortedImageData,n));
axis([0 90 0 2.2]);
title([ num2str(n) ' averaged']);

subplot(1,5,3);
n = 50; plot(meanNelements(sortedreRVectors,n),meanNelements(sortedImageData,n));
axis([0 90 0 2.2]);
title([ num2str(n) ' averaged']);

subplot(1,5,4);
n = 100; plot(meanNelements(sortedreRVectors,n),meanNelements(sortedImageData,n));
axis([0 90 0 2.2]);
title([ num2str(n) ' averaged']);

subplot(1,5,5);
n = 500; plot(meanNelements(sortedreRVectors,n),meanNelements(sortedImageData,n));
axis([0 90 0 2.2]);
title([ num2str(n) ' averaged']);
end

radialProfile(1,:) = sortedImageData;
radialProfile(2,:) = sortedreRVectors;

%radialProfile(1,:) = binResult(1,:);
%radialProfile(2,:) = binResult(2,:);

profileOutput = radialProfile;
center = cCenter;


end
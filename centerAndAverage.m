function [ outputImg ] = centerAndAverage( inputImgArray )
%Takes an array of images and finds the center of each using gaussian fits;
%then shifts each to be centered over a common center and finally averages.
%Returns a single image which is the averaged+centered image.
disp('Centering and Averging...');

imgNumber = length(inputImgArray(1,1,:));

centeredImgArray = centerImgArray(inputImgArray);

outputImg = mean(centeredImgArray,3);


end


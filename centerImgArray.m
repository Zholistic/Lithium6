function [ outputImgArray ] = centerImgArray( inputImgArray )
%Shift each image in the input image array to be over the center of the 
%first image. First find centers by fitting gaussians, then create a larger
%ROI around the image and center, then strip off the extra ROI. 

disp('Centering Image Array...');

%Find centers:
gcoefsX = []; gcoefsY = []; centers = [];
for i=1:length(inputImgArray(1,1,:))    
    gcoefsX(:,i) = gausFit1D(mean(inputImgArray(:,:,i),1)); %mean averages over y
    gcoefsY(:,i) = gausFit1D(mean(inputImgArray(:,:,i),2)); %mean averages over x
    
    centers(:,i) = [ceil(gcoefsX(2,i)), ceil(gcoefsY(2,i))]; %center = [x y]
    %centers(:,i)
end

%Spike: center to center on
%TODO convert this to median of the centers array
pin = centers(:,1);

%padbuffer
padbuffer = 100;

%Add buffer around each image:
originalSize = [length(inputImgArray(1,:,1)), length(inputImgArray(:,1,1))]; %[x y]

buffedImgArray = []; buffedImgArrayPre = [];
for i=1:length(inputImgArray(1,1,:))
    buffedImgArrayPre = [];
    shift(:,i) = pin - centers(:,i);
    %paddingpre = [30+shift(2,i) 30+shift(1,i)];%[y x] Padding depends on relative center
    %paddingpost = [30-shift(2,i) 30-shift(1,i)];
    paddingpre = [padbuffer+shift(2,i) padbuffer+shift(1,i)];%[y x] Padding depends on relative center
    paddingpost = [padbuffer-shift(2,i) padbuffer-shift(1,i)];
    buffedImgArrayPre(:,:) = padarray(inputImgArray(:,:,i),paddingpre,0,'pre'); %Adds a border to each image of zeros
    buffedImgArray(:,:,i) = padarray(buffedImgArrayPre(:,:),paddingpost,0,'post');
end

%Cut out final images:
for i=1:length(buffedImgArray(1,1,:))
    centeredImgArray(:,:,i) = buffedImgArray(padbuffer:padbuffer+originalSize(2)-1,padbuffer:padbuffer+originalSize(1)-1,i);  
end


outputImgArray = centeredImgArray;

end


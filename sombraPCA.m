directory = 'C:\Users\tpeppler\Downloads\Sombra\GoodImages\';
ImageFiles = dir([directory '/*.*']);
imageArray = [];
for i = 1:length(ImageFiles)-3
I = imread([directory ImageFiles(i+3).name]);
I2 = double(I);
%hi = imwrite(I);
imageArray(:,:,i) = I2(:,:,1);
end

imageArrayCrop = imageArray(55:290,111:527,:);

img = imread('C:\Users\tpeppler\Downloads\Sombra\GoodImages\img00177.png'); % Read image
red = img(:,:,1); % Red channel
green = img(:,:,2); % Green channel
blue = img(:,:,3); % Blue channel
a = zeros(size(img, 1), size(img, 2));
just_red = cat(3, red, a, a);
just_green = cat(3, a, green, a);
just_blue = cat(3, a, a, blue);
back_to_original_img = cat(3, red, green, blue);
figure, imagesc(img), title('Original image')
figure, imagesc(just_red), title('Red channel')
figure, imagesc(just_green), title('Green channel')
figure, imagesc(just_blue), title('Blue channel')
figure, imagesc(back_to_original_img), title('Back to original image')








        pcaEndPoint = 24; %Where to take the final holdtime array index

        imagesToReshape = [];
        %imagesToReshape = imageArrayC(:,:,1:35); %single images
        imagesToReshape = imageArrayCrop(:,:,1:end);
        %imagesToReshape = imagesAvgSet2;
        %imagesToReshape = imagesAvgSet3;
        %imagesToReshape = imagesAvgSet4;
        
        numberOfImages = length(imagesToReshape(1,1,:));
        vectorLength = length(imagesToReshape(:,1,1))*length(imagesToReshape(1,:,1));
        
        %reshape images into vectors:
        vectorImages = []; vectorSum = zeros(vectorLength,1);
        for i=1:length(imagesToReshape(1,1,:))
            vectorImages(:,i) = reshape(imagesToReshape(:,:,i),[],1);
            vectorSum = vectorSum + vectorImages(:,i);
        end
        
        %Mean Image:
        meanImageVector = (1/numberOfImages).*vectorSum;
        undividedImage = reshape(vectorSum,length(imagesToReshape(:,1,1)),length(imagesToReshape(1,:,1)));
        meanImage = reshape(meanImageVector,length(imagesToReshape(:,1,1)),length(imagesToReshape(1,:,1)));
        
        normImages = [];
        for i=1:length(imagesToReshape(1,1,:))
            normImages(:,i) = vectorImages(:,i) - meanImageVector;
        end
        
        %N x p Bmatrix (N = number of images, p = pixels per image)
        Bmatrix = normImages;
        
        %[coeff, score, latent] = pca(ingredients)
        %each column of score corresponds to one principal component
        %latent stores the variances of the N principal components
        [Y01,P01,E01,tsquared,percentV] = pca(Bmatrix);
              
        
        %Reconstruct images using desired principle components
        eigenVectorBmatrix = []; imagesFromEVector = [];
        for i=1:numberOfImages
            eigenVectorBmatrix = P01(:,i);
            imagesFromEVector(:,:,i) = reshape(eigenVectorBmatrix,length(imagesToReshape(:,1,1)),length(imagesToReshape(1,:,1)));
        end
        
    
                
        for i=1:5
            figure(i+100); imagesc(imagesFromEVector(:,:,i));
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
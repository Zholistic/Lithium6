function [beamReturnArray,atom1ReturnArray,atom2ReturnArray] = PullSPEBatch(fileloclist,Isat,OD,fringeRem)
%fileloclist is a cell of the image location names...

%SPE_File_Name = 'C:\Data\140311_Imaging_Parameter_Check_2D_Trap_1usROI1_5usROI2_10usROI4_10usNewBeatROI5_1us20IsatROI3_Image_Pulse\140310_2748.SPE';
atomImages1 = []; atomImages2 = []; beamImages1 = []; beamImages2 = [];
for i=1:length(fileloclist)
SPE_File_Name = fileloclist{i};

debug = 0;

warning('OFF');

xDim = 256;
yDim = 768;
NoFrames = 3;

Phi       =0.93; %Beam correction term determined for each beam.
PixelArea =(13e-6*83.0/400.0)^2; %Actual Pixel size area
%PixelArea =(2.84e-6)^2; %Macus given pxarea
Sigma     =3*((671e-9)^2)/(2*pi); %Resonant Cross Section

Delta     =0; %Detuning
Gamma     =5.9e6*(2*pi); %Natural Linewidth

[NumberOfImages,~] =size(SPE_File_Name);
%Read in SPE header:
fid    =fopen(SPE_File_Name);
header =fread(fid,4100/2,'int16'); 
fclose(fid);

%Mine header for information:
%NoFrames =header(741);
%xDim     =header(758);
%yDim     =header(761); 

%Inits:
ImMat    =zeros(xDim*yDim*NoFrames,1);
ZIm      =zeros(xDim,yDim,1);
ZBg      =zeros(xDim,yDim,1);
RawImage =zeros(xDim,yDim,NumberOfImages);
RawBackg =zeros(xDim,yDim,NumberOfImages);
RawImageNoBgSub =zeros(xDim,yDim,NumberOfImages);

%Open and allocate .SPE file data:

k=1;
    picx_name_spe =SPE_File_Name;
    fid           =fopen(picx_name_spe);
    header(:)     =fread(fid,4100/2,'int16'); %Read in first 4100 bytes
    ImMat(:)      =fread(fid,xDim*yDim*NoFrames,'int16'); %3 frames/image
    ZIm(:,:)      =reshape(ImMat((xDim*yDim)+1:2*(xDim*yDim)),xDim,yDim);
    if NoFrames==3
        ZBg(:,:)            =reshape(ImMat(2*(xDim*yDim)+1:end),xDim,yDim);
        %ZBg = medfilt2(ZBg,[3,3]); %medfilt the background. COULD BE SLOW!
        RawImage(:,:,k)     =double(ZIm(:,:)-ZBg(:,:)); %RawImage already has background subtraction.
        RawImageNoBgSub(:,:,k) =double(ZIm(:,:));
        RawBackg(:,:,k)     =double(ZBg(:,:));
    else
        disp(['Image ' SPE_File_Name ' does not contain 3 frames.']);
    end
    fclose(fid);
 
%disp(SPE_File_Name);

MeanBg=mean(RawBackg,3); %average of all backgrounds

%parfor i=1:size(RawImageNoBgSub,3)
%    RawImage(:,:,i) =RawImageNoBgSub(:,:,i)-MeanBg; %Subtract average of all backgrounds.
%end

%Cut out a better image:
ImagesCrop = [];
CropY = 35:251; CropX = 30:236;
crop1 = RawImage(35:241,30:236,1);
crop2 = RawImage(35:241,(30+256):(236+256),1);
crop3 = RawImage(35:241,(30+512):(236+512),1);
%crop1 = RawImage(15:245,25:248,1);
%crop2 = RawImage(15:245,(25+256):(248+256),1);
%crop3 = RawImage(15:245,(25+512):(248+512),1);
%ImagesCrop(:,1:175,1) = crop1;
%ImagesCrop(:,176:(175*2),1) = crop2;
%ImagesCrop(:,351:525,1) = crop3;

%figure(2)
%image(ImagesCrop,'CDataMapping','scaled');
%image(crop2(JunkTop:JunkBottom,JunkLeft:JunkRight,1),'CDataMapping','scaled');

BeamImg = crop1;
AtomImg1 = crop2;
AtomImg2 = crop3;

%Just use normal Raws:
%BeamImg = RawImage(:,1:yDim/NoFrames,1);
%AtomImg1 = RawImage(:,yDim/NoFrames + 1:(yDim/NoFrames)*2,1); 
%AtomImg2 = RawImage(:,(yDim/NoFrames)*2 + 1:end,1); 

%Re-order images from high to low, in order to find BeamScaleFactor:
SortBeamImg = sort(BeamImg(:,:,1));
SortAtomImg1 = sort(AtomImg1(:,:,1));
SortAtomImg2 = sort(AtomImg2(:,:,1));

SortAvgBeam = []; SortAvgAtomImg1 = []; SortAvgAtomImg2 = [];

%Sums over the x direction for each sorted image to give an average
%of largest and smallest values. 
for j=1:length(crop1(:,1,1))
    SortAvgBeam(j) = sum(SortBeamImg(j,:,1));
    SortAvgAtomImg1(j) = sum(SortAtomImg1(j,:,1));
    SortAvgAtomImg2(j) = sum(SortAtomImg2(j,:,1));
end

%BeamScaleFactor = (SortAvgBeam(end) - SortAvgBeam(1))/(SortAvgAtomImg1(end) - SortAvgAtomImg1(1));
BeamScaleFactor2 = (mean(SortAvgBeam(end-30:end-10))./mean(SortAvgAtomImg2(end-30:end-10)));
BeamScaleFactor1 = (mean(SortAvgBeam(end-30:end-10))./mean(SortAvgAtomImg1(end-30:end-10))); %Found on maximum intensities
%BeamScaleFactor = 0.93;
%disp(num2str(BeamScaleFactor2));
%disp(num2str(BeamScaleFactor1));

%Scale the beam:
BeamImg2 = BeamImg./BeamScaleFactor2;
BeamImg1 = BeamImg./BeamScaleFactor1; %Scale Beam Image down
%AtomImg1 = AtomImg1.*BeamScaleFactor1; %Scale Atom Image up
%AtomImg2 = AtomImg2.*BeamScaleFactor2;

%Beam scale factor off:
%BeamImg2 = BeamImg;
%BeamImg1 = BeamImg; %Scale Beam Image down

if(debug)
    ImagesBeforeOD = [];
    %ImagesBeforeOD(:,1:175,1) = BeamImg;
    %ImagesBeforeOD(:,176:(175*2),1) = AtomImg1;
    %ImagesBeforeOD(:,351:525,1) = AtomImg2;
    
    ImagesBeforeOD(:,1:yDim/NoFrames,1) = BeamImg;
    ImagesBeforeOD(:,yDim/NoFrames + 1:(yDim/NoFrames)*2,1) = AtomImg1;
    ImagesBeforeOD(:,(yDim/NoFrames)*2 + 1:yDim,1) = AtomImg2;
    
    figure(50)
    image(ImagesBeforeOD,'CDataMapping','scaled');
    
    DivImage2 = AtomImg2./(BeamImg);
    figure(1)
    imagesc(DivImage2);
    figure(2)
    imagesc(real(-log(DivImage2)));
    figure(100)
    imagesc(real(-log(DivImage2).*(PixelArea/Sigma)));
end

beamImages2(:,:,i) = BeamImg2;
beamImages1(:,:,i) = BeamImg1;
atomImages2(:,:,i) = AtomImg2;
atomImages1(:,:,i) = AtomImg1;
end

if(fringeRem)
    %Do fringe removal and return these:
backgroundMask = [];
for m = 1:207
    for n = 1:207
        if(m > 55 && m < 160 && n > 26 && n < 132)
        backgroundMask(n,m) = 0; %atoms region
        else
        backgroundMask(n,m) = 1; %background region
        end
    end
end


%[ odimages,optrefimages,avgimage,timer ] = fringeremoval2( absimages,refimages,bgmask,returnimgs)
%[ odimages,optrefimages,avgimage,timer ] = fringeremoval2( atomImages2,beamImages2,backgroundMask,returnimgs)
%[ odimages,optrefimages,avgimage,timer ] = fringeremoval2( atomImages2,beamImages2,backgroundMask)
if(0)
    close all;
    for j=1:length(odimages(1,1,:))
        if(mod(j,10) == 0)
            figure(j);
            imagesc(real(odimages(:,:,j)));
        end
    end
end
%-log(abs-dark)/(ref-dark)

beamReturnArray = beamf;
atom1ReturnArray = atom1f;
atom2ReturnArray = atom2f;

else
for i=1:length(fileloclist)
%Optical Density Calc-----------------------------------------------------%
%Division Term:
AtomsDiv1 = (1+4*(Delta^2)/Gamma^2).*real(-log(AtomImg1./(BeamImg1)));
AtomsDiv2 = (1+4*(Delta^2)/Gamma^2).*real(-log(AtomImg2./(BeamImg2)));

%Subtraction Term:
AtomsSub1 = -(AtomImg1 - (BeamImg1))./Isat;
AtomsSub2 = -(AtomImg2 - (BeamImg2))./Isat;

if(debug)
    figure(3)
    imagesc(AtomsSub2);
    figure(4)
    imagesc(AtomsDiv2);
end

%figure(1)
%image(AtomsSub2,'CDataMapping','scaled');

%Combine
%OD1 = (1/Sigma).*(AtomsDiv1 + AtomsSub1);
%OD2 = (1/Sigma).*(AtomsDiv2 + AtomsSub2);
%OD1 = (AtomsDiv1);
%OD2 = (AtomsDiv2);
OD1 = (AtomsDiv1 + AtomsSub1);
OD2 = (AtomsDiv2 + AtomsSub2);

%figure(2)
%image(OD2,'CDataMapping','scaled');

Atoms1ODnn = OD1;
Atoms2ODnn = OD2;


%Noise Correction --------------------------------------------------------%
%Use a region (not in the cloud) to determine the noise character:
JunkLeft = 20; JunkRight = 140; JunkTop = 140; JunkBottom = 180;
PxAreaNumber = (JunkRight - JunkLeft)*(JunkBottom - JunkTop); %Number of pixels in junk area.



AvgNoisePixel1 = sum(sum(OD1(JunkTop:JunkBottom,JunkLeft:JunkRight,1)))./PxAreaNumber; 
AvgNoisePixel2 = sum(sum(OD2(JunkTop:JunkBottom,JunkLeft:JunkRight,1)))./PxAreaNumber; 

%AvgNoisePixel2
%if(0) %Background subtraction off!
Atoms1ODnn = OD1 - AvgNoisePixel1;
Atoms2ODnn = OD2 - AvgNoisePixel2;
%end

if(debug)
    figure(5)
    imagesc(Atoms2ODnn);
end

%Gradient correction -----------------------------------------------------%

%Determine slice to fit over:
ToFit = []; xs = []; coefs = [];
ToFit(:,1) = mean(Atoms2ODnn(140:170,:,1));

%Fit:
fg = @(p,x)(p(1).*x + p(2)); %function to fit with
p0 = [0.0022 -0.5];
lb = [0 -5];
ub = [0.8 5];

curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
xs = 1:length(ToFit(:,1));
coefs = lsqcurvefit(fg,p0,xs(:),ToFit(:,1),lb,ub,curvefitoptions);

%if(coefs(1) > 0.35 || abs(coefs(2)) < 0.001)
%    disp('check gradient correction');
%    coefs
%end

%Now subtract out value based on the gradient:
%if(0) %Gradient correction off!
for k=1:length(Atoms2ODnn(:,1,1))
    Atoms2ODnn(k,:,1) = Atoms2ODnn(k,:,1) - fg(coefs,xs);
end
%end

    
if(debug)
    figure(10)
    hold on;
    plot(ToFit);
    plot(xs,fg(coefs,xs),'r');
    hold off;
end

%Without Optical Density:
%Atoms1ODnn = AtomImg1 - BeamImg;
%Atoms2ODnn = AtomImg2 - BeamImg;


%Invert so relevant values are positive:
%Atoms1ODnn = Atoms1ODnn.*(-1);
%Atoms2ODnn = Atoms2ODnn.*(-1);

if(debug)
    figure(6)
    imagesc(Atoms2ODnn);
end

%Convert from OD to real atoms:
Atoms1RN = Atoms1ODnn.*(PixelArea/Sigma);
Atoms2RN = Atoms2ODnn.*(PixelArea/Sigma);


%Vars to return from function:
%beamf = (BeamImg/BeamScaleFactor);


if(OD)
    beamf(:,:,i) = BeamImg;
    atom1f(:,:,i) = Atoms1ODnn;
    atom2f(:,:,i) = Atoms2ODnn;
else
    beamf(:,:,i) = BeamImg;
    atom1f(:,:,i) = Atoms1RN;
    atom2f(:,:,i) = Atoms2RN;   
end
end

beamReturnArray = beamf;
atom1ReturnArray = atom1f;
atom2ReturnArray = atom2f;

end

end


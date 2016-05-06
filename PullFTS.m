function [ftsImage] = PullFTS(fileloc,raw,Isat)
%'raw' uses the atom and beam images
Delta     =0; %Detuning
Gamma     =5.9e6*(2*pi); %Natural Linewidth

if(raw)
    %Rejig the strings:
    tempString = fileloc(1:end-7);
    number = fileloc(end-6:end-4);
    if(str2num(number(1)) == 0)

        if(str2num(number(2)) == 0)
        number = number(3:end);
        else
            number = number(2:end);
        end
                
    end
    atomloc = [tempString 'atom' number '.fts'];
    beamloc = [tempString 'beam' number '.fts'];
    
    newImageAtom = imread(atomloc);
    newImageBeam = imread(beamloc);
    
    %Isat Correction?
    newImage = (-1)*(double(newImageAtom) - double(newImageBeam));
    
    %Isat Correction%%%%%%%%%%
    %Division Term:
    AtomsDiv1 = (1+4*(Delta^2)/Gamma^2).*real(-log(double(newImageAtom)./(double(newImageBeam))));

    %Subtraction Term:
    AtomsSub1 = -(newImageAtom - (newImageBeam))./Isat;

    OD1 = (AtomsDiv1 + double(AtomsSub1));
    Atoms1ODnn = OD1;
    
    %Noise cancellation, using two large strips above and below the cloud:
    NROIy1 = 200:400; NROIy2 = 700:800;
    NROIx1 = 200:1000; NROIx2 = 200:1000;
    
    noiseVal1 = mean(mean(newImage(NROIy1,NROIx1)));
    noiseVal2 = mean(mean(newImage(NROIy2,NROIx2)));
    
    noiseVal = (noiseVal1 + noiseVal2)/2;
    
    %Subtraction of noise:
    finalImage = newImage - noiseVal;
    ftsImage = Atoms1ODnn;
    ftsImage = finalImage;
    
else
    
    %Read in image:
    newImage = imread(fileloc);
    
    %Noise cancellation, using two large strips above and below the cloud:
    NROIy1 = 200:400; NROIy2 = 700:800;
    NROIx1 = 200:1000; NROIx2 = 200:1000;
    
    noiseVal1 = mean(mean(newImage(NROIy1,NROIx1)));
    noiseVal2 = mean(mean(newImage(NROIy2,NROIx2)));
    
    noiseVal = (noiseVal1 + noiseVal2)/2;
    
    %Subtraction of noise:
    finalImage = newImage - noiseVal;
    
    
    ftsImage = finalImage;
    
end

end

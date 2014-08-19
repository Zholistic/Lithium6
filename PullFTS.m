function [ftsImage] = PullFTS(fileloc,raw)
%'raw' uses the atom and beam images

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
    newImage = (-1)*(newImageAtom - newImageBeam);
    
    %Noise cancellation, using two large strips above and below the cloud:
    NROIy1 = 200:400; NROIy2 = 700:800;
    NROIx1 = 200:1000; NROIx2 = 200:1000;
    
    noiseVal1 = mean(mean(newImage(NROIy1,NROIx1)));
    noiseVal2 = mean(mean(newImage(NROIy2,NROIx2)));
    
    noiseVal = (noiseVal1 + noiseVal2)/2;
    
    %Subtraction of noise:
    finalImage = newImage - noiseVal;
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

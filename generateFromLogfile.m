function [outputFileLocList,outputVarData] = generateFromLogfile(directory,datestring,varstring,camera)
%generateFileLocList: given a directory generates the list of filenames
%within it from the log file:

%Load Images: This section is filename format dependent.

sidecam = 0; topcam = 1; %Defaults to topcam
if(strcmp(camera,'top'))
    topcam = 1;
else
    sidecam = 1;
end

%Read in the log file:
disp('Reading in log file...');

logfilename = [directory datestring '_log_camera.txt'];
logfilename2 = [directory datestring '_log.txt'];
fid = fopen(logfilename,'rt');
if(fid == -1)
    fid = fopen(logfilename2,'rt');
end

C = textscan(fid, '%s', 'Delimiter','\t'); %tokenize into tab seperated tokens
C = C{1};
fclose(fid);

%Iterate through strings in the log file array 1-d C{:} to find the
%information we want:
varData = [];
TotalImages = 0;

%for loop to find the total number of images
for i=1:(length(C)-1)

    curr = C{i};
    next = C{i+1};
    
    if strcmp(curr,'MeasNr')
        TotalImages = TotalImages + 1;
    end
              
    
end


fileLocList = cell(TotalImages,1);
prev = 0;
j=0; k=0; picno = '0';
widthx = 0; widthy = 0; NsumROI = 0; centery = 0;
dataset = 1;

%for loop to get information out of the logfile and load it into arrays:
for i=1:(length(C)-1)
    
    if i>1
        prev = C{i-1};
    end
    curr = C{i};
    next = C{i+1};
    %ROImax = C{i+20};
    ROImax = '1';
    
    %load sidecam images
    if(sidecam)
        if strcmp(curr,'MeasNr')
            if str2num(ROImax) ~= 0
                j = j+1;
                filename = prev;
                picno = ['_' filename(end-2:end)];
                %picno = ['_' next];
                fullname = [directory filename picno '.fts'];
                %fullname = [directory filename '_001' '.fts'];
                %ftsImg = imread(fullname);
                fileLocList{j} = fullname;
                %imgArrayFresh(:,:,j) = ftsImg;
            else
                disp('fail on ROImax');
            end
        end
    end
    
    %load topcam images
    if(topcam)
        if strcmp(curr,'#WinView#')
            j = j+1;
            filename = next;
            fileLocList{j} = [directory filename];
            %[beamImage,atom1Image,atom2Image] = PullSPE(directory,filename);
            %imgAtom2Fresh(:,:,j) = atom2Image(:,:,1); %state 2 images in atom2Image
        end
    end
    
    %data set identifier
    %dataset = 1;
    if strcmp(curr,'width_mum_y1')
        dataset = 1;
    end
    if strcmp(curr,'width_mum_y2')
        dataset = 2;
    end
    if strcmp(curr,'width_mum_y3')
        dataset = 3;
    end
    if strcmp(curr,'width_mum_y4')
        dataset = 4;
    end
    if strcmp(curr,'width_mum_y5')
        dataset = 5;
    end
    if (strcmp(curr,'centery_m') || strcmp(curr,'COMy') || strcmp(curr,'centery_m14'))
        centery = next;
    end
    if strcmp(curr,'SIGMAx')
        widthx = next;
    end
    if strcmp(curr,'SIGMAy')
        widthy = next;
    end
    if strcmp(curr,'NsumROI')
        NsumROI = next;
    end
    if(topcam)
        if (strcmp(curr,'N_ROI') || strcmp(curr,'N_ROI1'))
            NsumROI = next;
        end
    end
    
    
    if strcmp(curr,varstring)
        k=k+1;
        %logholdtime(ceil(i/(length(C)/TotalImages))) = str2num(next);
        if(NsumROI)
        varData(k,5) = str2num(NsumROI); 
        else
          varData(k,5) = 0;   
        end
        
        if(widthy)
        varData(k,4) = str2num(widthy);
        else
         varData(k,4) = 0; 
        end
        
        if(widthx)
        varData(k,3) = str2num(widthx);
        else
         varData(k,3) = 0; 
        end
        
        if(centery)
            varData(k,6) = str2num(centery);
        else
            varData(k,6) = 0;
        end

        varData(k,2) = dataset;
        varData(k,1) = str2num(next);


        
    end
    
end

outputFileLocList = fileLocList;
outputVarData = varData;

end


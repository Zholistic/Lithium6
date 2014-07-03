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
fid = fopen(logfilename,'rt');
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
j=0; k=0;

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
    dataset = 1;
    if strcmp(curr,'N_ROI')
        dataset = 1;
    elseif strcmp(curr,'N_ROI2')
        dataset = 2;
    elseif strcmp(curr,'N_ROI3')
        dataset = 3;
    elseif strcmp(curr,'N_ROI4')
        dataset = 4;
    elseif strcmp(curr,'N_ROI5')
        dataset = 5;
    end
    
    if strcmp(curr,varstring)
        k=k+1;
        %logholdtime(ceil(i/(length(C)/TotalImages))) = str2num(next);
        varData(k,2) = dataset;
        varData(k,1) = str2num(next);
    end
    
end

outputFileLocList = fileLocList;
outputVarData = varData;

end

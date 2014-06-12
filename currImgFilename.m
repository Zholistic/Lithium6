
function [fn] = currImgFilename()

datetime=datestr(now);
datetime=strrep(datetime,':','_'); %Replace colon with underscore
datetime=strrep(datetime,'-','_');%Replace minus sign with underscore
datetime=strrep(datetime,' ','_');%Replace space with underscore

name = strcat('C:\MatlabImagesBackup\fig',datetime);

fn = name;
return
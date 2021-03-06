directory = 'C:\Users\tpeppler\Dropbox\PhD\2D_2016\GGTemp\Datas2\';
directorylist = dir(directory);

filenumber = length(directorylist)-2;

BetaEBCells = []; TempCells = []; ChemPotCells = [];
TempArray = []; ChemPotArray = []; BetaEBArray = [];
tempCount = 1; muCount = 1; betaEBCount = 1;
for i=1:length(directorylist)-2
    
filename = directorylist(i+2).name

%BetaEBCells{i} = str2num(filename(4:5))./100;

datafilename = [directory filename];

fid = fopen(datafilename,'rt');
%C = textscan(fid, '%s', 'Delimiter','\t'); %tokenize into tab seperated tokens
if i==1
C = textscan(fid,'%f %f %f %f');
elseif i==12
C = textscan(fid,'%f %f %f %f %f %f %f %f %f %f');   
elseif i>12
C = textscan(fid,'%f %f %f %f');
else
C = textscan(fid,'%f %f %f');
end
fclose(fid);


TempCells{i} = C{1};
for k=1:length(C{1})
   TempArray(tempCount) = C{1}(k);
   tempCount = tempCount + 1;
end

ChemPotCells{i} = C{3};
for k=1:length(C{3})
   ChemPotArray(muCount) = C{3}(k);
   muCount = muCount + 1;
end

for j=1:length(TempCells{i})
    BetaEBCells{i}(j) = str2num(filename(4:5))./100;
    BetaEBArray(betaEBCount) = str2num(filename(4:5))./100;
    betaEBCount = betaEBCount +1;
end

end



%plot3(x,y,z)
close all;
for i=1:filenumber
    
hold on;
%plot3(1./(2*pi*TempCells{i}),ChemPotCells{i}./TempCells{i},BetaEBCells{i},'.');
plot3(ChemPotCells{i}./TempCells{i},BetaEBCells{i},1./(2*pi*TempCells{i}),'.');
grid on;

end
xlabel('Beta mu') % x-axis label
ylabel('Beta Eb') % y-axis label
zlabel('n (1/ (2\pi T)') % y-axis label

%-------------------------------------------
x = []; y = []; z = []; xlin = []; ylin = []; Z = []; X = []; Y = [];

x = ChemPotArray./TempArray;
y = BetaEBArray;
z = 1./(2.*pi.*TempArray);

xlin = linspace(min(x),max(x),1000);
ylin = linspace(min(y),max(y),1000);
        

[X,Y] = meshgrid(xlin,ylin);

Z = griddata(x,y,z,X,Y,'cubic'); %v4 cubic linear



surf(X,Y,Z,'LineStyle','none','FaceColor','interp','FaceAlpha',0.6);
axis([0 5 0 10 0 2]); grid on;
xlabel('Beta mu') % x-axis label
ylabel('Beta Eb') % y-axis label
zlabel('n (1/ (2\pi T)') % y-axis label
view([1 0 0]);

if(0)
x = []; y = []; z = []; xlin = []; ylin = []; Z = [];
x = reshape(xdatas, 1, length(xdatas(:,1))*length(xdatas(1,:)));
y = reshape(ydatas, 1, length(ydatas(:,1))*length(ydatas(1,:)));
z = reshape(zdatas, 1, length(zdatas(:,1))*length(zdatas(1,:)));
%y = masterMatrix(2,:);
%z = masterMatrix(3,:).*2.*pixelLength./1e-6;

%reshape(logkfa2ds, 1, length(logkfa2ds(:,1))*length(logkfa2ds(1,:)));

xlin = linspace(min(x),max(x),80);
ylin = linspace(min(y),max(y),80);
        

[X,Y] = meshgrid(xlin,ylin);

%Z = griddata(x,y,z,X,Y,'cubic');
Z = griddata(x,y,z,X,Y,'cubic'); %v4 cubic linear
end


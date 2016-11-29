%Supergraph 3D surface plot:
close all;
open('C:\Users\tpeppler\Dropbox\PhD\2D_2016\EOS_Data\972_scatterGraph3D_Iso_indv clouds smoothed.fig');
handles(1) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2D_2016\EOS_Data\920_scatterGraph3D_Iso_indv clouds smoothed.fig');
handles(2) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2D_2016\EOS_Data\880_scatterGraph3D_Iso_indv clouds smoothed.fig');
handles(3) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2D_2016\EOS_Data\865_scatterGraph3D_Iso_indv clouds smoothed.fig');
handles(4) = gcf;

xdatas = []; ydatas = []; zdatas = [];
for i=1:4
    axesObjs = get(handles(i), 'Children');  %axes handles
    dataObjs = get(axesObjs(7), 'Children'); %handles to low-level graphics objects
    objType = get(dataObjs, 'Type');  %type of low-level graphics object
    xdatas{i} = get(dataObjs(3), 'XData');  %data from low-level grahics objects
    ydatas{i} = get(dataObjs(3), 'YData');
    zdatas{i} = get(dataObjs(3), 'ZData');
    %ldatas{i} = get(dataObjs,'Ldata');
end

%scatter3(logkfa2dsTrun,pTildeTrun,kappaTildeTrun);

x = []; y = []; z = []; xlin = []; ylin = []; Z = [];
%x = reshape(xdatas, 1, length(xdatas(:,1))*length(xdatas(1,:)));
%y = reshape(ydatas, 1, length(ydatas(:,1))*length(ydatas(1,:)));
%z = reshape(zdatas, 1, length(zdatas(:,1))*length(zdatas(1,:)));

%All data in one array:
x = [xdatas{1} xdatas{2} xdatas{3} xdatas{4}];
y = [ydatas{1} ydatas{2} ydatas{3} ydatas{4}];
z = [zdatas{1} zdatas{2} zdatas{3} zdatas{4}];  

if(0)
for i=1:4
    hold on;
    scatter3(xdatas{i},ydatas{i},zdatas{i});
end
end

xlin = linspace(min(x),max(x),80);
ylin = linspace(min(y),max(y),80);
        

[X,Y] = meshgrid(xlin,ylin);

%Z = griddata(x,y,z,X,Y,'cubic');
Z = griddata(x,y,z,X,Y,'cubic'); %v4 cubic linear


f= [1990 1998 2001 2004];
val = 2000; %value to find
tmp = abs(f-val);
[ival idx] = min(tmp); %index of closest value
closest = f(idx); %closest value

%To crop out along the maximal line:
%First convert the x-y surfaces to radial vectors from (0,0):


%Then find the closest vector to the point under consideration:


%Re-generate the surface, ignoring the point if the vector is greater:



close all;
figure(88);
hold on;
for i=1:4
    plot3(xdatas{i},ydatas{i},zdatas{i},'.');
    %plot3(xdatas{i},ydatas{i},zdatas{i});
end
surf(X,Y,Z,'LineStyle','none','FaceColor','interp','FaceAlpha',0.8);
hold off;
axis([0 5 0 10 0 2]);
grid on;

%Generate fancy surface subplots:
close all;
fa = figure(100);

subplot(2,2,1);
hold on;
for i=1:4
    plot3(xdatas{i},ydatas{i},zdatas{i},'.');
    %plot3(xdatas{i},ydatas{i},zdatas{i});
end
%surf(X,Y,Z,'LineStyle','none','FaceColor','interp','FaceAlpha',0.6);
axis([0 5 0 10 0 2]); grid on;
xlabel('log(kf*a2d)') % x-axis label
ylabel('pTilde') % y-axis label
zlabel('kappaTilde') % y-axis label
view(3);
%hmarkers = hs.MarkerHandle;
%hmarkers.EdgeColorData = uint8(255*[0;0;1;0.3]);

subplot(2,2,2);
hold on;
for i=1:4
    plot3(xdatas{i},ydatas{i},zdatas{i},'.');
    %plot3(xdatas{i},ydatas{i},zdatas{i});
end
%surf(X,Y,Z,'LineStyle','none','FaceColor','interp','FaceAlpha',0.6);
axis([0 5 0 10 0 2]); grid on;
xlabel('log(kf*a2d)') % x-axis label
ylabel('pTilde') % y-axis label
zlabel('kappaTilde') % y-axis label
view([1 0 0]);

subplot(2,2,3);
hold on;
for i=1:4
    plot3(xdatas{i},ydatas{i},zdatas{i},'.');
    %plot3(xdatas{i},ydatas{i},zdatas{i});
end
%surf(X,Y,Z,'LineStyle','none','FaceColor','interp','FaceAlpha',0.6);
axis([0 5 0 10 0 2]); grid on;
xlabel('log(kf*a2d)') % x-axis label
ylabel('pTilde') % y-axis label
zlabel('kappaTilde') % y-axis label
view([0 1 0]);

subplot(2,2,4);
hold on;
for i=1:4
    plot3(xdatas{i},ydatas{i},zdatas{i},'.');
    %plot3(xdatas{i},ydatas{i},zdatas{i});
end
%surf(X,Y,Z,'LineStyle','none','FaceColor','interp','FaceAlpha',0.6);
axis([0 5 0 10 0 2]); grid on;
xlabel('log(kf*a2d)') % x-axis label
ylabel('pTilde') % y-axis label
zlabel('kappaTilde') % y-axis label
view([0 0 1]);
hold off;


ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
titlestring = ['Kdata_scatter_indvClouds smoothed'];
text(0.5, 1,titlestring,'HorizontalAlignment','center','VerticalAlignment', 'top');
set(fa, 'Position', [100, 100, 1000, 800]);

figname = ['Kdata_scatter_indvClouds'];
figdirectory = 'C:\Users\tpeppler\Dropbox\PhD\2D_2016\EOS_Data\';
saveas(fa,[figdirectory figname '.fig'],'fig');
saveas(fa,[figdirectory figname '.png'],'png');






















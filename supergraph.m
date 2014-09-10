%Supergraph script:

handles = [];

if(0)
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\WidthvsAtom\CrossoverTopCam_850G_WidthvsAtom_140907.fig');
handles(1) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\WidthvsAtom\CrossoverTopCam_850G_WidthvsAtom_140905.fig');
handles(2) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\WidthvsAtom\CrossoverTopCam_880G_WidthvsAtom_140902.fig');
handles(3) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\WidthvsAtom\CrossoverTopCam_972G_WidthvsAtom_140901.fig');
handles(4) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\WidthvsAtom\CrossoverTopCam_920G_WidthvsAtom_140826.fig');
handles(5) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\WidthvsAtom\CrossoverTopCam_880G_WidthvsAtom_140819.fig');
handles(6) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\WidthvsAtom\CrossoverTopCam_972G_WidthvsAtom_140819.fig');
handles(7) = gcf;
end
if(1)
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\TonTFs\CrossoverTopCam_850G_TonTFvsAtom_140907.fig');
handles(1) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\TonTFs\CrossoverTopCam_850G_TonTFvsAtom_140905.fig');
handles(2) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\TonTFs\CrossoverTopCam_880G_TonTFvsAtom_140902.fig');
handles(3) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\TonTFs\CrossoverTopCam_972G_TonTFvsAtom_140901.fig');
handles(4) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\TonTFs\CrossoverTopCam_920G_TonTFvsAtom_140826.fig');
handles(5) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\TonTFs\CrossoverTopCam_880G_TonTFvsAtom_140819.fig');
handles(6) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover GraphDump 140908\All Figs\TonTFs\CrossoverTopCam_972G_TonTFvsAtom_140819.fig');
handles(7) = gcf;
end

xdatas = []; ydatas = [];

for i=1:7
axesObjs = get(handles(i), 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects
objType = get(dataObjs, 'Type');  %type of low-level graphics object
xdatas{i} = get(dataObjs, 'XData');  %data from low-level grahics objects
ydatas{i} = get(dataObjs, 'YData') + (i-1)*0.2;
%zdata = get(dataObjs, 'ZData');
end

figure(100);
hold on;
for j=1:7
    i=8-j;
    colourShift1 = i/7;
    %colourShift2 = (1-i)/7;
    plot(xdatas{i},ydatas{i},'MarkerSize',3,...
    'MarkerFaceColor',[colourShift1 colourShift1 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
end
hold off;
grid on;













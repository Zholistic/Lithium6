%Supergraph script:
close all;
handles = [];
imageNumber = 10;

magVector = [820 832.2 855 860 865 880 900 920 950 972];

if(0)
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\All Tight Binned Indv fig\CrossoverSideCam_855G_TightWidthvsAtom_140912_BinnedIndividualFits_withRaw.fig');
handles(1) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\All Tight Binned Indv fig\CrossoverSideCam_860G_TightWidthvsAtom_140912_BinnedIndividualFits_withRaw.fig');
handles(2) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\All Tight Binned Indv fig\CrossoverSideCam_865G_TightWidthvsAtom_140912_BinnedIndividualFits_withRaw.fig');
handles(3) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\All Tight Binned Indv fig\CrossoverSideCam_880G_TightWidthvsAtom_140911_BinnedIndividualFits_withRaw.fig');
handles(4) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\All Tight Binned Indv fig\CrossoverSideCam_900G_TightWidthvsAtom_140915_BinnedIndividualFits_withRaw.fig');
handles(5) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\All Tight Binned Indv fig\CrossoverSideCam_920G_TightWidthvsAtom_140914_BinnedIndividualFits_withRaw.fig');
handles(6) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\All Tight Binned Indv fig\CrossoverSideCam_972G_TightWidthvsAtom_140914_BinnedIndividualFits_withRaw.fig');
handles(7) = gcf;
end
if(1)
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_820G_Tight_50Bins.fig');
handles(1) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140916_sidecam_832p2G_Tight_50Bins.fig');
handles(2) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_855G_Tight_50Bins.fig');
handles(3) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_860G_Tight_50Bins.fig');
handles(4) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_865G_Tight_50Bins.fig');
handles(5) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140911_sidecam_880G_Tight_50Bins.fig');
handles(6) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_900G_Tight_50Bins.fig');
handles(7) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_920G_Tight_50Bins.fig');
handles(8) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_950G_Tight_50Bins.fig');
handles(9) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_972G_Tight_50Bins.fig');
handles(10) = gcf;
end

transverseImages = 3;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140909_sidecam_TransverseWidth_5000_PixelNumber.fig');
handlesTransverse(1) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140907_sidecam_TransverseWidth_10000_PixelNumber.fig');
handlesTransverse(2) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_TransverseWidth_13500_PixelNumber.fig');
handlesTransverse(3) = gcf;

xdatas = []; ydatas = [];
xdatasT = []; ydatasT = [];
for i=1:transverseImages
    axesObjsT = get(handlesTransverse(i), 'Children');  %axes handles
    dataObjsT = get(axesObjsT, 'Children'); %handles to low-level graphics objects
    objTypeT = get(dataObjsT, 'Type');  %type of low-level graphics object
    xdatasT{i} = get(dataObjsT, 'XData');  %data from low-level grahics objects
    ydatasT{i} = get(dataObjsT, 'YData');
end

transverseMagFieldIndexs = [13 14 16 17 18 19 21 23 26 28];

for i=1:imageNumber
    axesObjs = get(handles(i), 'Children');  %axes handles
    dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects
    objType = get(dataObjs, 'Type');  %type of low-level graphics object
    xdatas{i} = get(dataObjs, 'XData');  %data from low-level grahics objects
    ydatas{i} = get(dataObjs, 'YData');
    %zdata = get(dataObjs, 'ZData');
end

%Correct for TOF issue:
ydatas{4}{2} = ydatas{4}{2}.*1.0075;



figure(100);
hold on;
for j=1:imageNumber
    xToPlot = xdatas{j}{2};
    yToPlot = ydatas{j}{2}  + (j-1)*0.5;
    %i=8-j;
    colourShift1 = j/imageNumber;
    %colourShift2 = (1-i)/7;
    plot(xToPlot,yToPlot,'MarkerSize',3,...
    'MarkerFaceColor',[colourShift1 colourShift1 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
end
hold off;
grid on;

%zMat = [];
%zMatsC = [];
%for i=1:imageNumber
%    zMat = [];
%    %zdatas{i} = [ydatas{i}{2}, magVector]
%    for j=1:length(ydatas{i}{2}(:))
%        zMat(:,j) = [ydatas{i}{2}(j), magVector(i)];
%    end
%    zMatsC{i} = zMat;
%end

%longest ydata array?
longest = 0;
for i=1:imageNumber
    if( length(ydatas{i}{2}(:)) > longest)
        longest = length(ydatas{i}{2}(:));
    end
end

padAmount = 0;
zMatrix = [];
for i=1:imageNumber
    padAmount = longest - length(ydatas{i}{2}(:));
    zMatrix(:,i) = padarray(ydatas{i}{2}(:),padAmount,NaN,'post');
end

%For cftool:
prof820 = ydatas{1}{2}(:);
prof832 = ydatas{2}{2}(:);
prof855 = ydatas{3}{2}(:);
prof860 = ydatas{4}{2}(:);
prof865 = ydatas{5}{2}(:);
prof880 = ydatas{6}{2}(:);
prof900 = ydatas{7}{2}(:);
prof920 = ydatas{8}{2}(:);
prof950 = ydatas{9}{2}(:);
prof972 = ydatas{10}{2}(:);
profx = xdatas{2}{2};

elbowsN2D = [profx(1) profx(3) profx(8) profx(9) profx(10) profx(11) profx(14) profx(14) profx(18) profx(16)];
widthElbowsN2D = [prof820(1) prof832(3) prof855(8) prof860(9) prof865(10) prof880(11) prof900(14) prof920(14) prof950(18) prof972(16)];

%N2D scaling function:
N2D = @(x)(0.0484.*x.^2 - 130.72.*x + 116982);



figure(200);
hold on;

yToPlot = xdatas{2}{2};
xToPlot = magVector;
zToPlot = zMatrix;

surfc(xToPlot,yToPlot,zToPlot);
%shading interp;

hold off;
grid on;

fields = [];
for j=1:imageNumber
   for i=1:length(ydatas{j}{2})
   fields(i,j) = magVector(j);
   end 

end
%%%%%Fancy Surface Plot

figure(201);
hold on;

yToPlot = xdatas{4}{2};
xToPlot = magVector;
zToPlot = zMatrix;

yLocT1 = [5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000 5000];
yLocT2 = [10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000 10000];
yLocT3 = [13500 13500 13500 13500 13500 13500 13500 13500 13500 13500 13500 13500 13500 13500 13500 13500];
plot3(xdatasT{1}(13:end),yLocT1.*1.2,ydatasT{1}(13:end),'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
    'Marker','o',...
    'LineWidth',1,...
    'LineStyle',':');
plot3(xdatasT{2}(13:end),yLocT2.*1.2,ydatasT{2}(13:end),'MarkerFaceColor',[0 0.5 0],...
    'MarkerEdgeColor',[0.75 0.87 0.77],...
    'Marker','square',...
    'LineWidth',1,...
    'LineStyle','--',...
    'Color',[0 0.5 0]);
plot3(xdatasT{3}(13:end),yLocT3.*1.2,ydatasT{3}(13:end),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],...
    'Marker','diamond',...
    'LineWidth',1,...
    'LineStyle','--',...
    'Color',[1 0 0]);
for i=1:imageNumber
    plot3(fields(:,i),xdatas{i}{2},ydatas{i}{2},'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],...
    'Marker','.',...
    'LineWidth',1,...
    'LineStyle','--',...
    'Color',[0.5 0.5 0.5]);   
end

%Elbows plot:
plot3(magVector,elbowsN2D,widthElbowsN2D,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0],...
    'Marker','o',...
    'LineWidth',1.5,...
    'Color',[0 1 0]);

%N2D function plot:
plot3(magVector(1):magVector(end),N2D(magVector(1):magVector(end)),ones([length(magVector(1):magVector(end)) 1]).*7);
    
surfc(xToPlot,yToPlot,inpaint_nans(zToPlot),'LineStyle','none','FaceColor','interp','FaceAlpha',0.8);
%shading interp;
%surfc(xToPlot,yToPlot,zToPlot);
hold off;
grid on;
















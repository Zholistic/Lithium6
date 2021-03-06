%Supergraph script:
close all;
handles = [];
imageNumber = 10;
bins = 36;

massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
%pixel Lengths: 2.84 um topcam, topcam magnification = 4.58. 3.75 um
%sidecam, magnification 1.4 (? 1.33)
pixelLength = 3.75e-6 / 1.9; %2.84 um topcam, 1.9 magnification (old 1.4 mag)


magVector = [820 832.2 855 860 865 880 900 920 950 972];

if(1)
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

    if(0)
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_820G_Tight_50Bins.fig');
handles(1) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140916_sidecam_832p2G_Tight_50Bins.fig');
handles(2) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_855G_Tight_50Bins.fig');
handles(3) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_860G_Tight_50Bins.fig');
handles(4) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_865G_Tight_50Bins.fig');
handles(5) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140911_sidecam_880G_Tight_50Bins.fig');
handles(6) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_900G_Tight_50Bins.fig');
handles(7) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_920G_Tight_50Bins.fig');
handles(8) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_950G_Tight_50Bins.fig');
handles(9) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_972G_Tight_50Bins.fig');
handles(10) = gcf;
    end

if(0)
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_820G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(1) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140916_sidecam_832p2G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(2) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_855G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(3) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_860G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(4) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_865G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(5) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140911_sidecam_880G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(6) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_900G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(7) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_920G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(8) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_950G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(9) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_972G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(10) = gcf;

open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_820G_Radial_Log1050Bins.fig');
handlesLog(1) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140916_sidecam_832p2G_Radial_Log1050Bins.fig');
handlesLog(2) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_855G_Radial_Log1050Bins.fig');
handlesLog(3) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_860G_Radial_Log1050Bins.fig');
handlesLog(4) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_865G_Radial_Log1050Bins.fig');
handlesLog(5) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140911_sidecam_880G_Radial_Log1050Bins.fig');
handlesLog(6) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_900G_Radial_Log1050Bins.fig');
handlesLog(7) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_920G_Radial_Log1050Bins.fig');
handlesLog(8) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_950G_Radial_Log1050Bins.fig');
handlesLog(9) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_972G_Radial_Log1050Bins.fig');
handlesLog(10) = gcf;

open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_820G_Radial_50Bins.fig');
handlesRadial(1) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140916_sidecam_832p2G_Radial_50Bins.fig');
handlesRadial(2) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_855G_Radial_50Bins.fig');
handlesRadial(3) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_860G_Radial_50Bins.fig');
handlesRadial(4) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_865G_Radial_50Bins.fig');
handlesRadial(5) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140911_sidecam_880G_Radial_50Bins.fig');
handlesRadial(6) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_900G_Radial_50Bins.fig');
handlesRadial(7) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_920G_Radial_50Bins.fig');
handlesRadial(8) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_950G_Radial_50Bins.fig');
handlesRadial(9) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_972G_Radial_50Bins.fig');
handlesRadial(10) = gcf;

open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_820G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(1) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140916_sidecam_832p2G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(2) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_855G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(3) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_860G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(4) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_865G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(5) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140911_sidecam_880G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(6) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_900G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(7) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_920G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(8) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_950G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(9) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_972G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(10) = gcf;

transverseImages = 3;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140909_sidecam_TransverseWidth_5000_PixelNumber.fig');
handlesTransverse(1) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140907_sidecam_TransverseWidth_10000_PixelNumber.fig');
handlesTransverse(2) = gcf;
open('C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_TransverseWidth_13500_PixelNumber.fig');
handlesTransverse(3) = gcf;
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


%Radial is Tight in this dataset!!!!!!!! Whoops
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_820G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(1) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140916_sidecam_832p2G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(2) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_855G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(3) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_860G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(4) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_865G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(5) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140911_sidecam_880G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(6) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_900G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(7) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_920G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(8) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_950G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(9) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_972G_Tight_36Bins_BinAvg.fig');
handlesBin36CA(10) = gcf;

open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_820G_Radial_Log1050Bins.fig');
handlesLog(1) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140916_sidecam_832p2G_Radial_Log1050Bins.fig');
handlesLog(2) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_855G_Radial_Log1050Bins.fig');
handlesLog(3) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_860G_Radial_Log1050Bins.fig');
handlesLog(4) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_865G_Radial_Log1050Bins.fig');
handlesLog(5) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140911_sidecam_880G_Radial_Log1050Bins.fig');
handlesLog(6) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_900G_Radial_Log1050Bins.fig');
handlesLog(7) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_920G_Radial_Log1050Bins.fig');
handlesLog(8) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_950G_Radial_Log1050Bins.fig');
handlesLog(9) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_972G_Radial_Log1050Bins.fig');
handlesLog(10) = gcf;

open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_820G_Radial_50Bins.fig');
handlesRadial(1) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140916_sidecam_832p2G_Radial_50Bins.fig');
handlesRadial(2) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_855G_Radial_50Bins.fig');
handlesRadial(3) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_860G_Radial_50Bins.fig');
handlesRadial(4) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_865G_Radial_50Bins.fig');
handlesRadial(5) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140911_sidecam_880G_Radial_50Bins.fig');
handlesRadial(6) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_900G_Radial_50Bins.fig');
handlesRadial(7) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_920G_Radial_50Bins.fig');
handlesRadial(8) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_950G_Radial_50Bins.fig');
handlesRadial(9) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_972G_Radial_50Bins.fig');
handlesRadial(10) = gcf;

open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_820G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(1) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140916_sidecam_832p2G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(2) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_855G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(3) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_860G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(4) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140912_sidecam_865G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(5) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140911_sidecam_880G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(6) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_900G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(7) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_920G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(8) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_950G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(9) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140914_sidecam_972G_Radial_36Bins_BinAvg.fig');
handlesRadialBin36CA(10) = gcf;

transverseImages = 5;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140909_sidecam_TransverseWidth_5000_PixelNumber.fig');
handlesTransverse(1) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140907_sidecam_TransverseWidth_10000_PixelNumber.fig');
handlesTransverse(2) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140915_sidecam_TransverseWidth_13500_PixelNumber.fig');
handlesTransverse(3) = gcf;

open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140909_sidecam_TransverseWidth_5000_PixelNumber_CenterAvg.fig');
handlesTransverse(4) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Newbinned\140907_sidecam_TransverseWidth_10000_PixelNumber_CenterAvg.fig');
handlesTransverse(5) = gcf;

xdatas = []; ydatas = [];
xdatasT = []; ydatasT = [];
xdatasLog10 = []; ydatasLog10 = [];
xdatasR = []; ydatasR = [];

for i=1:transverseImages
    axesObjsT = get(handlesTransverse(i), 'Children');  %axes handles
    dataObjsT = get(axesObjsT, 'Children'); %handles to low-level graphics objects
    objTypeT = get(dataObjsT, 'Type');  %type of low-level graphics object
    xdatasT{i} = get(dataObjsT, 'XData');  %data from low-level grahics objects
    ydatasT{i} = get(dataObjsT, 'YData');
    ldatasT{i} = get(dataObjsT,'Ldata');
end

transverseMagFieldIndexs = [13 14 16 17 18 19 21 23 26 28];

for i=1:imageNumber
    axesObjs = get(handlesBin36CA(i), 'Children');  %axes handles
    %axesObjs = get(handles(i), 'Children');  %axes handles Original 50 bin 
    dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects
    objType = get(dataObjs, 'Type');  %type of low-level graphics object
    xdatas{i} = get(dataObjs, 'XData');  %data from low-level grahics objects
    ydatas{i} = get(dataObjs, 'YData');
    ldatas{i} = get(dataObjs(3),'Ldata');
    %zdata = get(dataObjs, 'ZData');
end

for i=1:imageNumber
    axesObjsR = get(handlesRadialBin36CA(i), 'Children');  %axes handles
    dataObjsR = get(axesObjsR, 'Children'); %handles to low-level graphics objects
    objTypeR = get(dataObjsR, 'Type');  %type of low-level graphics object
    xdatasR{i} = get(dataObjsR, 'XData');  %data from low-level grahics objects
    ydatasR{i} = get(dataObjsR, 'YData');
    %zdata = get(dataObjs, 'ZData');
end

for i=1:imageNumber
    axesObjsL = get(handlesLog(i), 'Children');  %axes handles
    dataObjsL = get(axesObjsL, 'Children'); %handles to low-level graphics objects
    objTypeL = get(dataObjsL, 'Type');  %type of low-level graphics object
    xdatasLog10{i} = get(dataObjsL, 'XData');  %data from low-level grahics objects
    ydatasLog10{i} = get(dataObjsL, 'YData');
    %zdata = get(dataObjs, 'ZData');
end

%Correct for TOF issue:
ydatas{4}{3} = ydatas{4}{3}.*1.0075;

figure(100);
hold on;
for j=1:imageNumber
    xToPlot = xdatas{j}{3};
    yToPlot = ydatas{j}{3}  + (j-1)*0.5;
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

figure(101);
hold on;
for j=1:imageNumber
    xToPlotL = xdatasLog10{j};
    yToPlotL = ydatasLog10{j}  + (j-1)*0.05;
    %i=8-j;
    colourShift1 = j/imageNumber;
    %colourShift2 = (1-i)/7;
    plot(xToPlotL,yToPlotL,'MarkerSize',3,...
    'MarkerFaceColor',[colourShift1 colourShift1 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
end
hold off;
grid on;

figure(102);
hold on;
for j=1:imageNumber
    xToPlotR = xdatasR{j};
    yToPlotR = ydatasR{j}  + (j-1)*5;
    %i=8-j;
    colourShift1 = j/imageNumber;
    %colourShift2 = (1-i)/7;
    plot(xToPlotR,yToPlotR,'MarkerSize',3,...
    'MarkerFaceColor',[colourShift1 colourShift1 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
end
hold off;
grid on;

%Aspect ratio plots
figure(300);
hold on;
for j=1:imageNumber
    %xToPlotAR = xdatasR{j}{2};
    %yToPlotAR = ydatasR{j}{2}./ydatas{j}{2}  + (j-1)*1;
    xToPlotAR = xdatasR{j};
    yToPlotAR = ydatasR{j}./ydatas{j}{3}  + (j-1)*1;
    %i=8-j;
    colourShift1 = j/imageNumber;
    %colourShift2 = (1-i)/7;
    plot(xToPlotAR,yToPlotAR,'MarkerSize',3,...
    'MarkerFaceColor',[colourShift1 colourShift1 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
end
hold off;
grid on;

%2015 Elbow work:
peakPoly = [36.5674, 53.2957, 45.5677, 54.4748, 54.6825, 62.1721, 53.5184, 64.859];
atomPoly = [7566, 12468, 10813, 15822, 18866, 18917, 20692, 18746];

figure(301);
hold on; coefsLF = []; coefsLFs = [];
for j=1:imageNumber
if(0)
    xToPlotAR = []; yToPlotAR = [];
    %xToPlotAR = xdatasR{j}{2};
    %yToPlotAR = ydatas{j}{2}./ydatasR{j}{2}  + (j-1)*0.05;
    xToPlotARt = xdatasR{j};
    yToPlotARt = ydatas{j}{3}./ydatasR{j}  + (j-1)*0.05;
    
    yToPlotAR = log10(yToPlotARt);
    xToPlotAR = log10(xToPlotARt);
    %i=8-j;
    colourShift1 = j/imageNumber;
    %colourShift2 = (1-i)/7;    
    
    %Line fit to > 30000:
    fgl = @(p,x)(p(1).*x + p(2)); 
    p0 = [0 0.3];
    lb = [-.1 0.15];
    ub = [.1 0.6];
    if(1)
        fgl = @(p,x)(0.*x + p(1));
        p0 = -0.3;
        lb = -1;
        ub = -0.1;
    end
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    %xs = 1:length(xToPlotAR);
    [coefsLF(:,j),resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fgl,p0,xToPlotAR(end-20:end),inpaint_nans(yToPlotAR(end-20:end)),lb,ub,curvefitoptions);

    
    %Line fit to < 10000:
    fgls = @(p,x)(p(1).*x + p(2)); 
    p0s = [-.1 -0.3];
    lbs = [-.3 -1];
    ubs = [0 0];
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    %xs = 1:length(xToPlotAR);
    [coefsLFs(:,j),resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fgls,p0s,xToPlotAR(1:6),inpaint_nans(yToPlotAR(1:6)),lbs,ubs,curvefitoptions);
       
    
    plot(xToPlotAR,inpaint_nans(yToPlotAR(:)),'MarkerSize',3,...
    'MarkerFaceColor',[colourShift1 colourShift1 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
    plot(xToPlotAR,fgl(coefsLF(:,j),xToPlotAR),'r');
    plot(xToPlotAR,fgls(coefsLFs(:,j),xToPlotAR),'g');
    end
end
hold off;
grid on;

fitsInitRange = [];
fitsInitRange{10} = 3:9;
fitsInitRange{9} = 3:10;
fitsInitRange{8} = 2:9;
fitsInitRange{7} = 2:8;
fitsInitRange{6} = 3:10;
fitsInitRange{5} = 1:5;
fitsInitRange{4} = 1:7;
fitsInitRange{3} = 1:4;
fitsInitRange{2} = 2:3;
fitsInitRange{1} = 1:2;

elem = 37;
fitsFinalRange = [];
fitsFinalRange{10} = [12:26 28:32];
fitsFinalRange{9} = 13:36;
fitsFinalRange{8} = [20 22:32 34:35];
fitsFinalRange{7} = [11:33 35:36];
fitsFinalRange{6} = [16:24 26:28 30:31 33:35];
fitsFinalRange{5} = [8:22];
fitsFinalRange{4} = [10 12 14 16 17 19 22 25 26 29 31];
fitsFinalRange{3} = [11 12 13 14 16 18 21 23 26 28 30 32];
fitsFinalRange{2} = [2:22 24:28 30:31 33:36];
fitsFinalRange{1} = [1:23 25 27:29 31:35];

%figure(301);
home = 0;
coefsLF = []; coefsLFs = [];
for j=1:imageNumber
    
    xToPlot = []; yToPlot = [];
    %xToPlotAR = xdatasR{j}{2};
    %yToPlotAR = ydatas{j}{2}./ydatasR{j}{2}  + (j-1)*0.05;
    %xToPlot = inpaint_nans(xdatas{j}{3});
    %yToPlot = inpaint_nans(ydatas{j}{3});
    xToPlot = xdatas{j}{3};
    yToPlot = ydatas{j}{3};
    errorY = ldatas{j};
    errorX = ydatas{j}{1};
    
    %Line fit to first points:
    %fgl = @(p,x)(p(1).*x + p(2));
    %p0 = [0 0.3];
    %lb = [-.1 0.15];
    %ub = [.1 0.6];
    if(1)
        fgl = @(p,x)(p(1) + p(2).*x);
        p0 = [4.5 0];
        lb = [3 -0.1];
        ub = [6 0.1];
    end
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    %xs = 1:length(xToPlotAR);
    [coefsLF(:,j),resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fgl,p0,xToPlot(fitsInitRange{j}),yToPlot(fitsInitRange{j}),lb,ub,curvefitoptions);
    
    
    %Line fit to latter part of data:
    fgls = @(p,x)(p(1).*x + p(2));
    p0s = [3/50000 5];
    lbs = [0 3];
    ubs = [1 6];
    curvefitoptions = optimset('MaxFunEvals',100000,'MaxIter',50000,'Display','off');
    %xs = 1:length(xToPlotAR);
    [coefsLFs(:,j),resnorm,r,exitflag,output,lambda,J] = lsqcurvefit(fgls,p0s,xToPlot(fitsFinalRange{j}),yToPlot(fitsFinalRange{j}),lbs,ubs,curvefitoptions);
    
    
    %figure(j)
    hold on;
    %plot(xToPlot,yToPlot,'MarkerSize',3,...
    %'MarkerFaceColor','b',...
    %'Marker','o',...
    %'LineStyle','--',...
    %'Color',[0 0 1]);
    %h = plotErrorLines(xToPlot,yToPlot,errorX,errorY,1,j, [0 0 1],'o');
    plot(2500:25000,fgl(coefsLF(:,j),2500:25000),'r');
    range = 15000:60000;
    if(j == 6)
        range = 12000:60000;
    end
    if(j == 5 || j == 4 || j == 3)
        range = 5000:60000;
    end
    if(j < 3)
        range = 1000:60000;
    end
    plot(range,fgls(coefsLFs(:,j),range),'g');
    line([xToPlot(fitsFinalRange{j}(1)) xToPlot(fitsFinalRange{j}(1))],[4.4 5.6], ...
        'Color',[0.7 0.7 0.7],'LineStyle','--')
    
    grid on;
    figname = [ 'ElbowGraph_Field_' num2str(magVector(j)) 'G'];
    if(home)
    figdirectory = 'C:\Users\Zholistic\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Elbows\';  
    else
    figdirectory = 'C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Elbows\';
    end
    %saveas(h,[figdirectory figname '.fig'],'fig');
    %saveas(h,[figdirectory figname '.png'],'png');
    hold off;
end


elbows = [];
elbows(:,10) = [18460 2000 0.03];
elbows(:,9) = [20000 2000 0.03];
elbows(:,8) = [19330 2000 0.03];
elbows(:,7) = [19500 3000 0.03];
elbows(:,6) = [16250 2000 0.06];
elbows(:,5) = [11530 1000 0.03];
elbows(:,4) = [12450 2500 0.05];
elbows(:,3) = [8000 500 0];
elbows(:,2) = [0 0 0];
elbows(:,1) = [0 0 0];


figure(102);
hold on;
for j=1:imageNumber
    %xToPlotR = xdatasR{j}{2};
    %yToPlotR = ydatasR{j}{2}  + (j-1)*0.5;
    xToPlotR = xdatasR{j};
    yToPlotR = ydatasR{j}  + (j-1)*0.5;
    %i=8-j;
    colourShift1 = j/imageNumber;
    %colourShift2 = (1-i)/7;
    plot(xToPlotR,yToPlotR,'MarkerSize',3,...
    'MarkerFaceColor',[colourShift1 colourShift1 1],...
    'Marker','o',...
    'LineStyle','--',...
    'Color',[0 0 1]);
end
hold off;
grid on;


%figure(105);
hold on;
k=1;
for j=[2 5 9]
    k = k+1;
    xToPlot = xdatas{j}{3};
    %yToPlot = ydatas{j}{2}  + (3-(k-1))*0.5;
    yToPlot = ydatas{j}{3};
    errorY = ldatas{j};
    errorX = ydatas{j}{1};
    %i=8-j;
    colourShift1 = j/imageNumber;
    %colourShift2 = (1-i)/7;
    if(0)
    %plot(xToPlot,yToPlot,'MarkerSize',3,...
    %'MarkerFaceColor',[colourShift1 colourShift1 1],...
    %'Marker','o',...
    %'LineStyle','--',...
    %'Color',[0 0 1]);
    end
    if( j == 2 )
        plotErrorLines(xToPlot,yToPlot + 0.3,errorX,errorY,1,1000, [1 0 0],'^');
    elseif( j == 5 )
        plotErrorLines(xToPlot,yToPlot,errorX,errorY,1,1000, [0 0.5 0],'square');
    else
        plotErrorLines(xToPlot,yToPlot - 0.7,errorX,errorY,1,1000, [0 0 1],'o');
    end
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
    if( length(ydatas{i}{3}(:)) > longest)
        longest = length(ydatas{i}{3}(:));
    end
end

padAmount = 0;
zMatrix = []; yMatrix = []; xMatrix = [];

for i=1:imageNumber
    %padAmount = longest - length(ydatas{i}{3}(:));
    zMatrix(:,i) = padarray(ydatas{i}{3}(:),padAmount,NaN,'post');
end

padAmount = 0;
zMatrixL = [];
for i=1:imageNumber
    padAmount = 0; %all same length...
    zMatrixL(:,i) = padarray(ydatasLog10{i}(:),padAmount,NaN,'post');
end

padAmount = 0;
zMatrixR = [];
for i=1:imageNumber
    padAmount = 0; %all same length...
    %zMatrixR(:,i) = padarray(ydatasR{i}{2}(:),padAmount,NaN,'post');
    zMatrixR(:,i) = padarray(ydatasR{i}(:),padAmount,NaN,'post');
end

%For cftool:
if(bins == 50)
end
prof820 = ydatas{1}{3}(:);
prof832 = ydatas{2}{3}(:);
prof855 = ydatas{3}{3}(:);
prof860 = ydatas{4}{3}(:);
prof865 = ydatas{5}{3}(:);
prof880 = ydatas{6}{3}(:);
prof900 = ydatas{7}{3}(:);
prof920 = ydatas{8}{3}(:);
prof950 = ydatas{9}{3}(:);
prof972 = ydatas{10}{3}(:);
profx = xdatas{2}{3};

elbowsN2D = [profx(1) profx(3) profx(8) profx(9) profx(10) profx(11) profx(14) profx(14) profx(18) profx(16)];
widthElbowsN2D = [prof820(1) prof832(3) prof855(8) prof860(9) prof865(10) prof880(11) prof900(14) prof920(14) prof950(18) prof972(16)];

%N2D scaling function:
N2D = @(x)(0.0512.*x.^2 - 138.49.*x + 124322);
%0.0512x2 - 138.49x + 124322


figure(200);
hold on;

yToPlot = xdatas{2}{3};
xToPlot = magVector;
zToPlot = zMatrix;

surfc(xToPlot,yToPlot,zToPlot);
%shading interp;

hold off;
grid on;

fields = [];
for j=1:imageNumber
   for i=1:length(ydatas{j}{3})
   fields(i,j) = magVector(j);
   end 

end


%%%%%Figure 2:


%plotErrorLines();
%Figure 2 a) Main: (12k atoms)
yToPlotTransverse = ydatasT{5};
xToPlotTransverse = xdatasT{5};
errorTransverse = ldatasT{5};
plotErrorLinesFig2aMain(xToPlotTransverse,yToPlotTransverse.*2.*pixelLength./1e-6,[],errorTransverse.*2.*pixelLength./1e-6,0,1,[0 0 1],'o');
csvwrite('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\fig2a_12kAtom.csv',[xToPlotTransverse yToPlotTransverse.*2.*pixelLength./1e-6 errorTransverse.*2.*pixelLength./1e-6]);
%Figure 2 a) Inset: (5k atoms)
yToPlotTransverse = ydatasT{4};
xToPlotTransverse = xdatasT{4};
errorTransverse = ldatasT{4};
plotErrorLinesFig2aInset(xToPlotTransverse,yToPlotTransverse.*2.*pixelLength./1e-6,[],errorTransverse.*2.*pixelLength./1e-6,0,2,[1 0 0],'square');
csvwrite('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\fig2a_6kAtom.csv',[xToPlotTransverse yToPlotTransverse.*2.*pixelLength./1e-6 errorTransverse.*2.*pixelLength./1e-6]);

%Figure 2 b)    
%hold on;
close all;
k=1;
for j=[2 5 9]
    k = k+1;
    xToPlot = xdatas{j}{3};
    %yToPlot = ydatas{j}{2}  + (3-(k-1))*0.5;
    yToPlot = ydatas{j}{3}.*2.*pixelLength./1e-6;
    errorY = ldatas{j}.*2.*pixelLength./1e-6;
    errorX = ydatas{j}{1};
    %i=8-j;
    %colourShift1 = j/imageNumber;
    %colourShift2 = (1-i)/7;
    range = 15000:60000;
    if(j == 6)
        range = 12000:60000;
    end
    if(j == 5 || j == 4 || j == 3)
        range = 5000:60000;
    end
    if(j < 3)
        range = 1000:60000;
    end
    
    if( j == 2 )
        hold off;
        
        h1 = plotErrorLines(xToPlot,yToPlot,errorX,errorY,1,1000, [1 0 0],[1 0.8 0.8],'^',[]);
        %plot(range,fgls(coefsLFs(:,j),range) + 0.3,'color',[0.5 0.5 0.5]);
        %plot(2500:25000,fgl(coefsLF(:,j),2500:25000)+ 0.3,'color',[0.5 0.5 0.5]);
    elseif( j == 5 )
        hold off;
        figure(1001);
        hold on;
        plot(18000:39000,fgls(coefsLFs(:,j),18000:39000).*2.*pixelLength./1e-6,'color',[0.4 0.4 0.4],'LineStyle','-');
        plot(3000:14000,fgl(coefsLF(:,j),3000:14000).*2.*pixelLength./1e-6,'color',[0.4 0.4 0.4],'LineStyle','-');
        plot(1:39000,fgls(coefsLFs(:,j),1:39000).*2.*pixelLength./1e-6,'color',[0.3 0.3 0.3],'LineStyle','--');
        plot(1:39000,fgl(coefsLF(:,j),1:39000).*2.*pixelLength./1e-6,'color',[0.3 0.3 0.3],'LineStyle','--');
        h2 = plotErrorLines(xToPlot,yToPlot,errorX,errorY,1,1001, [0 0.65 0], [0.6 0.83 0.6],'square',fitsInitRange{j},fitsFinalRange{j});
        
    else
        hold off;
        figure(1002);
        hold on;
        plot(17000:39000,fgls(coefsLFs(:,j),17000:39000).*2.*pixelLength./1e-6,'color',[0.4 0.4 0.4],'LineStyle','-');
        plot(3500:22000,fgl(coefsLF(:,j),3500:22000).*2.*pixelLength./1e-6,'color',[0.4 0.4 0.4],'LineStyle','-');
        plot(1:39000,fgls(coefsLFs(:,j),1:39000).*2.*pixelLength./1e-6,'color',[0.3 0.3 0.3],'LineStyle','--');
        plot(1:39000,fgl(coefsLF(:,j),1:39000).*2.*pixelLength./1e-6,'color',[0.3 0.3 0.3],'LineStyle','--');
        h3 = plotErrorLines(xToPlot,yToPlot,errorX,errorY,1,1002, [0 0 1], [0.7 0.7 1],'o',fitsInitRange{j},fitsFinalRange{j});       
    end
    
end
%hLegend = legend( ...
%    [h1, h2, h3], ...
%    '832.2' , ...
%    '865'      , ...
%    '950'       , ...
%    'location', 'NorthWest' );
%hold off;
%grid on;

%Wei Zhang Theory:
wei1 = xlsread('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\q2dbose.xlsx',1);
wei2 = xlsread('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\q2dbose.xlsx',2);
wei3 = xlsread('C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\q2dbose.xlsx',3);

h = figure(2); 
width = 3; % Width in inches
height = 2;
fsz = 8;
alw = 0.75;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]);
set(gca, 'FontSize', fsz, 'LineWidth', alw,'YMinorTick','on','XMinorTick','on');
plot(wei3(:,12),wei3(:,11),'LineWidth', 2,'LineStyle','--','color',[0 0 0]); xlim([680 990]);
hold on;
plot(wei2(1:end,7),wei2(1:end,6),'r','LineWidth', 2);
plot(wei1(1:end,7),wei1(1:end,6),'b','LineWidth', 2);
 line([832.2 832.2],[0 1],'Color',[0 0 0],'LineWidth', 1,'LineStyle','-.');  
hold off;

set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);

print(['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\figure2_wei_v2'],'-depsc2','-r300');
print(['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\figure2_wei_v2'],'-dpng','-r300');
saveas(h,['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure2\figure2_wei_v2'],'fig');


%%%%%Fancy Surface Plot
%%%%%%%%%%%%%%%%%%%%%%%%

%%TODO: remove inpaint_nans here...
%need to oversubscribe the points to account for atom number discrepancies...
%yToPlot = inpaint_nans(xdatas{4}{3});
yToPlotNums = (2000:(60000-2000)/36:60000)';
atomNumsShort = 2000:(60000-2000)/10:60000;
%xToPlot = magVector;
zToPlot = zMatrix;

%Scaling for dimensionless units:
omegaR = 2*pi*24.5;
omegaZ = 2*pi*5150;
a0 = 5.29e-11;
Lz = sqrt(hbar/(massL6*omegaZ));
%a2D/a0 (a0 is unit)
a2DVector = a0.*[15535.9 21729.5 40813.7 46471 52715.6 75487.2 116050 169722 277539 378012]; %A3D!!!! CARE
omegaRVector = 2.*pi.*[24.22100769 24.40201865 24.73680871 24.80963347 24.88224858 25.09885405 25.38483326 25.66767827 26.08628949 26.3891196];
omegaRLong = 2.*pi.*[25.221:(26.3891196-25.221)/36:26.3891196];
Lzona3d = [0.331976 0.00464385 -0.545481 -0.65335 -0.75747 -1.04901 -1.39474 -1.6987 -2.09034 -2.33704];

a2DonLz = a2DVector./Lz;

%yToPlot = sqrt(yToPlot.*2).*omegaR./omegaZ;
%Reparse each point, taking in its x&y and giving it new x&y:
xdatasEfOnHbarOmegaZ = []; logKfA2D = [];
for i=1:length(yToPlotNums)
    for j=1:length(magVector)
        if ~isnan(xdatas{j}{3}(i))
            %omegaR = N2D(magVector(j));
            xdatasEfOnHbarOmegaZ{j}(i) = (omegaRVector(j)/omegaZ).*sqrt(2.*xdatas{j}{3}(i));
            logKfA2D{j}(i) = log(sqrt((2.*massL6.*omegaRVector(j).*sqrt(2.*xdatas{j}{3}(i)))./hbar).*a2DVector(j));

            %sqrt((2.*massL6.*omegaRVector(j))./hbar).*((2.*xdatas{j}{3}(i)).^(1/4)).*a2DVector(j)
        else
            xdatasEfOnHbarOmegaZ{j}(i) = NaN;
            logKfA2D{j}(i) = NaN;
        end   
    end
end

%XYZ array:
masterMatrix = [];
k = 1;
for j=1:10
    for i=1:length(logKfA2D{1})
        if(isnan(logKfA2D{j}(i)) || isnan(xdatasEfOnHbarOmegaZ{j}(i)) || isnan(ydatas{j}{3}(i)))
        else
            %masterMatrix(:,k) = [logKfA2D{j}(i) xdatasEfOnHbarOmegaZ{j}(i) ydatas{j}{3}(i)];
            %masterMatrix(:,k) = [log(a2DonLz(j)) xdatasEfOnHbarOmegaZ{j}(i) ydatas{j}{3}(i)];
            masterMatrix(:,k) = [Lzona3d(j) xdatasEfOnHbarOmegaZ{j}(i) ydatas{j}{3}(i)];
            %logKfA2D{j}(i)
            k = k+1;
        end
    end
end

ElbowsN2DEf = (omegaRVector./omegaZ).*sqrt(2.*elbows(1,:));
ElbowsN2DEfError = (omegaRVector./omegaZ).*sqrt(2.*elbows(2,:));
%paintedNansMasterMatrix = inpaint_nans(masterMatrix);

close all;
h = figure(3000);
x = []; y = []; z = []; xlin = []; ylin = []; Z = [];
x = masterMatrix(1,:);
y = masterMatrix(2,:);
z = masterMatrix(3,:).*2.*pixelLength./1e-6;

xlin = linspace(min(x),max(x),1000);
ylin = linspace(min(y),max(y),1000);
[X,Y] = meshgrid(xlin,ylin);

%Z = griddata(x,y,z,X,Y,'cubic');
Z = griddata(x,y,z,X,Y,'cubic'); %v4 cubic linear
%zi = interp2(x, y, z, X, Y, 'cubic');

surf(X,Y,Z,'LineStyle','none','FaceColor','interp','FaceAlpha',1);
%contourf(X,Y,Z,15); 
%caxis([4.5,7]);
axis tight; hold on; grid off;

%Elbow error:
for i=1:length(Lzona3d)
    plot3([Lzona3d(i) Lzona3d(i)],[(ElbowsN2DEf(i) - ElbowsN2DEfError(i)./3) (ElbowsN2DEf(i) + ElbowsN2DEfError(i)./3)],ones([2 1]).*37,'LineStyle','-','LineWidth',1,'Color',[0.6 0.6 0.6]);
end
plot3(linspace(min(Lzona3d),max(Lzona3d),40),fittedmodelElbow4(linspace(min(Lzona3d),max(Lzona3d),40)),ones([length(linspace(min(Lzona3d),max(Lzona3d),40)) 1]).*37,'LineStyle','--','LineWidth',1,'Color',[0.8 0.8 0.8]);
plot3(Lzona3d(1:end),ElbowsN2DEf,ones([length(ElbowsN2DEf) 1]).*37,'MarkerFaceColor',[0.8 0.8 0.8],'MarkerEdgeColor',[0.6 0.6 0.6],...
    'Marker','o',...
    'MarkerSize',8,...
    'LineWidth',1,...
    'LineStyle','none',...
    'Color',[1 1 1]);
plot3(x,y,z + 1,'o','MarkerSize',2,'Color',[0.4 0.4 0.4]);
plot3(linspace(min(Lzona3d),max(Lzona3d),40),ones([40 1]),ones([40 1]).*37,'LineWidth',0.8,'Color',[1 1 1]);
view(2);
hold off;
set(gca,'xdir','reverse');
box on;
colorbar;

% Save the file as PNG
print(['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure3\figure3_y_Efonhbaromegaz_x_lzona3d_v6'],'-depsc2','-r300');
print(['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure3\figure3_y_Efonhbaromegaz_x_lzona3d_v6'],'-dpng','-r300');
saveas(h,['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure3\figure3_y_Efonhbaromegaz_x_lzona3d_v6'],'fig'); %(h,[figdirectory figname '.fig'],'fig')
csvwrite(['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure3\figure3_surface_interpolated.csv'],[X Y Z]);
csvwrite(['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure3\figure3_surface_raw.csv'],masterMatrix);
csvwrite(['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure3\figure3_elbows.csv'],ElbowsN2DEf);
csvwrite(['C:\Users\tpeppler\Dropbox\PhD\2DEOSandCrossover\Crossover Sidecam Sequence\Figure3\figure3_elbows_error.csv'],ElbowsN2DEfError./3);

for i=1:length(logKfA2D)
    minL = min(logKfA2D{i});
end
for i=1:length(yToPlotNums)
yToPlot(i) = sqrt(yToPlotNums(i).*2).*omegaRLong(i)./omegaZ;
end
for j=1:length(magVector)
xToPlot(j) = log(sqrt((2.*massL6.*omegaRVector(j).*sqrt(2.*atomNumsShort(j)))./hbar).*a2DVector(j));
end
figure(2001);
hold on;
for i=1:imageNumber
    plot3(logKfA2D{i},xdatasEfOnHbarOmegaZ{i},ydatas{i}{3},'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],...
    'Marker','.',...
    'LineWidth',1,...
    'LineStyle','--',...
    'Color',[0.5 0.5 0.5]);   
end
contour3(xToPlot,yToPlot,inpaint_nans(zToPlot));
%surfc(xToPlot,yToPlot,inpaint_nans(zToPlot),'LineStyle','none','FaceColor','interp','FaceAlpha',0.8);
%surfc(inpaint_nans(zToPlot),'LineStyle','none','FaceColor','interp','FaceAlpha',0.8);
%shading interp;
%surfc(xToPlot,yToPlot,zToPlot);
hold off;
grid on;

figure(201);
hold on;

yToPlot = []; xToPlot = [];
yToPlot = (2000:(60000-2000)/36:60000)';
xToPlot = magVector;
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
    plot3(fields(:,i),xdatas{i}{3},ydatas{i}{3},'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],...
    'Marker','.',...
    'LineWidth',1,...
    'LineStyle','--',...
    'Color',[0.5 0.5 0.5]);   
end

%Elbows plot:
if(bins == 50)
plot3(magVector,elbowsN2D,widthElbowsN2D,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0],...
    'Marker','o',...
    'LineWidth',1.5,...
    'Color',[0 1 0]);
end

%N2D function plot:
plot3(magVector(1):magVector(end),N2D(magVector(1):magVector(end))./2,ones([length(magVector(1):magVector(end)) 1]).*7);
    
surfc(xToPlot,yToPlot,inpaint_nans(zToPlot),'LineStyle','none','FaceColor','interp','FaceAlpha',0.8);
%shading interp;
%surfc(xToPlot,yToPlot,zToPlot);
hold off;
grid on;

figure(500)
hold on;
plot3(magVector,elbows(1,:),ones([length(magVector) 1]).*7);
for i=1:length(elbows(1,:))
    plot3([magVector(i) magVector(i)],[(elbows(1,i) - elbows(2,i)) (elbows(1,i) + elbows(2,i))],ones([2 1]).*7);
end
hold off;
grid on;


%%%%%%Log10 surface plot:
figure(202);

yToPlotL = xdatasLog10{4};
xToPlotL = magVector;
zToPlotL = zMatrixL;
surfc(xToPlotL,yToPlotL,inpaint_nans(zToPlotL),'LineStyle','none','FaceColor','interp','FaceAlpha',0.8);

grid on;

%%%%%%Radials surface plot:
figure(203);

%yToPlotR = xdatasR{4}{2};
yToPlotR = xdatasR{4};
xToPlotR = magVector;
zToPlotR = zMatrixR;
surfc(xToPlotR,yToPlotR,inpaint_nans(zToPlotR),'LineStyle','none','FaceColor','interp','FaceAlpha',0.8);

grid on;

%%%%%%Aspect Ratio surface plot:
figure(204);
hold on;
%yToPlotAR = xdatasR{4}{2};
yToPlotRR = xdatasR{4};
xToPlotAR = magVector;
zToPlotAR = zMatrixR./zMatrix;
surfc(xToPlotAR,yToPlotAR,inpaint_nans(zToPlotAR),'LineStyle','none','FaceColor','interp','FaceAlpha',0.8);
plot3(magVector(1):magVector(end),N2D(magVector(1):magVector(end))./2,ones([length(magVector(1):magVector(end)) 1]).*7);
hold off;
grid on;


%%%%%Delaunay plot:
figure(210);
hold on;
tri = delaunay(xMatrix,inpaint_nans(yMatrix));
[r,c] = size(tri);
disp(r)
h = trisurf(tri, xMatrix, inpaint_nans(yMatrix), inpaint_nans(zMatrix));
%axis vis3d
%l = light('Position',[-50 -15 29])
%set(gca,'CameraPosition',[208 -50 7687])
%lighting phong
shading interp
%colorbar EastOutside

hold off;

%%%%%%Griddata plot:
xMatrixRVector = reshape(xMatrix,[],1);
yMatrixRVector = reshape(inpaint_nans(yMatrix),[],1);
zMatrixRVector = reshape(inpaint_nans(zMatrix),[],1);

[xqS,yqS] = meshgrid(magVector, 2000:10:60000);

vqS = griddata(xMatrixRVector,yMatrixRVector,zMatrixRVector,xqS,yqS);









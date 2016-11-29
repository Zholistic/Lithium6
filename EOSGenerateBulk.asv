function [ ] = EOSGenerateBulk(inputDensityArray, inputRadiusArray, calcRegion, omegaR, field, a2d, smoothOn, zeroOn, saveOn, savestring)

massL6 = 9.988e-27; %9.988 x 10^27 kg
hbar = 1.05457e-34; %1.05457*10^-34 m^2 kg/s
kB = 1.38e-23; %m^2 kg s^-2 K^-1
a0 = 5.29e-11; %o.0
smoothAmount = 6;

pixelLength = 2.84e-6; %13 um topcam, topcam magnification = 4.58, ie 2.84um effective
kpixelLength = (13e-6*(83/400)); %kristian pixel length

if(0)
    for i=1:length(inputDensityArray)
        if(mod(i,10) == 0)
            figure(i);
            plot(inputRadiusArray(1,calcRegion,i),inputDensityArray(1,calcRegion,i),'r'); hold off;
        end
    end
end

inputDensityArrayS = []; shiftByThis = []; shiftByThisMin = [];
if(smoothOn)
    for i=1:length(inputDensityArray(1,1,:))
        inputDensityArrayS(1,calcRegion(1)-1:calcRegion(end)+1,i) = smooth(inputDensityArray(1,calcRegion(1)-1:calcRegion(end)+1,i),smoothAmount);
    end
    
    if(zeroOn)
        %zero
        for i=1:length(inputDensityArray(1,1,:))
            shiftByThis(i) = mean(inputDensityArrayS(1,calcRegion(end)-15:calcRegion(end),i));
            %shiftByThis(i) = min(inputDensityArrayS(1,calcRegion,i));
            %inputDensityArrayS(1,:,i) = inputDensityArrayS(1,:,i) - shiftByThis(i);
            shiftByThisMin(i) = min(inputDensityArrayS(1,calcRegion,i));
            inputDensityArrayS(1,:,i) = inputDensityArrayS(1,:,i) - shiftByThisMin(i);
        end
    end
else
    inputDensityArrayS = inputDensityArray;
    
    if(zeroOn)
        for i=1:length(inputDensityArray(1,1,:))
            shiftByThis(i) = mean(inputDensityArrayS(1,calcRegion(end)-5:calcRegion(end),i));
            inputDensityArrayS(1,:,i) = inputDensityArrayS(1,:,i) - shiftByThis(i);
        end        
    end
end

kappaTildeArray = []; pTildeArray = []; kfs = []; logkfa2ds = []; potentials = [];
for i=1:length(inputDensityArray(1,1,:)) %i is each cloud

%radproftoFityReal = radProfilesAvg(1,:,1)./(kpixelLength^2); %convert to density/m^-2
radProfilesDensityAdjusted = inputDensityArrayS(1,:,i); %density already in m^-2
radiusVector = inputRadiusArray(1,:,i); %x vector 
radiusVectorMeters = radiusVector.*kpixelLength; %convert x vector to m
potentials(:,i) = 0.5 .* massL6 .* omegaR^2 .* radiusVectorMeters.^2;
%potential_nk = (potentials(:,i)./kB).*(10^9);
%plot(potential_nk,radProfilesAvg(1,:,i));


%------------kappa and p calculations:
kappaTildeAvg = []; pTildeAvg = []; pAvg = []; dydx = []; 

%dydx = diff([eps radProfilesDensityAdjusted(calcRegion)])./diff([eps potentials(calcRegion,i)']);
%Default:
%kappaTildeAvg = (-1)*(pi * hbar^2 / massL6) * gradient(radProfilesDensityAdjusted(calcRegion),potentials(calcRegion,i));
%With Smooth:
%if(smoothOn)
%kappaTildeAvg = (-1)*(pi * hbar^2 / massL6) * gradient(smooth(radProfilesDensityAdjusted(calcRegion),8),potentials(calcRegion,i));
%else
kappaTildeAvg = (-1)*(pi * hbar^2 / massL6) * gradient(radProfilesDensityAdjusted(calcRegion),potentials(calcRegion,i));
%end
%kappaTildeAvg = (-1)*(pi * hbar^2 / massL6) .* dydx;

for j=1:length(calcRegion)
    %Default:
    %pAvg(j) = trapz(potentials(calcRegion(j):calcRegion(end)+1,i),radProfilesDensityAdjusted(calcRegion(j):calcRegion(end)+1));
    %With Smooth:
    %if(smoothOn)
    %pAvg(j) = trapz(potentials(calcRegion(j):calcRegion(end)+1,i),smooth(radProfilesDensityAdjusted(calcRegion(j):calcRegion(end)+1),8)); %integral from end to j point
    %else
    pAvg(j) = trapz(potentials(calcRegion(j):calcRegion(end)+1,i),radProfilesDensityAdjusted(calcRegion(j):calcRegion(end)+1));
    %end
    n = radProfilesDensityAdjusted(calcRegion(j));     %n, density on j'th point (varies across the potential)
    pTildeAvg(j) = ((2*massL6)/(n^2 * hbar^2 * pi)) * pAvg(j);
       
end
%-------------Beta_mu calc

betamu_initial = pTildeAvg(end);
for j=1:length(calcRegion)
    

    
    
    
    
end



kfs(:,i) = real(sqrt(2*pi*radProfilesDensityAdjusted(calcRegion(1):calcRegion(end))));
logkfa2ds(:,i) = log(kfs(:,i).*a2d);

kappaTildeArray(:,i) = kappaTildeAvg(1:end);
pTildeArray(:,i) = pTildeAvg(1:end);

end


    figure(2000);
    hold on;
    for i=1:length(kappaTildeArray(1,:))
        %plot(pTildeArray(3:end,i),kappaTildeArray(3:end,i),'.');
        plot(logkfa2ds(3:end,i),pTildeArray(3:end,i),'.');
    end
    %axis([-2 14 -2 8]);
    axis([0 5 0 8]);
    title(['ptilde vs log(kfa2d) 972G individual clouds']);
    
    kappaTildeAll = reshape(kappaTildeArray, 1, length(kappaTildeArray(:,1))*length(kappaTildeArray(1,:)));
    pTildeAll = reshape(pTildeArray, 1, length(pTildeArray(:,1))*length(pTildeArray(1,:)));    
    logkfa2dsAll = reshape(logkfa2ds, 1, length(logkfa2ds(:,1))*length(logkfa2ds(1,:)));
    kfsAll = reshape(kfs, 1, length(kfs(:,1))*length(kfs(1,:)));
    potentialsAll = reshape(potentials(calcRegion(1):calcRegion(end),:), 1, length(potentials(calcRegion(1):calcRegion(end),1))*length(potentials(1,:)));
    
    %figure(15); plot(pTildeAll,kappaTildeAll,'.'); axis([-2 14 -2 6]);
    %figure(11); plot(potentialsAll,kfsAll,'.');
    %figure(12); plot(potentialsAll,logkfa2dsAll,'.');
    %figure(13); plot(logkfa2dsAll,pTildeAll,'.'); axis([0 5 0 12]);
    %figure(14); scatter3(logkfa2dsAll,pTildeAll,kappaTildeAll); axis([0 5 0 12 -2 6]);
    
    %truncate pTildeAll kappaTildeAll for sensible points:
    j = 1; pTildeTrun = []; kappaTildeTrun = []; logkfa2dsTrun = [];
    for i=1:length(pTildeAll)
        if(pTildeAll(i) > 0)
            if(pTildeAll(i) < 20)
                if(logkfa2dsAll(i) > 0.1)
                    pTildeTrun(j) = pTildeAll(i);
                    kappaTildeTrun(j) = kappaTildeAll(i);
                    logkfa2dsTrun(j) = logkfa2dsAll(i);
                    j = j+1;
                end
            end
        end
    end
    
    [sortedpTilde,indexs] = sort(pTildeTrun,'ascend');
    sortedkappaTilde = kappaTildeTrun(indexs);
    
    n = 80;
    %pTildeElemAvg = meanNelements(sortedpTilde,n);
    %kappaTildeElemAvg = meanNelements(sortedkappaTilde,n);
    
    
    plot(pTildeTrun,kappaTildeTrun,'.'); axis([-2 14 -2 6]);
    
    close all;
    binResult = binMe(kappaTildeTrun,pTildeTrun,25);
    %figure(3); plot(binResult(2,:),binResult(1,:),'.'); axis([-2 14 -2 6]);
    %figure(4); errorbar(binResult(2,:),binResult(1,:),binResult(3,:)./2,'.'); axis([-3 15 -3 7]);
    %figure(5); plot(pTildeTrun,kappaTildeTrun,'.r'); hold on; errorbar(binResult(2,:),binResult(1,:),binResult(3,:)./2,'.'); axis([-3 15 -3 7]); grid on; hold off;
    %figure(6); plot(pTildeElemAvg,kappaTildeElemAvg,'.'); axis([-3 15 -3 7]); grid on;
    %binElemResult = binMe(kappaTildeElemAvg,pTildeElemAvg,25);
    %figure(7); plot(pTildeElemAvg,kappaTildeElemAvg,'.'); hold on; errorbar(binElemResult(2,:),binElemResult(1,:),binElemResult(3,:)./2,'.r'); axis([-3 15 -3 7]); grid on; hold off;
    %figure(15); scatter3(logkfa2dsTrun,pTildeTrun,kappaTildeTrun); axis([0 5 0 12 -2 6]);
    
close all; 

%figure(8); plot(potentials(20,calcRegion),kfs(20,calcRegion),'.');
%figure(9); plot(potentials(20,calcRegion),inputDensityArray(1,calcRegion,20),'.');
%figure(10); plot(potentials(20,calcRegion),logkfa2ds(20,calcRegion),'.');

%3D Binning:
k_p_Bin = binMe(kappaTildeTrun,pTildeTrun,25);
a2d_p_Bin = binMe(logkfa2dsTrun,pTildeTrun,25);

%%%%%%%%%%%%%%%%% BEGIN 4-way figure (iso + 3 perspectives)
close all;
fa = figure(100);

subplot(2,3,1);
hs = scatter3(logkfa2dsTrun,pTildeTrun,kappaTildeTrun,5); hold on;
scatter3(a2d_p_Bin(1,:),a2d_p_Bin(2,:),k_p_Bin(1,:)); 
plot3(a2d_p_Bin(1,:),a2d_p_Bin(2,:),k_p_Bin(1,:)); axis([0 5 0 12 -2 6]); grid on;
xlabel('log(kf*a2d)') % x-axis label
ylabel('pTilde') % y-axis label
zlabel('kappaTilde') % y-axis label
%hmarkers = hs.MarkerHandle;
%hmarkers.EdgeColorData = uint8(255*[0;0;1;0.3]);

subplot(2,3,2);
scatter3(logkfa2dsTrun,pTildeTrun,kappaTildeTrun,5); hold on;
scatter3(a2d_p_Bin(1,:),a2d_p_Bin(2,:),k_p_Bin(1,:)); 
plot3(a2d_p_Bin(1,:),a2d_p_Bin(2,:),k_p_Bin(1,:)); axis([0 5 0 12 -2 6]); grid on;
xlabel('log(kf*a2d)') % x-axis label
ylabel('pTilde') % y-axis label
zlabel('kappaTilde') % y-axis label
view([1 0 0]);

subplot(2,3,4);
scatter3(logkfa2dsTrun,pTildeTrun,kappaTildeTrun,5); hold on;
scatter3(a2d_p_Bin(1,:),a2d_p_Bin(2,:),k_p_Bin(1,:)); 
plot3(a2d_p_Bin(1,:),a2d_p_Bin(2,:),k_p_Bin(1,:)); axis([0 5 0 12 -2 6]); grid on;
xlabel('log(kf*a2d)') % x-axis label
ylabel('pTilde') % y-axis label
zlabel('kappaTilde') % y-axis label
view([0 1 0]);

subplot(2,3,5);
scatter3(logkfa2dsTrun,pTildeTrun,kappaTildeTrun,5); hold on;
scatter3(a2d_p_Bin(1,:),a2d_p_Bin(2,:),k_p_Bin(1,:)); 
plot3(a2d_p_Bin(1,:),a2d_p_Bin(2,:),k_p_Bin(1,:)); axis([0 5 0 12 -2 6]); grid on;
xlabel('log(kf*a2d)') % x-axis label
ylabel('pTilde') % y-axis label
zlabel('kappaTilde') % y-axis label
view([0 0 1]);
hold off;

subplot(2,3,3);
indvToPlot = 1;
if(smoothOn)
    plot(inputRadiusArray(1,calcRegion,indvToPlot),smooth(inputDensityArrayS(1,calcRegion,indvToPlot),smoothAmount)); hold on;
    if(zeroOn)
        plot(inputRadiusArray(1,calcRegion,indvToPlot),inputDensityArray(1,calcRegion,indvToPlot)-shiftByThisMin(indvToPlot));
    else
        plot(inputRadiusArray(1,calcRegion,indvToPlot),inputDensityArray(1,calcRegion,indvToPlot));
    end
else
    plot(inputRadiusArray(1,calcRegion,indvToPlot),inputDensityArrayS(1,calcRegion,indvToPlot));
end
grid on;


subplot(2,3,6);
hold on;
for i=1:length(kappaTildeArray(1,:))
    %plot(pTildeArray(3:end,i),kappaTildeArray(3:end,i),'.');
    plot(logkfa2ds(1:end,i),pTildeArray(1:end,i));
end
%axis([-2 14 -2 8]);
axis([0 5 0 8]);
grid on;



ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
%titlestring = [num2str(field) 'G each cloud data scatter with binned values'];
titlestring = [num2str(field) 'G' savestring];
text(0.5, 1,titlestring,'HorizontalAlignment','center','VerticalAlignment', 'top');
set(fa, 'Position', [100, 100, 1400, 800]);

figname = [num2str(field) '_scatterGraph3D_Iso_' savestring];
figdirectory = 'C:\Users\tpeppler\Dropbox\PhD\2D_2016\EOS_Data\';
if(saveOn)
saveas(fa,[figdirectory figname '.fig'],'fig');
saveas(fa,[figdirectory figname '.png'],'png');
end


%%%%%%%%%%%%%%%%% END 4-way figure (iso + 3 perspectives)

h = figure(5); plot(pTildeTrun,kappaTildeTrun,'.r'); hold on; errorbar(binResult(2,:),binResult(1,:),binResult(3,:)./2,'.'); axis([-3 15 -3 7]); grid on; hold off;
figname = [num2str(field) 'raw_and_binned_kappaTildevsPTilde' savestring];
figdirectory = 'C:\Users\tpeppler\Dropbox\PhD\2D_2016\EOS_Data\';
if(saveOn)
%saveas(h,[figdirectory figname '.fig'],'fig');
%saveas(h,[figdirectory figname '.png'],'png');
end

end
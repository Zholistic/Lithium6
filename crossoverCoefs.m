
crossData = [];
%[date magneticField coef2d coef3d]
crossData(:,1) = [140907 850 0.21 0.13];
crossData(:,2) = [140905 850 0.215 0.14];
crossData(:,3) = [140902 880 0.215 0.116];
crossData(:,4) = [140901 972 0.2 0.126];
crossData(:,5) = [140826 920 0.226 0.13];
crossData(:,6) = [140819 972 0.21 0.127];
crossData(:,7) = [140819 880 0.217 0.12];

plot(crossData(2,:),crossData(3,:),'MarkerFaceColor',[0.6 0.6 1],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[0 0 1]);
hold on;
plot(crossData(2,:),crossData(4,:),'MarkerFaceColor',[1 0 0],...
    'Marker','o',...
    'LineStyle','none',...
    'Color',[1 0 0]);
hold off;
grid on;
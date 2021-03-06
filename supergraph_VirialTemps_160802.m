%Supergraph 3D surface plot:
close all;
open('C:\Users\tpeppler\Dropbox\PhD\2D_2016\EOS_Data\865_TonTF_vs_logkFa2d_VirialFit.fig');
handles(1) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2D_2016\EOS_Data\880_TonTF_vs_logkFa2d_VirialFit.fig');
handles(2) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2D_2016\EOS_Data\920_TonTF_vs_logkFa2d_VirialFit.fig');
handles(3) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2D_2016\EOS_Data\972_TonTF_vs_logkFa2d_VirialFit.fig');
handles(4) = gcf;
open('C:\Users\tpeppler\Dropbox\PhD\2D_2016\EOS_Data\980_TonTF_vs_logkFa2d_VirialFit.fig');
handles(5) = gcf;


xdatas = []; ydatas = []; zdatas = [];
for i=1:5
    axesObjs = get(handles(i), 'Children');  %axes handles
    dataObjs = get(axesObjs(2), 'Children'); %handles to low-level graphics objects
    objType = get(dataObjs, 'Type');  %type of low-level graphics object
    xdatas(:,i) = get(dataObjs(1), 'XData');  %data from low-level grahics objects
    ydatas(:,i) = get(dataObjs(1), 'YData');
    %ldatas{i} = get(dataObjs,'Ldata');
end


figure(22);
hold on;
for i=1:5
    plot(xdatas(:,i),ydatas(:,i));
end
grid on;
axis([0 4 0 2]); title('T/TF vs log(kF a2d)');
legend('show');
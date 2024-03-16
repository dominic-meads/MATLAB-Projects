% pressure detect test
clear 
close all
clc
%% 
load ICP.mat
fi = PressureDetect(icp2,[],1);

fs = 125;
t = (0:length(icp2)-1)/fs;

figure('Color',[1 1 1]);
h = plot(t,icp2);
set(h,'LineWidth',1.2);
set(h,'Color',[0.2 0.73 1]);
hold on;
h = plot(fi/fs,icp2(fi),'r.');
set(h,'MarkerSize',18);
title('ICP Signal with Location of Percussion Peaks');

ylabel('ICP (mmHG)');
xlabel('Time (s)');
xlim([min(t) max(t)]);
hold on;
h = plot(dDT2/fs,icp2(dDT2),'r*');
set(h,'MarkerSize',10);
title('Expert Annotated vs Algorithm Detected Percussion Peaks')
legend('ICP', 'P1 (detected)', 'P1 (annotated)');
xlim([15 28])
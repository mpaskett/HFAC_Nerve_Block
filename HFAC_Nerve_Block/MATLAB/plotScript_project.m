% Generate most figures for project paper

clearvars
close all
axisFontSize = 36;
tickLabelFontSize = 26;
lgdFontSize = 20;
%% Generate AC Stim Plot (Myelinated)
threshArray = [ 180 053 058 073 ;   % 2 uM diameter
                156 046 047 056 ;   % 4 uM diameter
                134 029 028 031 ;   % 8 uM diameter
                160 050 046 048 ];  % 16 uM diameter    
            
freq        = [ 2 4 8 16 ];         % AC stim frequency kHz

figure(1)
semilogx(freq, threshArray','-o','LineWidth',1.4)
xlabel('Frequency (kHz)','FontSize',axisFontSize)
ylabel('Block Threshold (uA)','FontSize',axisFontSize)
axizer = gca;
axizer.FontSize = tickLabelFontSize;
axizer.YGrid = 'on';
leg = legend('2 (\mum)', '4 (\mum)', '8 (\mum)', '16 (\mum)');
title(leg,'Axon Diameter')
leg.FontSize = lgdFontSize;

%% Generate AC Stim Plot (Unmyelinated)
threshArray = [ 8.3 15.4 27.4 51.2; 
                4.6  8.2 14.9 28.5];
            
freq        = [ 2 4 8 16 ];         % AC stim frequency kHz

figure(2)
semilogx(freq, threshArray','-o','LineWidth',1.4)
xlabel('Frequency (kHz)','FontSize',axisFontSize)
ylabel('Block Threshold (mA)','FontSize',axisFontSize)
axizer = gca;
axizer.FontSize = tickLabelFontSize;
axizer.YGrid = 'on';
leg = legend('1 (\mum)', '2 (\mum)');
title(leg,'Axon Diameter')
leg.FontSize = lgdFontSize;

%% Generate DC Stim Plot

unMyThresh = [1350 1340];
unMyDiam   = [1.00 2.00];

myThresh   = [079 080 073 087];
myDiam     = [2.00 4.00 8.00 16.0];

figure(3)
loglog(unMyDiam, unMyThresh, '-o','LineWidth',1.4)
hold on
loglog(myDiam, myThresh,'-o','LineWidth',1.4)
xlabel('Axon Diameter (\mum)','FontSize',axisFontSize)
ylabel('Block Threshold (uA)','FontSize',axisFontSize)

xticks = get(gca,'XTick');
axizer = gca;
axizer.FontSize = tickLabelFontSize;
axizer.YGrid = 'on';
axizer.XTickLabel = xticks;
ylim([30 3000])
xlim([0 20])

leg = legend('Unmyelinated', 'Myelinated');
leg.FontSize = lgdFontSize;
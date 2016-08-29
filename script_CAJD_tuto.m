% example of convergence of Approximate Joint DIagonalization on data
% generated from a non-orthogonal bilinear model
% using three methods : AJD , BAJD and CAJD both with GPT optimization. See
% article below for details. 
%
% *** History: 09-Fev-2016
% *** Author: Florent BOUCHARD & Louis KORCZOWSKI, GIPSA-Lab, 2016
% *** Reference: MINING THE BILINEAR STRUCTURE OF DATA WITH APPROXIMATE
%                JOINT DIAGONALIZATION
%                L. Korczowski, F. Bouchard, C. Jutten, M. Congedo
% *** Contact: louis.korczowski@gmail.com
% *** Licence: GNU GPLv3
%
% Tested compatibility: Matlab 2015a
%
% Expected results ;
% If the data respect a bilinear model, AJD should perform worse than BAJD
% and CAJD. On the other hand, CAJD should perform better than BAJD on the
% estimation of the spatial mixing matrix thanks to the use of the
% diversity.

clear all;
close all;
addpath('./lib/')
%% define global parameters
%
nTest = 1; %number of simulations
N     = 8; %spatial matrix size
T     = 128; %temporal matrix size
K     = 100; %number of asymmetrical observations matrices
F     = 20; %number of simulated symmetric definite positive matrices
condA = 1; %true spatial mixing matrix conditioning
allCondE = {[50 100]}; % true range of temporal mixing matrix conditioning
allSNR   = {inf,20}; % signal-noise-ratio
outputDIR = '.\figures\';
%% perform Test on CAJD model
for tix=1:nTest
    disp(tix);
    for eix=1:length(allCondE)
        for six=1:length(allSNR)
            %% define current test parameters
            condE  = allCondE{eix};
            SNR    = allSNR(six);
            
            %% simulate data
            opSim = struct('N',N,'T',T,'K',K,'F',F,'condA',condA,'condE',condE,'SNR',SNR);
            [X,Cf,A,E] = simBAJD_dat(opSim);
            cond(E)
            % get all covariance matrices for AJD
            [C,Ct] = convertDat_BAJD(X,Cf);
            [U,S,V]=svd(mean(X,3)',0);
             V=V(:,1:N);
            %% define BSS parameters
            B0      = inv(V)';%eye(N);
            D0      = pinv(U)';%orth(randn(T,N));
            epsilon = 1e-18;
            itMax   = 50;
            
            %% perform CAJD
            opCAJD = struct('B0',B0,'D0',D0,'eps',epsilon,'itMax',itMax,'A',A,'E',E);
            [B_CAJD,D_CAJD,S_CAJD,C_CAJD,info_CAJD] = gp_CAJD_GPT(X,Cf,opCAJD);
            
            %% perform AJD
            opAJD = struct('B0',B0,'eps',epsilon,'itMax',itMax,'A',A);
            [B_AJD,C_AJD,info_AJD] = gp_AJD_GPT(C,opAJD);
            %% perform BAJD
            opBAJD = struct('B0',B0,'D0',D0,'eps',epsilon,'itMax',itMax,'A',A,'E',E);
            [B_BAJD,D_BAJD,S_BAJD,info_BAJD]= gp_BAJD_GPT(X,opBAJD);
            
            %% save results
            results.critA_CAJD{tix,eix,six} = [info_CAJD.critA];
            results.critE_CAJD{tix,eix,six} = [info_CAJD.critE];
            results.critA_BAJD{tix,eix,six} = [info_BAJD.critA];
            results.critE_BAJD{tix,eix,six} = [info_BAJD.critE];
            results.critA_AJD{tix,eix,six} =  [info_AJD.critA];
        end
    end
end
save( [outputDIR 'CAJD_sim'],'results');
% plot results for the Moreau-Macchi Criterion on B
close all
critA_CAJD = results.critA_CAJD;
critE_CAJD = results.critE_CAJD;
critA_BAJD = results.critA_BAJD;
critE_BAJD = results.critE_BAJD;
critA_AJD = results.critA_AJD;

limx = [1 39.5];
limy1 = [ -16 1];
limy2 = [ -3 1];
FontSize=12;
Paper=[0 0 15 10];
LineWidth=1.5;
MarkersSize=6;
Colors={[1 0 0],[0 0 1],[0.8 0.8 0],[0 0.8 .8]};
Legends={'AJD_{inf}','CAJD_{inf}','BAJD_{inf}','AJD_{100}','CAJD_{100}','BAJD_{100}'}
figure
subplot1(1,2,'YTickL','All')
for tix=1:nTest
    subplot1(1)
    plot(log10(critA_AJD{tix,1,1}),'x-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{1});
    hold on;
    plot(log10(critA_CAJD{tix,1,1}),'x-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{2});
    plot(log10(critA_BAJD{tix,1,1}),'x-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{3});
    
    
     ylabel('$I_{{M-M}}~~(\mathbf{B})$ ~(dB)', 'interpreter','latex','fontsize',FontSize);
       
           legend(Legends{1:3},'location','best');
     
    xlim(limx)
    ylim(limy1)
    xlabel('number of sweeps');
    title('(a)')
    
    subplot1(2)
    plot(log10(critA_AJD{tix,1,2}),'^-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{1});
    hold on;
    plot(log10(critA_CAJD{tix,1,2}),'^-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{2});
    plot(log10(critA_BAJD{tix,1,2}),'^-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{3});
    
    legend(Legends{4:6},'location','southwest');
    set(gca,'yaxislocation','right');
    %     set(gca,'yticklabel')
    ylim(limy2)
    xlabel('number of sweeps');
    title('(b)')
    
end
xlim(limx)
set(gcf,'paperposition',Paper)

% legend(Legends,'location','southwest');
set(gcf,'color',[1 1 1])
print(gcf, [outputDIR 'MoMaA'],'-dpng','-r450')


% plot results for the Moreau-Macchi Criterion on D

limx = [1 39.5];
limy1 = [ -16 1];
limy2 = [ -5 1];
FontSize=12;
Paper=[0 0 15 10];
LineWidth=1.5;
MarkersSize=6;
Colors={[1 0 0],[0 0 1],[0.8 0.8 0],[0 0.8 .8]};
Legends={'AJD_{inf}','CAJD_{inf}','BAJD_{inf}','AJD_{100}','CAJD_{100}','BAJD_{100}'};
figure
subplot1(1,2,'YTickL','All')
for tix=1:nTest
    subplot1(1)
    hold on;
    plot(log10(critE_CAJD{tix,1,1}),'x-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{2});
    plot(log10(critE_BAJD{tix,1,1}),'x-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{3});
    
    
    ylabel('$I_{{M-M}}~~(\mathbf{D})$ ~(dB)', 'interpreter','latex','fontsize',FontSize);
    legend(Legends{2:3},'location','southwest');
    xlim(limx)
    ylim(limy1)
    xlabel('number of sweeps');
    title('(a)')
    
    subplot1(2)
    hold on;
    plot(log10(critE_CAJD{tix,1,2}),'^-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{2});
    plot(log10(critE_BAJD{tix,1,2}),'^-','Markers',MarkersSize,'Linewidth',LineWidth,'color',Colors{3});
    
    legend(Legends{5:6},'location','southwest');
      set(gca,'yaxislocation','right');
     ylim(limy2)
     xlabel('number of sweeps');
    title('(b)')
    
end
xlim(limx)
set(gcf,'paperposition',Paper)

set(gcf,'color',[1 1 1])
print(gcf, [outputDIR '\MoMaE'],'-dpng','-r450')
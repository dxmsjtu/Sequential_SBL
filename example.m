% Version 1.0: (03/29/2021)
% written by Yongsung Park

% Yongsung Park, Florian Meyer, & Peter Gerstoft
% MPL/SIO/UCSD
% yongsungpark@ucsd.edu / flmeyer@ucsd.edu / gerstoft@ucsd.edu
% noiselab.ucsd.edu

% Citation
% Y. Park, F. Meyer, and P. Gerstoft, “Sequential sparse Bayesian learning for time-varying direction of arrival,” J. Acoust. Soc. Am. 149(3) (2021).
% https://doi.org/10.1121/10.0003802
% Y. Park, F. Meyer, and P. Gerstoft, “Sequential sparse Bayesian learning for DOA,” in Proc. Asilomar Conf. Signals Syst. Comput. (2020).
% https://doi.org/10.1109/IEEECONF51394.2020.9443444

%%
clear; clc;
close all;

% rngN = 0;
% rng(rngN)

dbstop if error;

% addpath([cd,'/_common'])

% SNRlist = [40 35 30 25 20 15 10 5 0];
% SNRlist = [40 35 30 25 20 15 10 5];
SNRlist = [20];
for nSNR = 1:length(SNRlist)
Nsim = 1;
for nsim=1:Nsim
Nrng=0; rng(Nrng+nsim)
% rng(nsim)
disp(['SNR',num2str(SNRlist(nSNR)),'_',num2str(nsim)])
% Environment parameters
c = 1500;       % speed of sound
f = 200;        % frequency
lambda = c/f;   % wavelength

% ULA-horizontal array configuration
Nsensor = 15;               % number of sensors
d = 1/2*lambda;             % intersensor spacing
q = (0:1:(Nsensor-1))';     % sensor numbering
xq = (q-(Nsensor-1)/2)*d;   % sensor locations

% sensor configuration structure
Sensor_s.Nsensor = Nsensor;
Sensor_s.lambda = lambda;
Sensor_s.d = d;
Sensor_s.q = q;
Sensor_s.xq = xq;

% signal generation parameters
SNR = SNRlist(nSNR);

% total number of snapshots
Nsnapshot = 50;

% range of angle space
thetalim = [-90 90];

theta_separation = 0.5;

% Angular search grid
theta = (thetalim(1):theta_separation:thetalim(2))';
Ntheta = length(theta);

% Design/steering matrix (Sensing matrix)
sin_theta = sind(theta);
sensingMatrix = exp(-1i*2*pi/lambda*xq*sin_theta.')/sqrt(Nsensor);

% Generate received signal
% anglesTrue = [-3; 2; 60]; % DOA of sources at first snapshot [deg]
% source_amp = [  7;   7;  7;  4;  4; 13];
anglesTrue = [-70; -55; -40; 35; 50; 65]; % DOA of sources at first snapshot [deg]
anglesTracks = repmat(anglesTrue,[1,Nsnapshot]);
anglesTracks(3,:) = anglesTracks(3,1) - 2*anglesTracks(3,1)./(1+exp(-.1*(-Nsnapshot/2:-Nsnapshot/2+Nsnapshot-1)));
anglesTracks(4,:) = anglesTracks(4,1) - 1.00*(0:Nsnapshot-1)';
sinAnglesTracks = sind(anglesTracks); 
Nsources = numel(anglesTrue);

receivedSignal = zeros(Nsensor,Nsnapshot);
source_amp = zeros(Nsources,Nsnapshot);

for snapshot = 1:Nsnapshot
    ActSrcInd{snapshot} = [1;2;3;4;5;6];
%     if snapshot <= 15
%         ActSrcInd{snapshot} = [1;2;3;4];
%     elseif snapshot <= 30
%         ActSrcInd{snapshot} = [1;2;5;6];
%     else
%         ActSrcInd{snapshot} = [3;4;5;6];
%     end
end

for snapshot = 1:Nsnapshot

    srcind = ActSrcInd{snapshot};
    
    % Source amplitude
%     source_amp(:,snapshot) = complex(randn(size(anglesTrue)),randn(size(anglesTrue)))/sqrt(2);
%     Xsource = source_amp(srcind,snapshot);
    source_amp = [  7;   7;  7;  4;  4; 13];
    Xsource = source_amp.*exp(1i*2*pi*rand(Nsources,1));    % random phase
    Xsource = Xsource(srcind);
    
    % Represenation matrix (steering matrix)
    transmitMatrix = exp( -1i*2*pi/lambda*xq*sinAnglesTracks(srcind,snapshot).' )/sqrt(Nsensor);
    
    % Received signal without noise
    receivedSignal(:,snapshot) = sum(transmitMatrix*diag(Xsource),2);
    
    % add noise to the signals
    rnl = 10^(-SNR/20)*norm(Xsource);
    nwhite = complex(randn(Nsensor,1),randn(Nsensor,1))/sqrt(2*Nsensor);
    
    e = nwhite * rnl;	% error vector
    receivedSignal(:,snapshot) = receivedSignal(:,snapshot) + e;
end

%% Conventional beamforming (CBF)
outputBeamformer = sensingMatrix' * receivedSignal;

for snapshot = 1:Nsnapshot
snapshot
srcind = ActSrcInd{snapshot};
%% Non-Sequential SBL
    options = SBLSet();
    options.convergence.error = 10^(-3);
    options.Nsource = ceil(Nsensor/2);
    
    threshold = .01;
    
    [gamma, reportSBL] = SBL_v3p12( sensingMatrix, receivedSignal(:,snapshot), options );
    [mu_val_SBL,peak_SBL] = findpeaks(gamma,'SORTSTR','descend','Npeaks', Nsensor);
    theta_filter_SBL = (theta(peak_SBL(abs(mu_val_SBL)>threshold*max(abs(mu_val_SBL)))));
    gamma_filter_SBL = abs( mu_val_SBL(abs(mu_val_SBL)>threshold*max(abs(mu_val_SBL))) );
    
    outputSBL = struct('gamma',gamma,'doas',theta_filter_SBL,'amps',gamma_filter_SBL);
    
    if nsim ==1 && snapshot ==1, outputsSBL = []; end
    outputsSBL = [outputsSBL;outputSBL];

%% Sequential SBL
%     options = SBLSet();
%     options.convergence.error = 10^(-3);
%     options.Nsource = ceil(Nsensor/2);
%     
%     threshold = .01;
    
    if snapshot==1
        [gammaSeq, reportSeqSBL] = SBL_v3p12( sensingMatrix, receivedSignal(:,snapshot), options );
    else
        [gammaSeq, reportSeqSBL] = SBLseq_v3p12( sensingMatrix, receivedSignal(:,snapshot), gammaSeq, options );
    end
    [mu_val_SBLseq,peak_SBLseq] = findpeaks(gammaSeq,'SORTSTR','descend','Npeaks', Nsensor);
    theta_filter_SBLseq = (theta(peak_SBLseq(abs(mu_val_SBLseq)>threshold*max(abs(mu_val_SBLseq)))));
    gamma_filter_SBLseq = abs( mu_val_SBLseq(abs(mu_val_SBLseq)>threshold*max(abs(mu_val_SBLseq))) );
    
    outputSBLseq = struct('gamma',gammaSeq,'doas',theta_filter_SBLseq,'amps',gamma_filter_SBLseq);
    
    if nsim ==1 && snapshot ==1, outputsSSBL = []; end
    outputsSSBL = [outputsSSBL;outputSBLseq];

end
end
% vars = who();
% varnames = vars(contains(vars, 'outputs'));
% save(['D_AD_CG6s_SNR',num2str(SNR),'_',num2str(nsim)],...
%     'anglesTracks',varnames{:});
end
% save('data_results_AD')

%% Plot: Non-sequential SBL
outputPlot = outputsSBL;
figure;
set(gcf,'position',[562,250,560,420]);
imagesc(1:Nsnapshot,theta,-inf);
caxis([-20 0])

rtIndex = []; rtTheta = []; rtAmp = [];
for index=1:Nsnapshot
    rTheta  = outputPlot(index).doas;
    rAmp    = outputPlot(index).amps;
    rAmp    = 10*log10( rAmp / max(rAmp) );
    
    rtIndex = [rtIndex;index*ones(size(rAmp))];
    rtTheta = [rtTheta;rTheta];
    rtAmp = [rtAmp;rAmp];
end

hold on; scatter(rtIndex,rtTheta,50,rtAmp,'filled','o','linewidth',.5,...
    'MarkerEdgeColor','k'); hold off;

title('Non-sequential SBL','interpreter','latex')
xlabel('Time step','interpreter','latex')
ylabel('DOA~[$^\circ$]','interpreter','latex')
box on
set(gca,'fontsize',18,'YDir','normal','TickLabelInterpreter','latex','YTick',-80:40:80)
axis([.5 Nsnapshot+.5 -90 90])
set(gca,'YTickLabel','')

%% Plot: Sequential SBL
outputPlot = outputsSSBL;
figure;
set(gcf,'position',[1124,250,560,420]);
imagesc(1:Nsnapshot,theta,-inf);
caxis([-20 0])

rtIndex = []; rtTheta = []; rtAmp = [];
for index=1:Nsnapshot
    rTheta  = outputPlot(index).doas;
    rAmp    = outputPlot(index).amps;
    rAmp    = 10*log10( rAmp / max(rAmp) );
    
    rtIndex = [rtIndex;index*ones(size(rAmp))];
    rtTheta = [rtTheta;rTheta];
    rtAmp = [rtAmp;rAmp];
end

hold on; scatter(rtIndex,rtTheta,50,rtAmp,'filled','o','linewidth',.5,...
    'MarkerEdgeColor','k'); hold off;

title('Sequential SBL','interpreter','latex')
xlabel('Time step','interpreter','latex')
ylabel('DOA~[$^\circ$]','interpreter','latex')
box on
set(gca,'fontsize',18,'YDir','normal','TickLabelInterpreter','latex','YTick',-80:40:80)
axis([.5 Nsnapshot+.5 -90 90])
set(gca,'YTickLabel','')

%% Conventional beamforming (CBF)
figure;
set(gcf,'position',[1,250,560,420]);
imagesc(1:Nsnapshot,theta,10*log10( (abs(outputBeamformer).^2) ./ max((abs(outputBeamformer).^2),[],1)));
caxis([-20 0])

title('Conventional beamforming (CBF)','interpreter','latex')
xlabel('Time step','interpreter','latex')
ylabel('DOA~[$^\circ$]','interpreter','latex')
box on
set(gca,'fontsize',18,'YDir','normal','TickLabelInterpreter','latex','YTick',-80:40:80)
axis([.5 Nsnapshot+.5 -90 90])

% Lexy von Diezmann, 2023-2024. Released under the GNU GPL v3.

function [kymo,options] = brownianMotionInSC_plotTracks(path) % qObs,D,frameTime, options,stickTime)

load(path,'qObs','D','frameTime','options','stickTime','nameStem');

% qObs is nSteps x nMols of 1D positions
% D is either a scalar or has length qObs
% frameTime is scalar in seconds
% option is "1" (start at point position) or "2" (start uniform, bind)
% stickTime is nMols x 1 of index where enter binding region (opt 2 only)

% name=[nameStem '_' num2str(setNum)];
% load(name,'frameTime','L','options','nMols','D');

%% conserved parameters

    frameTime = 0.2; % seconds. resolution of qObs.
    dt = 0.05; % seconds - sampling rate of raw simulation
    dispN = 100; % display progress in increments of dispN mols simulated
    nDims = 1; % number of dimensions (1 to 3)
    L = 3; % half-length in um
    r = 0.05; % radius of curvature (can be empty) - um
    er = 0; % isotropic localization error in um %40e-3;


    name = [nameStem '_opt_' num2str(options)];
%% define color space

colors = [0, 0.4470, 0.7410;...
          0.8500, 0.3250, 0.0980;...
          0.9290, 0.6940, 0.1250;...
          0.4940, 0.1840, 0.5560;...
          0.4660, 0.6740, 0.1880;...
          0.3010, 0.7450, 0.9330;...
          0.6350, 0.0780, 0.1840];
numtimepoints = size(qObs,1);



%% set up cell coordinates and plot examples nicely (1D)
    tracksToPlot = 17;%[2 3 6]; %     tracksToPlot = [10 4 8];

    h1=figure; hold on;
    for i = 1:length(tracksToPlot)
        mol = tracksToPlot(i);
        cval = 1+mod(i+1,size(colors,1));
        plot(qObs(:,mol),(1:size(qObs,1))*frameTime,'color',colors(cval,:))
        scatter(qObs(1,mol),1,50,colors(cval,:),'filled');
        scatter(qObs(end,mol),size(qObs,1)*frameTime,50,colors(cval,:),'filled');
        axis ij
    end
    xlim([-L,L]);
    xlim([-L-0.05,L+0.05]);
    title(path);
saveas(h1,[name '_example.fig'],'fig')

%% kymographs

binSize = 0.02; % um
nTimes = size(qObs,1);
edges = (-L-binSize/2):binSize:(L+binSize/2);
kymo = nan(size(qObs,1),length(edges)-1); % nTimes x Edges


for i = 1:nTimes
    kymo(i,:) = histcounts(qObs(i,:),edges);
end

h2 = figure('Position', [10 10 560 160]);
imagesc(kymo); colormap parula
title(path)

% add Gaussian filter to improve visual inspection of results
filt = [1 1 1 1 1]' * [0.0219    0.0983    0.1621    0.0983    0.0219];
filt = filt/sum(filt);
h4 = figure('Position', [10 10 560 160]);
imagesc(imfilter(kymo,filt));
saveas(h4,[name '_blurkymo.fig'],'fig')

% h3 = figure('Position', [10 10 560 160]);
% imagesc(log(kymo+1)); colormap gray
% title(path)ol

saveas(h2,[name '_kymo.fig'],'fig')
% saveas(h3,[name '_logkymo.fig'],'fig')
% saveas(h2,[nameStem '_example'],'epsc')
% figure; histogram(qObs(i,:),edges)

%%s
% 
% 
% region = nan(numtimepoints,3);
% 
% 
% switch options
%     case 1
%         lowerlim=-L/3;
%         upperlim=L/3;
%         for i = 1:numtimepoints
%             region(i,1) = sum(qObs(i,:)<lowerlim);
%             region(i,2) = sum(qObs(i,:)>=lowerlim&qObs(i,:) <upperlim);
%             region(i,3) = sum(qObs(i,:)>upperlim);
%         end
%     case 2
%         lowerlim=-L/3;
%         upperlim=L/3;
%         for i = 1:numtimepoints
%             region(i,1) = sum(qObs(i,:)<lowerlim);
%             region(i,2) = sum(qObs(i,:)>=lowerlim&qObs(i,:) <upperlim);
%             region(i,3) = sum(qObs(i,:)>upperlim);
%         end
%     otherwise
%         warning
% end
% h=figure; hold on; 
% for i = 1:size(region,2)
%     plot(frameTime*(1:numtimepoints),region(:,i));
% end
% ylim([0 nMols]);
% if length(D)==1
% title(['redistribution for molecules with D = ' num2str(D)]);
% else
%     title(['redistribution for molecules with mean D = ' num2str(mean(D)) 'w s.d. ' num2str(std(D))]);
% end
% saveas(h,['redistribution_' name],'fig')
% % legend({['all molecules -3 to ' num2str(threshold) ' µm'],['all molecules ' num2str(threshold) ' to 3 µm']});

%% plot profile

inRegion = sum(kymo(:,248:254),2);
outRegion = sum(kymo(:,53:59),2);
figure; plot(inRegion); hold on; plot(outRegion);

end
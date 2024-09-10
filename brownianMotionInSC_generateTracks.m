% Lexy von Diezmann, 2023-2024. Released under the GNU GPL v3.

function [name,qObs] = brownianMotionInSC_generateTracks(options,nMols,D,trajDur,unbindRate)
%% 1 to 3-dimensional simulations of Brownian diffusion with localization error and confinement

% Original code - Lexy von Diezmann - Moerner lab - 2018
% This version, LvD 2023-2024 for SYP-3 / ZHP-3 paper
% 
% Generates trajectories calling from a cell of all options
% Allows generation from a single point or 
% This version adds the ability to sample a vector of D vals - 20231126

% brownianMotionInSC_generateTracks(1,1000,'MSD_analyze_HTP3_LP.mat',600,0)

doPlot = 0;
dispInfo = 0; % verbose for some diagnostics

currentTime = datetime('now','Format','yMMdd-HHmm');


if nargin==0
options = 1; % 1: escape from initial spot. 2: sticky patch.
nMols = 10000; % number of simulated molecules
D = 0.0032; % 0.01176; % isotropic diffusion coefficient in um2/s
trajDur = 7; % duration of each trajectory in seconds
unbindRate = 0; % rate in seconds. set to 0 if no unbinding
end
    
    frameTime = 0.1; % seconds. resolution of qObs.
    dt = 0.05; % seconds - sampling rate of raw simulation
    dispN = 500; % display progress in increments of dispN mols simulated
    nDims = 1; % number of dimensions (1 to 3)
    L = 3; % half-length in um
    r = 0.05; % radius of curvature (can be empty) - um
    er = 15; %20e-3; % isotropic localization error in um %40e-3;

if ischar(D)
    load(D,'trackInfo');
    nameStem = D(1:end-4);
    fitD1 = [trackInfo(:).fitD1]; % using principal component
    numLocs = [trackInfo(:).numLocs]; goodTracks = numLocs>=10;
    % allLengths = [trackInfo(:).length]; goodTracks = allLengths>=10;
    
    % for 20230502 S1R6, numLocs>=10, global loc filter = 0, 19 of 210 tracks have negative D
    goodTracks = goodTracks & fitD1>0;
    D = [];    
    for i = find(goodTracks)
        D = vertcat(D,repmat(fitD1(i),numLocs(i),1));
    end
elseif length(D) == 1
    nameStem = ['D_' num2str(D,2)];
else
    nameStem = ['D_min_' num2str(min(D),2) '_max_' num2str(max(D),2)];
end

name = [nameStem '_N_' num2str(nMols) '_dur_' num2str(trajDur) '_opt_' num2str(options) '_unb_' num2str(unbindRate,3)];

% are walls absorptive or reflective? rescales magnitude of bounce \in[eps,1]
% see e.g. comment by Saxton 1995 on reflectivity at viscous bilayer
bounceMag = 3e-3; % basically absorptive but avoids getting completely stuck

switch options
    case 1
        diameter = 0.3; % 300 nm
        edges = [-L -L/3]; % a block of the chromosome
        center = mean(edges); % center position within the block
        window = center+ [-diameter/2, diameter/2];

        % D = 0.0056; %0.022; % SYP-3 sure why not
        startPos = unifrnd(window(1),window(2),1,nMols); % 1/3rd along, d = 0.3
        stickyZone = [10,11];
        % name = 'panel1';
        stickTime = nan(nMols,1);
    case 2

        startPos = unifrnd(-L,L,1,nMols);
        % D = 0.0056; %0.022;
        diameter = 0.15; % 150 nm
        edges = [L/3 L]; % a block of the chromosome
        center = mean(edges); % center position within the block
        stickyZone = center+ [-diameter/2, diameter/2];
        stickTime = inf(nMols,1);
        % name = 'panel2';
    otherwise
% startPos = [-0.9*L*ones(1,nMols);zeros(2,nMols)];
% startPos = [unifrnd(-(L*0.9),(-L*0.8),1,nMols);zeros(2,nMols)];
% % startPos = [unifrnd(-(L),(L),1,nMols);zeros(2,nMols)];
% startPos = startPos(1:nDims,:);
end

maxLim = 6; % axis limit for plotting cdfs

numDs = length(D);
numQ = round(trajDur/dt);

%% Run simulation

% Trajectory of mobile molecles

    % D = 0.001; % TESTING
    if length(D)==1
        dq = vertcat(...
                startPos,...
                sqrt(2*D*dt)*randn(numQ-1,nMols) );
    else
        numDs = length(D);
        scrambledDs = D(randi([1 numDs],1,nMols))';
        dq = vertcat(...
                startPos,...
                sqrt(2*dt)*randn(numQ-1,nMols).*sqrt(scrambledDs) );
    end

    %% sum up to get total q
    q = cumsum(dq,1);
    clear dq

    % add jump points as NaNs to be resolved next section

    qJumpProb = unbindRate*dt; % rate in unit steps
    q(rand(size(q))<qJumpProb)= nan;


    %% introduce constraints
    % case 1: walls in 1D + switching
    if ~isempty(L) && (isempty(r) || nDims == 1)
        constrained = false(1,nMols);
        for mol = 1:nMols
            
            if mol==dispN*round(mol/dispN)
                disp(mol)
            end
            while ~constrained(mol)
                % k1: introduce constraint and find first 'illegal' movement
                % k2: reposition localization if flagged as nan at time t
                molPos = q(:,mol);
                k1 = find(molPos.^2>L^2,1);
                k2 = find(isnan([0; molPos(2:end)]),1); % do not move on first step
                if ~isempty(k1) & ~isempty(k2)
                    k = min(k1,k2);
                else
                    k = [k1 k2];
                end
                
                if isempty(k) % if (now) fully within constraints
                    constrained(mol)=true;
                    continue
                elseif k==k1
                    % recreate trajectory reflecting 'illegal' step
                    k_step = q(k,mol);
                    a=sign(k_step); %dot(k_step,l)./(norm(l)*norm(k_step));
                    rmag = abs(q(k,mol))-L; % magnitude to reflect
%                     rmag = rmag*bounceMag; % user-scaling reflects wall viscosity, etc.
                    qr = -a*(1+bounceMag)*rmag;
                    q(k:end,mol) = q(k:end,mol)+repmat(qr,numQ-k+1,1);
                elseif k==k2
                    % assign target, then add difference vector
                    newPos = unifrnd(-L+eps,L-eps,1);
                    q(k,mol) = newPos;
                    qr = newPos-q(k-1,mol);
                    q(k+1:end,mol) = q(k+1:end,mol)+repmat(qr,numQ-k,1);
                end
            end
        end
        clear constrained
        % useful plots:
        % plot(q(:,1,mol));
        % plot(q(:,1,:));
    end
    % case 2: sticky area
    for mol = 1:nMols
        molQ = q(:,mol);
        k = find(molQ>stickyZone(1)&molQ<stickyZone(2),1); % find first case of sticking
        q(k:end,mol) = unifrnd(stickyZone(1),stickyZone(2),1);% q(k,mol);
        if ~isempty(k)
            stickTime(mol) = k*dt; % in units of seconds
        end
    end


    %% Reshape those into equal length frames
    % trajectory position / within each frame / dimension / realization
    qFr = (reshape(q,uint16(frameTime/dt),[],nMols,nDims));
    
    % Take the mean position of each frame to get the measured positions and add noise to each frame
    % Each row is a frame, each col is a molecule
    qObs = reshape(mean(qFr,1),[],nMols,nDims) + er*randn(trajDur/frameTime,nMols,nDims);
    clear qFr

    %% save and exit
% allQs{setNum} = qObs; %TEMP
% allStickTime(setNum,:) = stickTime;
% 
% simInfo = struct();
% simInfo.D = scrambledDs;
% for i = 1
% % simInfo.Q = qObs()  % TODONEXT
% simInfo.stickTime = stickTime;
% 
% save([name_sims .mat],'simInfo');
% 
% % save and transition to plotting
% clearvars -except name options qObs allQs startPos maxFromStart distFromStart D r L nDims trajDur nRuns nMols frameTime dt Px er maxLim xyObsHist xyPerfHist bounceMag numDs switchTypes settingsCells setNum nameStem doPlot dispInfo stickTime allStickTime
% save([name '.mat']); % redo to save specific variables

% clearvars -except allQs settingsCells setNum nameStem doPlot dispInfo allStickTime
save([name '.mat'],'qObs','frameTime','options','nMols','D','trajDur','unbindRate','stickTime','nameStem');


end

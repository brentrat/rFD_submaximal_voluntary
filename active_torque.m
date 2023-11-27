%% Analysis of individual participants from rFD experiment (DFG experiment #1)
% Run passive_torque.m before this script!
% Initialise workspace
close all; clear

%% Required functions
% Wfilt Wcu rmsDC Wvel vline hline plot_areaerrorbar cbrewer2 colorbrewer.mat

%% Inputs
% Note: 0 = no; 1 = yes
pname = 'C:\Users\Brent\Desktop\Research\DFG\MVC_TA\16\';
num = str2double(regexp(pname,'\d*','Match')); % subject number
fasData = 1;    % analyse fascicle data
zz = 1;         % muscle compartment to analyse: 1=superficial; 2=deep
check = 0;      % check individual trials before continuing
workCheck = 0;  % check work calculation
dur = 12.5;     % trial duration in sec
tqCU = 20;      % cut-off frequency for low-pass filter of torque
angCU = 6;      % cut-off frequency for low-pass filter of angle
avg_plot = 0;   % create plots of mean time traces
zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % returns zero-crossing indices

%% Select files to analyse
if fasData == 0
    cd(pname);
    if ~exist('files','var') || ischar(files) || ~iscell(files) && files == 0
        clear files
        [files(:,1), pname] = uigetfile('*.mat','MultiSelect','On');
    end
elseif fasData == 1
    load([pname 'Tracked\aFilesUSE.mat']);
    files = fAVG;
    clear fAVG
    load([pname 'Tracked\analyseds.mat'],'FLref');
    FLrest = mean(FLref(:,end));
    clear FLref
end

%% Catch if one file is selected
if ~iscell(files)
    error('Error: Select files to analyse.')
end

%% Load necessary files to calculate active torque and active force
% Load passive torque-angle fit:
% P refers to steady state torque;
% P1 to torque during passive lengthening;
% P2 to torque during passive shortening
cd(pname)
load('TQangFit.mat','P','P1','P2');
% Load TA moment arm-angle fit derived from Maganaris et al. 1999 data
load('C:\Users\brent\Documents\MATLAB\GitHub\DFG\maganaris1999_new.mat','MArest');
MA = MArest; clear MArest
% Load EMG data from MVCs for normalisation
load('C:\Users\Brent\Desktop\Research\DFG\MVC_TA\mvcNORM_10Hz.mat')
% Load participant height in cm
load('C:\Users\brent\Documents\MATLAB\GitHub\DFG\height.mat','ht');

if fasData > 0
    cd('Tracked') % trials with fascicle data need to be combined first
end

%% Initialise separate counters for dynamic and fixed-end reference trials
cDyn = 0;
cRef = 0;

%% Loop through files
for jj = 1:length(files)
    if isfile(['T',files{jj}])
        load(['T',files{jj}]) % load combined fascicle data
    else
        load(files{jj}) % load Spike2 data only
    end

    %% Define analysis windows
    % 3 s to ramp, 2 s ramp (5 s), 1 s wait (6 s), 2.5 s hold (8.5 s),
    % 2.5 s hold (11 s)
    s1 = 1/Torque.interval;         % 1 sec
    s2 = s1+s1;                     % 2 sec
    s3 = s2+s1;                     % 3 sec
    h0 = s3+s2+s1+s1-500;           % 6.75 sec
    % Hold phase 1: 8-9 s
    h1 = h0+s1;                     % 7.75 sec
    h2 = h1+s1-1;                   % 8.7495 sec
    % Hold phase 2: 9-10 s
    h3 = h2+1;                      % 8.75 sec
    h4 = h3+s1-1;                   % 9.7495 sec
    % Hold phase 3: 10-11 s
    h5 = h4+1;                      % 9.75 sec
    h6 = h5+s1-1;                   % 10.7495 sec

    %% Align trial starting times
    if isempty(Keyboard.times) % in the case of no keypress (very rare)
        rmsWIN = 125; % length of RMS window in ms
        emgx = rmsDC(TA.values,rmsWIN,1/TA.interval); % determine EMG rise
        [emgMax,emgMaxI] = max(emgx);
        [emgMin,emgMinI] = min(emgx(1:emgMaxI));
        startI = find(emgx(1:emgMaxI)<emgMax/20,1,'last');
        emgMean = mean(emgx(startI-1/Torque.interval*2:startI-1/Torque.interval));
        emgStd = std(emgx(startI-1/Torque.interval*2:startI-1/Torque.interval));
        emgThr = emgMean+emgStd*3; % define EMG start based on 3 std
        startx = find(emgx(1:startI)<=emgThr,1,'last');
        start = startx-1/Torque.interval*3; % crop trial 3 s before EMG rise
    else
        [~,start] = min(abs(Torque.times-Keyboard.times));
    end

    %% Synchronise data
    last = 1/Torque.interval*dur; % make all trials have a common end
    time(:,1) = Torque.times(start:end); % make all trials have a common start

    %% Filter torque and angle data with a dual-pass 2nd-order filter
    tqx(:,1) = Wfilt(Torque.values,tqCU,'low',1/Torque.interval);
    tq = interp1(Torque.times,tqx,time,'linear','extrap');
    angx(:,1) = Wfilt(Angle.values,angCU,'low',1/Angle.interval);
    ang = interp1(Angle.times,angx,time,'linear','extrap');

    %% Check filters
    if check == 1 && jj == 1
        hF = figure(1); 
        tqRaw = interp1(Torque.times,Torque.values,time,'linear','extrap');
        angRaw = interp1(Angle.times,Angle.values,time,'linear','extrap');
        figure(1); hold on; plot(time,tqRaw); plot(time,tq);
        figure(2); hold on; plot(time,angRaw); plot(time,ang);
        FFT(tqRaw,1/Torque.interval);
        waitfor(hF);
    end

    %% Filter the rectified EMG signal with a dual-pass 2nd-order filter
    emg(:,1) = interp1(TA.times,TA.values,time,'linear','extrap');
    % Bandpass software filter the raw EMG signal (10-500 Hz, 4th order)
    [b,a] = butter(4,[10 500]/(1/TA.interval/2));
    emg = filtfilt(b,a,emg);
    emgAbs = abs(emg-mean(emg)); % DC offset
    emgRMS = Wfilt(emgAbs,10,'low',1/TA.interval); % RMS was used initially
    % but this obviously caused timeshift, hence emgRMS = emgAbsFilt

    %% Logical indices for lengthening or shortening trials
    k = strfind(files{jj},'len');
    kk = strfind(files{jj},'sho');

    %% Loop through dynamic trials
    if ~isempty(k) || ~isempty(kk)

        cDyn = cDyn+1; % increment counter, which started at 0
        fDyn{cDyn,1} = files{jj}; % store filename

        %% Calculate active torque
        if ~isempty(k) % lengthening trial
            typeDyn(cDyn,1) = 1;
            tqPassive = polyval(P,ang); %P1 to fit lengthening correctly
            tqActive = tq-tqPassive;
        else
            typeDyn(cDyn,1) = 2; % shortening trial
            tqPassive = polyval(P,ang); %P2 to fit shortening correctly
            tqActive = tq-tqPassive;
        end

        %% Store if there was a preload or not
        if contains(files{jj},'p')
            pre(cDyn,1) = 1;
        else
            pre(cDyn,1) = 0;
        end

        %% Calculate start and end of rotation based on angular velocity
        angVel = Wvel(ang,time);
        % Find zero-crossing indices
        idx = zci(angVel);
        % Find biggest time difference in angular velocity
        [~,maxDiff] = max(diff(idx));
        % First zero-crossing before biggest time diff = start of rotation
        angStart = idx(maxDiff);
        % Next zero-crossing = end of rotation
        angEnd = idx(maxDiff+1);

        %% Determine mean velocity when velocity was constant at 5 deg/s
        if contains(files{jj},'len') % lengthening trials = + velocity

            velFirst = find(angVel>5,1,'first');
            velLast = find(angVel>5,1,'last');
            angVelMean(cDyn,1) = mean(angVel(velFirst:velLast));

        elseif contains(files{jj},'sho') % shortening trials = - velocity

            velFirst = find(angVel<-5,1,'first');
            velLast = find(angVel<-5,1,'last');
            angVelMean(cDyn,1) = mean(angVel(velFirst:velLast));

        end

        %% Store filtered torque, angle, and EMG data over a common time
        tqDyn(cDyn,1:last) = zeros(1,last);
        angDyn(cDyn,1:last) = zeros(1,last);
        emgDyn(cDyn,1:last) = zeros(1,last);

        if length(tqActive) > last

            tqDyn(cDyn,1:last) = tqActive(1:last);
            angDyn(cDyn,1:last) = ang(1:last);
            emgDyn(cDyn,1:last) = emgRMS(1:last);

        else

            tqDyn(cDyn,1:length(time)) = tqActive;
            angDyn(cDyn,1:length(time)) = ang;
            emgDyn(cDyn,1:length(time)) = emgRMS;

        end

        %% Calculate mean torque, angle, and EMG data over the steady state
        hold1(cDyn,1) = mean(tqActive(h1:h2));
        hold1(cDyn,2) = mean(emgRMS(h1:h2));
        hold1(cDyn,3) = mean(ang(h1:h2));
        hold2(cDyn,1) = mean(tqActive(h3:h4));
        hold2(cDyn,2) = mean(emgRMS(h3:h4));
        hold2(cDyn,3) = mean(ang(h3:h4));
        hold3(cDyn,1) = mean(tqActive(h5:h6));
        hold3(cDyn,2) = mean(emgRMS(h5:h6));
        hold3(cDyn,3) = mean(ang(h5:h6));

        %% Calculate the rotation amplitude and initial angle
        rot(cDyn,1) = ang(angEnd)-ang(angStart);
        rot_start(cDyn,1) = ang(angStart);

        %% Analyse fascicle data
        tqDynCh = tqActive; % store active torque as a new variable
        tqStart = 5501; % start time for work calculation

        if exist('TVDdata','var')

            % Load tracked fascicle data
            FLx(:,1) = Fdat.Region(zz).FL;
            FAx(:,1) = rad2deg(Fdat.Region(zz).PEN);

            % Upsample fascicle data rate to 2000 Hz
            FL = interp1(timeUS,FLx,time,'linear','extrap');
            FA = interp1(timeUS,FAx,time,'linear','extrap');

            if exist('FLrest','var')
                FL = FL-(FL(end-55)-FLrest); % Remove fascicle drift
            end

            %% Store fascicle data over a common time
            FLdyn(cDyn,1:last) = zeros(1,last);
            FAdyn(cDyn,1:last) = zeros(1,last);

            if length(FL) > last

                FLdyn(cDyn,1:last) = FL(1:last);
                FAdyn(cDyn,1:last) = FA(1:last);

            else

                FLdyn(cDyn,1:length(time)) = FL;
                FAdyn(cDyn,1:length(time)) = FA;

            end

            %% Calculate net fascicle work
            % Determine fascicle displacement relative to common start time
            FLdynCh = (FL(tqStart:end)-FL(tqStart))/1000;

            % Determine literature-based moment arms based on crank-arm
            % angle (0 deg = footplate perpendicular to tibia)
            MAdyn = polyval(MA,ang(tqStart:end))/100; % cm to m
            MArot = polyval(MA,ang(angStart:angEnd))/100; % cm to m

            % Determine TA tendon force
            FdynLon = tqDynCh(tqStart:end)./MAdyn*0.5;
            FdynLonRot = tqDynCh(angStart:angEnd)./MArot;

            % Determine TA fascicle force
            Fdyn = FdynLon./cosd(FA(tqStart:end));

            % Determine TA fascicle force over whole trial & store
            Fdyn3 = vertcat(tqDyn(1:tqStart-1)',Fdyn);

            ForceDyn(cDyn,1:last) = zeros(1,last);

            if length(Fdyn3) > last

                ForceDyn(cDyn,1:last) = Fdyn3(1:last);

            else

                ForceDyn(cDyn,1:length(time)) = Fdyn3;

            end

            %% Calculate mean TA fascicle force over the steady state
            hold1(cDyn,6) = mean(Fdyn3(h1:h2));
            hold2(cDyn,6) = mean(Fdyn3(h3:h4));
            hold3(cDyn,6) = mean(Fdyn3(h5:h6));

            %% Calculate TA fascicle work from common start to end time
            use = h4-tqStart; % end time

            % Use trapzoidal integration to calculate net fascicle work
            WorkDyn{cDyn}(:,1:2) = [FLdynCh(1:use),Fdyn(1:use)];
            [WDyn(cDyn,1) WDyn(cDyn,2) WDyn(cDyn,3)] = work_calc(-FLdynCh(1:use),Fdyn(1:use));

            %% Calculate TA MTU work from common start to end time
            % Determine TA MTU lengths
            MTU = 0.715+(-0.00130*gnegate(ang(tqStart:end)));

            % Scale MTU length based on participant height
            MTU = MTU*((ht(num)*0.25)/100);

            % Calculate MTU displacement
            MTU = MTU-MTU(1,1);

            % Use trapzoidal integration to calculate net MTU work
            [WDyn(cDyn,4) WDyn(cDyn,5) WDyn(cDyn,6)] = work_calc(-MTU(1:use),FdynLon(1:use));

            %% Calculate net ankle joint work from common start to end time
            % Calculate angular displacement
            angRad = deg2rad(ang);
            angDis = angRad(tqStart:h4)-angRad(tqStart);

            % Use trapzoidal integration to calculate net joint work
            [WDyn(cDyn,7) WDyn(cDyn,8) WDyn(cDyn,9)] = work_calc(-angDis(1:use),tqDynCh(tqStart:h4));

            %% Check work calculation by plotting
            if workCheck

                h = figure(100);
                set(h,'Name',files{jj})

                % Fascicle work
                subplot(211); plot(tqDynCh); hold on; plot(tqStart,tqDynCh(tqStart),'ro'); plot(use+tqStart-1,tqDynCh(use+tqStart-1),'ro');

                % Angular rotation
                subplot(212); plot(ang); hold on; plot(angStart,ang(angStart),'bo'); plot(angEnd,ang(angEnd),'bo');

                % Catch errors in the trial
                if ~contains(files{jj},'p') && angEnd < 7500

                    error('Rotation occurred at the wrong time.')

                elseif contains(files{jj},'p') && angEnd < 12500
                    error('Rotation occurred at the wrong time.')

                end

                waitfor(h);

            end

            %% Calculate fascicle stretch during rotation
            % Use the actual fascicle length trace to calculate stretch
            [~,timeUSstart] = min(abs(timeUS-time(start)));
            FLold = FLx;
            FLx = FLx(timeUSstart:end);

            % Downsample crank-arm angle data to US frame rate
            angUS = interp1(time,ang,timeUS(timeUSstart:end),'linear');

            % Determine angular velocity
            angVelUS = Wvel(angUS,timeUS);

            % Find zero-crossing points
            idx = zci(angVelUS);

            % Determine maximum time difference in zero-crossing points
            [~,maxDiff] = max(diff(idx));

            % Determine start and end of rotation
            angStartUS = idx(maxDiff);
            angEndUS = idx(maxDiff+1);

            % Calculate local and net fascicle stretch
            fasStretch = fas_stretch(FLx(angStartUS:angEndUS));

            %% Calculate fascicle shortening
            % Find start of last steady state
            [~,FLstartShoX] = min(abs(time(h4)-timeUS));

            % Find start of fascicle shortening
            FV = Wvel(FLold(1:FLstartShoX),time(1:FLstartShoX));

            % Find the maximum shortening velocity 
            FVy = min(FV);
            FVx = find(FV<=FVy*0.1,1,'first');

            % Determine the zero-crossing points until max. sho. vel.
            idx = zci(FV(1:FVx));

            % The last zero-crossing point before max. sho. vel. should
            % roughly correspond to the instant of fascicle shortening
            FLstartSho = idx(end-1);

            % Store fascicle length at instant of fascicle shortening
            FLmean = FLold(FLstartSho);

            % Calculate fascicle shortening
            fasSho = min(FLold(1:FLstartShoX))-FLmean;

            %% Store fascicle stretch data
            % L = lengthening
            LDyn(cDyn,1) = fasStretch(1,1); % local fascicle stretch
            LDyn(cDyn,2) = fasStretch(1,2); % cumulative fascicle stretch
            LDyn(cDyn,3) = fasStretch(1,3); % net fascicle stretch

            %% Store fascicle shortening data
            % H = (steady-state) hold phase
            HDyn(cDyn,1) = fasSho; %max. fascicle shortening

            %% Calculate mean TA fascicle length over the steady state
            HDyn(cDyn,2) = mean(FL(h1:h2)); % mean FL
            HDyn(cDyn,3) = mean(FL(h3:h4));
            HDyn(cDyn,4) = mean(FL(h5:h6));

            %% Store fascicle length at contraction onset
            HDyn(cDyn,5) = FLmean; % resting FL

            %% Store mean fascicle length after contraction (1 s)
            HDyn(cDyn,6) = mean(FLx(end-round(1/median(diff(timeUS)),0):end)); % end FL

            %% Store minimum fascicle length during the contraction
            HDyn(cDyn,7) = min(FLold(1:FLstartShoX)); % min FL

            %% Store net fascicle length change from onset to steady state
            HDyn(cDyn,8) = FL(h4)-FLmean; % FL change from rest to start of steady state

            %% Store TA EMG data
            % E = EMG
            EDyn(cDyn,1) = trapz(emgRMS(tqStart:h4)); % EMG area until steady state
            EDyn(cDyn,2) = mean(emgRMS(tqStart:h4)); % mean EMG until steady state
            EDyn(cDyn,3) = mean(emgRMS(tqStart:angEnd)); % mean EMG until end of rot.
            EDyn(cDyn,4) = min(emgRMS(tqStart:angEnd)); % min. EMG until end of rot.
            EDyn(cDyn,5) = max(emgRMS(tqStart:angEnd)); % max. EMG until end of rot.
            EDyn(cDyn,6) = mean(emgRMS(angStart:angEnd)); % mean EMG during rot.

            %% Store TA fascicle force / velocity data
            % F = force
            FDyn(cDyn,1) = mean(Fdyn(1:use)); % mean force until steady state
            FDyn(cDyn,2) = mean(FV(angStartUS:angEndUS));  % mean fasicle velocity during the rotation
            FDyn(cDyn,3) = Fdyn3(angEnd); % fascicle force at the end of the rotation
            FDyn(cDyn,4) = Fdyn3(angStart); % fascicle force at the start of the rotation
            FDyn(cDyn,5) = mean(Fdyn3(tqStart:angEnd)); % mean force until end of rotation
            FDyn(cDyn,6) = mean(Fdyn3(angStart:angEnd)); % mean force during rotation
            FDyn(cDyn,7) = round(mean(FdynLon(1:use)),0); % mean MTU force until steady state
            FDyn(cDyn,8) = round(mean(FdynLonRot),0); % mean MTU force during rotation
            FDyn(cDyn,9) = round((max(MTU(1:use))-min(MTU(1:use)))*1000,1); % MTU length change

        end

    else
        %% Loop through reference trials
        cRef = cRef+1; % increment counter, which started at 0
        fRef{cRef,1} = files{jj}; % store filenames

        %% Calculate active torque
        tqPassive = polyval(P,ang);
        tqActive = tq-tqPassive;

        %% Store filtered torque, angle, and EMG data over a common time
        tqRef(cRef,1:last) = zeros(1,last);
        angRef(cRef,1:last) = zeros(1,last);
        emgRef(cRef,1:last) = zeros(1,last);

        if length(tqActive) > last

            tqRef(cRef,1:last) = tqActive(1:last);
            angRef(cRef,1:last) = ang(1:last);
            emgRef(cRef,1:last) = emgRMS(1:last);

        else

            tqRef(cRef,1:length(time)) = tqActive;
            angRef(cRef,1:length(time)) = ang;
            emgRef(cRef,1:length(time)) = emgRMS;

        end

        %% Calculate mean torque, angle, and EMG data over the steady state
        hold1r(cRef,1) = mean(tqActive(h1:h2));
        hold1r(cRef,2) = mean(emgRMS(h1:h2));
        hold1r(cRef,3) = mean(ang(h1:h2));
        hold2r(cRef,1) = mean(tqActive(h3:h4));
        hold2r(cRef,2) = mean(emgRMS(h3:h4));
        hold2r(cRef,3) = mean(ang(h3:h4));
        hold3r(cRef,1) = mean(tqActive(h5:h6));
        hold3r(cRef,2) = mean(emgRMS(h5:h6));
        hold3r(cRef,3) = mean(ang(h5:h6));

        %% Analyse fascicle data
        tqRefCh = tqActive; % store active torque as a new variable
        tqStart = 5501; % start time for work calculation

        if exist('TVDdata','var')

            % Load tracked fascicle data
            FLx(:,1) = Fdat.Region(zz).FL;
            FAx(:,1) = rad2deg(Fdat.Region(zz).PEN);

            % Upsample fascicle data rate to 2000 Hz
            FL = interp1(timeUS,FLx,time,'linear');
            FA = interp1(timeUS,FAx,time,'linear');

            %% Store fascicle data over a common time
            FLref(cRef,:) = zeros(1,last);
            FAref(cRef,:) = zeros(1,last);

            if length(FL) > last

                FLref(cRef,1:last) = FL(1:last);
                FAref(cRef,1:last) = FA(1:last);

            else

                FLref(cRef,1:length(time)) = FL;
                FAref(cRef,1:length(time)) = FA;

            end

            %% Calculate net fascicle work
            % Determine fascicle displacement relative to common start time
            FLrefCh = (FL(tqStart:end)-FL(tqStart))/1000;

            % Determine literature-based moment arms based on crank-arm
            % angle (0 deg = footplate perpendicular to tibia)
            MAref = polyval(MA,ang(tqStart:end))/100; % cm to m

            % Determine TA tendon force
            FrefLon = tqRefCh(tqStart:end)./MAref*0.5;

            % Determine TA fascicle force
            Fref = FrefLon./cosd(FA(tqStart:end));

            % Determine TA fascicle force over whole trial & store
            Fref3 = vertcat(tqRef(1:tqStart-1)',Fref);

            ForceRef(cRef,1:last) = zeros(1,last);

            if length(Fref3) > last

                ForceRef(cRef,1:last) = Fref3(1:last);

            else

                ForceRef(cRef,1:length(time)) = Fref3;

            end

            %% Calculate mean TA fascicle force over the steady state
            hold1r(cRef,6) = mean(Fref3(h1:h2));
            hold2r(cRef,6) = mean(Fref3(h3:h4));
            hold3r(cRef,6) = mean(Fref3(h5:h6));

            %% Calculate TA fascicle work from common start to end time
            use = h4-tqStart; % end time

            % Use trapzoidal integration to calculate net fascicle work
            WorkRef{cRef}(:,1:2) = [FLrefCh(1:use),Fref(1:use)];
            [WRef(cRef,1) WRef(cRef,2) WRef(cRef,3)] = work_calc(-FLrefCh(1:use),Fref(1:use));

            %% Calculate TA MTU work from common start to end time
            % Determine TA MTU lengths
            MTU = 0.715+(-0.00130*gnegate(ang(tqStart:end)));

            % Scale MTU length based on participant height
            MTU = MTU*((ht(num)*0.25)/100);

            % Calculate MTU displacement
            MTU = MTU-MTU(1,1);

            % Use trapzoidal integration to calculate net MTU work
            [WRef(cRef,4) WRef(cRef,5) WRef(cRef,6)] = work_calc(-MTU(1:use),FrefLon(1:use));
            
            %% Calculate net ankle joint work from common start to end time
            % Calculate angular displacement
            angDis = ang(tqStart:h4)-ang(tqStart);

            % Use trapzoidal integration to calculate net joint work
            [WRef(cRef,7) WRef(cRef,8) WRef(cRef,9)] = work_calc(-angDis(1:use),tqRefCh(tqStart:h4));

            %% Check work calculation by plotting
            if workCheck

                h = figure(100);
                set(h,'Name',files{jj})

                % Fascicle work
                plot(tqRefCh); hold on; plot(tqStart,tqRefCh(tqStart),'ro'); plot(use+tqStart-1,tqRefCh(use+tqStart-1),'ro');
                waitfor(h);

            end

            %% Calculate fascicle shortening
            FLold = FLx;

            % Find start of last steady state
            [~,FLstartShoX] = min(abs(time(h4)-timeUS));

            % Find start of fascicle shortening
            FV = Wvel(FLold(1:FLstartShoX),time(1:FLstartShoX));

            % Find the maximum shortening velocity
            FVy = min(FV);
            FVx = find(FV<=FVy*0.1,1,'first');

            % Determine the zero-crossing points until max. sho. vel.
            idx = zci(FV(1:FVx));

            % The last zero-crossing point before max. sho. vel. should
            % roughly correspond to the instant of fascicle shortening
            FLstartSho = idx(end-1);

            % Store fascicle length at instant of fascicle shortening
            FLmean = FLold(FLstartSho);

            % Calculate fascicle shortening
            fasSho = min(FLold(1:FLstartShoX))-FLmean;

            %% Calculate fascicle lengthening during ramp
            % Use the actual fascicle length trace to calculate lengthening
            [~,timeUSstart] = min(abs(timeUS-time(start)));
            FLx = FLx(timeUSstart:end);

            % Downsample torque data to US frame rate to ensure min.
            % fascicle velocity is calculated during the contraction
            tqUS = interp1(time,tqActive,timeUS(timeUSstart:end),'linear');

            % Find max torque
            [tqMax,tqMaxX] = max(tqUS);

            % Stop FL trace at 95% of maximum torque
            FLendUS = find(tqUS>=tqMax*0.95,1,'first');

            % Start FL trace at start of fascicle shortening
            % Find zero-crossing points in fascicle velocity signal until
            % maximum torque
            timeUSx = timeUS(timeUSstart:end);
            FVx = Wvel(FLx(1:tqMaxX),timeUSx(1:tqMaxX)); 
            [~,FVshoX] = min(FVx);
            idx = zci(FVx(1:FVshoX));
            idx_end = idx(end);
            FLstartUS = idx(end-1);

            if FLstartUS-FLendUS < 22
                FLstartUS = 1;
            end

            % Check fascicle ramp
            if check
                figure; plot(FLx); vline(FLstartUS); vline(FLendUS)
            end

            % Calculate local and net fascicle stretch
            fasStretch = fas_stretch(FLx(FLstartUS:FLendUS));

            %% Store fascicle stretch data
            % L = lengthening
            LRef(cRef,1) = fasStretch(1,1); % local fascicle stretch
            LRef(cRef,2) = fasStretch(1,2); % cumulative fascicle stretch
            LRef(cRef,3) = fasStretch(1,3); % net fascicle stretch

            %% Store fascicle shortening data
            % H = (steady-state) hold phase
            HRef(cRef,1) = fasSho; %max. fascicle shortening

            %% Calculate mean TA fascicle length over the steady state
            HRef(cRef,2) = mean(FL(h1:h2)); % mean FL
            HRef(cRef,3) = mean(FL(h3:h4));
            HRef(cRef,4) = mean(FL(h5:h6));

            %% Store fascicle length at contraction onset
            HRef(cRef,5) = FLmean; % resting FL

            %% Store mean fascicle length after contraction (1 s)
            HRef(cRef,6) = mean(FLx(end-round(1/median(diff(timeUS)),0):end)); % end FL

            %% Store minimum fascicle length during the contraction
            HRef(cRef,7) = min(FLold(1:FLstartShoX)); % min FL

            %% Store net fascicle length change from onset to steady state
            HRef(cRef,8) = FL(h4)-FLmean; % FL change from rest to start of steady state

            %% Store TA EMG data
            % E = EMG
            ERef(cRef,1) = trapz(emgRMS(tqStart:h4)); % EMG area until steady state
            ERef(cRef,2) = mean(emgRMS(tqStart:h4)); % mean EMG until steady state

            %% Store TA fascicle force / velocity data
            % F = force
            FRef(cRef,1) = mean(Fref(1:use)); % mean force until steady state
            FRef(cRef,2) = mean(FV(FLstartUS:FLendUS));  % mean fasicle velocity during the ramp

        end
    end

    %% Plot reference trials
    if check && isempty(k) && isempty(kk)

        figure(1); subplot(311); hold on; plot(tqActive);
        vline(h1); vline(h3); vline(h5); vline(h6);
        subplot(312); hold on; plot(emgRMS);
        vline(h1); vline(h3); vline(h5); vline(h6);
        subplot(313); hold on; plot(ang);
        vline(h1); vline(h3); vline(h5); vline(h6);

    elseif check

        %% Plot dynamic trials
        figure(2); subplot(311); hold on; plot(tqActive);
        vline(h1); vline(h3); vline(h5); vline(h6);
        subplot(312); hold on; plot(emgRMS);
        vline(h1); vline(h3); vline(h5); vline(h6);
        subplot(313); hold on; plot(ang);
        vline(h1); vline(h3); vline(h5); vline(h6);

        if exist('angVel','var')
            figure(3); hold on; plot(angVel);
        end

    end

    clearvars -except ForceRef ForceDyn FDyn FRef FLrest fasData zz zci ht avg_plot files pname dur check workCheck tqCU angCU mvc* num P P1 P2 MA angVelMean cDyn fDyn fRef typeDyn tqDyn angDyn emgDyn hold* pre rot rot_start FLdyn FAdyn WorkDyn WDyn LDyn LRef HDyn EDyn ERef cRef tqRef angRef emgRef FLref FAref WorkRef WRef HRef h*
end

%% Normalize active torque and EMG steady-state data to MVC
T = mvcTQ(num,1);
E = mvcEMG(num,1);

if exist('hold1','var')
    hold1(:,4) = round(hold1(:,1)/T*100,2);
    hold2(:,4) = round(hold2(:,1)/T*100,2);
    hold3(:,4) = round(hold3(:,1)/T*100,2);
    hold1(:,5) = round(hold1(:,2)/E*100,2);
    hold2(:,5) = round(hold2(:,2)/E*100,2);
    hold3(:,5) = round(hold3(:,2)/E*100,2);
end

if exist('hold1r','var')
    hold1r(:,4) = round(hold1r(:,1)/T*100,2);
    hold2r(:,4) = round(hold2r(:,1)/T*100,2);
    hold3r(:,4) = round(hold3r(:,1)/T*100,2);
    hold1r(:,5) = round(hold1r(:,2)/E*100,2);
    hold2r(:,5) = round(hold2r(:,2)/E*100,2);
    hold3r(:,5) = round(hold3r(:,2)/E*100,2);
end

%% Check EMG matching between contraction conditions
if exist('hold1r','var') && exist('hold1','var')

    % Find minimum EMG amp. difference between dynamic and reference trials
    for aa = 1:size(hold1r,1)

        EMGdiff = [round((hold1(:,2)-hold1r(aa,2))/hold1r(aa,2)*100,0),...
            round((hold2(:,2)-hold2r(aa,2))/hold2r(aa,2)*100,0),...
            round((hold3(:,2)-hold3r(aa,2))/hold3r(aa,2)*100,0)];

        for bb = 1:size(EMGdiff,1)

            if pre(bb) == 1

                EMGdiffValid = find(abs(EMGdiff(bb,3))<=5);
                EMGvalid(bb,1) = bb;

            else

                EMGdiffValid = find(abs(EMGdiff(bb,2:3))<=5);

                if size(EMGdiffValid,2) == 2
                    EMGvalid(bb,1) = bb;
                end
            end

        end

        if exist('EMGvalid','var')
            EMGvalid(EMGvalid==0) = [];
            EMGcommon{aa,1} = EMGvalid;
        else
            EMGcommon{aa,1} = 0;
        end
        clear EMGdiff* EMGvalid
    end

    %% Check EMG matching between reference trials
    for aa = 1:size(EMGcommon,1)

        if size(EMGcommon{aa,1},1) == 1 && EMGcommon{aa,1} == 0
            error(['Delete reference trial ' num2str(aa) '.'])
        end

    end

end

%% Plot trial comparisons
c = 1;

if avg_plot

    if exist('FLref','var')

        for x = 1:size(tqDyn,1)

            % Change colors with color brewer
            c = c+1;
            if c == 9
                c = 2;
            end

            figure(199+x);
            subplot(411); hold on;
            plot_areaerrorbar(tqRef,1); hold on; plot_areaerrorbar(tqDyn(x,:),c);
            subplot(412); hold on;
            plot_areaerrorbar(emgRef,1); hold on; plot_areaerrorbar(emgDyn(x,:),c);
            subplot(413); hold on;
            plot_areaerrorbar(angRef,1); hold on; plot_areaerrorbar(angDyn(x,:),c);
            subplot(414); hold on;
            plot_areaerrorbar(FLref,1); hold on; plot_areaerrorbar(FLdyn(x,:),c);

            % Add vertical lines to plot
            subplot(411);
            vline(h1); vline(h3); vline(h5); vline(h6);
            subplot(412);
            vline(h1); vline(h3); vline(h5); vline(h6);
            subplot(413);
            vline(h1); vline(h3); vline(h5); vline(h6);
            subplot(414);
            vline(h1); vline(h3); vline(h5); vline(h6);

        end

    else

        for x = 1:size(tqDyn,1)

            % Change colors with color brewer
            c = c+1;
            if c == 9
                c = 2;
            end

            figure(249+x);
            subplot(311); hold on;
            plot_areaerrorbar(tqRefV,1); hold on; plot_areaerrorbar(tqDynV(x,:),c);
            subplot(312); hold on;
            plot_areaerrorbar(emgRefV,1); hold on; plot_areaerrorbar(emgDynV(x,:),c);
            subplot(313); hold on;
            plot_areaerrorbar(angRefV,1); hold on; plot_areaerrorbar(angDynV(x,:),c);

            % Add vertical lines to plot
            subplot(311);
            vline(h1); vline(h3); vline(h5); vline(h6);
            subplot(312);
            vline(h1); vline(h3); vline(h5); vline(h6);
            subplot(313);
            vline(h1); vline(h3); vline(h5); vline(h6);

        end

    end

end

%% Create new variables to store lengthening and shortening data separately
% Lengthening trials
lenTrials = find(typeDyn==1);
hold1L = hold1(lenTrials,:);
hold2L = hold2(lenTrials,:);
hold3L = hold3(lenTrials,:);

% Shortening trials
shoTrials = find(typeDyn==2);
hold1S = hold1(shoTrials,:);
hold2S = hold2(shoTrials,:);
hold3S = hold3(shoTrials,:);

rotMean = round(mean(rot),1);
rotStart = round(mean(rot_start),1);

%% Show if crank-arm kinematics were similar between dynamic trials
if exist('hold2','var')
    rotAmp_rotVel = [rot angVelMean];
end

%% Round data
if exist('HDyn','var')
    % Round and average data
    % Reference data first
    hold1ref = mean(hold1r,1);
    hold2ref = mean(hold2r,1);
    hold3ref = mean(hold3r,1);

    % M = mean
    HRefM = round(mean(HRef),1);
    WRefM = round(mean(WRef),3);

    % Column 1 doesn't need to be normalised
    EDyn(:,2:end) = EDyn(:,2:end)/E*100;
    ERef(:,2:end) = ERef(:,2:end)/E*100;

    % Calculate torque, EMG, and angle differences
    tqDiff = hold3(:,1)-hold3ref(:,1);
    emgDiff = hold3(:,2)-hold3ref(:,2);
    angDiff = hold3(:,3)-hold3ref(:,3);

    % Calcualte fascicle length differences
    fasStart(:,1) = HDyn(:,5)-HRefM(:,5);
    fasEnd(:,1) = HDyn(:,6)-HRefM(:,6);
    fasMin(:,1) = HDyn(:,7)-HRefM(:,7);

else

    tqDiff = hold3(:,1)-hold3ref(:,1);
    emgDiff = hold3(:,2)-hold3ref(:,2);
    angDiff = hold3(:,3)-hold3ref(:,3);

end

if exist('WDyn','var')
    % Max. fascicle shortening difference
    fasSho(:,1) = HDyn(:,1)-HRefM(:,1);

    % Absolute fascicle length difference during hold phases 1, 2, and 3
    fasLength(:,1) = HDyn(:,2)-HRefM(:,2);
    fasLength(:,2) = HDyn(:,3)-HRefM(:,3);
    fasLength(:,3) = HDyn(:,4)-HRefM(:,4);

    % Fascicle work difference
    fasWorkDiff(:,1) = WDyn(:,1)-WRefM(:,1);

    % MTU work difference
    MTUworkDiff(:,1) = WDyn(:,4)-WRefM(:,4);

    % Angular work difference
    angWorkDiff(:,1) = WDyn(:,7)-WRefM(:,7);

    % Maximum local stretch
    fasStrLoc(:,1) = LDyn(:,1);

    % Fascicle lengthening between min FL during rotation and end of rotation
    fasStrRot(:,1) = LDyn(:,3);

    % Store above variables in easy-to-read table
    AfasDiff = table(fDyn,fasStrLoc,fasStrRot,tqDiff, emgDiff, angDiff,fasSho,fasLength,fasWorkDiff,MTUworkDiff,angWorkDiff,fasStart,fasEnd,fasMin);
end

%% Store data in easy-to-read table
% Reference condition first
tableRows(1,1) = {'ref'};

% No rotation
ROT(1,1) = 0;

% Store steady-state torque, force, EMG, angle, and normalized torque and EMG (3 hold phases)
TQ(1,:) = [round(hold1ref(1,1),1) round(hold2ref(1,1),1) round(hold3ref(1,1),1)];
FORCE(1,:) = [round(hold1ref(1,6),0) round(hold2ref(1,6),0) round(hold3ref(1,6),0)];
EMG(1,:) = [hold1ref(1,2) hold2ref(1,2) hold3ref(1,2)];
ANG(1,:) = [round(hold1ref(1,3),1) round(hold2ref(1,3),1) round(hold3ref(1,3),1)];
TQmvc(1,:) = [round(hold1ref(1,4),2) round(hold2ref(1,4),2) round(hold3ref(1,4),2)];
EMGmvc(1,:) = [round(hold1ref(1,5),2) round(hold2ref(1,5),2) round(hold3ref(1,5),2)];

% Store EMG during ramp (called rot as I care about this data more from the
% dynamic conditions)
EMGrot(1,:) = [mean(ERef(:,1)) round(mean(ERef(:,2)),2) 0 0 0 0];

% Fascicle data
if exist('HRef','var')

    % Store fascicle lengths during the 3 hold phases
    FAS(1,:) = [round(mean(HRef(:,2)),1) round(mean(HRef(:,3)),1) round(mean(HRef(:,4)),1)];

    % Store max. fascicle shortening and net fascicle shortening to steady
    % state
    FASsho(1,:) = [round(mean(HRef(:,1)),1) round(mean(HRef(:,8)),1)];

    % Store fascicle, MTU, and net ankle joint work data (net, pos., neg.)
    FASwork(1,:) = [WRefM(1,1) WRefM(1,2) WRefM(1,3)];
    MTUwork(1,:) = [WRefM(1,4) WRefM(1,5) WRefM(1,6)];
    ANGwork(1,:) = [WRefM(1,7) WRefM(1,8) WRefM(1,9)];
    
    % Store fascicle stretch
    FASstr(1,:) = [round(mean(LRef(:,1)),1) round(mean(LRef(:,2)),1) round(mean(LRef(:,3)),1)];

    % Store force during ramp (called rot as I care about this data more from the
    % dynamic conditions)
    FORrot(1,:) = [round(mean(FRef(:,1)),1) round(mean(FRef(:,2)),1) 0 0 0 0 0 0 0];

end

%% Determine which lengthening condition is which
% lenS; lenL=lenM; lenLP
rotL = rot;
fDummy = fDyn;
fDummy(rotL<0) = [];
rotL(rotL<0) = [];

% lenS < 15 deg rotation
lenS = find(rotL<15);

% lenM >= 15 deg rotation
lenLx = find(rotL>=15);

% Preload or not
lenLP = find(contains(fDummy,'lenlp'));
if ~isempty(lenLP)
    lenL(:,1) = setdiff(lenLx,lenLP);
else
    lenL = lenLx;
end

% Store the same data as the reference condition here and throughout
% It is clearly NOT smart to code like this so avoid it in the future!
% lenS condition second
if ~isempty(lenS)

    tableRows(2,1) = {'lens'};
    ROT(2,1) = round(mean(rot(lenS)),1);
    TQ(2,:) = [round(mean(hold1L(lenS,1)),1) round(mean(hold2L(lenS,1)),1) round(mean(hold3L(lenS,1)),1)];
    FORCE(2,:) = [round(mean(hold1L(lenS,6)),0) round(mean(hold2L(lenS,6)),0) round(mean(hold3L(lenS,6)),0)];
    EMG(2,:) = [mean(hold1L(lenS,2)) mean(hold2L(lenS,2)) mean(hold3L(lenS,2))];
    ANG(2,:) = [round(mean(hold1L(lenS,3)),1) round(mean(hold2L(lenS,3)),1) round(mean(hold3L(lenS,3)),1)];
    TQmvc(2,:) = [round(mean(hold1L(lenS,4)),2) round(mean(hold2L(lenS,4)),2) round(mean(hold3L(lenS,4)),2)];
    EMGmvc(2,:) = [round(mean(hold1L(lenS,5)),2) round(mean(hold2L(lenS,5)),2) round(mean(hold3L(lenS,5)),2)];
    EMGrot(2,:) = [mean(EDyn(lenS,1)) round(mean(EDyn(lenS,2)),2) round(mean(EDyn(lenS,3)),2) round(mean(EDyn(lenS,4)),2) round(mean(EDyn(lenS,5)),2) round(mean(EDyn(lenS,6)),2)];

    if exist('HDyn','var')

        FAS(2,:) = [round(mean(HDyn(lenS,2)),1) round(mean(HDyn(lenS,3)),1) round(mean(HDyn(lenS,4)),1)];
        FASsho(2,:) = [round(mean(HDyn(lenS,1)),1) round(mean(HDyn(lenS,8)),1)];
        FASwork(2,:) = [round(mean(WDyn(lenS,1)),3) round(mean(WDyn(lenS,2)),3) round(mean(WDyn(lenS,3)),3)];
        MTUwork(2,:) = [round(mean(WDyn(lenS,4)),3) round(mean(WDyn(lenS,5)),3) round(mean(WDyn(lenS,6)),3)];
        ANGwork(2,:) = [round(mean(WDyn(lenS,7)),3) round(mean(WDyn(lenS,8)),3) round(mean(WDyn(lenS,9)),3)];
        FASstr(2,:) = [round(mean(LDyn(lenS,1)),1) round(mean(LDyn(lenS,2)),1) round(mean(LDyn(lenS,3)),1)];
        FORrot(2,:) = [round(mean(FDyn(lenS,1)),1) round(mean(FDyn(lenS,2)),1) round(mean(FDyn(lenS,3)),1) round(mean(FDyn(lenS,4)),1) round(mean(FDyn(lenS,5)),1) round(mean(FDyn(lenS,6)),1) round(mean(FDyn(lenS,7)),1) round(mean(FDyn(lenS,8)),1) round(mean(FDyn(lenS,9)),1)];

    end

end

%lenL condition third
if ~isempty(lenL)

    tableRows(3,1) = {'lenl'};
    ROT(3,1) = round(mean(rot(lenL)),1);
    TQ(3,:) = [round(mean(hold1L(lenL,1)),1) round(mean(hold2L(lenL,1)),1) round(mean(hold3L(lenL,1)),1)];
    FORCE(3,:) = [round(mean(hold1L(lenL,6)),0) round(mean(hold2L(lenL,6)),0) round(mean(hold3L(lenL,6)),0)];
    EMG(3,:) = [mean(hold1L(lenL,2)) mean(hold2L(lenL,2)) mean(hold3L(lenL,2))];
    ANG(3,:) = [round(mean(hold1L(lenL,3)),1) round(mean(hold2L(lenL,3)),1) round(mean(hold3L(lenL,3)),1)];
    TQmvc(3,:) = [round(mean(hold1L(lenL,4)),2) round(mean(hold2L(lenL,4)),2) round(mean(hold3L(lenL,4)),2)];
    EMGmvc(3,:) = [round(mean(hold1L(lenL,5)),2) round(mean(hold2L(lenL,5)),2) round(mean(hold3L(lenL,5)),2)];
    EMGrot(3,:) = [mean(EDyn(lenL,1)) round(mean(EDyn(lenL,2)),2) round(mean(EDyn(lenL,3)),2) round(mean(EDyn(lenL,4)),2) round(mean(EDyn(lenL,5)),2) round(mean(EDyn(lenL,6)),2)];

    if exist('HDyn','var')

        FAS(3,:) = [round(mean(HDyn(lenL,2)),1) round(mean(HDyn(lenL,3)),1) round(mean(HDyn(lenL,4)),1)];
        FASsho(3,:) = [round(mean(HDyn(lenL,1)),1) round(mean(HDyn(lenL,8)),1)];
        FASwork(3,:) = [round(mean(WDyn(lenL,1)),3) round(mean(WDyn(lenL,2)),3) round(mean(WDyn(lenL,3)),3)];
        MTUwork(3,:) = [round(mean(WDyn(lenL,4)),3) round(mean(WDyn(lenL,5)),3) round(mean(WDyn(lenL,6)),3)];
        ANGwork(3,:) = [round(mean(WDyn(lenL,7)),3) round(mean(WDyn(lenL,8)),3) round(mean(WDyn(lenL,9)),3)];
        FASstr(3,:) = [round(mean(LDyn(lenL,1)),1) round(mean(LDyn(lenL,2)),1) round(mean(LDyn(lenL,3)),1)];
        FORrot(3,:) = [round(mean(FDyn(lenL,1)),1) round(mean(FDyn(lenL,2)),1) round(mean(FDyn(lenL,3)),1) round(mean(FDyn(lenL,4)),1) round(mean(FDyn(lenL,5)),1) round(mean(FDyn(lenL,6)),1) round(mean(FDyn(lenL,7)),1) round(mean(FDyn(lenL,8)),1) round(mean(FDyn(lenL,9)),1)];

    end

end

% lenLP condition fourth
if ~isempty(lenLP)

    tableRows(4,1) = {'lenlp'};
    ROT(4,1) = round(mean(rot(lenLP)),1);
    TQ(4,:) = [round(mean(hold1L(lenLP,1)),1) round(mean(hold2L(lenLP,1)),1) round(mean(hold3L(lenLP,1)),1)];
    FORCE(4,:) = [round(mean(hold1L(lenLP,6)),0) round(mean(hold2L(lenLP,6)),0) round(mean(hold3L(lenLP,6)),0)];
    EMG(4,:) = [mean(hold1L(lenLP,2)) mean(hold2L(lenLP,2)) mean(hold3L(lenLP,2))];
    ANG(4,:) = [round(mean(hold1L(lenLP,3)),1) round(mean(hold2L(lenLP,3)),1) round(mean(hold3L(lenLP,3)),1)];
    TQmvc(4,:) = [round(mean(hold1L(lenLP,4)),2) round(mean(hold2L(lenLP,4)),2) round(mean(hold3L(lenLP,4)),2)];
    EMGmvc(4,:) = [round(mean(hold1L(lenLP,5)),2) round(mean(hold2L(lenLP,5)),2) round(mean(hold3L(lenLP,5)),2)];
    EMGrot(4,:) = [mean(EDyn(lenLP,1)) round(mean(EDyn(lenLP,2)),2) round(mean(EDyn(lenLP,3)),2) round(mean(EDyn(lenLP,4)),2) round(mean(EDyn(lenLP,5)),2) round(mean(EDyn(lenLP,6)),2)];

    if exist('HDyn','var')

        FAS(4,:) = [round(mean(HDyn(lenLP,2)),1) round(mean(HDyn(lenLP,3)),1) round(mean(HDyn(lenLP,4)),1)];
        FASsho(4,:) = [round(mean(HDyn(lenLP,1)),1) round(mean(HDyn(lenLP,8)),1)];
        FASwork(4,:) = [round(mean(WDyn(lenLP,1)),3) round(mean(WDyn(lenLP,2)),3) round(mean(WDyn(lenLP,3)),3)];
        MTUwork(4,:) = [round(mean(WDyn(lenLP,4)),3) round(mean(WDyn(lenLP,5)),3) round(mean(WDyn(lenLP,6)),3)];
        ANGwork(4,:) = [round(mean(WDyn(lenLP,7)),3) round(mean(WDyn(lenLP,8)),3) round(mean(WDyn(lenLP,9)),3)];
        FASstr(4,:) = [round(mean(LDyn(lenLP,1)),1) round(mean(LDyn(lenLP,2)),1) round(mean(LDyn(lenLP,3)),1)];
        FORrot(4,:) = [round(mean(FDyn(lenLP,1)),1) round(mean(FDyn(lenLP,2)),1) round(mean(FDyn(lenLP,3)),1) round(mean(FDyn(lenLP,4)),1) round(mean(FDyn(lenLP,5)),1) round(mean(FDyn(lenLP,6)),1) round(mean(FDyn(lenLP,7)),1) round(mean(FDyn(lenLP,8)),1) round(mean(FDyn(lenLP,9)),1)];

    end

end

%% Determine which shortening condition is which
% shoS; shoSP; shoL; shoLP
v = max(vertcat(lenS,lenL,lenLP));
rotS = rot;
fDummy = fDyn;
fDummy(rotS>0) = [];
rotS(rotS>0) = [];

% shoS < 15 deg rotation
shoSx = find(rotS>-15);

% shoL >= 15 deg rotation
shoLx = find(rotS<=-15);

% Preload or not
shoSP = find(contains(fDummy,'shosp'));
if ~isempty(shoSP)
    shoS(:,1) = setdiff(shoSx,shoSP);
else
    shoS = shoSx;
end

shoLPx = find(contains(fDummy,'sholp'));
if ~isempty(shoLPx)
    shoLxx(:,1) = setdiff(shoLx,shoLPx);
end

% Additional condition with max activation during shortening
shoLPM = find(contains(fDummy,'sholpm'));
if ~isempty(shoLPM)
    shoLP(:,1) = setdiff(shoLPx,shoLPM);
    shoL(:,1) = setdiff(shoLxx,shoLPM);
else
    shoLP = shoLPx;
    if exist('shoLxx','var')
        shoL = shoLxx;
    else
        shoL = shoLx;
    end
end

% shoS condition fifth
if ~isempty(shoS)

    tableRows(5,1) = {'shos'};
    ROT(5,1) = round(mean(rot(v+shoS)),1);
    TQ(5,:) = [round(mean(hold1S(shoS,1)),1) round(mean(hold2S(shoS,1)),1) round(mean(hold3S(shoS,1)),1)];
    FORCE(5,:) = [round(mean(hold1S(shoS,6)),0) round(mean(hold2S(shoS,6)),0) round(mean(hold3S(shoS,6)),0)];
    EMG(5,:) = [mean(hold1S(shoS,2)) mean(hold2S(shoS,2)) mean(hold3S(shoS,2))];
    ANG(5,:) = [round(mean(hold1S(shoS,3)),1) round(mean(hold2S(shoS,3)),1) round(mean(hold3S(shoS,3)),1)];
    TQmvc(5,:) = [round(mean(hold1S(shoS,4)),2) round(mean(hold2S(shoS,4)),2) round(mean(hold3S(shoS,4)),2)];
    EMGmvc(5,:) = [round(mean(hold1S(shoS,5)),2) round(mean(hold2S(shoS,5)),2) round(mean(hold3S(shoS,5)),2)];
    EMGrot(5,:) = [mean(EDyn(v+shoS,1)) round(mean(EDyn(v+shoS,2)),2) round(mean(EDyn(v+shoS,3)),2) round(mean(EDyn(v+shoS,4)),2) round(mean(EDyn(v+shoS,5)),2) round(mean(EDyn(v+shoS,6)),2)];

    if exist('HDyn','var')

        FAS(5,:) = [round(mean(HDyn(v+shoS,2)),1) round(mean(HDyn(v+shoS,3)),1) round(mean(HDyn(v+shoS,4)),1)];
        FASsho(5,:) = [round(mean(HDyn(v+shoS,1)),1) round(mean(HDyn(v+shoS,8)),1)];
        FASwork(5,:) = [round(mean(WDyn(v+shoS,1)),3) round(mean(WDyn(v+shoS,2)),3) round(mean(WDyn(v+shoS,3)),3)];
        MTUwork(5,:) = [round(mean(WDyn(v+shoS,4)),3) round(mean(WDyn(v+shoS,5)),3) round(mean(WDyn(v+shoS,6)),3)];
        ANGwork(5,:) = [round(mean(WDyn(v+shoS,7)),3) round(mean(WDyn(v+shoS,8)),3) round(mean(WDyn(v+shoS,9)),3)];
        FASstr(5,:) = [round(mean(LDyn(v+shoS,1)),1) round(mean(LDyn(v+shoS,2)),1) round(mean(LDyn(v+shoS,3)),1)];
        FORrot(5,:) = [round(mean(FDyn(v+shoS,1)),1) round(mean(FDyn(v+shoS,2)),1) round(mean(FDyn(v+shoS,3)),1) round(mean(FDyn(v+shoS,4)),1) round(mean(FDyn(v+shoS,5)),1) round(mean(FDyn(v+shoS,6)),1) round(mean(FDyn(v+shoS,7)),1) round(mean(FDyn(v+shoS,8)),1) round(mean(FDyn(v+shoS,9)),1)];

    end

end

% shoSP condition sixth
if ~isempty(shoSP)

    tableRows(6,1) = {'shosp'};
    ROT(6,1) = round(mean(rot(v+shoSP)),1);
    TQ(6,:) = [round(mean(hold1S(shoSP,1)),1) round(mean(hold2S(shoSP,1)),1) round(mean(hold3S(shoSP,1)),1)];
    FORCE(6,:) = [round(mean(hold1S(shoSP,6)),0) round(mean(hold2S(shoSP,6)),0) round(mean(hold3S(shoSP,6)),0)];
    EMG(6,:) = [mean(hold1S(shoSP,2)) mean(hold2S(shoSP,2)) mean(hold3S(shoSP,2))];
    ANG(6,:) = [round(mean(hold1S(shoSP,3)),1) round(mean(hold2S(shoSP,3)),1) round(mean(hold3S(shoSP,3)),1)];
    TQmvc(6,:) = [round(mean(hold1S(shoSP,4)),2) round(mean(hold2S(shoSP,4)),2) round(mean(hold3S(shoSP,4)),2)];
    EMGmvc(6,:) = [round(mean(hold1S(shoSP,5)),2) round(mean(hold2S(shoSP,5)),2) round(mean(hold3S(shoSP,5)),2)];
    EMGrot(6,:) = [mean(EDyn(v+shoSP,1)) round(mean(EDyn(v+shoSP,2)),2) round(mean(EDyn(v+shoSP,3)),2) round(mean(EDyn(v+shoSP,4)),2) round(mean(EDyn(v+shoSP,5)),2) round(mean(EDyn(v+shoSP,6)),2)];

    if exist('HDyn','var')

        FAS(6,:) = [round(mean(HDyn(v+shoSP,2)),1) round(mean(HDyn(v+shoSP,3)),1) round(mean(HDyn(v+shoSP,4)),1)];
        FASsho(6,:) = [round(mean(HDyn(v+shoSP,1)),1) round(mean(HDyn(v+shoSP,8)),1)];
        FASwork(6,:) = [round(mean(WDyn(v+shoSP,1)),3) round(mean(WDyn(v+shoSP,2)),3) round(mean(WDyn(v+shoSP,3)),3)];
        MTUwork(6,:) = [round(mean(WDyn(v+shoSP,4)),3) round(mean(WDyn(v+shoSP,5)),3) round(mean(WDyn(v+shoSP,6)),3)];
        ANGwork(6,:) = [round(mean(WDyn(v+shoSP,7)),3) round(mean(WDyn(v+shoSP,8)),3) round(mean(WDyn(v+shoSP,9)),3)];
        FASstr(6,:) = [round(mean(LDyn(v+shoSP,1)),1) round(mean(LDyn(v+shoSP,2)),1) round(mean(LDyn(v+shoSP,3)),1)];
        FORrot(6,:) = [round(mean(FDyn(v+shoSP,1)),1) round(mean(FDyn(v+shoSP,2)),1) round(mean(FDyn(v+shoSP,3)),1) round(mean(FDyn(v+shoSP,4)),1) round(mean(FDyn(v+shoSP,5)),1) round(mean(FDyn(v+shoSP,6)),1) round(mean(FDyn(v+shoSP,7)),1) round(mean(FDyn(v+shoSP,8)),1) round(mean(FDyn(v+shoSP,9)),1)];

    end

end

% shoL condition seventh
if ~isempty(shoL)

    tableRows(7,1) = {'shol'};
    ROT(7,1) = round(mean(rot(v+shoL)),1);
    TQ(7,:) = [round(mean(hold1S(shoL,1)),1) round(mean(hold2S(shoL,1)),1) round(mean(hold3S(shoL,1)),1)];
    FORCE(7,:) = [round(mean(hold1S(shoL,6)),0) round(mean(hold2S(shoL,6)),0) round(mean(hold3S(shoL,6)),0)];
    EMG(7,:) = [mean(hold1S(shoL,2)) mean(hold2S(shoL,2)) mean(hold3S(shoL,2))];
    ANG(7,:) = [round(mean(hold1S(shoL,3)),1) round(mean(hold2S(shoL,3)),1) round(mean(hold3S(shoL,3)),1)];
    TQmvc(7,:) = [round(mean(hold1S(shoL,4)),2) round(mean(hold2S(shoL,4)),2) round(mean(hold3S(shoL,4)),2)];
    EMGmvc(7,:) = [round(mean(hold1S(shoL,5)),2) round(mean(hold2S(shoL,5)),2) round(mean(hold3S(shoL,5)),2)];
    EMGrot(7,:) = [mean(EDyn(v+shoL,1)) round(mean(EDyn(v+shoL,2)),2) round(mean(EDyn(v+shoL,3)),2) round(mean(EDyn(v+shoL,4)),2) round(mean(EDyn(v+shoL,5)),2) round(mean(EDyn(v+shoL,6)),2)];

    if exist('HDyn','var')

        FAS(7,:) = [round(mean(HDyn(v+shoL,2)),1) round(mean(HDyn(v+shoL,3)),1) round(mean(HDyn(v+shoL,4)),1)];
        FASsho(7,:) = [round(mean(HDyn(v+shoL,1)),1) round(mean(HDyn(v+shoL,8)),1)];
        FASwork(7,:) = [round(mean(WDyn(v+shoL,1)),3) round(mean(WDyn(v+shoL,2)),3) round(mean(WDyn(v+shoL,3)),3)];
        MTUwork(7,:) = [round(mean(WDyn(v+shoL,4)),3) round(mean(WDyn(v+shoL,5)),3) round(mean(WDyn(v+shoL,6)),3)];
        ANGwork(7,:) = [round(mean(WDyn(v+shoL,7)),3) round(mean(WDyn(v+shoL,8)),3) round(mean(WDyn(v+shoL,9)),3)];
        FASstr(7,:) = [round(mean(LDyn(v+shoL,1)),1) round(mean(LDyn(v+shoL,2)),1) round(mean(LDyn(v+shoL,3)),1)];
        FORrot(7,:) = [round(mean(FDyn(v+shoL,1)),1) round(mean(FDyn(v+shoL,2)),1) round(mean(FDyn(v+shoL,3)),1) round(mean(FDyn(v+shoL,4)),1) round(mean(FDyn(v+shoL,5)),1) round(mean(FDyn(v+shoL,6)),1) round(mean(FDyn(v+shoL,7)),1) round(mean(FDyn(v+shoL,8)),1) round(mean(FDyn(v+shoL,9)),1)];

    end

end

% shoLP condition eighth
if ~isempty(shoLP)

    tableRows(8,1) = {'sholp'};
    ROT(8,1) = round(mean(rot(v+shoLP)),1);
    TQ(8,:) = [round(mean(hold1S(shoLP,1)),1) round(mean(hold2S(shoLP,1)),1) round(mean(hold3S(shoLP,1)),1)];
    FORCE(8,:) = [round(mean(hold1S(shoLP,6)),0) round(mean(hold2S(shoLP,6)),0) round(mean(hold3S(shoLP,6)),0)];
    EMG(8,:) = [mean(hold1S(shoLP,2)) mean(hold2S(shoLP,2)) mean(hold3S(shoLP,2))];
    ANG(8,:) = [round(mean(hold1S(shoLP,3)),1) round(mean(hold2S(shoLP,3)),1) round(mean(hold3S(shoLP,3)),1)];
    TQmvc(8,:) = [round(mean(hold1S(shoLP,4)),2) round(mean(hold2S(shoLP,4)),2) round(mean(hold3S(shoLP,4)),2)];
    EMGmvc(8,:) = [round(mean(hold1S(shoLP,5)),2) round(mean(hold2S(shoLP,5)),2) round(mean(hold3S(shoLP,5)),2)];
    EMGrot(8,:) = [mean(EDyn(v+shoLP,1)) round(mean(EDyn(v+shoLP,2)),2) round(mean(EDyn(v+shoLP,3)),2) round(mean(EDyn(v+shoLP,4)),2) round(mean(EDyn(v+shoLP,5)),2) round(mean(EDyn(v+shoLP,6)),2)];

    if exist('HDyn','var')

        FAS(8,:) = [round(mean(HDyn(v+shoLP,2)),1) round(mean(HDyn(v+shoLP,3)),1) round(mean(HDyn(v+shoLP,4)),1)];
        FASsho(8,:) = [round(mean(HDyn(v+shoLP,1)),1) round(mean(HDyn(v+shoLP,8)),1)];
        FASwork(8,:) = [round(mean(WDyn(v+shoLP,1)),3) round(mean(WDyn(v+shoLP,2)),3) round(mean(WDyn(v+shoLP,3)),3)];
        MTUwork(8,:) = [round(mean(WDyn(v+shoLP,4)),3) round(mean(WDyn(v+shoLP,5)),3) round(mean(WDyn(v+shoLP,6)),3)];
        ANGwork(8,:) = [round(mean(WDyn(v+shoLP,7)),3) round(mean(WDyn(v+shoLP,8)),3) round(mean(WDyn(v+shoLP,9)),3)];
        FASstr(8,:) = [round(mean(LDyn(v+shoLP,1)),1) round(mean(LDyn(v+shoLP,2)),1) round(mean(LDyn(v+shoLP,3)),1)];
        FORrot(8,:) = [round(mean(FDyn(v+shoLP,1)),1) round(mean(FDyn(v+shoLP,2)),1) round(mean(FDyn(v+shoLP,3)),1) round(mean(FDyn(v+shoLP,4)),1) round(mean(FDyn(v+shoLP,5)),1) round(mean(FDyn(v+shoLP,6)),1) round(mean(FDyn(v+shoLP,7)),1) round(mean(FDyn(v+shoLP,8)),1) round(mean(FDyn(v+shoLP,9)),1)];

    end

end

% shoLPM condition ninth
if ~isempty(shoLPM)

    tableRows(9,1) = {'sholpm'};
    ROT(9,1) = round(mean(rot(v+shoLPM)),1);
    TQ(9,:) = [round(mean(hold1S(shoLPM,1)),1) round(mean(hold2S(shoLPM,1)),1) round(mean(hold3S(shoLPM,1)),1)];
    FORCE(9,:) = [round(mean(hold1S(shoLPM,6)),0) round(mean(hold2S(shoLPM,6)),0) round(mean(hold3S(shoLPM,6)),0)];
    EMG(9,:) = [mean(hold1S(shoLPM,2)) mean(hold2S(shoLPM,2)) mean(hold3S(shoLPM,2))];
    ANG(9,:) = [round(mean(hold1S(shoLPM,3)),1) round(mean(hold2S(shoLPM,3)),1) round(mean(hold3S(shoLPM,3)),1)];
    TQmvc(9,:) = [round(mean(hold1S(shoLPM,4)),2) round(mean(hold2S(shoLPM,4)),2) round(mean(hold3S(shoLPM,4)),2)];
    EMGmvc(9,:) = [round(mean(hold1S(shoLPM,5)),2) round(mean(hold2S(shoLPM,5)),2) round(mean(hold3S(shoLPM,5)),2)];
    EMGrot(9,:) = [mean(EDyn(v+shoLPM,1)) round(mean(EDyn(v+shoLPM,2)),2) round(mean(EDyn(v+shoLPM,3)),2) round(mean(EDyn(v+shoLPM,4)),2) round(mean(EDyn(v+shoLPM,5)),2) round(mean(EDyn(v+shoLPM,5)),2) round(mean(EDyn(v+shoLPM,6)),2)];

    if exist('HDyn','var')

        FAS(9,:) = [round(mean(HDyn(v+shoLPM,2)),1) round(mean(HDyn(v+shoLPM,3)),1) round(mean(HDyn(v+shoLPM,4)),1)];
        FASsho(9,:) = [round(mean(HDyn(v+shoLPM,1)),1) round(mean(HDyn(v+shoLPM,8)),1)];
        FASwork(9,:) = [round(mean(WDyn(v+shoLPM,1)),3) round(mean(WDyn(v+shoLPM,2)),3) round(mean(WDyn(v+shoLPM,3)),3)];
        MTUwork(9,:) = [round(mean(WDyn(v+shoLPM,4)),3) round(mean(WDyn(v+shoLPM,5)),3) round(mean(WDyn(v+shoLPM,6)),3)];
        ANGwork(9,:) = [round(mean(WDyn(v+shoLPM,7)),3) round(mean(WDyn(v+shoLPM,8)),3) round(mean(WDyn(v+shoLPM,9)),3)];
        FASstr(9,:) = [round(mean(LDyn(v+shoLPM,1)),1) round(mean(LDyn(v+shoLPM,2)),1) round(mean(LDyn(v+shoLPM,3)),1)];
        FORrot(9,:) = [round(mean(FDyn(v+shoLPM,1)),1) round(mean(FDyn(v+shoLPM,2)),1) round(mean(FDyn(v+shoLPM,3)),1) round(mean(FDyn(v+shoLPM,4)),1) round(mean(FDyn(v+shoLPM,5)),1) round(mean(FDyn(v+shoLPM,6)),1) round(mean(FDyn(v+shoLPM,7)),1) round(mean(FDyn(v+shoLPM,8)),1) round(mean(FDyn(v+shoLPM,9)),1)];

    end

end

% Create a table with the above data
if exist('HDyn','var')

    A = table(tableRows,ROT,TQ,FORCE,EMG,ANG,TQmvc,EMGmvc,EMGrot,FAS,FASsho,FASwork,MTUwork,ANGwork,FASstr,FORrot);

else

    % If fascicle data not analysed
    A = table(tableRows,ROT,TQ,EMG,EMGrot,ANG,TQmvc,EMGmvc);

end

% Create another table with dynamic conditions only
Adiff = table(tableRows(2:end),round((TQ(2:end,2)-TQ(1,2))/TQ(1,2)*100,2),round((TQ(2:end,3)-TQ(1,3))/TQ(1,3)*100,2),round((EMG(2:end,2)-EMG(1,2))/EMG(1,2)*100,2),round((EMG(2:end,3)-EMG(1,3))/EMG(1,3)*100,2),round((TQmvc(2:end,2)-TQmvc(1,2))/TQmvc(1,2)*100,2),round((TQmvc(2:end,3)-TQmvc(1,3))/TQmvc(1,3)*100,2),round((EMGmvc(2:end,2)-EMGmvc(1,2))/EMGmvc(1,2)*100,2),round((EMGmvc(2:end,3)-EMGmvc(1,3))/EMGmvc(1,3)*100,2));

%% Save data
save('analyzedJP.mat')
%[min(rot_start) max(rot_start)]
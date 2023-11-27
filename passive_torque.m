%% Passive torque analysis to create a passive torque-angle fits
% Initialise workspace without reselecting files of interest
close all
clearvars -except files

%% Required functions
% Wfilt Wcu rmsDC Wvel vline hline

%% Inputs
% Note: 0 = no; 1 = yes
usePas = 0;     % use passive data in steady-state passive torque-angle fit
velDesired = 5; % constant velocity of the passive rotation
df = -30;       % maximum dorsiflexion angle to perform the fit to *-10
pf = 10;        % maximum plantar flexion angle to perform the fit to *50
order = 3;      % order of the polynomial fit
flip = 0;       % make negative angles positive
check = 0;      % check individual trials before continuing
se = 1;         % average start & end values
duration = 1000;% duration for average where 2000 = 1 second
tqCU = 20;      % cut-off frequency for LP filter of torque
angCU = 6;      % cut-off frequency for LP filter of angle
rmsWIN = 125;   % duration of RMS window in ms

%% Select files to analyse
if ~exist('files','var') || ischar(files) || ~iscell(files) && files == 0
    clear files
    [files(:,1), pname] = uigetfile('*.mat','MultiSelect','On');
end

if ~iscell(files) && files == 0
    error('Select files to analyse.')
end

files(contains(files,'TQangFit')) = [];

%% Loop through files
for jj = 1:length(files)

    load(files{jj})

    if jj == 1 && duration ~= 1/Torque.interval
        error('Duration is incorrect.')
    end

    %% Synchronise data
    if isfield(Torque,'times')
        time(:,1) = Torque.times;
    else
        error('Export the time vectors with the data vectors.')
    end

    %% Filter torque and angle data with a dual-pass 2nd-order filter
    tq(:,1) = Wfilt(Torque.values,tqCU,'low',1/Torque.interval);
    angx(:,1) = Wfilt(Angle.values,angCU,'low',1/Angle.interval);
    ang = interp1(Angle.times,angx,time,'linear','extrap');

    %% Calculate RMS amplitudes from the raw EMG signal
    if exist('TA','var')
        emg(:,1) = interp1(TA.times,TA.values,time,'linear','extrap');
        emgRMS = rmsDC(emg,rmsWIN,1,1/TA.interval);

    elseif exist('EMG_TA','var')
        emg(:,1) = interp1(EMG_TA.times,EMG_TA.values,time,'linear','extrap');
        emgRMS = rmsDC(emg,rmsWIN,1,1/EMG_TA.interval);

    elseif exist('TA_EMG','var')
        emg(:,1) = interp1(TA_EMG.times,TA_EMG.values,time,'linear','extrap');
        emgRMS = rmsDC(emg,rmsWIN,1,1/TA_EMG.interval);

    end

    %% Passive rotation trials
    % Determine start and end of rotation where velocity is constant
    % based on angular velocity vector
    if contains(files{jj},'pas')

        angVel = Wvel(ang,time);

        % Passive lengthening trials
        if contains(files{jj},'len')

            if velDesired > 0
                s = find(angVel>velDesired,1,'first');
                e = find(angVel>velDesired,1,'last');
            else
                s = find(angVel<velDesired,1,'first');
                e = find(angVel<velDesired,1,'last');
            end
            col = 1;

            % Passive shortening trials
        elseif contains(files{jj},'sho')

            if  velDesired*-1 < 0
                s = find(angVel<velDesired*-1,1,'first');
                e = find(angVel<velDesired*-1,1,'last');
            else
                s = find(angVel>velDesired*-1,1,'first');
                e = find(angVel>velDesired*-1,1,'last');
            end
            col = 2;

        end

        angStart = s;
        angEnd = e;

        %% Store torque and angle during the rotation
        tqFitPas{col} = tq(angStart:angEnd);
        angFitPas{col} = ang(angStart:angEnd);

        %% Store torque and angle before and after the rotation
        % Note: B = Beginning; E = End
        zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % returns zero-crossing indices
        idx = zci(angVel);
        % Find the last zero-crossing index before the rotation starts
        [~,max_diff] = max(diff(idx));
        % Start and end of rotation
        s = idx(max_diff);
        e = idx(max_diff+1);
        % chPts = findchangepts(angVel,'MaxNumChanges',2,'Statistic','mean');
        s2 = s;
        e1 = e;

        % Determine the steady state before the rotation
        if s < duration
            s1 = 1;
        else
            s1 = s-duration;
        end

        % Store steady-state angle and torque before the rotation
        angFitB(jj,1) = mean(ang(s1:s2));
        tqFitB(jj,1) = mean(tq(s1:s2));

        % Determine the steady state after the rotation
        if e + duration > length(time)
            e2 = length(time);
        else
            e2 = e+duration-1;
        end

        % Store steady-state angle and torque after the rotation
        angFitE(jj,1) = mean(ang(e1:e2));
        tqFitE(jj,1) = mean(tq(e1:e2));

        %% Plot passive rotation trials
        h = figure(jj); set(h,'Name',files{jj});
        subplot(311); plot(tq);
        ylabel('Torque (Nm)')
        vline(angStart,'b'); vline(angEnd,'b'); vline(s1); vline(s2); vline(e1); vline(e2);
        subplot(312); plot(ang);
        ylabel('Angle (deg)')
        vline(angStart,'b'); vline(angEnd,'b'); vline(s1); vline(s2); vline(e1); vline(e2);
        if exist('emgRMS','var')
            subplot(313); plot(emgRMS);
        end
        ylabel('EMG RMS (V)')
        vline(angStart,'b'); vline(angEnd,'b'); vline(s1); vline(s2); vline(e1); vline(e2);

        if check
            waitfor(h);
        end
    else
        %% All other trials
        if flip
            ang = gnegate(ang);
        end

        %% Find minimum EMG value before contraction
        if exist('emgRMS','var')
        [~,emgMaxI] = max(emgRMS);
        [emgMin,emgMinI] = min(emgRMS(1:emgMaxI));

        % Find the passive torque before the contraction
        if ~isempty(Keyboard.times)
            [~,start] = min(abs(time-Keyboard.times(end)));
        else
            % Determine passive torque until minimum EMG value
            start = emgMinI-duration;
        end
        end

        if exist('start','var') && start < 1 % check non-negative 1st index
            start = duration*5; % assuming that contraction starts at 6 s
        elseif ~exist('start','var')
            start = 1;
        end

        % Determine steady state before the contraction starts
        s1 = start;
        s2 = start+duration-1;

        % Change passive torque calculation if maximum EMG > 2*minimum
        if exist('emgRMS','var')
        if max(emgRMS(s1:s2)) > emgMin*2
            s2 = s1;
            s1 = s2-duration;

            if s1<1 % check non-negative 1st index
                s1 = 1;
                s2 = s1+duration-1;
            end

            % Check again (3)
            if max(emgRMS(s1:s2)) > emgMin*2
                s2 = s1;
                s1 = s2-duration;

                if s1<1 % check non-negative 1st index
                    s1 = 1;
                    s2 = s1+duration-1;
                end

                % Check again (4)
                if  max(emgRMS(s1:s2)) > emgMin*2
                    s2 = s1;
                    s1 = s2-duration;

                    if s1<1 % check non-negative 1st index
                        s1 = 1;
                        s2 = s1+duration-1;
                    end

                    % Check again (5) - 5 sec before contraction
                    if max(emgRMS(s1:s2)) > emgMin*2
                        h = figure(jj); set(h,'Name',files{jj});
                        subplot(311); plot(tq);
                        vline(s1); vline(s2);
                        subplot(312); plot(emgRMS);
                        vline(s1); vline(s2); hline(emgMin*2)
                        subplot(313); plot(ang);
                        vline(s1); vline(s2);
                        error('Muscle was activated more than 2x the minimum EMG RMS during the start window.')
                    end
                end
            end
        end
        end

        %% Find minimum EMG value after contraction
        e2 = length(time);
        e1 = e2-duration/4;

        if exist('emgRMS','var')
        % Change passive torque calculation if maximum EMG > 2*minimum
        if se && max(emgRMS(e1:e2)) > emgMin*2
            e2 = e1;
            e1 = e2-duration;

            % Check again (3)
            if max(emgRMS(e1:e2)) > emgMin*2
                e2 = e1;
                e1 = e2-duration;

                % Check again (4)
                if  max(emgRMS(e1:e2)) > emgMin*2
                    h = figure(jj); set(h,'Name',files{jj});
                    subplot(3,1,1); plot(tq);
                    vline(e1); vline(e2);
                    subplot(3,1,2); plot(ang);
                    vline(e1); vline(e2);
                    subplot(3,1,3); plot(emgRMS);
                    vline(e1); vline(e2); hline(emgMin*2)
                    error('Muscle was activated more than 2x the minimum EMG RMS during the start window.')
                end
            end
        end
        end

        % Store steady-state angle and torque before and after the rotation
        angFitB(jj,1) = mean(ang(s1:s2));
        tqFitB(jj,1) = mean(tq(s1:s2));
        angFitE(jj,1) = mean(ang(e1:e2));
        tqFitE(jj,1) = mean(tq(e1:e2));
    end

    if check && ~contains(files{jj},'pas')
        %% Plot all other trials
        h = figure(jj); set(h,'Name',files{jj});
        subplot(3,1,1); plot(tq);
        vline(s1); vline(s2); vline(e1); vline(e2);
        subplot(3,1,2); plot(ang);
        vline(s1); vline(s2); vline(e1); vline(e2);
        subplot(3,1,3); plot(emg);
        vline(s1); vline(s2); vline(e1); vline(e2);

        waitfor(h);
    end

    clearvars -except files pname order flip check se duration tqCU angCU rmsWIN angFit* tqFit* df pf velDesired usePas col
end

%% Create angle vs. torque fits from DF to PF
x(:,1) = df:0.01:pf;
figure(length(files)+1); hold on
xlabel('Plantar flexion angle (deg)','fontweight','bold')
ylabel('Net ankle joint torque (Nm)','fontweight','bold')

%% Passive rotation fits first - lengthening then shortening
if exist('tqFitPas','var')

    if ~isempty(angFitPas{1}) && ~isempty(angFitPas{2})
        P1 = polyfit(angFitPas{1},tqFitPas{1},order);
        y1 = polyval(P1,x);
        plot(x,y1,'r-');
        P2 = polyfit(angFitPas{2},tqFitPas{2},order);
        y2 = polyval(P2,x);
        plot(x,y2,'b-');

    elseif ~isempty(angFitPas{1})
        P1 = polyfit(angFitPas{1},tqFitPas{1},order);
        y1 = polyval(P1,x);
        plot(x,y1,'r-');

    elseif ~isempty(angFitPas{2})
        P2 = polyfit(angFitPas{2},tqFitPas{2},order);
        y2 = polyval(P2,x);
        plot(x,y2,'b-');
    end
end

%% Fit all other trials, including start and end of passive rotation trials
if exist('tqFitB','var')

    %% Fit only start values not including passive rotations
    if se ~= 1 && usePas == 0
        ang = angFitB(:,1);
        tq = tqFitB(:,1);
        P = polyfit(ang,tq,order);
        y = polyval(P,x);
        plot(x,y,'k:');

        %% Fit only start values including passive rotations
    elseif se ~= 1 && usePas == 1
        ang = [angFitB(:,1); angFitPas{col}];
        tq = [tqFitB(:,1); tqFitPas{col}];
        P = polyfit(ang,tq,order);
        y = polyval(P,x);
        plot(x,y,'k:');

        %% Fit start and end values not including passive rotations
    elseif se == 1 && usePas ~= 1
        ang = vertcat(angFitB(:,1),angFitE(:,1));
        tq = vertcat(tqFitB(:,1),tqFitE(:,1));
        P = polyfit(ang,tq,order);
        y = polyval(P,x);
        plot(x,y,'k:');

    elseif se == 1 && usePas == 1

        %% Fit start and end values including passive rotations
        ang = vertcat(angFitB(:,1),angFitE(:,1),angFitPas{col});
        tq = vertcat(tqFitB(:,1),tqFitE(:,1),tqFitPas{col});
        P = polyfit(ang,tq,order);
        y = polyval(P,x);
        plot(x,y,'k:');

    end
end

%% Show mean passive torque values and remove outliers
figure(length(files)+2); hold on
%ang = gnegate(ang);
hPts = plot(ang,tq,'o');
xlabel('Plantar flexion angle (deg)')
ylabel('Net ankle joint torque (Nm)')
brush on
title('Select points to keep by highlighting them and clicking confirm')
hPos = get(gcf);
h = uicontrol('String', 'Confirm', 'Position', [hPos.OuterPosition(1,3)-120 10 100 50], ...
    'Callback', 'uiresume(gcbf)');
uiwait(gcf);
ptsSelected = logical(hPts.BrushData.');
data = [hPts.XData(ptsSelected).' ...
    hPts.YData(ptsSelected).'];
names = {'use'};
assignin('base',names{1},data)
close(gcf);

figure(length(files)+1); hold on

% Plot fits depending on inputs
if se ~= 1
    plot(angFitB,tqFitB,'ro');
    plot(data(:,1),data(:,2),'ko');
    P = polyfit(data(:,1),data(:,2),order);
    y = polyval(P,x);
    plot(x,y,'k.');
    plot(angFitE,tqFitE,'rx')

else
    plot(ang,tq,'ro');
    plot(data(:,1),data(:,2),'ko');
    P = polyfit(data(:,1),data(:,2),order);
    y = polyval(P,x);
    plot(x,y,'k.');
end

ang = data(:,1); tq = data(:,2); clear data

%% Save analysed data
save('TQangFit.mat')
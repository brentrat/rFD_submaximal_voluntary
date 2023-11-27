%% Summary of participants from rFD experiment (1st DFG funding)
% Run passive_torque.m before this script!
% Run active_torque.m before this script!
% Initialise workspace
clear; close all

%% Inputs
pathName = 'C:\Users\Brent\Desktop\Research\DFG\MVC_TA';
RMtext = 0; % Save textfile for linear repeated-measures correlation

%% Load necessary file to normalize EMG data
load('C:\Users\Brent\Desktop\Research\DFG\MVC_TA\mvcNORM_10Hz.mat')

%% Loop through files
cd(pathName)

for jj =  1:16

    cd(pathName)

    if jj < 10
        folder = ['0' num2str(jj)];
    else
        folder = num2str(jj);
    end

    cd(folder)
    cd('Tracked')

    %% Load variables of interest
    active_all{jj,1} = load('analyzedJP.mat','A','ANG','fRef','lenS','lenL','lenLP','shoS','shoSP','shoL','shoLP','AfasDiff','hold3r','HRef','WRef');
    load('analyzedJP.mat','ForceRef','ForceDyn','emgRef','emgDyn','angDyn','angRef','tqDyn','tqRef','FLdyn','FLref','fRef','lenS','lenL','lenLP','shoS','shoSP','shoL','shoLP','v','rot_start','fDyn');

    %% Determine number of submaximal contractions per participant
    % Determine number of recorded files
    cd ..
    files = dir;
    files = {files.name};
    cond = {'ref','lens','lenl','lenlp','shos','shosp','shol','sholp'};

    % Determine number of analyzed trials
    condv = [length(fRef) length(lenS) length(lenL) length(lenLP) length(shoS) length(shoSP) length(shoL) length(shoLP)];
    for bb = 1:8
        C{bb}(jj,1) = sum(contains(files,cond{bb}));
        C{bb}(jj,2) = condv(bb);
    end

    %% Determine start of angular rotation
    ANGstart{jj,1} = rot_start;
    ANGstart{jj,2} = fDyn;

    %% Make a nice figure
    E = mvcEMG(jj,1); % variable to normalize individual EMG data

    % C1 = REF, C2 = LENsmall etc.
    % 1st EMG, 2nd angle, 3rd torque, 4th fascicle length, 5th fascicle
    % force
    if size(emgRef,1) == length(fRef)

        C1a(jj,:) = mean(emgRef/E*100,1,'omitnan');
        C1b(jj,:) = mean(angRef,1,'omitnan');
        C1c(jj,:) = mean(tqRef,1,'omitnan');
        C1d(jj,:) = mean(FLref,1,'omitnan');
        C1e(jj,:) = mean(ForceRef,1,'omitnan');

    end

    %% I should be shot for this code
    C2a(jj,:) = mean(emgDyn(lenS,:)/E*100,1,'omitnan');
    C2b(jj,:) = mean(angDyn(lenS,:),1,'omitnan');
    C2c(jj,:) = mean(tqDyn(lenS,:),1,'omitnan');
    C2d(jj,:) = mean(FLdyn(lenS,:),1,'omitnan');
    C2e(jj,:) = mean(ForceDyn(lenS,:),1,'omitnan');

    C3a(jj,:) = mean(emgDyn(lenL,:)/E*100,1,'omitnan');
    C3b(jj,:) = mean(angDyn(lenL,:),1,'omitnan');
    C3c(jj,:) = mean(tqDyn(lenL,:),1,'omitnan');
    C3d(jj,:) = mean(FLdyn(lenL,:),1,'omitnan');
    C3e(jj,:) = mean(ForceDyn(lenL,:),1,'omitnan');

    C4a(jj,:) = mean(emgDyn(lenLP,:)/E*100,1,'omitnan');
    C4b(jj,:) = mean(angDyn(lenLP,:),1,'omitnan');
    C4c(jj,:) = mean(tqDyn(lenLP,:),1,'omitnan');
    C4d(jj,:) = mean(FLdyn(lenLP,:),1,'omitnan');
    C4e(jj,:) = mean(ForceDyn(lenLP,:),1,'omitnan');

    C5a(jj,:) = mean(emgDyn(v+shoS,:)/E*100,1,'omitnan');
    C5b(jj,:) = mean(angDyn(v+shoS,:),1,'omitnan');
    C5c(jj,:) = mean(tqDyn(v+shoS,:),1,'omitnan');
    C5d(jj,:) = mean(FLdyn(v+shoS,:),1,'omitnan');
    C5e(jj,:) = mean(ForceDyn(v+shoS,:),1,'omitnan');

    C6a(jj,:) = mean(emgDyn(v+shoSP,:)/E*100,1,'omitnan');
    C6b(jj,:) = mean(angDyn(v+shoSP,:),1,'omitnan');
    C6c(jj,:) = mean(tqDyn(v+shoSP,:),1,'omitnan');
    C6d(jj,:) = mean(FLdyn(v+shoSP,:),1,'omitnan');
    C6e(jj,:) = mean(ForceDyn(v+shoSP,:),1,'omitnan');

    C7a(jj,:) = mean(emgDyn(v+shoL,:)/E*100,1,'omitnan');
    C7b(jj,:) = mean(angDyn(v+shoL,:),1,'omitnan');
    C7c(jj,:) = mean(tqDyn(v+shoL,:),1,'omitnan');
    C7d(jj,:) = mean(FLdyn(v+shoL,:),1,'omitnan');
    C7e(jj,:) = mean(ForceDyn(v+shoL,:),1,'omitnan');

    C8a(jj,:) = mean(emgDyn(v+shoLP,:)/E*100,1,'omitnan');
    C8b(jj,:) = mean(angDyn(v+shoLP,:),1,'omitnan');
    C8c(jj,:) = mean(tqDyn(v+shoLP,:),1,'omitnan');
    C8d(jj,:) = mean(FLdyn(v+shoLP,:),1,'omitnan');
    C8e(jj,:) = mean(ForceDyn(v+shoLP,:),1,'omitnan');

    %% Determine fatigue across fixed-end reference trials
    [~,ndx] = natsort(active_all{jj,1}.fRef);

    % Active torque
    Fat(jj,1) = active_all{jj,1}.hold3r(ndx(1),1);
    Fat(jj,2) = active_all{jj,1}.hold3r(ndx(end),1);

    % Absolute EMG
    Fat(jj,3) = active_all{jj,1}.hold3r(ndx(1),2);
    Fat(jj,4) = active_all{jj,1}.hold3r(ndx(end),2);

    % Normalized active torque
    Fat(jj,5) = active_all{jj,1}.hold3r(ndx(1),4);
    Fat(jj,6) = active_all{jj,1}.hold3r(ndx(end),4);

    % Normalized EMG
    Fat(jj,7) = active_all{jj,1}.hold3r(ndx(1),5);
    Fat(jj,8) = active_all{jj,1}.hold3r(ndx(end),5);

    % Net fascicle work
    Fat(jj,9) = active_all{jj,1}.WRef(ndx(1),1);
    Fat(jj,10) = active_all{jj,1}.WRef(ndx(end),1);

    % Positive fascicle work
    Fat(jj,11) = active_all{jj,1}.WRef(ndx(1),2);
    Fat(jj,12) = active_all{jj,1}.WRef(ndx(end),2);

    % Fascicle shortening
    Fat(jj,13) = active_all{jj,1}.HRef(ndx(1),1);
    Fat(jj,14) = active_all{jj,1}.HRef(ndx(end),1);

    % Steady-state fascicle length (3rd hold phase)
    Fat(jj,15) = active_all{jj,1}.HRef(ndx(1),4);
    Fat(jj,16) = active_all{jj,1}.HRef(ndx(end),4);

    %% Store data for later statistical analysis
    % Active torque
    TQ_all(1:8,jj) = round(mean(active_all{jj,1}.A.TQ(1:8,3),2),1);
    TQ_mvc(1:8,jj) = round(mean(active_all{jj,1}.A.TQmvc(1:8,3),2),1);

    % Active fascicle force
    FORCE_all(1:8,jj) = round(mean(active_all{jj,1}.A.FORCE(1:8,3),2),1);

    % Normalized EMG
    EMG_all(1:8,jj) = mean(active_all{jj,1}.A.EMGmvc(1:8,3),2);

    % Rotation amplitudes
    ROT_all(1:8,jj) = active_all{jj,1}.A.ROT(1:8,1);

    % Absolute EMG area until the steady state
    EMG_area(1:8,jj) = active_all{jj,1}.A.EMGrot(1:8,1);

    % Normalized mean EMG until the steady state
    EMG_meanSS(1:8,jj) = active_all{jj,1}.A.EMGrot(1:8,2);

    % Normalized mean EMG until the end of the rotation
    EMG_mean(1:8,jj) = active_all{jj,1}.A.EMGrot(1:8,3);

    % Normalized min. and max. EMG until the end of the rotation
    EMG_min(1:8,jj) = active_all{jj,1}.A.EMGrot(1:8,4);
    EMG_max(1:8,jj) = active_all{jj,1}.A.EMGrot(1:8,5);

    % Normalized mean EMG during rotation
    EMG_meanRot(1:8,jj) = active_all{jj,1}.A.EMGrot(1:8,6);

    % Mean active fascicle force until the steady state
    Force_meanSS(1:8,jj) = active_all{jj,1}.A.FORrot(1:8,1);

    % Mean fascicle velocity
    FV_mean(1:8,jj) = active_all{jj,1}.A.FORrot(1:8,2);

    % Min. and max. active fascicle force until the end of the rotation
    Force_end(1:8,jj) = active_all{jj,1}.A.FORrot(1:8,3);
    Force_start(1:8,jj) = active_all{jj,1}.A.FORrot(1:8,4);

    % Mean active fascicle force until the end of the rotation
    Force_mean(1:8,jj) = active_all{jj,1}.A.FORrot(1:8,5);

    % Mean active MTU force and MTU shortening amplitude
    FAS_forceRot(1:8,jj) = active_all{jj,1}.A.FORrot(1:8,6);
    MTU_force(1:8,jj) = active_all{jj,1}.A.FORrot(1:8,7);
    MTU_forceRot(1:8,jj) = active_all{jj,1}.A.FORrot(1:8,8);
    MTU_amp(1:8,jj) = active_all{jj,1}.A.FORrot(1:8,9);

    if sum(strcmp('FAS',active_all{jj, 1}.A.Properties.VariableNames)) == 1

        % Fascicle length during the steady-state (3rd hold phase)
        FAS_hold(1:8,jj) = active_all{jj,1}.A.FAS(1:8,3);

        % Fascicle shortening amplitude
        FAS_sho(1:8,jj) = active_all{jj,1}.A.FASsho(1:8,1);

        % Fascicle stretch amplitude
        FAS_str(1:8,jj) = max([active_all{jj,1}.A.FASstr(1:8,1) active_all{jj,1}.A.FASstr(1:8,3)],[],2);

        % Fascicle work
        FAS_work(1:8,jj) = active_all{jj,1}.A.FASwork(1:8,1);
        FAS_workPos(1:8,jj) = active_all{jj,1}.A.FASwork(1:8,2);
        FAS_workNeg(1:8,jj) = active_all{jj,1}.A.FASwork(1:8,3);

        % MTU work
        MTU_work(1:8,jj) = active_all{jj,1}.A.MTUwork(1:8,1);
        MTU_workPos(1:8,jj) = active_all{jj,1}.A.MTUwork(1:8,2);
        MTU_workNeg(1:8,jj) = active_all{jj,1}.A.MTUwork(1:8,3);

        % Joint work
        ANG_work(1:8,jj) = active_all{jj,1}.A.ANGwork(1:8,1);
        ANG_workPos(1:8,jj) = active_all{jj,1}.A.ANGwork(1:8,2);
        ANG_workNeg(1:8,jj) = active_all{jj,1}.A.ANGwork(1:8,3);

    end

    % Common ankle angle
    ANG_tested(jj,1) = round(mean(active_all{jj,1}.ANG(:,3)),1);

    % Rotation amplitude
    ANG_rot(jj,:) = round(active_all{jj,1}.A.ROT,1);

end

% Number of recorded and valid trials per condition
% rows: mean, min, max
% columns: ref([all valid]), lens+lenm([all valid]), 
% lenlp+shos+shosp+shol+sholp([all valid])
trials(1,1:2) = round(mean(C{1}),0);
trials(2,1:2) = round(min(C{1}),0);
trials(3,1:2) = round(max(C{1}),0);
lensm = vertcat(C{2},C{3});
trials(1,3:4) = round(mean(lensm),0);
trials(2,3:4) = round(min(lensm),0);
trials(3,3:4) = round(max(lensm),0);
lenlpshospsholp = vertcat(C{4},C{5},C{6},C{7},C{8});
trials(1,5:6) = round(mean(lenlpshospsholp),0);
trials(2,5:6) = round(min(lenlpshospsholp),0);
trials(3,5:6) = round(max(lenlpshospsholp),0);
reflensmlpshospsholp = horzcat(C{1}(:,1),C{2}(:,1),C{3}(:,1),C{4}(:,1),C{5}(:,1),C{6}(:,1),C{7}(:,1),C{8}(:,1));
trialsAll(1,1) = round(mean(sum(reflensmlpshospsholp,2)));
trialsAll(2,1) = round(min(sum(reflensmlpshospsholp,2)));
trialsAll(3,1) = round(max(sum(reflensmlpshospsholp,2)));
reflensmlpshospsholpV = horzcat(C{1}(:,2),C{2}(:,2),C{3}(:,2),C{4}(:,2),C{5}(:,2),C{6}(:,2),C{7}(:,2),C{8}(:,2));
trialsAll(1,2) = round(mean(sum(reflensmlpshospsholpV,2)));
trialsAll(2,2) = round(min(sum(reflensmlpshospsholpV,2)));
trialsAll(3,2) = round(max(sum(reflensmlpshospsholpV,2)));

% Remove zeros in matrix
EMG_all(EMG_all==0) = NaN;
TQ_all(TQ_all==0) = NaN;
TQ_mvc(TQ_mvc==0) = NaN;
ROT_all(ROT_all==0) = NaN;
ROT_all = ROT_all'; % transpose
ROT_all(:,1) = 0;

if exist('FAS_hold','var')

    % Remove zeros in matrix
    FORCE_all(FORCE_all==0) = NaN;
    FAS_hold(FAS_hold==0) = NaN;
    FAS_sho(FAS_sho==0) = NaN;
    FAS_sho = FAS_sho'*-1; % transpose and make negative values positive
    FAS_str(2,7) = NaN;
    FAS_str(3,11) = NaN;
    FAS_str = FAS_str'; % transpose
    FAS_work = FAS_work'; % transpose
    FAS_work(7,2) = NaN;
    FAS_work(11,3) = NaN;
    FAS_workPos(2,7) = NaN;
    FAS_workPos(3,11) = NaN;
    FAS_workPos = FAS_workPos'; % transpose
    FAS_workNeg(2,7) = NaN;
    FAS_workNeg(3,11) = NaN;
    FAS_workNeg = FAS_workNeg'; % transpose
    MTU_work(2,7) = NaN;
    MTU_work(3,11) = NaN;
    MTU_work = MTU_work'; % transpose
    MTU_workPos(2,7) = NaN;
    MTU_workPos(3,11) = NaN;
    MTU_workPos = MTU_workPos'; % transpose
    MTU_workNeg(2,7) = NaN;
    MTU_workNeg(3,11) = NaN;
    MTU_workNeg = MTU_workNeg'; % transpose
    FAS_forceRot(2,7) = NaN;
    FAS_forceRot(3,11) = NaN;
    FAS_forceRot = FAS_forceRot'; % transpose
    MTU_force(2,7) = NaN;
    MTU_force(3,11) = NaN;
    MTU_force = MTU_force'; % transpose
    MTU_forceRot(2,7) = NaN;
    MTU_forceRot(3,11) = NaN;
    MTU_forceRot = MTU_forceRot'; % transpose
    MTU_amp(2,7) = NaN;
    MTU_amp(3,11) = NaN;
    MTU_amp = MTU_amp';% transpose

end

% Transpose a bunch of data
TQ_all = TQ_all';
TQ_mvc = TQ_mvc';
FORCE_all = FORCE_all';
FAS_hold = FAS_hold';
EMG = EMG_all';
%EMG = round(EMG./mvcEMG*100,2);
% No need to remove zeros here as other variable in correlation has NaNs in
% their place
EMG_area = EMG_area';
EMG_meanSS = EMG_meanSS';
EMG_mean = EMG_mean';
EMG_min = EMG_min';
EMG_max = EMG_max';
EMG_meanRot = EMG_meanRot';
Force_meanSS = Force_meanSS';
Force_mean = Force_mean';
FV_mean(2,7) = NaN;
FV_mean(3,11) = NaN;
FV_mean = FV_mean';
Force_end = Force_end';
Force_start = Force_start';
ANG_rot(ANG_rot==0) = NaN;
ANG_rot(:,1) = 0;
ANG_work(2,7) = NaN;
ANG_work(3,11) = NaN;
ANG_work = round(ANG_work,1)'; % transpose
ANG_workPos(2,7) = NaN;
ANG_workPos(3,11) = NaN;
ANG_workPos = round(ANG_workPos,1)'; % transpose
ANG_workNeg(2,7) = NaN;
ANG_workNeg(3,11) = NaN;
ANG_workNeg = round(ANG_workNeg,1)'; % transpose


%% Summary rFD data
% Fatigue
Fat = round(Fat(:,1:end),2);
% rFD active torque (%)
TQ_diff = round((TQ_all(:,1:end)-TQ_all(:,1))./TQ_all(:,1)*100,2);

% rFD active fascicle force (%)
FORCE_diff = round((FORCE_all(:,1:end)-FORCE_all(:,1))./FORCE_all(:,1)*100,2);

% rFE active fascicle force (N)
FORCE_diffAbs = round((FORCE_all(:,2:end)-FORCE_all(:,1)),2);

%% Go back to main directory
cd(pathName)

%% Reliability
counter = 0;

for aa = 1:2:size(Fat,2)

    counter = counter+1;

    var1 = Fat(:,aa);
    var2 = Fat(:,aa+1);

    p(:,1) = 1:1:length(var1);
    t = table(p,var1,var2,'VariableNames',{'Participant';'Var1';'Var2'});
    meas = table([1 2]','VariableNames',{'Measurements'});
    rm = fitrm(t,'Var1-Var2 ~ Participant','WithinDesign',meas);
    ranovatbl = ranova(rm);
    SEM(counter,1) = sqrt(ranovatbl{3,3});
    tCrit = tinv(1-0.05/2, length(var1)-1);
    MD(counter,1) = std(var1-var2)*tCrit;
    MD2(counter,1) = SEM(counter,1)*sqrt(2)*tCrit;
    ICC_all = icc([var1 var2],0.05);
    ICC(counter,:) = [ICC_all{1,3}.est ICC_all{1,3}.confInterval];

end

Rel = table(round(ICC,2),round(SEM,1),round(MD,1),round(MD2,1));
Rel.Properties.RowNames = {'Torque (Nm)','EMG (V)','Norm. torque (%)','Norm. EMG (%)','Net work (J)','Pos. work (J)', 'Fas. sho. (mm)', 'Fas. length (mm)'};
clear ICC SEM MD* p

%% Prism statistical analysis
Prism.fatigueTorqueEmgPercent = [Fat(:,5) Fat(:,6) Fat(:,7) Fat(:,8)];
Prism.emg = EMG;
Prism.torque = TQ_all;
Prism.torquePercent = round((TQ_all-TQ_all(:,1))./TQ_all(:,1)*100,2);
Prism.force = FORCE_all;
Prism.forcePercent = round((FORCE_all-FORCE_all(:,1))./FORCE_all(:,1)*100,2);

% Grouped dataset
torquePercentSummary = [mean(Prism.torquePercent,'omitnan'); std(Prism.torquePercent,'omitnan'); sum(~isnan(Prism.torquePercent))];
forcePercentSummary = [mean(Prism.forcePercent,'omitnan'); std(Prism.forcePercent,'omitnan'); sum(~isnan(Prism.forcePercent))];

for jj = 2:8
    Prism.torqueForcePercentSummary(jj-1,1:6) = round([torquePercentSummary(:,jj)' forcePercentSummary(:,jj)'],2);
end

Prism.fascicleLength = FAS_hold;
Prism.fascicleStretch = FAS_str*-1;
Prism.fascicleShortening = FAS_sho;
Prism.fascicleNetWork = round(FAS_work,1);
Prism.mtuNetWork = round(MTU_work,1);
Prism.angularWork = round(ANG_work,1);

fascicleShorteningSummary = [mean(Prism.fascicleShortening,'omitnan'); std(Prism.fascicleShortening,'omitnan'); sum(~isnan(Prism.fascicleShortening))];
fascicleStretchSummary = [mean(Prism.fascicleStretch,'omitnan'); std(Prism.fascicleStretch,'omitnan'); sum(~isnan(Prism.fascicleStretch))];

for jj = 1:8
    Prism.fascicleLengthChangeSummary(jj,1:6) = round([fascicleShorteningSummary(:,jj)' fascicleStretchSummary(:,jj)'],1);
end

FASnet = [mean(FAS_work,'omitnan'); std(FAS_work,'omitnan'); sum(~isnan(FAS_work))];
FASpos = [mean(FAS_workPos,'omitnan'); std(FAS_workPos,'omitnan'); sum(~isnan(FAS_workPos))];
FASneg = [mean(FAS_workNeg,'omitnan'); std(FAS_workNeg,'omitnan'); sum(~isnan(FAS_workNeg))];

for jj = 1:8
    Prism.fascicleWorkAll(jj,1:9) = round([FASnet(:,jj)' FASpos(:,jj)' FASneg(:,jj)'],1);
end

MTUnet = [mean(MTU_work,'omitnan'); std(MTU_work,'omitnan'); sum(~isnan(MTU_work))];
MTUpos = [mean(MTU_workPos,'omitnan'); std(MTU_workPos,'omitnan'); sum(~isnan(MTU_workPos))];
MTUneg = [mean(MTU_workNeg,'omitnan'); std(MTU_workNeg,'omitnan'); sum(~isnan(MTU_workNeg))];

for jj = 1:8
    Prism.mtuWorkAll(jj,1:9) = round([MTUnet(:,jj)' MTUpos(:,jj)' MTUneg(:,jj)'],1);
end

ANGnet = [mean(ANG_work,'omitnan'); std(ANG_work,'omitnan'); sum(~isnan(ANG_work))];
ANGpos = [mean(ANG_workPos,'omitnan'); std(ANG_workPos,'omitnan'); sum(~isnan(ANG_workPos))];
ANGneg = [mean(ANG_workNeg,'omitnan'); std(ANG_workNeg,'omitnan'); sum(~isnan(ANG_workNeg))];

for jj = 1:8
    Prism.angularWorkAll(jj,1:9) = round([ANGnet(:,jj)' ANGpos(:,jj)' ANGneg(:,jj)'],1);
end

%% Prism repeated-measures scatterplot figure
% A = Net fascicle work
FASWORK = FAS_work(:,5:8);
RFD = Prism.forcePercent(:,5:8)*-1;

c = 1;
for jj = 1:16
    if jj == 1
        FASfig(c:c+3,1) = FASWORK(jj,:)';
        RFDfig(1:4,jj) = RFD(jj,:)';
    else
        FASfig(1:c+3,1) = vertcat(FASfig,FASWORK(jj,:)');
        RFDfig(1:4,jj) = RFD(jj,:)';
    end
    c = c+4;
end

FASfig = vertcat(FASfig, mean(FASWORK,2));
FASfig = round(FASfig,1);
RFDfig = round(RFDfig,2);
RFDfigMean = round(mean(RFD,2),2);

% B - Net MTU work
MTUWORK = MTU_work(:,5:8);

c = 1;
for jj = 1:16
    if jj == 1
        MTUfig(c:c+3,1) = MTUWORK(jj,:)';
    else
        MTUfig(1:c+3,1) = vertcat(MTUfig,MTUWORK(jj,:)');
    end
    c = c+4;
end

MTUfig = vertcat(MTUfig, mean(MTUWORK,2));
MTUfig = round(MTUfig,1);

%% Prism time-trace figure
REF(:,2) = mean(C1a,1,'omitnan'); REF(:,4) = mean(C1b,1,'omitnan'); REF(:,1) = mean(C1c,1,'omitnan'); REF(:,3) = mean(C1d,1,'omitnan');
LENS(:,2) = mean(C2a,1,'omitnan'); LENS(:,4) = mean(C2b,1,'omitnan'); LENS(:,1) = mean(C2c,1,'omitnan'); LENS(:,3) = mean(C2d,1,'omitnan');
LENL(:,2) = mean(C3a,1,'omitnan'); LENL(:,4) = mean(C3b,1,'omitnan'); LENL(:,1) = mean(C3c,1,'omitnan'); LENL(:,3) = mean(C3d,1,'omitnan');
LENLP(:,2) = mean(C4a,1,'omitnan'); LENLP(:,4) = mean(C4b,1,'omitnan'); LENLP(:,1) = mean(C4c,1,'omitnan'); LENLP(:,3) = mean(C4d,1,'omitnan');
SHOS(:,2) = mean(C5a,1,'omitnan'); SHOS(:,4) = mean(C5b,1,'omitnan'); SHOS(:,1) = mean(C5c,1,'omitnan'); SHOS(:,3) = mean(C5d,1,'omitnan');
SHOSP(:,2) = mean(C6a,1,'omitnan'); SHOSP(:,4) = mean(C6b,1,'omitnan'); SHOSP(:,1) = mean(C6c,1,'omitnan'); SHOSP(:,3) = mean(C6d,1,'omitnan');
SHOL(:,2) = mean(C7a,1,'omitnan'); SHOL(:,4) = mean(C7b,1,'omitnan'); SHOL(:,1) = mean(C7c,1,'omitnan'); SHOL(:,3) = mean(C7d,1,'omitnan');
SHOLP(:,2) = mean(C8a,1,'omitnan'); SHOLP(:,4) = mean(C8b,1,'omitnan'); SHOLP(:,1) = mean(C8c,1,'omitnan'); SHOLP(:,3) = mean(C8d,1,'omitnan');
FIG = table(REF,LENS,LENL,LENLP,SHOS,SHOSP,SHOL,SHOLP);

FORCE(:,1) = mean(C1e,1,'omitnan'); FORCE(:,2) = mean(C2e,1,'omitnan'); FORCE(:,3) = mean(C3e,1,'omitnan'); FORCE(:,4) = mean(C4e,1,'omitnan');
FORCE(:,5) = mean(C5e,1,'omitnan'); FORCE(:,6) = mean(C6e,1,'omitnan'); FORCE(:,7) = mean(C7e,1,'omitnan'); FORCE(:,8) = mean(C8e,1,'omitnan');

%% Linear repeated-measures correlation
%https://lmarusich.shinyapps.io/shiny_rmcorr/
if RMtext == 1

    % Create textfile with x and y variables and the number of repeated
    % measurements
    x = [];
    y = [];
    n = [];

    % Define x and y variables to be correlated
    for jj = 1:16
%         x = vertcat(x,(FORCE_diff(jj,[2:4]))'*-1);
%         y = vertcat(y,(FAS_stretch(jj,[2:4]))'*-1); 
%         x = vertcat(x,(FORCE_diff(jj,[5:8]))'*-1);
%         y = vertcat(y,(FAS_work(jj,[5:8]))'); 
         x = vertcat(x,(FORCE_diff(jj,[5:8]))'*-1);
         y = vertcat(y,(MTU_work(jj,[5:8]))'); 
    end

    % Define number of repeated measurements
    for jj = 1:16
%         n = vertcat(n,ones(3,1)*jj);
        n = vertcat(n,ones(4,1)*jj);
    end

    % Create table
    clear T
    T = table(n,x,y, 'VariableNames', { 'Subject', 'rFD', 'Work' });

    % Store textfile in new location
    cd('C:\Users\Brent\Desktop\Research\DFG\MVC_TA\RM\New')

    % Save table as textfile
    writetable(T, ['RM_rFD_vs_MTUWork_sho.txt']); 

end

save('analyzedJP.mat')

fn = fieldnames(Prism);
fn([7 14:17]) = [];
T = table(Prism.(fn{1}),Prism.(fn{2}),Prism.(fn{3}),Prism.(fn{4}),Prism.(fn{5}),Prism.(fn{6}),Prism.(fn{7}),Prism.(fn{8}),Prism.(fn{9}),Prism.(fn{10}),Prism.(fn{11}),Prism.(fn{12}));
T.Properties.VariableNames = fn;
writetable(T,'analyzedJP.xlsx');
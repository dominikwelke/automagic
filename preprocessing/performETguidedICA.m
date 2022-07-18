function [EEG, ET] = performETguidedICA(EEG, params)

% For details see the underlying publication: Dimigen, 2020, NeuroImage


% Copyright (C) 2017  Amirreza Bahreini, methlabuzh@gmail.com
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Check and unzip if EYE_EEG() does not exist
addEYE_EEG()

%% params
et_namePattern = params.fileext_edit;
eeg_fileName = EEG.fileName;
x = split(et_namePattern, '.');
et_fileExt = x{end};
et_datafolder = params.datafolder_edit; 
% dataFolder = EEG.dataFolder;
dataFolder = et_datafolder;

L_GAZE_X = params.l_gaze_x_edit;
L_GAZE_Y = params.l_gaze_y_edit;
R_GAZE_X = params.r_gaze_x_edit;
R_GAZE_Y = params.r_gaze_y_edit;


% build a structure (can migrate to defaults later?)
ETparams = struct;
% ET setup
ETparams.setup.SCREEN_RES_X       = str2double(params.screenWidth_edit);  % in px
ETparams.setup.SCREEN_RES_Y       = str2double(params.screenHeight_edit); % in px
ETparams.setup.VIEW_DIST          = 720;   % in mm
ETparams.setup.SCREEN_WIDTH       = 530;   % in mm
ETparams.setup.DEG_PER_PIXEL      = px2deg(ETparams.setup.SCREEN_WIDTH,ETparams.setup.SCREEN_RES_X,ETparams.setup.VIEW_DIST); % calculate transfer from screen px to deg vis. angle
% syncronization
ETparams.sync.importEyeEvents     = false; % dont import them, recompute saccades
ETparams.sync.PLOTFIG             = 0;     % plot and save the sync figure
% ET based data rejection
ETparams.markbad.REJECTMODE       = 2;     % don't reject data, add extra "bad_ET" events to EEG.event
% saccade detection
ETparams.detection.THRESH         = 6;     % eye velocity threshold (in median-based SDs)
ETparams.detection.MINDUR         = 4;     % minimum saccade duration (samples)
ETparams.detection.SMOOTH         = 1;     % smooth eye velocities? (recommended if SR > 250 Hz)
ETparams.detection.PLOTFIG        = 1;     % plot and save the saccade detection figure
ETparams.detection.WRITESAC       = 1;     % add saccades as events to EEG.event?
ETparams.detection.WRITEFIX       = 1;     % add fixations as events to EEG.event?
% OPTICAT proper
ETparams.OPTICAT.OW_PROPORTION    = 1.0;   % overweighting proportion
ETparams.OPTICAT.SACCADE_WINDOW   = [str2double(params.from_edit) str2double(params.to_edit)];  % time window to overweight (-20 to 10 ms is default)
ETparams.OPTICAT.REMOVE_EPOCHMEAN = true;  % subtract mean from overweighted epochs? (recommended)


%% find first and last trigger, if not provided

if isempty(params.startTrigger_edit) & isempty(params.endTrigger_edit)
    ETparams.sync.startTrigger = str2double(EEG.event(1).type);
    ETparams.sync.endTrigger = str2double(EEG.event(end).type);
else
    ETparams.sync.startTrigger = str2double(params.startTrigger_edit);
    ETparams.sync.endTrigger = str2double(params.endTrigger_edit);
end
%ETparams.sync.endTrigger          = 56;    % overwriting it here

%% find corresponding et file - tricky, if more ET files in a folder
d = dir(fullfile(dataFolder, ['*' , et_namePattern]));

if length(d) == 1
    et_fileName = d(1).name;   
else % assume that that the filenames for EEG and ET are identical up to "_eeg/_EEG" vs et_namePattern
    i = regexp(lower(eeg_fileName), '_eeg');
    patt = eeg_fileName(1:i); % remove "eeg" from the name
    if strcmp(et_namePattern(1), '_')
        et_fileName = [patt(1:end-1), et_namePattern];
    else
        et_fileName = [patt, et_namePattern];
    end
end


%% if .txt, convert to .mat and save as a mat file
if strcmp(et_fileExt, 'txt')
    ET = parsesmi(fullfile(dataFolder, et_fileName), dataFolder);
elseif strcmp(et_fileExt, 'mat')
    ET = load(fullfile(dataFolder, et_fileName));
end

%% import & synchronize ET data
ETparams.sync.importColumns = 2:length(ET.colheader); % exclude time trace
ETparams.sync.newLabels = ET.colheader(2:end);

EEG = pop_importeyetracker(EEG, fullfile(dataFolder, et_fileName), ...
    [ETparams.sync.startTrigger, ETparams.sync.endTrigger], ETparams.sync.importColumns, ETparams.sync.newLabels, ETparams.sync.importEyeEvents,1,0,ETparams.sync.PLOTFIG,4);

% save the figure
if ETparams.sync.PLOTFIG
    ETparams.sync.fname = [et_datafolder et_fileName(1:end-4) '_sync.fig'];
    savefig(gcf,ETparams.sync.fname)
    close gcf
end

%% Mark intervals with bad eye tracking data
% important, so these intervals will not influence our saccade detection
% This function is also useful to objectively reject intervals during
% which the participant blinked or did not look at the stimulus

% eye position channels
ETrecmode = [];
if any(strcmp({EEG.chanlocs(:).labels},L_GAZE_X)) & ~any(strcmp({EEG.chanlocs(:).labels}, R_GAZE_X))
    LX = find(strcmp({EEG.chanlocs(:).labels}, L_GAZE_X));
    LY = find(strcmp({EEG.chanlocs(:).labels}, L_GAZE_Y));
    ETparams.markbad.rej_chans = [LX LY];
    ETparams.markbad.rej_minvals = [1 1];
    ETparams.markbad.rej_maxvals = [ETparams.setup.SCREEN_RES_X ETparams.setup.SCREEN_RES_Y];
    ETrecmode = 'L';
    left_eye_xy = [LX LY];
    right_eye_xy = [];
elseif any(strcmp({EEG.chanlocs(:).labels}, R_GAZE_X)) & ~any(strcmp({EEG.chanlocs(:).labels}, L_GAZE_X))
    RX = find(strcmp({EEG.chanlocs(:).labels}, R_GAZE_X));
    RY = find(strcmp({EEG.chanlocs(:).labels}, R_GAZE_Y));
    ETparams.markbad.rej_chans = [RX RY];
    ETparams.markbad.rej_minvals = [1 1];
    ETparams.markbad.rej_maxvals = [ETparams.setup.SCREEN_RES_X ETparams.setup.SCREEN_RES_Y];
    ETrecmode = 'R';
    left_eye_xy = [];
    right_eye_xy = [RX RY];
elseif any(strcmp({EEG.chanlocs(:).labels}, L_GAZE_X)) & any(strcmp({EEG.chanlocs(:).labels}, R_GAZE_X))
    LX = find(strcmp({EEG.chanlocs(:).labels}, L_GAZE_X));
    LY = find(strcmp({EEG.chanlocs(:).labels}, L_GAZE_Y));
    RX = find(strcmp({EEG.chanlocs(:).labels}, R_GAZE_X));
    RY = find(strcmp({EEG.chanlocs(:).labels}, R_GAZE_Y));
    ETparams.markbad.rej_chans = [LX LY RX RY];
    ETparams.markbad.rej_minvals = [1 1 1 1];
    ETparams.markbad.rej_maxvals = [ETparams.setup.SCREEN_RES_X ETparams.setup.SCREEN_RES_Y ETparams.setup.SCREEN_RES_X ETparams.setup.SCREEN_RES_Y];
    ETrecmode = 'B';
    left_eye_xy = [LX LY];
    right_eye_xy = [RX RY];
end

EEG = pop_rej_eyecontin(EEG, ETparams.markbad.rej_chans, ETparams.markbad.rej_minvals, ETparams.markbad.rej_maxvals, 25, ETparams.markbad.REJECTMODE);

%% Detect (micro)saccades & fixations (Engbert & Kliegl, 2003)

% ### GUI: "Eyetracker" > "Detect saccades & fixations"
% % see "help pop_detecteyemovements" to see all options
% % 

EEG = pop_detecteyemovements(EEG,left_eye_xy,right_eye_xy,ETparams.detection.THRESH,ETparams.detection.MINDUR,ETparams.setup.DEG_PER_PIXEL,ETparams.detection.SMOOTH,0,25,2,ETparams.detection.PLOTFIG,ETparams.detection.WRITESAC,ETparams.detection.WRITEFIX);

% save the figure
if ETparams.detection.PLOTFIG
    ETparams.detection.fname = [et_datafolder et_fileName(1:end-4) '_movements.fig'];
    savefig(gcf,ETparams.detection.fname)
    close gcf
end

%% keep ET data for reattaching it at the end of the preprocessing stream
nbchan = length(ETparams.sync.importColumns);
ET = EEG;
ET.data = ET.data(end-nbchan+1:end, :); % keep only ET data
ET.nbchan = nbchan;
ET.chanlocs = ET.chanlocs(end-nbchan+1:end);

%% Create optimized data for ICA training (OPTICAT, Dimigen, 2018)

% find name of saccade event 
i = find(~cellfun(@isempty, regexp({EEG.event(:).type}, 'sac')));
for idx = i
    EEG.event(idx).type = 'saccade';
end

% find name of fixation event 
i = find(~cellfun(@isempty, regexp({EEG.event(:).type}, 'fix')));
for idx = i
    EEG.event(idx).type = 'fixation';
end

% Overweight saccade intervals (containing spike potential)
EEG = pop_overweightevents(EEG,'saccade',ETparams.OPTICAT.SACCADE_WINDOW,ETparams.OPTICAT.OW_PROPORTION,ETparams.OPTICAT.REMOVE_EPOCHMEAN);

% Run ICA on optimized training data
fprintf('\nTraining ICA on the optimized data ...')


%% nested function for transformation of screen px to deg
function alpha_per_pix = px2deg(SCREEN_WIDTH_MM, SCREEN_RES_X, VIEW_DIST_MM)
    mm_per_pix    = SCREEN_WIDTH_MM/SCREEN_RES_X;
    alpha_per_pix = 180/pi*(2*atan(mm_per_pix/VIEW_DIST_MM/2));
end
end
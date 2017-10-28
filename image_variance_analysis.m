% image_variance_analysis
% Find regions of image data in which the time-variance of brightness is
% high enough to indicate significant movement
% EJR 2017
% License: CC-BY%
% 
% Method
% 1. Read in video data
% 2. Get user-defined time range for analysis
% 3. Compute time-variance of image data
% 4. Threshold the time-variance data to identify regions of movement


% 0.0 CONTROL PARAMETERS FOR THIS SCRIPT
filenameMP4 = 'Nov_07_Exp_5_H130D45.MP4'; % Name of file to process

% 0.1 PARTICLE TRACK ANALYSIS - SET TIME RANGE OF INTEREST
tInit   = 50;   % Timestamp in MP4 data (seconds) for first frame to analyse
tStep   = 2;    % Time step to get next frame to evaluate particle position
nSteps  = 6;   % Number of steps to consider 

% 0.1.1 PARTICLE TRACK ANALYSIS - GET USER CONFIRMATION OF TIME RANGE
prompt = {'time start (s)','time increment (s)','number of time steps'};
dlg_title = 'Please confirm time range for analysis';
num_lines = 1;
defaultans = { num2str(tInit),num2str(tStep),num2str(nSteps) };
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

tInit   = str2num(answer{1});
tStep   = str2num(answer{2});
nSteps  = str2num(answer{3});


% 1. INPUT
%    Create video reader object: 
v = VideoReader(filenameMP4);
vid_number_frames = v.Duration.*v.FrameRate;
vid_width = v.Width;
vid_height = v.Height;

% 2. ANALYSIS
% 2.1 Pre-allocate memory for image analysis
imDatGraySum        = zeros(vid_height,vid_width);
imDatGraySumSquares = zeros(vid_height,vid_width);
imDatGraySumMBG     = zeros(vid_height,vid_width);
imDatGraySumSquaresMBG = zeros(vid_height,vid_width);

for lpImDat =1:nSteps % For the required number of frames of MP4 data

v.CurrentTime = tInit + (lpImDat-1)*tStep;  % Find  desired timepoint
imDat     = readFrame(v);                   % Read this frame
imDatGray = mean(imDat, 3);

imDatGraySum        = imDatGraySum + imDatGray;
imDatGraySumSquares = imDatGraySumSquares + imDatGray.^2;

end

imDatGrayMean = imDatGraySum / vid_number_frames;

for lpImDat =1:nSteps % For the required number of frames of MP4 data

v.CurrentTime = tInit + (lpImDat-1)*tStep;  % Find  desired timepoint
imDat     = readFrame(v);                   % Read this frame
imDatGray = mean(imDat, 3) - imDatGrayMean;

imDatGraySumMBG        = imDatGraySumMBG + imDatGray;
imDatGraySumSquaresMBG = imDatGraySumSquaresMBG + imDatGray.^2;

end

% 
% imDat_time_variance = imDatGraySumSquares/vid_number_frames  - (imDatGraySum/vid_number_frames ).^2;
% 
% imDat_time_variance_n = imDat_time_variance ./ imDatGraySum;

figure(1)
imagesc(imDatGraySumSquaresMBG)
colorbar


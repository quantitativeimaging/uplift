% trackuplift
% Applies tracking methods to image data from uplift rig
% EJR 2016
% License: CC-BY
%
% Notes
% 1. Try to get uniform illumination - non-flat brightness affects pkfnd
%
% 2. The camera pointing direction seems to drift during aquisition, 
%    which introduced a constant translation. Needs immobilising. 
%    Or use a fiducial mark to counter drift - prefer immobilising.
%
% 3. Remove red lines on glass. These create detections at thfrac=0.4
% 
% 4. Note vertical stretch in image - measure graph paper to identify
%    this ratio (1.29 in some sample data), due to non-horizontal optic
%    axis, most likely.
%    2016/11/7: included YXratio to 'correct for' this - calculated
%    displacements in listDisps are now in 'horizontal pixel widths'
%
% 5. Co-ordinates
%    X is horizontal position (column number: leftmost =1)
%    Y is vertical position (row number: topmost = 1)
%
% 6. This script currently does not check whether the 
%    Time ranges defined under 0. (PARTICLE TRACK ANALYSIS) actually exist
%    in the MP4 video image data. If they don't, it will fail. 
%
% Method
% 0. This script uses a CONTROL / INPUT / ANALYIS / OUTPUT (visualisation)
% layout
% 1. Set scale manually for each video (in case camera angle shifts)
%    Use graph paper in image to check scale is uniform. 
%    Try to align camera so that one camera pixel width corresponds to 
%   the same distance horizontally and vertically on the testing rig. 
%   If this is not the case, then make sure you know the calibration.
% 2. 
% 

% ADD PATH FOR TRACKING FUNCTIONS WRITTEN BY JCC/DLB
addpath([pwd,'/tracking']);

% 0. CONTROL PARAMETERS FOR THIS SCRIPT
filenameMP4 = 'Nov_07_Exp_5_H130D45.MP4'; % Name of file to process

thfrac = 0.95;  % Fractional value to set threshold for particle finding
                % Note 'median' might be more robust than 0.25*(max-bg)

scaleX   = 30/108 ; % mm distance per horizontal pixel width
scaleY   = 30/144 ; % mm distance per vertical pixel width. 

YXratio = scaleX/scaleY; % How many vertical pixel widths equal one
                % horizontal width. Should be 1.00 if camera axis
                % horizontal, but in practice it may vary! Calibrate. 

% REGION OF INTEREST SELECTION
flagROIauto = 0;  % Set to 1 to bypass manual roi selection
roiBorder   = 50; % Width of border (pixels) where we will not try to interpolate displacement

flagShowAllTracks = 0;

% PARTICLE IDENTIFICATION
szPKFND = 11;    % Set to bigger than particle diameter.
szCNTRD = 11;    % Diameter of window for centroid finding. Avoid other particles.
lobjectBPASS = 5;% object diameter for bandpass filter. Try 5~diameter

flagBGgs = 1;   % Set to 1 to apply gaussian-blurred background subtration
radgauss = 20;  % Std deviation value for background subtraction

% PARTICLE TRACK IDENTIFICATION
maxdisp = 4;    % Maximum particle displacement per frame 
                % (distance is in pixel widths - somewhat wrongly assumed identical in X and Y) 

% PARTICLE TRACK ANALYSIS
tInit   = 42;   % Timestamp in MP4 data (seconds) for first frame to analyse
tStep   = 1;    % Time step to get next frame to evaluate particle position
nSteps  = 40;   % Number of steps to consider 

% MOVEMENT / STATIC region analysis
threshMove = 0.025*nSteps; % Threshold for identifying movement. 
                  % Is applied to dispacements in millimetres 
                  % Seems sensible to define as (speed X nSteps)

% 1. INPUT

v = VideoReader(filenameMP4);
numFrames = v.Duration.*v.FrameRate;

pos = [];       % Empty array to store identified positions

v.CurrentTime = tInit;
imDat = readFrame(v); % skip a few frames ahead
if(flagROIauto)
  roiRect = [480, 356, 550, 618]; % xmin, ymin, width, height
else
  figure(1)
  imagesc(imDat)
  title('Please select rectangular region of interest')
  roiRect = floor(getrect(1));
end


% 2. ANALYSIS
lpIm = 0;
for lpIm =1:nSteps % For the required number of frames of MP4 data

v.CurrentTime = tInit + (lpIm-1)*tStep;  % Find frame at desired timepoint
imDat = readFrame(v);                    % Read this frame

% The region of interest, imROI, should only contain material to be tracked
% Try finding mean RGB value level (could use just red, or rgb2gray)
imROI = mean(imcrop(imDat,roiRect) ,3);
% figure(2)
% imagesc(imROI)
% colormap(gray)

if(flagBGgs)
  imBG = imgaussfilt(imROI, radgauss);
  im1 = double(imBG) - double(imROI);
else
  im1 = 255-imROI;
end
    
figure(4)
imagesc(im1)
  colormap(gray)
  title('Unfiltered region of interest')

im2 = bpass(im1, 1, lobjectBPASS);
% figure(3)
% imagesc(im2)
%   colormap(gray)
%   title('Filtered region of interest')

if lpIm == 1 % set threshold using data in frame 1
  % thresh = min(im2(:)) + thfrac*(max(im2(:)) - min(im2(:)));
  thresh = quantile(im2(:), thfrac);
end

% APPLY JCC / DLB peak- and centroid-finding
pk  = pkfnd(im2,thresh,11);
cnt = cntrd(im2,pk,szCNTRD);

figure(5)
imagesc(im2)
colormap('gray')
hold on
 scatter(cnt(:,1),cnt(:,2), 'r')
hold off

pause(0.1)
% listThresh = [listThresh;thresh];
newpos = [cnt(:,1), cnt(:,2), lpIm*ones(size(cnt,1),1)];

pos = [pos;newpos];
end

%% Apply tracking method of JCC/DLB
res = track(pos,maxdisp);

% Plot all identified tracks
if(flagShowAllTracks)
 figure(3)
 imagesc(imROI)
 colormap('gray')
 hold on
 for lpT = 1:res(end,4)
  resA = res(res(:,4)==lpT,:);
  plot(resA(:,1),resA(:,2) ,'r');
  scatter(resA(1,1), resA(1,2), 'c');
  scatter(resA(end,1), resA(end,2), 'r');
 end
 legend('track','start','end')
 hold off
end

%% 3. OUTPUT / VISUALISATION
% --------------------------------------------
% POST PROCESSING FOR VISUALISATION
% % Now find initial and displaced positions,
% Find particles that were tracked from frame 1 up to nSteps (all frames)
numTracks = (res(end,4));
posInit = -ones(numTracks,2);
posDisp = -ones(numTracks,2);
for lpTrack = 1:numTracks
   
    myCoords = res((res(:,4)==lpTrack),[1:4]);
    if( ((sum(myCoords(:,3)==1))+(sum(myCoords(:,3)==nSteps)))==2 )
    posInit(lpTrack,:) = myCoords((myCoords(:,3)==1),1:2);
    posDisp(lpTrack,:) = myCoords((myCoords(:,3)==nSteps),1:2);
    end
end

invalid = (posInit(:,1)==-1) | (posDisp(:,1)==-1);

posInit(invalid,:) = [];
posDisp(invalid,:) = [];

figure(10)
imagesc(imROI)
colormap(gray)
hold on
  scatter(posInit(:,1), posInit(:,2), 'c');
  scatter(posDisp(:,1), posDisp(:,2), 'r');
hold off
title('Identified positions of real particles');
legend('frame 1', ['frame', int2str(nSteps)])
axis equal
xlabel('X-position, pixels')
ylabel('Y-position, pixels')
set(gca, 'fontSize', 14)

% ---------------------------- CP 2 TFORM
mytform = cp2tform(posDisp, posInit, 'lwm');
% 'piecewise linear' is stricter than local weighted mean 'lwm'

[XX,YY] = meshgrid([roiBorder:25:(roiRect(3)-roiBorder)], ...
                   [roiBorder:25:(roiRect(4)-roiBorder)]);

listXinit = XX(:);
listYinit = YY(:);

[listXfinal, listYfinal] = tforminv(mytform, listXinit, listYinit);

listDisps = sqrt((listXfinal - listXinit).^2 +...
            ((listYfinal - listYinit)/YXratio).^2); % YXratio for vstretch
matrDisps = reshape(listDisps, size(XX));

% % Scatterplot just showing initial and displaced positions, no image
% figure(11)
% scatter(listXinit, -listYinit, 'r')
% hold on
%   scatter(listXfinal, -listYfinal, 'c')
% hold off

figure(12)
imagesc(imROI)
colormap(gray)
hold on
  scatter(listXinit, listYinit, 'r')
  scatter(listXfinal, listYfinal, 'c')
hold off
axis equal
legend('frame 1', ['frame', int2str(nSteps)])
title('Interpolated positions');
xlabel('X-position, pixels')
ylabel('Y-position, pixels')
set(gca, 'fontSize', 14)

figure(13)
imagesc(imROI)
colormap(gray)
hold on
for lpPt = 1:size(listXinit,1)
   plot([listXinit(lpPt);listXfinal(lpPt)],...
        [listYinit(lpPt);listYfinal(lpPt)], 'r', 'lineWidth', 1) 
end
hold off
axis equal
title(['Displacement during ', int2str(nSteps), ' frames']);
xlabel('X-position, pixels')
ylabel('Y-position, pixels')
set(gca, 'fontSize', 14)

figure(14)
mesh(XX, YY, matrDisps)
colormap(jet)
xlabel('X position, pixels', 'fontSize', 14)
ylabel('Y position, pixels', 'fontSize', 14)
zlabel('speed, pixel per 9 frames', 'fontSize', 14)
set(gca, 'fontSize', 12)

% IDENTIFY MOVEMENT / NON-MOVED REGIONS
% % Next section of analysis - threshold the displacement field. 
% Try overlaying thresholded slip field to identiy regions where
% speed exceeds some specific value
% First calculate displacement field at 1:1 scale
[XX,YY] = meshgrid([roiBorder:1:(roiRect(3)-roiBorder)], ...
                   [roiBorder:1:(roiRect(4)-roiBorder)]);
listXinit = XX(:);
listYinit = YY(:);
[listXfinal, listYfinal] = tforminv(mytform, listXinit, listYinit);
listDisps = sqrt( (scaleX*(listXfinal - listXinit)).^2 +...
                  (scaleY*(listYfinal - listYinit)).^2); % In millimetres
matrDisps = reshape(listDisps, size(XX));

% Get displacements separately in X- and Y-
listDispX = listXfinal - listXinit;
listDispY = (listYfinal - listYinit);

maskSlip = (matrDisps>threshMove);

figure(15)
  imagesc(maskSlip)
  title(['Region of movement > ', num2str(threshMove), ' mm'])
  xlabel('X position, pixels', 'fontSize', 14)
  ylabel('Y position, pixels', 'fontSize', 14)

imOverlay = double(imDat)/255;
imOverlay( roiRect(2)-0+[roiBorder:1:(roiRect(4)-roiBorder)], ...
           roiRect(1)-0+[roiBorder:1:(roiRect(3)-roiBorder)],2) = ...
    0.5*imOverlay( roiRect(2)-0+[roiBorder:1:(roiRect(4)-roiBorder)], ...
                   roiRect(1)-0+[roiBorder:1:(roiRect(3)-roiBorder)],2)+...
    0.5*double(maskSlip);
figure(16)
imagesc(imOverlay)

% NEED TO FIX PROFILE PLOTS:
% figure(17)
% plot(listXinit(listYinit==500)*scaleX, listDisps(listYinit==500)*scaleX,...
%      'b','lineWidth', 2);
% hold on
%   plot(listXinit(listYinit==100)*scaleX, listDisps(listYinit==100)*scaleX,...
%        'r', 'lineWidth', 2);
%   
%   plot([roiXplate(1), roiXplate(1)]*scaleX, [0,4], 'k--', 'lineWidth', 1)
%   plot([roiXplate(2), roiXplate(2)]*scaleX, [0,4], 'k--', 'lineWidth', 1)
% hold off
% legend('Near base','Near top', 'Plate edges')
% set(gca,'fontSize', 14)
% xlabel('Horizontal position / mm')
% ylabel('Displacement during 18 seconds / mm')
% grid on

% % Try integrating speedss - change this to integrating vertical velocities
% figure(18)
% matrDispY = reshape(listDispY, size(YY));
% plot(-sum(matrDispY,2))
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
%    Or use a fiducial mark to counter drift - although less ideal.
%
% 3. Remove red lines on glass. These create detections at thfrac=0.4
% 
% 4. Note vertical stretch in image - measure graph paper to identify
%    this ratio (1.29 in some sample data), due to non-horizontal optic
%    axis, most likely.
%    2016/11/7: included YXratio to 'correct for' this - calculated
%    displacements in listDisps are now in 'horizontal pixel widths'
%

addpath([pwd,'/tracking']);

thfrac = 0.50;

flagBGgs = 1;
radgauss = 20;

YXratio = 1.29; % From image data, how many vertical pixel widths equal one
                % horizontal width. Should be 1.00 if camera axis
                % horizontal, but in practice it may vary! Calibrate. 
hScale  = 40/149 ;   % mm per horizontal pixel width
                
% 1. Input

v = VideoReader('Exp1.MP4');

numFrames = v.Duration.*v.FrameRate;

pos = [];
% listThresh =[];
v.CurrentTime = 20;
lpIm = 0;
while lpIm < 12 %hasFrame(v)
lpIm = lpIm + 1;

v.CurrentTime = 20 + lpIm*2;

imDat = readFrame(v); % skip a few frames ahead
% figure(1)
% imshow(imDat)

% Select a region containing only particles
% Try finding mean grey level (could use red channel)
imROI = mean( imDat(330:(330+689), 370:(370+759), :), 3); 
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

im2 = bpass(im1, 1, 5);
figure(5)
imagesc(im2)

if lpIm == 1 % set threshold on frame 1
  thresh = min(im2(:)) + thfrac*(max(im2(:)) - min(im2(:)));
end

pk = pkfnd(im2,thresh,11);

figure(5)
imagesc(im2)
colormap('gray')
hold on
 scatter(pk(:,1),pk(:,2), 'r')
hold off

cnt = cntrd(im2,pk,15);

pause(1)
% listThresh = [listThresh;thresh];
newpos = [cnt(:,1), cnt(:,2), lpIm*ones(size(cnt,1),1)];

pos = [pos;newpos];
end

% Apply tracking method:
res = track(pos,6);

% plot
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


% --------------------------------------------
% % Now find initial and displaced positions,
% Try arbitrary frames (1 and 10) at first:
numTracks = (res(end,4));
posInit = -ones(numTracks,2);
posDisp = -ones(numTracks,2);
for lpTrack = 2:numTracks
   
    myCoords = res((res(:,4)==lpTrack),[1:4]);
    if( ((sum(myCoords(:,3)==1))+(sum(myCoords(:,3)==10)))==2 )
    posInit(lpTrack,:) = myCoords((myCoords(:,3)==1),1:2);
    posDisp(lpTrack,:) = myCoords((myCoords(:,3)==10),1:2);
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
legend('frame 1', 'frame 10')
axis equal

% ---------------------------- CP 2 TFORM
mytform = cp2tform(posDisp, posInit, 'lwm');
% 'piecewise linear' is stricter than local weighted mean 'lwm'

[XX,YY] = meshgrid([100:25:600], [100:25:600]);

listXinit = XX(:);
listYinit = YY(:);

[listXfinal, listYfinal] = tforminv(mytform, listXinit, listYinit);

listDisps = sqrt((listXfinal - listXinit).^2 +...
            ((listYfinal - listYinit)/YXratio).^2); % YXratio for vstretch
matrDisps = reshape(listDisps, size(XX));

figure(11)
scatter(listXinit, -listYinit, 'r')
hold on
  scatter(listXfinal, -listYfinal, 'c')
hold off

figure(12)
imagesc(imROI)
hold on
  scatter(listXinit, listYinit, 'r')
  scatter(listXfinal, listYfinal, 'c')
hold off
axis equal

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

figure(14)
mesh(XX, YY, matrDisps)
colormap(jet)
xlabel('X position, pixels', 'fontSize', 14)
ylabel('Y position, pixels', 'fontSize', 14)
zlabel('speed, pixel per 9 frames', 'fontSize', 14)
set(gca, 'fontSize', 12)

% % Next section of analysis - threshold the displacement field. 
% Try overlaying thresholded slip field to identiy regions where
% speed exceeds some specific value
% First calculate displacement field at 1:1 scale
[XX,YY] = meshgrid([100:1:600], [100:1:600]);
listXinit = XX(:);
listYinit = YY(:);
[listXfinal, listYfinal] = tforminv(mytform, listXinit, listYinit);
listDisps = sqrt((listXfinal - listXinit).^2 +...
            ((listYfinal - listYinit)/YXratio).^2); % YXratio for vstretch
matrDisps = reshape(listDisps, size(XX));

% Get displacements separately in X- and Y-
listDispX = listXfinal - listXinit;
listDispY = (listYfinal - listYinit);

figure(15)
maskSlip = (matrDisps>1);
imagesc(maskSlip)

imOverlay = double(imDat)/255;
imOverlay(430:(430+500), 470:(470+500),2) = ...
    0.5*(imOverlay(430:(430+500), 470:(470+500),2)) + 0.5*double(maskSlip);
figure(16)
imagesc(imOverlay)

figure(17)
plot(listXinit(listYinit==500)*hScale, listDisps(listYinit==500)*hScale,...
     'b','lineWidth', 2);
hold on
  plot(listXinit(listYinit==100)*hScale, listDisps(listYinit==100)*hScale,...
       'r', 'lineWidth', 2);
  
  plot([(560-370), (560-370)]*hScale, [0,4], 'k--', 'lineWidth', 1)
  plot([(795-370), (795-370)]*hScale, [0,4], 'k--', 'lineWidth', 1)
hold off
legend('Near base','Near top', 'Plate edges')
set(gca,'fontSize', 14)
xlabel('Horizontal position / mm')
ylabel('Displacement during 18 seconds / mm')
grid on

% Try integrating speedss - change this to integrating vertical velocities
figure(18)
matrDispY = reshape(listDispY, size(YY));
plot(-sum(matrDispY,2))
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

addpath([pwd,'/tracking']);

thfrac = 0.50;

flagBGgs = 1;
radgauss = 20;

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

% Select a region containing only particles
% Try finding mean grey level (could use red channel)
imROI = mean( imDat(330:(330+689), 370:(370+759), :), 3); 
% imagesc(imROI)

if(flagBGgs)
  imBG = imgaussfilt(imROI, radgauss);
  im1 = double(imBG) - double(imROI);
else
  im1 = 255-imROI;
end
    
figure(4)
imagesc(im1)

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


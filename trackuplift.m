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


addpath([pwd,'/tracking']);

thfrac = 0.70;

% 1. Input

v = VideoReader('Exp1.MP4');

numFrames = v.Duration.*v.FrameRate;

pos = [];
% listThresh =[];
lpIm = 0;
while lpIm < 25 %hasFrame(v)
lpIm = lpIm + 1;

v.CurrentTime = lpIm*2;

imDat = readFrame(v); % skip a few frames ahead

% Select a region containing only particles
% select the red channel so the lines on the glass are fainter
imROI = imDat(330:(330+689), 370:(370+759) , 1); 

figure(3)
imagesc(imROI)

im1 = 255-imROI;

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
hold on
 for lpT = 1:res(end,4)
 plot(res(res(:,4)==lpT,1),res(res(:,4)==lpT,2) ,'r')
 end
hold off


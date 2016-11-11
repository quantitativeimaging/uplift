% timeresolved_uplift
% EJR 2016 CC-BY
%
% Notes
%   This is a draft method for time-resolving particle uplift.
%   First run trackuplift_v2
%   Then run this script

outDirP = '\outVidP\';
outDirO = '\outVidO\';

framesPerSlice = 5;
nSlices = floor( nSteps / framesPerSlice ) ;

threshMove = 0.025*framesPerSlice;

for lpSlice = 1:nSlices
frameStart = 1+(lpSlice-1)*framesPerSlice;
frameEnd   = lpSlice*framesPerSlice;

v.CurrentTime = tInit + (frameEnd-1)*tStep;  % Find frame at desired timepoint
imDat = readFrame(v);                    % Read this frame
imROI = mean(imcrop(imDat,roiRect) ,3);


% Identify initial and displaced points from reliable tracks
numTracks = (res(end,4));
posInit = -ones(numTracks,2);
posDisp = -ones(numTracks,2);
for lpTrack = 1:numTracks
    myCoords = res((res(:,4)==lpTrack),[1:4]);
    if( ((sum(myCoords(:,3)==1))+ ...
         (sum(myCoords(:,3)==frameStart))+ ...
         (sum(myCoords(:,3)==frameEnd))) ==3 )
    posInit(lpTrack,:) = myCoords((myCoords(:,3)==frameStart),1:2);
    posDisp(lpTrack,:) = myCoords((myCoords(:,3)==frameEnd),1:2);
    end
end
invalid = (posInit(:,1)==-1) | (posDisp(:,1)==-1);
posInit(invalid,:) = [];
posDisp(invalid,:) = [];

% Plot identified displacements in current slice
figure(10)
imagesc(imROI)
colormap(gray)
hold on
  scatter(posInit(:,1), posInit(:,2), 'c');
  scatter(posDisp(:,1), posDisp(:,2), 'r');
hold off
title(['Identified positions of real particles in slice',num2str(lpSlice)]);
legend(['frame ', int2str(frameStart)], ['frame', int2str(frameEnd)])
axis equal
xlabel('X-position, pixels')
ylabel('Y-position, pixels')
set(gca, 'fontSize', 14)

myFrame = getframe(gcf);
myFrDat = myFrame.cdata;
imwrite(myFrDat, ['D:\EJR_GIT\reverse_hopper\outVidP\part',int2str(lpSlice),'.png']);

% Identify current displacement field in this slice
% ---------------------------- CP 2 TFORM
mytform = cp2tform(posDisp, posInit, 'lwm');
% 'piecewise linear' is stricter than local weighted mean 'lwm'


% Generate slip/static overlay image:
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

imOverlay = double(imDat)/255;
imOverlay( roiRect(2)-1+[roiBorder:1:(roiRect(4)-roiBorder)], ...
           roiRect(1)-1+[roiBorder:1:(roiRect(3)-roiBorder)],2) = ...
    0.5*imOverlay( roiRect(2)-1+[roiBorder:1:(roiRect(4)-roiBorder)], ...
                   roiRect(1)-1+[roiBorder:1:(roiRect(3)-roiBorder)],2)+...
    0.5*double(maskSlip);
figure(16)
imagesc(imOverlay)
    
drawnow;
myFrame = getframe(gcf);
myFrDat = myFrame.cdata;
imwrite(myFrDat, ['D:\EJR_GIT\reverse_hopper\outVidO\over',int2str(lpSlice),'.png']);

end


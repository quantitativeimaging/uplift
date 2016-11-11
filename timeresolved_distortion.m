% timeresolved_distortion
% EJR 2016, CC-BY
%
% Generates a map showing final state of an initially-square grid

outDirO = '\outVidM\';

framesPerSlice = 5;
nSlices = floor( nSteps / framesPerSlice ) ;

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
         (sum(myCoords(:,3)==frameEnd))) == 2 )
    posInit(lpTrack,:) = myCoords((myCoords(:,3)==1),1:2);
    posDisp(lpTrack,:) = myCoords((myCoords(:,3)==frameEnd),1:2);
    end
  end
  invalid = (posInit(:,1)==-1) | (posDisp(:,1)==-1);
  posInit(invalid,:) = [];
  posDisp(invalid,:) = [];
  
  
  % Identify current displacement field in this slice
  % ---------------------------- CP 2 TFORM
  mytform = cp2tform(posDisp, posInit, 'lwm');
  % 'piecewise linear' is stricter than local weighted mean 'lwm'


  % Generate slip/static overlay image:
  [XX,YY] = meshgrid([roiBorder:25:(roiRect(3)-roiBorder)], ...
                   [roiBorder:25:(roiRect(4)-roiBorder)]);
  listXinit = XX(:);
  listYinit = YY(:);
  [listXfinal, listYfinal] = tforminv(mytform, listXinit, listYinit);
  listDisps = sqrt( (scaleX*(listXfinal - listXinit)).^2 +...
                    (scaleY*(listYfinal - listYinit)).^2); % In millimetres

  
  XXfinal = reshape(listXfinal, size(XX));
  YYfinal = reshape(listYfinal, size(YY));

  figure(20)
  plot(XX,YY,'c');
  set(gca,'YDir','reverse');
  hold on
    plot(XX',YY','c')
    plot(XXfinal, YYfinal, 'r')
    plot(XXfinal',YYfinal','r')
  hold off
  title(['Distortion up to slice ',num2str(lpSlice)]);
  axis equal
  xlabel('X-position, pixels')
  ylabel('Y-position, pixels')
  set(gca, 'fontSize', 14)
  xlim([0 700])
  ylim([0 600])
  drawnow()
  
  
  myFrame = getframe(gcf);
  myFrDat = myFrame.cdata;
  imwrite(myFrDat, ['D:\EJR_GIT\reverse_hopper\outVidM\mesh',int2str(lpSlice),'.png']);

end
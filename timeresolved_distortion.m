% timeresolved_distortion
% EJR 2016, CC-BY
%
% NOTES
% -Generates a map showing final state of an initially-square grid
%
% -Also tries to calculate area strain at each (non-edge) vertex of mesh

outName = 'D:\EJR_GIT\reverse_hopper\outVidM\mesh'; % File for output

framesPerSlice = 5;
nSlices = floor( nSteps / framesPerSlice ) ;

[XX,YY] = meshgrid([roiBorder:25:(roiRect(3)-roiBorder)], ...
                   [roiBorder:25:(roiRect(4)-roiBorder)]);
arrayA = zeros([size(XX)-2,nSlices]);

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
  
  
  % Identify current displacement field at end of this slice
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
  
  % Capture and save the distored mesh
  myFrame = getframe(gcf);
  myFrDat = myFrame.cdata;
  imwrite(myFrDat, [outName,int2str(lpSlice),'.png']);

  % Try to calculate area strain at each vertex of current mesh
  % (a) Evaluate two current displacement vectors at each vertex, 
  %      for vectors that were initially 'right' and 'down' in start grid
  %      size = (rows, cols, 2) - add third direction (Dz = 0) for 'cross'
  vectH = -ones([size(XXfinal),3]); % Residual -ones should be ignored
  vectH(1:end,(1:end-1),1) = XXfinal(1:end,2:end)-XXfinal(1:end,(1:end-1));
  vectH(1:end,(1:end-1),2) = YYfinal(1:end,2:end)-YYfinal(1:end,(1:end-1));
  vectH(:,:,3) = 0; % Rel. displacement in z-direction is taken to be zero
  % Apply pixel width to get physical displacements
  vectH(:,:,1) = vectH(:,:,1)*scaleX;
  vectH(:,:,2) = vectH(:,:,2)*scaleY;
   
  vectV = -ones([size(YYfinal),3]);
  vectV(1:(end-1),1:end,1) = XXfinal(2:end,1:end)-XXfinal(1:(end-1),1:end);
  vectV(1:(end-1),1:end,2) = YYfinal(2:end,1:end)-YYfinal(1:(end-1),1:end);
  vectV(:,:,3) = 0;
  % Apply pixel width scale to get physical displacements
  vectV(:,:,1) = vectV(:,:,1)*scaleX;
  vectV(:,:,2) = vectV(:,:,2)*scaleY;
  
  % Calculate areas of elements
  vectAreas = abs(cross(vectH,vectV,3));
  AA = vectAreas(:,:,3);
  
  matrArea = zeros(size(XXfinal));
  matrArea(2:(end-1),2:(end-1)) = 0.5.*( AA(2:(end-1),2:(end-1) ) + ...
                                        AA(1:(end-2),1:(end-2)) + ...
                                        AA(1:(end-2),2:(end-1)    ) + ... 
                                        AA(2:(end-1),1:(end-2)) );
  if(lpSlice ==1)
    matrAreaInit = matrArea;
  end
  arrayA(:,:,lpSlice) = matrArea(2:(end-1),2:(end-1));
  
  figure(21)
  % imagesc(matrArea(2:(end-1),2:(end-1)))
  imagesc(matrArea(2:(end-1),2:(end-1))./matrAreaInit(2:(end-1),2:(end-1)))
  % view([180 90 ])
  caxis([0.9 1.1])
  colorbar
   title(['Area strain of elements up to slice ',num2str(lpSlice)]);
  axis equal
  xlabel('X-position, grid points')
  ylabel('Y-position, grid points')
  set(gca, 'fontSize', 14)
  
  % Calculate aspect ratios of distorted elements
end

% Display (smoothed) area strain at the end of some time-slice
BB = arrayA(:,:,4)./arrayA(:,:,1); % area strain
CC = imgaussfilt(BB-1, 1) +1; % slight smoothing - just for visualisation 
figure(22)
  imagesc(CC)
  caxis([1 1.05])
xlabel('X-position, grid points')
ylabel('Y-position, grid points')
set(gca, 'fontSize', 14)
colorbar
title('(Smoothed) area strain... dilation in shear zone')

% Try integrating energy corresponding to uplift and expansion... 

  
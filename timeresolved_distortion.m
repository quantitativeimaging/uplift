% timeresolved_distortion
% EJR 2016, CC-BY
%
% NOTES. This script:
% 1. Post-processes particle tracking (velocimetry) results from
% trackuplift_v2
%
% 2. Generates a map showing final state of an initially-square grid
%
% 3. Also tries to calculate area strain at each (non-edge) vertex of mesh
%
% 4. Then try to calculate potential energy for
% -> gravitational potential (uplift)
% -> volume expantion (try initially as: ideal gas?)
%
% 5. Also tries to calculate the shear rate (as a tensor?) for each area
% element in the initial mesh
%
% 6. Hypothesis: area strain is proportional to shear rate, providing
% enough info to make a Lagrange mechanics model of uplift assuming PdV
% work and gravitational potential are the key factors. (Also: that area
% strain is all vertical due to confinement.)
% 
% 7: Warning. Using zero border for mesh interpolation may fail due to lack
% of well-tracked particles (moving outside ROI) at large time points.
%
% 8. Warning. Pressure interpolation is codged from grid size, and should 
% be re-done with an exact initial position for the top of the ballotini.
% Also, the grid spacing should be changed so that the entire area of 
% distored material is always captured exactly. 
%
% 9. Try: particle tracking in larger ROI, with throw-away for tracks that
%  start outside the area of ballotini in frame 1 (i.e. ignore air).
%
% 10. Try tracking video frames at 10 Hz and taking overlapping 1 sec
% slices to get repeat estimates of area strain for checking / smoothing
%
% 11. Warning. Grid distortion measurement at 50 px (1/3 plate) is 
% sensitive to x/y edges. Need small enough cells to capture plate edge,
% but big enough to be free of Nyquist (?) like problems. 

outName = 'D:\EJR_GIT\reverse_hopper\outVidM\mesh'; % File for output

roiBorder = 8; % Unfortunately, strain is important at very grid bottom!
gridSpacing = 25;
% gridSpacing = (roiRect(4)-2*roiBorder) /12; % Capture full vertical range

initialHeight = 345; % Row in imDat corresponding to top of ballotini.

framesPerSlice = 1;
nSlices = 20;
% nSlices = floor( nSteps / framesPerSlice ) ;

[XX,YY] = meshgrid([roiBorder:gridSpacing:(roiRect(3)-roiBorder)], ...
                   [roiBorder:gridSpacing:(roiRect(4)-roiBorder)]);

% Array to store areas corresponding to non-boundary vertices in XX,YY
arrayA = zeros([size(XX)-1,nSlices]);

% Array of hydrostatic pressures estimated at centre of each 'square'
arrayPre = ( roiRect(2) - initialHeight + gridSpacing/2 + ...
            YY(1:(end-1),1:(end-1)) ).*scaleY*0.001*9.81*1800;


% Array to store estimated gravitational and strain energies
arrayGPE = zeros(size(arrayA)); % Let slice 1 = zero for each element
arrayEPE = zeros(size(arrayA)); % Let slice 1 = zero for each element
listZs = zeros(nSlices, 1);

for lpSlice = 1:nSlices
  frameStart = 1+(lpSlice-1)*framesPerSlice;
  frameEnd   = lpSlice*framesPerSlice;
  frameMax   = nSlices*framesPerSlice;
    
  v.CurrentTime = tInit + (frameEnd-1)*tStep;  % Find frame at desired timepoint
  imDat = readFrame(v);                    % Read this frame
  imROI = mean(imcrop(imDat,roiRect) ,3);
  % imagesc(imDat)

  % Identify initial and displaced points from reliable tracks
  numTracks = (res(end,4));
  posInit = -ones(numTracks,2);
  posDisp = -ones(numTracks,2);
  for lpTrack = 1:numTracks
    myCoords = res((res(:,4)==lpTrack),[1:4]);
    if( ((sum(myCoords(:,3)==1))+ ...
         (sum(myCoords(:,3)==frameEnd)) + ...
         (sum(myCoords(:,3)==frameMax))) == 3 )
    posInit(lpTrack,:) = myCoords((myCoords(:,3)==1),1:2);
    posDisp(lpTrack,:) = myCoords((myCoords(:,3)==frameEnd),1:2);
    end
  end
  invalid = (posInit(:,1)==-1) | (posDisp(:,1)==-1);
  posInit(invalid,:) = [];
  posDisp(invalid,:) = [];
  
  
  % Identify current displacement field at end of this slice
  % ---------------------------- CP 2 TFORM
  % mytform = cp2tform(posDisp, posInit, 'lwm');
  % 'piecewise linear' is stricter than local weighted mean 'lwm'
  % 2016-11-23: change to use fitgeotrans
  %   Note n can be as small as 6, but risks ill-conditioning
  %   In cp2tform, n defaults to 12. Higher values somewhat help to avoid
  % bad extrapolations at grid edge.
  n = 14;
  mytform = fitgeotrans(posDisp,posInit,'lwm',n);
  % mytform = fitgeotrans(posDisp,posInit,'pwl');
  
  % Generate slip/static overlay image:
  [XX,YY] = meshgrid([roiBorder:gridSpacing:(roiRect(3)-roiBorder)], ...
                   [roiBorder:gridSpacing:(roiRect(4)-roiBorder)]);
  listXinit = XX(:);
  listYinit = YY(:);
  % [listXfinal, listYfinal] = tforminv(mytform, listXinit, listYinit); %
  % for co2tform based method
  [listXfinal, listYfinal] = transformPointsInverse(mytform, listXinit, listYinit);
  listDisps = sqrt( (scaleX*(listXfinal - listXinit)).^2 +...
                    (scaleY*(listYfinal - listYinit)).^2); % In millimetres
  % matrDisps = reshape(listDisps, size(XX));
  
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
  xlim([-50 700])
  ylim([-50 650])
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
  matrArea(1:(end-1), 1:(end-1)) = 1.0 * AA(1:(end-1), 1:(end-1));
%   matrArea(2:(end-1),2:(end-1)) = 0.25.*( AA(2:(end-1),2:(end-1) ) + ...
%                                         AA(1:(end-2),1:(end-2)) + ...
%                                         AA(1:(end-2),2:(end-1)    ) + ... 
%                                         AA(2:(end-1),1:(end-2)) );
%                                     % 0.25 gets area nearest each vertex?
  if(lpSlice ==1)
    matrAreaInit = matrArea;
  end
  arrayA(:,:,lpSlice) = matrArea(1:(end-1), 1:(end-1));
  
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
  
  % Optional: calculate energies of the interpolated grid
  volOfEle = gridSpacing*scaleX*gridSpacing*scaleY*0.001*0.001*0.024;
  arrayGPE(:,:,lpSlice) = -( YYfinal(1:(end-1),1:(end-1)) - YY(1:(end-1),1:(end-1)) )...
                          * scaleY * 0.001 * 9.81 * 1800 * volOfEle; 
         % mm per px * 0.001, g, 1800 kg/m^3, 24 mm depth, width
         % negative as rows are counted downwards
  arrayAst = arrayA(:,:,lpSlice)./arrayA(:,:,1) - 1;
  arrayAst(arrayAst<0)=0;
  arrayEPE(:,:,lpSlice) = volOfEle * arrayAst.*arrayPre;
  
  listZs(lpSlice) = -YYfinal(end,floor(end/2))*0.001*scaleY;

  figure(27) % Overlay inferred mesh on image data
  imagesc(imDat)
  hold on
    plot(XX+roiRect(1),YY+roiRect(2),'c');
    plot(XX'+roiRect(1),YY'+roiRect(2),'c')
    plot(XXfinal+roiRect(1), YYfinal+roiRect(2), 'r')
    plot(XXfinal'+roiRect(1),YYfinal'+roiRect(2),'r')
  hold off

end

% % Display (smoothed) area strain at the end of some time-slice
BB = arrayA(:,:,10)./arrayA(:,:,1); % area strain
CC = imgaussfilt(BB-1, 1) +1; % slight smoothing - just for visualisation 
figure(22)
  imagesc(BB)
  caxis([1 1.05])
xlabel('X-position, grid points')
ylabel('Y-position, grid points')
set(gca, 'fontSize', 14)
colorbar
title('(Smoothed) area strain... dilation in shear zone')

% Try including the following in the loop:
% Try integrating energy corresponding to uplift and expansion... 
% Assume:
%  Slice 1 has zero GPE (height energy) and EPE (strain) in each element
%  Only calculate energy from increasing height (ignore falling)
%  Only estimate energy from expansion Vs hydrostatic pressure ('ideal
%  gas')
%  Neglect energy of compressive strain, if there is any compression
%  Try to neglect the 'flow around' region. (Small volume - may be OK.) 
% arrayGPE = zeros(size(arrayA)); % Let slice 1 = zero for each element
% arrayEPE = zeros(size(arrayA)); % Let slice 1 = zero for each element
figure(24)
imagesc(arrayGPE(:,:,15) - arrayGPE(:,:,3))
 caxis([0 60E-5])
xlabel('X-position, grid points')
ylabel('Y-position, grid points')
set(gca, 'fontSize', 14)
colorbar
title('Observed GPE change / J per element , XX tsec')

figure(25)
imagesc(arrayEPE(:,:,15) - arrayEPE(:,:,3))
 caxis([0 60E-5])
xlabel('X-position, grid points')
ylabel('Y-position, grid points')
set(gca, 'fontSize', 14)
colorbar
title('Estimated p dV work/ J per element , XX sec')

% Try summing up the estimated energies.
listGPEs = squeeze(sum(sum(arrayGPE)));
listEPEs = squeeze(sum(sum(arrayEPE)));
listTEs  = listGPEs + listEPEs;

listW1s = listGPEs(2:end) - listGPEs(1:(end-1));
listW1s(listW1s<0) = 0; % assume any relaxation is wasted as heat
% listW2s = listEPEs(2:end) - listEPEs(1:(end-1));
% alternative
deltEPE =  arrayEPE(:,:,(2:end))-arrayEPE(:,:,1:(end-1));
%deltEPE(deltEPE<0)=0; % make choice to integrate only +ve PdV work?
                       % But this would sum up +ve values from noise
listW2s = squeeze(sum(sum(deltEPE,1),2)) ;
listW2s(listW2s<0) = 0; % assume any relaxation is wasted as heat
listWs = listW1s + listW2s;

% listDzs = ones(size(listZs));
listDzs = listZs(2:end) - listZs(1:(end-1));

% From visual analysis:
% Alternative: for 42 s to 62 s in sample data: inspect pixel posit of plt:
listZplate=[5,5,5,5,5,6,7,8,9.5,11,13,15,18,20,22.5,25,28,31,34,37]'*scaleY*0.001;
listDzs = listZplate(2:end) - listZplate(1:(end-1));

listFs = listWs./listDzs ; % at 1.2 mm / 5s

figure(26)
plot(listFs)
  hold on
    plot(listW1s./listDzs)
    plot(listW2s./listDzs)
    plot([0 20], [2.97, 2.97], 'k--')
  hold off
title('Force corresponding to imaged energy');
legend('total force', 'mg dz work','PdV work', 'column weight')
xlim([5 20])
xlabel('time / seconds')
ylabel('Force / N')
set(gca, 'fontSize', 14)



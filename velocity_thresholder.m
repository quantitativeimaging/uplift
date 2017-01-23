threshMove = 0.025*7*1.8; % Threshold for identifying movement. 

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
% extrapolate_cone_apex
%
% Method
%   Run trackuplift_v2
%   Then run this script. It will try to extrapolate cone apex position.
%
% Note: dimensions are in pixels
%       So: to get alpha, you need 
%       alpha = atan( scaleX*(D/2) / (scaleY * Y_A)  )
%       where
% scaleX   =   mm distance per horizontal pixel width
% scaleY   =   mm distance per vertical pixel width. 
%
% Suggest you: 
%  1. Read value of D for the experiment (e.g. 45 mm)
%  2. Measure y-coordinate of plate (e.g. 620)
%  3. Read y-co-ordinate fitted to cone apex (e.g. apex_position(2) )
%     apex_position(2)
%  4. Get distance of cone apex below plate
%     Y_A = ( apex_position(2) - y_coordinate_plate ) * scaleY
%  5. Get alpha using:
%     atan( (D / 2) /Y_A )  (radians)
%     atan( (D / 2) /Y_A ) * 180/pi

% The following code plots interpolated positions
[XX,YY] = meshgrid([roiBorder:25:(roiRect(3)-roiBorder)], ...
                   [roiBorder:25:(roiRect(4)-roiBorder)]);
listXinit = XX(:);
listYinit = YY(:);

[listXfinal, listYfinal] = tforminv(mytform, listXinit, listYinit);
listDisps = sqrt((listXfinal - listXinit).^2 +...
            ((listYfinal - listYinit)/YXratio).^2); % YXratio for vstretch
matrDisps = reshape(listDisps, size(XX));

listDX = listXfinal - listXinit;
listDY = listYfinal - listYinit;

% Plot initial and final positions
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

% Find points with significant movement (threshold of one pixel moved).
listMoved = listDisps > 2;
listGood = listMoved & (listYinit > 225) & (listYinit < 575);
numberMoved = sum(listGood);

% 
xI = listXinit(listGood);
yI = listYinit(listGood);
xF = listXfinal(listGood);
yF = listYfinal(listGood);


% Plot initial and final positions
figure(31)
imagesc(imROI)
colormap(gray)
hold on
  scatter(xI, yI, 'r')
  scatter(xF, yF, 'c')
hold off
axis equal
legend('frame 1', ['frame', int2str(nSteps)])
title('Interpolated positions');
xlabel('X-position, pixels')
ylabel('Y-position, pixels')
set(gca, 'fontSize', 14)

a = yI - yF;
b = xF - xI;
c = yI.*(xF - xI) + (yI - yF).*xI;

A = [a,b];

apex_position = A\c


x_result    = apex_position(1);
y_predicted = yI + (-a./b).*(x_result - xI);


%% Now work out which lines are far from the median result
median_y = median(y_predicted);
std_y    = std(y_predicted); % Try inter-quartile range (some outliers)

listVeryGood = zeros(size(xI));
listVeryGood( abs(y_predicted - median_y) < 1*std_y ) = 1;

xI_vg = xI(listVeryGood==1);
yI_vg = yI(listVeryGood==1);
xF_vg = xF(listVeryGood==1);
yF_vg = yF(listVeryGood==1);


% Plot initial and final positions
figure(31)
imagesc(imROI)
colormap(gray)
hold on
  scatter(xI_vg, yI_vg, 'r')
  scatter(xF_vg, yF_vg, 'c')
hold off
axis equal
legend('frame 1', ['frame', int2str(nSteps)])
title('Interpolated positions');
xlabel('X-position, pixels')
ylabel('Y-position, pixels')
set(gca, 'fontSize', 14)

a = yI_vg - yF_vg;
b = xF_vg - xI_vg;
c = yI_vg.*(xF_vg - xI_vg) + (yI_vg - yF_vg).*xI_vg;

A = [a,b];

apex_position = A\c

% 

figure(31)
hold on
  scatter(apex_position(1), apex_position(2), 'r+')
	scatter(xI(listVeryGood==0), yI(listVeryGood==0),'gx')
hold off
ylim([0 apex_position(2)*1.1])

% Plot initial and final positions
figure(12)
hold on
  scatter(apex_position(1), apex_position(2), 'r+')
	scatter(xI(listVeryGood==0), yI(listVeryGood==0),200,'gx')
hold off
legend('frame 1', ['frame', int2str(nSteps)], 'apex position')
ylim([0 apex_position(2)*1.1])


% Rough work for uplift rig image processing. 
%
% INSTRUCTIONS
%   Set working directory to contain the sample data, 'test_stack.tif'
%   Run this script.
%   Use the cpselect() tool to select matching points
%   The second call of cpselect can add extra pairs
%   Generate overlays with tracked particle displacements
%   Infer displacement field, and interpolate a grid of displacements

% Input

imRaw1 = imread('test_stack.tif',1); % Both RGB uint8
imRaw2 = imread('test_stack.tif',2);

% In cpselect(in,base), "in" is to be warped to "base"
[input_points, base_points] = cpselect(imRaw1, imRaw2, ...
                                       'Wait',true);
% % Add new points:
[input_points, base_points] = cpselect(imRaw1, imRaw2, ...
                              input_points, base_points, 'Wait',true);
                                   
figure(2)
imagesc(imRaw2)
hold on
scatter(base_points(:,1), base_points(:,2), 'c')
scatter(input_points(:,1), input_points(:,2), 'r')
legend('peak','initial')

figure(3)
imagesc(imRaw2)
hold on
for lpPt = 1:size(base_points,1)
   plot([base_points(lpPt,1);input_points(lpPt,1)],...
        [base_points(lpPt,2);input_points(lpPt,2)], 'b', 'lineWidth', 2) 
end

mytform = cp2tform(base_points, input_points, 'piecewise linear');
% 'piecewise linear' is stricter than local weighted mean 'lwm'


[XX,YY] = meshgrid([500:25:1200], [450:25:850]);

listXinit = XX(:);
listYinit = YY(:);

[listXfinal, listYfinal] = tforminv(mytform, listXinit, listYinit);

figure(5)
scatter(listXinit, -listYinit, 'r')
hold on
  scatter(listXfinal, -listYfinal, 'c')
hold off

figure(6)
imagesc(imRaw2)
hold on
  scatter(listXinit, listYinit, 'r')
  scatter(listXfinal, listYfinal, 'c')
hold off


figure(7)
imagesc(imRaw2)
hold on
for lpPt = 1:size(listXinit,1)
   plot([listXinit(lpPt);listXfinal(lpPt)],...
        [listYinit(lpPt);listYfinal(lpPt)], 'b', 'lineWidth', 1) 
end
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 1: CORRELATION OF OPTICAL IMAGES OF SPECKLED SAMPLES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ex1_OpticImages_SoftwoodGrowthRings.m

%% Moisture-induced swelling gradients in growth rings at hygroscopic equilibrium
%% are determined from optical images of speckled softwood samples
%% Reference image: RH = 0%, Test image: RH = 95%

% Reference: [Sanabria S J, Lanvermann C, Michel F, Mannes D, Niemz P (2014) Adaptive neutron radiography correlation for simultaneous 
% imaging of moisture transport and deformation in hygroscopic materials]

addpath('../src');

visualizeResults = 1; % 1: only visualize stored results, 0: compute again

if (visualizeResults ~= 1)

% Load test images
input_ref  = double(imread('../data/Ex1_ImageRef.tif'));  input_ref = input_ref/max(abs(input_ref(:)));
input_test = double(imread('../data/Ex1_ImageTest.tif')); input_test = input_test/max(abs(input_test(:)));

% Define ANRC deformation computation parameters
parameters_deformation = struct('WR', [71 71], ...
                    'WS', [3 3], ...
                    'WL', [5 5], ...
                    'subset', [9 9],... 
                    'ov',     [1 1], ... 
                    'u_max',  [71 71], ...
                    'u_max_increment', [3 3], ...
                    'threshold_R', 0.7, ...
                    'NEXT', [10 10], ...
                    'NERR', 5, ...
                    'monitoring_loop', [0 1], ...
                    'tracking_on', 1, ...
                    'raster_scan_modus', 1, ...
                    'decimate_output', 1, ...
                    'precorrect_image', 0, ...
                    'independent_error', 1, ...
                    'avoid_edges', [1 1], ...
                    'S', 10, ...
                    'RTpauseTime', 0.1);
% Define ANRC strain computation parameters
parameters_strain = struct('latresx_ux', 9/3,  ...         
                           'latresy_ux', 9/3, ...
                           'latresx_uy', 9/3, ...
                           'latresy_uy', 9/3, ...
                           'fixnans', 0);
                
% Compute deformation fields with ANRC
output_def = ADIC_deformation(input_test, input_ref, parameters_deformation);

% Compute strain fields with ANRC
output_str = ADIC_strain(output_def, parameters_strain);


% Calculate average profiles
epsXX = squeeze(output_str(:,:, 6))'; epsYY = squeeze(output_str(:,:, 7))'; epsXY = squeeze(output_str(:,:, 8))';
whereNaNs_epsXX = isnan(epsXX);      whereNaNs_epsYY = isnan(epsYY);      whereNaNs_epsXY = isnan(epsXY);
epsXX(whereNaNs_epsXX == 1) = 0;     epsYY(whereNaNs_epsYY == 1) = 0;     epsXY(whereNaNs_epsXY == 1) = 0;
sumNaNs_epsXX = sum(whereNaNs_epsXX == 0, 1); sumNaNs_epsYY = sum(whereNaNs_epsYY == 0, 1); sumNaNs_epsXY = sum(whereNaNs_epsXY == 0, 1);
epsXX_mean = mean(epsXX, 1).*size(epsXX, 1)./sumNaNs_epsXX;
epsYY_mean = mean(epsYY, 1).*size(epsYY, 1)./sumNaNs_epsYY;
epsXY_mean = mean(epsXY, 1).*size(epsXY, 1)./sumNaNs_epsXY;



% Store results
output_def = single(output_def); output_str = single(output_str);
epsXX_mean = single(epsXX_mean); epsYY_mean = single(epsYY_mean); 
epsXY_mean = single(epsXY_mean);

% In Ex1, we decimate the results by a factor of 3 to reduce storage space - This line can
% be removed
output_def = output_def(1:3:end,1:3:end, 1:end);
output_str = output_str(1:3:end,1:3:end, 1:end);
epsXX_mean = epsXX_mean(1:3:end);epsYY_mean = epsYY_mean(1:3:end);epsXY_mean = epsXY_mean(1:3:end);

save('../results/Ex1_results.mat', 'output_def', 'output_str', 'epsXX_mean', 'epsYY_mean', 'epsXY_mean');

end

% Load results
load('../results/Ex1_results.mat');

% Display deformation fields
fh1 = figure( 'Name', 'Deformation fields', 'NumberTitle', 'off', 'Visible', 'on', 'Color', [1,1,1]);
subplot(3,3,1); imshow(double(output_def(:,:, 1)), []); colormap bone;colorbar; daspect('auto'); title('Reference image');
subplot(3,3,2); imshow(double(output_def(:,:, 2)), []); colormap bone;colorbar; daspect('auto'); title('Test image');
subplot(3,3,3); imshow(double(output_def(:,:, 3)), []); colormap bone;colorbar; daspect('auto'); title('Corrected image'); 
subplot(3,3,4); imshow(double(output_def(:,:, 4)), []); colormap jet; colorbar;daspect('auto'); title('Deformation u_X (px)');
subplot(3,3,5); imshow(double(output_def(:,:, 5)), []); colormap jet; colorbar;daspect('auto'); title('Deformation u_Y (px)');
subplot(3,3,7); imshow(double(output_def(:,:, 6)), []); colormap jet; colorbar;daspect('auto'); title('Error u_X');
subplot(3,3,8); imshow(double(output_def(:,:, 7)), []); colormap jet; colorbar;daspect('auto'); title('Error u_Y');
subplot(3,3,9); imshow(double(output_def(:,:, 8)), []); colormap jet; colorbar; daspect('auto'); title('Error code');

% Display strain fields
fh2 = figure( 'Name', 'Strain fields', 'NumberTitle', 'off', 'Visible', 'on', 'Color', [1,1,1]);
subplot(2,3,1); imshow(double(output_str(:,:, 4)), []); colormap jet; colorbar; daspect('auto'); title('Spline fit deformation u_X (px)');
subplot(2,3,2); imshow(double(output_str(:,:, 5)), []); colormap jet; colorbar; daspect('auto'); title('Spline fit deformation u_Y (px)');
subplot(2,3,4); imshow(double(output_str(:,:, 6))*100, []); colormap jet; colorbar; daspect('auto'); title('Strain {\epsilon}_{XX} (%)');
subplot(2,3,5); imshow(double(output_str(:,:, 7))*100, []); colormap jet; colorbar; daspect('auto'); title('Strain {\epsilon}_{YY} (%)');
subplot(2,3,6); imshow(double(output_str(:,:, 8))*100, []); colormap jet; colorbar; daspect('auto'); title('Strain {\epsilon}_{XY} (%)');

% Display average profiles
fh3 = figure( 'Name', 'Mean strain profiles in X', 'NumberTitle', 'off', 'Visible', 'on', 'Color', [1,1,1]);
subplot(3,1,1); plot(double(epsXX_mean)*100, 'k'); ylabel('Strain {\epsilon}_{XX} (%)');
subplot(3,1,2); plot(double(epsYY_mean)*100, 'k'); ylabel('Strain {\epsilon}_{YY} (%)');
subplot(3,1,3); plot(double(epsXY_mean)*100, 'k'); ylabel('Strain {\epsilon}_{XY} (%)'); xlabel('Y position (px)');



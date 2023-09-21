function   output = ADIC_strain(input, parameters)
% function output = ADIC_strain(input, parameters)
% ANRC : Adaptive Neutron Radiography Correlation. 
%
% Computes strain fields by fitting discontinuous deformation fields to
% continuous cubic smoothing splines and differentiating along the spatial coordinates
% Strains are defined as: eps_xx = dux/dx; eps_yy = duy/dy; eps_xy = 0.5*(dux/dy + duy/dx)
%
% Copyright: Sergio Sanabria, ETH Zürich, 03.09.2013 (ssanabria@ethz.ch)
%
% Reference for citation: [Sanabria S J, Lanvermann C, Michel F, Mannes D, Niemz P (2013) Adaptive neutron radiography correlation for simultaneous imaging
% of hygroscopic moisture transport and swelling in biological composites]
%
% input (N1 X N2 X 8): deformation fields calculated with ANRC140206_deformation
% parameters: struct containing
%   latresx_ux: Lateral resolution in X of deformation field u_X
%   latresy_ux: Lateral resolution in Y of deformation field u_X
%   latresx_uy: Lateral resolution in X of deformation field u_Y
%   latresy_uy: Lateral resolution in Y of deformation field u_Y
%   fixnans:    If 1, fix NaN values in the deformation fields to
%               interpolated spline values, otherwise set these values to NaN
% output is a (N1 x N2 x X) matrix containing:
%   output(1:end, 1:end, 1) = test image compensated for deformation with respect to corrected def. fields
%   output(1:end, 1:end, 2) = corrected u_x (if fixnans = 1,  NaN values set to u_x_spline)
%   output(1:end, 1:end, 3) = corrected u_y (if fixnans = 1,  NaN values set to u_y_spline)
%   output(1:end, 1:end, 4) = spline fit of u_x (u_x_spline)
%   output(1:end, 1:end, 5) = spline fit of u_y (u_y_spline)
%   output(1:end, 1:end, 6) = strain field eps_xx
%   output(1:end, 1:end, 7) = strain field eps_yy
%   output(1:end, 1:end, 8) = strain field eps_xy

fix_Nans =  parameters.fixnans;
latres1_u = parameters.latresx_ux; latres2_u = parameters.latresy_ux;
latres1_v = parameters.latresx_uy; latres2_v = parameters.latresy_uy;
u = input(1:end,1:end, 4); error_u =  1 - input(1:end, 1:end, 6);
v = input(1:end,1:end, 5); error_v =  1 - input(1:end, 1:end, 7);
yd = input(1:end, 1:end, 2);
ox1 = 1:size(input, 1); x1 = ox1;
ox2 = 1:size(input, 2); x2 = ox2;
output = zeros(size(input, 1), size(input, 2), 8);

% Resolution parameters for smoothing function
p1u = 1/(1 + latres1_u^3/6); p2u = 1/(1 + latres2_u^3/6); p1v = 1/(1 + latres1_v^3/6); p2v = 1/(1 + latres2_v^3/6); % Opt. p-wert w.r.t. lateral resolution

% Weights of error function         
whereNaNs_u = isnan(error_u); whereNaNs_v = isnan(error_v); error_u(whereNaNs_u == 1) = 0; error_v(whereNaNs_v == 1) = 0;
sumNaNs_u_x1 = sum(whereNaNs_u == 0, 2); sumNaNs_u_x2 = sum(whereNaNs_u == 0, 1); sumNaNs_v_x1 = sum(whereNaNs_v == 0, 2); sumNaNs_v_x2 = sum(whereNaNs_v == 0, 1);
error_u_x1 = mean(error_u, 2).*size(error_u, 2)./sumNaNs_u_x1; error_u_x2 = mean(error_u, 1).*size(error_u, 1)./sumNaNs_u_x2; 
error_v_x1 = mean(error_v, 2).*size(error_v, 2)./sumNaNs_v_x1; error_v_x2 = mean(error_v, 1).*size(error_v, 1)./sumNaNs_v_x2;
error_u_x1(isnan(error_u_x1) == 1) = 0; error_u_x2(isnan(error_u_x2) == 1) = 0; error_v_x1(isnan(error_v_x1) == 1) = 0;  error_v_x2(isnan(error_v_x2) == 1) = 0;    

% Nan cleaning for spline interpolation 
u_clean = u; v_clean = v; ox2_nonan = []; ox1_nonan = ox1;
       
% Get rid of lines with only NaN values - at least two data sites
for index2 = 1:length(ox2), if((numel(find(isnan(u(1:end, index2)) == 0)) >= 2) && (numel(find(isnan(v(1:end, index2)) == 0)) >=2)); ox2_nonan  = [ox2_nonan, index2]; end; end
for index2 = ox2_nonan, % Substitute NaN by spline values
    u1 = u(ox1_nonan, index2); cs_u1 = csaps(x1(ox1_nonan), u1, p1u, [], error_u_x1(ox1_nonan)); u_spline1 = fnval(cs_u1, x1(ox1_nonan));
    u1(isnan(u1) == 1) = u_spline1(isnan(u1) == 1); u_clean(ox1_nonan, index2) = u1; 
    v1 = v(ox1_nonan, index2); cs_v1 = csaps(x1(ox1_nonan), v1, p1v, [], error_v_x1(ox1_nonan)); v_spline1 = fnval(cs_v1, x1(ox1_nonan));
    v1(isnan(v1) == 1) = v_spline1(isnan(v1) == 1); v_clean(ox1_nonan, index2) = v1; end

% Cubic smoothing splines
cs_u = csaps({x1(ox1_nonan), x2(ox2_nonan)},u_clean(ox1_nonan, ox2_nonan), {p1u,p2u}, [], {error_u_x1(ox1_nonan), error_u_x2(ox2_nonan)}); cs_u = fnxtr(cs_u);% cubic smoothing spline
cs_v = csaps({x1(ox1_nonan), x2(ox2_nonan)},v_clean(ox1_nonan, ox2_nonan), {p1v, p2v}, [], {error_v_x1(ox1_nonan), error_v_x2(ox2_nonan)}); cs_v = fnxtr(cs_v);% cubic smoothing spline
u_spline = fnval(cs_u,{x1, x2}); v_spline = fnval(cs_v,{x1, x2});
% p = 0, least-squares straight line fit, variational or 'natural' cubic spline interpolate
    
% First derivate (strain), fit to splines
csder_ux = fnder(cs_u, [1 , 0]); csder_uy = fnder(cs_u, [0 , 1]); csder_vx = fnder(cs_v, [1,  0]); csder_vy = fnder(cs_v, [0, 1]);
csder_ux_cs= csaps({x1, x2}, fnval(csder_ux,{x1, x2}), {p1u,p2u}, [], {error_u_x1, error_u_x2});csder_ux_cs = fnxtr(csder_ux_cs);
csder_uy_cs= csaps({x1, x2}, fnval(csder_uy,{x1, x2}), {p1u,p2u}, [], {error_u_x1, error_u_x2});csder_uy_cs = fnxtr(csder_uy_cs);
csder_vx_cs= csaps({x1, x2}, fnval(csder_vx,{x1, x2}), {p1v,p2v}, [], {error_v_x1, error_v_x2});csder_vx_cs = fnxtr(csder_vx_cs);
csder_vy_cs= csaps({x1, x2}, fnval(csder_vy,{x1, x2}), {p1v,p2v}, [], {error_v_x1, error_v_x2});csder_vy_cs = fnxtr(csder_vy_cs);
exx_spline = fnval(csder_ux_cs, {x1, x2}); eyy_spline = fnval(csder_vy_cs, {x1, x2});
exy_spline = 0.5*(fnval(csder_uy_cs, {x1, x2}) + fnval(csder_vx_cs, {x1, x2}));
 
% Fix NaNs?
u_clean = u; v_clean = v; u_clean(isnan(u) == 1) = u_spline(isnan(u) == 1); v_clean(isnan(v) == 1) = v_spline(isnan(v) == 1);  
if (fix_Nans == 1); u_out = u_clean; v_out = v_clean;
else u_out = u; v_out = v; u_spline(isnan(u) == 1) = NaN; v_spline(isnan(v) == 1) = NaN; exx_spline(isnan(u) == 1) = NaN; eyy_spline(isnan(v) == 1) = NaN; exy_spline((isnan(u) == 1) | (isnan(v) == 1)) = NaN;  end
% Interpolate deformed image with NaN-corrected deformation fields
NaN_cancel_active = 1; % Cancel NaNs obtained at image interpolation with discontinuous deformation fields(may slightly corrupt statistics at edge values)   
[ox2g, ox1g] = meshgrid(ox2,ox1); yd_hR_corr = ANRCAlgorithmNaNFreeInterp2(yd, ox2g + v_out, ox1g + u_out, 'spline', NaN_cancel_active);   

% Output
output(1:end,1:end, 1) = yd_hR_corr; output(1:end,1:end, 2) = u_out; output(1:end,1:end, 3) = v_out;
output(1:end,1:end, 4) = u_spline; output(1:end,1:end, 5) = v_spline;
output(1:end,1:end, 6) = exx_spline; output(1:end,1:end, 7) = eyy_spline; output(1:end,1:end, 8) = exy_spline;    
end
   
function [inputi] = ANRCAlgorithmNaNFreeInterp2(input, x2_values, x1_values, interpol_mode, activate)
% function [inputi] = ANRCAlgorithmNaNFreeInterp2(input, x2_values, x1_values, interpol_mode, activate)
% activate = 1: Performs interpolation of successively low order to cancel NaN values. If NaNs in higher order interp disappear in lower interp step, then take
% values of lower interp step. 
% Can lead to some loss of accuracy due to lower order interpolation
% One can deactivate lower order interpolation, set flag activate = 0
% if interpol_mode = 'spline' -> do spline, then linear, then nearest.
warning('off', 'all');
if (activate == 1)
    if (strcmp(interpol_mode, 'spline') == 1)
        inputi = interp2(input, x2_values, x1_values, 'spline'); inputi_nonan = interp2(input, x2_values, x1_values, 'linear');
        inputi((isnan(inputi) == 1) & (isnan(inputi_nonan) == 0)) = inputi_nonan((isnan(inputi) == 1) & (isnan(inputi_nonan) == 0)); inputi_nonan = interp2(input, x2_values, x1_values, 'nearest');
        inputi((isnan(inputi) == 1) & (isnan(inputi_nonan) == 0)) = inputi_nonan((isnan(inputi) == 1) & (isnan(inputi_nonan) == 0));
    elseif (strcmp(interpol_mode, 'linear') == 1)
        inputi = interp2(input, x2_values, x1_values, 'linear'); inputi_nonan = interp2(input, x2_values, x1_values, 'nearest');
        inputi((isnan(inputi) == 1) & (isnan(inputi_nonan) == 0)) = inputi_nonan((isnan(inputi) == 1) & (isnan(inputi_nonan) == 0));
    else inputi = interp2(input, x2_values, x1_values, 'nearest'); end
else inputi = interp2(input, x2_values, x1_values, interpol_mode); end
warning('on', 'all');
end



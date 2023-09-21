function   output = ADIC_deformation(input1, input2, parameters)
% function output = ADIC_deformation(input1, input2,parameters)
% ANRC : Adaptive Neutron Radiography Correlation. 
% 
% Computes deformation fields out of image natural texture by performing digital image correlation 
% between test and reference images. The ANRC algorithm is a texture correlation (TC) algorithm, which is optimized and expanded to robustly extract
% deformation information from neutron radiographies while ignoring "shadow regions", where only poor correlation statistics are available.
% The inner TC core performs a zero-order search (rigid subset) with a zero-normalized cross-correlation calculated at integer positions and
% refined with bicubic interpolation. The center position and size of the search window are adaptively adjusted. Error control ensure correlation,
% uniqueness and continuity of deformation estimates.
%
% Copyright: Sergio Sanabria, ETH Zürich, 03.09.2013 (ssanabria@ethz.ch)
%
% Reference for citation: [Sanabria S J, Lanvermann C, Michel F, Mannes D, Niemz P (2013) Adaptive neutron radiography correlation for simultaneous imaging
% of hygroscopic moisture transport and swelling in biological composites]
%
% input1 (N1 x N2): test image 
% input2 (N1 x N2): reference image
% parameters: struct containing
%   WR (2 x 1): Reset search window size (dimension 1: dim1, dimension 2: dim2)
%   WS (2 x 1): Short search window size
%   WL (2 x 1): Long search window size
%   subset (2 x1): Subset size
%   ov (2 x 1):    Step size
%   u_max (2 x 1): Maximum allowed deformation (typically = WR)
%   u_max_increment (2 x 1): Maximum allowed deformation increment (typically = WS)
%   avoid_edges (2 x 1): If 1, do not include image edges in correlation
%   threshold_R:  Minimum accepted correlation coefficient
%   NEXT (1 x 1): No. values used   for extrapolation (20)
%   NERR (1 x 1): No. errors before reset (typically NERR <= NEXT) (10)
%   monitoring_loop (2 x 1): If 1, show monitoring info in dim1, dim2
%   tracking_on (1 x1): If 1, use adaptive operation, otherwise don't use
%   raster_scan_modus (1x 1): If 1, explore image in raster fashion, otherwise treat lines independently
%   decimate_output (1 x 1):  If 1, provide output decimated with ov, if 0  provide output in same grid as input1, input2
%   precorrect_image (1 x 1): If 1, previous to ANRC eliminate reproducible patterns (assuming symmetry along dim2) 
%   independent_error (1 x 1): If 1, treat errors in dim1 and dim2 separately
%   calculation. If 0, include all available pixels within search range window, excluding pixels out of edges and NaN values
%   S: interpolation factor for correlation function (typically 100)
%   RTpauseTime: pause time after each monitoring snapshot in s (-1: wait for user to press keyboard)
%
%   with decimate_output = 0, output is a (N1 x N2 x 8) matrix containing:
%       output(1:end, 1:end, 1) =  input1 % test image 
%       output(1:end, 1:end, 2) =  input2 % reference image (reference
%       output(1:end, 1:end, 3) =  input2 compensated for deformation w.r.t. to reference
%       output(1:end, 1:end, 4) =  u1: deformation field in dim1
%       output(1:end, 1:end, 5) =  u2: deformation field in dim2
%       output(1:end, 1:end, 6) =  error1: error in dim1 (=1: error,  <1: 1 - regression_coefficient)
%       output(1:end, 1:end, 7) =  error2: error in dim2 (=1: error,  <1: 1 - regression_coefficient)
%       output(1:end, 1:end, 8) =  error_code (for error1 or error2 = 1, gives identified error cause):
%               EVERYTHING_OK = 0;                   % 0: no error
%               OUT_OF_LIMITS = 1;                   % 1: deformation vector out of image boundaries
%               UNDER_THRESHOLD_R = 2;               % 2: correlation below threshold R
%               OUTSIDE_CORR_SEARCH = 3;             % 3: correlation maximum at edge of search window 
%               DEF_TOO_LARGE = 4;                   % 4: u or v go over max. deformation limits 
%               INCR_TOO_LARGE = 5;                  % 5: too large deformation increment
%               LARGEWINDOW_UNCONSISTENCY = 6;       % 6: deformation fields change by more than one pixel when enlarging search window
%
%   with decimate_output = 1, output is decimated with respect to input2(i1, i2) as follows: 
%       i1d = [1 + floor((ov(1) - 1)/2)]: ov(1): [end - floor(ov(1)/2];
%       i2d = [1 + floor((ov(2) - 1)/2)]: ov(2): [end - floor(ov(2)/2];
% For correct operation, the following points need to be observed:
%   1) input1, input2: Same size. 
%   2) subset: Odd values
%   3) at image edges, subset values outside image are not included in correlation computation, thus the correlation at image edges rely on a smaller number of points and may be more noisy. 
%      this uncertainty region corresponds for zero deformation to a  distance subset/2 from the edges. The noisy region can be avoided in the calculation by
%      setting avoid_edges = 1 in the desired dimensions
%   4) for swelling deformation, the size of the usable area in the deformed image may be larger than in the reference image. It is possible to avoid performing correlation in non-usable pixels of the reference 
%      image by setting them to NaN

    NDIMENSIONS = 2; S = parameters.S; WR = parameters.WR; WS = parameters.WS; WL = parameters.WL; u_max = parameters.u_max; u_max_increment = parameters.u_max_increment; 
    subset = parameters.subset; ov = parameters.ov; NEXT = parameters.NEXT;  monitoring_loop = parameters.monitoring_loop; avoid_edges =parameters.avoid_edges;
    WR = reshape(WR, [length(WR), 1]); WS = reshape(WS, [length(WS), 1]); WL = reshape(WL, [length(WL), 1]); u_max = reshape(u_max, [length(u_max), 1]); avoid_edges = reshape(avoid_edges, [length(avoid_edges), 1]); 
    u_max_increment = reshape(u_max_increment, [length(u_max_increment), 1]); subset = reshape(subset, [length(subset), 1]); ov = reshape(ov, [length(ov), 1]); NEXT = reshape(NEXT, [length(NEXT), 1]);
    monitoring_loop = reshape(monitoring_loop, [length(monitoring_loop), 1]); 
    
    threshold_R = parameters.threshold_R; tracking_on = parameters.tracking_on; independent_error = parameters.independent_error; raster_scan_modus = parameters.raster_scan_modus; 
    decimate_output = parameters.decimate_output; precorrect_image = parameters.precorrect_image; NERR = parameters.NERR;
    RTpauseTime = parameters.RTpauseTime;
    
    EVERYTHING_OK = 0;                   % 0: no error
    OUT_OF_LIMITS = 1;                   % 1: window out of sample limits
    UNDER_THRESHOLD_R = 2;               % 2: correlation below threshold R
    OUTSIDE_CORR_SEARCH = 3;             % 3: correlation maximum at edge of search window
    DEF_TOO_LARGE = 4;                   % 4: u or v go over max. deformation limits
    INCR_TOO_LARGE = 5;                  % 5: too large deformation increment
    LARGEWINDOW_UNCONSISTENCY = 6;       % 6: deformation fields change by more than one pixel when enlarging search window
    
    % Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ydcorr  = input1; % Deformed dataset
    yucorr  = input2; % Non-deformed dataset to be corrected   

    
    dimsize = zeros(1, NDIMENSIONS); 
    for index=1:NDIMENSIONS, dimsize(index) = size(input1, index) + 2*(WR(index) + (subset(index) - 1)/2 ); end;
    yu = NaN*ones(dimsize); yd = yu;               
    rn = struct('x1', []); % X1 coordinates
    rn_start = rn; rn_end = rn; rn_starteff = rn; rn_endeff = rn; rn_crop = rn; 
    for index = 1:NDIMENSIONS, 
        rn.(['x', num2str(index)]) = transp(1:size(yu, index)); 
        rn_start.(['x', num2str(index)]) =  (WR(index) + (subset(index) - 1)/2 ) + 1; 
        rn_end.(['x', num2str(index)]) =  length(rn.(['x', num2str(index)])) - ((subset(index) - 1)/2 + WR(index));   % Enough pixels for healthy correlation   
        if (avoid_edges(index) == 1); rn_starteff.(['x', num2str(index)]) = rn_start.(['x', num2str(index)]) + floor((ov(index) - 1)/2) + (WR(index) + (subset(index) - 1)/2 ) ;
        else rn_starteff.(['x', num2str(index)]) = rn_start.(['x', num2str(index)]) + floor((ov(index) - 1)/2); end
        if (avoid_edges(index) == 1); rn_endeff.(['x', num2str(index)])   = rn_end.(['x', num2str(index)]) - ((subset(index) - 1)/2 + WR(index));
        else rn_endeff.(['x', num2str(index)])   = rn_end.(['x', num2str(index)]); end
        rn_crop.(['x', num2str(index)]) = transp(rn_starteff.(['x', num2str(index)]):ov(index):rn_endeff.(['x', num2str(index)])); % X1, X2 coordinates of calculated deformation sequence       
        if (numel(rn_crop.(['x', num2str(index)])) == 0); disps(['Too few pixels in dimension ', num2str(index)]); output = []; return; end
    end % X2, ... XN    
    
  
    
    dimsize= zeros(1, NDIMENSIONS); for index =1:NDIMENSIONS, dimsize(index) = length( rn_crop.(['x', num2str(index)])); end
    u = NaN*ones([dimsize, NDIMENSIONS]);          % Deformation field
    u_track = zeros([dimsize, NDIMENSIONS]);       % Tracking vector (actual iteration)
    W = zeros([dimsize, NDIMENSIONS]);             % Search window
    u_cont  = zeros([dimsize, NDIMENSIONS]);       % Continuity vector (actual iteration)
    error_u = zeros([dimsize, NDIMENSIONS]);     % Error measure (1 - r)
    error_code = zeros([dimsize, NDIMENSIONS]);  % Error code
    count_error = zeros(dimsize);                % Error count
    count_reset = zeros(dimsize);                % Reset count % 0 - no reset, 1 - first consecutive reset, 2 - second consecutive reset, etc.
    count_nonused = zeros(dimsize);              % 1 if pixel is not used for correlation, 0 otherwise
    count_indexloop = zeros([dimsize, NDIMENSIONS]); % Records looping index used at each image pixel, useful for debugging 
    if (decimate_output ~= 1); output = zeros(size(input1, 1), size(input1, 2), 8);
    else                       output = zeros(size(error_u, 1), size(error_u, 2), 8); end
    
    fhmonitor_loop = struct('x1', []); % Figure handles for monitor in each dimension
    
    for index = 1:NDIMENSIONS,  if (monitoring_loop(index) == 1),
            fhmonitor_loop.(['x', num2str(index)]) = figure( 'Name', 'Online monitor', 'NumberTitle', 'off', ...
                'Visible', 'on', 'Color', [1,1,1]); end; end;
        
    %%%%%%%%%%% Image precorrection, eliminate reproducible patterns
    if (precorrect_image == 1);        [yucorr, ydcorr] = ANRCAlgorithmPrecorr(yucorr, ydcorr, NDIMENSIONS);    end   
    
    yu(rn_start.x1:rn_end.x1, rn_start.x2:rn_end.x2) = yucorr(:,:); % Make images larger to avoid cropping of results
    yd(rn_start.x1:rn_end.x1, rn_start.x2:rn_end.x2) = ydcorr(:,:);
        
    %%%%%%%%%%% Raster loop definition
    numeloop = 1; for index = 1:NDIMENSIONS, numeloop = numeloop*length(rn_crop.(['x', num2str(index)])); end
    if ((tracking_on == 1) && (raster_scan_modus == 1))
        index1_loop_l = 1:numeloop; rn_loop = zeros(numeloop, NDIMENSIONS); off_index = 1;
        for index2 = 1:length(rn_crop.x2),             
            if (mod(index2, 2) == 1), rn_loop(off_index:(off_index + length(rn_crop.x1) - 1), 1) = 1:length(rn_crop.x1); 
            else rn_loop(off_index:(off_index + length(rn_crop.x1) - 1), 1) = length(rn_crop.x1):-1:1; end
            rn_loop(off_index:(off_index + length(rn_crop.x1) - 1), 2) = index2;
            off_index = off_index + length(rn_crop.x1);
        end
        index2_loop_l = 1;
    else
        index1_loop_l = 1:length(rn_crop.x1); rn_loop = zeros(length(rn_crop.x1), NDIMENSIONS); 
        rn_loop(1:length(rn_crop.x1), 1) = 1:length(rn_crop.x1);
        rn_loop(1:length(rn_crop.x2), 2) = 1;
        index2_loop_l = 1:length(rn_crop.x2);
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ANRC LOOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for index2_loop = index2_loop_l,        
        for index1_loop = index1_loop_l,  
            index1 = rn_loop(index1_loop, 1); % Current index in 2D coordinates
            index2 = max(rn_loop(index1_loop, 2), index2_loop); % Allows managing raster_scan_modus = 1 and raster_scan_modus = 0
            count_indexloop(index1, index2, 1) = index1_loop; count_indexloop(index1, index2, 2) = index2_loop;
            
            %if ((index2_loop == 1) && (index1_loop == 1333)) % For debugging
            %      pause(0.1); end
            
            if (index1_loop > 1);  
                index1p = rn_loop(index1_loop - 1, 1); index2p = max(rn_loop(index1_loop - 1, 2), index2_loop);
                if ((tracking_on == 1) && (raster_scan_modus == 1))                
                    if ((index2 > index2p) && (monitoring_loop(2) == 1)); 
                        ANRCAlgorithmLoop2Monitor(yd, rn_crop, index2p, u, fhmonitor_loop, yu, error_code, error_u, count_error, count_reset, u_track, W, count_indexloop, RTpauseTime); 
                    end
                else
                    if ((index2 > 1) && (index1 == 1)&& (monitoring_loop(2) == 1));
                        ANRCAlgorithmLoop2Monitor(yd, rn_crop, index2 - 1, u, fhmonitor_loop, yu, error_code, error_u, count_error, count_reset, u_track, W, count_indexloop, RTpauseTime); 
                    end
                end
            else first_loopindex = 1; end % First index, this flag is necessary to filter out NaN values in reference image, for which the looping does not count
            
            
            corrvalue = 0;                       
            % *********************************************************************************
            % RESET CONTROL
            % *********************************************************************************
            if ( (first_loopindex == 1) || (tracking_on == 0)) % Start iteration
                W(index1, index2, :) = WR(:); u_track(index1, index2, :) = 0; u_cont(index1, index2, :) = 0;
            else
                if (count_reset(index1p, index2p) >= 1) 
                    if (numel(find(error_u(index1p, index2p, :) == 1)) == 0); count_reset(index1, index2) = 0; else count_reset(index1, index2) = count_reset(index1p, index2p); end; end
                    % If we were in reset mode in previous it. and there is an error in any dimension keep same reset mode, otherwise cancel
                if (count_error(index1p, index2p) >= NERR) % We have exceeded the consecutive error count before reset      
                    % Starvation control - Get NERR + 2 pixels
                    [NSTARV_EFF, NSTARV_O] = ANRCAlgorithmGetNusedPixels(NERR + 2, index1_loop, index2_loop, count_nonused, rn_loop);
                    index1_starv = rn_loop(index1_loop - NSTARV_EFF, 1); % Index at starvation index
                    index2_starv = max(rn_loop(index1_loop - NSTARV_EFF, 2), index2_loop); 
                    if (NSTARV_O == NERR + 2) % Check for starvation situation detected (reset index = 1, NERR + 1 pixels ago, which was removed in this pixel given rise to error since, so that count_reset(index1p, index2p) = 0)
                         count_reset(index1, index2) = max(count_reset(index1_starv, index2_starv) + 1, count_reset(index1p, index2p) + 1); 
                    else count_reset(index1, index2) = count_reset(index1p, index2p) + 1; end
                    count_error(index1, index2) = 0; end % Allows detecting more than one consecutive reset               
                if (count_reset(index1, index2) >= 1); W(index1, index2, :) = WR(:); else W(index1, index2, :) = WS(:); end; end
                        
            if (monitoring_loop(1) == 1); disps(['Row ', num2str(rn_crop.x1(index1)), ' of ', num2str(length(rn.x1))]); end
           
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % DIC ALGORITHM
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Define center search positions and subset sizes
            ru = transp([rn_crop.x1(rn_loop(index1_loop, 1)), rn_crop.x2(max(rn_loop(index1_loop, 2), index2_loop))]); subsetu = (subset - 1)/2; 
            if (count_reset(index1, index2) > 0); rd = ru; else rd = ru(:) + squeeze(u_track(index1, index2, :)); end; subsetd = subsetu + squeeze(W(index1, index2, :));
            if (isnan(yu(ru(1), ru(2))) == 1); count_nonused(index1, index2) = 1; end; % Non-used pixel
            
            % If there is NaN in the undeformed image skip iteration 
            if (count_nonused(index1, index2) == 1)
                if (index1_loop > 1); error_u(index1, index2, :) = error_u(index1p, index2p, :); error_code(index1, index2, :) = error_code(index1p, index2p, :); count_error(index1, index2) = count_error(index1p, index2p); end
            else
                % Check for out of limits error
                out_of_limits_start = zeros(1, NDIMENSIONS); out_of_limits_end = zeros(1, NDIMENSIONS);
                for indexdim = 1:NDIMENSIONS,
                    if (avoid_edges(indexdim) == 1); out_of_limits_start(indexdim) = rd(indexdim) - subsetd(indexdim) - subsetu(indexdim) - WR(indexdim);
                                                     out_of_limits_end(indexdim)   = rd(indexdim) + subsetd(indexdim) + subsetu(indexdim) + WR(indexdim);
                    else  out_of_limits_start(indexdim) = rd(indexdim) - subsetd(indexdim);out_of_limits_end(indexdim)   = rd(indexdim) + subsetd(indexdim); end; end

                if ((ANRCAlgorithmOutOfLimits(out_of_limits_start, yd) == 1) || (ANRCAlgorithmOutOfLimits(out_of_limits_end, yd) == 1))
                    u(index1, index2, :) = NaN; error_u(index1, index2, :) = 1; error_code(index1, index2, :) = OUT_OF_LIMITS;
                else
                    % Here one can extend to higher order DIC algorithms
                    % DIC core in separate function % mode = 1; % ZERO-ORDER DIC without DATASET INTERPOLATION
                    mode = 1; [corrpos, corrvalue, p] = ANRCAlgorithmDIC2D(yu, yd, ru, rd, subsetu, subsetd, mode, monitoring_loop(1), fhmonitor_loop.x1, RTpauseTime);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % ERROR CONTROL
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if (corrvalue < threshold_R)
                        u(index1, index2, :) = NaN; error_u(index1, index2,:) = 1; error_code(index1, index2, :) = UNDER_THRESHOLD_R; % If below threshold_R we cancel search
                    else
                        if (independent_error == 1)
                            for index = 1:NDIMENSIONS,
                                if(p.atedge(index) == 1); error_u(index1, index2, index) = 1; u(index1, index2, index) = NaN; error_code(index1, index2, index) = OUTSIDE_CORR_SEARCH;
                                else if (abs(corrpos(index)) > u_max(index)); error_u(index1, index2, index) = 1; u(index1, index2, index) = NaN; error_code(index1, index2, index) = DEF_TOO_LARGE; end; end; end
                        else
                            if (sum(p.atedge) >= 1); u(index1, index2, :) = NaN; error_u(index1, index2,:) = 1; error_code(index1, index2, :) = OUTSIDE_CORR_SEARCH;
                            else for index = 1:NDIMENSIONS, if (abs(corrpos(index)) > u_max(index)); error_u(index1, index2, :) = 1; u(index1, index2, :) = NaN; error_code(index1, index2, :) = DEF_TOO_LARGE; break; end; end; end
                        end;
                    end
                    
                    % Long window routine
                    if ((tracking_on == 1) && (count_reset(index1, index2) == 0))
                    select_errorcode = 0; for index = 1:NDIMENSIONS, 
                    if ((error_code(index1, index2, index) == UNDER_THRESHOLD_R) ||(error_code(index1, index2, index) == OUTSIDE_CORR_SEARCH) ||(error_code(index1, index2, index) == DEF_TOO_LARGE)); select_errorcode = 1; break; end; end;
                    if (select_errorcode == 1) % We give it a second try with a larger tracking window
                        error_u_old = squeeze(error_u(index1, index2, :)); error_u(index1, index2, :) = 0; corrpos_old = corrpos; % Backup data
                        
                        W(index1, index2, :) = WL(:); % Make window larger and try again
                        subsetd = subsetu + squeeze(W(index1, index2, :));
                                                
                        % Check for out of limits error
                        out_of_limits_start = zeros(1, NDIMENSIONS); out_of_limits_end = zeros(1, NDIMENSIONS);
                        for indexdim = 1:NDIMENSIONS,
                            if (avoid_edges(indexdim) == 1); out_of_limits_start(indexdim) = rd(indexdim) - subsetd(indexdim) - subsetu(indexdim) - WR(indexdim);
                                out_of_limits_end(indexdim)   = rd(indexdim) + subsetd(indexdim) + subsetu(indexdim) + WR(indexdim);
                            else  out_of_limits_start(indexdim) = rd(indexdim) - subsetd(indexdim);out_of_limits_end(indexdim)   = rd(indexdim) + subsetd(indexdim); end; end
                        
                        if ((ANRCAlgorithmOutOfLimits(out_of_limits_start, yd) == 1) || (ANRCAlgorithmOutOfLimits(out_of_limits_end, yd) == 1))
                            u(index1, index2, :) = NaN; error_u(index1, index2, :) = NaN; error_code(index1, index2, :) = OUT_OF_LIMITS;
                        else
                            % DIC core in separate function % mode = 1; % ZERO-ORDER DIC without DATASET INTERPOLATION
                            mode = 1; [corrpos, corrvalue, p] = ANRCAlgorithmDIC2D(yu, yd, ru, rd, subsetu, subsetd, mode, monitoring_loop(1), fhmonitor_loop.x1, RTpauseTime);
                            if (corrvalue < threshold_R);
                                u(index1, index2, :) = NaN; error_u(index1, index2,:) = 1; error_code(index1, index2, :) = UNDER_THRESHOLD_R; % If below threshold_R we cancel search
                            else
                                if (independent_error == 1)
                                    for index = 1:NDIMENSIONS,
                                        if(p.atedge(index) == 1); error_u(index1, index2, index) = 1; u(index1, index2, index) = NaN; error_code(index1, index2, index) = OUTSIDE_CORR_SEARCH;
                                        else if (abs(corrpos(index)) > u_max(index)); error_u(index1, index2, index) = 1; u(index1, index2, index) = NaN; error_code(index1, index2, index) =DEF_TOO_LARGE; end; end; end
                                else
                                    if (sum(p.atedge) >= 1); u(index1, index2, :) = NaN; error_u(index1, index2,:) = 1; error_code(index1, index2, :) = OUTSIDE_CORR_SEARCH;
                                    else for index = 1:NDIMENSIONS, if (abs(corrpos(index)) > u_max(index)); error_u(index1, index2, :) = 1; u(index1, index2, :) = NaN; error_code(index1, index2, :) = DEF_TOO_LARGE; break; end; end; end
                                end
                            end
                        end 
                        % Check consistency with short window results
                        if ( (numel(find(error_u(index1, index2, :) == 1)) >= 1) && (numel(find(error_u(index1, index2, :) == 1)) < NDIMENSIONS) ); % Check than non-erroneous variables are consistent
                            for index = 1:NDIMENSIONS, if ((error_u_old(index) ~= 1) && (abs(corrpos_old(index) - corrpos(index)) > 1)); 
                                    u(index1, index2, index) = NaN; error_u(index1, index2, index) = 1; error_code(index1, index2, index) = LARGEWINDOW_UNCONSISTENCY; end; end; end
                        % If no errors, update error_code
                        for index = 1:NDIMENSIONS, if (error_u(index1, index2, index) ~= 1); error_code(index1, index2, index) = EVERYTHING_OK; end; end
                    end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % INTERPOLATION OF CORRELATION FUNCTION  % Correlation mode DIC orden 0
                if ((mode == 1) && (numel (find(error_u(index1, index2, :) == 1)) < NDIMENSIONS)) % If both with errors it does not make sense to peform this step
                           [corrpos, corrvalue, p] = ANRCAlgorithmDIC2DInt(S, p);
                end
                
                %% Continuity error check
                if (independent_error == 1)
                    for index = 1:NDIMENSIONS,  if ((error_u(index1, index2, index) ~= 1) && (first_loopindex == 0) && (count_reset(index1, index2) == 0) && (tracking_on == 1))
                                if (abs(corrpos(index) - u_cont(index1, index2, index)) > u_max_increment(index)) % No reset, no start, no error
                            error_code(index1, index2, index) = INCR_TOO_LARGE; u(index1, index2, index) = NaN;  error_u(index1, index2, index) = 1; end; end; end
                else
                    for index = 1:NDIMENSIONS,  if ((error_u(index1, index2, index) ~= 1) && (first_loopindex == 0) && (count_reset(index1, index2) == 0) && (tracking_on == 1))
                                if (abs(corrpos(index) - u_cont(index1, index2, index)) > u_max_increment(index)) % No reset, no start, no error
                            error_code(index1, index2, :) = INCR_TOO_LARGE; u(index1, index2, :) = NaN;  error_u(index1, index2, :) = 1; break; end; end; end; end
                
                %% Finish, assign variables
                % Error count                
                if (first_loopindex == 1); count_error_old = 0; else count_error_old = count_error(index1p, index2p);  if (count_error_old >= NERR); count_error_old = 0; end; end                
                if (numel (find(error_u(index1, index2, :) == 1)) > 0); count_error(index1, index2) = count_error_old + 1; else count_error(index1, index2) = max(count_error_old - 1, 0); end
                
                % Deformation fields
                for index = 1:NDIMENSIONS, if (error_u(index1, index2, index) ~= 1); u(index1, index2, index) = corrpos(index); error_u(index1, index2, index) = 1 - corrvalue; end; end;
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% FEEDBACK VARIABLE UPDATE                                                              %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if ((index1_loop ~= length(index1_loop_l)) && (tracking_on == 1)) % In last iteration, nothing to actualize
                index1n = rn_loop(index1_loop + 1, 1); % index for next iteration
                index2n = max(rn_loop(index1_loop + 1, 2), index2_loop); % index for next iteration
                
                for index = 1:NDIMENSIONS, % Independent dimension management                   
                    %%%%%%%%%%%%%%%%%%%%%%% u_track
                    % Create track lists, do to the possibility of non-used yu values (NaNs wo. errors) we need an iterative procedure
                    [NEXT_EFF, notused] = ANRCAlgorithmGetNusedPixels(NEXT(index), index1_loop, index2_loop, count_nonused, rn_loop);                                        
                    % Takes pixels including actual pixel until we have as much as we need NEXT(index), ignoring non-used pixels                    
                    index1_loop_trkl = transp((index1_loop - NEXT_EFF):index1_loop);
                    index1_trkl = rn_loop((index1_loop - NEXT_EFF):index1_loop, 1); % create track lists
                    index2_trkl = max(rn_loop((index1_loop - NEXT_EFF):index1_loop, 2), index2_loop);
                    lin_index_trkl_u = sub2ind(size(u), index1_trkl, index2_trkl, index*ones(size(index1_trkl))); lin_index_trkl_count_reset = sub2ind(size(count_reset), index1_trkl, index2_trkl);
                    u_trkl      = u(lin_index_trkl_u); error_u_trkl = error_u(lin_index_trkl_u); count_reset_trkl = count_reset(lin_index_trkl_count_reset); count_nonused_trkl = count_nonused(lin_index_trkl_count_reset);
                    % Select valid tracking values
                    validvalues_trkl = (error_u_trkl ~= 1) & (count_reset_trkl ~= 1) & (count_nonused_trkl ~= 1); % Used indexes, non-error values, with reset ~= 1(larger resets are considered to be good data)
                    index1_loop_trkl_valid = index1_loop_trkl(validvalues_trkl); u_trkl_valid = u_trkl(validvalues_trkl);
                    % Extrapolate next tracking value
                    numel_valid = numel(index1_loop_trkl_valid);
                    if (numel_valid >= 3); % Sufficient points for linear extrapolation
                        warning('off', 'all'); pol = polyfit(index1_loop_trkl_valid, u_trkl_valid, 1); warning('on', 'all');
                        u_track(index1n, index2n, index) = round((index1_loop + 1)*pol(1) + pol(2)); % Linear extrapolation (higher order splines could be fitted, avoid overfitting...)
                    else if (validvalues_trkl(end) == 0); % Actual value not valid
                            u_track(index1n, index2n, index) = u_track(index1, index2, index);                     % If current deformation fields not valid take last tracking value (similar to u_cont)
                        else if (numel_valid == 2); u_track(index1n, index2n, index) = round(mean(u_trkl_valid));  % No error: Take mean value (polynom of order 0)
                            elseif (numel_valid == 1);
                                if (first_loopindex == 1); u_track(index1n, index2n, index) = round(u_trkl_valid(1)); % Start or more than one reset (optionally introduce reset count > 2 here -> less smoothing:  "|| (count_reset(index1, index2) > 1")
                                else u_track(index1n, index2n, index) = round((u_track(index1, index2, index) + u_trkl_valid(1))*0.5); end % No reset or reset count > 2, smooth with current value to avoid discontinuous behaviour
                            end; end; end
                    
                    %%%%%%%%%%%%%%%%%%%%%%% u_cont
                    if (error_u(index1, index2, index) == 1);
                        u_cont(index1n, index2n, index) = u_cont(index1, index2, index);    % With error continue searching for same position as stored previously (we should never reach to a situation in which data is used once a reset comes)
                        % - Upon reset no more searchs are performed                        % - Reset is only cancelled, when both u and v have no errors
                    else if (count_reset(index1, index2) == 0);             u_cont(index1n, index2n, index) = u(index1, index2, index);       % No error and no reset take actual position (holds as well for first pixel)
                        else if (error_u(index1p, index2p, index) == 0);    u_cont(index1n, index2n, index) = u(index1, index2, index);       % No error, neither in previous position take actual position
                            else if (count_reset(index1, index2) == 1);     u_cont(index1n, index2n, index) = u_track(index1n, index2n, index); % Reset = 1 (first pixel after error), take tracking index for actual pixel
                                else                                        u_cont(index1n, index2n, index) = u_track(index1n, index2n, index); % Reset > 1 (first pixel after error), eliminate any past dependency
                                end; end; end; end
                end
                
                if(count_nonused(index1, index2) == 0); first_loopindex = 0; end % Not anymore the first loop index                
            end
       
            
        end
    end           
    if (monitoring_loop(2) == 1); ANRCAlgorithmLoop2Monitor(yd, rn_crop, index2, u, fhmonitor_loop, yu, error_code, error_u, count_error, count_reset, u_track, W, count_indexloop, RTpauseTime); end % For last iteration
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ENDING TASKS                                                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Clean deformation fields from isolated pixels
    for index = 1:NDIMENSIONS,        
        for index2_loop = index2_loop_l, 
            for index1_loop = index1_loop_l,
                index1 = rn_loop(index1_loop, 1); index2 = max(rn_loop(index1_loop, 2), index2_loop);
                if (index1_loop > 1); index1p = rn_loop(index1_loop - 1, 1); index2p = max(rn_loop(index1_loop - 1, 2), index2_loop); end
                if (index1_loop < length(index1_loop_l)); index1n = rn_loop(index1_loop + 1, 1); index2n = max(rn_loop(index1_loop + 1, 2), index2_loop); end
                % Clear first reset pixels if they don't achieve continuity requirements               
                if (index1_loop >  1); if(abs(u(index1, index2, index) - u(index1p, index2p, index)) > u_max_increment(index)); u(index1p, index2p, index) = NaN; error_u(index1p, index2p, index) = 1;  end; end
                % Clear lonely pixels to avoid issues in strain computation
                if (error_u(index1, index2, index) ~= 1);
                    if (index1_loop == 1); if (error_u(index1n, index2n, index) == 1); u(index1, index2, index) = NaN; error_u(index1, index2, index) = 1; end 
                    elseif (index1_loop == length(index1_loop_l)); if(error_u(index1p, index2p, index) == 1); u(index1, index2, index) = NaN; error_u(index1, index2, index) = 1; end 
                    else if((error_u(index1p, index2p, index) == 1) && (error_u(index1n, index2n, index) == 1)); u(index1, index2, index) = NaN; error_u(index1, index2, index) = 1; end; end; end    
            end; end; end
    
    % Optimize information representation - 2D function
    NaN_cancel_active = 1; % Cancel NaNs obtained at image interpolation with discontinuous deformation fields(may slightly corrupt statistics at edge values)   
    rn_nodec = struct('x1', []); for index = 1:NDIMENSIONS, rn_nodec.(['x', num2str(index)]) = rn_start.(['x', num2str(index)]):rn_end.(['x', num2str(index)]); end
    [rn_nodec_x2gi, rn_nodec_x1gi] = meshgrid(rn_nodec.x2, rn_nodec.x1);
    rni = struct('x1', []); for index = 1:NDIMENSIONS, 
            rni.(['x', num2str(index)]) = transp(((1:size(input2, index)) + (rn_start.(['x', num2str(index)]) - rn_starteff.(['x', num2str(index)]) - 1))/ov(index) + 1);       end
    [u_exp] = ANRCAlgorithmExpandArray(u, NaN_cancel_active); ui = zeros(length(rni.x1), length(rni.x2), NDIMENSIONS); 
    for index = 1:NDIMENSIONS, ui(1:end, 1:end, index) = ANRCAlgorithmNaNFreeInterp2(u_exp(1:end,1:end, index), rni.x2 + NaN_cancel_active, transp(rni.x1) + NaN_cancel_active, 'linear', NaN_cancel_active);end
    yd_hR_corr = ANRCAlgorithmNaNFreeInterp2(input1, rn_nodec_x2gi + ui(1:end,1:end, 2) - rn_start.x2 + 1, rn_nodec_x1gi + ui(1:end,1:end, 1) - rn_start.x1 + 1, 'spline', NaN_cancel_active);   

    if (decimate_output ~= 1) % Interpolate datasets to fit grid of original images
        [error_u_exp] =  ANRCAlgorithmExpandArray(error_u, NaN_cancel_active); [error_code_exp] =  ANRCAlgorithmExpandArray(error_code, NaN_cancel_active); error_ui = zeros(size(ui)); error_codei = zeros(size(ui));
        for index = 1:NDIMENSIONS, error_ui(1:end, 1:end, index) = ANRCAlgorithmNaNFreeInterp2(error_u_exp(1:end,1:end, index), rni.x2 + NaN_cancel_active, transp(rni.x1) + NaN_cancel_active, 'nearest', NaN_cancel_active);
            error_codei(1:end, 1:end, index) = ANRCAlgorithmNaNFreeInterp2(error_code_exp(1:end,1:end, index), rni.x2 + NaN_cancel_active, transp(rni.x1) + NaN_cancel_active, 'nearest', NaN_cancel_active); end
        output(1:end, 1:end, 1) = input2; output(1:end, 1:end, 2) = input1; output(1:end, 1:end, 3) = yd_hR_corr; output(1:end,1:end, 4) = ui(1:end, 1:end, 1); output(1:end, 1:end, 5) = ui(1:end, 1:end, 2);
        output(1:end, 1:end, 6) = error_ui(1:end, 1:end, 1); output(1:end, 1:end, 7) = error_ui(1:end, 1:end, 2); output(1:end, 1:end, 8) = max(error_codei, [], 3);    
    else % Decimate output datasets according to computation step (perform overlapped averaging to reduce noise)
        yu_crop = input2; yd_crop = input1;
        yu_crop = yu_crop(rn_starteff.x1 - rn_start.x1 + 1:end - (rn_end.x1 - rn_endeff.x1), rn_starteff.x2 - rn_start.x2 + 1:end - (rn_end.x2 - rn_endeff.x2)); 
        yd_crop = yd_crop(rn_starteff.x1 - rn_start.x1 + 1:end - (rn_end.x1 - rn_endeff.x1), rn_starteff.x2 - rn_start.x2 + 1:end - (rn_end.x2 - rn_endeff.x2));
        yd_hR_corr_crop = yd_hR_corr(rn_starteff.x1 - rn_start.x1 + 1:end - (rn_end.x1 - rn_endeff.x1), rn_starteff.x2 - rn_start.x2 + 1:end - (rn_end.x2 - rn_endeff.x2));       
        fun = @(block_struct) mean2(block_struct.data);
        yu_dec = blockproc(yu_crop, [ov(1) ov(2)], fun); yu_dec_nonan = yu_crop(1:ov(1):end, 1:ov(2):end);
        if (NaN_cancel_active == 1); yu_dec((isnan(yu_dec) == 1) & (isnan(yu_dec_nonan) == 0)) = yu_dec_nonan((isnan(yu_dec) == 1) & (isnan(yu_dec_nonan) == 0)); end % At the edge pixels, no overlapped averaging, just decimation, in order to avoid NaNs
        yd_dec = blockproc(yd_crop, [ov(1) ov(2)], fun); yd_dec_nonan = yd_crop(1:ov(1):end, 1:ov(2):end);
        if (NaN_cancel_active == 1); yd_dec((isnan(yd_dec) == 1) & (isnan(yd_dec_nonan) == 0)) = yd_dec_nonan((isnan(yd_dec) == 1) & (isnan(yd_dec_nonan) == 0)); end
        yd_hR_corr_dec = blockproc(yd_hR_corr_crop, [ov(1) ov(2)], fun); yd_hR_corr_dec_nonan = yd_hR_corr_crop(1:ov(1):end, 1:ov(2):end);
        if (NaN_cancel_active == 1); yd_hR_corr_dec((isnan(yd_hR_corr_dec) == 1) & (isnan(yd_hR_corr_dec_nonan) == 0)) = yd_hR_corr_dec_nonan((isnan(yd_hR_corr_dec) == 1) & (isnan(yd_hR_corr_dec_nonan) == 0));  end
        output(1:end,1:end, 1) = yu_dec; output(1:end,1:end, 2) = yd_dec; output(1:end,1:end, 3) = yd_hR_corr_dec; output(1:end,1:end, 4) = u(1:end, 1:end, 1); output(1:end,1:end, 5) = u(1:end, 1:end, 2);
        output(1:end,1:end, 6) = error_u(1:end, 1:end, 1); output(1:end,1:end, 7) = error_u(1:end, 1:end, 2); output(1:end,1:end, 8) = max(error_code, [], 3);
    end
    
    % Close monitoring windows
    for index = 1:NDIMENSIONS, if (monitoring_loop(index) == 1); close(fhmonitor_loop.(['x', num2str(index)])); end; end
end


function [yucorr, ydcorr] = ANRCAlgorithmPrecorr(yucorr, ydcorr, NDIMENSIONS)  
% function [yucorr, ydcorr] = ANRCAlgorithmPrecorr(yucorr, ydcorr, NDIMENSIONS)  
% Randomizes images by substracting mean along dim2 and normalizing to std deviation along dim 2
   if(NDIMENSIONS == 2);       
       meanD = NaN*ones([size(ydcorr, 1), 1]); stdD = meanD; meanU = meanD; stdU = meanD;
       for index = 1:size(ydcorr, 1),
           ydcorr1 = ydcorr(index, 1:end); ydcorr1_nonan = ydcorr1(isnan(ydcorr1) == 0); yucorr1 = yucorr(index, 1:end); yucorr1_nonan = yucorr1(isnan(yucorr1) == 0); 
           if(numel(ydcorr1_nonan) > 0); meanD(index) = mean(ydcorr1_nonan); stdD(index) = std(ydcorr1_nonan); end;
           if(numel(yucorr1_nonan) > 0); meanU(index) = mean(yucorr1_nonan); stdU(index) = std(yucorr1_nonan); end; end           
       meanD_matrix = meanD*ones(1, size(ydcorr,2)); stdD_matrix = stdD*ones(1, size(ydcorr, 2)); ydcorr = (ydcorr - meanD_matrix)./stdD_matrix; ydcorr(stdD_matrix == 0) = 0;       
       meanU_matrix = meanU*ones(1, size(yucorr,2)); stdU_matrix = stdU*ones(1, size(yucorr, 2)); yucorr = (yucorr - meanU_matrix)./stdU_matrix; yucorr(stdU_matrix == 0) = 0; end
end
   
function outoflimits = ANRCAlgorithmOutOfLimits(vector, matrix)
% function outoflimits = ANRCAlgorithmOutOfLimits(vector, matrix)
% Checks whether the components of a vector are outside the index range of a matrix
    size_matrix = size(matrix); outoflimits = 0;    
    for index = 1:length(vector),if ((vector(index) < 1) || (vector(index) > size_matrix(index))); outoflimits = 1; break;  end; end
end

function  [corrpos, corrvalue, p] = ANRCAlgorithmDIC2D(yu, yd, ru, rd, subsetu, subsetd, mode, monitoring, fh_monitoring, RTpauseTime)
% function  [corrpos, corrvalue, p] = ANRCAlgorithmDIC2D(yu, yd, ru, rd, subsetu, subsetd, mode, monitoring, fh_monitoring, RTpauseTime)
% Performs Digital Image Correlation to find deformation vector at specified subsets
% yu, yd: image data (undeformed, deformed)
% ru, rd: center of search windows (2x1)
% subsetu, subsetd: length of subsets (2x1)
% monitoring (1): activated, show monitor info in figure fh_monitoring
% mode: Algorithm used -> 1: DIC order 0, not interpolate datasets
% corrpos: deformation vector, corrvalue: correlation coefficient
% p: correlation parameters 
%   for mode = 1 -> p.atedge (2x1): is correlation at edge of search interval in dimensions 1, 2 (1: yes, 0: no), p.c_c: correlation function
W = subsetd - subsetu; u_track = rd - ru;
if (mode == 1) % DIC ORDEN 0 - NOT INTERPOLATE DATASETS
    % Subset coordinates
    x1u_subset = ru(1) - subsetu(1):ru(1) + subsetu(1); x2u_subset = ru(2) - subsetu(2):ru(2) + subsetu(2);
    x1d_subset = rd(1) - subsetd(1):rd(1) + subsetd(1); x2d_subset = rd(2) - subsetd(2):rd(2) + subsetd(2);
    % Get data
    yu_subset = yu(x1u_subset, x2u_subset); yd_subset = yd(x1d_subset, x2d_subset);    
    % Perform correlation
    c = ANRCAlgorithmNormxcorr2(yu_subset, yd_subset); % Search for yu_subset in yd_subset   
    % Search maximum in interest area
    lagso1 = (size(c,1) - 1)/2 + 1; lagso2 = (size(c,2) - 1)/2 + 1;    
    c_c = c((lagso1 - W(1)):(lagso1+ W(1)), (lagso2 - W(2)):(lagso2 + W(2)));    
    lagso_cc1 = W(1) + 1; lagso_cc2 = W(2) + 1;    
    [lagsox2, lagsox1] = meshgrid((lagso_cc2 - W(2)):(lagso_cc2 + W(2)), (lagso_cc1 - W(1)):(lagso_cc1 + W(1)) );    
    [last_corr, maxcorr] = max(c_c(:)); [maxcorrx1, maxcorrx2] = ind2sub(size(c_c), maxcorr(1));    
    maxcorrx1_2 = lagsox1(maxcorrx1, maxcorrx2) -lagso_cc1 + u_track(1); maxcorrx2_2 = lagsox2(maxcorrx1, maxcorrx2) -lagso_cc2 + u_track(2);  
    corrpos  = [maxcorrx1_2, maxcorrx2_2]; corrvalue = last_corr;
    
    if (monitoring == 1); figure(fh_monitoring); subplot(1,3,1); imshow(yu_subset, []); colormap bone; title('undeformed'); %caxis([0 1]);
        subplot(1,3,2); imshow(yd_subset, []); colormap jet; title('deformed');%caxis([0 1]);
        subplot(1,3,3); imshow(c_c, []); colormap jet; title({['r = ', num2str(last_corr)]; ['(', num2str(maxcorrx1_2), ',' num2str(maxcorrx2_2),')']}); ANRCAlgorithmPause(RTpauseTime); end
    
    p = struct('atedge', [], 'c_c', [], 'maxcorrx1', [], 'maxcorrx2', [], 'mode', 1, 'u_track', u_track); 
    p.atedge = [0, 0]; if ((maxcorrx1 == 1) || (maxcorrx1 == size(c_c, 1))); p.atedge(1) = 1; end   
    if ((maxcorrx2 == 1) || (maxcorrx2 == size(c_c, 2))); p.atedge(2) = 1; end % At edge from search window    
    p.c_c = c_c; p.maxcorrx1 = maxcorrx1; p.maxcorrx2 = maxcorrx2; 
else
    corrpos = [NaN, NaN]; corrvalue = 0;  p = [];
end
end

function  [corrposi, corrvaluei, pi] = ANRCAlgorithmDIC2DInt(S, p)
% function  [corrposi, corrvaluei, pi] = ANRCAlgorithmDIC2DInt(S, p)
% Performs interpolation of correlation function by factor S 
% S: interpolation factor
% p: correlation parameters obtained from ANRCAlgorithmDIC2D
% corrposi: interpolated deformation vector, corrvaluei: interpolated correlation coefficient
% pi: output interpolation parameters
%          for mode = 1 -> p.c_ci: interpolated correlation function
    u_track = p.u_track;

    if (p.mode == 1)
        base = sqrt(S); if (abs(round(base) - base) == 0); exponentS = log(S)/log(base); % Here one could define iterative refining by chosing smaller base
        else base = S; exponentS = log(S)/log(base); end % To be improved
        lagso_cc1 = (size(p.c_c, 1) - 1)/2 + 1; lagso_cc2 = (size(p.c_c, 2) - 1)/2 + 1;
        lags_maxcorrx2 = p.maxcorrx2; lags_maxcorrx1 = p.maxcorrx1;
        for indexS = 1:exponentS, % Iterative refining of corr. function
            lagsox2i = (lags_maxcorrx2 - 1/(base^(indexS - 1))):1/(base^(indexS)):(lags_maxcorrx2 + 1/(base^(indexS - 1)));
            lagsox1i = (lags_maxcorrx1 - 1/(base^(indexS - 1))):1/(base^(indexS)):(lags_maxcorrx1 + 1/(base^(indexS - 1)));
            [lagsox2i_g, lagsox1i_g] = meshgrid(lagsox2i, lagsox1i);
            c_ci = interp2(p.c_c, lagsox2i_g, lagsox1i_g, 'spline');
            [last_corri, maxcorri] = max(c_ci(:)); [maxcorrx1i, maxcorrx2i] = ind2sub(size(c_ci), maxcorri(1));
            lags_maxcorrx2 = lagsox2i(maxcorrx2i);  lags_maxcorrx1 = lagsox1i(maxcorrx1i);
        end
        maxcorrx1_2i = lags_maxcorrx1 -lagso_cc1 + u_track(1);  maxcorrx2_2i = lags_maxcorrx2 -lagso_cc2 + u_track(2);
        
        pi = struct('c_ci', []); pi.c_ci = c_ci; corrposi  = [maxcorrx1_2i, maxcorrx2_2i]; corrvaluei = last_corri;        
    else
        corrposi = [NaN, NaN]; corrvaluei = 0; pi = [];
    end
end

function ANRCAlgorithmPause(time)
% function ANRCAlgorithmPause(time)
% Slight modification of the pause time of Matlab. If time is positive,
% waits for that time. If time is negative stops and waits for keyboard
% input
    if (time >= 0); pause(time); else  pause; end; end

function c = ANRCAlgorithmNormxcorr2(yu, yd)
% function c = ANRCAlgorithmNormxcorr2(yu, yd)
% normxcorr2 including edge management. Out of sample values are not considered in correlation. Out of sample values are defined in yu (undeformed) and yd
% (deformed) as "NaN" values    
    numel0yu = numel(yu(isnan(yu) == 1)); numel0yd = numel(yd(isnan(yd) == 1)); % Check for NaN values    
    % yu = 2L + 1 % yd = 2M + 1 % c is a matrix of size (L + M)*2 + 1
    if ((numel0yu > 0) || (numel0yd > 0)) % Out of edge values       
        L1 = (size(yu, 1) - 1)/2; L2 = (size(yu, 2) - 1)/2; M1 = (size(yd, 1) - 1)/2; M2 = (size(yd, 2) - 1)/2; mid_c1 = (L1 + M1) + 1; mid_c2 = (L2 + M2) + 1;
        c = zeros((L1 + M1)*2 + 1, (L2 + M2)*2 + 1);  
        for dispu = -(M1 - L1):(M1 - L1),
            for dispv = -(M2 - L2):(M2 - L2),
                yd_c = yd((M1 + 1) + dispu - L1:(M1 + 1) + dispu + L1, (M2 + 1) + dispv - L2:(M2 + 1) + dispv + L2); yu_c = yu;
                yd_c(isnan(yu_c) == 1) = NaN; yu_c(isnan(yd_c) == 1) = NaN; % Set to same number of 0s in yd_c and yu_c 
                sum_yd_c = sum(sum(yd_c(isnan(yd_c) == 0))); sum_yu_c = sum(sum(yu_c(isnan(yu_c) == 0)));
                if ((numel(yd_c(isnan(yd_c) == 0)) == 0) || (numel(yu_c(isnan(yu_c) == 0)) == 0)); c(mid_c1 + dispu, mid_c2 + dispv) = 0; % No info available
                else numel_yd_c_NoNaN = numel(yd_c(isnan(yd_c) == 0)); numel_yu_c_NoNaN = numel(yu_c(isnan(yu_c) == 0));
                    yd_c(isnan(yd_c) == 1) = 0; yu_c(isnan(yu_c) == 1) = 0; % Set all NaNs to 0                    
                    yd_c_mean = sum_yd_c/numel_yd_c_NoNaN; yu_c_mean = sum_yu_c/numel_yu_c_NoNaN; yd_c_energy = sum(sum((yd_c - yd_c_mean).^2));
                    yu_c_energy = sum(sum((yu_c - yu_c_mean).^2)); yd_yu_cross = sum(sum((yd_c - yd_c_mean).*(yu_c - yu_c_mean)));
                    c(mid_c1 + dispu, mid_c2 + dispv) = yd_yu_cross/sqrt(yd_c_energy*yu_c_energy);
                end; end; end;                
    else c = normxcorr2(yu, yd); end % Use built-in function, fastest
end


function  [NUSED_EFF, NUSED_O] = ANRCAlgorithmGetNusedPixels(NUSED, index1_loop, index2_loop, count_nonused, rn_loop)
% function  [NUSED_EFF, NUSED_O] = ANRCAlgorithmGetNusedPixels(NUSED, index1_loop, index2_loop, count_nonused, rn_loop)
% Gets NUSED pixels including actual and previous loop iterations selecting only usable pixels (count_nonused = 1 in non-usable pixels).
% NUSED_EFF gives the loop index where the last pixel is taken from (index1_loop - NUSED_EFF)
% NUSED_O is the number of usable pixels returned, if sufficient pixels in loop NUSED_O = NUSED, otherwise NUSED_O < NUSED
NUSED_EFF = -1; NUSED_O = 0;
while (NUSED_O < NUSED);
    NUSED_EFF = NUSED_EFF + 1; if (index1_loop - NUSED_EFF <= 1); break; end; % If there are no more pixels go out
    index1t = rn_loop(index1_loop - NUSED_EFF, 1); index2t = max(rn_loop(index1_loop - NUSED_EFF, 2), index2_loop); if (count_nonused(index1t, index2t) == 0); NUSED_O = NUSED_O + 1; end;
end; 
end

function [input_exp] = ANRCAlgorithmExpandArray(input, activate)
% function [input_exp] = ANRCAlgorithmExpandArray(input, activate)    
% Expands input by one pixels on all edges, by copying actual edge values.
% Avoids NaNs in interpolation functions, at the cost of small corruption of statistics at edge values (activate = 1)
% If one does want to preserve the NaNs, set flag activate = 0
if (activate == 1)
    input_exp = zeros(size(input, 1) + 2, size(input, 2) + 2, size(input, 3));
    input_exp(2:end-1,2:end-1,:) = input(1:end,1:end,:); input_exp(end,2:end-1,:) = input_exp(end-1,2:end-1,:); input_exp(2:end-1,end,:) = input_exp(2:end-1,end-1,:); input_exp(end, end,:) = input_exp(end-1, end-1,:);
    input_exp(1,2:end-1, :) = input_exp(2,2:end-1, :); input_exp(2:end-1,1,:) = input_exp(2:end-1,2,:); input_exp(end, end,:) = input_exp(end-1, end-1,:); input_exp(1,1,:) = input_exp(2,2,:);
    input_exp(end,1,:) = input_exp(end - 1,2,:); input_exp(1,end,:) = input_exp(2,end - 1,:);
else    input_exp = input; end
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

function  ANRCAlgorithmLoop2Monitor(yd, rn_crop, index2p, u, fhmonitor_loop, yu, error_code, error_u, count_error, count_reset, u_track, W, count_indexloop, RTpauseTime)
% function  ANRCAlgorithmLoop2Monitor(yd, rn_crop, index2p, u, fhmonitor_loop, yu, error_code, error_u, count_error, count_reset, u_track, W, count_indexloop, RTpauseTime)
% Displays monitor info for dimension 2 loop     
yd_corr_temp = interp2(yd,  rn_crop.x2(index2p) + u(1:end, index2p, 2), rn_crop.x1(1:end) + u(1:end, index2p) , 'linear');
%yd_corr_temp = interp1(yd(1:end, index2p), transp(rn_loop.x1(1:end) + u(1:end, index2p)));
figure(fhmonitor_loop.x2); subplot(13,1,1); plot(rn_crop.x1 - rn_crop.x1(1) + 1, yu(rn_crop.x1(1:end), rn_crop.x2(index2p)), 'k'); ...
    hold on; plot(rn_crop.x1 - rn_crop.x1(1) + 1, yd(rn_crop.x1(1:end), rn_crop.x2(index2p)), 'b');
hold on; plot(rn_crop.x1 - rn_crop.x1(1) + 1, yd_corr_temp, 'r'); hold off; xlim([1, rn_crop.x1(end) - rn_crop.x1(1) + 1]);
title({['Colum ', num2str(index2p), ' of ', num2str(length(rn_crop.x2))]; 'Black: Undeformed, Blue: Deformed, Red: corrected'});
subplot(13,1,2); plot(rn_crop.x1 - rn_crop.x1(1) + 1, u(1:end, index2p, 1), 'k');  ylabel('u (pixels)'); hold off; xlim([1, rn_crop.x1(end) - rn_crop.x1(1) + 1]);
subplot(13,1,3); plot(rn_crop.x1 - rn_crop.x1(1) + 1, u(1:end, index2p, 2), 'k');  ylabel('v (pixels)'); hold off; xlim([1, rn_crop.x1(end) - rn_crop.x1(1) + 1]);
subplot(13,1,4); plot(rn_crop.x1 - rn_crop.x1(1) + 1, error_code(1:end, index2p, 1), 'k');  ylabel('Error code u (#)'); hold off; xlim([1, rn_crop.x1(end) - rn_crop.x1(1) + 1]);
subplot(13,1,5); plot(rn_crop.x1 - rn_crop.x1(1) + 1, error_code(1:end, index2p, 2), 'k');  ylabel('Error code v (#)'); hold off;xlim([1, rn_crop.x1(end) - rn_crop.x1(1) + 1]);
subplot(13,1,6); plot(rn_crop.x1 - rn_crop.x1(1) + 1, 1 - error_u(1:end, index2p, 1), 'k');  ylabel('1- Error_u (#)'); hold off; xlim([1, rn_crop.x1(end) - rn_crop.x1(1) + 1]);
subplot(13,1,7); plot(rn_crop.x1 - rn_crop.x1(1) + 1, 1 - error_u(1:end, index2p, 2), 'k');  ylabel('1- Error_v (#)'); hold off; xlim([1, rn_crop.x1(end) - rn_crop.x1(1) + 1]);
subplot(13,1,8); plot(rn_crop.x1 - rn_crop.x1(1) + 1, count_error(1:end, index2p), 'k'); ylabel('Error count (#)'); hold off; xlim([1, rn_crop.x1(end) - rn_crop.x1(1) + 1]);
subplot(13,1,9); plot(rn_crop.x1 - rn_crop.x1(1) + 1, count_reset(1:end, index2p), 'k');  ylabel('Reset count (#)'); hold off; xlim([1, rn_crop.x1(end) - rn_crop.x1(1) + 1]);
subplot(13,1,10); plot(rn_crop.x1 - rn_crop.x1(1) + 1, u_track(1:end, index2p, 1), 'k'); ylabel('TRK u (#)'); hold off; xlim([1, rn_crop.x1(end) - rn_crop.x1(1) + 1]);
subplot(13,1,11); plot(rn_crop.x1 - rn_crop.x1(1) + 1, u_track(1:end, index2p, 2), 'k');  ylabel('TRK v (#)'); hold off; xlim([1, rn_crop.x1(end) - rn_crop.x1(1) + 1]);
subplot(13,1,12); plot(rn_crop.x1 - rn_crop.x1(1) + 1, W(1:end, index2p, 1), 'k');  ylabel('W u (#)'); hold off; xlim([1, rn_crop.x1(end) - rn_crop.x1(1) + 1]);
subplot(13,1,13); plot(rn_crop.x1 - rn_crop.x1(1) + 1, count_indexloop(1:end, index2p, 1), 'k'); xlabel('Sample index (x1)'); ylabel('Loop 1'); hold off; xlim([1, rn_crop.x1(end) - rn_crop.x1(1) + 1]);
ANRCAlgorithmPause(RTpauseTime); disps(['Colum ', num2str(index2p), ' of  ', num2str(length(rn_crop.x2))]); 

end
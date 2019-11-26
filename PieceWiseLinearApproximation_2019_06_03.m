%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin constructing a Continuous Piece-wise Linear Approximation to the input MT-length history data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "Loop_Thru_Inputs.m" calls this file.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtain significant local extrema from the MT Length Data via findpeaks().
% Input parameter for the Matlab function 'findpeaks' to identify local max. 
% Use negative version data with 'findpeaks' to identify local min.

MIN_PROMINENCE_FOR_MAJOR_PEAKS = ERROR_TOLERANCE_LEVEL; % Use linear approximation error tolerance as peak prominence
MIN_PEAK_HEIGHT = NUC_HEIGHT_THRESHOLD; % min peak height. Any peaks below this height are ignored.

[major_peak_heights, major_peak_indices, major_peak_widths, major_peak_prominences] ...
    = findpeaks(MT_length, 'MinPeakHeight', MIN_PEAK_HEIGHT, 'MinPeakProminence', MIN_PROMINENCE_FOR_MAJOR_PEAKS );

[major_valley_heights, major_valley_indices, major_valley_widths, major_valley_prominences] ...
    = findpeaks(-MT_length, 'MinPeakProminence', MIN_PROMINENCE_FOR_MAJOR_PEAKS);
major_valley_heights = -major_valley_heights;

% Put Peaks and Valleys in chronological order in which they occur
peakcounter = 1;
valleycounter = 1;
peakvalleycounter = 1;
nuc_valley1 = 0;
nuc_valley2 = 0;
nuc_peak = 0;
PeakValley_indices = [];

while(peakcounter<=length(major_peak_indices) && valleycounter<length(major_valley_indices))
    %Save the first valley as the starting point
    PeakValley_indices(peakvalleycounter) = major_valley_indices(valleycounter);
    peakvalleycounter = peakvalleycounter + 1; % Seek next peak/valley
    
    %Find the last valley before next peak, and store it 
    nextvalley = valleycounter + 1;
    
    while(major_peak_indices(peakcounter)> major_valley_indices(valleycounter + 1))
        valleycounter = valleycounter + 1;
        if(valleycounter==length(major_valley_indices)) % prevent index out of bounds error
            break
        end
    end
    if (nextvalley~=valleycounter)
        PeakValley_indices(peakvalleycounter) = major_valley_indices(valleycounter);
        peakvalleycounter = peakvalleycounter + 1;
        valleycounter = valleycounter + 1;
    end 
    
    %Store all peaks prior to next valley
    while(major_peak_indices(peakcounter) < major_valley_indices(valleycounter) && peakcounter<length(major_peak_indices))
        PeakValley_indices(peakvalleycounter) = major_peak_indices(peakcounter);
        peakvalleycounter = peakvalleycounter + 1;
        peakcounter = peakcounter + 1;
    end
end
% Add in last peak and valley indices, and the very first and last data points
PeakValley_indices = [1, PeakValley_indices, ...
    min(major_peak_indices(end),major_valley_indices(end)), ...
    max(major_peak_indices(end),major_valley_indices(end)),...
    length(MT_length)];
% Add in first/last data points between stitched data, and stitch vertices
StitchedData_StartEnd_indicies = 1:length(MT_length);
StitchedData_StartEnd_indicies = StitchedData_StartEnd_indicies(MT_length<0);
Stitch_indicies = StitchedData_StartEnd_indicies;
StitchedData_StartEnd_indicies(1:2:length(StitchedData_StartEnd_indicies)) = StitchedData_StartEnd_indicies(1:2:length(StitchedData_StartEnd_indicies)) - 1;
StitchedData_StartEnd_indicies(2:2:length(StitchedData_StartEnd_indicies)) = StitchedData_StartEnd_indicies(2:2:length(StitchedData_StartEnd_indicies)) + 1;
PeakValley_indices = sort([PeakValley_indices,Stitch_indicies,StitchedData_StartEnd_indicies], 'ascend');


% Remove redundant vertices (i.e. same slope before and after vertex, and repeated indices)
PeakValley_indices = unique(PeakValley_indices);
delta_h = diff(MT_length(PeakValley_indices));
diff_delta_h = diff(delta_h);
Bad_PeakValley_list = (diff_delta_h==0);
while(sum(Bad_PeakValley_list)>0)
    Bad_PeakValley_list = [0;Bad_PeakValley_list(1:end)];
    PeakValley_indices(Bad_PeakValley_list==1) = []; % Remove bad peak/valley points
    delta_h = diff(MT_length(PeakValley_indices)); % Double-check for redundant points
    diff_delta_h = diff(delta_h);
    Bad_PeakValley_list = (diff_delta_h==0);
end
delta_h = 0; %flush delta_h variable

% Keep only Peaks/Valleys points that are more than MIN_TIME_STEP apart
PeakValleyTimeSteps = diff(Time(PeakValley_indices));
while(min(PeakValleyTimeSteps)<MIN_TIME_STEP)
    if(PeakValleyTimeSteps(end)<MIN_TIME_STEP)
        PeakValley_indices(end-1) = [];
        PeakValleyTimeSteps = diff(Time(PeakValley_indices));
    end
    Bad_PeakValley_list = PeakValleyTimeSteps<MIN_TIME_STEP; % Find Peak/Valley pairs that are too close
    Bad_PeakValley_list = [0;Bad_PeakValley_list(1:end-1)]; % Identify the latter of points that are too close
    PeakValley_indices(Bad_PeakValley_list==1) = []; % Remove bad peak/valley points
    PeakValleyTimeSteps = diff(Time(PeakValley_indices)); % Double-check for close points
end


% Identify new veritices for entry into/exit out of nucleation levels. 
% These vertices can over-ride peak/vallies with regard to min-time steps
Num_PeakValley_indices = length(PeakValley_indices);
GoodPeakValley_indicies = PeakValley_indices;
for(i = 1:Num_PeakValley_indices-1)
    start_index = PeakValley_indices(i);
    end_index = PeakValley_indices(i+1);
    vertex_now = [start_index,end_index];
    nuc_vertex_index = 1;
        
    % If ONLY start OR end vertex is below nucleation level, then find nucleation entry/exit vertex    
    if(MT_length(start_index)<=NUC_HEIGHT_THRESHOLD && MT_length(end_index)>NUC_HEIGHT_THRESHOLD)
        MT_length_start2end = MT_length(start_index:end_index);
        [nuc_vertex_hvalue,nuc_vertex_index] = max(~(MT_length_start2end<NUC_HEIGHT_THRESHOLD));
    elseif(MT_length(start_index)>NUC_HEIGHT_THRESHOLD && MT_length(end_index)<=NUC_HEIGHT_THRESHOLD)
        MT_length_start2end = MT_length(start_index:end_index);
        [nuc_vertex_hvalue,nuc_vertex_index] = max(MT_length_start2end<=NUC_HEIGHT_THRESHOLD);
    end
    
    % Check min-time step criteria if nucleation entry/exit vertex was added
    if(nuc_vertex_index~=1)
        if(Time(nuc_vertex_index+start_index-1) - Time(start_index) < MIN_TIME_STEP)
            vertex_now = [(nuc_vertex_index+start_index-1),end_index];
        elseif(Time(end_index) - Time(nuc_vertex_index) < MIN_TIME_STEP)
            vertex_now = [start_index,(nuc_vertex_index+start_index-1)];
        else
            vertex_now = [start_index,(nuc_vertex_index+start_index-1),end_index];
        end
    end

    % Update vertex list
    GoodPeakValley_indicies = [GoodPeakValley_indicies(GoodPeakValley_indicies<start_index),...
        vertex_now, GoodPeakValley_indicies(GoodPeakValley_indicies>end_index)];
end
% Reload Variable for peak valley indices, now with nucleation
PeakValley_indices = GoodPeakValley_indicies;


display('-----Start iterations in PW Linear Approximation');

% Start iterations of piecewise linear approximations:
% Begin by considering line segments connecting consecutive Peak/Valley
% points as the initial piecewise linear function to approximate
% to the MT length data. Next, identify nucleation entry/exit vertices 
% and avoid error minimization in these regions (behavior below nucleation 
% levels is not analyzed). To minimize errors including into the approximation 
% those data points that incur the highest error within the segment of 
% interest. Track and save the "vertex" points where the largest errors 
% are iteratively zeroed out, and their corresponding slope changes. These 
% will serve as possible separation/boundary points between DI phases.
vertex_all = PeakValley_indices(1);

Num_MT_length = length(MT_length);
linear_approx = zeros(Num_MT_length,1);
new_err = ERROR_TOLERANCE_LEVEL*100*ones(Num_MT_length,1);

Num_PeakValley_indices = length(PeakValley_indices);
for(i = 1:Num_PeakValley_indices-1)
    start_index = PeakValley_indices(i);
    end_index = PeakValley_indices(i+1);
    vertex_now = 0;
    max_err = ERROR_TOLERANCE_LEVEL*100;
    prior_err = ERROR_TOLERANCE_LEVEL*100;
    repeat_err_count = 0;

     while(max_err>ERROR_TOLERANCE_LEVEL) % repeat until error is less than tolerance
        vertex_now = [start_index, vertex_now(vertex_now>start_index & vertex_now<end_index),end_index];
        for(j=1:(length(vertex_now)-1))
            % Create linear segment using the points available until now
            m = (MT_length(vertex_now(j+1)) - MT_length(vertex_now(j)))/(Time(vertex_now(j+1)) - Time(vertex_now(j)));
            linear_approx(vertex_now(j):vertex_now(j+1)) = m*(Time(vertex_now(j):vertex_now(j+1)) - Time(vertex_now(j))) + MT_length(vertex_now(j));
            new_err(vertex_now(j):vertex_now(j+1)) = linear_approx(vertex_now(j):vertex_now(j+1)) - MT_length(vertex_now(j):vertex_now(j+1));
        end
        % Identify max error values and indices
        [max_err, maxerr_index] = max(abs(new_err(start_index:end_index)));
        
        if(max_err > ERROR_TOLERANCE_LEVEL) % If error criteria not met, identify new vertex for the next iteration
            new_vertex_index = maxerr_index + start_index - 1;
            
            % Check if new data vertex is at least the min time step away from the vertices already included
            new_time = Time(new_vertex_index);
            TempTime = Time(vertex_now);
            [last_time, last_time_index] = max(TempTime(TempTime<new_time));
            [next_time, next_time_index] = min(TempTime(TempTime>new_time));
            % Treat new vertices with time issue
            if(new_time - last_time < MIN_TIME_STEP && next_time - min(Time(Time>=last_time + MIN_TIME_STEP)) >= MIN_TIME_STEP)
                [new_vertex_time,new_vertex_index] = min(Time(Time>=last_time + MIN_TIME_STEP));
                display(['New vertex: Nudge to the right at time = ', num2str(new_vertex_time)])
            elseif(next_time - new_time < MIN_TIME_STEP && max(Time(Time<=next_time - MIN_TIME_STEP)) - last_time >= MIN_TIME_STEP)
                [new_vertex_time,new_vertex_index] = max(Time(Time<=next_time - MIN_TIME_STEP));
                display(['New vertex: Nudge to the left at time = ', num2str(new_vertex_time)])
            elseif(new_time - last_time < MIN_TIME_STEP && next_time - new_time < MIN_TIME_STEP)
                new_vertex_index = []; % Don't include problematic vertices
                display(['Problems with new vertex near time = ', num2str(new_time)])
            end
            
            % Incorporate new vertex into total vertex list
            vertex_now = sort([new_vertex_index, vertex_now(vertex_now>=start_index & vertex_now<=end_index)]);
            vertex_now = unique(vertex_now);
        end
        
        if(repeat_err_count == 5) % Catch repeated errors loops due to impossible tolerance level
            disp('Max error repeated 5 times');
            ErrMsg = ['i=', int2str(i), '; Error = ', int2str(max_err), '; index = ', int2str(maxerr_index), '; Time = ', num2str(new_time)];
            disp(ErrMsg);
            break;
        else
            if(prior_err == max_err)
                repeat_err_count = repeat_err_count + 1;
            else
                prior_err = max_err;
            end
        end
    end %end while
    
    %Update all the vertices from this section
    if(i==1)
        vertex_all = vertex_now(vertex_now>=start_index);
    else
        vertex_all = [vertex_all, vertex_now(vertex_now>start_index)];
    end
end % end iteration for loop 

% Plot MT length, the major peaks and valleys, and the Piece-wise Linear Approximation
if PLOT_FIG_1 == 1;
    fig1handle = figure;
    set(0,'units','pixels') % Resolution and screen size info needed for PNG image files
    Pixels= get(0,'screensize');
    set(0,'units','inches')
    Inches= get(0,'screensize');
    Res = Pixels/Inches;
    Papersize = [0 0 Pixels(3) 0.55*Pixels(4)];
    set(fig1handle, 'PaperUnits', 'inches', 'PaperPosition', Papersize/Res);   
    set(fig1handle, 'units','normalized','Position', [0 1 1 0.55]) 
    hold on
    % Omit artificial MT lengths from multi-sample stitch segments
    for(i=1:(length(vertex_all)-1))
        % Plot dummy lines & peaks/valleys to create legend first without issues from multiple plot commands
        if(i==1)
            plot(Time(1), MT_length(1), 'r-', 'LineWidth', 3)
            plot(Time(1), MT_length(1), 'b-', 'LineWidth', 3, 'Marker', '.', 'MarkerSize', 30)
            plot(Time(major_peak_indices),   major_peak_heights,   'd', 'LineWidth', 3, 'MarkerSize', 15, 'Color', [1.0, 0.6, 0.2])
            plot(Time(major_valley_indices), major_valley_heights, 's', 'LineWidth', 3, 'MarkerSize', 15, 'Color', [0.0, 1.0, 0.0])
            [hleg, hobj] = legend('Actual MT length', 'Approximate MT length', 'Major Peaks', 'Major Valleys');
            textobj = findobj(hobj, 'type', 'text');
            set(textobj,'fontsize', 16);
            set(hleg,'position',[0.2 0.6 0.2 0.25]);
        end
        if(MT_length(vertex_all(i))>=0 && MT_length(vertex_all(i+1))>=0)
            plot([Time(vertex_all(i)),Time(vertex_all(i+1))], [MT_length(vertex_all(i)),MT_length(vertex_all(i+1))], 'b-', 'LineWidth', 4, 'Marker', '.', 'MarkerSize', 30)
            plot(Time(Time>=Time(vertex_all(i)) & Time<=Time(vertex_all(i+1))), MT_length(Time>=Time(vertex_all(i)) & Time<=Time(vertex_all(i+1))), 'r-', 'LineWidth', 2)
        end
    end
    xlabel('Time (seconds)')
    ylabel('Microtubule Length (dimers)')
    set(gca, 'Fontsize', 16);
    axis([FIG_WINDOW_START,FIG_WINDOW_END,0,1.1*max(MT_length)])
    title({'Piece-wise Linear Approximation to MT Length History Plot'; ['Error Threshold = ',int2str(ERROR_TOLERANCE_LEVEL)]}, 'FontSize', 20);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    At this point, the Linear Approximation stage is complete    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

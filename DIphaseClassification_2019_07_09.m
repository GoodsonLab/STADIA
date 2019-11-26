%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for Classifying the Segments in the Approx. MT-length data into DI Phases %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Input_and_Run.m is the file that you run, and it calls the file   %
% Loop_Thru_Inputs.m, which calls this file.                        %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Library of DI phase identifiers:
    %  -12 = Long Shortening .. (Light Red)
    %  -11 = Brief Shortening . (Dark red)
    %   -1 = Down Stutter ..... (Purple)
    %    0 = Flat Stutter ..... (Blue)
    %    1 = Up Stutter ....... (Purple)
    %   11 = Brief Growth ..... (Dark Green)
    %   12 = Long Growth ...... (Light Green)
    %  NaN = Nucleation ....... (Gray)
    % -inf = Sample Stiches ... (Black)
    
Sh2_color = [1 0.3 0.3];    % Light Red for Long Shortening
Sh1_color = [0.8 0 0];      % Dark Red for Brief Shortening
DnSt_color = [0.4 0 1];     % Purple for Down Stutter
FlSt_color = [0.1 0.1 1];   % Blue for Flat Stutter
UpSt_color = [0 0.5 1];     % Aqua for Up Stutter
Gr1_color = [0 0.6 0];      % Light Green for Brief Growth
Gr2_color = [0.2 0.8 0.2];  % Dark Green for Long Growth
Nuc_color = [0.5 0.5 0.5];  % Gray for Nucleation
Stch_color = [1, 1, 1];     % Stitching between Multiple Samples 
Phase_Colors = [Gr2_color;Gr1_color;...
    UpSt_color;FlSt_color;DnSt_color;...
    Sh1_color;Sh2_color;...
    Nuc_color; Stch_color];


% Determine changes in height, segment duration, and slope of each segment 
% as defined between vertices. This is the 3D data used from this point on.
deltah_all = 0;
deltat_all = 0;
slope_all = 0;
for(i=1:length(vertex_all)-1)
    deltah_all(i) = (MT_length(vertex_all(i+1)) - MT_length(vertex_all(i)));
    deltat_all(i) = Time(vertex_all(i+1)) - Time(vertex_all(i));
    slope_all(i) = deltah_all(i)/deltat_all(i);
end
X_all = [deltat_all; deltah_all; slope_all]';

% User Input thresholds for Identifying Nucleation and Flat Stutter segments
min_deltah = FLATSTUT_MAX_SEGMENTHEIGHT_THRESHOLD; % Cutoff delta-h for Flat Stutters
min_slope = FLATSTUT_MAX_SEGMENTSLOPE_THRESHOLD; % Cutoff slope for Flat Stutters
Min_MT_Height = NUC_HEIGHT_THRESHOLD; % Min height for MT to be out of nucleation

% Label Nucleation and Flat Stutters in the "DI_phases" term/array.
% Also label segments used to stitch/separate multiple samples when
% more than one input file/column is used for the length history data.
DI_phases = 0.1*ones(length(X_all),1); % Initialize with value different from id values
for(i=1:length(DI_phases))
    if( MT_length(vertex_all(i))<0 || MT_length(vertex_all(i+1))<0 )
        DI_phases(i) = -inf;
    elseif( max(MT_length(vertex_all(i)),MT_length(vertex_all(i+1))) <= Min_MT_Height)
        DI_phases(i) = NaN;
    elseif( abs(deltah_all(i)) <= min_deltah || abs(slope_all(i)) <= min_slope)
        DI_phases(i) = 0;
    end
end

%%% ONLY CONSIDER NON-CLASSIFIED SEGMENTS (EXCLUDING NUCLEATION, FLAT STUTTER PHASES, and STITCH SEGMENTS)
exclusion_idx = ~isnan(DI_phases) & DI_phases>0;
deltah = deltah_all(exclusion_idx);
slope = slope_all(exclusion_idx);
deltat = deltat_all(exclusion_idx);

%%% Explore 3D Scatter-plot Excluding Nucleation and Flat Stutters
% Raw data
X = [deltat;deltah;slope]'; % Define Data: segment time, height changes, and segment slopes
dist2o = (X_all(:,1).^2 + X_all(:,2).^2 + X_all(:,3).^2).^0.5;

if(KMEANS_DIAGNOSTIC_AllFLAG == 1)
    DataFig = figure;
    scatter3(X_all(:,1),X_all(:,2),X_all(:,3),5,dist2o,'ko')
    caxis([0,1000])
	xlabel('Time')
    ylabel('Height')
    zlabel('Slope')
    set(gca, 'Fontsize', 16);
    title({'DI Segment Data'}, 'Fontsize', 24);
    colormap(copper)
    StndrdDataFig = figure;
    Y = sign(X).*log(abs(X)+1);
    Z = zscore(Y); %Just the values
    plot3(Z(:,1),Z(:,2),Z(:,3),'.')
    title({'Log Transformed then Standardized Data'}, 'Fontsize', 24);
    
    % Gap Statistics for k-means using ALL segments (to help choose best k value)
    ALL_eval = evalclusters(Z,'kmeans','Gap','KList',[1:12], 'SearchMethod', 'firstMaxSE');
    alldiagnosticfigure = figure;
    set(alldiagnosticfigure, 'units','normalized','Position', [0 1 1 0.55]) 
    subplot(1,2,1)
    plot(ALL_eval)
    set(gca, 'Fontsize', 16);
    title({'Gap Statistics for k-means on'; 'Log Scaled & Standardized Data'},'Fontsize',24)
    suggested_optimal_k_all = ALL_eval.OptimalK;
    
%    % Use User defined number of k-s to cluster All DI Segments
%    k_all=KMEANS_NumClust_All;
    
    % Classify the Log Transformed Data using k-means Clustering
%    [idx,C,~] = kmeans(Z,k_all,'Replicates',101);
    [idx,C,~] = kmeans(Z,suggested_optimal_k_all,'Replicates',101);
    figure(alldiagnosticfigure)
    subplot(1,2,2)
    hold all
    for(i=1:suggested_optimal_k_all)
        plot3(Z(idx==i,1),Z(idx==i,2),Z(idx==i,3),'.','MarkerSize',12)
        plot3(C(i,1),C(i,2),C(i,3),'kx','MarkerSize',15,'LineWidth',3)
    end
    xlabel('Time')
    ylabel('Height')
    zlabel('Slope')
    set(gca, 'Fontsize', 16);
    titletxt = [{['k-means Classification using k=',int2str(suggested_optimal_k_all)]; 'on Log Scaled & Standardized Data'}];
    title(titletxt, 'Fontsize', 24);
    grid on
    
    
    % Save Results of Diagnostic Mode for All Segments
    original_dir = pwd;
    
    % Remove filename suffix
    filename = sprintf('%s',char(FILE_NAME));
    filename_prefix = split(filename,".");
    filename_prefix = filename_prefix{1};
    % Remove directories
    filename_prefix = split(filename_prefix, "/"); % for Non-Windows style directory
    Dir_Suffix = filename_prefix{length(filename_prefix)};
    filename_prefix = split(filename_prefix, "\"); % for Windows style directory
    Dir_Suffix = filename_prefix{length(filename_prefix)};
    
    % Create & Move to new directory
    Dir_name = "Diagnostic_Mode_Outputs_for_" + Dir_Suffix;
    if ~exist(Dir_name, 'dir')
        mkdir(Dir_name);
    end
    cd(Dir_name);
    
    % Use date & time to save files
    this_datetime = datetime;
    DateTimeString = sprintf('%.0f_%02dh%.0fmin%.0fsec', ...
        yyyymmdd(this_datetime), hour(this_datetime), ...
        minute(this_datetime), round(second(this_datetime)));
    
    % Save diagnostic figure
    diagnosticALL_figure_filename = sprintf('Diagnostic_Mode_All_Segments_Figure_%s', DateTimeString);
    savefig(alldiagnosticfigure ,diagnosticALL_figure_filename)
    print(alldiagnosticfigure,'-dpng','-r0', diagnosticALL_figure_filename);
    close(alldiagnosticfigure);
    
    % Save Gap Statistic Results
    diagnosticALL_table_filename = sprintf('Diagnostic_Mode_All_Segments_Table_%s.txt', DateTimeString);
    k_values = ALL_eval.InspectedK';
    gap_values = ALL_eval.CriterionValues';
    StdErrors = ALL_eval.SE';
    diagnosticALL_table = table(k_values, gap_values, StdErrors);
    writetable(diagnosticALL_table, diagnosticALL_table_filename);
    
    % Return to original directory
    cd(original_dir);
    
    disp(' ');
    disp('***###***###***###***###***###***###***###***###***###***###***###***###***');
    disp('Diagnostic Mode ON for All Slope Segments');
    disp('Suggested k-value: ' + string(suggested_optimal_k_all));
    disp('Corresponding output folder: ' + Dir_name);
    disp('Corresponding output file time stamp: ' + string(DateTimeString));
    disp('***###***###***###***###***###***###***###***###***###***###***###***###***');
    disp(' ');
    
    % Generate error message for diagnostic mode
    error('Diagnostic Mode on for All Segments')
end

% Define Positive sloped subgroups for clustering application
% Use Log transform values for "plus"
intindx = [1:length(DI_phases)]'; 
plusindx = intindx(exclusion_idx & X_all(:,2)>0);
plus_idx_mat = repmat((X(:,2)>0),1,3);

Xplus = X_all(plusindx,:); % Positive slope segments only
Yplus = plus_idx_mat.*sign(X).*log(abs(X)+1); % Log Transform
Yplusvec = Yplus(Yplus~=0);
Yplusmat = [Yplusvec(1:length(Xplus)),Yplusvec(length(Xplus)+1:2*length(Xplus)),Yplusvec(2*length(Xplus)+1:3*length(Xplus))];
MU_PLUS_TRANSFORM = mean(Yplusmat); 
SIGMA_PLUS_TRANSFORM = std(Yplusmat);
Zplus = zscore(Yplus(plus_idx_mat(:,1),:)); % Standardized

if(KMEANS_DIAGNOSTIC_PosSlopeFLAG == 1)
    % Positive Slope data: Observe and Apply Log Transform
    RawPosDataFig = figure;
    scatter3(Xplus(:,1),Xplus(:,2),Xplus(:,3),5)
    title('Raw Positive Slope Data', 'Fontsize', 24);
    PosDataFig = figure;
    plot3(Zplus(:,1),Zplus(:,2),Zplus(:,3),'.')
    grid on
    title({'Log Scaled & Standardized'; 'Positive Slope Data'}, 'Fontsize', 24);
    
    % Gap Statistics for k-means (to help choose best k value)
    clustFcn = @(X,K) kmeans(X, K, 'Replicates',101);
    POS_eval = evalclusters(Zplus,'kmeans','Gap','KList',[1:12], 'SearchMethod', 'firstMaxSE');
    posdiagnosticfigure = figure;
    set(posdiagnosticfigure, 'units','normalized','Position', [0 1 1 0.55]) 
    subplot(1,2,1)
    plot(POS_eval)
    set(gca, 'Fontsize', 16);
    title({'Gap Statistics for k-means Classification'; 'on Log Scaled & Standardized'; 'Positive Slope Data'}, 'Fontsize', 24);
    suggested_optimal_k_pos = POS_eval.OptimalK;
end

% Assign number of k-s to cluster Positive Slope Segments 
if(KMEANS_DIAGNOSTIC_PosSlopeFLAG == 1)
    % Use sugested value when in diagnostic mode
    k_pos = suggested_optimal_k_pos;
else
    % Use User defined k-value
    k_pos=KMEANS_NumClust_PosSlope;
end

% Classify the Log Transformed Positive Slope data using k-means Clustering
[idx,C_pos,~] = kmeans(Zplus,k_pos,'Replicates',201);

if(KMEANS_DIAGNOSTIC_PosSlopeFLAG == 1)
    figure(posdiagnosticfigure)
    subplot(1,2,2)
elseif(PLOT_FIG_2 == 1)
    fig2handle = figure;
    set(0,'units','pixels') % Resolution and screen size info needed for PNG image files
    Pixels= get(0,'screensize');
    set(0,'units','inches')
    Inches= get(0,'screensize');
    Res = Pixels/Inches;
    Papersize = [0 0 Pixels(3) 0.55*Pixels(4)];
    set(fig2handle, 'PaperUnits', 'inches', 'PaperPosition', Papersize/Res);  
    set(fig2handle, 'units','normalized','Position', [0 1 1 0.55])
    subplot(1,2,1)
end
if(KMEANS_DIAGNOSTIC_PosSlopeFLAG == 1 || PLOT_FIG_2 == 1)
    hold all
    for(i=1:k_pos)
        plot3(Zplus(idx==i,1),Zplus(idx==i,2),Zplus(idx==i,3),'.','MarkerSize',12)
        plot3(C_pos(i,1),C_pos(i,2),C_pos(i,3),'kx','MarkerSize',15,'LineWidth',3)
    end
    xlabel('Time')
    ylabel('Height')
    zlabel('Slope')
    set(gca, 'Fontsize', 16);
    titletxt = [{['k-means Classification using k=',int2str(k_pos)]; 'for Log Scaled & Standardized'; 'Positive Slope Data'}];
    title(titletxt, 'Fontsize', 24);
    grid on
end

if(KMEANS_DIAGNOSTIC_PosSlopeFLAG == 1)
    % Save Results of Diagnostic Mode for All Segments
    original_dir = pwd;
    
    % Remove filename suffix
    filename = sprintf('%s',char(FILE_NAME));
    filename_prefix = split(filename,".");
    filename_prefix = filename_prefix{1};
    % Remove directories
    filename_prefix = split(filename_prefix, "/"); % for Non-Windows style directory
    Dir_Suffix = filename_prefix{length(filename_prefix)};
    filename_prefix = split(filename_prefix, "\"); % for Windows style directory
    Dir_Suffix = filename_prefix{length(filename_prefix)};
    
    % Create & Move to new directory
    Dir_name = "Diagnostic_Mode_Outputs_for_" + Dir_Suffix;
    if ~exist(Dir_name, 'dir')
        mkdir(Dir_name);
    end
    cd(Dir_name);
    
    % Use date & time to save files
    this_datetime = datetime;
    DateTimeString = sprintf('%.0f_%02dh%.0fmin%.0fsec', ...
        yyyymmdd(this_datetime), hour(this_datetime), ...
        minute(this_datetime), round(second(this_datetime)));
    
    % Save diagnostic figure
    diagnosticPOS_figure_filename = sprintf('Diagnostic_Mode_PositiveSlope_Segments_Figure_%s', DateTimeString);
    savefig(posdiagnosticfigure ,diagnosticPOS_figure_filename)
    print(posdiagnosticfigure,'-dpng','-r0', diagnosticPOS_figure_filename);
    close(posdiagnosticfigure);
    
    % Save Gap Statistic Results
    diagnosticPOS_table_filename = sprintf('Diagnostic_Mode_PositiveSlope_Segments_Table_%s.txt', DateTimeString);
    k_values = POS_eval.InspectedK';
    gap_values = POS_eval.CriterionValues';
    StdErrors = POS_eval.SE';
    diagnosticPOS_table = table(k_values, gap_values, StdErrors);
    writetable(diagnosticPOS_table, diagnosticPOS_table_filename);
    
    % Return to original directory
    cd(original_dir);
    
    disp(' ');
    disp('***###***###***###***###***###***###***###***###***###***###***###***###***');
    disp('Diagnostic Mode ON for Positive Slope Segments');
    disp('Suggested k-value: ' + string(suggested_optimal_k_pos));
    disp('Corresponding output folder: ' + Dir_name);
    disp('Corresponding output file time stamp: ' + string(DateTimeString));
    disp('***###***###***###***###***###***###***###***###***###***###***###***###***');
    disp(' ');
    
    error('Diagnostic Mode on for Positive Slope Segments')
end

% Update DI-phases for k=k_pos:
posindx = [1:k_pos];
pos_phase_idx = zeros(size(idx));
if(k_pos == 1) % if k_pos = 1: Identify 1 Growth Phase
    [gr_Cdeltah,gr_idx] = max(C_pos(:,3));
    pos_phase_idx(idx==gr_idx) = 12;
elseif(k_pos == 2) % if k_pos = 2: Identify (A) 2 Growth OR (B) 1 Up-Stutters and 1 Growth Phases
    if(KMEANS_Pos2_Option == 'A') % 1st option (2 Growth Phases, but no Up-Stutters)
        % Separate by segment time durations
        [~,gr1_idx] = min(C_pos(:,1));
        [~,gr2_idx] = max(C_pos(:,1)); 
        pos_phase_idx(idx==gr1_idx) = 11;
        pos_phase_idx(idx==gr2_idx) = 12;
    elseif(KMEANS_Pos2_Option == 'B') % 2nd option (1 Growth Phase, and 1 Up-Stutters Phase)
        % Separate by segment slope steepness (Recall that greater positive slope segment values are steeper)
        %[~,upst_idx] = min(C_pos(:,3));
        %[~,gr_idx] = max(C_pos(:,3)); 
        % Separate by segment time durations
        [~,upst_idx] = min(C_pos(:,1));
        [~,gr_idx] = max(C_pos(:,1)); 
        pos_phase_idx(idx==gr_idx) = 12;
        pos_phase_idx(idx==upst_idx) = 1;
    else
        error('Invalid option for k=2 for positive slope segments. Please use a valid choice for KMEANS_Pos2_Option');
    end
elseif(k_pos == 3) % if k_pos = 3: Identify 2 Growth Phases and Up-Stutters
    [gr2_Ctime, gr2_idx] = max(C_pos(:,1)); % First, identify long time durations cluster as Long Growth
    [upst_Cslope, ~] = min(C_pos((C_pos(:,1)~=gr2_Ctime),3)); % Second, choose the smaller of the remaining cluster slopes as stutter
    [~, upst_idx] = max(C_pos(:,3)==upst_Cslope);
    [gr1_Ctime, ~] = min(C_pos((C_pos(:,3)~=upst_Cslope & C_pos(:,1)~=gr2_Ctime),1)); % Finally, the remaining cluster is Breif Growth
    gr1_idx = posindx(C_pos(:,1) == gr1_Ctime);
    pos_phase_idx(idx==gr2_idx) = 12;
    pos_phase_idx(idx==gr1_idx) = 11;
    pos_phase_idx(idx==upst_idx) = 1;
else
    error('Invalid number of k-s used to cluster positive slope data. Please limit KMEANS_NumClust_PosSlope to 1, 2, or 3.');
end
DI_phases(plusindx) = pos_phase_idx;  %The phase ids indexed according to the complete data set

display('-----Positive Slope Segment Classification Complete');
toc


% Define Negative sloped subgroups for clustering application
% Use square root transform for "minus" labels
intindx = [1:length(DI_phases)]'; 
minusindx = intindx(exclusion_idx & X_all(:,2)<0);
minus_idx_mat = repmat((X(:,2)<0),1,3);

Xminus = X_all(minusindx,:); % Negative slope segments only
Yminus = minus_idx_mat.*sign(X).*log(abs(X)+1); % Log Transform
Yminusvec = Yminus(Yminus~=0);
Yminusmat = [Yminusvec(1:length(Xminus)),Yminusvec(length(Xminus)+1:2*length(Xminus)),Yminusvec(2*length(Xminus)+1:3*length(Xminus))];
MU_MINUS_TRANSFORM = mean(Yminusmat); 
SIGMA_MINUS_TRANSFORM = std(Yminusmat);
Zminus = zscore(Yminus(minus_idx_mat(:,1),:)); % Standardized

if(KMEANS_DIAGNOSTIC_NegSlopeFLAG == 1)
    % Negative Slope data: Observe Raw and Log Transform
    RawNegDataFig = figure;
    scatter3(Xminus(:,1),Xminus(:,2),Xminus(:,3),5)
    title('Raw Negative Slope Data', 'Fontsize', 24);
    NegDataFig = figure;
    plot3(Zminus(:,1),Zminus(:,2),Zminus(:,3),'.');
    grid on
    title({'Log Transformed then Standardized'; 'Negative Slope Data'}, 'Fontsize', 24);
    % Gap Statistics for k-means (to help choose best k value)
    clustFcn = @(X,K) kmeans(X, K, 'Replicates',101);
    NEG_eval = evalclusters(Zminus,'kmeans','Gap','KList',[1:12], 'SearchMethod', 'firstMaxSE');
    negdiagnosticfigure = figure;
    set(negdiagnosticfigure, 'units','normalized','Position', [0 1 1 0.55]) 
    subplot(1,2,1)
    plot(NEG_eval)
    set(gca, 'Fontsize', 16);
    title({'Gap Statistics for k-means Classification'; 'on Log Scaled & Standardized'; 'Negative Slope Data'}, 'Fontsize', 24);
    suggested_optimal_k_neg = NEG_eval.OptimalK;
end

% Assign number of k-s to cluster Negative Slope Segments 
if(KMEANS_DIAGNOSTIC_NegSlopeFLAG == 1)
    % Use sugested value when in diagnostic mode
    k_neg = suggested_optimal_k_neg;
else
    % Use User defined k-value
    k_neg=KMEANS_NumClust_NegSlope;
end

% Classify the Log Transformed Negative Slope data using k-means Clustering
[idx,C_neg,~] = kmeans(Zminus,k_neg,'Replicates',201);

if(KMEANS_DIAGNOSTIC_NegSlopeFLAG == 1)
    figure(negdiagnosticfigure)
    subplot(1,2,2)
elseif(PLOT_FIG_2 == 1)
    figure(fig2handle)
    subplot(1,2,2)
end
if(KMEANS_DIAGNOSTIC_NegSlopeFLAG == 1 || PLOT_FIG_2 == 1)
    hold all
    for(i=1:k_neg)
        plot3(Zminus(idx==i,1),Zminus(idx==i,2),Zminus(idx==i,3),'.','MarkerSize',12)
        plot3(C_neg(i,1),C_neg(i,2),C_neg(i,3),'kx','MarkerSize',15,'LineWidth',3)
    end
    xlabel('Time')
    ylabel('Height')
    zlabel('Slope')
    set(gca, 'Fontsize', 16);
    titletxt = [{['k-means Classification using k=',int2str(k_neg)]; 'for Log Scaled & Standardized'; 'Negative Slope Data'}];
    title(titletxt, 'Fontsize', 24);
    grid on
end
if(KMEANS_DIAGNOSTIC_NegSlopeFLAG == 1)
    % Save Results of Diagnostic Mode for All Segments
    original_dir = pwd;
    
    % Remove filename suffix
    filename = sprintf('%s',char(FILE_NAME));
    filename_prefix = split(filename,".");
    filename_prefix = filename_prefix{1};
    % Remove directories
    filename_prefix = split(filename_prefix, "/"); % for Non-Windows style directory
    Dir_Suffix = filename_prefix{length(filename_prefix)};
    filename_prefix = split(filename_prefix, "\"); % for Windows style directory
    Dir_Suffix = filename_prefix{length(filename_prefix)};
    
    % Create & Move to new directory
    Dir_name = "Diagnostic_Mode_Outputs_for_" + Dir_Suffix;
    if ~exist(Dir_name, 'dir')
        mkdir(Dir_name);
    end
    cd(Dir_name);
    
    % Use date & time to save files
    this_datetime = datetime;
    DateTimeString = sprintf('%.0f_%02dh%.0fmin%.0fsec', ...
        yyyymmdd(this_datetime), hour(this_datetime), ...
        minute(this_datetime), round(second(this_datetime)));
    
    % Save diagnostic figure
    diagnosticNEG_figure_filename = sprintf('Diagnostic_Mode_NegativeSlope_Segments_Figure_%s', DateTimeString);
    savefig(negdiagnosticfigure ,diagnosticNEG_figure_filename)
    print(negdiagnosticfigure,'-dpng','-r0', diagnosticNEG_figure_filename);
    close(negdiagnosticfigure);
    
    % Save Gap Statistic Results
    diagnosticNEG_table_filename = sprintf('Diagnostic_Mode_NegativeSlope_Segments_Table_%s.txt', DateTimeString);
    k_values = NEG_eval.InspectedK';
    gap_values = NEG_eval.CriterionValues';
    StdErrors = NEG_eval.SE';
    diagnosticNEG_table = table(k_values, gap_values, StdErrors);
    writetable(diagnosticNEG_table, diagnosticNEG_table_filename);
    
    % Return to original directory
    cd(original_dir);
    
    disp(' ');
    disp('***###***###***###***###***###***###***###***###***###***###***###***###***');
    disp('Diagnostic Mode ON for Negative Slope Segments');
    disp('Suggested k-value: ' + string(suggested_optimal_k_neg));
    disp('Corresponding output folder: ' + Dir_name);
    disp('Corresponding output file time stamp: ' + string(DateTimeString));
    disp('***###***###***###***###***###***###***###***###***###***###***###***###***');
    disp(' ');
    
    error('Diagnostic Mode on for Negative Slope Segments')
end

% Update DI-phases for k=k_neg:
% Identify Down-Stutters and 2 Shortening Phases
negindx = [1:k_neg];
neg_phase_idx = zeros(size(idx));
if(k_neg == 1) % if k_neg = 1: Identify 1 Shortening Phase
    [sh_Cdeltah,sh_idx] = max(C_neg(:,3));
    neg_phase_idx(idx==sh_idx) = -12;
elseif(k_neg == 2) % if k_neg = 2: Identify (A) 2 Shortening OR (B) 1 Down-Stutter and 1 Growth Phases
    if(KMEANS_Neg2_Option == 'A') % 1st option (choose if 2 Growth Phases, but no Up-Stutters)
        % Separate by segment time durations
        [~,sh1_idx] = min(C_neg(:,1));
        [~,sh2_idx] = max(C_neg(:,1)); 
        neg_phase_idx(idx==sh1_idx) = -11;
        neg_phase_idx(idx==sh2_idx) = -12;
    elseif(KMEANS_Neg2_Option == 'B') % 2nd option (choose if 1 Growth Phase, and Up-Stutters too)
        % Separate by segment slope steepness 
        % (Recall that lower negative slope segment values are steeper)
        [~,dnst_idx] = max(C_neg(:,3));
        [~,sh_idx] = min(C_neg(:,3)); 
        neg_phase_idx(idx==sh_idx) = -12;
        neg_phase_idx(idx==dnst_idx) = -1;
    else
        error('Invalid option for k=2 for negative slope segments. Please use a valid choice for KMEANS_Neg2_Option');
    end
elseif(k_neg == 3) % if k_neg = 3: Identify 2 Shortening Phases, and 1 Down-Stutter Phase
    [dnst_Cslope,dnst_idx] = max(C_neg(:,3)); % First, identify maximum slope closest to actual 0-slope as down stutter
    [sh2_Ctime,~] = max(C_neg((C_neg(:,3)~=dnst_Cslope),1)); % Second, identify long time duration of the remaining clusters as Long Shortening 
    [~,sh2_idx] = max(C_neg(:,1)==sh2_Ctime);
    [sh1_Ctime, ~] = min(C_neg((C_neg(:,3)~=dnst_Cslope & C_neg(:,1)~=sh2_Ctime),1)); % The remaining cluster is Breif Shortening
    [~,sh1_idx] = max(C_neg(:,1)==sh1_Ctime);
%     [sh2_Ctime,sh2_idx] = max(C_neg(:,1));
%     [dnst_Cheight,dnst_idx] = max(C_neg((C_neg(:,1)~=sh2_Ctime),2));
%     dnst_idx = negindx(C_neg(:,2) == dnst_Cheight);
%     [sh1_Cheight, ~] = min(C_neg((C_neg(:,2)~=dnst_Cheight & C_neg(:,1)~=sh2_Ctime),2));
%     sh1_idx = negindx(C_neg(:,2) == sh1_Cheight);
    neg_phase_idx(idx==sh2_idx) = -12;
    neg_phase_idx(idx==sh1_idx) = -11;
    neg_phase_idx(idx==dnst_idx) = -1;
else
    error('Invalid number of k-s used to cluster negative slope data. Please limit KMEANS_NumClust_NegSlope to 1, 2, or 3.');
end

DI_phases(minusindx,:) = neg_phase_idx;  %The phase ids indexed according to the complete data set

display('-----Negative Slope Segment Classification Complete');
toc


if(PLOT_FIG_3 == 1)
    display('---Start Figure 3');
    % Scatter plot of classified points
    fig3handle = figure;
    set(0,'units','pixels') % Resolution and screen size info needed for PNG image files
    Pixels= get(0,'screensize');
    set(0,'units','inches')
    Inches= get(0,'screensize');
    Res = Pixels/Inches;
    Papersize = [0 0 0.65*Pixels(3) Pixels(4)];
    set(fig3handle, 'PaperUnits', 'inches', 'PaperPosition', Papersize/Res);  
    set(fig3handle, 'units','normalized','Position', [0 1 0.65 1])
    hold all
    plot3(deltat_all(DI_phases==12),deltah_all(DI_phases==12),slope_all(DI_phases==12),'.','MarkerSize',12, 'color', Gr2_color)
    plot3(deltat_all(DI_phases==11),deltah_all(DI_phases==11),slope_all(DI_phases==11),'.','MarkerSize',12, 'color', Gr1_color)
    plot3(deltat_all(DI_phases==1),deltah_all(DI_phases==1),slope_all(DI_phases==1),'.','MarkerSize',12, 'color', UpSt_color)
    plot3(deltat_all(DI_phases==-12),deltah_all(DI_phases==-12),slope_all(DI_phases==-12),'.','MarkerSize',12, 'color', Sh2_color)
    plot3(deltat_all(DI_phases==-11),deltah_all(DI_phases==-11),slope_all(DI_phases==-11),'.','MarkerSize',12, 'color', Sh1_color)
    plot3(deltat_all(DI_phases==-1),deltah_all(DI_phases==-1),slope_all(DI_phases==-1),'.','MarkerSize',12, 'color', DnSt_color)
    plot3(deltat_all(DI_phases==0),deltah_all(DI_phases==0),slope_all(DI_phases==0),'.','MarkerSize',12, 'color', FlSt_color)
    title('');
    grid on
end

if(SAVE_FIG_3_movie == 1)
    figure(fig3handle)
    hold on
    
    % Movie pics filename
    moviefilename = sprintf('%s',char(FILE_NAME));
    moviefolder = sprintf('%s_DIPhaseClassification_movie',moviefilename(1:end-4));
    mkdir(moviefolder);
    moviefilename = sprintf('%s/3DplotMovie',moviefolder);
    
    % Spin Plot for 3D movie
    pixcounter = 0;
    RotDeg = 5; % Degree increments for rotation
    az = -45;
    el1 = 15;
    el2 = -10;
    el = el1;
    view([az,el]);
    axis vis3d
    
    print(fig3handle, [moviefilename '000'],'-dpng');
    pixcounter = pixcounter+1;
    for(i=1:floor(360*2/RotDeg))
        az = az + RotDeg;
        view([az,el])
        print(fig3handle, [moviefilename num2str(pixcounter,'%03.0f')], '-dpng');
        pixcounter = pixcounter+1;
        if(mod(i,floor(360/RotDeg))==0)
            RotVals = [el1:-5:el2];
            if(el == el2) RotVals = fliplr(RotVals); end
            for(j=RotVals)
                el = j;
                view([az,el])
                print(fig3handle, [moviefilename num2str(pixcounter,'%03.0f')], '-dpng');
                pixcounter = pixcounter+1;
            end
        end
    end
end

if(PLOT_FIG_3 == 1)
    figure(fig3handle)
    view([50,30]);
    xlabel('Time')
    ylabel('Height')
    zlabel('Slope')
    set(gca,'Fontsize', 16)
    title('DI Segment Phase Classification', 'Fontsize', 24)
end

% MT Length History Plot with DI-Phases Labeled for each segment
if PLOT_FIG_4 == 1;
    display('---Start Figure 4');
    fig4handle = figure;
    set(0,'units','pixels') % Resolution and screen size info needed for PNG image files
    Pixels= get(0,'screensize');
    set(0,'units','inches')
    Inches= get(0,'screensize');
    Res = Pixels/Inches;
    Papersize = [0 0 Pixels(3) 0.55*Pixels(4)];
    set(fig4handle, 'PaperUnits', 'inches', 'PaperPosition', Papersize/Res);  
    set(fig4handle, 'units','normalized','Position', [0 1 1 0.55]) 
    hold on
    
    % Limit plot window with user-defined values
    if( FIG_WINDOW_START<min(Time) | FIG_WINDOW_END>max(Time) )
        error('Figure Window values outside data time range')
    end
      
    % Collect vertices associated with times ranging the figure window size
    limited_vertices_all = vertex_all(Time(vertex_all)>=FIG_WINDOW_START & Time(vertex_all)<=FIG_WINDOW_END);
    while( Time(min(limited_vertices_all)) > FIG_WINDOW_START )
        limited_vertices_all = [vertex_all(find(vertex_all == min(limited_vertices_all)) - 1), limited_vertices_all];
    end
    while( Time(max(limited_vertices_all)) < FIG_WINDOW_END )
        limited_vertices_all = [limited_vertices_all, vertex_all(find(vertex_all == max(limited_vertices_all)) + 1)];
    end
    
    limited_indices = 1:length(vertex_all);
    limited_indices = limited_indices(vertex_all == min(limited_vertices_all)):limited_indices(vertex_all == max(limited_vertices_all));
    
    % Assign DI Phase Classifications color to each Segment
    %for(i=1:length(DI_phases(1:1000)))
    for(i=limited_indices)
        if( isnan(DI_phases(i)) )
            patch_color = Nuc_color; % Gray for Nucleation
        elseif(DI_phases(i)==12)
            patch_color = Gr2_color; % Light Green for Growth2
        elseif(DI_phases(i)==11)
            patch_color = Gr1_color; % Dark Green for Growth1
        elseif(DI_phases(i)==1)
            patch_color = UpSt_color; % Aqua for Up Stutter
        elseif(DI_phases(i)==0)
            patch_color = FlSt_color; % Blue for Flat Stutter
        elseif(DI_phases(i)==-1)
            patch_color = DnSt_color; % Purple for Down Stutter
        elseif(DI_phases(i)==-11)
            patch_color = Sh1_color; % Dark Red for Shortening 1
        elseif(DI_phases(i)==-12)
            patch_color = Sh2_color; % Light Red for Shortening 2
        elseif(DI_phases(i)==-Inf)
            patch_color = Stch_color; % for Stitches between multiple Length History data samples
        else
            error('DI_phase assignment error')
        end
        xpts = [Time(vertex_all(i)) Time(vertex_all(i+1)) Time(vertex_all(i+1)) Time(vertex_all(i))];
        ypts = [0 0 1.1*max(MT_length) 1.1*max(MT_length)];
        % Color segments according to DI label in Plot #3
        figure(fig4handle)
        patch(xpts, ypts, patch_color, 'LineStyle', 'none');
        % Plot approximation only in non-stitch segments
        if(DI_phases(i)~=-Inf)
            plot([Time(vertex_all(i)),Time(vertex_all(i+1))], [linear_approx(vertex_all(i)),linear_approx(vertex_all(i+1))], 'LineWidth', 3, 'color', [0 0 0]);
        end
    end
    
    % Complete the MT Length History with DI Phases labeled for each segment
    figure(fig4handle)
    plot(Time(min(limited_vertices_all):max(limited_vertices_all)),MT_length(min(limited_vertices_all):max(limited_vertices_all)), 'LineWidth', 0.5, 'color', [1 1 1]);
    xlabel('Time (seconds)')
    ylabel('Microtubule Length (dimers)')
    set(gca, 'Fontsize', 16);
    title('Labeled DI Phase Classes Based on MT Length History', 'Fontsize', 24)
    axis([Time(min(limited_vertices_all)),Time(max(limited_vertices_all)),0,1.1*max(MT_length)])
end % End Figure 4: MT Length History with DI Phases



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    At this point, the DI Phase Classification stage is complete    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create variables for file saving purposes

% Table of DI Phase Segment data
% Start_Time   MT-Length    Delta-t     Delta-h     Slope   DI-Phase-Label
format long g
StartTimes = round(Time(vertex_all(1:end-1)),4);
StartMTlength = round(MT_length(vertex_all(1:end-1)),4);
SegmentTime = round(X_all(:,1),4);
SegmentHeight = round(X_all(:,2),4);
SegmentSlope = round(X_all(:,3),4);
Export_Table1 = table(StartTimes,StartMTlength,SegmentTime,SegmentHeight,SegmentSlope,DI_phases);


if (SET_BENCHMARK_CLASSES == 1) % Create file of Benchmark class info 
    % Organize positive slope segment info
    if(KMEANS_NumClust_PosSlope == 1)
        PosRownames = {'PosSlope_Mean','PosSlope_StDev','C_LngGr'}';
        PosSegTime = [MU_PLUS_TRANSFORM(1); SIGMA_PLUS_TRANSFORM(1);C_pos(gr_idx,1)];
        PosSegHeight = [MU_PLUS_TRANSFORM(2); SIGMA_PLUS_TRANSFORM(2);C_pos(gr_idx,2)];
        PosSegSlope = [MU_PLUS_TRANSFORM(3); SIGMA_PLUS_TRANSFORM(3);C_pos(gr_idx,3)];
    elseif(KMEANS_NumClust_PosSlope == 2)
        if(KMEANS_Pos2_Option == 'A')
            PosRownames = {'PosSlope_Mean','PosSlope_StDev','C_LngGr','C_BrfGr'}';
            PosSegTime = [MU_PLUS_TRANSFORM(1); SIGMA_PLUS_TRANSFORM(1);C_pos(gr2_idx,1);C_pos(gr1_idx,1)];
            PosSegHeight = [MU_PLUS_TRANSFORM(2); SIGMA_PLUS_TRANSFORM(2);C_pos(gr2_idx,2);C_pos(gr1_idx,2)];
            PosSegSlope = [MU_PLUS_TRANSFORM(3); SIGMA_PLUS_TRANSFORM(3);C_pos(gr2_idx,3);C_pos(gr1_idx,3)];
        elseif(KMEANS_Pos2_Option == 'B')
            PosRownames = {'PosSlope_Mean','PosSlope_StDev','C_LngGr','C_UpSt'}';
            PosSegTime = [MU_PLUS_TRANSFORM(1); SIGMA_PLUS_TRANSFORM(1);C_pos(gr_idx,1);C_pos(upst_idx,1)];
            PosSegHeight = [MU_PLUS_TRANSFORM(2); SIGMA_PLUS_TRANSFORM(2);C_pos(gr_idx,2);C_pos(upst_idx,2)];
            PosSegSlope = [MU_PLUS_TRANSFORM(3); SIGMA_PLUS_TRANSFORM(3);C_pos(gr_idx,3);C_pos(upst_idx,3)];
        end        
    elseif(KMEANS_NumClust_PosSlope == 3)
        PosRownames = {'PosSlope_Mean','PosSlope_StDev','C_LngGr','C_BrfGr','C_UpSt'}';
        PosSegTime = [MU_PLUS_TRANSFORM(1); SIGMA_PLUS_TRANSFORM(1);C_pos(gr2_idx,1);C_pos(gr1_idx,1);C_pos(upst_idx,1)];
        PosSegHeight = [MU_PLUS_TRANSFORM(2); SIGMA_PLUS_TRANSFORM(2);C_pos(gr2_idx,2);C_pos(gr1_idx,2);C_pos(upst_idx,2)];
        PosSegSlope = [MU_PLUS_TRANSFORM(3); SIGMA_PLUS_TRANSFORM(3);C_pos(gr2_idx,3);C_pos(gr1_idx,3);C_pos(upst_idx,3)];
    end
    
    % Organize negative slope segment info
    if(KMEANS_NumClust_NegSlope == 1)
        NegRownames = {'NegSlope_Mean','NegSlope_StDev','C_LngSh'}';
        NegSegTime = [MU_MINUS_TRANSFORM(1); SIGMA_MINUS_TRANSFORM(1);C_neg(sh_idx,1)];
        NegPosSegHeight = [MU_MINUS_TRANSFORM(2); SIGMA_MINUS_TRANSFORM(2);C_neg(sh_idx,2)];
        NegPosSegSlope = [MU_MINUS_TRANSFORM(3); SIGMA_MINUS_TRANSFORM(3);C_neg(sh_idx,3)];
    elseif(KMEANS_NumClust_NegSlope == 2)
        if(KMEANS_Neg2_Option == 'A')
            NegRownames = {'NegSlope_Mean','NegSlope_StDev','C_LngSh','C_BrfSh'}';
            NegSegTime = [MU_MINUS_TRANSFORM(1); SIGMA_MINUS_TRANSFORM(1);C_neg(sh2_idx,1);C_neg(sh1_idx,1);C_neg(dnst_idx,1)];
            NegPosSegHeight = [MU_MINUS_TRANSFORM(2); SIGMA_MINUS_TRANSFORM(2);C_neg(sh2_idx,2);C_neg(sh1_idx,2);C_neg(dnst_idx,2)];
            NegPosSegSlope = [MU_MINUS_TRANSFORM(3); SIGMA_MINUS_TRANSFORM(3);C_neg(sh2_idx,3);C_neg(sh1_idx,3);C_neg(dnst_idx,3)];
        elseif(KMEANS_Neg2_Option == 'B')
            NegRownames = {'NegSlope_Mean','NegSlope_StDev','C_LngSh','C_DnSt'}';
            NegSegTime = [MU_MINUS_TRANSFORM(1); SIGMA_MINUS_TRANSFORM(1);C_neg(sh_idx,1);C_neg(dnst_idx,1)];
            NegPosSegHeight = [MU_MINUS_TRANSFORM(2); SIGMA_MINUS_TRANSFORM(2);C_neg(sh_idx,2);C_neg(dnst_idx,2)];
            NegPosSegSlope = [MU_MINUS_TRANSFORM(3); SIGMA_MINUS_TRANSFORM(3);C_neg(sh_idx,3);C_neg(dnst_idx,3)];
        end        
    elseif(KMEANS_NumClust_NegSlope == 3)
        NegRownames = {'NegSlope_Mean','NegSlope_StDev','C_LngSh','C_BrfSh','C_DnSt'}';
        NegSegTime = [MU_MINUS_TRANSFORM(1); SIGMA_MINUS_TRANSFORM(1);C_neg(sh2_idx,1);C_neg(sh1_idx,1);C_neg(dnst_idx,1)];
        NegPosSegHeight = [MU_MINUS_TRANSFORM(2); SIGMA_MINUS_TRANSFORM(2);C_neg(sh2_idx,2);C_neg(sh1_idx,2);C_neg(dnst_idx,2)];
        NegPosSegSlope = [MU_MINUS_TRANSFORM(3); SIGMA_MINUS_TRANSFORM(3);C_neg(sh2_idx,3);C_neg(sh1_idx,3);C_neg(dnst_idx,3)];
    end
    
    % Combine positive and negative info into Benchmark Classes table
    Rownames = [PosRownames;NegRownames];
    SegTime = [PosSegTime;NegSegTime];
    SegHeight = [PosSegHeight;NegPosSegHeight];
    SegSlope = [PosSegSlope;NegPosSegSlope];
    Benchmark_Classes_Info_table = table(Rownames,SegTime,SegHeight,SegSlope);
end

display('--End DI Segment Phase Classification');
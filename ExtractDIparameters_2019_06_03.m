%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for Extracting DI Parameters from Phase Identified in Input Data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%% Organize Number of Segments Data by DI Phases
NumSegs_in_Gr2 = sum(DI_phases==12);
NumSegs_in_Gr1 = sum(DI_phases==11);
NumSegs_in_UpSt = sum(DI_phases==1);
NumSegs_in_FlSt = sum(DI_phases==0);
NumSegs_in_DnSt = sum(DI_phases==-1);
NumSegs_in_Sh1 = sum(DI_phases==-11);
NumSegs_in_Sh2 = sum(DI_phases==-12);
NumSegs_in_Nuc = sum(isnan(DI_phases));
TotNumSegs = [NumSegs_in_Gr2, NumSegs_in_Gr1, ...
    NumSegs_in_UpSt, NumSegs_in_FlSt, NumSegs_in_DnSt, ...
    NumSegs_in_Sh1, NumSegs_in_Sh2, NumSegs_in_Nuc];

%%% Organize Time Duration Data by DI Phases
% Total segment times for each phase:
TotTime_in_Gr2 = sum(deltat_all(DI_phases==12));
TotTime_in_Gr1 = sum(deltat_all(DI_phases==11));
TotTime_in_UpSt = sum(deltat_all(DI_phases==1));
TotTime_in_FlSt = sum(deltat_all(DI_phases==0));
TotTime_in_DnSt = sum(deltat_all(DI_phases==-1));
TotTime_in_Sh1 = sum(deltat_all(DI_phases==-11));
TotTime_in_Sh2 = sum(deltat_all(DI_phases==-12));
TotTime_in_Nuc = sum(deltat_all(isnan(DI_phases)));
TotSegTimes_in_minutes = round([TotTime_in_Gr2, TotTime_in_Gr1, ...
    TotTime_in_UpSt, TotTime_in_FlSt, TotTime_in_DnSt, ...
    TotTime_in_Sh1, TotTime_in_Sh2, TotTime_in_Nuc]/60,4);
% Percentage of time spent in each phase:
SegTimes_Perc = TotSegTimes_in_minutes/sum(TotSegTimes_in_minutes);
PercTime_in_Gr2 = SegTimes_Perc(1);
PercTime_in_Gr1 = SegTimes_Perc(2);
PercTime_in_UpSt = SegTimes_Perc(3);
PercTime_in_FlSt = SegTimes_Perc(4);
PercTime_in_DnSt = SegTimes_Perc(5);
PercTime_in_Sh1 = SegTimes_Perc(6);
PercTime_in_Sh2 = SegTimes_Perc(7);
PercTime_in_Nuc = SegTimes_Perc(8);

%%% Organize Height Change Data by DI Phases
% Total segment height changes for each phase:
TotHeight_in_Gr2 = sum(deltah_all(DI_phases==12));
TotHeight_in_Gr1 = sum(deltah_all(DI_phases==11));
TotHeight_in_UpSt = sum(deltah_all(DI_phases==1));
TotHeight_in_FlSt = sum(deltah_all(DI_phases==0));
TotHeight_in_DnSt = sum(deltah_all(DI_phases==-1));
TotHeight_in_Sh1 = sum(deltah_all(DI_phases==-11));
TotHeight_in_Sh2 = sum(deltah_all(DI_phases==-12));
TotHeight_in_Nuc = sum(deltah_all(isnan(DI_phases)));
%TotHeight_all = sum(abs(deltah_all));
SegHeights = [TotHeight_in_Gr2, TotHeight_in_Gr1, ...
    TotHeight_in_UpSt, TotHeight_in_FlSt, TotHeight_in_DnSt, ...
    TotHeight_in_Sh1, TotHeight_in_Sh2, TotHeight_in_Nuc];
TotHeight_all = sum(abs(SegHeights));
% Percentage of height change for each phase:
SegHeights_Perc = abs(SegHeights)/TotHeight_all;
PercHeight_in_Gr2 = SegHeights_Perc(1);
PercHeight_in_Gr1 = SegHeights_Perc(2);
PercHeight_in_UpSt = SegHeights_Perc(3);
PercHeight_in_FlSt = SegHeights_Perc(4);
PercHeight_in_DnSt = SegHeights_Perc(5);
PercHeight_in_Sh1 = SegHeights_Perc(6);
PercHeight_in_Sh2 = SegHeights_Perc(7);
PercHeight_in_Nuc = SegHeights_Perc(8);


%%% Organize Slope Data by DI Phases
% Median segment slopes for each phase
MedSlope_in_Gr2 = median(slope_all(DI_phases==12));
MedSlope_in_Gr1 = median(slope_all(DI_phases==11));
MedSlope_in_UpSt = median(slope_all(DI_phases==1));
MedSlope_in_FlSt = median(slope_all(DI_phases==0));
MedSlope_in_DnSt = median(slope_all(DI_phases==-1));
MedSlope_in_Sh1 = median(slope_all(DI_phases==-11));
MedSlope_in_Sh2 = median(slope_all(DI_phases==-12));
MedSlope_in_Nuc = median(slope_all(isnan(DI_phases)));
MedSlope_all = median(slope_all);
MedSegSlopes = round([MedSlope_in_Gr2, MedSlope_in_Gr1,...
    MedSlope_in_UpSt, MedSlope_in_FlSt, MedSlope_in_DnSt,...
    MedSlope_in_Sh1, MedSlope_in_Sh2, MedSlope_in_Nuc]',4);
% Average segment slopes for each phase:
AvgSlope_in_Gr2 = mean(slope_all(DI_phases==12));
AvgSlope_in_Gr1 = mean(slope_all(DI_phases==11));
AvgSlope_in_UpSt = mean(slope_all(DI_phases==1));
AvgSlope_in_FlSt = mean(slope_all(DI_phases==0));
AvgSlope_in_DnSt = mean(slope_all(DI_phases==-1));
AvgSlope_in_Sh1 = mean(slope_all(DI_phases==-11));
AvgSlope_in_Sh2 = mean(slope_all(DI_phases==-12));
AvgSlope_in_Nuc = mean(slope_all(isnan(DI_phases)));
AvgSlope_all = mean(slope_all);
AvgSegSlopes = round([AvgSlope_in_Gr2, AvgSlope_in_Gr1,...
    AvgSlope_in_UpSt, AvgSlope_in_FlSt, AvgSlope_in_DnSt,...
    AvgSlope_in_Sh1, AvgSlope_in_Sh2, AvgSlope_in_Nuc]',4);
% St.Dev. of slopes for each phase:
StDevSlope_in_Gr2 = std(slope_all(DI_phases==12));
StDevSlope_in_Gr1 = std(slope_all(DI_phases==11));
StDevSlope_in_UpSt = std(slope_all(DI_phases==1));
StDevSlope_in_FlSt = std(slope_all(DI_phases==0));
StDevSlope_in_DnSt = std(slope_all(DI_phases==-1));
StDevSlope_in_Sh1 = std(slope_all(DI_phases==-11));
StDevSlope_in_Sh2 = std(slope_all(DI_phases==-12));
StDevSlope_in_Nuc = std(slope_all(isnan(DI_phases)));
StDevSlope_all = mean(slope_all);
StDevSlopes = round([StDevSlope_in_Gr2, StDevSlope_in_Gr1, ...
    StDevSlope_in_UpSt, StDevSlope_in_FlSt, StDevSlope_in_DnSt,...
    StDevSlope_in_Sh1, StDevSlope_in_Sh2, StDevSlope_in_Nuc]',4);
% Weighted Average Slopes for each phase:
WeightedAvgSlope_in_Gr2 = sum(slope_all(DI_phases==12) .* deltat_all(DI_phases==12)) / sum(deltat_all(DI_phases==12));
WeightedAvgSlope_in_Gr1 = sum(slope_all(DI_phases==11) .* deltat_all(DI_phases==11)) / sum(deltat_all(DI_phases==11));
WeightedAvgSlope_in_UpSt = sum(slope_all(DI_phases==1) .* deltat_all(DI_phases==1)) / sum(deltat_all(DI_phases==1));
WeightedAvgSlope_in_FlSt = sum(slope_all(DI_phases==0) .* deltat_all(DI_phases==0)) / sum(deltat_all(DI_phases==0));
WeightedAvgSlope_in_DnSt = sum(slope_all(DI_phases==-1) .* deltat_all(DI_phases==-1)) / sum(deltat_all(DI_phases==-1));
WeightedAvgSlope_in_Sh1 = sum(slope_all(DI_phases==-11) .* deltat_all(DI_phases==-11)) / sum(deltat_all(DI_phases==-11));
WeightedAvgSlope_in_Sh2 = sum(slope_all(DI_phases==-12) .* deltat_all(DI_phases==-12)) / sum(deltat_all(DI_phases==-12));
WeightedAvgSlope_in_Nuc = sum(slope_all(isnan(DI_phases)) .* deltat_all(isnan(DI_phases))) / sum(deltat_all(isnan(DI_phases)));
WeightedAvgSegSlopes = round([WeightedAvgSlope_in_Gr2, WeightedAvgSlope_in_Gr1,...
    WeightedAvgSlope_in_UpSt, WeightedAvgSlope_in_FlSt, WeightedAvgSlope_in_DnSt,...
    WeightedAvgSlope_in_Sh1, WeightedAvgSlope_in_Sh2, WeightedAvgSlope_in_Nuc]',4);
% Create a matrix consolidating all slope information
slopetable = [WeightedAvgSegSlopes,MedSegSlopes,AvgSegSlopes,StDevSlopes];


if(PLOT_FIG_5 == 1)
    % Make plots for Average DI Phase Measurements
    
    % MAKE BOXPLOTS FOR ALL 3 VARIABLES
    % Create Boxplots for each segment time duration
    % Vector of Segment time durations
    segtime_all_vector = [deltat_all(DI_phases==12),deltat_all(DI_phases==11),...
        deltat_all(DI_phases==1),deltat_all(DI_phases==0),deltat_all(DI_phases==-1),...
        deltat_all(DI_phases==-11),deltat_all(DI_phases==-12)]';
    % Vector of Segment Height changes
    segheight_all_vector = [deltah_all(DI_phases==12),deltah_all(DI_phases==11),...
        deltah_all(DI_phases==1),deltah_all(DI_phases==0),deltah_all(DI_phases==-1),...
        deltah_all(DI_phases==-11),deltah_all(DI_phases==-12)]';
    % Vector of Slopes
    slope_all_vector = [slope_all(DI_phases==12),slope_all(DI_phases==11),...
        slope_all(DI_phases==1),slope_all(DI_phases==0),slope_all(DI_phases==-1),...
        slope_all(DI_phases==-11),slope_all(DI_phases==-12)]';
    % Vector of DI-Phases
    DIphaseClass = {'LngGr','BrfGr','UpSt','FltSt','DnSt','BrfSh','LngSh','Nuc'}';
    phaselabel_values = [12,11,1,0,-1,-11,-12];
    DIphaseClass_vector=[];
    for(i=1:length(phaselabel_values))
        for(j=1:length(deltat_all(DI_phases==phaselabel_values(i))))
            DIphaseClass_vector = [DIphaseClass_vector;DIphaseClass(i)];
        end
    end
    
    % Make plots for DI Phase Measurements
    fig5handle = figure;
    set(0,'units','pixels') % Resolution and screen size info needed for PNG image files
    Pixels= get(0,'screensize');
    set(0,'units','inches')
    Inches= get(0,'screensize');
    Res = Pixels/Inches;
    Papersize = [0 0 Pixels(3) 0.52*Pixels(4)];
    set(fig5handle, 'PaperUnits', 'inches', 'PaperPosition', Papersize/Res);  
    set(fig5handle, 'units','normalized','Position', [0 1 1 0.52])
    
    % Create Figure of Tables for DI Phase #'s, & Mean and StDev for Time
    % Duration, Height Change, and Slope
    for(i=1:7)
        N_Segs(i) = sum(DI_phases==phaselabel_values(i));
        MeanTime(i) = round(mean(deltat_all(DI_phases==phaselabel_values(i))),1);
        StDevTime(i) = round(std(deltat_all(DI_phases==phaselabel_values(i))),1);
        MeanHeight(i) = round(mean(deltah_all(DI_phases==phaselabel_values(i))),1);
        StDevHeight(i) = round(std(deltah_all(DI_phases==phaselabel_values(i))),1);
        MeanSlope(i) = round(mean(slope_all(DI_phases==phaselabel_values(i))),1);
        StDevSlope(i) = round(std(slope_all(DI_phases==phaselabel_values(i))),1);
    end
    DIPhaseInfoTable = [N_Segs',MeanTime',StDevTime',MeanHeight',StDevHeight',MeanSlope',StDevSlope'];
    cnames = {'# of Segments',...
        'Time Mean','Time StDev',...
        'Height Mean','Height StDev',...
        'Slope Mean','Slope StDev'};
    rnames = {'LngGr','BrfGr','UpSt','FltSt','DnSt','BrfSh','LngSh'}';
    subplot(2,6,[2:5])
    TableSize = size(DIPhaseInfoTable);
    for(i=1:TableSize(1))
        for(j=1:TableSize(2))
            if(j>1)
                strTable(i,j) = {sprintf('%.1f',DIPhaseInfoTable(i,j))};
            else
                strTable(i,j) = {sprintf('%d',DIPhaseInfoTable(i,j))};
            end
        end
    end
    t = uitable(fig5handle,'Data',strTable,...
        'ColumnName',cnames,...
        'RowName',rnames,...
        'ColumnWidth',{83});
    pos = get(subplot(2,6,[2:5]),'position');
    set(subplot(2,6,[2:5]),'yTick',[])
    set(subplot(2,6,[2:5]),'xTick',[])
    set(t,'units','normalized')
    set(t,'position',pos)
    set(t,'FontSize',14)
    text('Position',[0.15 1.08],'string','Segment Statistics for each DI Phase','FontSize',24,'FontWeight','bold')
    % For centering data in "uitable" cells
    jscrollpane = findjobj(t);
    jTable = jscrollpane.getViewport.getView;
    cellStyle = jTable.getCellStyleAt(0,0);
    cellStyle.setHorizontalAlignment(cellStyle.CENTER);
    jTable.repaint;
    
    % Create Boxplots for 3 variables for each DI phase
    h = subplot(2,6,[7:8]);
    p = get(h, 'position');
    p(1) = p(1)/2;
    set(h, 'position',p);
    boxplot(segtime_all_vector,DIphaseClass_vector);
    set(gca,'FontSize',12,'fontweight','bold');
    title('Time Durations','FontSize',18,'FontWeight','bold')
    
    subplot(2,6,[9:10]);
    boxplot(segheight_all_vector,DIphaseClass_vector)
    set(gca,'FontSize',12,'fontweight','bold');
    title('Height Changes','FontSize',18,'FontWeight','bold')
    
    h = subplot(2,6,[11:12]);
    p = get(h, 'position');
    p(1) = p(1) + (1-(p(1)+p(3)))/2;
    set(h, 'position',p);
    boxplot(slope_all_vector,DIphaseClass_vector)
    set(gca,'FontSize',12,'fontweight','bold');
    title('Slopes','FontSize',18,'FontWeight','bold')
end


if(PLOT_FIG_6 == 1)
   
    % Make plots for Total DI Phase Measurements
    fig6handle = figure;
    set(0,'units','pixels') % Resolution and screen size info needed for PNG image files
    Pixels= get(0,'screensize');
    set(0,'units','inches')
    Inches= get(0,'screensize');
    Res = Pixels/Inches;
    Papersize = [0 0 Pixels(3) Pixels(4)];
    set(fig6handle, 'PaperUnits', 'inches', 'PaperPosition', Papersize/Res);  
    set(fig6handle, 'units','normalized','Position', [0 1 1 1])

    % Pie chart with Total # of Segments for each Phase
    subplot(2,2,1);
    blanklabels = {'','','','','','','',''};
    TotNumSegs_Modified = TotNumSegs;
    TotNumSegs_Modified(TotNumSegs==0) = 0.000001;
    SegNumpiefig = pie(TotNumSegs_Modified, blanklabels);
    % Change colors to match DI Phase Color Scheme
    SegNumpiefig(1).FaceColor = Gr2_color;
    SegNumpiefig(3).FaceColor = Gr1_color;
    SegNumpiefig(5).FaceColor = UpSt_color;
    SegNumpiefig(7).FaceColor = FlSt_color;
    SegNumpiefig(9).FaceColor = DnSt_color;
    SegNumpiefig(11).FaceColor = Sh1_color;
    SegNumpiefig(13).FaceColor = Sh2_color;
    SegNumpiefig(15).FaceColor = Nuc_color;
    % Add Title and Legend with DI Phase Labels, and Percent Time Spent
    LabelNames = {'Long Growth:       ';'Brief Growth:        ';'Up Stutter:            ';...
        'Flat Stutter:          '; 'Down Stutter:       ';'Brief Shortening:  ';...
        'Long Shortening: ';'Nucleation:           '};
    LabelFull = {};
    for(i=1:length(LabelNames))
        LabelFull(i) = strcat(LabelNames(i),[' ',num2str(TotNumSegs(i))]); % strings of values
    end
    legend(LabelFull,'Location','eastoutside','Orientation','vertical','FontSize',18);
    text('Position',[-1 1.4],'string','Number of Segments for each DI Phase','FontSize',24,'FontWeight','bold')
    
    % Pie chart with Total Time Duration for each Phase
    subplot(2,2,2);
    TotSegTimes_in_minutes_Modified = TotSegTimes_in_minutes;
    TotSegTimes_in_minutes_Modified(TotSegTimes_in_minutes==0) = 0.000001;
    timepiefig = pie(TotSegTimes_in_minutes_Modified, blanklabels);
    % Change colors to match DI Phase Color Scheme
    timepiefig(1).FaceColor = Gr2_color;
    timepiefig(3).FaceColor = Gr1_color;
    timepiefig(5).FaceColor = UpSt_color;
    timepiefig(7).FaceColor = FlSt_color;
    timepiefig(9).FaceColor = DnSt_color;
    timepiefig(11).FaceColor = Sh1_color;
    timepiefig(13).FaceColor = Sh2_color;
    timepiefig(15).FaceColor = Nuc_color;
    % Add Title and Legend with DI Phase Labels, and Percent Time Spent
    LabelNames = {'Long Growth:       ';'Brief Growth:       ';'Up Stutter:           ';...
        'Flat Stutter:          '; 'Down Stutter:      ';'Brief Shortening: ';...
        'Long Shortening: ';'Nucleation:          '};
    LabelFull = {};
    for(i=1:length(LabelNames))
        LabelFull(i) = strcat(LabelNames(i),[' ',num2str(100*round(SegTimes_Perc(i),3)),'%']); % strings and percent values
    end
    legend(LabelFull,'Location','eastoutside','Orientation','vertical','FontSize',18);
    text('Position',[-0.45 1.4],'string','Percent Time Spent in each DI Phase','FontSize',24,'FontWeight','bold')
    
    % Pie chart with Total Height Change for each Phase
    subplot(2,2,4); 
    SegHeights_Modified = abs(SegHeights);
    SegHeights_Modified(SegHeights==0) = 0.000001;
    heightpiefig = pie(SegHeights_Modified, blanklabels);
    % Change colors to match DI Phase Color Scheme
    heightpiefig(1).FaceColor = Gr2_color;
    heightpiefig(3).FaceColor = Gr1_color;
    heightpiefig(5).FaceColor = UpSt_color;
    heightpiefig(7).FaceColor = FlSt_color;
    heightpiefig(9).FaceColor = DnSt_color;
    heightpiefig(11).FaceColor = Sh1_color;
    heightpiefig(13).FaceColor = Sh2_color;
    heightpiefig(15).FaceColor = Nuc_color;
    % Add Title and Legend with DI Phase Labels, and Percent Time Spent
    LabelNames = {'Long Growth:       ';'Brief Growth:       ';'Up Stutter:           ';...
        'Flat Stutter:          '; 'Down Stutter:      ';'Brief Shortening: ';...
        'Long Shortening: ';'Nucleation:          '};
    LabelFull = {};
    for(i=1:length(LabelNames))
        LabelFull(i) = strcat(LabelNames(i),[' ',num2str(100*round(SegHeights_Perc(i),3)),'%']); % strings and percent values
    end
    legend(LabelFull,'Location','eastoutside','Orientation','vertical','FontSize',18);
    text('Position',[-1 1.4],'string','Percent Height Change for each DI Phase','FontSize',24,'FontWeight','bold')
    
    % Bar plots for different slope calculations
    subplot(2,2,3)
    cnames = {'Weighted Slope','Slope Median','Slope Mean'};
    rnames = {'Long Growth','Brief Growth',...
        'Up Stutter','Flat Stutter','Down Stutter',...
        'Brief Shortening','Long Shortening','Nucleation'}';
    SlopeBarfig = bar(slopetable(1:(end-1),1:(end-1)));
    set(SlopeBarfig(1), 'FaceColor',[0.6 0.6 0])
    set(SlopeBarfig(2), 'FaceColor',[0.8 0.8 0.5])
    set(SlopeBarfig(3), 'FaceColor',[1 1 0.2])
    legend(cnames,'Location','southwest');
    xposC = 1:length(phaselabel_values);
    xposL = xposC - 0.25;
    xposR = xposC + 0.25;
    max_y_axis = 1.2*max(max(slopetable(:,1:3)));
    min_y_axis = 1.1*min(min(slopetable(:,1:3)));
    
    for(i=1:length(phaselabel_values))
        % Background color per DI Phase
        patch([xposC(i)-0.5,xposC(i)+0.5,xposC(i)+0.5,xposC(i)-0.5],[min_y_axis,min_y_axis,max_y_axis,max_y_axis],[Phase_Colors(i,1),Phase_Colors(i,2),Phase_Colors(i,3)], ...
            'HandleVisibility','off');
        % Plot Weighted Avg. Slopes
        WtdAvgSlope = sum(slope_all(DI_phases==phaselabel_values(i)) .* deltat_all(DI_phases==phaselabel_values(i))) / sum(deltat_all(DI_phases==phaselabel_values(i)));
        patch([xposL(i)-0.1,xposL(i)+0.1,xposL(i)+0.1,xposL(i)-0.1],...
            [0,0,WtdAvgSlope,WtdAvgSlope],...
            [0.6 0.6 0], 'HandleVisibility','off');
        % Plot Slope Medians
        patch([xposC(i)-0.1,xposC(i)+0.1,xposC(i)+0.1,xposC(i)-0.1],...
            [0,0,median(slope_all(DI_phases==phaselabel_values(i))),median(slope_all(DI_phases==phaselabel_values(i)))],...
            [0.8 0.8 0.5], 'HandleVisibility','off');
        % Plot Slope Means
        patch([xposR(i)-0.1,xposR(i)+0.1,xposR(i)+0.1,xposR(i)-0.1],...
            [0,0,mean(slope_all(DI_phases==phaselabel_values(i))),mean(slope_all(DI_phases==phaselabel_values(i)))],...
            [1 1 0.2], 'HandleVisibility','off');
    end
    %add x-axis line
    hold on
    plot([min(xposC)-0.5,max(xposC)+0.5],[0,0],'k','LineWidth',1.5, 'HandleVisibility','off')
    axis([min(xposC)-0.5,max(xposC)+0.5,min_y_axis,max_y_axis])
    set(gca,'XTickLabel',DIphaseClass,'FontSize',12,'fontweight','bold');
    text('Position',[0 (0.1*(max_y_axis-min_y_axis)+max_y_axis)],'string','Segment Slope Averages for each DI Phase','FontSize',24,'FontWeight','bold')
end


% Bundle Growth, Stutter, and Shortening Phases in new variable
% Bundled Labels:
% Growth = 10;
% Stutter = 0;
% Shortening = -10;
phase_count = 1;
bundle_phase_count = 0;
DI_phases_bundled = 0;
Nuc_Start_Index = [];
Nuc_End_Index = [];
Nuc_Start_Time = [];
Nuc_End_Time = [];
Nuc_Start_Height = [];
Nuc_End_Height = [];
Growth_Start_Index = [];
Growth_End_Index = [];
Growth_Start_Time = [];
Growth_End_Time = [];
Growth_Start_Height = [];
Growth_End_Height = [];
Shortening_Start_Index = [];
Shortening_End_Index = [];
Shortening_Start_Time = [];
Shortening_End_Time = [];
Shortening_Start_Height = [];
Shortening_End_Height = [];
Stutter_Start_Index = [];
Stutter_End_Index = [];
Stutter_Start_Time = [];
Stutter_End_Time = [];
Stutter_Start_Height = [];
Stutter_End_Height = [];

while(phase_count<length(DI_phases))
    bundle_phase_count = bundle_phase_count + 1;
    Bundle_Phase_Start_Index(bundle_phase_count) = vertex_all(phase_count);
    Bundle_Phase_Start_Time(bundle_phase_count) = Time(vertex_all(phase_count));
    Bundle_Phase_Start_Height(bundle_phase_count) = MT_length(vertex_all(phase_count));
    if(isnan(DI_phases(phase_count))) % Nucleation
        DI_phases_bundled(bundle_phase_count) = NaN;
        Nuc_Start_Index = [Nuc_Start_Index,vertex_all(phase_count)];
        Nuc_Start_Time = [Nuc_Start_Time,Time(vertex_all(phase_count))];
        Nuc_Start_Height = [Nuc_Start_Height,MT_length(vertex_all(phase_count))];
        while(isnan(DI_phases(phase_count+1)) && (phase_count+1)<length(DI_phases))
            phase_count = phase_count + 1;
        end
        Nuc_End_Index = [Nuc_End_Index,vertex_all(phase_count+1)];
        Nuc_End_Time = [Nuc_End_Time,Time(vertex_all(phase_count+1))];
        Nuc_End_Height = [Nuc_End_Height,MT_length(vertex_all(phase_count+1))];        
    elseif(DI_phases(phase_count) > 10) % Growth
        DI_phases_bundled(bundle_phase_count) = 10;
        Growth_Start_Index = [Growth_Start_Index,vertex_all(phase_count)];
        Growth_Start_Time = [Growth_Start_Time,Time(vertex_all(phase_count))];
        Growth_Start_Height = [Growth_Start_Height,MT_length(vertex_all(phase_count))];        
        while(DI_phases(phase_count+1) > 10 && (phase_count+1)<length(DI_phases))
            phase_count = phase_count + 1;
        end
        Growth_End_Index = [Growth_End_Index,vertex_all(phase_count+1)];
        Growth_End_Time = [Growth_End_Time,Time(vertex_all(phase_count+1))];
        Growth_End_Height = [Growth_End_Height,MT_length(vertex_all(phase_count+1))];   
    elseif(DI_phases(phase_count) < -10 && DI_phases(phase_count) > -inf) % Shortening
        DI_phases_bundled(bundle_phase_count) = -10;
        Shortening_Start_Index = [Shortening_Start_Index,vertex_all(phase_count)];
        Shortening_Start_Time = [Shortening_Start_Time,Time(vertex_all(phase_count))];
        Shortening_Start_Height = [Shortening_Start_Height,MT_length(vertex_all(phase_count))];
        while(DI_phases(phase_count+1) < -10 && DI_phases(phase_count+1) > -inf && (phase_count+1)<length(DI_phases))
            phase_count = phase_count + 1;
        end
        Shortening_End_Index = [Shortening_End_Index,vertex_all(phase_count+1)];
        Shortening_End_Time = [Shortening_End_Time,Time(vertex_all(phase_count+1))];
        Shortening_End_Height = [Shortening_End_Height,MT_length(vertex_all(phase_count+1))];   
    elseif(DI_phases(phase_count) > -10 && DI_phases(phase_count) < 10) % Stutter
        DI_phases_bundled(bundle_phase_count) = 0; 
        Stutter_Start_Index = [Stutter_Start_Index,vertex_all(phase_count)];
        Stutter_Start_Time = [Stutter_Start_Time,Time(vertex_all(phase_count))];
        Stutter_Start_Height = [Stutter_Start_Height,MT_length(vertex_all(phase_count))];
        while(DI_phases(phase_count+1) > -10 && DI_phases(phase_count+1) < 10 && (phase_count+1)<length(DI_phases))
            phase_count = phase_count + 1;
        end
        Stutter_End_Index = [Stutter_End_Index,vertex_all(phase_count+1)];
        Stutter_End_Time = [Stutter_End_Time,Time(vertex_all(phase_count+1))];
        Stutter_End_Height = [Stutter_End_Height,MT_length(vertex_all(phase_count+1))];
    elseif(DI_phases(phase_count) == -Inf) % Stitch between multiple length history data samples
        DI_phases_bundled(bundle_phase_count) = -Inf;
        while(DI_phases(phase_count+1) == -Inf && (phase_count+1)<length(DI_phases))
            phase_count = phase_count + 1;
        end
    end
    phase_count = phase_count+1;
end % End While: Bundling phases (except for last phase segment)
%Check if last DI Phase Segment is the same as the last bundled phase
%IF yes, include last segment into last bundled phase
%ELSE, add a new last bundled phase
if(isnan(DI_phases(phase_count))) % Nucleation
    if(DI_phases_bundled(bundle_phase_count) == NaN)
        Nuc_End_Index(end) = vertex_all(phase_count+1);
        Nuc_End_Time(end) = Time(vertex_all(phase_count+1));
        Nuc_End_Height(end) = MT_length(vertex_all(phase_count+1));
    else
        bundle_phase_count = bundle_phase_count + 1;
        DI_phases_bundled(bundle_phase_count) = NaN;
        Bundle_Phase_Start_Index(bundle_phase_count) = vertex_all(phase_count);
        Bundle_Phase_Start_Time(bundle_phase_count) = Time(vertex_all(phase_count));
        Bundle_Phase_Start_Height(bundle_phase_count) = MT_length(vertex_all(phase_count));
        Nuc_Start_Index = [Nuc_Start_Index,vertex_all(phase_count)];
        Nuc_Start_Time = [Nuc_Start_Time,Time(vertex_all(phase_count))];
        Nuc_Start_Height = [Nuc_Start_Height,MT_length(vertex_all(phase_count))];
        Nuc_End_Index = [Nuc_End_Index,vertex_all(phase_count+1)];
        Nuc_End_Time = [Nuc_End_Time,Time(vertex_all(phase_count+1))];
        Nuc_End_Height = [Nuc_End_Height,MT_length(vertex_all(phase_count+1))];
    end
elseif(DI_phases(phase_count) > 10) % Growth
    if(DI_phases_bundled(bundle_phase_count) == 10)
        Growth_End_Index(end) = vertex_all(phase_count+1);
        Growth_End_Time(end) = Time(vertex_all(phase_count+1));
        Growth_End_Height(end) = MT_length(vertex_all(phase_count+1));
    else
        bundle_phase_count = bundle_phase_count + 1;
        DI_phases_bundled(bundle_phase_count) = 10;
        Bundle_Phase_Start_Index(bundle_phase_count) = vertex_all(phase_count);
        Bundle_Phase_Start_Time(bundle_phase_count) = Time(vertex_all(phase_count));
        Bundle_Phase_Start_Height(bundle_phase_count) = MT_length(vertex_all(phase_count));
        Growth_Start_Index = [Growth_Start_Index,vertex_all(phase_count)];
        Growth_Start_Time = [Growth_Start_Time,Time(vertex_all(phase_count))];
        Growth_Start_Height = [Growth_Start_Height,MT_length(vertex_all(phase_count))];        
        Growth_End_Index = [Growth_End_Index,vertex_all(phase_count+1)];
        Growth_End_Time = [Growth_End_Time,Time(vertex_all(phase_count+1))];
        Growth_End_Height = [Growth_End_Height,MT_length(vertex_all(phase_count+1))];   
    end
elseif(DI_phases(phase_count) < -10 && DI_phases(phase_count) > -inf) % Shortening
    if(DI_phases_bundled(bundle_phase_count) == -10)
        Shortening_End_Index(end) = vertex_all(phase_count+1);
        Shortening_End_Time(end) = Time(vertex_all(phase_count+1));
        Shortening_End_Height(end) = MT_length(vertex_all(phase_count+1));
    else
        bundle_phase_count = bundle_phase_count + 1;
        DI_phases_bundled(bundle_phase_count) = -10;
        Bundle_Phase_Start_Index(bundle_phase_count) = vertex_all(phase_count);
        Bundle_Phase_Start_Time(bundle_phase_count) = Time(vertex_all(phase_count));
        Bundle_Phase_Start_Height(bundle_phase_count) = MT_length(vertex_all(phase_count));
        Shortening_Start_Index = [Shortening_Start_Index,vertex_all(phase_count)];
        Shortening_Start_Time = [Shortening_Start_Time,Time(vertex_all(phase_count))];
        Shortening_Start_Height = [Shortening_Start_Height,MT_length(vertex_all(phase_count))];
        Shortening_End_Index = [Shortening_End_Index,vertex_all(phase_count+1)];
        Shortening_End_Time = [Shortening_End_Time,Time(vertex_all(phase_count+1))];
        Shortening_End_Height = [Shortening_End_Height,MT_length(vertex_all(phase_count+1))];   
    end
elseif(DI_phases(phase_count) > -10 && DI_phases(phase_count) < 10) % Stutter
    if(DI_phases_bundled(bundle_phase_count) == 0)
        Stutter_End_Index(end) = vertex_all(phase_count+1);
        Stutter_End_Time(end) = Time(vertex_all(phase_count+1));
        Stutter_End_Height(end) = MT_length(vertex_all(phase_count+1));
    else
        bundle_phase_count = bundle_phase_count + 1;
        DI_phases_bundled(bundle_phase_count) = 0; 
        Bundle_Phase_Start_Index(bundle_phase_count) = vertex_all(phase_count);
        Bundle_Phase_Start_Time(bundle_phase_count) = Time(vertex_all(phase_count));
        Bundle_Phase_Start_Height(bundle_phase_count) = MT_length(vertex_all(phase_count));
        Stutter_Start_Index = [Stutter_Start_Index,vertex_all(phase_count)];
        Stutter_Start_Time = [Stutter_Start_Time,Time(vertex_all(phase_count))];
        Stutter_Start_Height = [Stutter_Start_Height,MT_length(vertex_all(phase_count))];
        Stutter_End_Index = [Stutter_End_Index,vertex_all(phase_count+1)];
        Stutter_End_Time = [Stutter_End_Time,Time(vertex_all(phase_count+1))];
        Stutter_End_Height = [Stutter_End_Height,MT_length(vertex_all(phase_count+1))];
    end
elseif(DI_phases(phase_count) == -Inf) % Stitch between multiple length history data samples
    if(DI_phases_bundled(bundle_phase_count) > -Inf)
        bundle_phase_count = bundle_phase_count + 1;
        DI_phases_bundled(bundle_phase_count) = -Inf;
        Bundle_Phase_Start_Index(bundle_phase_count) = vertex_all(phase_count);
        Bundle_Phase_Start_Time(bundle_phase_count) = Time(vertex_all(phase_count));
        Bundle_Phase_Start_Height(bundle_phase_count) = MT_length(vertex_all(phase_count));
    end
end %DONE BUNDLING PHASES


Bundled_Nuc_Table = table(Nuc_Start_Index',Nuc_End_Index',Nuc_Start_Time',Nuc_End_Time',Nuc_Start_Height',Nuc_End_Height');
Bundled_Nuc_Table.Properties.VariableNames = {'Nuc_Start_Index','Nuc_End_Index', ...
    'Nuc_Start_Time', 'Nuc_End_Time', ...
    'Nuc_Start_Height', 'Nuc_End_Height'};

Bundled_Growth_Table = table(Growth_Start_Index',Growth_End_Index',Growth_Start_Time',Growth_End_Time',Growth_Start_Height',Growth_End_Height');
Bundled_Growth_Table.Properties.VariableNames = {'Growth_Start_Index','Growth_End_Index', ...
'Growth_Start_Time','Growth_End_Time', ...
'Growth_Start_Height','Growth_End_Height'};

Bundled_Shortening_Table = table(Shortening_Start_Index',Shortening_End_Index',Shortening_Start_Time',Shortening_End_Time',Shortening_Start_Height',Shortening_End_Height');
Bundled_Shortening_Table.Properties.VariableNames = {'Shortening_Start_Index','Shortening_End_Index',...
'Shortening_Start_Time','Shortening_End_Time',...
'Shortening_Start_Height','Shortening_End_Height'};

Bundled_Stutter_Table = table(Stutter_Start_Index',Stutter_End_Index',Stutter_Start_Time',Stutter_End_Time',Stutter_Start_Height',Stutter_End_Height');
Bundled_Stutter_Table.Properties.VariableNames = {'Stutter_Start_Index','Stutter_End_Index',...
'Stutter_Start_Time','Stutter_End_Time',...
'Stutter_Start_Height','Stutter_End_Height'};



Bundled_Nuc_Array = [Nuc_Start_Index;Nuc_End_Index;Nuc_Start_Time;Nuc_End_Time;Nuc_Start_Height;Nuc_End_Height]';
Bundled_Growth_Array = [Growth_Start_Index;Growth_End_Index;Growth_Start_Time;Growth_End_Time;Growth_Start_Height;Growth_End_Height]';
Bundled_Shortening_Array = [Shortening_Start_Index;Shortening_End_Index;Shortening_Start_Time;Shortening_End_Time;Shortening_Start_Height;Shortening_End_Height]';
Bundled_Stutter_Array = [Stutter_Start_Index;Stutter_End_Index;Stutter_Start_Time;Stutter_End_Time;Stutter_Start_Height;Stutter_End_Height]';



% Count the different phase changes:
Number_of_Abrupt_Catastrophe = 0;
Abrupt_Catastrophe_Start_Index = [];
Abrupt_Catastrophe_Start_Time = [];
Abrupt_Catastrophe_Start_Height = [];

Number_of_Abrupt_Rescue = 0;
Abrupt_Rescue_Start_Index = [];
Abrupt_Rescue_Start_Time = [];
Abrupt_Rescue_Start_Height = [];

Number_of_Transitional_Catastrophe = 0;
Transitional_Catastrophe_Stutter_Start_Index = [];
Transitional_Catastrophe_Shortening_Start_Index = [];
Transitional_Catastrophe_Stutter_Start_Time = [];
Transitional_Catastrophe_Shortening_Start_Time = [];
Transitional_Catastrophe_Stutter_Start_Height = [];
Transitional_Catastrophe_Shortening_Start_Height = [];

Number_of_Transitional_Rescue = 0;
Transitional_Rescue_Stutter_Start_Index = [];
Transitional_Rescue_Growth_Start_Index = [];
Transitional_Rescue_Stutter_Start_Time = [];
Transitional_Rescue_Growth_Start_Time = [];
Transitional_Rescue_Stutter_Start_Height = [];
Transitional_Rescue_Growth_Start_Height = [];

Number_of_Interupted_Growth = 0;
Interupted_Growth_Stutter_Start_Index = [];
Interupted_Growth_Growth_Start_Index = [];
Interupted_Growth_Stutter_Start_Time = [];
Interupted_Growth_Growth_Start_Time = [];
Interupted_Growth_Stutter_Start_Height = [];
Interupted_Growth_Growth_Start_Height = [];

Number_of_Interupted_Shortening = 0;
Interupted_Shortening_Stutter_Start_Index = [];
Interupted_Shortening_Shortening_Start_Index = [];
Interupted_Shortening_Stutter_Start_Time = [];
Interupted_Shortening_Shortening_Start_Time = [];
Interupted_Shortening_Stutter_Start_Height = [];
Interupted_Shortening_Shortening_Start_Height = [];

    
for(i=1:(length(DI_phases_bundled)-2))
    SegI = DI_phases_bundled(i);
    SegII = DI_phases_bundled(i+1);
    SegIII = DI_phases_bundled(i+2);
    if( SegI==10 && SegII==-10 ) % Abrupt Catastrophe: Gr-Sh
        Number_of_Abrupt_Catastrophe = Number_of_Abrupt_Catastrophe + 1;
        Abrupt_Catastrophe_Start_Index = [Abrupt_Catastrophe_Start_Index,Bundle_Phase_Start_Index(i+1)];
        Abrupt_Catastrophe_Start_Time = [Abrupt_Catastrophe_Start_Time,Bundle_Phase_Start_Time(i+1)];
        Abrupt_Catastrophe_Start_Height = [Abrupt_Catastrophe_Start_Height,Bundle_Phase_Start_Height(i+1)];
    elseif( SegI==-10 && SegII==10 ) % Abrupt Rescue: Sh-Gr
        Number_of_Abrupt_Rescue = Number_of_Abrupt_Rescue + 1;
        Abrupt_Rescue_Start_Index = [Abrupt_Rescue_Start_Index,Bundle_Phase_Start_Index(i+1)];
        Abrupt_Rescue_Start_Time = [Abrupt_Rescue_Start_Time,Bundle_Phase_Start_Time(i+1)];
        Abrupt_Rescue_Start_Height = [Abrupt_Rescue_Start_Height,Bundle_Phase_Start_Height(i+1)];
    elseif( SegI==10 && SegII==0 && SegIII==-10 ) % Transitional Catastrophe: Gr-St-Sh
        Number_of_Transitional_Catastrophe = Number_of_Transitional_Catastrophe + 1;
        Transitional_Catastrophe_Stutter_Start_Index = [Transitional_Catastrophe_Stutter_Start_Index,Bundle_Phase_Start_Index(i+1)];
        Transitional_Catastrophe_Shortening_Start_Index = [Transitional_Catastrophe_Shortening_Start_Index,Bundle_Phase_Start_Index(i+2)];
        Transitional_Catastrophe_Stutter_Start_Time = [Transitional_Catastrophe_Stutter_Start_Time,Bundle_Phase_Start_Time(i+1)];
        Transitional_Catastrophe_Shortening_Start_Time = [Transitional_Catastrophe_Shortening_Start_Time,Bundle_Phase_Start_Time(i+2)];
        Transitional_Catastrophe_Stutter_Start_Height = [Transitional_Catastrophe_Stutter_Start_Height,Bundle_Phase_Start_Height(i+1)];
        Transitional_Catastrophe_Shortening_Start_Height = [Transitional_Catastrophe_Shortening_Start_Height,Bundle_Phase_Start_Height(i+2)];        
    elseif( SegI==-10 && SegII==0 && SegIII==10 ) % Transitional Rescue: Sh-St-Gr
        Number_of_Transitional_Rescue = Number_of_Transitional_Rescue + 1;
        Transitional_Rescue_Stutter_Start_Index = [Transitional_Rescue_Stutter_Start_Index,Bundle_Phase_Start_Index(i+1)];
        Transitional_Rescue_Growth_Start_Index = [Transitional_Rescue_Growth_Start_Index,Bundle_Phase_Start_Index(i+2)];
        Transitional_Rescue_Stutter_Start_Time = [Transitional_Rescue_Stutter_Start_Time,Bundle_Phase_Start_Time(i+1)];
        Transitional_Rescue_Growth_Start_Time = [Transitional_Rescue_Growth_Start_Time,Bundle_Phase_Start_Time(i+2)];
        Transitional_Rescue_Stutter_Start_Height = [Transitional_Rescue_Stutter_Start_Height,Bundle_Phase_Start_Height(i+1)];
        Transitional_Rescue_Growth_Start_Height = [Transitional_Rescue_Growth_Start_Height,Bundle_Phase_Start_Height(i+2)];        
    elseif( SegI==10 && SegII==0 && SegIII==10 ) % Interupted Growth: Gr-St-Gr
        Number_of_Interupted_Growth = Number_of_Interupted_Growth + 1;
        Interupted_Growth_Stutter_Start_Index = [Interupted_Growth_Stutter_Start_Index,Bundle_Phase_Start_Index(i+1)];
        Interupted_Growth_Growth_Start_Index = [Interupted_Growth_Growth_Start_Index,Bundle_Phase_Start_Index(i+2)];
        Interupted_Growth_Stutter_Start_Time = [Interupted_Growth_Stutter_Start_Time,Bundle_Phase_Start_Time(i+1)];
        Interupted_Growth_Growth_Start_Time = [Interupted_Growth_Growth_Start_Time,Bundle_Phase_Start_Time(i+2)];
        Interupted_Growth_Stutter_Start_Height = [Interupted_Growth_Stutter_Start_Height,Bundle_Phase_Start_Height(i+1)];
        Interupted_Growth_Growth_Start_Height = [Interupted_Growth_Growth_Start_Height,Bundle_Phase_Start_Height(i+2)];                
    elseif( SegI==-10 && SegII==0 && SegIII==-10 ) % Interupted Shortening: Sh-St-Sh
        Number_of_Interupted_Shortening = Number_of_Interupted_Shortening + 1;
        Interupted_Shortening_Stutter_Start_Index = [Interupted_Shortening_Stutter_Start_Index,Bundle_Phase_Start_Index(i+1)];
        Interupted_Shortening_Shortening_Start_Index = [Interupted_Shortening_Shortening_Start_Index,Bundle_Phase_Start_Index(i+2)];
        Interupted_Shortening_Stutter_Start_Time = [Interupted_Shortening_Stutter_Start_Time,Bundle_Phase_Start_Time(i+1)];
        Interupted_Shortening_Shortening_Start_Time = [Interupted_Shortening_Shortening_Start_Time,Bundle_Phase_Start_Time(i+2)];
        Interupted_Shortening_Stutter_Start_Height = [Interupted_Shortening_Stutter_Start_Height,Bundle_Phase_Start_Height(i+1)];
        Interupted_Shortening_Shortening_Start_Height = [Interupted_Shortening_Shortening_Start_Height,Bundle_Phase_Start_Height(i+2)];                
    end
end

Frequency_of_Abrupt_Catastrophe = Number_of_Abrupt_Catastrophe/(TotTime_in_Gr2+TotTime_in_Gr1);
Frequency_of_Transitional_Catastrophe = Number_of_Transitional_Catastrophe/(TotTime_in_Gr2+TotTime_in_Gr1);
Frequency_of_Abrupt_Rescue = Number_of_Abrupt_Rescue/(TotTime_in_Sh2+TotTime_in_Sh1);
Frequency_of_Transitional_Rescue = Number_of_Transitional_Rescue/(TotTime_in_Sh2+TotTime_in_Sh1);
Frequency_of_Interupted_Growth = Number_of_Interupted_Growth/(TotTime_in_Gr2+TotTime_in_Gr1);
Frequency_of_Interupted_Shortening = Number_of_Interupted_Shortening/(TotTime_in_Sh2+TotTime_in_Sh1);

DIPhaseChangeNumbers = [Number_of_Abrupt_Catastrophe, Number_of_Transitional_Catastrophe, ...
    Number_of_Abrupt_Rescue, Number_of_Transitional_Rescue,...
    Number_of_Interupted_Growth, Number_of_Interupted_Shortening]';
DIPhaseChangeFrequencies = [Frequency_of_Abrupt_Catastrophe, Frequency_of_Transitional_Catastrophe,...
    Frequency_of_Abrupt_Rescue, Frequency_of_Transitional_Rescue, ...
	Frequency_of_Interupted_Growth, Frequency_of_Interupted_Shortening]';


% Count the different Catastrophe+Rescue events:
% ACat = Abrupt Catastrophe;    TCat = Transitional Catastrophe
% ARes = Abrupt Rescue;         TRes = Transitional Rescue
% IntGr = Interupted Growth;    IntSh = Interupted Shortening
Number_of_ACat_ARes = 0; %
Number_of_TCat_ARes = 0; %
Number_of_ACat_TRes = 0; %
Number_of_TCat_TRes = 0; %
Number_of_IntGr_ACat = 0; %
Number_of_IntGr_TCat = 0; %
Number_of_IntSh_ARes = 0; %
Number_of_IntSh_TRes = 0; %
Number_of_ACat_IntSh = 0; %
Number_of_TCat_IntSh = 0; %
Number_of_ARes_IntGr = 0; %
Number_of_TRes_IntGr = 0; %
for(i=1:(length(DI_phases_bundled)-3))
    SegI = DI_phases_bundled(i);
    SegII = DI_phases_bundled(i+1);
    SegIII = DI_phases_bundled(i+2);
    SegIV = DI_phases_bundled(i+3);
    if(i<=(length(DI_phases_bundled)-4)) 
        SegV = DI_phases_bundled(i+4);
    end
    if( SegI==10 && SegII==-10 && SegIII==10 ) % ACat + ARes
        Number_of_ACat_ARes = Number_of_ACat_ARes + 1;
    elseif( SegI==10 && SegII==0 && SegIII==-10 && SegIV==10 ) % TCat + ARes
        Number_of_TCat_ARes = Number_of_TCat_ARes + 1;
    elseif( SegI==10 && SegII==-10 && SegIII==0 && SegIV==10 ) % ACat + TRes
        Number_of_ACat_TRes = Number_of_ACat_TRes + 1;
    elseif( SegI==10 && SegII==0 && SegIII==10 && SegIV==-10 ) % IntGr + ACat
        Number_of_IntGr_ACat = Number_of_IntGr_ACat + 1;
    elseif( SegI==-10 && SegII==0 && SegIII==-10 && SegIV==10 ) % IntSh + ARes
        Number_of_IntSh_ARes = Number_of_IntSh_ARes + 1;
    elseif( SegI==10 && SegII==-10 && SegIII==0 && SegIV==-10 ) % ACat + IntSh
        Number_of_ACat_IntSh = Number_of_ACat_IntSh + 1;
    elseif( SegI==10 && SegII==-10 && SegIII==0 && SegIV==-10 ) % ARes + IntGr
        Number_of_ARes_IntGr = Number_of_ARes_IntGr + 1;
    elseif( i<=(length(DI_phases_bundled)-4))
        if( SegI==10 && SegII==0 && SegIII==-10 && SegIV==0 && SegV==10 ) % TCat + TRes
            Number_of_TCat_TRes = Number_of_TCat_TRes + 1;
        elseif( SegI==10 && SegII==0 && SegIII==10 && SegIV==0 && SegV==-10 ) % IntGr + TCat
            Number_of_IntGr_TCat = Number_of_IntGr_TCat + 1;
        elseif( SegI==-10 && SegII==0 && SegIII==-10 && SegIV==0 && SegV==10 ) % IntSh + TRes
            Number_of_IntSh_TRes = Number_of_IntSh_TRes + 1;
        elseif( SegI==10 && SegII==0 && SegIII==-10 && SegIV==0 && SegV==-10 ) % TCat + IntSh
            Number_of_TCat_IntSh = Number_of_TCat_IntSh + 1;
        elseif( SegI==-10 && SegII==0 && SegIII==10 && SegIV==0 && SegV==10 ) % TRes + IntGr
            Number_of_TRes_IntGr = Number_of_TRes_IntGr + 1;
        end
    end
end
% table(Number_of_ACat_ARes,Number_of_TCat_ARes,Number_of_ACat_TRes,Number_of_TCat_TRes,...
%    Number_of_IntGr_ACat, Number_of_IntGr_TCat, Number_of_IntSh_ARes, Number_of_IntSh_TRes)

CatPlusResVals = [Number_of_ACat_ARes,Number_of_TCat_ARes,Number_of_ACat_TRes,Number_of_TCat_TRes]';
IntGrShCombosVals = [Number_of_IntGr_ACat, Number_of_IntGr_TCat, Number_of_IntSh_ARes, Number_of_IntSh_TRes,...
    Number_of_ACat_IntSh, Number_of_TCat_IntSh, Number_of_ARes_IntGr, Number_of_TRes_IntGr]';

% 
% 
% % Create Figure of Tables for DI Phase Change Results
% DIPhaseChangeTableFig = figure;
% set(DIPhaseChangeTableFig, 'units','normalized','Position', [0 0.4 0.475 0.6]) 
% DIPhaseChangeTable = [DIPhaseChangeNumbers,DIPhaseChangeFrequencies];
% subplot(3,1,1)
% cnames = {'Number','Frequency'};
% rnames = {'Abrupt Catastrophes','Transitional Catastrophes',...
%     'Abrupt Rescues','Transitional Rescues',...
%     'Interupted_Growths', 'Interupted_Shortenings'}';
%  t = uitable(DIPhaseChangeTableFig,'Data',DIPhaseChangeTable,...
%              'ColumnName',cnames,... 
%              'RowName',rnames,...
%              'ColumnWidth',{100});
% pos = get(subplot(2,1,1),'position');
% set(subplot(2,1,1),'yTick',[])
% set(subplot(2,1,1),'xTick',[])
% set(t,'units','normalized')
% set(t,'position',pos)
% text('Position',[0.12 1.05],'string','DI Phase Change Measurements','FontSize',24,'FontWeight','bold')
% 
% subplot(3,1,2)
% CatResCombosVals = [Number_of_ACat_ARes,Number_of_TCat_ARes,Number_of_ACat_TRes,Number_of_TCat_TRes,...
%     Number_of_IntGr_ACat, Number_of_IntGr_TCat, Number_of_IntSh_ARes, Number_of_IntSh_TRes]';
% rnames = {'Abrupt Cat + Abrupt Resc','Transitional Cat + Abrupt Resc',...
%     'Abrupt Cat + Transitional Resc','Transitional Cat + Transitional Resc',...
%     'Interupted Short + Abrupt Cat', 'Interupted Short + Transitional Cat',...
%     'Interupted Growth + Abrupt Resc', 'Interupted Growth + Transitional Resc',};
%  t = uitable(DIPhaseChangeTableFig,'Data',CatResCombosVals,...
%              'ColumnName',{''},... 
%              'RowName',rnames,...
%              'ColumnWidth',{100});
% pos = get(subplot(2,1,2),'position');
% set(subplot(2,1,2),'yTick',[])
% set(subplot(2,1,2),'xTick',[])
% set(t,'units','normalized')
% set(t,'position',pos)
% text('Position',[0.12 1.05],'string','Catastrophe & Rescue Combinations','FontSize',24,'FontWeight','bold')
% 
% 


if(PLOT_FIG_7 == 1)
    % Figures for DI Phase Change Measurements
    fig7handle = figure;
    set(0,'units','pixels') % Resolution and screen size info needed for PNG image files
    Pixels= get(0,'screensize');
    set(0,'units','inches')
    Inches= get(0,'screensize');
    Res = Pixels/Inches;
    Papersize = [0 0 Pixels(3) Pixels(4)];
    set(fig7handle, 'PaperUnits', 'inches', 'PaperPosition', Papersize/Res);  
    set(fig7handle, 'units','normalized','Position', [0 1 1 1])
    
    % Pie Chart: Abrupt vs. Transitional Catastrophes
    subplot(6,5,[1:2,6:7]);
    blanklabels = {'',''};
    CompareCatVals = [Number_of_Abrupt_Catastrophe,Number_of_Transitional_Catastrophe];
    CompareCatVals_Modified = CompareCatVals;
    CompareCatVals_Modified(CompareCatVals==0) = 0.000001;
    CompareCatfig = pie(CompareCatVals_Modified, blanklabels);
    % Change colors to match DI Phase Color Scheme
    colorchoice(1,:) = [0 0 0];
    colorchoice(2,:) = [1 1 1];
    for(i=1:floor(length(CompareCatfig)/2))
        CompareCatfig(2*i-1).FaceColor = colorchoice(i,:);
    end
    % Add Title and Legend with Number of Occurrences
    LabelNames = {'Abrupt           ','Transitional   '};
    LabelFull = {};
    for(i=1:length(LabelNames))
        LabelFull(i) = strcat(LabelNames(i),' ',int2str(CompareCatVals(i))); % strings and values
    end
    legend(LabelFull,'Location','eastoutside','Orientation','vertical','FontSize',18);
    text('Position',[-0.75 1.2],'string','Abrupt vs. Transitional Catastrophes','FontSize',18,'FontWeight','bold')
    
    % Abrupt vs. Transitional Rescues
    subplot(6,5,[11:12,16:17]);
    blanklabels = {'',''};
    CompareRescVals = [Number_of_Abrupt_Rescue,Number_of_Transitional_Rescue];
    CompareRescVals_Modified = CompareRescVals;
    CompareRescVals_Modified(CompareRescVals==0) = 0.000001;
    CompareRescfig = pie(CompareRescVals_Modified, blanklabels);
    % Change colors to match DI Phase Color Scheme
    colorchoice(1,:) = [0 0 0];
    colorchoice(2,:) = [1 1 1];
    for(i=1:floor(length(CompareRescfig)/2))
        CompareRescfig(2*i-1).FaceColor = colorchoice(i,:);
    end
    % Add Title and Legend with Number of Occurrences
    LabelNames = {'Abrupt           ','Transitional   '};
    LabelFull = {};
    for(i=1:length(LabelNames))
        LabelFull(i) = strcat(LabelNames(i),' ',int2str(CompareRescVals(i))); % strings and values
    end
    legend(LabelFull,'Location','eastoutside','Orientation','vertical','FontSize',18);
    text('Position',[-0.5 1.2],'string','Abrupt vs. Transitional Rescues','FontSize',18,'FontWeight','bold')
    
    
    % Create Figure of Tables for DI Phase Change Results
    DIPhaseChangeTable = [DIPhaseChangeNumbers,DIPhaseChangeFrequencies];
    subplot(6,5,[21:22,26:27]);
    cnames = {'Number','Frequency'};
    rnames = {'Abrupt Cat.','Trans. Cat.',...
        'Abrupt Rescues','Trans. Rescues',...
        'Interupted Growths', 'Interupted Short'}';
    t = uitable(fig7handle,'Data',DIPhaseChangeTable,...
        'ColumnName',cnames,...
        'RowName',rnames,...
        'ColumnWidth',{90});
    pos = get(subplot(6,5,[21:22,26:27]),'position');
    set(subplot(6,5,[21:22,26:27]),'yTick',[])
    set(subplot(6,5,[21:22,26:27]),'xTick',[])
    set(t,'units','normalized')
    set(t,'position',pos)
    set(t,'FontSize',18)
    text('Position',[0.15 1.075],'string','DI Phase Change Measurements','FontSize',18,'FontWeight','bold')
    
    
    % DI Change Patterns with Stutters
    subplot(6,5,[3:5,8:10]);
    blanklabels = {'','','',''};
    CompareStutterVals = [Number_of_Transitional_Catastrophe,Number_of_Transitional_Rescue,Number_of_Interupted_Growth,Number_of_Interupted_Shortening];
    CompareStutterVals_Modified = CompareStutterVals;
    CompareStutterVals_Modified(CompareStutterVals==0) = 0.000001;
    CompareStutterfig = pie(CompareStutterVals_Modified, blanklabels);
    % Change colors to match DI Phase Color Scheme
    CompareStutterfig(1).FaceColor = [0 0 0];
    CompareStutterfig(3).FaceColor = [1 1 0];
    CompareStutterfig(5).FaceColor = [0.5 0.5 0.5];
    CompareStutterfig(7).FaceColor = [1 1 1];
    % Add Title and Legend with Number of Occurrences
    LabelNames = {'Trans. Catastrophe     ','Trans. Rescue             ',...
        'Interupted Growth       ','Interupted Shortening '};
    LabelFull = {};
    for(i=1:length(LabelNames))
        LabelFull(i) = strcat(LabelNames(i),' ',int2str(CompareStutterVals(i))); % strings and values
    end
    legend(LabelFull,'Location','eastoutside','Orientation','vertical','FontSize',18);
    text('Position',[0.7 1.1],'string','DI Changes with Stutters','FontSize',18,'FontWeight','bold')
   
    
    % All Combinations of DI Change Patterns
    subplot(6,5,[13:15,18:20,23:25,28:30]);
    blanklabels = {'','','','','','','','','','','',''};
    AllComboVals = [CatPlusResVals;IntGrShCombosVals];
    AllComboVals(AllComboVals==0) = 0.000001;
    IntGrShCombosfig = pie(AllComboVals, blanklabels);
    % Change colors to match DI Phase Color Scheme
    IntGrShCombosfig(1).FaceColor = [0 0 0];
    IntGrShCombosfig(3).FaceColor = [0.25 0.25 0.25];
    IntGrShCombosfig(5).FaceColor = [0.5 0.5 0.5];
    IntGrShCombosfig(7).FaceColor = [0.6 0.6 0.25];
    IntGrShCombosfig(9).FaceColor = [0.7 0.7 0.1];
    IntGrShCombosfig(11).FaceColor = [0.8 0.8 0.1];
    IntGrShCombosfig(13).FaceColor = [1 1 0];
    IntGrShCombosfig(15).FaceColor = [1 1 0.3];
    IntGrShCombosfig(17).FaceColor = [0.9 0.9 0.4];
    IntGrShCombosfig(19).FaceColor = [0.85 0.85 0.6];
    IntGrShCombosfig(21).FaceColor = [0.8 0.8 0.8];
    IntGrShCombosfig(23).FaceColor = [1 1 1];
    % Add Title and Legend with Number of Occurrences
    LabelNames = {'A.Cat + A.Resc   ','T.Cat + A.Resc   ','A.Cat + T.Resc   ','T.Cat + T.Resc   ',...
        'Int.Gr + A.Cat     ','Int.Gr + T.Cat     ','Int.Sh + A.Resc  ','Int.Sh + T.Resc  ',...
        'A.Cat + Int.Sh    ','T.Cat + Int.Sh    ','A.Res + Int.Gr    ','T.Res + Int.Gr    '};
    LabelFull = {};
    for(i=1:length(LabelNames))
        LabelFull(i) = strcat(LabelNames(i),' ',num2str(round(AllComboVals(i),1))); % strings and values
    end
    legend(LabelFull,'Location','eastoutside','Orientation','vertical','FontSize',18);
    text('Position',[-0.3 1.1],'string','Pairwise Combinations of DI Changes','FontSize',18,'FontWeight','bold')
end

Export_Table2 = table( {'',             [];
    'FILE_NAME',                        FILE_NAME;
    'FIRST_DATA_ROW',                   FIRST_DATA_ROW;
    'TIME_COLUMN',                      TIME_COLUMN;
    'TIME_CONVERSION_FACTOR',           TIME_CONVERSION_FACTOR;
    'MT_LENGTH_COLUMN_INDEX',           MT_LENGTH_COLUMN_INDICES;
    '',                                 [];
    'MIN_TIME_STEP',                    MIN_TIME_STEP;
    'ERROR_TOLERANCE_LEVEL',            ERROR_TOLERANCE_LEVEL;
    '',                                 [];
    'FLATSTUT_MAX_SEGMENTHEIGHT_THRESHOLD_INPUT'	FLATSTUT_MAX_SEGMENTHEIGHT_THRESHOLD_INPUT;
    'FLATSTUT_MAX_SEGMENTSLOPE_THRESHOLD'           FLATSTUT_MAX_SEGMENTSLOPE_THRESHOLD;
    'NUC_HEIGHT_THRESHOLD'              NUC_HEIGHT_THRESHOLD;
    '',                                 [];
    'KMEANS_NumClust_All'               KMEANS_NumClust_All;
    'KMEANS_NumClust_PosSlope'          KMEANS_NumClust_PosSlope;
    'KMEANS_NumClust_NegSlope'          KMEANS_NumClust_NegSlope
    '',                                 [];
    'MIN_PROMINENCE_FOR_MAJOR_PEAKS',   MIN_PROMINENCE_FOR_MAJOR_PEAKS;
    'MIN_PEAK_HEIGHT',                  MIN_PEAK_HEIGHT;
    '',                                 [];
	'WeightedAvgSlope_in_Gr2',          WeightedAvgSlope_in_Gr2;
    'WeightedAvgSlope_in_Gr1',          WeightedAvgSlope_in_Gr1;
    'WeightedAvgSlope_in_UpSt',         WeightedAvgSlope_in_UpSt;
    'WeightedAvgSlope_in_FlSt',         WeightedAvgSlope_in_FlSt;
    'WeightedAvgSlope_in_DnSt',         WeightedAvgSlope_in_DnSt;
    'WeightedAvgSlope_in_Sh1',          WeightedAvgSlope_in_Sh1;
    'WeightedAvgSlope_in_Sh2',          WeightedAvgSlope_in_Sh2;
    'WeightedAvgSlope_in_Nuc',          WeightedAvgSlope_in_Nuc;
    '',                                 [];
    'MedSlope_in_Gr2',                  MedSlope_in_Gr2;
    'MedSlope_in_Gr1',                  MedSlope_in_Gr1;
    'MedSlope_in_UpSt',                 MedSlope_in_UpSt;
    'MedSlope_in_FlSt',                 MedSlope_in_FlSt;
    'MedSlope_in_DnSt',                 MedSlope_in_DnSt;
    'MedSlope_in_Sh1',                  MedSlope_in_Sh1;
    'MedSlope_in_Sh2',                  MedSlope_in_Sh2;
    'MedSlope_in_Nuc',                  MedSlope_in_Nuc;
    '',                                 [];
    'AvgSlope_in_Gr2',                  AvgSlope_in_Gr2;
    'AvgSlope_in_Gr1',                  AvgSlope_in_Gr1;
    'AvgSlope_in_UpSt',                 AvgSlope_in_UpSt;
    'AvgSlope_in_FlSt',                 AvgSlope_in_FlSt;
    'AvgSlope_in_DnSt',                 AvgSlope_in_DnSt;
    'AvgSlope_in_Sh1',                  AvgSlope_in_Sh1;
    'AvgSlope_in_Sh2',                  AvgSlope_in_Sh2;
    'AvgSlope_in_Nuc',                  AvgSlope_in_Nuc;
    '',                                 [];
    'StDevSlope_in_Gr2',                StDevSlope_in_Gr2;
    'StDevSlope_in_Gr1',                StDevSlope_in_Gr1;
    'StDevSlope_in_UpSt',               StDevSlope_in_UpSt;
    'StDevSlope_in_FlSt',               StDevSlope_in_FlSt;
    'StDevSlope_in_DnSt',               StDevSlope_in_DnSt;
    'StDevSlope_in_Sh1',                StDevSlope_in_Sh1;
    'StDevSlope_in_Sh2',                StDevSlope_in_Sh2;
    'StDevSlope_in_Nuc',                StDevSlope_in_Nuc;
    '',                                 [];
    'TotTime_in_Gr2_in_minutes',        TotSegTimes_in_minutes(1);
    'TotTime_in_Gr1_in_minutes',        TotSegTimes_in_minutes(2);
    'TotTime_in_UpSt_in_minutes',       TotSegTimes_in_minutes(3);
    'TotTime_in_FlSt_in_minutes',       TotSegTimes_in_minutes(4);
    'TotTime_in_DnSt_in_minutes',       TotSegTimes_in_minutes(5);
    'TotTime_in_Sh1_in_minutes',        TotSegTimes_in_minutes(6);
    'TotTime_in_Sh2_in_minutes',        TotSegTimes_in_minutes(7);
    'TotTime_in_Nuc_in_minutes',        TotSegTimes_in_minutes(8);
    '',                                 [];
    'Number_of_Abrupt_Catastrophe',     Number_of_Abrupt_Catastrophe;
    'Number_of_Transitional_Catastrophe',   Number_of_Transitional_Catastrophe;
    'Number_of_Abrupt_Rescue',          Number_of_Abrupt_Rescue;
    'Number_of_Transitional_Rescue',    Number_of_Transitional_Rescue;
    'Number_of_Interupted_Growth',      Number_of_Interupted_Growth;
    'Number_of_Interupted_Shortening',  Number_of_Interupted_Shortening;
    '',                                 [];
    'Frequency_of_Abrupt_Catastrophe',  Frequency_of_Abrupt_Catastrophe;
    'Frequency_of_Transitional_Catastrophe',    Frequency_of_Transitional_Catastrophe;
    'Frequency_of_Abrupt_Rescue',       Frequency_of_Abrupt_Rescue;
    'Frequency_of_Transitional_Rescue', Frequency_of_Transitional_Rescue;
    'Frequency_of_Interupted_Growth',   Frequency_of_Interupted_Growth;
    'Frequency_of_Interupted_Shortening',   Frequency_of_Interupted_Shortening;
    '',                                 [];
    'Number_of_ACat_ARes',              Number_of_ACat_ARes;
    'Number_of_TCat_ARes',              Number_of_TCat_ARes;
    'Number_of_ACat_TRes',              Number_of_ACat_TRes;
    'Number_of_TCat_TRes',              Number_of_TCat_TRes;
    'Number_of_IntGr_ACat',             Number_of_IntGr_ACat;
    'Number_of_IntGr_TCat',             Number_of_IntGr_TCat;
    'Number_of_IntSh_ARes',             Number_of_IntSh_ARes;
    'Number_of_IntSh_TRes',             Number_of_IntSh_TRes;
    'Number_of_ACat_IntSh',             Number_of_ACat_IntSh;
    'Number_of_TCat_IntSh',             Number_of_TCat_IntSh;
    'Number_of_ARes_IntGr',             Number_of_ARes_IntGr;
    'Number_of_TRes_IntGr',             Number_of_TRes_IntGr});


Export_Table3 = table(Abrupt_Catastrophe_Start_Index', ...
    Abrupt_Catastrophe_Start_Time', ...
    Abrupt_Catastrophe_Start_Height');
Export_Table3.Properties.VariableNames = {'Abrupt_Catastrophe_Start_Index', ...
    'Abrupt_Catastrophe_Start_Time', ...
    'Abrupt_Catastrophe_Start_Height'};

Export_Table4 = table(Abrupt_Rescue_Start_Index',...
    Abrupt_Rescue_Start_Time',...
    Abrupt_Rescue_Start_Height');
Export_Table4.Properties.VariableNames = {'Abrupt_Rescue_Start_Index',...
    'Abrupt_Rescue_Start_Time',...
    'Abrupt_Rescue_Start_Height'};

Export_Table5 = table(Transitional_Catastrophe_Stutter_Start_Index',...
    Transitional_Catastrophe_Shortening_Start_Index',...
    Transitional_Catastrophe_Stutter_Start_Time',...
    Transitional_Catastrophe_Shortening_Start_Time',...
    Transitional_Catastrophe_Stutter_Start_Height',...
    Transitional_Catastrophe_Shortening_Start_Height');
Export_Table5.Properties.VariableNames = {'Transitional_Catastrophe_Stutter_Start_Index',...
    'Transitional_Catastrophe_Shortening_Start_Index',...
    'Transitional_Catastrophe_Stutter_Start_Time',...
    'Transitional_Catastrophe_Shortening_Start_Time',...
    'Transitional_Catastrophe_Stutter_Start_Height',...
    'Transitional_Catastrophe_Shortening_Start_Height'};

Export_Table6 = table(Transitional_Rescue_Stutter_Start_Index',...
    Transitional_Rescue_Growth_Start_Index',...
    Transitional_Rescue_Stutter_Start_Time',...
    Transitional_Rescue_Growth_Start_Time',...
    Transitional_Rescue_Stutter_Start_Height',...
    Transitional_Rescue_Growth_Start_Height');
Export_Table6.Properties.VariableNames = {'Transitional_Rescue_Stutter_Start_Index',...
    'Transitional_Rescue_Growth_Start_Index',...
    'Transitional_Rescue_Stutter_Start_Time',...
    'Transitional_Rescue_Growth_Start_Time',...
    'Transitional_Rescue_Stutter_Start_Height',...
    'Transitional_Rescue_Growth_Start_Height'};

Export_Table7 = table(Interupted_Growth_Stutter_Start_Index',...
    Interupted_Growth_Growth_Start_Index',...
    Interupted_Growth_Stutter_Start_Time',...
    Interupted_Growth_Growth_Start_Time',...
    Interupted_Growth_Stutter_Start_Height',...
    Interupted_Growth_Growth_Start_Height');
Export_Table7.Properties.VariableNames = {'Interupted_Growth_Stutter_Start_Index',...
    'Interupted_Growth_Growth_Start_Index',...
    'Interupted_Growth_Stutter_Start_Time',...
    'Interupted_Growth_Growth_Start_Time',...
    'Interupted_Growth_Stutter_Start_Height',...
    'Interupted_Growth_Growth_Start_Height'};

Export_Table8 = table(Interupted_Shortening_Stutter_Start_Index',...
    Interupted_Shortening_Shortening_Start_Index',...
    Interupted_Shortening_Stutter_Start_Time',...
    Interupted_Shortening_Shortening_Start_Time',...
    Interupted_Shortening_Stutter_Start_Height',...
    Interupted_Shortening_Shortening_Start_Height');
Export_Table8.Properties.VariableNames = {'Interupted_Shortening_Stutter_Start_Index',...
    'Interupted_Shortening_Shortening_Start_Index',...
    'Interupted_Shortening_Stutter_Start_Time',...
    'Interupted_Shortening_Shortening_Start_Time',...
    'Interupted_Shortening_Stutter_Start_Height',...
    'Interupted_Shortening_Shortening_Start_Height'};
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







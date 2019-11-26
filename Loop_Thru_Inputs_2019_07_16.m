% Input_and_Run.m is the file that you run, and it calls this file,
% which calls the following Matlab scripts:
%         1) PieceWiseLinearApproximation.m
%         2) DIphaseClassification.m  or  PredefinedDIphaseClassification.m
%   and   3) ExtractDIparameters.m
% Then, all files are saved



% Read input data files
tic
if(SKIP_FILE_READ==1)
    display('-Skipped Read Data File.');
else
    for file_index = 1:length(FILE_NAME_INPUT)
        
        FILE_NAME = FILE_NAME_INPUT{file_index};
        display(['-Begin Read Data File ', FILE_NAME]);
 
        M_new = dlmread(FILE_NAME, INPUT_FILE_DELIMETER, FIRST_DATA_ROW - 1, 0);
        Time_new = M_new(:, TIME_COLUMN) * TIME_CONVERSION_FACTOR;
        
        % Structure of length history data, M:
        % 1st column = Time values, 2nd column = Length values
        for( LengthColIndex = MT_LENGTH_COLUMN_INDICES )
            if( file_index==1 && LengthColIndex==MT_LENGTH_COLUMN_INDICES(1) ) 
                % For first file & first length data column only
                M = [Time_new, M_new(:, LengthColIndex)];
                Time = M(:, 1);
            else
                % Stitch length histories from all subsequent data columns and files to M
                Time_end = Time(end);
                This_Time_new = Time_new + Time_end + 11*MIN_TIME_STEP_INPUT; %Add time shift for next length history times
                %Add buffer between sampled length history data (negative length values for separations)
                M = [M;...
                    [(Time_end + MIN_TIME_STEP_INPUT), -ERROR_TOLERANCE_LEVEL];...
                    [(Time_end + 10*MIN_TIME_STEP_INPUT), -ERROR_TOLERANCE_LEVEL];...
                    [This_Time_new, M_new(:, LengthColIndex)]];
                Time = M(:, 1);
            end
        end
        
        display('--File Read Completed.');
    end
end
toc
    
%for MT_LENGTH_COLUMN_INDEX = MT_LENGTH_COLUMN_INDICES
if(SKIP_FILE_READ==0)
    MT_length = M(:, 2);
end

for MIN_TIME_STEP = MIN_TIME_STEP_INPUT
    for FLATSTUT_MAX_SEGMENTHEIGHT_THRESHOLD = FLATSTUT_MAX_SEGMENTHEIGHT_THRESHOLD_INPUT
        for FLATSTUT_MAX_SEGMENTSLOPE_THRESHOLD = FLATSTUT_MAX_SEGMENTSLOPE_THRESHOLD_INPUT
            for NUC_HEIGHT_THRESHOLD = NUC_HEIGHT_THRESHOLD_INPUT

                % Segment Length History Plot in DI data
                % Segmentation via a piecewise linear approximation
                display('-Begin Length History Plot Segmentation and Linear Approximation.');
                PieceWiseLinearApproximation_2019_06_03
                display('--Length History Plot Segmentation & Approximation Completed.');
                toc

                % Classify each segment into DI Phases
                display('-Begin DI Segment Phase Classification.');
                % Classify DI phases using k-means
                    DIphaseClassification_2019_07_09

%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %
%  %  % Transfer Learning Options to be added in next STADIA Version %  %  %  
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %   
%                 if(SET_BENCHMARK_CLASSES==1)
%                     % Classify DI phases using k-means
%                     DIphaseClassification_2019_07_09
%                 else
%                     % Classify DI phases using benchmark info
%                     PredefinedDIphaseClassification_2017_07_29
%                 end
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  
                display('--DI Phase Classification Complete');
                toc


                % Extract DI parameters
                display('-Begin DI Parameter Extraction.');
                ExtractDIparameters_2019_06_03
                display('--DI Parameter Extraction Complete.');
                toc


                % Save output
                
                % Remove filename suffix
                filename = sprintf('%s',char(FILE_NAME));
                filename_prefix = split(filename,".");
                filename_prefix = filename_prefix{1};
                % Remove directories
                filename_prefix = split(filename_prefix, "/"); % for Non-Windows style directory
                FILE_NAME = filename_prefix{length(filename_prefix)};
                filename_prefix = split(filename_prefix, "\"); % for Windows style directory
                FILE_NAME = filename_prefix{length(filename_prefix)};
                
                
                this_datetime = datetime;
                DateTimeString = sprintf('_%.0f_%02dh%.0fmin%.0fsec', ...
                    yyyymmdd(this_datetime), hour(this_datetime), ...
                    minute(this_datetime), round(second(this_datetime)));
                
                Dir_name = "Output_for_" + FILE_NAME + DateTimeString;
                mkdir(Dir_name);
                original_dir = pwd;
                cd(Dir_name);
                
                writetable(Export_Table1,...
                    [FILE_NAME '_DIsegmentPhaseData.txt'],...
                    'Delimiter','\t');
                if(SET_BENCHMARK_CLASSES == 1)
                    writetable(Benchmark_Classes_Info_table,...
                        [FILE_NAME '_BenchmarkClassificationInfo.dat'],...
                        'Delimiter',',');
                end
                writetable(Export_Table2, ...
                    [FILE_NAME '_DI_parameters_output.txt'], ...
                    'Delimiter','\t')
                writetable(Export_Table3, ...
                    [FILE_NAME '_Abrupt_Catastrophe_output.txt'], ...
                    'Delimiter','\t')
                writetable(Export_Table4, ...
                    [FILE_NAME '_Abrupt_Rescue_output.txt'], ...
                    'Delimiter','\t')
                writetable(Export_Table5, ...
                    [FILE_NAME '_Transitional_Catastrophe_output.txt'], ...
                    'Delimiter','\t')
                writetable(Export_Table6, ...
                    [FILE_NAME '_Transitional_Rescue_output.txt'], ...
                    'Delimiter','\t')
                writetable(Export_Table7, ...
                    [FILE_NAME '_Interupted_Growth_output.txt'], ...
                    'Delimiter','\t')
                writetable(Export_Table8, ...
                    [FILE_NAME '_Interupted_Shortening_output.txt'], ...
                    'Delimiter','\t')

                if SAVE_FIG_1 == 1 & PLOT_FIG_1 == 1
                    fig1filename = sprintf('%s_PWLinearFit_Tol%d',FILE_NAME,ERROR_TOLERANCE_LEVEL);
                    savefig(fig1handle ,fig1filename)
                    print(fig1handle,'-dpng','-r0', fig1filename);
                    close(fig1handle);
                end
                if SAVE_FIG_2 == 1 & PLOT_FIG_3 == 1
                    fig2filename = sprintf('%s_PosNegSlopeClusterResults',FILE_NAME);
                    savefig(fig2handle ,fig2filename)
                    print(fig2handle,'-dpng','-r0', fig2filename);
                    close(fig2handle);
                end
                if SAVE_FIG_3 == 1 & PLOT_FIG_3 == 1
                    fig3filename = sprintf('%s_DIPhaseClassification',FILE_NAME);
                    savefig(fig3handle ,fig3filename)
                    print(fig3handle,'-dpng','-r0', fig3filename);
                    close(fig3handle);
                end
                if SAVE_FIG_4 == 1 & PLOT_FIG_4 == 1
                    fig4filename = sprintf('%s_DIPhasesOnLengthHistory',FILE_NAME);
                    savefig(fig4handle ,fig4filename)
                    print(fig4handle,'-dpng','-r0', fig4filename);
                    close(fig4handle);
                end
                if SAVE_FIG_5 == 1 & PLOT_FIG_5 == 1
                    fig5filename = sprintf('%s_AvgDIMeasurements',FILE_NAME);
                    savefig(fig5handle ,fig5filename)
                    print(fig5handle,'-dpng','-r0', fig5filename);
                    close(fig5handle);
                end
                if SAVE_FIG_6 == 1 & PLOT_FIG_6 == 1
                    fig6filename = sprintf('%s_TotDIMeasurements',FILE_NAME);
                    savefig(fig6handle ,fig6filename)
                    print(fig6handle,'-dpng','-r0', fig6filename);
                    close(fig6handle);
                end
                if SAVE_FIG_7 == 1 & PLOT_FIG_7 == 1
                    fig7filename = sprintf('%s_DIPhaseChangeResults',FILE_NAME);
                    savefig(fig7handle ,fig7filename)
                    print(fig7handle,'-dpng','-r0', fig7filename);
                    close(fig7handle);
                end
                
                % Move back to the previous directory
                cd(original_dir);


                if SAVE_MATLAB_WORKSPACE == 1
                    %    close all
                    save(['DI_parameters_output_' FILE_NAME DateTimeString], ...
                        '-regexp', '^(?!(M|Time|MT_length)$).')
                end

            end
        end
    end
end
%end


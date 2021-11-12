%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RDE aft-end image processing call script
%Written by: Zach Ayers
%Created: 5/25/2021
%Last updated: 11/12/2021

%%%Inputs
    %impath     -- For .tif or .tiff files, path to folder containing
                    %images. For .cine files, path to file including file
                    %name.
    %filetype   -- Image filetype given as '.tif', '.tiff', or '.cine'
    %framerate  -- Frame rate of the camera in Hz
    %d_ann      -- RDE annulus mean diameter in meters
    %imrange    -- Bounds of the frame range to be processed, given as 
                    %[start, end]
    %batchsize  -- Number of frames to be processed in each batch -
                    %normally determined by how much the annulus moves in 
                    %the frame and how much memory is available
    %rad_range  -- A range that the annulus radius falls in (in pixels),
                    %only for use with the 'hough' circle detection option
    %polbins    -- Number of polar bins to establish when creating the DSP
    %ann_detect -- Annulus detection method (given as string):
                    %'freq' - use frequency filtering and Taubin circle fit
                        %(more time consuming, more accurate)
                    %'hough' - use a Hough transform to quickly detect a
                        %full annulus
                    
    %dsp_plot   -- Boolean option to plot the detonation surface plot
    %mode_freq_plot -- Boolean option to plot the wave mode vs. frequency
                        %plot
    %wave_stat_plots  -- Boolean option to plot wave mode and
                        %batch-averaged wave velocity
 
%%%Returns:
    %All requested plots
    %TD (testdata class instance) - contains all reduced data
    
%Developer notes:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
%% Inputs
impath          = strcat('C:\Users\zacha\OneDrive - purdue.edu\',...
    'Projects\NSTGRO\Data\20211025\Exhaust_Imaging\Test1.mraw');
filetype        = '.mraw';
framerate       = 1e5;
d_ann           = mean([5.34,4.74])*.0254;
imrange         = [1 999];
batchsize       = 100;
rad_range       = [25 55];
polbins         = 200;
ann_detect      = 'freq';

dsp_plot            = 1;
mode_freq_plot      = 1;
wave_stat_plots     = 1;

print_summary       = 0;
testnum             = 1;
summary_file        = 'example.xlsx';
%% Create testdata instance and perform batch processing
TD = testdata(impath,filetype,framerate,d_ann,rad_range,imrange);
TD.im_process(batchsize,polbins,ann_detect);

%% Plotting
if wave_stat_plots
    TD.wavestats_plot()
end

if mode_freq_plot
    TD.plot_mfp(true);
end

if dsp_plot
    TD.plot_dsp([1,999],true);
end

%% Write summary data to test summary spreadsheet
if print_summary
    det_frac = sum((~isnan(TD.num_cw)).*(~isnan(TD.num_ccw)))/...
        length(TD.num_cw);          %Fraction of the test showing detonations
    
    num_cw = mode(TD.num_cw(TD.num_cw>0));
    num_ccw = mode(TD.num_ccw(TD.num_cw>0));
    if isnan(num_cw)
        num_cw = 0;
    end
    if isnan(num_ccw)
        num_ccw = 0;
    end
    
    if num_cw >0
        v_cw = mean(TD.velos_cw(TD.num_cw==num_cw));
    else
        v_cw = '';
    end
    if num_ccw>0
        v_ccw = mean(TD.velos_ccw(TD.num_ccw==num_cw));
    else
        v_ccw = '';
    end
    
    T = table(...
        det_frac,...
        num_cw,...
        num_ccw,...
        v_cw,...
        v_ccw...
    );

    writerow = testnum + 2;
    rangecode = sprintf('J%0.0d',writerow);
    writetable(T,summary_file,'Sheet',1,'Range',rangecode,...
        'WriteVariableNames',false)
end
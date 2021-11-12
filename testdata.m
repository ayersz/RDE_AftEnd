classdef testdata < handle
   properties 
       impath
       filetype
       imrange
       framerate
       d_ann
       rad_range
       batch_time
       polbins
       dsp
       mfp
       velos_cw
       velos_ccw
       num_cw
       num_ccw
       ccw_fpk
       cw_fpk
       cntr
       radius
   end
   
   methods
       %% Constructor method to initialize class values
       function TD = testdata(impath,filetype,framerate,d_ann,rad_range,...
               imrange)
           %If no imrange is given, use all the images
            if nargin < 5
                tag = strcat('*',filetype);
                filenames = dir(fullfile(strcat(impath,'\',tag)));
                TD.imrange = [1,length(filenames)];
            end 
            
            %Populate attributes
            TD.impath = impath;
            TD.filetype = filetype;
            
            %Check that the Phantom SDK is on the path and add it if not
            if strcmp(TD.filetype,'.cine')
                try
                    LoadPhantomLibraries();
                catch
                    ph_path = 'U:\My Documents\Phantom\PhSDK792';
                    pathCell = regexp(path,pathsep,'split');
                    if ~any(strcmpi(ph_path, pathCell))
                        addpath(genpath(ph_path));
                    end
                end
                
                try
                    LoadPhantomLibraries();
                catch
                    error('Failed to load phantom libraries. Make sure that the Phantom SDK is installed and the libraries are on the Matlab path.');
                end
                RegisterPhantom(true);
            end
            
            TD.framerate = framerate;
            TD.d_ann = d_ann;
            if isempty(imrange)
                TD.imrange = TD.numframes();
            else
                TD.imrange = imrange;
            end
            TD.rad_range = rad_range;
            
       end
       
       
        %% Method to run image processing
        function im_process(TD,batchsize,polbins,ann_method)
            TD.polbins = polbins;
            batches = floor((TD.imrange(2)-TD.imrange(1))/batchsize)+1;
            TD.batch_time = (1:batches)*1000*(batchsize/TD.framerate)...
                -1000*(round(batchsize/2)/TD.framerate);
            
            
            for batch = 1:batches
                fprintf('Processing batch %0.0f of %0.0f\n',...
                    [batch,batches])
        
                %If this is the final batch, only use the images that are
                %left
                if batch == batches
                    framebords = [batchsize*(batch-1)+TD.imrange(1),...
                        TD.imrange(2)];
                else
                    framebords = [batchsize*(batch-1)+TD.imrange(1),...
                        batchsize*batch + TD.imrange(1) - 1];
                end

                %Load images
                imgs = TD.img_import(framebords);
                
                %Background subtraction
                imgs = TD.bsub(imgs);
                
                %Circle identification
                if strcmp(ann_method,'freq')
                    try
                        [center, rad] = TD.findcirc_freq(imgs,3,false);
                    catch
                        center = NaN;
                        rad = NaN;
                    end
                elseif strcmp(ann_method,'hough')
                    [center, rad] = TD.findcirc_hough(imgs);
                end
                
                if isnan(rad)
                    if ~isempty(TD.cntr)
                        center = TD.cntr(batch-1,:);
                        rad = TD.radius(batch-1);
                        TD.cntr(batch,:) = center;
                        TD.radius(batch) = rad;
                        %Create DSP
                        dsp_temp = TD.make_dsp(imgs,new_center);
                        if batch == 1
                            TD.dsp = dsp_temp;
                        else
                            TD.dsp = cat(2,TD.dsp,dsp_temp);
                        end
                    else
                        dsp_temp = zeros(TD.polbins,size(imgs,3));
                        if batch == 1 
                            TD.dsp = dsp_temp;
                        else
                            TD.dsp = cat(2,TD.dsp,dsp_temp);
                        end
                    end
                    
                    TD.num_cw(batch) = NaN;
                    TD.num_ccw(batch) = NaN;
                    TD.velos_cw(batch) = NaN;
                    TD.velos_ccw(batch) = NaN;
                    TD.cw_fpk(batch) = NaN;
                    TD.ccw_fpk(batch) = NaN;
                else
                    if isempty(center)
                        fprintf('No annulus detected using Hough transform\n')
                        [center,rad] = TD.findcirc_Taubin(imgs);
                    end

                    if isempty(center)
                        fprintf('No annulus detected using Taubin fit\n')
                    end

                    %Image cropping
                    [imgs,new_center] = TD.crop_ims(imgs,center,rad);

                    %Create DSP
                    dsp_temp = TD.make_dsp(imgs,new_center);
                    if batch == 1
                        TD.dsp = dsp_temp;
                    else
                        TD.dsp = cat(2,TD.dsp,dsp_temp);
                    end

                    %Calculate operating statistics
                    [cw,ccw] = TD.av_wavestats(dsp_temp);
                    TD.num_cw(batch) = cw.num;
                    TD.num_ccw(batch) = ccw.num;
                    TD.velos_cw(batch) = cw.v;
                    TD.velos_ccw(batch) = ccw.v;
                    TD.cw_fpk(batch) = cw.fpk;
                    TD.ccw_fpk(batch) = ccw.fpk;
                    TD.cntr(batch,:) = center;
                    TD.radius(batch,:) = rad;
                end
                
            end
        end
        
        
        %% Method for plotting wave mode and speed
        function [] = wavestats_plot(TD)
            %change zero wave numbers to nan so they don't get plotted
            TD.num_cw(TD.num_cw==0) = nan;
            TD.num_ccw(TD.num_ccw==0) = nan;
            
            %Wave mode plot
            figure
            plot(TD.batch_time,-TD.num_cw,'bx','MarkerSize',10)
            hold on
            plot(TD.batch_time,TD.num_ccw,'ro','MarkerSize',10)
            ylim([-12,12])
            yticks(-12:3:12)
            set(gca,'FontSize',14,'FontName','Arial');
            patch([min(xlim) max(xlim) max(xlim) min(xlim)],...
                [mean(ylim) mean(ylim) max(ylim) max(ylim)],...
                [0 0 0],'FaceAlpha',0.4,'EdgeColor','none')
            plot(TD.batch_time,-TD.num_cw,'bx','MarkerSize',10)
            plot(TD.batch_time,TD.num_ccw,'ro','MarkerSize',10)
            text(min(xlim)+0.05*(max(xlim)-min(xlim)),10,'CCW')
            text(min(xlim)+0.05*(max(xlim)-min(xlim)),-10,'CW')
            xlabel('Time [ms]')
            ylabel('Number of Waves')
            hold off
            
            %Wave velocity plot
            figure
            plot(TD.batch_time,TD.velos_cw,'bx','MarkerSize',10)
            hold on
            plot(TD.batch_time,TD.velos_ccw,'ro','MarkerSize',10)
            xlabel('Time [ms]')
            ylabel('Wave Speed [m/s]')
            legend('Clockwise','Counter-Clockwise')
            set(gca,'FontSize',14,'FontName','Arial');
            hold off
            
        end
       
        %% Method for finding the frequency peaks, wave modes, and average
        %   velocities in a dsp
        function [cw,ccw] = av_wavestats(TD,dsp)
            %Create frequency power surface
            Y = fft2(dsp);
            L = size(dsp,2);
            P2 = abs(Y/L);
            mfp_now = P2(:,1:round(L/2));
            mfp_now(:,2:end-1) = 2*(mfp_now(:,2:end-1));
            mfp_now = cat(1,flipud(mfp_now(1:size(mfp_now,1)/2,:)),...
                flipud(mfp_now(size(mfp_now,1)/2+1:end,:)));
            mfp_now = mfp_now/max(mfp_now(:));
            
            %Create the axes for the surface
            f = TD.framerate*(0:(L/2))/L/1000;
            w = -round(size(mfp_now,1))/2:round(size(mfp_now,1)/2);
            w = w(2:end);
            
            [F,W] = meshgrid(f,w);
            
            %Extract the max amplitudes for cw and ccw modes
            p = reshape(mfp_now,[1,numel(mfp_now)]);
            Fvec = reshape(F,[1,numel(F)]);
            Wvec = reshape(W,[1,numel(W)]);
            
            [amp,locs] = findpeaks(p,'MinPeakHeight',0.3);
            [a,b] = sort(amp,'descend');
            locs = locs(b);
            
            use_w = Wvec(locs);     
            use_f = Fvec(locs);
            
            ipos = find(use_w>0 & use_f>1 & use_f<f(end)*0.9);
            if ~isempty(ipos)
                %   Prevent overtones from being selected as the dominant mode
                if length(ipos) > 1
                    if a(ipos(2)) > 0.3 && use_w(ipos(2)) == use_w(ipos(1))/2
                        ccw_pk = ipos(2);
                    else
                        ccw_pk = ipos(1);
                    end
                else
                    ccw_pk = ipos(1);
                end
            
                ccw.num = use_w(ccw_pk);
                ccw.fpk = use_f(ccw_pk);
                ccw.v = pi*TD.d_ann*(1000*ccw.fpk)/ccw.num;
            else
                ccw.num = 0;
                ccw.fpk = nan;
                ccw.v = nan;
            end
            
            ineg = find(use_w<0 & use_f>1 & use_f<f(end)*0.9);
            if ~isempty(ineg)
                %   Prevent overtones from being selected as the dominant mode
                if length(ineg) > 1
                    if a(ineg(2)) > 0.3 && a(ineg(2)) == a(ineg(1))/2
                        cw_pk = ineg(2);
                    else
                        cw_pk = ineg(1);
                    end
                else
                    cw_pk = ineg(1);
                end
                
                cw.num = abs(use_w(cw_pk));
                cw.fpk = use_f(cw_pk);
                cw.v = pi*TD.d_ann*(1000*cw.fpk)/cw.num;
            else
                cw.num = 0;
                cw.fpk = nan;
                cw.v = nan;
            end

        end
        
        %% Method for plotting the mode-frequency plot (MFP)
        function [f,w,mfp] = plot_mfp(TD,showplot,range)
            %if no range is given, perform frequency analysis on the entire
            %detonation surface plot (TD.dsp)
            if nargin < 3
                range = [1 size(TD.dsp,2)];
            end
            
            dsp_use = TD.dsp(:,range(1):range(2));
            
            %Create frequency amplitude surface
            Y = fft2(dsp_use);
            L = size(dsp_use,2);
            P2 = abs(Y/L);
            TD.mfp = P2(:,1:round(L/2));
            TD.mfp(:,2:end-1) = 2*(TD.mfp(:,2:end-1));
            TD.mfp = cat(1,flipud(TD.mfp(1:size(TD.mfp,1)/2,:)),...
                flipud(TD.mfp(size(TD.mfp,1)/2+1:end,:)));
            TD.mfp = TD.mfp/max(TD.mfp(:));
            mfp = TD.mfp;
            
            %Create MFP axis units
            f = TD.framerate*(0:(L/2))/L/1000;
            w = -round(size(TD.mfp,1))/2:round(size(TD.mfp,1)/2);
            w = w(2:end);
            
            %Create the plot
            if showplot
                figure
                imagesc(f,w,TD.mfp)
                ax = gca;
                ax.YDir = 'normal';
                h = colorbar;
                ylabel(h, 'Normalized Amplitude [a.u.]')
                ylim([-10,10])
                xlim([0,50]);
                yticks(-10:2:10)
                xlabel('Operating Frequency, f_d_e_t, [kHz]')
                ylabel('Number of Waves, m')

                set(findall(gcf,'-property','FontSize'),'FontSize',12) 
                set(findall(gcf,'-property','Font'),'Font','Arial')
            end
        end
        
        
        %% Method to use a Hough transform to detect the annulus
        function [center, rad] = findcirc_hough(TD,imgs)
            %Apply a threshold to each image
            imgs(imgs<0.7) = 0;

            %Average the thresholded images
            avim = mean(imgs,3)./max(mean(imgs,3),[],'all');
            
            %If no radius range was given, prompt the user to draw a
            %diameter line across the annulus to approximate the diameter
            if isempty(TD.rad_range)
                f = figure;
                imshow(avim,[])
                f.Position = [550 250 600 500];
                title(strcat('Click and hold to draw a diameter line',...
                    ' across the annulus'))
                d = drawline;
                pos = d.Position;
                close(gcf)
                diffPos = diff(pos);
                diameter = hypot(diffPos(1),diffPos(2));
                rad = round(diameter/2);
                rmin = rad-15;
                rmax = rad+15;
                TD.rad_range = [rmin,rmax];
            end
            
            %Detect the annulus
            [center, rad] = imfindcircles(mean(avim,3),TD.rad_range,...
                'Sensitivity',0.9);
        end
        
        
        %% Method for detecting the annulus pixels using frequency
        %   filtering
        function [cntr,rad] = findcirc_freq(TD,imgs,binsize,show_circfit)
            if nargin < 3
                show_circfit = false;
            end
            
            %perform image binning
            numrows = round(size(imgs,1)/binsize);
            numcols = round(size(imgs,2)/binsize);
            bims = zeros(numrows,numcols,size(imgs,3));
            high_bims = zeros(numrows,numcols,size(imgs,3));
            for frame = 1:size(imgs,3)
                bims(:,:,frame) = imresize(imgs(:,:,frame),...
                    [numrows numcols]);
                bim = bims(:,:,frame);
                bim(bim<0.1*max(imgs(:))) = 0;
                bim(bim~=0) = 1;
                high_bims(:,:,frame) = bim;
            end
            
            %create a binary mask where pixels in the mean image above 10%
            %of the max value are equal to 1
            mask = mean(high_bims,3);
            mask(mask<0.1*max(mask)) = 0;
            mask(mask~=0) = 1;
            
            %frequency analysis
            L = size(bims,3);
            count = 1;
            x = zeros(size(bims,1)*size(bims,2),1);
            y = zeros(size(bims,1)*size(bims,2),1);
            for row = 1:numrows
                for col = 1:numcols
                    Y = fft(bims(row,col,:));
                    Y = squeeze(Y(1,1,:));
                    P2 = abs(Y/L);
                    P1 = P2(1:round(L/2)+1);
                    P1(2:end-1) = 2*P1(2:end-1);
                    f = TD.framerate*(0:(round(L/2)))/L;
                    
                    [pk,loc] = max(P1(2:end));
                    pk_freq = f(loc+1);
                    
                    %Freq must be btwn 2 kHz & 50 kHz and pixel is 1 in the
                    %mask and freq. power is > 50% of the zero freq. power
                    if pk_freq > 2e3 && pk_freq < 5e4 &&...
                            mask(row,col) == 1 && pk > 0.5*P1(1)
                       x(count,1) = col;
                       y(count,1) = row;
                    end
                    count = count + 1;
                end
            end
            
            %get shortlist of pixel locations to be used for circle fitting
            x = x(x~=0);
            y = y(y~=0);
            XY = [x y];
            
            %calculate binned center
            [c,r] = TD.Taubin(XY);
            
            if show_circfit
                %show the mean binned image with the pixels used for circle
                %fitting indicated by red dots
                figure
                imshow(mean(high_bims,3),[]);
                hold on
                plot(x,y,'r.')
                plot(c(1),c(2),'go')
                plot(size(bims,1)/2,size(bims,2)/2,'r+')
            end
            
            %calculate unbinned center
            cntr = c*binsize;
            rad = r*binsize;
        end

        %% Method for plotting a detonation surface plot
        function [time,dsp_plot] = plot_dsp(TD,framerange,showdsp,showlines)
            if nargin < 4
                showlines = 0;
            end
            
            dsp_plot = TD.dsp(:,framerange(1):framerange(2));
            time = [framerange(1),framerange(2)].*1e3*(1/TD.framerate);
            
            if showdsp
                figure
                image(time,[0,359],dsp_plot,'CDataMapping','scaled')
                ax = gca;
                ax.YDir = 'normal';
                colormap gray
                ylabel('\theta [deg.]')
                xlabel('Time [ms]')
                title('DSP')
                hold on

                if showlines
                    [~,~,~,lines] = TD.hough_lines(dsp_plot);

                    for k = 1:length(lines)
                        xy = [lines(k).point1; lines(k).point2];
                        xy(:,2) = xy(:,2)*(360/TD.polbins);
                        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

                        % Plot beginnings and ends of lines
                        plot(xy(1,1),xy(1,2)-1,'x','LineWidth',2,'Color','yellow');
                        plot(xy(2,1),xy(2,2)-1,'x','LineWidth',2,'Color','red');
                    end
                end
                hold off
            end
        end
        
        %% Method for finding the wave mode and wave speeds
       function [velos,num_cw,num_ccw,lines] = hough_lines(TD,dsp)
            %Identify lines in the DSP
            %Sample different thetas and find the number of waves to look for
            timevec = dsp(round(2*size(dsp,1)/5),:);
            pks1 = findpeaks(timevec,'MinPeakHeight',...
                0.2*max(timevec));
            timevec = dsp(round(size(dsp,1)/3),:);
            pks2 = findpeaks(timevec,'MinPeakHeight',...
                0.2*max(timevec));
            timevec = dsp(round(2*size(dsp,1)/3),:);
            pks3 = findpeaks(timevec,'MinPeakHeight',...
                0.2*max(timevec));
            numlines = median([length(pks1),length(pks2),length(pks3)])-1;

            %Binarize DSP and use Hough transform to detect lines
            temp = dsp;
            temp(temp<0.2*max(temp)) = 0;
            temp(temp>=0.2*max(temp)) = 1;
            [H,theta,rho] = hough(temp);
            P = houghpeaks(H,numlines,'threshold',ceil(0.3*max(H(:))));
            lines = houghlines(temp,theta,rho,P,'MinLength',100);

            %Convert line points to time and calculate velocities
            velos = zeros(length(lines),1);
            for k = 1:length(lines)
                lines(k).point1(1) = lines(k).point1(1)...
                    *(1000/TD.framerate);
                lines(k).point2(1) = lines(k).point2(1)...
                    *(1000/TD.framerate);

                %Get m and b for the equation of each line (y=m*x+b)
                lines(k).m = (lines(k).point2(2)-lines(k).point1(2))...
                    /(lines(k).point2(1)-lines(k).point1(1));
                lines(k).b = lines(k).point1(2) -...
                    lines(k).m*lines(k).point1(1);

                %Calculate velocities
                c_ann = pi*TD.d_ann;   %Annulus mean circumference
                velos(k) = 1000*(c_ann/TD.polbins)*lines(k).m;
            end

            %Identify the wave mode
            time = (1:(size(dsp,2)))*(1000/TD.framerate);
            ws = zeros(1,size(dsp,2));
            for col = 1:size(dsp,2)
               w = 0;
               for k = 1:length(lines)
                  y = lines(k).m*time(col)+lines(k).b;
                  if y < 360 && y > 0
                      w = w + 1;
                  end              
               end   
               ws(col) = w;
            end
            numwaves = round(mean(ws));
            num_cw = round(numwaves*(sum(velos<0)/length(velos)));
            num_ccw = round(numwaves*(sum(velos>0)/length(velos)));

       end
       
       %% Method for getting the number of frames in a .cine file
       function [range] = numframes(TD)
            switch TD.filetype
                case '.cine'
                   %Get the frame numbers of the .cine file
                    [HRES, cineHandle] = PhNewCineFromFile(TD.impath);
                    if (HRES<0)
                        [message] = PhGetErrorMessage( HRES );
                        error(['Cine handle creation error: ' message]);
                    end
        
                    pFirstIm = libpointer('int32Ptr',0);
                    PhGetCineInfo(cineHandle, PhFileConst.GCI_FIRSTIMAGENO,...
                        pFirstIm);
                    firstIm = pFirstIm.Value;
                    pImCount = libpointer('uint32Ptr',0);
                    PhGetCineInfo(cineHandle, PhFileConst.GCI_IMAGECOUNT,...
                        pImCount);
                    lastIm = int32(double(firstIm) + double(pImCount.Value)...
                        - 1);
                    
                    range = double([1,lastIm+1-firstIm]);
                case '.mraw'
                    fid1 = fopen(sprintf('%s.cih',TD.impath(1:end-5)),'r');
                    if fid1 < 0
                        fid1 = fopen(sprintf('%s.cihx',TD.impath(1:end-5)),'r');
                        cihx = true;
                    else
                        cihx = false;
                    end
                    fid2 = fopen(sprintf('%s',TD.impath),'r');
                    if fid1 < 1 || fid2 < 1
                        error(['Could not locate .CIH or .CIHX header for file: ''' filename '''']);
                    end
                    Header = textscan(fid1,'%s','delimiter',':');
                    Header = Header{1};
         
                    frame_ind = find(contains(Header, 'Total Frame')) + 1;
                    Total_Frames = str2double(cell2mat(Header(frame_ind(1))));
                    range = [1,Total_Frames];
            end
       end
        
        %% Method to load images
        function [imgs,numframes] = img_import(TD,imrange,resolution)
            if nargin < 3
                resolution = 1;
            end
            
            %.tif and .tiff file import
            if strcmp(TD.filetype,'.tif') || strcmp(TD.filetype,'.tiff')
                tag = strcat('*',TD.filetype);
                filenames = dir(fullfile(strcat(TD.impath,'\',tag)));

                %If no imrange is given, import all the images
                if nargin < 2
                    imrange = [1,length(filenames)];
                end

                %Start for loop to import each image
                testim = imread(strcat(TD.impath,'\',filenames(1).name));
                imgs = zeros(size(testim,1),size(testim,2),...
                    imrange(2)-imrange(1)+1);
                count = 1;
                for i = imrange(1):resolution:imrange(2)
                    current = filenames(i).name;
                    imgs(:,:,count) = imread(strcat(TD.impath,'\',current));
                    count = count + 1;
                end        
            
            %.cine file import
            elseif strcmp(TD.filetype,'.cine')
                %Get the frame numbers of the .cine file
                [HRES, cineHandle] = PhNewCineFromFile(TD.impath);
                if (HRES<0)
                    [message] = PhGetErrorMessage( HRES );
                    error(['Cine handle creation error: ' message]);
                end
           
                pFirstIm = libpointer('int32Ptr',0);
                PhGetCineInfo(cineHandle, PhFileConst.GCI_FIRSTIMAGENO,...
                    pFirstIm);
                firstIm = pFirstIm.Value;
                pImCount = libpointer('uint32Ptr',0);
                PhGetCineInfo(cineHandle, PhFileConst.GCI_IMAGECOUNT,...
                    pImCount);
                lastIm = int32(double(firstIm) + double(pImCount.Value)...
                    - 1);
                numframes = lastIm-firstIm + 1;
                
                %Make sure the requested range is allowable
                if imrange(2)-imrange(1)+1 > numframes
                    error('Requested images must be in the range [1,%d]',...
                        numframes);
                end
                
                %Get the size of the images
                [testim,~] = ReadCineFileImage(TD.impath,firstIm,false);
                
                %Test if the image is an RGB color image or not
                if length(size(testim)) == 3
                    iscolor = true;
                else
                    iscolor = false;
                end
                
                %Preallocate memory
                if resolution ~= 1 && iscolor
                    imgs = zeros(size(testim,1),size(testim,2),3,...
                        ceil((imrange(2)-imrange(1))/resolution));
                elseif resolution ~= 1 && ~iscolor
                    imgs = zeros(size(testim,1),size(testim,2),...
                        ceil((imrange(2)-imrange(1))/resolution));
                elseif resolution == 1 && iscolor
                    imgs = zeros(size(testim,1),size(testim,2),3,...
                        imrange(2)-imrange(1)+1);
                elseif resolution ==1 && ~iscolor
                    imgs = zeros(size(testim,1),size(testim,2),...
                        imrange(2)-imrange(1)+1);
                end
                
                %Load all the images
                count = 1;
                for i = imrange(1)-1+firstIm:resolution:imrange(2)-1+firstIm
                    if iscolor
                        [imgs(:,:,:,count),~] =...
                           ReadCineFileImage(TD.impath,i,false); 
                    else
                        [imgs(:,:,count),~] =...
                            ReadCineFileImage(TD.impath,i,false); 
                    end
                   count = count + 1;
                end
                
                %Convert RGB images to black and white
                if iscolor
                    imgs = mean(imgs,3);
                    imgs = squeeze(imgs(:,:,1,:));
                end

            %.mraw file import from Photron cameras
            elseif strcmp(TD.filetype,'.mraw')
                if resolution == 1
                    [imgs, numframes] = readmraw(TD.impath,imrange);
                    imgs = double(imgs);
                else
                    disp('Multi-image resolution for .mraw files not yet programmed.')
                end
            
            %.mat file import
            elseif strcmp(TD.filetype,'.mat')
                all_imgs = load(TD.impath);
                if nargin < 2
                    imgs = all_imgs.imgs;
                else
                    imgs = all_imgs.imgs(:,:,...
                        imrange(1):resolution:imrange(2));
                end

            end
        end

        %% Method to make a detonation surface plot
        function [dsp] = make_dsp(TD,imgs,center)
            %Find the annulus center in the cropped image
            avim = mean(imgs,3);
            avim = avim/max(avim(:));
            avim(avim<0.5) = 0;

            %Convert the average image to polar coordinates and find the bounds
            %of the annulus
            pol_avim = TD.cart2pol2d(avim,center);        

            pol_vec = mean(pol_avim,2);
            [mx,mx_idx] = max(pol_vec);
            low_r = find(pol_vec(1:mx_idx)<mx/2,1,'last');
            high_r = find(pol_vec(mx_idx:end)<mx/2,1,'first')+mx_idx;

            %Convert all images to polar, radially average across the annulus,
            %and concatenate theta vectors in time to create DSP
            imgs = TD.cart2pol3d(imgs,center);
            dsp = mean(imgs(low_r:high_r,:,:),1);
            dsp = squeeze(dsp(1,:,:));
        end
        
        %% Method for converting a 2D image from cartesian to polar
        function [polim] = cart2pol2d(TD,cartim,center)
            maxrad = floor(sqrt(center(1)^2+center(2)^2));
            polim = zeros(maxrad,TD.polbins);
            for row = 1:size(cartim,1)
                for col = 1:size(cartim,2)
                    dx = col - center(1);
                    dy = center(2) - row;
                    r = round(sqrt(dx^2+dy^2));
                    if dx > 0 && dy > 0
                        theta = round(atand(dy/dx));
                    elseif (dx < 0 && dy > 0) || (dx < 0 && dy < 0)
                        theta = 180 + round(atand(dy/dx));
                    elseif dx > 0 && dy < 0
                        theta = 360 + round(atand(dy/dx));
                    end

                    polim(r+1,ceil(TD.polbins*(theta+1)/360)) =...
                        cartim(row,col);
                end
            end

        end

        %% Method for converting a 3D image array from cartesian to polar
        function [polim] = cart2pol3d(TD,cartims,center)
            maxrad = floor(sqrt(center(1)^2+center(2)^2));
            polim = zeros(maxrad,TD.polbins,size(cartims,3));
            for row = 1:size(cartims,1)
                for col = 1:size(cartims,2)
                    dx = col - center(1);
                    dy = center(2) - row;
                    r = round(sqrt(dx^2+dy^2));
                    if dx > 0 && dy > 0
                        theta = round(atand(dy/dx));
                    elseif (dx < 0 && dy > 0) || (dx < 0 && dy < 0)
                        theta = 180 + round(atand(dy/dx));
                    elseif dx > 0 && dy < 0
                        theta = 360 + round(atand(dy/dx));
                        if theta == 360
                            theta = 0;
                        end
                    end
                    polim(r+1,ceil(TD.polbins*(theta+1)/360),:) =...
                        cartims(row,col,:);
                end
            end

        end
   end
   
   methods (Static)
       %% Background subtraction method
        function [bsub_imgs] = bsub(imgs)
            avim = mean(imgs,3);

            bsub_imgs = imgs - avim;
            bsub_imgs(bsub_imgs<0) = 0;

            %Perform image normalization - normalize to max intensity value
            %in the image set
            for frame = 1:size(bsub_imgs,3)
                bsub_imgs(:,:,frame) = bsub_imgs(:,:,frame)...
                    /max(bsub_imgs,[],'all');
            end
        end


        %% Method for cropping images about a center with a certain radius
        function [c_imgs,new_center] = crop_ims(imgs,center,rad)
            %Crop images about the detected annulus
            mult = 1.2;     %Multiplier to choose how many radii to include
                            %in the cropped image
            rowbounds = [floor(center(2)-mult*rad),...
                ceil(center(2)+mult*rad)];
            if rowbounds(1) < 1
                rowbounds(1) = 1;
            end
            if rowbounds(2) > size(imgs,1)
                rowbounds(2) = size(imgs,1);
            end

            colbounds = [floor(center(1)-mult*rad),...
                ceil(center(1)+mult*rad)];
            if colbounds(1) < 1
                colbounds(1) = 1;
            end
            if colbounds(2) > size(imgs,2)
                colbounds(2) = size(imgs,2);
            end

            c_imgs = imgs(rowbounds(1):rowbounds(2),colbounds(1):colbounds(2),:);
            new_center = [center(1,1)-colbounds(1),center(1,2)-rowbounds(1)];
        end
        
        %% Method to use a Taubin fit to detect the annulus
        %   adapted from MATLAB file exchange function:
        %   https://www.mathworks.com/matlabcentral/fileexchange/
        %   22678-circle-fit-taubin-method
        function [cntr,rad] = Taubin(XY)
                        
            n = size(XY,1);      % number of data points

            centroid = mean(XY);   % the centroid of the data set

            %     computing moments (note: all moments will be normed, i.e.
            %       divided by n)

            Mxx = 0; Myy = 0; Mxy = 0; Mxz = 0; Myz = 0; Mzz = 0;

            for i=1:n
                Xi = XY(i,1) - centroid(1);  %  centering data
                Yi = XY(i,2) - centroid(2);  %  centering data
                Zi = Xi*Xi + Yi*Yi;
                Mxy = Mxy + Xi*Yi;
                Mxx = Mxx + Xi*Xi;
                Myy = Myy + Yi*Yi;
                Mxz = Mxz + Xi*Zi;
                Myz = Myz + Yi*Zi;
                Mzz = Mzz + Zi*Zi;
            end

            Mxx = Mxx/n;
            Myy = Myy/n;
            Mxy = Mxy/n;
            Mxz = Mxz/n;
            Myz = Myz/n;
            Mzz = Mzz/n;

            %computing the coefficients of the characteristic polynomial

            Mz = Mxx + Myy;
            Cov_xy = Mxx*Myy - Mxy*Mxy;
            A3 = 4*Mz;
            A2 = -3*Mz*Mz - Mzz;
            A1 = Mzz*Mz + 4*Cov_xy*Mz - Mxz*Mxz - Myz*Myz - Mz*Mz*Mz;
            A0 = Mxz*Mxz*Myy + Myz*Myz*Mxx - Mzz*Cov_xy - 2*Mxz*Myz*Mxy...
                + Mz*Mz*Cov_xy;
            A22 = A2 + A2;
            A33 = A3 + A3 + A3;

            xnew = 0;
            ynew = 1e+20;
            epsilon = 1e-12;
            IterMax = 20;

            % Newton's method starting at x=0

            for iter=1:IterMax
                yold = ynew;
                ynew = A0 + xnew*(A1 + xnew*(A2 + xnew*A3));
                if abs(ynew) > abs(yold)
                   disp(strcat('Newton-Taubin goes wrong direction:',...
                       '|ynew| > |yold|'));
                   xnew = 0;
                   break;
                end
                Dy = A1 + xnew*(A22 + xnew*A33);
                xold = xnew;
                xnew = xold - ynew/Dy;
                if (abs((xnew-xold)/xnew) < epsilon), break, end
                if (iter >= IterMax)
                    disp('Newton-Taubin will not converge');
                    xnew = 0;
                end
                if (xnew<0.)
                    fprintf(1,'Newton-Taubin negative root:  x=%f\n',xnew);
                    xnew = 0;
                end
            end

            %  computing the circle parameters

            DET = xnew*xnew - xnew*Mz + Cov_xy;
            Center = [Mxz*(Myy-xnew)-Myz*Mxy ,...
                Myz*(Mxx-xnew)-Mxz*Mxy]/DET/2;
            
            rad = sqrt(Center*Center'+Mz);
            cntr = Center+centroid;

        end
        
        %% Method to detect the annulus using thresholding and Taubin 
        % fitting
        function [cntr,rad] = findcirc_threshold(imgs)
            %Apply a threshold to each image
            imgs(imgs<0.7) = 0;

            %Average the thresholded images
            avim = mean(imgs,3)./max(mean(imgs,3),[],'all');
            
            %Binarize the average image
            avim(avim>=0.2) = 1;
            avim(avim<0.2) = 0;
            
            %Get the [x,y] locations of each non-zero pixel
            [y,x] = ind2sub(size(avim),find(avim>0));
            XY = [x,y];
            
            [cntr,rad] = Taubin(XY);
        end
   end
   
end




function [imgs, Total_Frames] = readmraw(filename, numImgs)
%   readmraw.m 
%   READMRAW Read Photron MRAW files into MATLAB
%   
%   imgs = READMRAW('C:\Photron\Filename.mraw', [a,b]) loads images a
%   through b from 'C:\Photron\Filename.mraw' into matrix imgs. 
%
%   Remarks
%   -------
%   This function must be handed the common *.cih(x) and *.mraw file name 
%   and the range of images to be loaded (MATLAB may not handle the entire
%   image range for large files). 
%   NOTE: Both the *.cih(x) file and the *.mraw file are utilized
%   Autor: SEP                                Creation Date: June 20,2013
%   Editor: Phil Kreth                    Modification Date: Aug 31, 2021
%
%   Examples
%   --------
%   % Load all images
%   imgs = readmraw('C:\Photron\Moviefile.mraw', 0); 
% 
%   % Load images 10 through 50
%   imgs = readmraw('C:\Photron\Moviefile.mraw', [10,50]);
% 
%   % Load image 10
%   imgs = readmraw('C:\Photron\Moviefile.mraw', 10);
% 
fid1 = fopen(sprintf('%s.cih',filename(1:end-5)),'r');
if fid1 < 0
    fid1 = fopen(sprintf('%s.cihx',filename(1:end-5)),'r');
    cihx = true;
else
    cihx = false;
end
fid2 = fopen(sprintf('%s',filename),'r');
if fid1 < 1 || fid2 < 1
    error(['Could not locate .CIH or .CIHX header for file: ''' filename '''']);
end
if ~cihx % CIH FILE
    % Read Header Information
    Header = textscan(fid1,'%s','delimiter',':');
    Header = Header{1};
    color_ind = find(contains(Header, 'Color Type')) + 1;
    if strcmp(cell2mat(Header(color_ind(1))), 'Color')
        color = true;
    else
        color = false;
    end
    bit_ind = find(contains(Header, 'Color Bit')) + 1;
    bits = str2double(cell2mat(Header(bit_ind(1))));
    if color
        bits = bits/3;
    end
    bit_depth = sprintf('ubit%d', bits);
    frame_ind = find(contains(Header, 'Total Frame')) + 1;
    Total_Frames = str2double(cell2mat(Header(frame_ind(1))));
    width_ind = find(contains(Header, 'Image Width')) + 1;
    Width = str2double(cell2mat(Header(width_ind(1))));
    height_ind = find(contains(Header, 'Image Height')) + 1;
    Height = str2double(cell2mat(Header(height_ind(1))));
    % fps_ind = find(contains(Header, 'Record Rate(fps)')) + 1;
    % fps = str2double(cell2mat(Header(fps_ind(1))));
    
else % CIHX FILE
    % Read Header Information
    Header = textscan(fid1,'%s','delimiter',{'<','>'});
    Header = Header{1};
    color_ind = find(contains(Header, 'type')) + 1;
    if strcmp(cell2mat(Header(color_ind(1))), 'Color')
        color = true;
    else
        color = false;
    end
    bit_ind = find(contains(Header, 'bit')) + 1;
    bits = str2double(cell2mat(Header(bit_ind(1))));
    if color
        bits = bits/3;
    end
    bit_depth = sprintf('ubit%d', bits);
    frame_ind = find(contains(Header, 'totalFrame')) + 1;
    Total_Frames = str2double(cell2mat(Header(frame_ind(1))));
    width_ind = find(contains(Header, 'width')) + 1;
    Width = str2double(cell2mat(Header(width_ind(1))));
    height_ind = find(contains(Header, 'height')) + 1;
    Height = str2double(cell2mat(Header(height_ind(1))));
    % fps_ind = find(contains(Header, 'recordRate')) + 1;
    % fps = str2double(cell2mat(Header(fps_ind(1))));
end
Pixels = Width*Height;
fclose(fid1);
% Define Image Range
if numImgs == 0             % load all the images 
    first_frame = 1;
    frames = Total_Frames;
elseif length(numImgs) == 1 % load a single image
    first_frame = numImgs;
    frames = 1;
else                        % load a specified range of images
    first_frame = numImgs(1,1);
    last_frame = numImgs(1,2);
    frames = last_frame-first_frame+1;
end
% Load Images
bytes_offset = (first_frame-1)*Pixels*bits/8;
if color
    bytes_offset = bytes_offset*3;
end
fseek(fid2, bytes_offset, 'bof');
if bits > 8
	data_fmt = 'uint16';
else
	data_fmt = 'uint8';
end
if color
    imgs = zeros(Pixels*3, frames, data_fmt);
    for n = 1:frames
        imgs(:,n) = fread(fid2, Pixels*3, bit_depth, 0, 'b');
    end
    imgs = [imgs(1:3:end,:); imgs(2:3:end,:); imgs(3:3:end,:)]; % separate color channels
    imgs = reshape(imgs, [Width*Height 3 frames]); % reshape to separate color channels
    N = [Width Height 3 frames];
    imgs = permute(reshape(imgs, N), [2 1 3 4]); % standard reshape and permute
else 
    imgs = zeros(Pixels, frames, data_fmt);
    try
        for n = 1:frames
            imgs(:,n) = fread(fid2, Pixels, bit_depth, 0, 'b');
        end
        N = [Width Height frames];
        imgs = permute(reshape(imgs, N), [2 1 3]);
    catch
        imgs = zeros(Height,Width,frames);
    end
end
fclose(fid2);
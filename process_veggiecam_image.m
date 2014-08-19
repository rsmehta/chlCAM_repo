% process_veggiecam_image.m

% Takes a series of infrared and color images, the first of which are
% assumed to be "base" images.  Produces new images in which infrared and
% color images have been combined such that areas of high chlorophyll are
% very bright, and areas of low chlorophyll are very dark.  All pixel
% values are normalized to the base image to account for changes in light
% when the photographs were taken.

% Inputs:

% Data for images.  Open each image in ImageJ and select Image -> Type
% -> RGB stack.  This will open all your images in Red, Green and Blue
% frames.  Outline the black, grey and white white balance cards in each
% image using the "Polygon Selections" tool.  Select Analyze -> Measure to
% obtain the mean and standard deviation pixel values for each card in each
% image.  Save these data in "WB Data.txt".  Column names in "WB Data.txt"
% are: [Year, Month, Day, Site#, color/black, stdev, color/grey, stdev,
% color/white, stdev, ir/black, stdev, ir/grey, stdev, ir/white, stdev]
% where the year, month and day give the date the images were taken, and
% "color/black" indicates the average pixel value on the black card in a
% particular color image.  

% Original images - images from field work.

% Outputs:

% WB_fit_data.txt - contains the slope and R-squared value for the linear
% fits which standardize light exposure between images.  These linear fits
% are also plotted.

% New images - new images are saved which have been normalized, converted
% to reflectence, and then merged to produce the final chlCAM image.

clear
clc

%% LOAD FILES

drive = 'C:\Users\Diana\Documents\work\Stanford\LimpetCAM\Field\Chlorophyll CAM\Analyzing Images\'; % Ender
% address of file containing images

cd(drive);
% change the working directory to drive

% load the color images
color_images = dir('*_color.jpg'); 
num_color_images = length(color_images);
color_data = cell(1, num_color_images);

for i = 1:num_color_images 
  color_data{i} = imread(color_images(i).name); 
end

% load the ir images
ir_images = dir('*_ir.jpg'); 
num_ir_images = length(ir_images);
ir_data = cell(1, num_ir_images);

for i = 1:num_ir_images 
    ir_data{i} = imread(ir_images(i).name); 
end

% load the wb data
name = [drive 'WB Data.txt'];
wb_data = load(name);

%% PLOT BASE IMAGE VALUES BY SUBJECT IMAGE VALUES; SAVE BEST FIT LINES

% Make color plots

color_fit_data = zeros(num_color_images - 1, 7);
color_fit_data(:,1) = wb_data(2:end,1);
color_fit_data(:,2) = wb_data(2:end,2);
color_fit_data(:,3) = wb_data(2:end,3);
color_fit_data(:,4) = wb_data(2:end,4);

for i = 2:length(wb_data(:,1))
    
    p = polyfit(wb_data(i, [5,7,9]), wb_data(1, [5,7,9]), 1);
    m = p(1);
    b = p(2);
    yfit = polyval(p, wb_data(i, [5,7,9]));
    yresid = wb_data(1, [5,7,9]) - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(wb_data(1, [5,7,9])-1) * var(wb_data(1, [5,7,9])));
    rsq = 1 - SSresid/SStotal;
    
    image_name = sprintf('Site %u', wb_data(i,4));
    image_date = sprintf('%u', wb_data(i, 1:3));

    figure()

    hold on
    scatter(wb_data(i, [5,7,9]), wb_data(1, [5,7,9]), 'r')
    y = m*wb_data(i, [5,7,9]) + b;
    plot(wb_data(i, [5,7,9]), y, 'b')
    figure_title = ['Normalizing to Base Image: ' image_name ', ', ...
        image_date];
    title(figure_title)
    xlabel('Subject Image (color) [Pixel Value]')
    ylabel('Base Image (color) [Pixel Value]')
    text(200, 80, sprintf('m = %1.2f', m)) 
    text(200, 70, sprintf('b = %1.2f', b))
    text(200, 60, sprintf('R^2 = %1.2f', rsq))
    hold off
    
    color_fit_data(i-1, 5) = m;
    color_fit_data(i-1, 6) = b;
    color_fit_data(i-1, 7) = rsq;
    
end

% Make ir plots

ir_fit_data = zeros(num_color_images - 1, 7);
ir_fit_data(:,1) = wb_data(2:end,1);
ir_fit_data(:,2) = wb_data(2:end,2);
ir_fit_data(:,3) = wb_data(2:end,3);
ir_fit_data(:,4) = wb_data(2:end,4);

for i = 2:length(wb_data(:,1))
    
    p = polyfit(wb_data(i, [11,13,15]), wb_data(1, [11,13,15]), 1);
    m = p(1);
    b = p(2);
    yfit = polyval(p, wb_data(i, [11,13,15]));
    yresid = wb_data(1, [11,13,15]) - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(wb_data(1, [11,13,15])-1) * var(wb_data(1, [11,13,15])));
    rsq = 1 - SSresid/SStotal;
    
    image_name = sprintf('Site %u', wb_data(i,4));
    image_date = sprintf('%u', wb_data(i, 1:3));

    figure()

    hold on
    scatter(wb_data(i, [11,13,15]), wb_data(1, [11,13,15]), 'r')
    y = m*wb_data(i, [11,13,15]) + b;
    plot(wb_data(i, [11,13,15]), y, 'b')
    figure_title = ['Normalizing to Base Image: ' image_name ', ', ...
        image_date];
    title(figure_title)
    xlabel('Subject Image (ir) [Pixel Value]')
    ylabel('Base Image (ir) [Pixel Value]')
    text(160, 80, sprintf('m = %1.2f', m)) 
    text(160, 70, sprintf('b = %1.2f', b))
    text(160, 60, sprintf('R^2 = %1.2f', rsq))
    hold off
    
    ir_fit_data(i-1, 5) = m;
    ir_fit_data(i-1, 6) = b;
    ir_fit_data(i-1, 7) = rsq;
    
end

%% STANDARDIZE IMAGES TO BASE IMAGE

output_drive = 'C:\Users\Diana\Documents\work\Stanford\LimpetCAM\Field\Chlorophyll CAM\Analyzing Images\Corrected_Images\';

% Standardize color images

dims = size(color_data{1}(:,:,1));
num_rows = dims(1);
num_cols = dims(2);

std_red_layers = color_data{1}(:,:,1);
% initialize std_red_layers with the red layer from the base color image

%***********************************
% color and ir contain the color image and infrared image
% respectively.  After using imread(), MATLAB contains
% two arrays (class uint 8) that are 1944x2592x3.  The size
% of the images is 1944 by 2592 pixels.  The 3 tells us that
% they are RGB.  Dimension 1 is red, dimension 2 is green
% and dimension 3 is blue.
%***********************************

for k = 2:length(wb_data(:,1))

    red_layer = color_data{k}(:,:,1);
    b = color_fit_data(k-1,6);
    m = color_fit_data(k-1,5);

    new_red_layer = zeros(num_rows, num_cols);

    for i = 1:num_rows

        for j = 1:num_cols        
            new_red_layer(i,j) = b + m*red_layer(i,j);
        end

    end
    
   std_red_layers(:,:,k) = new_red_layer;
    
end

% Standardize ir images

dims = size(ir_data{1}(:,:,1));
num_rows = dims(1);
num_cols = dims(2);

std_green_layers = color_data{1}(:,:,2);
% initialize std_green_layers with the green layer from the base ir image

for k = 2:length(wb_data(:,1))

    green_layer = ir_data{k}(:,:,2);
    b = ir_fit_data(k-1,6);
    m = ir_fit_data(k-1,5);

    new_green_layer = zeros(num_rows, num_cols);

    for i = 1:num_rows

        for j = 1:num_cols        
            new_green_layer(i,j) = b + m*green_layer(i,j);
        end

    end
    
    std_green_layers(:,:,k) = new_green_layer;
    
end

%% CONVERT TO REFLECTANCE

% Convert color images

dims = size(color_data{1}(:,:,1));
num_rows = dims(1);
num_cols = dims(2);

ref_standard = 0.18;

% ref_red_layers = zeros(num_rows, num_cols);

for k = 1:length(wb_data(:,1))
    
    std_red_layer = std_red_layers(:,:,k);
    ref_red_layer = zeros(num_rows, num_cols);

    for i = 1:num_rows

        for j = 1:num_cols        
            ref_red_layer(i,j) = (std_red_layer(i,j)*ref_standard)/wb_data(k,7);
        end

    end
    
    ref_red_layers(:,:,k) = ref_red_layer;
    
end

% Convert ir images

dims = size(ir_data{1}(:,:,1));
num_rows = dims(1);
num_cols = dims(2);

for k = 1:length(wb_data(:,1))
    
    std_green_layer = std_green_layers(:,:,k);
    ref_green_layer = zeros(num_rows, num_cols);

    for i = 1:num_rows

        for j = 1:num_cols        
            ref_green_layer(i,j) = (std_green_layer(i,j)*ref_standard)/wb_data(k,7);
        end

    end
    
    ref_green_layers(:,:,k) = ref_green_layer;
    
end

%% RECOMBINE CHANNELS

for i = 1:length(wb_data(:,1))
    
    color = ref_red_layers(:,:,i);
    ir = ref_green_layers(:,:,i);
    
    chl_images(:,:,i) = ir./color;
    
end

% this operation produces bright pixels where chlorophyll is high

%% AVERAGE PIXELS

for k = 1:length(wb_data(:,1))
    
    current_image = chl_images(:,:,k);
    
    for i = 1:8:1936
        % iterate through each row

        for j = 1:8:2584
            % iterate through each column

            selection = current_image(i:i+7, j:j+7);
            % select an 8x8 pixel square
            average = mean2(selection);
            % average pixel values within that square
            current_image(i:i+7, j:j+7) = average;
            % replace the image values with averaged values

        end

    end

end

%% SAVE AND DISPLAY FINAL IMAGE

for i = 1:length(wb_data(:,1))
    figure()
    imshow(chl_images(:,:,i))
end

%%
veggieCAM_image_adjusted = imadjust(veggieCAM_image);
% adjusts image so that pixel intensities are normalized across 1 to 255

subplot(1,2,1), subimage(veggieCAM_image_adjusted)
subplot(1,2,2), subimage(color)
% displays final image next to original color image

output_name = 'test1-chla';
% enter desired name for output file
output_filename = [drive output_name '.JPG'];
imwrite(veggieCAM_image_adjusted, output_filename)
% write the final image to file

disp('DONE')












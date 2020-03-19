% scGelDockSwarmingQuantification --- Quantification of swarming motility
clear, close all; 

%%  DEFINE THE BASE DIRECTORY AND THE VARIABLES
% % set number
% setNumber = 22;

% what is the background of the imported file? 
% 1 = white/gray
% 2 = black
bgValue = 1;
% in which format is the file?
fileType = 'tif';
% erosion constant, signal amplification
erosionConstant = 87;
maskAmp = 20;   % background value
swarmAmp = 3;
swarmFillstroke = 1000;
plateFillStroke = 5000;
% % set type
% groups = [1 1 11 11 11 11 11 12 12 12 12 12];
% ignore = [1 1 0 0 0 0 0 0 0 0 0 0];
% groups = [groups; ignore]';

% base directory
myfolder = '030520';
basedir = ['/Users/santamag/Desktop/GUILLEM/pics/' myfolder '/'];
cd /Users/santamag/Desktop/GUILLEM/pics
% get plate file names
platePics = dir(fullfile(basedir, '*.tif'));
platePicsName={platePics.name};
noPlates = length(platePicsName);
%%  CALCULATE THROUGH THE IMAGES
swarmData = array2table(zeros(noPlates, 7), ...
                            'VariableNames', {'Area' 'AreaPercentage' 'Circularity' 'Perimeter' 'MaxLength' 'skel' 'ecc'});
swarmNames = array2table(repmat({''}, noPlates, 1), ...
                            'VariableNames', {'Strains'});
for i = 1:noPlates 
% for i = 1
    swarmNames.Strains{i} = platePicsName{i};
    im = imread([basedir platePicsName{i}], fileType);
%     im = im2uint8(im);
%     im= round(im/max(im(:))*(2^12-1));
    
%     size(im)
    %%%%%%%%%% generate mask
    if bgValue == 1
%         im  = im - min(im(:));
%         im = 255-im;
        im = max(im)-im;
    end
%     figure(i); imshow(im*maskAmp);title(platePicsName{i});
    mask = maskAmp * im;
    mask = im2bw(mask, graythresh(mask));   %graythresh() computes a global threshold T from grayscale image I, using Otsu's method;
                                            %im2bw() generate a binary image, with both outline of plates and swarm, plates should be closed circle
    mask = bwareaopen(mask, plateFillStroke);          %bwareaopen() removes all connected components (objects) that have fewer than P pixels from the binary image BW, producing another binary image;
                                            %here we open the binary image and remove area size < 5000
    mask = imfill(mask, 'holes');           %fills holes in the input binary image BW. In this syntax, a hole is a set of background pixels that cannot 
                                            %be reached by filling in the background from the edge of the image
    mask = imerode(mask, strel('octagon', erosionConstant)); %get ride of the outmost layer of plates
%     figure(i + 100); imshow(mask); title(['mask', platePicsName{i}]);     % display the mask in image
    
    %%%%%%%Measure properties of image regions (plate)
    platePixels  = regionprops(mask, 'Area');   

%     'MajorAxisLength'
%     'MinorAxisLength'
%     'MaxFeretProperties'
%     'MinFeretProperties'
%     figure; imshow(plateConverIm.ConvexImage);

%%%% process swarm

    swarm = swarmAmp * im;   
    swarm = im2bw(swarm, graythresh(swarm));    %this step should generate the image closest to the swarming area, but with petridish
    swarm = swarm.* mask;       % get ride of petridish
    swarm = bwareaopen(swarm, swarmFillstroke); 
    swarm = imfill(swarm, 'holes');
    figure(i + 200); imshow(swarm);
    
    [B,L] = bwboundaries(swarm, 'noholes');
    [B2, L2]=bwboundaries(mask, 'noholes');
%     figure (i+800); imshow(label2rgb(L, @jet, [.5 .5 .5]));
%     hold on;
        

        [m, n]  = size(swarm);
        x       = 1:n;
        y       = 1:m;
        [x, y]  = meshgrid(x, y);
        maxX    = max(x(swarm));
        minX    = min(x(swarm));
        maxY = max(y(swarm));
        minY = min(y(swarm));
%         [indexOfMaxX, ~] = find(x(swarm) == max(x(swarm)));
%         [indexOfminX, ~] = find(x(swarm) == min(x(swarm)));
        
        distanceToCenter = sqrt((maxX(1) - minX(1)).^2 ...
                + (maxY(1) - minY(1)).^2); 
            
        yCenter = mean([y(and(x == max(maxX), swarm==1));...
                        y(and(x == min(minX), swarm==1))]);
        xCenter = mean([maxX, minX]);
%         r       = maxX - xCenter;
%         withoutPlate = mask;
%         withoutPlate(distanceToCenter > 0.79*r) = 0;

        
%         xLength = max(x(swarm))...
%             - min(x(swarm));
%         yLength = max(y(swarm))...
%             - min(y(swarm));
%         maxLength = max(xLength, yLength)* 90/(2*r);
    skeltn = bwskel(swarm);
    figure
    imshow(skeltn)
    for k = 1:length(B)
        boundary        = B{k};
        plateboundary =B2{1};
        figure(i+200)
        imshow(swarmAmp * im)
%         title(['swarmAmp-' platePicsName{i}]);
        hold on
        plot (boundary(:, 2),...
            boundary(:, 1),...
            'r-',...
            'LineWidth', 2);        
        hold on;
        plot(plateboundary(:, 2),...
            plateboundary(:, 1),...
            'g-',...
            'LineWidth', 2);        
        hold on;
         plot([maxX...
                minX],...
                [maxY...
                minY],...
                'LineWidth', 1.5,...
                'Color', [1 0 1]);   
        hold off
        title(platePicsName{i}, 'fontSize', 14);
%         stats                  = regionprops (L, 'Area', 'Centroid', 'Circularity', 'Perimeter');
        stats                  = regionprops (L, 'all');
        

        swarmData.skel(i) = sum(skeltn(:));
        swarmData.ecc(i) = stats(k).Eccentricity;
%         threshold              = 0.94;
%         delta_sq               = diff(boundary).^2;
%         perimeter              = sum(sqrt(sum(delta_sq, 2)));
        swarmData.Area(i)       = stats(k).Area;
        swarmData.AreaPercentage(i)        = 100*(stats(k).Area/platePixels.Area);
        swarmData.Circularity(i)= stats(k).Circularity;
        swarmData.Perimeter(i)        = stats(k).Perimeter;
        swarmData.MaxLength(i)        = distanceToCenter;
%         metric                 = 4 * pi * area/perimeter^2;
%         metric_string          = sprintf('%2.2f',metric);
    end
    clear k stats threshold delta_sq perimeter area metric metric_string...
        boundary B L;
end
clear i bgValue;
swarmData = [swarmNames swarmData cell2table(repmat({myfolder}, noPlates,1), 'VariableNames', {'exp'})]
writetable(swarmData, ['swarmData_' myfolder '.xlsx'], 'writeVariableNames', true);
% %% SUPPLEMENT THE DATA WITH THE GROUPING FOR LATER ON
% swarmData = [swarmData groups];
% 
% fid = fopen(strcat(basedir, '../set', num2str(setNumber), '.txt'), 'w');
% fprintf(fid, [repmat('%5.4f ', 1, size(swarmData, 2) - 1), '%5.4f\n'], swarmData');
% fclose(fid);
% 
% clear fid erosionConstant groups maskAmp setNumber swarmAmp swarm mask im...
%     fileType noPlates;
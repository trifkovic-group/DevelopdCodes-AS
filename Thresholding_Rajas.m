%*Code thresholding Rajas' images
%* Original version written by Edna MB.
% Modifications: Allen Shi
%*Last mod.:July 25th;
%*This code has to be run in MATLAB R2018b or newer versions due to the
%imflatfield function

close all hidden
clear 
clc

cd 'F:\Allen\Rajas\PP-PEO\20180731\20180810';

% Files=dir('*s021*ch00.tif'); %Read all of the images that end with ch00.tif
Files=dir('*ch00.tif'); %Read all of the images that end with ch00.tif
ak=1:1:length(Files);%Number of images in file
initial_layer=0;% initial z layer under consideration
final_layer=40;%final z layer under consideration

timestep =30; %timestep under consideration. Should be modified if many time steps are under consideration.

for i=timestep:timestep; %time steps under consideration

time_s=i;%time_step we're looking at

%------------------------------------------------------------------------------------------------------------------------------%
%Used to convert time step in a string compatible with the name of the
%images
m = floor(log10(time_s)); 
D = mod(floor(time_s ./ 10 .^ (m:-1:0)), 10);
time_step=zeros(1,2);

if length(D)==1
    time_step(1,2)=D;
else if length(D)==2
        time_step(1,1:2)=D;
%     else if length(D)==3
%         time_step=D;
%     end
    end
end
time_step=strjoin(string(time_step));
time_step=regexprep(time_step, '\s+', '');

%------------------------------------------------------------------------------------------------------------------------------%

initialname=strcat('50-50.lif_01_t',time_step,'_z0',num2str(initial_layer),'_ch00.tif'); %name of initial z layer at time t
finalname=strcat('50-50.lif_01_t',time_step,'_z',num2str(final_layer),'_ch00.tif'); %name of final z layer at time t
initial_index = find(strcmp({Files.name}, initialname)==1); % determine where in File those layers are located
final_index = find(strcmp({Files.name}, finalname)==1);

aux=0;

for i=initial_index:final_index
    %    
    imch1=imread((Files(i).name));  %read Clay image at z=i
    
    % Cleaning and Filtering
    image2=histeq(histeq(histeq(imch1))); %enhancement of contrast
  
    n = 9; %size of kernel iused in the average box filter
    %rf = 0.95; %thresholding factor
    so = 25; %size threshold to remove of small pores
    soo = 25; %size threshold to remove of small objects
    %-->
    FiltIm = imfilter(image2,ones(n)/n^2,'symmetric'); %Apply an average box filter, symmetric means that 'Input array values outside the bounds of the array are computed by mirror-reflecting the array across the array border.'
    for z = 1:2 %the filter is applied 2 time smore to blur the image and obtain smoother edges   
        FiltIm = imfilter(FiltIm,ones(n)/n^2,'symmetric');
    end
    


%threshold factors based on the image you're evaluating    
    if aux >=15 && aux <=34
        rf = 1.05;
        BW1 = imbinarize(FiltIm,graythresh(FiltIm)*rf);
    elseif aux>=35 & aux <=41
        rf = 0.95;
        BW1 = imbinarize(FiltIm,graythresh(FiltIm)*rf);
        BW2= im2bw(FiltIm,0.55);
    end
    BW = BW1.* BW2;
    %-->
    %small objects and pores
    BW = bwareaopen(~BW,so); %Remove small pores
    BW = bwareaopen(~BW,soo); %Remove small objects
    
    cd 'F:\Allen\Rajas\PP-PEO\20180731\binary'
    %change 'time step' in the imwrite for the correct numbering on second series
    %num2str(str2num(time_step)+297)
    imwrite(BW,strcat('Binary_step_t',time_step,'z_',num2str(aux),'.tiff')); %write the binary images as tiff
    aux=aux+1;
    
    cd 'F:\Allen\Rajas\PP-PEO\20180731\20180810';
end

end

binaryPath='F:\Allen\Rajas\PP-PEO\20180731\binary\';
rawPath='F:\Allen\Rajas\PP-PEO\20180731\20180810\';
for i=1:41
    raw(:,:,i)=imread(strcat(rawPath,'50-50.lif_01_t',num2str(timestep,'%02u'),'_z',num2str(i-1,'%02u'),'_ch00.tif'));
    %raw1(:,:,i)=raw(20:340,20:340,i);
    bw(:,:,i)=imread(strcat(binaryPath,'Binary_step_t',num2str(timestep,'%02u'),'z_',num2str(i-1),'.tiff'));
    outline(:,:,i)=bwperim(bw(:,:,i));
    outlineImg(:,:,i)=raw(:,:,i);
    outlineImg(outline)=255;
    
end
handle1 = implay(raw); 
handle1.Parent.Position = [300 200 550 550];
handle2 = implay(outlineImg);
handle2.Parent.Position = [900 200 550 550];
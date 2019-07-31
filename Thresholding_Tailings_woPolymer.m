% Code for Rahul's images
%Last modification: June 25th-2019
%Code by; Allen Shi

close all hidden
clear all
clc
warning off


cd 'F:\Allen\Rahul_20190619\optimal1' %Folder where tif images are located
Files=dir('*ch01.tif'); %Read names and locations of all of the images that end with ch01.tif (Usually Clay channel) 
Files2=dir('*ch00.tif');%Read names and locations of all of the images that end with ch00.tif (Usually Bitumen channel)

Z=length(Files); %total numbe rof layers

for i=1:Z %binarization of each layer
imch1=imread(strcat(Files(i).name));  %read Clay image at z=i
imch2=imread(strcat(Files2(i).name)); %read Bitumen image at z=i
    
%% Bitumen treatment
%-->
  n=4;
  FiltIm = imfilter(imch2,ones(n)/n^2,'symmetric'); %Apply an average box filter, symmetric means that 'Input array values outside the bounds of the array are computed by mirror-reflecting the array across the array border.'
    for k = 1:2 %the filter is applied 2 time smore to blur the image and obtain smoother edges   
        FiltIm = imfilter(FiltIm,ones(n)/n^2,'symmetric');
    end    
    %-->
   n2=0.25;
    BW= im2bw(FiltIm,n2); %As you increase the numerical value, only brighter values are considered.
   %--> 
   n3=2;
    bw_bitumen(:,:,i)=bwareaopen(BW,n3);
     
    %% This part separates touching objects
%     D = -bwdist(~BW);
%     imshow(D,[])
%     
%     mask = imextendedmin(D,2);
%     D2 = imimposemin(D,mask);
%     Ld2 = watershed(D2);
%     bw3 = BW;
%     bw3(Ld2 == 0) = 0;    

    %% Clay
 
    % Cleaning and Filtering
    image22=imfilter(imch1,ones(3)/3^2); % Averaging filter. Other filter that can be used is the median filter 
    image22=imfilter(image22,ones(3)/3^2);
    image22_i(:,:,i)=image22;
    %-----------------------------------------%-------------------------------------------------%
    se = strel('disk',5);
    Ie = imerode(image22,se); %Erosion of image 
    Iobr = imreconstruct(Ie,image22); %Used to identify high-intensity objects in the mask 
    % External edges. Thresholding after opening-closing by reconstruction
    Iobrd = imdilate(Iobr,se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    %Threshholding
    %-->
    a=0.12;b=0.6;c =0.6;d=0.85;
    if i >=36 && i <=51
         c = 0.52;d=0.8;
% %         bw_clay1(:,:,i)=im2bw(Iobrcbr,a);
%     elseif i ==31
%         a=0.08; c=0.54;
%     bw_clay1(:,:,i)=imbinarize(Iobrcbr,graythresh(Iobrcbr)*c);
    end
    bw_clay1(:,:,i)=imbinarize(Iobrcbr,graythresh(Iobrcbr)*c);
    %bw_clay1(:,:,i)=im2bw(Iobrcbr,a);%% as you increase the second value [0 1], only the brightest pixels are considered 
    %bw_clay2(:,:,i)=imbinarize(Iobrcbr,'adaptive','ForegroundPolarity','bright','Sensitivity',b);%As you reduce the 'b' value, only brighter values are considered.
    bw_clay2(:,:,i)=imbinarize(Iobrcbr,adaptthresh(Iobrcbr,d,'ForegroundPolarity','bright'));
    bw_clay(:,:,i)= bw_clay1(:,:,i).*bw_clay2(:,:,i);
    
    %-----------------------------------------%-------------------------------------------------%
end
cd 'F:\Allen\Rahul_20190619'
[bitumensq]=Bitumen_squish(Z,bw_bitumen);%Bitumen squishing function

cd 'F:\Allen\Rahul_20190619\optimal1\binary\';

for i =1:Z
    bw_clay(:,:,i) = im2double(imread(strcat('Clay',num2str(i),'.tiff')));
    %bitumensq(:,:,i) = im2double(imread(strcat('Bitumen_sq',num2str(i),'.tiff')));
end

clay_bitsq=bw_clay-bitumensq;

%writes 2d images
for i =1:Z
    water(:,:,i)=imcomplement(clay_bitsq(:,:,i)+bitumensq(:,:,i));
    imwrite(water(:,:,i),strcat('Water',num2str(i),'.tiff'));
    imwrite((bitumensq(:,:,i)),strcat('Bitumen_sq',num2str(i),'.tiff'));
    imwrite((clay_bitsq(:,:,i)),strcat('Clay_bitsq',num2str(i),'.tiff'));
    imwrite((bw_bitumen(:,:,i)),strcat('Bitumen_bw',num2str(i),'.tiff'));
end

binaryPath='F:\Allen\Rahul_20190619\optimal1\binary\';
rawPath='F:\Allen\Rahul_20190619\optimal1\';

%video display
for i=1:Z
    
    raw_bitumen(:,:,i)=imread(strcat(rawPath,Files2(i).name)); 
    outline_bitumen(:,:,i)=bwperim(bw_bitumen(:,:,i)); 
    outlineImg_bitumen(:,:,i)=raw_bitumen(:,:,i);
    outlineImg_bitumen(outline_bitumen)=255;
end

handle2 = implay(outlineImg_bitumen);
handle2.Parent.Position = [625 200 550 550];%[1100 200 550 550]
% handle3 = implay(bw_clay2(:,:,i)); 
% handle3.Parent.Position = [2200 200 550 550];
% handle4 = implay(outlineImg_bitumen,2);
% handle4.Parent.Position = [500 200 550 550];%[2800 200 550 550]
%     outline=bwperim(bw_clay(:,:,end));
%     outlineImg=imch1;
%     outlineImg(outline)=255;
%     figure,imshow(outlineImg)
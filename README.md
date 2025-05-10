# Chromatin-heterogeneity-and-phase-behavior

MATLAB code for nuclear volume analysis





function AreaVolume_Sep(folderin,Areavoxel,z_increment,outputfile_name,select,imagedirection,threshold_correction_factor)

%%Code to calculate the area and volume of the cell, Only one cell in one image
%input variable
%folderin='C:/.....'%without slash at the end %put all the tif file for one type of sample in same folder
%Areavoxel  %unit: um^2
%z_increment %z increment in um
%select
%select==0 for cell,  select==1 for nucleus
%imagedirection should be the same as confocal imaging(used to calculate cell height)
%imagedirection=0 means slice 1-->end  ==  cell bottom --> cell top
%imagedirection=1 means slice 1-->end  ==  cell top --> cell bottom
%top_area_thresh_um2 %area threshold to define top of the cell

top_area_thresh_um2=10
top_area_thresh=top_area_thresh_um2/Areavoxel;

%clear all

cd(folderin)
d = dir( '*.tif' );
filenumber=numel(d);
d.name
filename_list=cell(filenumber,1); %initialization file name list

%loop through all the tif file
for k=1:filenumber %Loop over all images

    filename = d(k).name % get individual tif file
    filename_list(k,1)={filename};  %get file name list for output

    filepath=strcat(folderin,'\',filename); % get tif file path
 
    ImInf = imfinfo(filepath);%imageinfo(filename) creates an Image Information tool containing image metadata from the graphics file filename
    zslices =  max(size(ImInf)); %How many slices are in the image
    imstk=readstk(filepath,1,zslices); %Create a variable where all data is loaded into a matrix
    
    %use max z stack projection to calculate threshold
    maxproject=max(imstk,[],3); %Take maximum along the z-direction

    gry=gry*threshold_correction_factor

    bwim=im2bw(maxproject,gry); %Convert image to black and white. Use graythresh to automatically find threshold for conversion.


   
    % create the figure for saturated region, use color RED
    img_saturate=(maxproject==255);
    RGB = zeros(512, 512, 3);
    RGB(:, :, 1) = img_saturate;

    %plot figure
    figure;
    imshow(RGB);
    hold on
    
    %overlay the RED image with(oringinal image + threshold image)
    img_overlay=max(maxproject,uint8(bwim)*255);  % do overlay of grey image with threshold image
    hhh=imshow(img_overlay);
    set( hhh, 'AlphaData', .65 )  % set transparency to 50%
    saveas( gcf, strcat( d(k).name, '_CodeSep.jpg' ) );
   
    regiontest=regionprops(bwim); %Measure all connected regions.
    [maxArea, index]=max([regiontest.Area]); %Find the maximum area of all the connected regions measured.
   
    grylist(k)=gry;
    CellArea(k)=maxArea*Areavoxel; %Convert number of pixels to um^2.
   
    %calculate the cell volume
    for j=1:zslices
        bwimz(:,:,j)=im2bw(imstk(:,:,j),gry);
    end
    
    %find cell height
    imageinten=sum(sum(bwimz,1),2);
    [max_val,max_loc]=max(imageinten);
    max_loc
    if imagedirection==1
        jj=1;
        for j=max_loc:(-1):1
            if imageinten(:,:,j)<top_area_thresh
                jj=j;
                break;
            end
        end
        CellHeight(k)=(max_loc-jj)*z_increment;
    else
        jj=zslices;
        for j=max_loc:zslices
            if imageinten(:,:,j)<top_area_thresh
                jj=j;
                break;
            end
        end
        CellHeight(k)=(jj-max_loc)*z_increment;
    end
    jj
   
    if select==0
        %calculate the volume above the max intensity plane
        %CellVolume(k)=sum(sum(sum(bwimz(:,:,jj:max_loc))))*Areavoxel*z_increment;
       
        %calculate all the pixels whether above or under the max intensity plane
        CellVolume(k)=sum(sum(sum(bwimz)))*Areavoxel*z_increment;
    else
        CellVolume(k)=sum(sum(sum(bwimz)))*Areavoxel*z_increment;
    end

end
   
writefile=folderin;
writefileout_sum=strcat(writefile,'\',outputfile_name,'_CodeSep.xls');

output={};
if select==0
output(1,1)={'Filename'};
output(2:(1+filenumber),1)=filename_list;
output(1,2)={'cell area'};
output(2:(1+filenumber),2)=num2cell(CellArea);
output(1,3)={'cell volume'};
output(2:(1+filenumber),3)=num2cell(CellVolume);
output(1,4)={'threshold list'};
output(2:(1+filenumber),4)=num2cell(grylist);
output(1,5)={'CellHeight'};
output(2:(1+filenumber),5)=num2cell(CellHeight);
output(1,6)={'cell area mean'};
output(2,6)=num2cell(mean(CellArea));
output(1,7)={'cell area std'};
output(2,7)=num2cell(std(CellArea));
output(1,8)={'cell volume mean'};
output(2,8)=num2cell(mean(CellVolume));
output(1,9)={'cell volume std'};
output(2,9)=num2cell(std(CellVolume));
output(1,10)={'CellHeight mean'};
output(2,10)=num2cell(mean(CellHeight));
output(1,11)={'CellHeight std'};
output(2,11)=num2cell(std(CellHeight));
output(1,12)={'threshold mean'};
output(2,12)=num2cell(mean(grylist));
output(1,13)={'threshold std'};
output(2,13)=num2cell(std(grylist));
else
    output(1,1)={'Filename'};
    output(2:(1+filenumber),1)=filename_list;
    output(1,2)={'nucleus crosssection area'};
    output(2:(1+filenumber),2)=num2cell(CellArea);
    output(1,3)={'nucleus volume'};
    output(2:(1+filenumber),3)=num2cell(CellVolume);
    output(1,4)={'threshold list'};
    output(2:(1+filenumber),4)=num2cell(grylist);
    output(1,5)={'nucleus crosssection area mean'};
    output(2,5)=num2cell(mean(CellArea));
    output(1,6)={'nucleus crosssection area std'};
    output(2,6)=num2cell(std(CellArea));
    output(1,7)={'nucleus volume mean'};
    output(2,7)=num2cell(mean(CellVolume));
    output(1,8)={'nucleus volume std'};
    output(2,8)=num2cell(std(CellVolume));
    output(1,9)={'threshold mean'};
    output(2,9)=num2cell(mean(grylist));
    output(1,10)={'threshold std'};
    output(2,10)=num2cell(std(grylist));
end


   
xlswrite(writefileout_sum,output); %Save areas to excel file.


end





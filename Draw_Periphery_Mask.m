function []=Draw_Periphery_Mask(TotalImages)

%If "D_or_Sigma=1" then it will calculate the subsection results for D otherwise it does the analysis based on Sigma
%TotalImages is the number of Images in your folder.
%This code should read all the images in your folder. It is assumed that
%the structure of the folder is the same as the data produces by
%PWSsoftware (Version Jan 30 2022)

D_or_Sigma=0;
%If "D_or_Sigma=1" then it will calculate the subsection results for D otherwise it does the analysis based on Sigma
Noise=0.03; %PWS Noise level

clc;
close all;  % Close all figures (except those of imtool.)
imtool close all;
clearvars -except Noise D_or_Sigma TotalImages Noise % Make sure the workspace panel is showing.

fontSize=22; 
%

% Check if user has installed the Image Processing Toolbox.
hasIPT = license('test', 'image_toolbox');
if ~hasIPT
    % the Image Processing Toolbox is not installed.
    message = sprintf('Sorry, but you do not seem to have the Image Processing Toolbox.\nDo you want to try to continue anyway?');
    reply = questdlg(message, 'Toolbox missing', 'Yes', 'No', 'Yes');
    if strcmpi(reply, 'No')
        % User said No, so exit.
        return;
    end
end

%[baseFileName,folder] = uigetfile
%fullFileName = fullfile(folder, baseFileName);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Have user browse for a file, from a specified "starting folder."
% For convenience in browsing, set a starting folder from which to browse
startingFolder = '\\backmanlabnas.myqnapcloud.com\Public\Jane_Emily_Vasundhara Share';
if ~exist(startingFolder, 'dir')
    % If that folder doesn't exist, just start in the current folder.
    startingFolder = pwd;
end
% Get the name of the file that the user wants to use.
defaultFileName = fullfile(startingFolder, '*.*');
[baseFileName, folder] = uigetfile(defaultFileName, 'Select a file');
if baseFileName == 0
    % User clicked the Cancel button.
    return;
end
fullFileName = fullfile(folder, baseFileName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=1;








%%%----------------------------------------------------------------------

TotalImages = 10;
for v=1:TotalImages  %
    clearvars -except Cell_Marker_Array Noise adaptiveorFix Nf thickIn system_config nuSys liveCellRI sMax sMin ncCenterLambda sigma D_or_Sigma TotalImages CellMarkerArray MaskRow DatasetRow MaskfullName AnalysisName v fullFileName k folder baseFileName kk
    
    %Search within the folder we selected our analysis file and deterimin the
    %Rms file in that folder.
    %We are using the phares "cell" and "PWS" to automatically determine the location of
    %Rms file. So make sure you have not used these two phrases for naming your
    %folders. Otherwise the program may get confused and do not find the Rms
    %image.
    %----------------
    idx1 = strfind(fullFileName,"Cell");
    idx2 = strfind(fullFileName,"PWS");
    fullFileName(idx1+4:idx2-2)=[];
    CellNumber=num2str(v);
    fullFileName=insertAfter(fullFileName,"Cell",CellNumber);
    
    
    idx1 = strfind(folder,"Cell");
    idx2 = strfind(folder,"PWS");
    folder(idx1+4:idx2-2)=[];
    folder=insertAfter(folder,"Cell",CellNumber);
    %-------------------
    
    %a = h5read('C:\Users\ali_n\Documents\Northwestern PC\Projects\PWS Project\Projects\SOP\95% Ethanol Fixation _Time Variant Analysis\Hella Cell_Practice\liveorig\Cell E\Cell52\PWS\analyses\analysisResults_fixed.h5', '/rms')
    NW_PWSImage = h5read(fullFileName, '/rms');
    noise=Noise;% You have to calculate the noise level by drawing Roi in your background and average it over all your images.
    
    %Your noise level will be higher if you use lower exposure time. For my
    %Images at were 0.037
    true_rms = sqrt(abs(NW_PWSImage.^2 - noise.^2));% Calculate the True Rms by subtracting the noise from Rms signal produced by our software
    NW_PWSImage=true_rms;% We will analyze the True Rms for D calculation. TrueRms will be more different that produced Rms for lower exposure time
    CellImage=NW_PWSImage;
    
    
    %%%%---------------------------------------------------------------------
  
    
    
    
    f=folder;
    f1=folder(1:end-13);
    d=dir(f1);
    d_folder=d.folder;
    dd=d(end);
    kk=1;
    
    %Now we go and find Roi for each analysis. We assume there ara no more than
    %3 analysis
    for ii=1:3 %Find ROI 01,   ROI 02 and ROI 03.
        dd=d(end+1-ii);
        hh=dd.name(end-1:end);
        if hh=='h5'
            maskname{ii}=dd.name;
            maskfolder{ii}=dd.folder;
            maskfullFileName{ii}=fullfile(dd.folder, dd.name);
            
            info=h5info(maskfullFileName{ii});
            
            %info.Groups.Name
            GroupNames=struct2cell(info.Groups);
            [m n]=size(GroupNames);
            %n is the number of Rois we have
            
            for vv=1:n
                
                %GroupName=struct2cell(info.Groups.Name)
                GroupName=char(GroupNames(1,vv));
                DatasetName=char({'/mask/'});
                as=([GroupName DatasetName]);
                msk = h5read(maskfullFileName{ii}, as);
                %m4 = h5read(maskfullFileName{1,1}, '/4/mask/');
                Mask(ii,vv,:,:)=msk;
                
                MaskRow(k,:,:)=msk;
                DatasetRow{k}=GroupName;
                MaskfullName{k}=fullfile(dd.folder, dd.name);
                AnalysisName{k}=f;
                
                Mask_Row(v,kk,:,:)=msk;
                Dataset_Row{v,kk}=GroupName;
                Mask_fullName{v,kk}=fullfile(dd.folder, dd.name);
                Analysis_Name{v,kk}=f;
                
                

                NucMask=double(msk);
                degree=2;
                [Center,Boundry1,Boundry2,Boundry3,Boundry4,Boundry5,Boundry6]=Creat_Periphery_Mask(NucMask,degree);
                
                MaskfullNameChar=char(Mask_fullName{v,kk});
                idx1 = strfind(MaskfullNameChar,".h5");

                Tagged_Text=MaskfullNameChar(idx1-20:idx1-1);
                Periphery_Name='_Cell';
                Lable_Boundry1=[Periphery_Name,num2str(v),'_ROI',num2str(vv),'_Boundry1']
                Lable_Boundry2=[Periphery_Name,num2str(v),'_ROI',num2str(vv),'_Boundry2']
                Lable_Boundry3=[Periphery_Name,num2str(v),'_ROI',num2str(vv),'_Boundry3']
                Lable_Boundry4=[Periphery_Name,num2str(v),'_ROI',num2str(vv),'_Boundry4']
                Lable_Boundry5=[Periphery_Name,num2str(v),'_ROI',num2str(vv),'_Boundry5']
                Lable_Boundry6=[Periphery_Name,num2str(v),'_ROI',num2str(vv),'_Boundry6']
                Lable_Center=[Periphery_Name,num2str(v),'_ROI',num2str(vv),'_Center']

                
                fullFileNameBoundry1=insertAfter(MaskfullNameChar,Tagged_Text,Lable_Boundry1);
                fullFileNameBoundry2=insertAfter(MaskfullNameChar,Tagged_Text,Lable_Boundry2);
                fullFileNameBoundry3=insertAfter(MaskfullNameChar,Tagged_Text,Lable_Boundry3);
                fullFileNameBoundry4=insertAfter(MaskfullNameChar,Tagged_Text,Lable_Boundry4);
                fullFileNameBoundry5=insertAfter(MaskfullNameChar,Tagged_Text,Lable_Boundry5);
                fullFileNameBoundry6=insertAfter(MaskfullNameChar,Tagged_Text,Lable_Boundry6);
                fullFileNameCenter=insertAfter(MaskfullNameChar,Tagged_Text,Lable_Center);
                
                fullFileNameBoundry1(end-2:end)=[];
                fullFileNameBoundry1=[fullFileNameBoundry1, '.tif']
                imwrite(Boundry1,fullFileNameBoundry1,'tiff');
                %imshow(Boundry1)
                %saveas(gcf, fullFileNameBoundry1,'tiffn');%Save Tiff with no compression ;)

                fullFileNameBoundry2(end-2:end)=[];
                fullFileNameBoundry2=[fullFileNameBoundry2, '.tif'];
                imwrite(Boundry2,fullFileNameBoundry2,'tiff');
                %imshow(Boundry2)
                %saveas(gcf, fullFileNameBoundry2,'tiffn');%Save Tiff with no compression ;)

                fullFileNameBoundry3(end-2:end)=[];
                fullFileNameBoundry3=[fullFileNameBoundry3, '.tif'];
                imwrite(Boundry3,fullFileNameBoundry3,'tiff');
                %imshow(Boundry3)
                %saveas(gcf, fullFileNameBoundry3,'tiffn');%Save Tiff with no compression ;)

                fullFileNameBoundry4(end-2:end)=[];
                fullFileNameBoundry4=[fullFileNameBoundry4, '.tif'];
                imwrite(Boundry4,fullFileNameBoundry4,'tiff');
                %imshow(Boundry4)
                %saveas(gcf, fullFileNameBoundry4,'tiffn');%Save Tiff with no compression ;)
                
              
                fullFileNameBoundry5(end-2:end)=[];
                fullFileNameBoundry5=[fullFileNameBoundry5, '.tif'];
                imwrite(Boundry5,fullFileNameBoundry5,'tiff');
                %imshow(Boundry5)
                %saveas(gcf, fullFileNameBoundry5,'tiffn');%Save Tiff with no compression ;)

                fullFileNameBoundry6(end-2:end)=[];
                fullFileNameBoundry6=[fullFileNameBoundry6, '.tif'];
                imwrite(Boundry6,fullFileNameBoundry6,'tiff');
                %imshow(Boundry6)
                %saveas(gcf, fullFileNameBoundry6,'tiffn');%Save Tiff with no compression ;)

                fullFileNameCenter(end-2:end)=[];
                fullFileNameCenter=[fullFileNameCenter, '.tif'];
                imwrite(Center,fullFileNameCenter,'tiff');
                %imshow(Center)
                %saveas(gcf, fullFileNameCenter,'tiffn');%Save Tiff with no compression ;)

                k=k+1;
                kk=kk+1;
                
            end
            
            
        end
    end
end




%Ali Daneshkhah, PhD
%Research Associate Fellow
%Backman Photonics Laboratory
%Biomedical Engineering Department
%McCormick School of Engineering & Applied Sciences
%Northwestern University
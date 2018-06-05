%cript for converting meteorological data into desired format 
%Not a generalized scripts, only to convert "SondeKEF_2010_04_14_00Z.txt"
%to "SondeKEF_2010_04_14_00Z_converted.txt"
%For any new data, users have to do converting by themselfve. 
%Should be simply to do either manually or making use of this script.
clc
clear
close all

%Open and read file 
filename='/gpfs/scratch/zhixuanc/IceApp/Plume_SPH_Model/run/Meteo_data/SondeKEF_2010_04_14_00Z.txt';
fileID = fopen(filename);
data1 = textscan(fileID, '%d%s%s%d%f%d%f%f%d%f', 2321, 'headerLines', 1);
fclose(fileID);

C=zeros(length(data1{1}), 7);
% In Plume-SPH, it is assumed that the meteorological data are store in the
% following order:
% 0 height in km
% 1 density  (kg/m^3)
% 2 pressure (Bar = 100 Pa)
% 3 temperature (K)
% 4 specific humidity
% 5 wind velocity West->East
% 6 wind velocity North->South
C(:,1)=double(data1{1, 4})./1000.0; %convert m to km
C(:,2)=data1{1, 10};
C(:,3)=data1{1, 8};
C(:,4)=data1{1, 5}+273.15;

% The left was store in order:  "D"     "F"  "RH" 
C(:,5)=data1{1, 6};
C(:,6)=data1{1, 7};
C(:,7)=data1{1, 9};

%Now write data into a new text file: 
savefile='/gpfs/scratch/zhixuanc/IceApp/Plume_SPH_Model/run/Meteo_data/SondeKEF_2010_04_14_00Z_converted.txt';
fileID = fopen(savefile,'w');
fprintf(fileID,'%14.4f %18.4f %18.4f %14.4f %14.0f %14.4f %14.0f \r\n', C');
fclose(fileID);


% We only need data from 0 ~ 7000 m  ---> do not need to save all data
C=zeros(2167, 7);
% In Plume-SPH, it is assumed that the meteorological data are store in the
% following order:
% 0 height km
% 1 density  (kg/m^3)
% 2 pressure (Bar = 100 Pa)
% 3 temperature (K)
% 4 specific humidity
% 5 wind velocity West->East
% 6 wind velocity North->South
C(1:2167,1)=double(data1{1, 4}(134:2300))./1000.0; % convert to km
C(1:2167,2)=data1{1, 10}(134:2300);
C(1:2167,3)=data1{1, 8}(134:2300);
C(1:2167,4)=data1{1, 5}(134:2300)+273.15;

% The left was store in order:  "D"     "F"  "RH" 
C(1:2167,5)=data1{1, 6}(134:2300);
C(1:2167,6)=data1{1, 7}(134:2300);
C(1:2167,7)=data1{1, 9}(134:2300);

%Now write data into a new text file: 
savefile='/gpfs/scratch/zhixuanc/IceApp/Plume_SPH_Model/run/Meteo_data/SondeKEF_2010_04_14_00Z_converted.txt';
fileID = fopen(savefile,'w');
fprintf(fileID,'%4.4f %6.4f %8.4f %8.4f %3.0f %7.4f %2.0f \r\n', C');
fclose(fileID);

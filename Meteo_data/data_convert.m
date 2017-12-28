%cript for converting meteorological data into desired format 
%Not a generalized scripts, only to convert "SondeKEF_2010_04_14_00Z.txt"
%to "SondeKEF_2010_04_14_00Z_converted.txt"
%For any new data, users have to do converting by themselfve. 
%Should be simply to do either manually or making use of this script.
clc
clear
close all

%Open and read file 
filename='/rohit1/data/users/zhixuanc/Documents/C_plusplus/Plume_SPH_Model/Meteo_data/SondeKEF_2010_04_14_00Z.txt';
fileID = fopen(filename);
data1 = textscan(fileID, '%d%s%s%d%f%d%f%f%d%f', 2321, 'headerLines', 1);
fclose(fileID);

C=zeros(length(data1{1}), 7);
% In Plume-SPH, it is assumed that the meteorological data are store in the
% following order:
% 0 height
% 1 density  (kg/m^3)
% 2 pressure (Bar = 100 Pa)
% 3 temperature (K)
% 4 specific humidity
% 5 wind velocity West->East
% 6 wind velocity North->South
C(:,1)=data1{1, 4};
C(:,2)=data1{1, 10};
C(:,3)=data1{1, 8};
C(:,4)=data1{1, 5}+273.15;

% The left was store in order:  "D"     "F"  "RH" 
C(:,5)=data1{1, 6};
C(:,6)=data1{1, 7};
C(:,7)=data1{1, 9};

%Now write data into a new text file: 
savefile='/rohit1/data/users/zhixuanc/Documents/C_plusplus/Plume_SPH_Model/Meteo_data/SondeKEF_2010_04_14_00Z_converted.txt';
fileID = fopen(savefile,'w');
fprintf(fileID,'%14.0f %18.4f %18.4f %14.4f %14.0f %14.4f %14.0f \r\n', C');
fclose(fileID);


% We only need data from 0 ~ 7000 m  ---> do not need to save all data
C=zeros(434, 7);
% In Plume-SPH, it is assumed that the meteorological data are store in the
% following order:
% 0 height
% 1 density  (kg/m^3)
% 2 pressure (Bar = 100 Pa)
% 3 temperature (K)
% 4 specific humidity
% 5 wind velocity West->East
% 6 wind velocity North->South
C(1:434,1)=data1{1, 4}(134:567);
C(1:434,2)=data1{1, 10}(134:567);
C(1:434,3)=data1{1, 8}(134:567);
C(1:434,4)=data1{1, 5}(134:567)+273.15;

% The left was store in order:  "D"     "F"  "RH" 
C(1:434,5)=data1{1, 6}(134:567);
C(1:434,6)=data1{1, 7}(134:567);
C(1:434,7)=data1{1, 9}(134:567);

%Now write data into a new text file: 
savefile='/rohit1/data/users/zhixuanc/Documents/C_plusplus/Plume_SPH_Model/Meteo_data/SondeKEF_2010_04_14_00Z_converted.txt';
fileID = fopen(savefile,'w');
fprintf(fileID,'%4.0f %6.4f %8.4f %8.4f %3.0f %7.4f %2.0f \r\n', C');
fclose(fileID);

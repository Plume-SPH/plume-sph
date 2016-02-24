
%function that uesed to visualize results
function vis_results

%Simulation information
nproc = 2;

%get file information
filename = 'pvplot000.h5part';
hdfile_info (filename);

info = h5info(filename);
size = info.Groups(1).Datasets(2);
%read data
for (i=1:nproc)
    
end
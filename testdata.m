
%Data and reference(background)files should be imported and placed into the
%'alldata.mat' file prior to running fitting script. The format of columns
%should go as follows: 
%wavelength // intensity data (all) // reference data (all) // std error
%-----------//----------------------//----------------------//-----------
%where the total number of columns is up to the user, the main script will
%work out the rest. There should always be an equal number of data and
%reference columns. There should only be a single error column. 
dat_1 = importdata('userdata.csv');
dat_2 = importdata('userdata.csv');
dat_3 = importdata('userdata.csv');
refdat_1 = importdata('ref_userdata.csv');
refdat_2 = importdata('ref_userdata.csv');
refdat_3 = importdata('ref_userdata.csv');
%%
wavelength = dat_1.data(2,2:end);
int_1 = dat_1.data(3,2:end);
int_2 = dat_2.data(3,2:end);
int_3 = dat_3.data(3,2:end);
int_all = vertcat(int_1,int_2,int_3);

ref_1 = refdat_1.data(3,2:end);
ref_2 = refdat_2.data(3,2:end);
ref_3 = refdat_3.data(3,2:end);

stderr = zeros(length(wavelength),1);
stderr(:,1) = std(int_all);
for i=1:length(stderr)
    if stderr(i) == 0
        stderr(i) = mean(stderr);
    end
end

alldata(:,1) = wavelength;
alldata(:,2) = int_1;
alldata(:,3) = int_2;
alldata(:,4) = int_3;
alldata(:,5) = ref_1;
alldata(:,6) = ref_2;
alldata(:,7) = ref_3;
alldata(:,8) = stderr;

save testdata.mat alldata %'testdata.mat' will be the input file for the
% fitting script

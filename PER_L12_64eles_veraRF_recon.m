%% Configure local data paths before running this script.
% Keep private dataset locations outside source control. Set VRUSI_DATA_ROOT
% to a directory containing dataset_01 ... dataset_07, or replace the
% generic entries below with paths on your own system.
data_root = getenv('VRUSI_DATA_ROOT');
if isempty(data_root)
    error('Set the VRUSI_DATA_ROOT environment variable before running.');
end

pathdata = {fullfile(data_root,'dataset_01','rf_processed');...
    fullfile(data_root,'dataset_01');...
    fullfile(data_root,'dataset_02');...
    fullfile(data_root,'dataset_03','rf_processed');...
    fullfile(data_root,'dataset_04','rf_processed');...
    fullfile(data_root,'dataset_04');...
    fullfile(data_root,'dataset_05','rf_processed')};

rfdatap = {'RFdata_dataset_01';...
    'RFdata_dataset_01';...
    'RFdata_dataset_02';...
    'RFdata_dataset_03';...
    'RFdata_dataset_04';...
    'RFdata_dataset_04';...
    'RFdata_dataset_05'};
for temp_datap = 5:5
clc;
% clearvars -EXCEPT rf_all;
clearvars -EXCEPT simData_all rfdatap pathdata temp_datap;
% clear
eval('field_init(0)','1;')

L12_64eles_veraRF_recon_continue1;%param defination and tx , measure range defination
% simData_all = [];

frames_num = 1000; %500 frames

for num_rcv = 5:5
rf_all = [];rf_batch1=[];rf_batch2=[];rf_batch3=[];rf_batch4=[];
Rf_name = [pathdata{temp_datap},rfdatap{temp_datap}];
L12_64eles_veraRF_recon_continue2; %rf fft
L12_64eles_veraRF_recon_continue3; %simdata updata

TWIST_vera_L12_64eles;

[x,z] = meshgrid(mx,mz);% in mm
bm_kk = reshape(abs(x_twist(:,1)),Nmpz,Nmpx);
figure(12)
 bm_kk2 = rescale(bm_kk,0,255);
imagesc(uint16(bm_kk2))
% pcolor(x*1000,z*1000,bm_kk)
colormap(gray)
title('Optimization image')
xlabel('[mm]'), ylabel('[mm]')
shading interp, axis equal ij tight

mkdir([pathdata{temp_datap},'optirecon\opti tau0.02 c1540 pdataregion Ir2wyConv CUSI_Regulations_OnlyL1_1angles_tau0_02']);
temp_bm_kk = zeros(Nmpz,Nmpx,frames_num);
temp_bm_kk_complex = zeros(Nmpz,Nmpx,frames_num);
for frame_temp = 1:frames_num

    bm_kk = reshape(abs(x_twist(:,frame_temp)),Nmpz,Nmpx);
    bm_kk = single(bm_kk);
    temp_bm_kk(:,:,frame_temp)=bm_kk;
    temp_bm_kk_complex(:,:,frame_temp) = reshape(x_twist(:,frame_temp),Nmpz,Nmpx);
end
outpic=[pathdata{temp_datap},'optirecon\opti tau0.02 c1540 pdataregion Ir2wyConv CUSI_Regulations_OnlyL1_1angles_tau0_02\','img_',sprintf('%04d',num_rcv),'.tif'] ;
tifwrite(single(temp_bm_kk),outpic) ;  
outpic2=[pathdata{temp_datap},'optirecon\opti tau0.02 c1540 pdataregion Ir2wyConv CUSI_Regulations_OnlyL1_1angles_tau0_02\','x_twist_',sprintf('%04d',num_rcv),'.mat'] ;
save(outpic2,"x_twist");
end
end








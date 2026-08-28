%% make sure that field_init has been called 

pathdata = {'O:\OPTI\20250217 rabbit kid lim\20250217 opti rabbit kid07 rfsvd0.02 svdnumframes1000\';...
    'O:\OPTI\20250217 rabbit kid lim\20250217 opti rabbit kid07\';...
    'O:\OPTI\20250310 mbs freeflow\';...
    'O:\OPTI\20250311 musclembs flow rfsvd0.02\';...
    'I:\OPTI\20250103 rabbitlowerlimb\20250103 opti 5v2cycs lowMBs rabbit02 rfsvd0.05 svdnumframes1000\';...
    'I:\OPTI\20250103 rabbitlowerlimb\20250103 opti 5v2cycs lowMBs rabbit02\';...
    'I:\OPTI\20250103 rabbitlowerlimb\20250103 opti 5v2cycs highMBs rabbit01 rfsvd0.02 svdnumframes1000\';};



rfdatap={'RFdata_17-February-2025_15-13-50';...
    'RFdata_17-February-2025_15-13-50';...
    'RFdata_10-March-2025_13-58-31';...
    'RFdata_11-March-2025_11-06-44';...
    'RFdata_03-January-2025_11-36-32';...
    'RFdata_03-January-2025_11-36-32';...
    'RFdata_03-January-2025_11-34-19';};
for temp_datap = 5:5
clc;
% clearvars -EXCEPT rf_all;
clearvars -EXCEPT simData_all rfdatap pathdata temp_datap;
% clear
% addpath(genpath('C:\Program Files\Polyspace\R2020a\toolbox\Field_II_windows'));
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








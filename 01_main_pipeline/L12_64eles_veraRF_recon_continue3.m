%% simdata all pixels 5 angels and 64 elements
trans_times = 5;% angels nums
simData_all = [];
tic()
for trans_num_now  =  3:3

simdata_one_angle = [];
for re_ele = 1:N_elements
    temp_path2 = getenv('VRUSI_SIMULATION_ROOT');
    if isempty(temp_path2)
        error('Set the VRUSI_SIMULATION_ROOT environment variable before running.');
    end
    simdata_name = fullfile(temp_path2,[num2str(trans_num_now,'%02d'), ...
        'angle_',num2str(re_ele,'%03d'),'elements.mat']);
    simData2 = load(simdata_name);simData2 = simData2.temp_angle_temp_ele_pulseecho_allregion;
    %% fft p-e r
    simData2_com = zeros(nfft,size(simData2,2));
    for i = 1:Nmpx*Nmpz
        simData2_com(:,i) = fftshift(fft(simData2(:,i),nfft));
    end
    %% norm f-p-e r
    simData2_com_norm = simData2_com;
    simData2_com_norm(isnan(simData2_com_norm) ==1) = 0 ;
    simData2_com_norm2 = simData2_com_norm(w_start_per:w_end_per,:);
    %normalize
     %simData2_com_norm2 = normalize_complex_golabel(simData2_com_norm2);
    simdata_one_angle = cat(1,simdata_one_angle,simData2_com_norm2);

%     fprintf('element %d \n simdata len %d \n rf len %d \n',re_ele,size(simData2_com_norm2,1),size(rf_data_com_norm2,1));
    fprintf('angels %d, element %d\n',trans_num_now,re_ele);
    %% plot and consider Nstart & Nend
%     figure(9)
%     subplot(321)
%     x_range_per = fs/size(simData2_com_norm,1)*(-size(simData2_com_norm,1)/2:size(simData2_com_norm,1)/2-1);
%     plot(x_range_per,abs(simData2_com_norm(:,Nmpz*(Nmpx/2-1)+Nmpx/2)));% f-p-e r
%     subplot(323)
%     x_range = (fs)/size(rf_data_com_norm,1)*(-size(rf_data_com_norm,1)/2:size(rf_data_com_norm,1)/2-1);
%     plot(x_range,abs(rf_data_com_norm(:,24)));% f-rf
%     subplot(325)
%     plot(abs(rf_data_com_norm(:,24)));% f-rf
%     subplot(322)
%     plot(abs(simData2_com_norm2(:,Nmpz*(Nmpx/2-1)+Nmpx/2)));% f-p-e r
%     subplot(324)
%     plot(abs(rf_data_com_norm2(:,24)));% f-rf
end
tic()
simData_all = cat(1,simData_all,simdata_one_angle);
fprintf('cat simdata of %d th angle, using %d s\n',trans_num_now,toc());
end
fprintf('simdata for all  over,used %d s\n',toc());
%%

%% verasonics data
% Rf_name = [data_path_temp(6,:),Rf_nameall(6,:)];

Rf_name2 = [Rf_name,'.mat'];
load(Rf_name2);
temp_trans_rf_len = Receive(1).endSample;
% trans_time_each_angle = 2;temp_trans_rf_len = Receive(1).endSample/trans_time_each_angle;
%% frequence range determine parameters
N_len = 1200;f00 = f0;
nfft = temp_trans_rf_len*2;
x_range_rf = (fs)/nfft*(-nfft/2:nfft/2-1);
x_range_per = (fs)/nfft*(-nfft/2:nfft/2-1);
n_w_mid_rf = find(abs((x_range_rf-f00)) == min(abs((x_range_rf-f00))));
n_w_mid_per = find(abs((x_range_per-f00)) == min(abs((x_range_per-f00))));
if nfft > N_len
    %     n_w_mid = size(simData2_com_norm,1)/2;
    w_start_rf = fix(n_w_mid_rf-N_len/2+1);
    if w_start_rf <= 0
        w_start_rf = 1;
    end
    w_end_rf = fix(w_start_rf + N_len);
    if w_end_rf > nfft
        w_end_rf = nfft;
    end
    
    w_start_per = fix(n_w_mid_per-N_len/2+1);
    if w_start_per <= 0
        w_start_per = 1;
    end
    w_end_per = fix(w_start_per + N_len);
    if w_end_per > nfft
        w_end_per = nfft;
    end
end 
%% rf all frames in 5 angels
rf_cat_batch_num = 1000;
batch_num = max([Receive.framenum])/rf_cat_batch_num

batch_num = 1;
for temp_batch = 1:batch_num
temp_co = ['rf_batch',num2str(temp_batch),'=[];'];
eval(temp_co);
for temp_frame_id = 1:1:rf_cat_batch_num                                                               % control frame-rate down-sampling
frame_id = temp_frame_id + (temp_batch-1)*rf_cat_batch_num;
rf_name = [Rf_name,'-',num2str(num_rcv),'-',num2str(frame_id),'.mat'];
load(rf_name);
rf_data_one_frame_five_angles = [];
for angle_id = 3:3
rf_data = RcvDataFrame(((angle_id-1)*temp_trans_rf_len+1):angle_id*temp_trans_rf_len,65:128);
%% norm RF
v_temp = double(rf_data);
v_temp(isnan(v_temp) ==1) = 0 ;
% rf_data_norm = normc(v_temp);
rf_data_norm = v_temp;
rf_data_norm(isnan(rf_data_norm) ==1) = 0 ;
% rf_data_norm = v_temp;
%% fft RF
rf_data_com = zeros(nfft,N_elements);
for i = 1:size(rf_data_norm,2)
    rf_data_com(:,i) = fftshift(fft(rf_data_norm(:,i),nfft));
end
rf_data_com_norm = rf_data_com;
rf_data_com_norm(isnan(rf_data_com_norm) ==1) = 0 ;
rf_data_com_norm2 = rf_data_com_norm(w_start_rf:w_end_rf,:);
rf_data_one_frame_five_angles = cat(1,reshape(rf_data_one_frame_five_angles,size(rf_data_one_frame_five_angles,1)*size(rf_data_one_frame_five_angles,2),1),...
    reshape(rf_data_com_norm2,size(rf_data_com_norm2,1)*size(rf_data_com_norm2,2),1)); % cat one frame 5 angle's frmaes
fprintf('now batch %d rf is %d frame %d angle\n',temp_batch,frame_id,angle_id);
end
temp_co2 = ['rf_batch',num2str(temp_batch),' = cat(2',',','rf_batch',num2str(temp_batch),',rf_data_one_frame_five_angles);']; % cat 2000 frames rf data, 2000 cols
eval(temp_co2);
end
end
% num is 2000frames batch = 4
rf_all = cat(2,rf_batch1,rf_batch2);











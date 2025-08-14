%% Generate the simulated signal datasets
clc
clear
close all

addpath(genpath(pwd))

%% Physical parameter ranges
min_T1 = 50;  % unit: ms
max_T1 = 3000;
min_T2 = 5;  % unit: ms
max_T2 = 250;
min_B1 = 0.62;  % unit: %
max_B1 = 0.62;

% Sequence parameters
xDoc = load('./seq_para.mat').xDoc;

%% Training dataset
save_path = './Simulated_Signals/train.mat';

min_HR_mean = 20;  % unit: beat/min
max_HR_mean = 120;
min_RR_mean = 60000/max_HR_mean;  % unit: ms
max_RR_mean = 60000/min_HR_mean;
min_RR_var = 0.00;  % unit: RR
max_RR_var = 0.50;

RR_mean_num = 100;
RR_var_num = 50;
sig_num_per_RR = 2500;

sig_record = zeros(RR_mean_num*RR_var_num*sig_num_per_RR,10);
RR_record = zeros(RR_mean_num*RR_var_num*sig_num_per_RR,9);
para_record = zeros(RR_mean_num*RR_var_num*sig_num_per_RR,3);

RR_mean_record = linspace(min_RR_mean,max_RR_mean,RR_mean_num);
RR_var_record = linspace(min_RR_var,max_RR_var,RR_var_num);
RR_counter = 0;
for RR_mean_ind = 1:RR_mean_num

    RR_mean = RR_mean_record(RR_mean_ind);

    for RR_var_ind = 1:RR_var_num

        RR_var = RR_var_record(RR_var_ind);

        RR = RR_var*RR_mean*randn(9,1)+RR_mean;
        RR(RR < 0) = (1-RR_var)*RR_mean;

        new_T1_record = (max_T1-min_T1)*rand(sig_num_per_RR,1)+min_T1;
        new_T2_record = min(new_T1_record,max_T2);
        new_T2_record = (new_T2_record-min_T2).*rand(sig_num_per_RR,1)+min_T2;
        new_B1_record = (max_B1-min_B1)*rand(sig_num_per_RR,1)+min_B1;

        D = generate_multimapping_dictionary_sim(new_T1_record,new_T2_record,new_B1_record,RR,xDoc);

        RR_counter = RR_counter+1;
        sig_record(sig_num_per_RR*(RR_counter-1)+1:sig_num_per_RR*RR_counter,:) = D.';
        RR_record(sig_num_per_RR*(RR_counter-1)+1:sig_num_per_RR*RR_counter,:) = repmat(RR.',[sig_num_per_RR,1]);
        para_record(sig_num_per_RR*(RR_counter-1)+1:sig_num_per_RR*RR_counter,:) = [new_T1_record,new_T2_record,new_B1_record];

    end

end

sig_record = single(sig_record);
RR_record = single(RR_record);
para_record = single(para_record);

save(save_path,'sig_record','RR_record','para_record')

%% Validation dataset
save_path = './Simulated_Signals/valid.mat';

min_HR_mean = 20;  % unit: beat/min
max_HR_mean = 120;
min_RR_mean = 60000/max_HR_mean;  % unit: ms
max_RR_mean = 60000/min_HR_mean;
min_RR_var = 0.05;  % unit: RR
max_RR_var = 0.45;

RR_mean_num = 50;
RR_var_num = 25;
sig_num_per_RR = 500;

sig_record = zeros(RR_mean_num*RR_var_num*sig_num_per_RR,10);
RR_record = zeros(RR_mean_num*RR_var_num*sig_num_per_RR,9);
para_record = zeros(RR_mean_num*RR_var_num*sig_num_per_RR,3);

RR_mean_record = linspace(min_RR_mean,max_RR_mean,RR_mean_num);
RR_var_record = linspace(min_RR_var,max_RR_var,RR_var_num);
RR_counter = 0;
for RR_mean_ind = 1:RR_mean_num

    RR_mean = RR_mean_record(RR_mean_ind);

    for RR_var_ind = 1:RR_var_num

        RR_var = RR_var_record(RR_var_ind);

        RR = RR_var*RR_mean*randn(9,1)+RR_mean;
        RR(RR < 0) = (1-RR_var)*RR_mean;

        new_T1_record = (max_T1-min_T1)*rand(sig_num_per_RR,1)+min_T1;
        new_T2_record = min(new_T1_record,max_T2);
        new_T2_record = (new_T2_record-min_T2).*rand(sig_num_per_RR,1)+min_T2;
        new_B1_record = (max_B1-min_B1)*rand(sig_num_per_RR,1)+min_B1;

        D = generate_multimapping_dictionary_sim(new_T1_record,new_T2_record,new_B1_record,RR,xDoc);

        RR_counter = RR_counter+1;
        sig_record(sig_num_per_RR*(RR_counter-1)+1:sig_num_per_RR*RR_counter,:) = D.';
        RR_record(sig_num_per_RR*(RR_counter-1)+1:sig_num_per_RR*RR_counter,:) = repmat(RR.',[sig_num_per_RR,1]);
        para_record(sig_num_per_RR*(RR_counter-1)+1:sig_num_per_RR*RR_counter,:) = [new_T1_record,new_T2_record,new_B1_record];

    end

end

sig_record = single(sig_record);
RR_record = single(RR_record);
para_record = single(para_record);

save(save_path,'sig_record','RR_record','para_record')

%% Testing dataset
save_path = './Simulated_Signals/test.mat';

min_HR_mean = 20;  % unit: beat/min
max_HR_mean = 120;
min_RR_mean = 60000/max_HR_mean;  % unit: ms
max_RR_mean = 60000/min_HR_mean;
min_RR_var = 0.05;  % unit: RR
max_RR_var = 0.45;

RR_mean_num = 50;
RR_var_num = 25;
sig_num_per_RR = 500;

sig_record = zeros(RR_mean_num*RR_var_num*sig_num_per_RR,10);
RR_record = zeros(RR_mean_num*RR_var_num*sig_num_per_RR,9);
para_record = zeros(RR_mean_num*RR_var_num*sig_num_per_RR,3);

RR_mean_record = linspace(min_RR_mean,max_RR_mean,RR_mean_num);
RR_var_record = linspace(min_RR_var,max_RR_var,RR_var_num);
RR_counter = 0;
for RR_mean_ind = 1:RR_mean_num

    RR_mean = RR_mean_record(RR_mean_ind);

    for RR_var_ind = 1:RR_var_num

        RR_var = RR_var_record(RR_var_ind);

        RR = RR_var*RR_mean*randn(9,1)+RR_mean;
        RR(RR < 0) = (1-RR_var)*RR_mean;

        new_T1_record = (max_T1-min_T1)*rand(sig_num_per_RR,1)+min_T1;
        new_T2_record = min(new_T1_record,max_T2);
        new_T2_record = (new_T2_record-min_T2).*rand(sig_num_per_RR,1)+min_T2;
        new_B1_record = (max_B1-min_B1)*rand(sig_num_per_RR,1)+min_B1;

        D = generate_multimapping_dictionary_sim(new_T1_record,new_T2_record,new_B1_record,RR,xDoc);

        RR_counter = RR_counter+1;
        sig_record(sig_num_per_RR*(RR_counter-1)+1:sig_num_per_RR*RR_counter,:) = D.';
        RR_record(sig_num_per_RR*(RR_counter-1)+1:sig_num_per_RR*RR_counter,:) = repmat(RR.',[sig_num_per_RR,1]);
        para_record(sig_num_per_RR*(RR_counter-1)+1:sig_num_per_RR*RR_counter,:) = [new_T1_record,new_T2_record,new_B1_record];

    end

end

sig_record = single(sig_record);
RR_record = single(RR_record);
para_record = single(para_record);

save(save_path,'sig_record','RR_record','para_record')


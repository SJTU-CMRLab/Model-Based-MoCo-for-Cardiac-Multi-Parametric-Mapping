# %%
import torch
import torch.nn as nn
import os
from scipy.io import loadmat, savemat
from Models import MLP

# Select the visible GPUs
os.environ['CUDA_VISIBLE_DEVICES'] = '0'
print('Let\'s use', torch.cuda.device_count(), 'GPU!')
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

# Define the signal simulation network
sig_dim_num = 10
RR_dim_num = 9
para_dim_num = 2
ss_net_input_dim_num = para_dim_num + RR_dim_num
ss_net_output_dim_num = sig_dim_num
ss_net_width = 200
ss_net = MLP(ss_net_input_dim_num, ss_net_output_dim_num, ss_net_width)  # Signal simulation network
ss_net = nn.DataParallel(ss_net)
ss_net.to(device)

# Generate the parameter combinations
T1_record = range(200, 2201, 10)
T2_record = range(20, 171, 5)
T1_num = len(T1_record)
T2_num = len(T2_record)
para_record = torch.zeros(T1_num * T2_num, para_dim_num)
para_num = 0
for T1 in T1_record:
    for T2 in T2_record:
        if T1 > T2:
            para_record[para_num, 0] = T1
            para_record[para_num, 1] = T2
            para_num = para_num + 1
para_record = para_record[0:para_num, :]

# %% Dictionary generation
if __name__ == '__main__':
    sub_num = 5
    slice_num = 3

    para_load_path = f'./Parameters/ss_net.pth'
    ss_net.load_state_dict(torch.load(para_load_path))
    ss_net.eval()
    with torch.no_grad():
        for sub_ind in range(1, sub_num + 1):
            for slice_ind in range(1, slice_num + 1):
                data_path = f'./data/in_vivo/subject{sub_ind:d}/slice{slice_ind:d}/data.mat'
                data = loadmat(data_path)
                RR = data['RR']

                RR = torch.tensor(RR, dtype=torch.float32)
                RR = torch.permute(RR, [1, 0])
                RR = RR.repeat(para_num, 1)

                time_scaling = 1000
                x = torch.concatenate([para_record, RR], dim=1) / time_scaling
                D = ss_net(x.to(device))

                D = torch.permute(D, [1, 0])
                D_norm = torch.sum(D ** 2, dim=0, keepdim=True) ** 0.5
                D = D / D_norm  # Normalization
                D = D.cpu().numpy()

                new_T1_record = para_record[:, 0:1].numpy()
                new_T2_record = para_record[:, 1:2].numpy()

                save_path = f'./data/in_vivo/subject{sub_ind:d}/slice{slice_ind:d}/dict.mat'
                save_dict = {'D': D, 'new_T1_record': new_T1_record, 'new_T2_record': new_T2_record}
                savemat(save_path, save_dict)

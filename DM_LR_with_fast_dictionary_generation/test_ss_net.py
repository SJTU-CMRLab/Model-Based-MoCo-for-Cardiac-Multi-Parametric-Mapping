# %%
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
from Models import MLP, MSE

# Select the visible GPUs
os.environ['CUDA_VISIBLE_DEVICES'] = '0'
print('Let\'s use', torch.cuda.device_count(), 'GPU!')
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

# Define the network
sig_dim_num = 10
RR_dim_num = 9
para_dim_num = 2
ss_net_input_dim_num = para_dim_num + RR_dim_num
ss_net_output_dim_num = sig_dim_num
ss_net_width = 200
ss_net = MLP(ss_net_input_dim_num, ss_net_output_dim_num, ss_net_width)  # Signal simulation network
ss_net = nn.DataParallel(ss_net)
ss_net.to(device)

# Define the loss function
mse = MSE()  # Mean squared error

# Create the testing dataset
time_scaling = 1000
batch_size = 1024

test_data_path = f'./data/simulation/Simulated_Signals/test.mat'
test_data = loadmat(test_data_path)
test_sig_record = test_data['sig_record']
test_RR_record = test_data['RR_record'] / time_scaling
test_para_record = test_data['para_record'][:, 0:2] / time_scaling
test_sig_record, test_RR_record, test_para_record \
    = map(lambda x: torch.tensor(x, dtype=torch.float32), (test_sig_record, test_RR_record, test_para_record))
test_X = torch.concatenate((test_para_record, test_RR_record), dim=1)
test_Y = test_sig_record
test_ds = TensorDataset(test_X, test_Y)
test_dl = DataLoader(test_ds, batch_size=batch_size, shuffle=False)

# %% Perform testing
if __name__ == '__main__':
    para_load_path = f'./Parameters/ss_net.pth'
    ss_net.load_state_dict(torch.load(para_load_path))
    ss_net.eval()
    with torch.no_grad():
        test_loss = [mse(ss_net(test_xb.to(device)), test_yb.to(device)).item()
                     for test_xb, test_yb in test_dl]
        test_loss = sum(test_loss) / len(test_loss)

    print(f'Test loss = {test_loss ** 0.5:.4f}')

# %% Show some examples in the testing dataset
if __name__ == '__main__':
    test_dl = DataLoader(test_ds, batch_size=4, shuffle=True)
    test_xb, test_yb = next(iter(test_dl))
    with torch.no_grad():
        test_pred = ss_net(test_xb.to(device))
    test_pred = test_pred.cpu().numpy()
    test_xb, test_yb = test_xb.numpy(), test_yb.numpy()
    test_xb[:, 0:2] = test_xb[:, 0:2] * time_scaling

    fig, axs = plt.subplots(2, 2, figsize=(8, 6))
    plt.subplots_adjust(bottom=0.05, top=0.95, left=0.05, right=0.99)

    t = np.arange(1, 11, 1)

    axs[0][0].plot(t, test_pred[0, :], label='pred')
    axs[0][0].plot(t, test_yb[0, :], label='true')
    axs[0][0].set_xlim((0.5, 10.5))
    axs[0][0].set_ylim((0, 0.6))
    axs[0][0].legend(loc='best')
    axs[0][0].set_title(f'T1={test_xb[0][0]:.0f}, T2={test_xb[0][1]:.0f}')

    axs[0][1].plot(t, test_pred[1, :], label='pred')
    axs[0][1].plot(t, test_yb[1, :], label='true')
    axs[0][1].set_xlim((0.5, 10.5))
    axs[0][1].set_ylim((0, 0.6))
    axs[0][1].legend(loc='best')
    axs[0][1].set_title(f'T1={test_xb[1][0]:.0f}, T2={test_xb[1][1]:.0f}')

    axs[1][0].plot(t, test_pred[2, :], label='pred')
    axs[1][0].plot(t, test_yb[2, :], label='true')
    axs[1][0].set_xlim((0.5, 10.5))
    axs[1][0].set_ylim((0, 0.6))
    axs[1][0].legend(loc='best')
    axs[1][0].set_title(f'T1={test_xb[2][0]:.0f}, T2={test_xb[2][1]:.0f}')

    axs[1][1].plot(t, test_pred[3, :], label='pred')
    axs[1][1].plot(t, test_yb[3, :], label='true')
    axs[1][1].set_xlim((0.5, 10.5))
    axs[1][1].set_ylim((0, 0.6))
    axs[1][1].legend(loc='best')
    axs[1][1].set_title(f'T1={test_xb[3][0]:.0f}, T2={test_xb[3][1]:.0f}')

    plt.show()


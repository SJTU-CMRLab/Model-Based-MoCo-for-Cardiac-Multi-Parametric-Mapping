# %%
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader
from torch.utils.tensorboard import SummaryWriter
import os
import shutil
from scipy.io import loadmat
from Models import MLP, MSE

torch.set_num_threads(4)

# Select the visible GPUs
os.environ['CUDA_VISIBLE_DEVICES'] = '0'
print('Let\'s use', torch.cuda.device_count(), 'GPU!')
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

# Define the network
sig_dim_num = 10
RR_dim_num = 9
para_dim_num = 2
input_dim_num = para_dim_num + RR_dim_num
output_dim_num = sig_dim_num
width = 200
ss_net = MLP(input_dim_num, output_dim_num, width)  # Signal simulation network
ss_net = nn.DataParallel(ss_net)
ss_net.to(device)

# Define the loss function
mse = MSE()  # Mean squared error

# Define the optimizer
learning_rate = 1e-3
optimizer = optim.Adam(ss_net.parameters(), lr=learning_rate)
scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=5, gamma=0.5)  # Learning rate decay

# Create the training and validation datasets
time_scaling = 1000
batch_size = 1024

train_data_path = f'./data/simulation/Simulated_Signals/train.mat'
train_data = loadmat(train_data_path)
train_sig_record = train_data['sig_record']
train_RR_record = train_data['RR_record'] / time_scaling
train_para_record = train_data['para_record'][:, 0:2] / time_scaling
train_sig_record, train_RR_record, train_para_record \
    = map(lambda x: torch.tensor(x, dtype=torch.float32), (train_sig_record, train_RR_record, train_para_record))
train_X = torch.concatenate((train_para_record, train_RR_record), dim=1)
train_Y = train_sig_record
train_ds = TensorDataset(train_X, train_Y)
train_dl = DataLoader(train_ds, batch_size=batch_size, shuffle=True)

valid_data_path = f'./data/simulation/Simulated_Signals/valid.mat'
valid_data = loadmat(valid_data_path)
valid_sig_record = valid_data['sig_record']
valid_RR_record = valid_data['RR_record'] / time_scaling
valid_para_record = valid_data['para_record'][:, 0:2] / time_scaling
valid_sig_record, valid_RR_record, valid_para_record \
    = map(lambda x: torch.tensor(x, dtype=torch.float32), (valid_sig_record, valid_RR_record, valid_para_record))
valid_X = torch.concatenate((valid_para_record, valid_RR_record), dim=1)
valid_Y = valid_sig_record
valid_ds = TensorDataset(valid_X, valid_Y)
valid_dl = DataLoader(valid_ds, batch_size=batch_size, shuffle=False)

# %% Perform training
if __name__ == '__main__':
    # Create an empty log directory
    log_dir = './Runs'
    shutil.rmtree(log_dir)
    os.mkdir(log_dir)
    writer = SummaryWriter(log_dir, flush_secs=100)

    para_save_path = f'./Parameters/ss_net.pth'
    epoch_num = 50
    iter_ind = 0
    min_valid_loss = torch.inf
    ss_net.train()
    for epoch_ind in range(epoch_num):

        for train_xb, train_yb in train_dl:

            train_xb, train_yb = train_xb.to(device), train_yb.to(device)

            train_xb_pred = ss_net(train_xb)
            train_xb_loss = mse(train_xb_pred, train_yb)

            train_xb_loss.backward()
            optimizer.step()
            optimizer.zero_grad()

            if iter_ind % 1000 == 0:
                ss_net.eval()
                with torch.no_grad():
                    valid_loss = [mse(ss_net(valid_xb.to(device)), valid_yb.to(device)).item()
                                  for valid_xb, valid_yb in valid_dl]
                    valid_loss = sum(valid_loss) / len(valid_loss)
                ss_net.train()

                if valid_loss < min_valid_loss:
                    torch.save(ss_net.state_dict(), para_save_path)
                    min_valid_loss = valid_loss
                    min_valid_loss_epoch_ind = epoch_ind
                    min_valid_loss_iter_ind = iter_ind

                writer.add_scalars('Root train loss vs. Root valid loss',
                                   {'Train': train_xb_loss ** 0.5, 'Valid': valid_loss ** 0.5}, iter_ind)

                print(f'Epoch {epoch_ind:d}, iter {iter_ind:d}: '
                      f'root train loss = {train_xb_loss ** 0.5:.4f}, root valid loss = {valid_loss ** 0.5:.4f}, '
                      f'learn rate = {scheduler.get_last_lr()[0]:g}')

            iter_ind += 1

        scheduler.step()

    print(f'\nEpoch {min_valid_loss_epoch_ind:d}, iter {min_valid_loss_iter_ind:d}: '
          f'minimal root valid loss = {min_valid_loss ** 0.5:.4f}')

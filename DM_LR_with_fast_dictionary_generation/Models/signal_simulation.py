import torch
import torch.nn as nn


class MLP(nn.Module):
    # Plain MLP, 4 layers
    def __init__(self, input_dim, output_dim, width):
        super().__init__()

        self.mlp_block1 = nn.Sequential(
            nn.Linear(in_features=input_dim, out_features=width),
            nn.BatchNorm1d(num_features=width),
            nn.LeakyReLU()
        )
        self.mlp_block2 = nn.Sequential(
            nn.Linear(in_features=width, out_features=width),
            nn.BatchNorm1d(num_features=width),
            nn.LeakyReLU()
        )
        self.mlp_block3 = nn.Sequential(
            nn.Linear(in_features=width, out_features=width),
            nn.BatchNorm1d(num_features=width),
            nn.LeakyReLU()
        )
        self.mlp_block4 = nn.Sequential(
            nn.Linear(in_features=width, out_features=output_dim)
        )

    def forward(self, x):
        x = self.mlp_block1(x)
        x = self.mlp_block2(x)
        x = self.mlp_block3(x)
        x = self.mlp_block4(x)
        return x

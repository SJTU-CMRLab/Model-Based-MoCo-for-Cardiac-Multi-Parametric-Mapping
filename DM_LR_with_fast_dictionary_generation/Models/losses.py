import torch
import torch.nn as nn
import torch.nn.functional as F
import math


class NMAE(nn.Module):
    # Normalized mean absolute error (2D or 3D)
    def __init__(self):
        super().__init__()

    def forward(self, x, y):
        return torch.mean(torch.abs((x - y) / (y + 1e-5)))


class MSE(nn.Module):
    # Mean squared error (2D or 3D)
    def __init__(self):
        super().__init__()

    def forward(self, x, y):
        return torch.mean((x - y) ** 2)

import torch
import torch.nn as nn


class ResNet(nn.Module):
    def __init__(self,
                 input_channels: int=4,
                 output_channels: int=32,
                 kernel_size: int=11,
                 stride: int=1,
                 padding: int=5,
                 dilation: int=1
                 ):
        super(ResNet, self).__init__()
        self.conv1 = nn.Conv1d(in_channels=output_channels, out_channels=output_channels, kernel_size=kernel_size,
                               stride=stride, padding=padding, dilation=dilation)
        self.relu = nn.ReLU()
        self.BN1 = nn.BatchNorm1d(num_features=output_channels)

    def forward(self, x):
        y = self.conv1(self.relu(self.BN1(x)))
        y = y + x
        return y


class VariantSeqClassify(nn.Module):
    def __init__(self,
                 input_channels: int=4,
                 output_channels: int=32,
                 kernel_size: int=11,
                 stride: int=1,
                 padding: int=5,
                 dilation: int=1,
                 reduction1: int=4096,
                 reduction2: int=2048,
                 reduction3: int=2048
                 ):
        super(VariantSeqClassify, self).__init__()
        self.conv1 = nn.Conv1d(in_channels=input_channels, out_channels=output_channels, kernel_size=1, stride=1,
                               padding=0, dilation=1)
        self.ResBlock1 = ResNet(input_channels=output_channels, output_channels=output_channels, kernel_size=kernel_size,
                               stride=stride, padding=padding, dilation=dilation)
        self.ResBlock2 = ResNet(input_channels=output_channels, output_channels=output_channels, kernel_size=kernel_size,
                                stride=stride, padding=padding, dilation=dilation)
        self.relu = nn.ReLU()
        self.BN1 = nn.BatchNorm1d(num_features=output_channels)
        self.conv3 = nn.Conv1d(in_channels=output_channels, out_channels=3, kernel_size=1, stride=1, padding=0,
                               dilation=1)
        self.pool1 = nn.AdaptiveMaxPool1d(reduction1)
        self.pool2 = nn.AdaptiveAvgPool1d(reduction2)
        self.flatten = nn.Flatten()
        self.pool3 = nn.AdaptiveAvgPool1d(reduction3)
        self.linear1 = nn.Linear(in_features=reduction3, out_features=1024)
        self.BN2 = nn.BatchNorm1d(num_features=1024)
        self.linear2 = nn.Linear(in_features=1024, out_features=256)
        self.BN3 = nn.BatchNorm1d(num_features=256)
        self.linear3 = nn.Linear(in_features=256, out_features=64)
        self.linear4 = nn.Linear(in_features=64, out_features=4)
        self.dropout = nn.Dropout(p=0.2)

    def forward(self, x):
        x = self.conv1(x)
        y = self.ResBlock1(x)
        y = self.ResBlock2(y)
        y = self.relu(self.pool2(self.relu(self.pool1(self.conv3(y)))))
        y = self.relu(self.pool3(self.flatten(y)))
        y = self.dropout(self.relu(self.BN2(self.linear1(y))))
        y = self.dropout(self.relu(self.BN3(self.linear2(y))))
        y = self.linear3(y)
        y = self.linear4(y)
        y = torch.softmax(y, dim=1)
        return y
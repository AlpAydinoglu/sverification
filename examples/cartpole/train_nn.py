###Loads the samples (data_training.mat) and trains a NN with ReLU activations. 
###Converts the NN into the LCP format (Section 3.1, Lemma 2). 
###Saves the LCP parameterization of the NN as LCP_param.mat.

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
import scipy.io

#set the seed for repeatability
torch.manual_seed(2)

#load the data and set up the training samples
inputs = np.load('inputs_ARE.npy');
X = torch.from_numpy(np.load('states_ARE.npy')).float()
Y = torch.from_numpy(inputs.reshape([inputs.shape[0], 1])).float()
num_samples = X.shape[0]
num_train_samples = int(np.ceil(num_samples*0.8))
num_test_samples = num_samples - num_train_samples
Xtrain,Xtest = torch.split(X,[num_train_samples,num_test_samples])
Ytrain,Ytest = torch.split(Y,[num_train_samples,num_test_samples])
trainset = torch.utils.data.TensorDataset(torch.Tensor(Xtrain),torch.Tensor(Ytrain))

##neural network structure
sz = 10; #number of neurons on each layer
## neural net
net = nn.Sequential(
        nn.Linear(4,sz,bias=True),
        nn.ReLU(),
        nn.Linear(sz,sz,bias=True),
        nn.ReLU(),
        nn.Linear(sz,1,bias=False),
        )


## optimization setup
criterion = nn.MSELoss(size_average=None, reduce=None, reduction='mean')
optimizer = optim.Adam(net.parameters(), lr=1e-4) 
verbose = False
train_batch_size = 100
trainloader = torch.utils.data.DataLoader(trainset, batch_size= train_batch_size, shuffle=True, num_workers=0)
epoch = 4000
net.train()

# train
for t in range(epoch):
    for i, (X,Y) in enumerate(trainloader):
        batch_size = X.shape[0]
        out = net(X)
        loss = nn.MSELoss()(out,Y)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
    if(np.mod(t,100)==0):
        print('epoch: ',t,'MSE loss: ',loss.item())
        

#obtain the neural network parameters
params = list(net.parameters())
n_layer = 2
n_input = 1
n_states = 4
n_contact = sz*n_layer

#convert the neural network into the LCP format
F = np.identity(( sz*n_layer ))
F[sz:sz+sz,0:sz] = -params[2].data.numpy()
c = np.zeros((sz*n_layer))
c[0:sz] = -params[1].data.numpy()
c[sz:sz+sz] = -params[3].data.numpy()
D = np.zeros((n_input,sz*n_layer))
D[0,sz:sz+sz] = params[4].data.numpy()
k = 0
E   = np.zeros((sz*n_layer, n_states))
E[0:sz,:] = -params[0].data.numpy()

#save the neural network parameters in the LCP format
mdic = {"Fc": F, "c": c, "D": D, "k": k, "Ec": E }
scipy.io.savemat('LCP_param.mat', mdic)


import numpy as np
import copy
import multiprocessing

import torchvision
import torchvision.transforms as transforms

import os

def threshold_func(d,lamda):
    d = d-lamda
    d[d<0] = 0
    return d

def shrinkage_func(M,lamda):
    try:
        u,d,vh = np.linalg.svd(M,full_matrices=False)
        M_ = np.matmul(np.matmul(u,np.diag(threshold_func(d,lamda))),vh) 
    except:
        try: 
            u,d,vh = np.linalg.svd(M+1e-3,full_matrices=False)
            M_ = np.matmul(np.matmul(u,np.diag(threshold_func(d,lamda))),vh) 
        except:
            M_ = M
    return M_

def truncation_func(M,mu):
    return np.sign(M) * threshold_func(np.abs(M),mu)

def rpca_func(process_idx,ret_que,M_orig,mu=None,lamda=None,itertime=30):
    lens = M_orig.shape[0]
    p = M_orig.shape[1]
    q = M_orig.shape[2]
    ret_L = np.zeros_like(M_orig)
    ret_S = np.zeros_like(M_orig)
    for i in range(lens):
        if (i//50 *50==i):
            print('process',process_idx,'iteration',i)
        for j in range(3):
            # initial
            M = M_orig[i,:,:,j]
            S = np.zeros_like(M)
            Y = np.zeros_like(M)
    
            if mu is None:
                mu = p * q / 4 / np.linalg.norm(M,ord=1)
            if lamda is None:
                lamda = 1/(max(p,q))**0.5
        
            for k in range(itertime):
                L = shrinkage_func(M - S + 1 / mu * Y, 1 / mu)
                S = truncation_func(M - L + 1/mu*Y,lamda/mu)
                Y = Y + mu*(M-L-S)
            ret_S[i,:,:,j] = S
            ret_L[i,:,:,j] = L
    print('finished prceosses',process_idx)
    ret_que.put({'L':ret_L,'S':ret_S,'idx':process_idx})
    
        

def compress_func(M,x_index,y_index,size):
    # M: array batchsize*32*32*3 need to compress
    M_compress = np.zeros_like(M)
    batch_size = M.shape[0]
    for i in range(batch_size):
        for j in range(3):
            M_compress[i,x_index[i,:,j],y_index[i,:,j],j] = copy.deepcopy(M[i,x_index[i,:,j],y_index[i,:,j],j]) + np.random.randn(size)*20
    return np.float32(M_compress)

def Lipschitz_func(x_index,y_index):
    L_mat = np.zeros(shape=(32,32))
    L_mat[x_index,y_index] = 1
    return 2*np.linalg.norm(L_mat,'fro')
    
def inner_pro_func(x_index,y_index,M):
    return M[x_index,y_index]

def grad_func(M0,inner_pro,x_index,y_index):
    grad_mat = np.zeros(shape=(32,32))
    tem_mat = inner_pro - inner_pro_func(x_index,y_index,M0)
    grad_mat[x_index,y_index] = tem_mat
    return 2*grad_mat

def obj_func(M0,inner_pro,x_index,y_index):
    return np.sum((M0[x_index,y_index]-inner_pro)**2)
    
def recovery_func(process_idx,save_que,M_compress,x_index_total,y_index_total,lamda,itertime=10000,tor=10,seed=20230520):
    lens = M_compress.shape[0]
    # sample index i and color index j
    ret_M = np.zeros_like(M_compress)
    for i in range(lens):
        if (i//10 *10 == i):
            print('process',process_idx,'iteration',i)
        for j in range(3):
            # initial
            t = 1
            x_index = x_index_total[i,:,j]
            y_index = y_index_total[i,:,j]
            L = Lipschitz_func(x_index,y_index)
            M0 = M_compress[i,:,:,j]
            M = copy.deepcopy(M0) + 100
            N = copy.deepcopy(M)
            M_old = copy.deepcopy(M)
            obj_value_before = np.inf
            # iteration
            for iter in range(itertime):
                inner_pro = inner_pro_func(x_index,y_index,N)
                grad_N = grad_func(M0,inner_pro,x_index,y_index)
                try:
                    u,d,vh = np.linalg.svd(N - 1/L*grad_N,full_matrices=False)
                    M_ = np.matmul(np.matmul(u,np.diag(threshold_func(d,lamda/L))),vh) 
                except:
                    try: 
                        u,d,vh = np.linalg.svd(N - 1/L*grad_N/10,full_matrices=False)
                        M_ = np.matmul(np.matmul(u,np.diag(threshold_func(d,lamda/L))),vh) 
                    except:
                        M_ = N - 1/L*grad_N/10
                t_ = (1 + (1 + 4*(t**2))**0.5)/2
                N = copy.deepcopy(M + (t-1)/t_*(M_ - M))
                M = copy.deepcopy(M_)
                t = t_
                if (iter//20 * 20 == iter):
                    # print('current iteration time is ',iter)
                    obj_value = obj_func(M0,inner_pro,x_index,y_index)
                    # print('loss',obj_value)
                    if abs(obj_value - obj_value_before) < tor and iter >= 40:
                        break
                    obj_value_before = obj_value
                    M_old = copy.deepcopy(M)
            ret_M[i,:,:,j] = copy.deepcopy(M_old)
    print('finished processes',process_idx)   
    save_que.put({'M':ret_M,'idx':process_idx})
                    
     
def compress_and_recovery(M,x_index,y_index,size,recovery=False,lamda =1000,k=500,l=100):
    # do robust pca first to split lowrank and sparse part
    processes_compress = []
    q_compress = multiprocessing.Queue()
    for dx in range(l):
        process_compress = multiprocessing.Process(target=rpca_func,args = (dx,q_compress,M[(dx*k):(dx*k + k),:,:,:]),kwargs={'mu':0.00005,'lamda':0.13,'itertime':50})
        process_compress.start()
        processes_compress.append(process_compress)
    L = np.ndarray(shape=(l*k,32,32,3))
    S = np.ndarray(shape=(l*k,32,32,3))
    for i in range(l):
        ret = q_compress.get()
        idx_ = ret['idx']
        L_part = ret['L']
        S_part = ret['S']
        print('collect process',idx_)
        L[(idx_*k):(idx_*k+k),:,:,:] = L_part
        S[(idx_*k):(idx_*k+k),:,:,:] = S_part
        
    S[abs(S)<5] = 0
    M_compress = compress_func(L,x_index,y_index,size)                        
    if recovery:
        processes = []
        q = multiprocessing.Queue()
        for idx in range(l):
            process = multiprocessing.Process(target=recovery_func, args=(idx,q,M_compress[(idx*k):(idx*k + k),:,:,:],x_index[(idx*k):(idx*k + k),:,:],y_index[(idx*k):(idx*k + k),:,:]),kwargs={'lamda':lamda,'tor':1})
            process.start()
            processes.append(process)
        #print('joint processes...')
        # Wait for all processes to complete
        #for process in processes:
        #    process.join()
        L_rec = np.ndarray(shape=(l*k,32,32,3))
        for i in range(l):
            ret_ = q.get()
            id = ret_['idx']
            dta = ret_['M']
            print('collect process',id)
            L_rec[(id*k):(id*k+k),:,:,:] = dta
        ret_M = L_rec+ S
        ret_M[ret_M>255]= 255
        ret_M[ret_M<0] = 0
        return np.uint8(ret_M)
    ret_M = M_compress + S
    ret_M[ret_M>255]= 255
    ret_M[ret_M<0] = 0
    return np.uint8(ret_M)

# Data
transform_train = transforms.Compose([
    transforms.RandomCrop(32, padding=4),
    transforms.RandomHorizontalFlip(),
    transforms.ToTensor(),
    transforms.Normalize((0.4914, 0.4822, 0.4465), (0.2023, 0.1994, 0.2010)),
])

transform_test = transforms.Compose([
    transforms.ToTensor(),
    transforms.Normalize((0.4914, 0.4822, 0.4465), (0.2023, 0.1994, 0.2010)),
])

np.random.seed(20230527)
size = 32*32
recovery = True
lamda = 300


x_index_train = np.zeros((50000,size,3),dtype=np.int64)
y_index_train = np.zeros((50000,size,3),dtype=np.int64)
x_index_test = np.zeros((10000,size,3),dtype=np.int64)
y_index_test = np.zeros((10000,size,3),dtype=np.int64)
for i in range(50000):
    for j in range(3):
        index = np.random.choice(32*32,size,replace=False)
        x_index_train[i,:,j] = index//32
        y_index_train[i,:,j] = index - x_index_train[i,:,j]*32
        
for i in range(10000):
    for j in range(3):
        index = np.random.choice(32*32,size,replace=False)
        x_index_test[i,:,j] = index//32
        y_index_test[i,:,j] = index - x_index_test[i,:,j]*32
        


trainset = torchvision.datasets.CIFAR10(
    root='./data', train=True, download=True, transform=transform_train)

train_rec = compress_and_recovery(trainset.data,x_index_train,y_index_train,size,recovery=recovery,lamda=lamda,k=500,l=100)
np.save('train_recovery_'+str(recovery)+'_'+str(size) +'_'+ str(lamda) + '_rpca_.npy',train_rec)

testset = torchvision.datasets.CIFAR10(
    root='./data', train=False, download=True, transform=transform_test)
test_rec = compress_and_recovery(testset.data,x_index_test,y_index_test,size,recovery=recovery,lamda=lamda,k=100,l=100)
np.save('test_recovery_'+str(recovery)+'_'+str(size) +'_'+ str(lamda) + '_rpca_.npy',test_rec)



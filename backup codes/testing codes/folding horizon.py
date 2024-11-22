#folding horizon
#%%simulation integrated script 11.04.2020 backorder is optimized in simulation 
import xlrd 
import numpy as np
import time
import matplotlib.pyplot as plt
from RC_fold import RC_fold
from docplex.mp.model import Model
import copy
rnd = np.random
seed=1
rnd.seed(seed)
# index and data
ww = 2
kk = 3
jj = 1
rr = 1
tt = 5
iti = 100
CAP = np.zeros((jj))
aa= 1.5
fixed_cap=1
for zeta in ([-0.4]):
    q_p=np.zeros((jj,tt))
    for t in range (tt):
             q_p[0][t] = fixed_cap*180


    for j in range(jj):
        q_p[j][0]=q_p[j][0]*0.6
        q_p[j][tt-1]=q_p[j][tt-1]*0.6

    q_w=np.zeros((ww,tt))
    for w in range(ww):
        for t in range(tt): 
            q_w[w][t]=500
# =============================================================================
#     p=np.zeros((jj))
#     #level 1: 24 level2: 4   
#     for j in range (jj):
#         p[j]=np.random.randint(25,70)
# =============================================================================
    p=[15]    
    dis_per_unit= 10
    q=[1,1,2,2,1,1,1,1]
    transperunit = np.zeros((jj,rr))
    c_b=np.zeros((jj))
    h=np.zeros((kk,jj))
    #level 1: 24 level2: 4   
    for j in range (jj):
        c_b[j]=p[j]*0.25

        for r in range (rr):
            transperunit[j][r] = p[j]*0.10*q[r]

    #Demand
    dz=np.zeros((rr,jj,tt))
    for t in range (tt):
        dz [:,:,t]= [15]
    D_W=np.zeros((rr,jj,tt))
    dh=np.zeros((rr,jj,tt))
    for j in range (jj):   
        CAP[j]=sum(sum(dz[r][j][t] for r in range (rr))for t in (0,1))
    for r in range(rr):
        for j in range(jj):
             for t in range(tt): 
                 D_W[r][j][t]=dz[r][j][t]*(1+zeta)                
                 dh[r][j][t]=dz[r][j][t]*zeta
                 
    #Initial Inventory
    fixed_inv=1
    ini_inv=np.zeros((ww,jj,kk))
    ini_inv[0][0][1]=fixed_inv*sum(dz[r][0][0] for r in range (rr))


    #possiblity of opening aplant or warehouse
    o_w=[[0,1,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,1,1]]
    o_p=[[0,1,1],[0,0,1],[0,0,1],[0,0,1],[0,0,1],[0,1,1]]
    k_1=[[1,2],[2],[2],[2],[2],[2],[1,2]]
    #lead time
    q=[1,1,2,2,1,1,1,1]
    #setup
    sigma=[[0,400,500],[250,0,300],[280,260,0]]
    #opening cost of warehouse
    C_O_W=[0,800,1000]
    bigM=sum(sum(sum(D_W[r][j][t] for r in range (rr)) for j in range (jj)) for t in range (tt))*1.5
    total_demand=sum(sum(sum( dz[r][j][t] for r in range (rr)) for j in range (jj)) for t in range (tt))     
    SIMstart = time.time()
    
    from NewRCfunction import RC
    objdeter,Z_P_V_deter,Ydeter,Xdeter,Bdeter,Ideter,setup_det_first,setup_det,DETmodeling,DETtotal,backorder_det,prod_det,holding_det,transport_det,opening_det,number_of_DV_RC_model=RC(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,dz,o_p,o_w,q,sigma,bigM,C_O_W,ini_inv,k_1)
    zaman=open('random data BS.txt','a')
    zaman.write(str('************')+'     '+str(total_demand)+'     '+str(jj)+'     '+str(kk)+'     '+str(tt)+'     '+str(rr)+'     '+str(ww)+'     '+str(zeta)+'     '+str('************'))
    zaman.write('\n')
    zaman.write(str(DETmodeling)+'     '+str(DETtotal)+'     '+str(objdeter)+'     '+str(setup_det_first)+'     '+str(setup_det)+'     '+str(opening_det)+'     '+str(holding_det)+'     '+str(backorder_det)+'     '+str(transport_det)+'     '+str(prod_det))
    zaman.write('\n')
    zaman.close()
    realized_demand=copy.deepcopy(dz)
    fold_X=  np.zeros((rr,jj,tt))
    fold_Y=  np.zeros((jj,tt))
    fold_I=  np.zeros((ww,jj,kk,tt))
    fold_Z=  np.zeros((ww,kk,tt))


    for period in (list(range(0, tt))):
        realized_demand[:,:, period] = copy.deepcopy(D_W[:,:, period])
        objdeter_fold,American_value,disposal_valu,Z_P_V_fold,Yfold,Xfold,Bfold,Ifold,costdispose,costamrican,setup_det_first_fold,setup_fold,foldmodeling,foldtotal,backorder_fold,prod_fold,holding_fold,transport_fold,opening_fold,number_of_DV_RC_model=RC_fold(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,realized_demand,o_p,o_w,q,sigma,bigM,C_O_W,ini_inv,k_1,fold_X,fold_Y, fold_Z,fold_I,period,dis_per_unit)
        fold_Y[:, period]=Yfold[:, period]
        fold_X[:,:, period]=Xfold[:,:,period]
        fold_I[:,:,:, period]=Ifold[:,:,:, period]
        fold_Z[:,:, period]=Z_P_V_fold[:,:, period]
# =============================================================================
#     from detsim import detsim
#     American_value_d,disposal_amount_D,obj_cor_d,costamrican_c_d,costdispose_c_d,costbackorder_c_d,costransport_c_d,costsetupfirst_c_d,costsetupd_c_d,costholding_c_d,costopeningwarehouse_c_d=detsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,realized_demand,o_p,o_w,q,sigma,bigM,C_O_W,Ydeter,prod_det,dis_per_unit,ini_inv)                                                                                                                     
#     
#         
# =============================================================================
        

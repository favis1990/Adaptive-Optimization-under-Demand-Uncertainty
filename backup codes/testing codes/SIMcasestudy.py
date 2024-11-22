#%%simulation integrated script 11.04.2020 backorder is optimized in simulation 
#we do the simulation for casestudy sized instances
# dispose and penalty cost are input   31.05.2023
import xlrd 
import numpy as np
import time
import matplotlib.pyplot as plt
rnd = np.random
seed=1
rnd.seed(seed)
from docplex.mp.model import Model
import copy
# index and data
ww = 2
kk = 3
jj = 6
rr = 8
tt = 11
iti = 100
CAP = np.zeros((jj))
aa= 1.5
fixed_cap=1
for zeta in ([0.1, 0.2, 0.3, 0.4]):
    q_p=np.zeros((jj,tt))
    for t in (0,2,4,6,8,10):
             q_p[0][t] = fixed_cap*180
             q_p[1][t] = fixed_cap*130
             q_p[2][t] = fixed_cap*130
             q_p[3][t] = fixed_cap*130
             q_p[4][t] = fixed_cap*130
             q_p[5][t] = fixed_cap*130
    for t in (1,3,5,7,9):
             q_p[0][t] = fixed_cap*180
             q_p[1][t] = fixed_cap*130
             q_p[2][t] = fixed_cap*130
             q_p[3][t] = fixed_cap*130
             q_p[4][t] = fixed_cap*130
             q_p[5][t] = fixed_cap*130
    for j in range(jj):
        q_p[j][0]=q_p[j][0]*0.6
        q_p[j][tt-1]=q_p[j][tt-1]*0.6
    q_w=np.zeros((ww,tt))
    for w in range(ww):
        for t in range(tt): 
            q_w[w][t]=500
# =============================================================================
#     p = [[ 14, 12, 50, 30 , 45, 12],
#           [ 14*aa, 12*aa, 50*aa, 30*aa , 45*aa, 12*aa],
#         [ 14, 12, 50, 30 , 45, 12],
#           [ 14*aa, 12*aa, 50*aa, 30*aa , 45*aa, 12*aa],
#           [ 14, 12, 50, 30 , 45, 12],
#           [ 14*aa, 12*aa, 50*aa, 30*aa , 45*aa, 12*aa],
#           [ 14, 12, 50, 30 , 45, 12],
#           [ 14*aa, 12*aa, 50*aa, 30*aa , 45*aa, 12*aa],
#           [ 14, 12, 50, 30 , 45, 12],
#           [ 14*aa, 12*aa, 50*aa, 30*aa , 45*aa, 12*aa],
#           [ 14, 12, 50, 30 , 45, 12]]
# =============================================================================
    dis_per_unit= 10
    carpi=10
    p = [ 14, 12, 50, 30 , 45, 12]
    q=[1,1,2,2,1,1,1,1]
    transperunit = np.zeros((jj,rr))
    c_b=np.zeros((jj))
    h=np.zeros((kk,jj))
    #level 1: 24 level2: 4   
    for j in range (jj):
        c_b[j]=p[j]*0.25
        h[1][j]=p[j]*0.15
        h[2][j]=p[j]*0.20
        for r in range (rr):
            transperunit[j][r] = p[j]*0.10*q[r]

    #Demand
    dz=np.zeros((rr,jj,tt))
    for t in range (tt):
        dz [:,:,t]= [[15, 20, 5, 7, 3, 8],
                     [9,  12, 8, 6, 8, 8,],
                     [25, 0,  13, 10, 7, 8],
                     [20, 0,  10, 12, 6, 14],
                     [15, 15, 10, 10, 8, 12],
                     [10, 0,  5, 5, 3, 9],
                     [8,  8,  4, 6, 3, 9],
                     [8,  8,  3, 5, 4, 10]]
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
    ini_inv[1][1][2]=fixed_inv*sum(dz[r][1][0] for r in range (rr))
    ini_inv[1][2][2]=fixed_inv*sum(dz[r][2][0] for r in range (rr))
    ini_inv[1][3][2]=fixed_inv*sum(dz[r][3][0] for r in range (rr))
    ini_inv[1][4][2]=fixed_inv*sum(dz[r][4][0] for r in range (rr))
    ini_inv[0][5][1]=fixed_inv*sum(dz[r][5][0] for r in range (rr))
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
    farzad=3
    #%%roboust counterpart with DW
    objRC,Z_P_V_RC,YRC,XRC,BRC,IRC,setup_RC_first,setup_RC,RCmodeling,RCtotal,backorder_RC,prod_RC,holding_RC,transport_RC,opening_RC,number_of_DV_RC_model=RC(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,D_W,o_p,o_w,q,sigma,bigM,C_O_W,ini_inv,k_1) 
    zaman=open('random data BS.txt','a')
    zaman.write(str(RCmodeling)+'     '+str(RCtotal)+'     '+str(objRC)+'     '+str(setup_RC_first)+'     '+str(setup_RC)+'     '+str(opening_RC)+'     '+str(holding_RC)+'     '+str(backorder_RC)+'     '+str(transport_RC)+'     '+str(prod_RC))
    zaman.write('\n')
    zaman.close()
    openingbinary=open('random data ZPV.txt','a')
    openingbinary.write(str('#############################################################################'))
    openingbinary.write('\n')
    openingbinary.write(str(Z_P_V_RC))
    openingbinary.write('\n')
    openingbinary.close()
    #%% Adjustable robust counterpart
    from NewARCfunction import ARC  
    ARCobj,IARC,alphaARC,ARCprodcost,XARC,betaARC,ARCtranscost,BARC,gammaARC,ARCbackcost,ARCholdcost,ARCmodeling,ARCtotal=ARC(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,D_W,dz,dh,o_p,o_w,q,sigma,bigM,setup_RC,setup_RC_first,opening_RC,Z_P_V_RC,ini_inv,k_1) 
    zaman=open('random data BS.txt','a')
    zaman.write(str(ARCmodeling)+'     '+str(ARCtotal)+'     '+str(ARCobj)+'     '+str(setup_RC_first)+'     '+str(setup_RC)+'     '+str(opening_RC)+'     '+str(ARCholdcost)+'     '+str(ARCbackcost)+'     '+str(ARCtranscost)+'     '+str(ARCprodcost))
    zaman.write('\n')
    zaman.close()
    #%% reoptimization
    from NewREARCfunction import REARC
    (RE_obj,YB,alphaN,RE_prodcost,XB,betaN,RE_transcost,BB,gammaN,RE_backcost,RE_holdcost,RE_modeling,RE_total)=REARC (jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,D_W,dz,dh,o_p,o_w,q,sigma,bigM,setup_RC,setup_RC_first,opening_RC,Z_P_V_RC,ARCobj,ini_inv)                               
    zaman=open('random data BS.txt','a')
    zaman.write(str(RE_modeling)+'     '+str(RE_total)+'     '+str(RE_obj)+'     '+str(setup_RC_first)+'     '+str(setup_RC)+'     '+str(opening_RC)+'     '+str(RE_holdcost)+'     '+str(RE_backcost)+'     '+str(RE_transcost)+'     '+str(RE_prodcost))
    zaman.write('\n')
    zaman.close()
    #%% simulation starts
    benevis=open('SIMiulation case.txt','a')
    benevis.write(str('***************'))
    benevis.write(str(jj)+'     '+str(kk)+'     '+str(tt)+'     '+str(rr)+'     '+str(ww)+'     '+str(zeta))
    benevis.write(str('***************'))
    dh=np.zeros((rr,jj,tt))
    realized_demand=np.zeros((rr,jj,tt))
    worsecase_demand=sum(sum(sum((D_W[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt))
#here we start the simulation , in total we have 3 simulation, uniform, uniform right hand side and triangular    
#first simulation is uniform
    for iteration in range(iti):
        #calculate samplede demand for each iteration (realized_demand)
        #sample normaly    
   
        for r in range(rr):
            for j in range(jj):
                for t in range(tt): 
                    #dh[r][j][t]=dz[r][j][t]*(np.random.randint(-zeta*1000,zeta*1000)/1000) 
                    dh[r][j][t]=dz[r][j][t]*(np.random.randint(-zeta*1000,zeta*1000)/1000) 
                    #dh[r][j][t]=+dz[r][j][t]*0.2
                    realized_demand[r][j][t]=dh[r][j][t]+dz[r][j][t]
   
         
        total_realized_demand=sum(sum(sum( realized_demand[r][j][t] for r in range (rr)) for j in range (jj)) for t in range (tt))         
        #deterministic feasiblity check
        from detsim import detsim
        American_value_d,disposal_amount_D,obj_cor_d,costamrican_c_d,costdispose_c_d,costbackorder_c_d,costransport_c_d,costsetupfirst_c_d,costsetupd_c_d,costholding_c_d,costopeningwarehouse_c_d=detsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,realized_demand,o_p,o_w,q,sigma,bigM,C_O_W,Ydeter,prod_det,dis_per_unit,ini_inv,carpi)
        total_American_value_d=np.sum(np.sum(np.sum(np.sum(American_value_d))))
        total_disposal_amount_D=np.sum(np.sum(np.sum(np.sum(np.sum(disposal_amount_D)))))
        loss_ration_d=total_American_value_d/total_realized_demand
        #RC
        from RCsim import RCsim
        newobj,disposal_amount_rc,costdispose_c_rc,newbackordercost=RCsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,realized_demand,q,objRC,backorder_RC,XRC,dis_per_unit)                                                                                                                             
        #ARC
        from reARCsim import reARCsim
        disposal_value_Arc2,X_samp_Arc,obj_cor_Arc,disp_cost_Arc,prodcost_samp_Arc,dist_samp_Arc,costbackorder_c_Arc,costsetupfirst_c_Arc,costsetupd_c_Arc,costholding_c_Arc,opening_RC=reARCsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,q,realized_demand,setup_RC_first,setup_RC,opening_RC,IARC,alphaARC,XARC,betaARC,dis_per_unit,ini_inv)  
        total_disposal_value_Arc2=np.sum(np.sum(np.sum(np.sum(np.sum(disposal_value_Arc2)))))
        #REARC
        disposal_value_RE_2,X_samp_RE,obj_cor_RE,disp_cost_RE,prodcost_samp_RE,dist_samp_RE,costbackorder_c_RE,costsetupfirst_c_RE,costsetupd_c_RE,costholding_c_RE,opening_RC=reARCsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,q,realized_demand,setup_RC_first,setup_RC,opening_RC,YB,alphaN,XB,betaN,dis_per_unit,ini_inv)
        total_exxtra_value_RE2=np.sum(np.sum(np.sum(np.sum(np.sum(disposal_value_RE_2)))))
        #rollong horizon 
        from RC_fold import RC_fold
        #rolling for nominal
        fold_X=  np.zeros((rr,jj,tt))
        fold_Y=  np.zeros((jj,tt))
        fold_I=  np.zeros((ww,jj,kk,tt))
        fold_Z=  np.zeros((ww,kk,tt))

        nominal_folding_demand= copy.deepcopy(dz)
        for period in (list(range(0, tt))):
            nominal_folding_demand[:,:, period] = copy.deepcopy(realized_demand[:,:, period])
            objdeter_fold_nominal,American_value_nominal,disposal_valu_nominal,Z_P_V_fold,Yfold,Xfold,Bfold,Ifold,costdispose_nom,costamrican_nom,setup_det_first_fold,setup_fold,foldmodeling,foldtotal,backorder_fold,prod_fold,holding_fold,transport_fold,opening_fold,number_of_DV_RC_model=RC_fold(carpi,jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,nominal_folding_demand,o_p,o_w,q,sigma,bigM,C_O_W,ini_inv,k_1,fold_X,fold_Y, fold_Z,fold_I,period,dis_per_unit)
            fold_Y[:, period]=Yfold[:, period]
            fold_X[:,:, period]=Xfold[:,:,period]
            fold_I[:,:,:, period]=Ifold[:,:,:, period]
            fold_Z[:,:, period]=Z_P_V_fold[:,:, period]
        #rolling for robust    
        fold_X=  np.zeros((rr,jj,tt))
        fold_Y=  np.zeros((jj,tt))
        fold_I=  np.zeros((ww,jj,kk,tt))
        fold_Z=  np.zeros((ww,kk,tt))
        robust_folding_demand= copy.deepcopy(D_W)
        for period in (list(range(0, tt))):
            robust_folding_demand[:,:, period] = copy.deepcopy(realized_demand[:,:, period])
            objdeter_fold_rob,American_value_rob,disposal_valu_rob,Z_P_V_fold_r,Yfold_r,Xfold_r,Bfold,Ifold_r,costdispose_rob,costamrican_rob,setup_det_first_fold,setup_fold,foldmodeling,foldtotal,backorder_fold,prod_fold,holding_fold,transport_fold,opening_fold,number_of_DV_RC_model=RC_fold(carpi,jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,robust_folding_demand,o_p,o_w,q,sigma,bigM,C_O_W,ini_inv,k_1,fold_X,fold_Y, fold_Z,fold_I,period,dis_per_unit)
            fold_Y[:, period]=Yfold_r[:, period]
            fold_X[:,:, period]=Xfold_r[:,:,period]
            fold_I[:,:,:, period]=Ifold_r[:,:,:, period]
            fold_Z[:,:, period]=Z_P_V_fold_r[:,:, period]
        #recording the simulation result
        benevis.write('\n')
        benevis.write(str(total_realized_demand)+'  '+str(total_American_value_d)+'  '+str(total_disposal_amount_D)+'  '+str(loss_ration_d)+'  '+str(obj_cor_d)+'  '+str(costamrican_c_d)+'  '+str(costdispose_c_d)+'  '+str(costsetupfirst_c_d)+'  '+str(costsetupd_c_d)+'  '+str(costholding_c_d)+'  '+str(costopeningwarehouse_c_d)+'  '+str(costbackorder_c_d)+'  '+str(costransport_c_d)+'  '+str(prod_det)+'  '
                      +str(newobj)+'  '+str(disposal_amount_rc)+'  '+str(costdispose_c_rc)+'  '+str(newbackordercost)+'  '+str(backorder_RC)+'   '+str(transport_RC)+'  '+str(prod_RC)+'  '
                      +str(obj_cor_Arc)+'   '+str(disp_cost_Arc)+'   '+str(total_disposal_value_Arc2)+'  '+str(costbackorder_c_Arc)+'  '+str(dist_samp_Arc)+'  '+str(prodcost_samp_Arc)+'  '
                      +str(obj_cor_RE)+'   '+str(disp_cost_RE)+'   '+str(total_exxtra_value_RE2)+'  '+str(costbackorder_c_RE)+'  '+str(dist_samp_RE)+'  '+str(prodcost_samp_RE)+'  '
                      +str(objdeter_fold_nominal)+'   '+str(costamrican_nom)+'   '+str(costdispose_nom)+'  '+str(objdeter_fold_rob)+'  '+str(costamrican_rob)+'  '+str(costdispose_rob))
        
    benevis.write('\n')
    benevis.write(str('**********END OF THIS SET***********'))
    benevis.write('\n')
    benevis.write(str('#############################################################################'))
    benevis.write('\n')
    benevis.write('\n')
    benevis.write('\n')
    benevis.close()     
    end= time.time()
    
    
    benevis=open('SIMiulation case.txt','a')
    benevis.write(str('***************'))
    benevis.write(str(jj)+'     '+str(kk)+'     '+str(tt)+'     '+str(rr)+'     '+str(ww)+'     '+str(zeta))
    benevis.write(str('***************'))
    dh=np.zeros((rr,jj,tt))
    realized_demand=np.zeros((rr,jj,tt))
    worsecase_demand=sum(sum(sum((D_W[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt))
#second simulation is uniform  right hand side  
    for iteration in range(iti):
        #calculate samplede demand for each iteration (realized_demand)
        #sample normaly    
        for r in range(rr):
            for j in range(jj):
                for t in range(tt): 
                    dh[r][j][t]=dz[r][j][t]*(np.random.randint(0,zeta*1000)/1000) 
                    #dh[r][j][t]=dz[r][j][t]*(np.random.randint(0,zeta*1000)/1000) 
                    #dh[r][j][t]=+dz[r][j][t]*0.2
                    realized_demand[r][j][t]=dh[r][j][t]+dz[r][j][t]

        total_realized_demand=sum(sum(sum( realized_demand[r][j][t] for r in range (rr)) for j in range (jj)) for t in range (tt))         
        #deterministic feasiblity check
        from detsim import detsim
        American_value_d,disposal_amount_D,obj_cor_d,costamrican_c_d,costdispose_c_d,costbackorder_c_d,costransport_c_d,costsetupfirst_c_d,costsetupd_c_d,costholding_c_d,costopeningwarehouse_c_d=detsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,realized_demand,o_p,o_w,q,sigma,bigM,C_O_W,Ydeter,prod_det,dis_per_unit,ini_inv,carpi) 
        total_American_value_d=np.sum(np.sum(np.sum(np.sum(American_value_d))))
        total_disposal_amount_D=np.sum(np.sum(np.sum(np.sum(np.sum(disposal_amount_D)))))
        loss_ration_d=total_American_value_d/total_realized_demand
        #RC
        from RCsim import RCsim
        newobj,disposal_amount_rc,costdispose_c_rc,newbackordercost=RCsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,realized_demand,q,objRC,backorder_RC,XRC,dis_per_unit)                                                                                                                             
        #ARC
        from reARCsim import reARCsim
        disposal_value_Arc2,X_samp_Arc,obj_cor_Arc,disp_cost_Arc,prodcost_samp_Arc,dist_samp_Arc,costbackorder_c_Arc,costsetupfirst_c_Arc,costsetupd_c_Arc,costholding_c_Arc,opening_RC=reARCsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,q,realized_demand,setup_RC_first,setup_RC,opening_RC,IARC,alphaARC,XARC,betaARC,dis_per_unit,ini_inv)  
        total_disposal_value_Arc2=np.sum(np.sum(np.sum(np.sum(np.sum(disposal_value_Arc2)))))
        #REARC
        disposal_value_RE_2,X_samp_RE,obj_cor_RE,disp_cost_RE,prodcost_samp_RE,dist_samp_RE,costbackorder_c_RE,costsetupfirst_c_RE,costsetupd_c_RE,costholding_c_RE,opening_RC=reARCsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,q,realized_demand,setup_RC_first,setup_RC,opening_RC,YB,alphaN,XB,betaN,dis_per_unit,ini_inv)
        total_exxtra_value_RE2=np.sum(np.sum(np.sum(np.sum(np.sum(disposal_value_RE_2)))))
        #rollong horizon 
        from RC_fold import RC_fold
        #rolling for nominal
        fold_X=  np.zeros((rr,jj,tt))
        fold_Y=  np.zeros((jj,tt))
        fold_I=  np.zeros((ww,jj,kk,tt))
        fold_Z=  np.zeros((ww,kk,tt))

        nominal_folding_demand= copy.deepcopy(dz)
        for period in (list(range(0, tt))):
            nominal_folding_demand[:,:, period] = copy.deepcopy(realized_demand[:,:, period])
            objdeter_fold_nominal,American_value_nominal,disposal_valu_nominal,Z_P_V_fold,Yfold,Xfold,Bfold,Ifold,costdispose_nom,costamrican_nom,setup_det_first_fold,setup_fold,foldmodeling,foldtotal,backorder_fold,prod_fold,holding_fold,transport_fold,opening_fold,number_of_DV_RC_model=RC_fold(carpi,jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,nominal_folding_demand,o_p,o_w,q,sigma,bigM,C_O_W,ini_inv,k_1,fold_X,fold_Y, fold_Z,fold_I,period,dis_per_unit)
            fold_Y[:, period]=Yfold[:, period]
            fold_X[:,:, period]=Xfold[:,:,period]
            fold_I[:,:,:, period]=Ifold[:,:,:, period]
            fold_Z[:,:, period]=Z_P_V_fold[:,:, period]
        #rolling for robust    
        fold_X=  np.zeros((rr,jj,tt))
        fold_Y=  np.zeros((jj,tt))
        fold_I=  np.zeros((ww,jj,kk,tt))
        fold_Z=  np.zeros((ww,kk,tt))

        robust_folding_demand= copy.deepcopy(D_W)
        for period in (list(range(0, tt))):
            robust_folding_demand[:,:, period] = copy.deepcopy(realized_demand[:,:, period])
            objdeter_fold_rob,American_value_rob,disposal_valu_rob,Z_P_V_fold_r,Yfold_r,Xfold_r,Bfold,Ifold_r,costdispose_rob,costamrican_rob,setup_det_first_fold,setup_fold,foldmodeling,foldtotal,backorder_fold,prod_fold,holding_fold,transport_fold,opening_fold,number_of_DV_RC_model=RC_fold(carpi,jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,robust_folding_demand,o_p,o_w,q,sigma,bigM,C_O_W,ini_inv,k_1,fold_X,fold_Y, fold_Z,fold_I,period,dis_per_unit)
            fold_Y[:, period]=Yfold_r[:, period]
            fold_X[:,:, period]=Xfold_r[:,:,period]
            fold_I[:,:,:, period]=Ifold_r[:,:,:, period]
            fold_Z[:,:, period]=Z_P_V_fold_r[:,:, period]
        #recording the simulation result
        benevis.write('\n')
        benevis.write(str(total_realized_demand)+'  '+str(total_American_value_d)+'  '+str(total_disposal_amount_D)+'  '+str(loss_ration_d)+'  '+str(obj_cor_d)+'  '+str(costamrican_c_d)+'  '+str(costdispose_c_d)+'  '+str(costsetupfirst_c_d)+'  '+str(costsetupd_c_d)+'  '+str(costholding_c_d)+'  '+str(costopeningwarehouse_c_d)+'  '+str(costbackorder_c_d)+'  '+str(costransport_c_d)+'  '+str(prod_det)+'  '
                      +str(newobj)+'  '+str(disposal_amount_rc)+'  '+str(costdispose_c_rc)+'  '+str(newbackordercost)+'  '+str(backorder_RC)+'   '+str(transport_RC)+'  '+str(prod_RC)+'  '
                      +str(obj_cor_Arc)+'   '+str(disp_cost_Arc)+'   '+str(total_disposal_value_Arc2)+'  '+str(costbackorder_c_Arc)+'  '+str(dist_samp_Arc)+'  '+str(prodcost_samp_Arc)+'  '
                      +str(obj_cor_RE)+'   '+str(disp_cost_RE)+'   '+str(total_exxtra_value_RE2)+'  '+str(costbackorder_c_RE)+'  '+str(dist_samp_RE)+'  '+str(prodcost_samp_RE)+'  '
                      +str(objdeter_fold_nominal)+'   '+str(costamrican_nom)+'   '+str(costdispose_nom)+'  '+str(objdeter_fold_rob)+'  '+str(costamrican_rob)+'  '+str(costdispose_rob))


        
    benevis.write('\n')
    benevis.write(str('**********END OF THIS SET***********'))
    benevis.write('\n')
    benevis.write(str('#############################################################################'))
    benevis.write('\n')
    benevis.write('\n')
    benevis.write('\n')
    benevis.close()     
    end= time.time()
    
    
    benevis=open('SIMiulation case.txt','a')
    benevis.write(str('***************'))
    benevis.write(str(jj)+'     '+str(kk)+'     '+str(tt)+'     '+str(rr)+'     '+str(ww)+'     '+str(zeta))
    benevis.write(str('***************'))
    dh=np.zeros((rr,jj,tt))
    realized_demand=np.zeros((rr,jj,tt))
    worsecase_demand=sum(sum(sum((D_W[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt))
#third  simulation is trinagular
    for iteration in range(iti):
        #calculate samplede demand for each iteration (realized_demand)


        randnums=(np.random.triangular(-zeta*1000,zeta*500,zeta*1000,rr*jj*tt)/1000)
        hi = plt.hist(randnums,density=True)
        plt.show()
        f=0
        for r in range(rr):
            for j in range(jj):
               for t in range(tt): 
                    dh[r][j][t]=dz[r][j][t]*randnums[f]
                    realized_demand[r][j][t]=dh[r][j][t]+dz[r][j][t]
                    f=f+1  
    
                
        total_realized_demand=sum(sum(sum( realized_demand[r][j][t] for r in range (rr)) for j in range (jj)) for t in range (tt))         
        #deterministic feasiblity check
        from detsim import detsim
        American_value_d,disposal_amount_D,obj_cor_d,costamrican_c_d,costdispose_c_d,costbackorder_c_d,costransport_c_d,costsetupfirst_c_d,costsetupd_c_d,costholding_c_d,costopeningwarehouse_c_d=detsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,realized_demand,o_p,o_w,q,sigma,bigM,C_O_W,Ydeter,prod_det,dis_per_unit,ini_inv,carpi)
        total_American_value_d=np.sum(np.sum(np.sum(np.sum(American_value_d))))
        total_disposal_amount_D=np.sum(np.sum(np.sum(np.sum(np.sum(disposal_amount_D)))))
        loss_ration_d=total_American_value_d/total_realized_demand
        #RC
        from RCsim import RCsim
        newobj,disposal_amount_rc,costdispose_c_rc,newbackordercost=RCsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,realized_demand,q,objRC,backorder_RC,XRC,dis_per_unit)                                                                                                                             
        #ARC
        from reARCsim import reARCsim
        disposal_value_Arc2,X_samp_Arc,obj_cor_Arc,disp_cost_Arc,prodcost_samp_Arc,dist_samp_Arc,costbackorder_c_Arc,costsetupfirst_c_Arc,costsetupd_c_Arc,costholding_c_Arc,opening_RC=reARCsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,q,realized_demand,setup_RC_first,setup_RC,opening_RC,IARC,alphaARC,XARC,betaARC,dis_per_unit,ini_inv)  
        total_disposal_value_Arc2=np.sum(np.sum(np.sum(np.sum(np.sum(disposal_value_Arc2)))))
        #REARC
        disposal_value_RE_2,X_samp_RE,obj_cor_RE,disp_cost_RE,prodcost_samp_RE,dist_samp_RE,costbackorder_c_RE,costsetupfirst_c_RE,costsetupd_c_RE,costholding_c_RE,opening_RC=reARCsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,q,realized_demand,setup_RC_first,setup_RC,opening_RC,YB,alphaN,XB,betaN,dis_per_unit,ini_inv)
        total_exxtra_value_RE2=np.sum(np.sum(np.sum(np.sum(np.sum(disposal_value_RE_2)))))
         #rollong horizon 
        from RC_fold import RC_fold
        #rolling for nominal
        fold_X=  np.zeros((rr,jj,tt))
        fold_Y=  np.zeros((jj,tt))
        fold_I=  np.zeros((ww,jj,kk,tt))
        fold_Z=  np.zeros((ww,kk,tt))
 
        nominal_folding_demand= copy.deepcopy(dz)
        for period in (list(range(0, tt))):
            nominal_folding_demand[:,:, period] = copy.deepcopy(realized_demand[:,:, period])
            objdeter_fold_nominal,American_value_nominal,disposal_valu_nominal,Z_P_V_fold,Yfold,Xfold,Bfold,Ifold,costdispose_nom,costamrican_nom,setup_det_first_fold,setup_fold,foldmodeling,foldtotal,backorder_fold,prod_fold,holding_fold,transport_fold,opening_fold,number_of_DV_RC_model=RC_fold(carpi,jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,nominal_folding_demand,o_p,o_w,q,sigma,bigM,C_O_W,ini_inv,k_1,fold_X,fold_Y, fold_Z,fold_I,period,dis_per_unit)
            fold_Y[:, period]=Yfold[:, period]
            fold_X[:,:, period]=Xfold[:,:,period]
            fold_I[:,:,:, period]=Ifold[:,:,:, period]
            fold_Z[:,:, period]=Z_P_V_fold[:,:, period]
        #rolling for robust    
        fold_X=  np.zeros((rr,jj,tt))
        fold_Y=  np.zeros((jj,tt))
        fold_I=  np.zeros((ww,jj,kk,tt))
        fold_Z=  np.zeros((ww,kk,tt))
 
        robust_folding_demand= copy.deepcopy(D_W)
        for period in (list(range(0, tt))):
            robust_folding_demand[:,:, period] = copy.deepcopy(realized_demand[:,:, period])
            objdeter_fold_rob,American_value_rob,disposal_valu_rob,Z_P_V_fold_r,Yfold_r,Xfold_r,Bfold,Ifold_r,costdispose_rob,costamrican_rob,setup_det_first_fold,setup_fold,foldmodeling,foldtotal,backorder_fold,prod_fold,holding_fold,transport_fold,opening_fold,number_of_DV_RC_model=RC_fold(carpi,jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,robust_folding_demand,o_p,o_w,q,sigma,bigM,C_O_W,ini_inv,k_1,fold_X,fold_Y, fold_Z,fold_I,period,dis_per_unit)
            fold_Y[:, period]=Yfold_r[:, period]
            fold_X[:,:, period]=Xfold_r[:,:,period]
            fold_I[:,:,:, period]=Ifold_r[:,:,:, period]
            fold_Z[:,:, period]=Z_P_V_fold_r[:,:, period]
        #recording the simulation result
        benevis.write('\n')
        benevis.write(str(total_realized_demand)+'  '+str(total_American_value_d)+'  '+str(total_disposal_amount_D)+'  '+str(loss_ration_d)+'  '+str(obj_cor_d)+'  '+str(costamrican_c_d)+'  '+str(costdispose_c_d)+'  '+str(costsetupfirst_c_d)+'  '+str(costsetupd_c_d)+'  '+str(costholding_c_d)+'  '+str(costopeningwarehouse_c_d)+'  '+str(costbackorder_c_d)+'  '+str(costransport_c_d)+'  '+str(prod_det)+'  '
                      +str(newobj)+'  '+str(disposal_amount_rc)+'  '+str(costdispose_c_rc)+'  '+str(newbackordercost)+'  '+str(backorder_RC)+'   '+str(transport_RC)+'  '+str(prod_RC)+'  '
                      +str(obj_cor_Arc)+'   '+str(disp_cost_Arc)+'   '+str(total_disposal_value_Arc2)+'  '+str(costbackorder_c_Arc)+'  '+str(dist_samp_Arc)+'  '+str(prodcost_samp_Arc)+'  '
                      +str(obj_cor_RE)+'   '+str(disp_cost_RE)+'   '+str(total_exxtra_value_RE2)+'  '+str(costbackorder_c_RE)+'  '+str(dist_samp_RE)+'  '+str(prodcost_samp_RE)+'  '
                      +str(objdeter_fold_nominal)+'   '+str(costamrican_nom)+'   '+str(costdispose_nom)+'  '+str(objdeter_fold_rob)+'  '+str(costamrican_rob)+'  '+str(costdispose_rob))


    benevis.write('\n')
    benevis.write(str('**********END OF THIS SET***********'))
    benevis.write('\n')
    benevis.write(str('#############################################################################'))
    benevis.write('\n')
    benevis.write('\n')
    benevis.write('\n')
    benevis.close()     
    end= time.time()
    
    
    
    
    SIMtime=end-SIMstart
    zaman=open('random data BS.txt','a')
    zaman.write(str(Z_P_V_RC))
    zaman.write('\n')
    zaman.write(str('**********END OF THIS SET***********')+'      '+str(SIMtime)+'     '+str('**********END OF THIS SET***********'))
    zaman.write('\n')
    zaman.write(str('###############################################################################################################'))
    zaman.write('\n')
    zaman.write('\n')
    zaman.write('\n')
    zaman.close()
    
    #(w*r*j*k*t)+(w*r*j*k*t*r*j*t)+(w*j*k*t)+(w*j*k*t*r*j*t)+(r*j*t)+(r*j*t*r*j*t)+(r*j*t)+(j*t*r*j*t)+(w*t*r*j*t)+(r*j*t*r*j*t)+(w*k*t*r*j*t)+(j*k*r*j*t)+(j*k*r*j*t)+(w*r*j*k*t*r*j*t)+(r*j*t*r*j*t)+(w*j*k*t*r*j*t)

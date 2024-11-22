#%%simulation integrated script 
## you can run both synthatic and case study
##default is synthatic



import numpy as np
import time
import matplotlib.pyplot as plt
rnd = np.random
seed=1
rnd.seed(seed)
from docplex.mp.model import Model
# index and data
ww = 2
kk = 3
jj = 5
rr = 5
tt = 5
zeta = 0.2
seed = 1

iti = 100
CAP = np.zeros((jj))
aa= 1.5
fixed_cap=2
for zeta in ([0.2]):
    dis_per_unit = 1
    rnd = np.random    
    rnd.seed(seed)

    carpi =1
    CAP = np.zeros((jj))
    #COEF is the coeficent of the cost for sensitivity analysis
    #Demand
    dz=  demand_matrix = [[[np.random.randint(5, 26) for t in range(tt)] for j in range(jj)] for r in range(rr)]
    D_W= np.zeros((rr,jj,tt))
    dh= np.zeros((rr,jj,tt))
    
    CAP = [sum(sum(demand_matrix[r][j][t] for r in range(rr)) for t in (0, 1)) for j in range(jj)]
    D_W = [[[demand_matrix[r][j][t] * (1 + zeta) for t in range(tt)] for j in range(jj)] for r in range(rr)]
    dh = [[[demand_matrix[r][j][t] * zeta for t in range(tt)] for j in range(jj)] for r in range(rr)]

   
    q_p= np.zeros((jj,tt))#capacity of producing product j in period t 
    
    for j in range (jj):
        for t in range (tt):
            q_p[j][t] = np.random.randint(1, 5)*sum (dz[r][j][t]  for r in range (rr))
    q_w= np.zeros((ww,tt))
    for w in range(ww):
        for t in range(tt): 
            q_w[w][t]= np.random.randint(1, 5)*sum (sum (dz[r][j][t] for j in range(jj)) for r in range (rr))
    
   #production cost
    p =  [np.random.randint(20, 51)  for j in range(jj)] 
    
    #lead time
    q =  [np.random.randint(1, 3) for r in range(rr)]
    
    
    transperunit = np.zeros((jj,rr)) #transporation cost
    c_b=np.zeros((jj))# backorder cost
    h=np.zeros((kk,jj)) #holding cost
     
    for j in range (jj):
        c_b[j]= p[j]*np.random.randint(25, 35)/10
        
        for k in range (kk):
            h[k][j]= p[j]*np.random.randint(10, 20)/10

        
        for r in range (rr):
            transperunit[j][r] = p[j]*0.10*q[r]

   
                 
    #possiblity of opening a plant or warehouse
    o_w= [[np.random.randint(0, 2) for k in range(kk)]for j in range(jj)]
    for j in range (jj):   
        try :
            temp = o_w[j].index(1)
        except ValueError:
            o_w[j][1]=1
            
        
    
    
    #Initial Inventory
    ini_inv=np.zeros((ww,jj,kk))
    for j in range (jj):    
        ini_inv[o_w[j].index(1)][j][o_w[j].index(1)]= sum(dz[r][0][0] for r in range (rr))
        
    
    #setup
    sigma = [[np.random.randint(250, 500) for k in range(kk)] for k in range(kk)] 
    #opening cost of warehouse
    C_O_W = np.zeros((jj))
    C_O_W=  [np.random.randint(800, 1000)  for j in range(jj)]
    
    bigM= sum(sum(sum(D_W[r][j][t] for r in range (rr)) for j in range (jj)) for t in range (tt))*1.5
    total_demand= sum(sum(sum( dz[r][j][t] for r in range (rr)) for j in range (jj)) for t in range (tt)) 
    
    from NewRCfunction import RC
    objective_value, Z_P_Values, Ydeter, Xdeter, Bdeter, Ideter, setup_cost_first, setup_cost, DETmodeling, DETtotal, backorder_det, prod_det, holding_det, transport_det, opening_cost, number_of_DV_RC_model = RC(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, D_W, o_w, q, sigma, bigM, C_O_W, ini_inv)
    zaman=open('random data BS.txt','a')
    zaman.write(str('************')+'     '+str(total_demand)+'     '+str(jj)+'     '+str(kk)+'     '+str(tt)+'     '+str(rr)+'     '+str(ww)+'     '+str(zeta)+'     '+str('************'))
    zaman.write('\n')
    zaman.write(str(DETmodeling)+'     '+str(DETtotal)+'     '+str(objective_value)+'     '+str(setup_cost_first)+'     '+str(setup_cost)+'     '+str(opening_cost)+'     '+str(holding_det)+'     '+str(backorder_det)+'     '+str(transport_det)+'     '+str(prod_det))
    zaman.write('\n')
    zaman.close()


    #%%roboust counterpart with DW
    objRC, Z_P_V_RC, YRC, XRC , IRC, IRC, setup_RC_first, setup_RC, RCmodeling,   RCtotal, backorder_RC, prod_RC, holding_RC, transport_RC, opening_RC, number_of_DV_RC_model = RC(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, D_W, o_w, q, sigma, bigM, C_O_W, ini_inv)
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
    from ARPDP1 import ARC
    ARCobj,IARC,alphaARC,ARCprodcost,XARC,betaARC,ARCtranscost,BARC,gammaARC,ARCbackcost,ARCholdcost,ARCmodeling,ARCtotal = ARC(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, D_W, dz, dh, o_w,q, sigma, bigM, opening_cost, setup_cost_first, setup_cost, Z_P_Values, ini_inv)
    

    zaman=open('random data BS.txt','a')
    zaman.write(str(ARCmodeling)+'     '+str(ARCtotal)+'     '+str(ARCobj)+'     '+str(setup_RC_first)+'     '+str(setup_RC)+'     '+str(opening_RC)+'     '+str(ARCholdcost)+'     '+str(ARCbackcost)+'     '+str(ARCtranscost)+'     '+str(ARCprodcost))
    zaman.write('\n')
    zaman.close()
    #%% reoptimization
    from NewREARCfunction import REARC
    RE_obj, YB, alphaN, RE_prodcost, XB, betaN, RE_transcost, BB, gammaN, RE_backcost, RE_holdcost, RE_modeling, RE_total = REARC(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, D_W, dz, dh, o_w, q, sigma, bigM, opening_cost, setup_cost_first, setup_cost, Z_P_Values, objective_value, ini_inv)
    
    zaman=open('random data BS.txt','a')
    zaman.write(str(RE_modeling)+'     '+str(RE_total)+'     '+str(RE_obj)+'     '+str(setup_RC_first)+'     '+str(setup_RC)+'     '+str(opening_RC)+'     '+str(RE_holdcost)+'     '+str(RE_backcost)+'     '+str(RE_transcost)+'     '+str(RE_prodcost))
    zaman.write('\n')
    zaman.close()
    #%% simulation starts
    benevis=open('random data Sim.txt','a')
    benevis.write(str('***************'))
    benevis.write(str(jj)+'     '+str(kk)+'     '+str(tt)+'     '+str(rr)+'     '+str(ww)+'     '+str(zeta))
    benevis.write(str('***************'))
    dh=np.zeros((rr,jj,tt))
    realized_demand=np.zeros((rr,jj,tt))
    worsecase_demand=sum(sum(sum((D_W[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt))
    
    for iteration in range(iti):
        #calculate samplede demand for each iteration (realized_demand)
        #sample normaly    
    # =============================================================================
    #     for r in range(rr):
    #         for j in range(jj):
    #             for t in range(tt): 
    #                 #dh[r][j][t]=dz[r][j][t]*(np.random.randint(-zeta*1000,zeta*1000)/1000) 
    #                 dh[r][j][t]=dz[r][j][t]*(np.random.randint(0,zeta*1000)/1000) 
    #                 #dh[r][j][t]=+dz[r][j][t]*0.2
    #                 realized_demand[r][j][t]=dh[r][j][t]+dz[r][j][t]
    # =============================================================================
        #sample from a triangular 
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
        American_value_d,disposal_amount_D,obj_cor_d,costamrican_c_d,costdispose_c_d,costbackorder_c_d,costransport_c_d,costsetupfirst_c_d,costsetupd_c_d,costholding_c_d,costopeningwarehouse_c_d=detsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,realized_demand,o_w,q,sigma,bigM,C_O_W,Ydeter,prod_det,dis_per_unit,ini_inv, carpi)                                                                                                                     
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
        #recording the simulation result
        benevis.write('\n')
        benevis.write(str(total_realized_demand)+'  '+str(total_American_value_d)+'  '+str(total_disposal_amount_D)+'  '+str(loss_ration_d)+'  '+str(obj_cor_d)+'  '+str(costamrican_c_d)+'  '+str(costdispose_c_d)+'  '+str(costsetupfirst_c_d)+'  '+str(costsetupd_c_d)+'  '+str(costholding_c_d)+'  '+str(costopeningwarehouse_c_d)+'  '+str(costbackorder_c_d)+'  '+str(costransport_c_d)+'  '+str(prod_det)+'  '
                      +str(newobj)+'  '+str(disposal_amount_rc)+'  '+str(costdispose_c_rc)+'  '+str(newbackordercost)+'  '+str(backorder_RC)+'   '+str(transport_RC)+'  '+str(prod_RC)+'  '
                      +str(obj_cor_Arc)+'   '+str(disp_cost_Arc)+'   '+str(total_disposal_value_Arc2)+'  '+str(costbackorder_c_Arc)+'  '+str(dist_samp_Arc)+'  '+str(prodcost_samp_Arc)+'  '
                      +str(obj_cor_RE)+'   '+str(disp_cost_RE)+'   '+str(total_exxtra_value_RE2)+'  '+str(costbackorder_c_RE)+'  '+str(dist_samp_RE)+'  '+str(prodcost_samp_RE))
    benevis.write('\n')
    benevis.write(str('**********END OF THIS SET***********'))
    benevis.write('\n')
    benevis.write(str('#############################################################################'))
    benevis.write('\n')
    benevis.write('\n')
    benevis.write('\n')
    benevis.close()     
    end= time.time()
    
    
    benevis=open('random data Sim.txt','a')
    benevis.write(str('***************'))
    benevis.write(str(jj)+'     '+str(kk)+'     '+str(tt)+'     '+str(rr)+'     '+str(ww)+'     '+str(zeta))
    benevis.write(str('***************'))
    dh=np.zeros((rr,jj,tt))
    realized_demand=np.zeros((rr,jj,tt))
    worsecase_demand=sum(sum(sum((D_W[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt))
    
    for iteration in range(iti):
        #calculate samplede demand for each iteration (realized_demand)
        #sample normaly    
        for r in range(rr):
            for j in range(jj):
                for t in range(tt): 
                    dh[r][j][t]=dz[r][j][t]*(np.random.randint(-zeta*1000,zeta*1000)/1000) 
                    #dh[r][j][t]=dz[r][j][t]*(np.random.randint(0,zeta*1000)/1000) 
                    #dh[r][j][t]=+dz[r][j][t]*0.2
                    realized_demand[r][j][t]=dh[r][j][t]+dz[r][j][t]
        #sample from a triangular 
    # =============================================================================
    #     randnums=(np.random.triangular(-zeta*1000,zeta*500,zeta*1000,rr*jj*tt)/1000)
    #     hi = plt.hist(randnums,density=True)
    #     plt.show()
    #     f=0
    #     for r in range(rr):
    #         for j in range(jj):
    #             for t in range(tt): 
    #                 dh[r][j][t]=dz[r][j][t]*randnums[f]
    #                 realized_demand[r][j][t]=dh[r][j][t]+dz[r][j][t]
    #                 f=f+1                
    # =============================================================================
        total_realized_demand=sum(sum(sum( realized_demand[r][j][t] for r in range (rr)) for j in range (jj)) for t in range (tt))         
        #deterministic feasiblity check
        from detsim import detsim
        American_value_d,disposal_amount_D,obj_cor_d,costamrican_c_d,costdispose_c_d,costbackorder_c_d,costransport_c_d,costsetupfirst_c_d,costsetupd_c_d,costholding_c_d,costopeningwarehouse_c_d=detsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,realized_demand,o_w,q,sigma,bigM,C_O_W,Ydeter,prod_det,dis_per_unit,ini_inv,carpi)                                                                                                                     
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
        #recording the simulation result
        benevis.write('\n')
        benevis.write(str(total_realized_demand)+'  '+str(total_American_value_d)+'  '+str(total_disposal_amount_D)+'  '+str(loss_ration_d)+'  '+str(obj_cor_d)+'  '+str(costamrican_c_d)+'  '+str(costdispose_c_d)+'  '+str(costsetupfirst_c_d)+'  '+str(costsetupd_c_d)+'  '+str(costholding_c_d)+'  '+str(costopeningwarehouse_c_d)+'  '+str(costbackorder_c_d)+'  '+str(costransport_c_d)+'  '+str(prod_det)+'  '
                      +str(newobj)+'  '+str(disposal_amount_rc)+'  '+str(costdispose_c_rc)+'  '+str(newbackordercost)+'  '+str(backorder_RC)+'   '+str(transport_RC)+'  '+str(prod_RC)+'  '
                      +str(obj_cor_Arc)+'   '+str(disp_cost_Arc)+'   '+str(total_disposal_value_Arc2)+'  '+str(costbackorder_c_Arc)+'  '+str(dist_samp_Arc)+'  '+str(prodcost_samp_Arc)+'  '
                      +str(obj_cor_RE)+'   '+str(disp_cost_RE)+'   '+str(total_exxtra_value_RE2)+'  '+str(costbackorder_c_RE)+'  '+str(dist_samp_RE)+'  '+str(prodcost_samp_RE))
    benevis.write('\n')
    benevis.write(str('**********END OF THIS SET***********'))
    benevis.write('\n')
    benevis.write(str('#############################################################################'))
    benevis.write('\n')
    benevis.write('\n')
    benevis.write('\n')
    benevis.close()     
    end= time.time()
    
    
    benevis=open('random data Sim.txt','a')
    benevis.write(str('***************'))
    benevis.write(str(jj)+'     '+str(kk)+'     '+str(tt)+'     '+str(rr)+'     '+str(ww)+'     '+str(zeta))
    benevis.write(str('***************'))
    dh=np.zeros((rr,jj,tt))
    realized_demand=np.zeros((rr,jj,tt))
    worsecase_demand=sum(sum(sum((D_W[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt))
    for iteration in range(iti):
        #calculate samplede demand for each iteration (realized_demand)
        #sample normaly    upper bound
        for r in range(rr):
            for j in range(jj):
                for t in range(tt): 
                    #dh[r][j][t]=dz[r][j][t]*(np.random.randint(-zeta*1000,zeta*1000)/1000) 
                    dh[r][j][t]=dz[r][j][t]*(np.random.randint(0,zeta*1000)/1000) 
                    #dh[r][j][t]=+dz[r][j][t]*0.2
                    realized_demand[r][j][t]=dh[r][j][t]+dz[r][j][t]
        #sample from a triangular 
    # =============================================================================
    #     randnums=(np.random.triangular(-zeta*1000,zeta*500,zeta*1000,rr*jj*tt)/1000)
    #     hi = plt.hist(randnums,density=True)
    #     plt.show()
    #     f=0
    #     for r in range(rr):
    #         for j in range(jj):
    #             for t in range(tt): 
    #                 dh[r][j][t]=dz[r][j][t]*randnums[f]
    #                 realized_demand[r][j][t]=dh[r][j][t]+dz[r][j][t]
    #                 f=f+1                
    # =============================================================================
        total_realized_demand=sum(sum(sum( realized_demand[r][j][t] for r in range (rr)) for j in range (jj)) for t in range (tt))         
        #deterministic feasiblity check
        from detsim import detsim
        American_value_d,disposal_amount_D,obj_cor_d,costamrican_c_d,costdispose_c_d,costbackorder_c_d,costransport_c_d,costsetupfirst_c_d,costsetupd_c_d,costholding_c_d,costopeningwarehouse_c_d=detsim(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,transperunit,CAP,realized_demand,o_w,q,sigma,bigM,C_O_W,Ydeter,prod_det,dis_per_unit,ini_inv, carpi)                                                                                                                     
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
        #recording the simulation result
        benevis.write('\n')
        benevis.write(str(total_realized_demand)+'  '+str(total_American_value_d)+'  '+str(total_disposal_amount_D)+'  '+str(loss_ration_d)+'  '+str(obj_cor_d)+'  '+str(costamrican_c_d)+'  '+str(costdispose_c_d)+'  '+str(costsetupfirst_c_d)+'  '+str(costsetupd_c_d)+'  '+str(costholding_c_d)+'  '+str(costopeningwarehouse_c_d)+'  '+str(costbackorder_c_d)+'  '+str(costransport_c_d)+'  '+str(prod_det)+'  '
                      +str(newobj)+'  '+str(disposal_amount_rc)+'  '+str(costdispose_c_rc)+'  '+str(newbackordercost)+'  '+str(backorder_RC)+'   '+str(transport_RC)+'  '+str(prod_RC)+'  '
                      +str(obj_cor_Arc)+'   '+str(disp_cost_Arc)+'   '+str(total_disposal_value_Arc2)+'  '+str(costbackorder_c_Arc)+'  '+str(dist_samp_Arc)+'  '+str(prodcost_samp_Arc)+'  '
                      +str(obj_cor_RE)+'   '+str(disp_cost_RE)+'   '+str(total_exxtra_value_RE2)+'  '+str(costbackorder_c_RE)+'  '+str(dist_samp_RE)+'  '+str(prodcost_samp_RE))
    benevis.write('\n')
    benevis.write(str('**********END OF THIS SET***********'))
    benevis.write('\n')
    benevis.write(str('#############################################################################'))
    benevis.write('\n')
    benevis.write('\n')
    benevis.write('\n')
    benevis.close()     
    end= time.time()
    
    
    

    zaman=open('random data BS.txt','a')
    zaman.write(str(Z_P_V_RC))
    zaman.write('\n')
    zaman.write(str('**********END OF THIS SET***********')+'       '+str('**********END OF THIS SET***********'))
    zaman.write('\n')
    zaman.write(str('###############################################################################################################'))
    zaman.write('\n')
    zaman.write('\n')
    zaman.write('\n')
    zaman.close()
    

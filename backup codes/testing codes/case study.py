#%%simulation integrated script 11.04.2020 backorder is optimized in simulation
#%% data is for case study 
# 1 may
import xlrd 
import numpy as np
import time
rnd = np.random
rnd.seed(1)
from docplex.mp.model import Model
# index and data
jj=5
kk=3
tt=7
rr=8
ww=2
zeta=0.1
#opening cost of warehouse
C_O_W=[0,180000,135000,100000]
q_p=[[120,200,200,200,200,100,0],[80,120,120,120,120,60,0],[65,110,110,110,110,55,0],[50,80,80,80,80,40,0],[55,90,90,90,90,90,45,0]]
q_w=np.zeros((ww,tt))
for w in range(ww):
    for t in range(tt): 
        q_w[w][t]=200
#cost of back rder holding and production
c_b=np.zeros((jj))
p=[5500,12000,8500,6000,7000]
c_b=[1900,2100,2000,2000,2000]
h=[0,3000,2500]
# Cots of transportation T in matlab
c=[0,2400,2000]
#Demand
dz=[[[15,15,15,15,15,15,15],[5,5,5,5,5,5,5],[7,7,7,7,7,7,7],[3,3,3,3,3,3,3],[5,5,5,5,5,5,5]],
    [[9,9,9,9,9,9,9],[8,8,8,8,8,8,8],[6,6,6,6,6,6,6],[8,8,8,8,8,8,8],[5,5,5,5,5,5,5]],
    [[15,15,15,15,15,15,15],[13,13,13,13,13,13,13],[10,10,10,10,10,10,10],[7,7,7,7,7,7,7],[5,5,5,5,5,5,5]],
    [[20,20,20,20,20,20,20],[10,10,10,10,10,10,10],[12,12,12,12,12,12,12],[6,6,6,6,6,6,6],[8,8,8,8,8,8,8]],
    [[15,15,15,15,15,15,15],[10,10,10,10,10,10,10],[10,10,10,10,10,10,10],[8,8,8,8,8,8,8],[12,12,12,12,12,12,12]],
    [[10,10,10,10,10,10,10],[5,5,5,5,5,5,5],[5,5,5,5,5,5,5],[3,3,3,3,3,3,3],[5,5,5,5,5,5,5]],
    [[8,8,8,8,8,8,8],[4,4,4,4,4,4,4],[6,6,6,6,6,6,6],[3,3,3,3,3,3,3],[5,5,5,5,5,5,5]],
    [[8,8,8,8,8,8,8],[3,3,3,3,3,3,3],[5,5,5,5,5,5,5],[4,4,4,4,4,4,4],[6,6,6,6,6,6,6]]]
D_W=np.zeros((rr,jj,tt))
dh=np.zeros((rr,jj,tt))
for r in range(rr):
    for j in range(jj):
         for t in range(tt): 
             D_W[r][j][t]=dz[r][j][t]*(1+zeta)                
             dh[r][j][t]=dz[r][j][t]*zeta
#possiblity of opening aplant or warehouse
o_p=[[0,1,1],[0,1,0],[0,1,0],[0,1,0],[0,1,1]]
o_w=[[0,1,1],[0,1,0],[0,1,0],[0,1,0],[0,1,1]]
CAP=500
#lead time
# =============================================================================
# q=np.zeros((rr))
# for r in range ((rr)):
#     q[r]=1
# =============================================================================
q=[0,1,3,3,1,2,2,2]
#setup
# =============================================================================
# sigma=np.zeros((kk,kk))
# for k in range ((kk)):
#     for kprime in  ((min(k+1,kk-1),kk-1)):
#         sigma[k][kprime]=np.random.randint(150,200)
# sigma[kk-1][kk-1]=0
# =============================================================================
sigma=[[0,80000,60000],[0,0,0],[0,20000,0]]
bigM=sum(sum(sum(D_W[r][j][t] for r in range (rr)) for j in range (jj)) for t in range (tt))*1.5
SIMstart = time.time()
#%%deterministic function  we use RC wit nominal data input
from NewRCfunction import RC
objdeter,Z_P_V_deter,Ydeter,Xdeter,Bdeter,setupdet,setupdeter,DETmodeling,DETtotal,deterbackorder,detprodcost,detholding,dettransport,detcostopenin=RC(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,c,CAP,dz,o_p,o_w,q,sigma,bigM,C_O_W)
zaman=open('timingcase.txt','a')
zaman.write(str(DETmodeling)+'     '+str(DETtotal)+'     '+str(objdeter)+'     '+str(setupdeter)+'     '+str(setupdet)+'     '+str(detprodcost)+'     '+str(detholding)+'     '+str(deterbackorder)+'     '+str(dettransport)+'     '+str(detcostopenin))
zaman.write('\n')
zaman.close()
farzad=5
#%%roboust counterpart with DW
from NewRCfunction import RC 
objRC,Z_P_V_RC,YRC,XRC,Bdeter,setupupN,setupfirst,RCmodeling,RCtotal,RCrbackorder,RCprodcost,RCholding,RCtransport,RCcostopen=RC(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,c,CAP,D_W,o_p,o_w,q,sigma,bigM,C_O_W) 
zaman=open('timingcase.txt','a')
zaman.write(str(RCmodeling)+'     '+str(RCtotal)+'     '+str(objRC)+'     '+str(setupfirst)+'     '+str(setupupN)+'     '+str(RCprodcost)+'     '+str(RCholding)+'     '+str(RCrbackorder)+'     '+str(RCtransport)+'     '+str(RCcostopen))
zaman.write('\n')
zaman.close()
#%% Adjustable robust counterpart
from NewARCfunction import ARC  
ARCobj,YARC,alphaARC,ARCprodcost,XARC,betaARC,ARCtranscost,BARC,gammaARC,ARCbackcost,ARCholdcost,ARCmodeling,ARCtotal=ARC(jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,c,CAP,D_W,dz,dh,o_p,o_w,q,sigma,bigM,setupupN,setupfirst,RCcostopen,Z_P_V_RC) 
zaman=open('timingcase.txt','a')
zaman.write(str(ARCmodeling)+'     '+str(ARCtotal)+'     '+str(ARCobj)+'     '+str(setupfirst)+'     '+str(setupupN)+'     '+str(ARCprodcost)+'     '+str(ARCholdcost)+'     '+str(ARCbackcost)+'     '+str(ARCtranscost)+'     '+str(RCcostopen))
zaman.write('\n')
zaman.close()
#%% reoptimization
from NewREARCfunction import REARC
(reoptobj,YB,alphaN,Reprodcost,XB,betaN,REtranscost,BB,gammaN,REbackcost,REholdcost,REARCmodeling,REARCtotal)=REARC (jj,kk,tt,rr,ww,q_p,q_w,c_b,h,p,c,CAP,D_W,dz,dh,o_p,o_w,q,sigma,bigM,setupupN,setupfirst,RCcostopen,Z_P_V_RC,ARCobj)                               
zaman=open('timingcase.txt','a')
zaman.write(str(REARCmodeling)+'     '+str(REARCtotal)+'     '+str(reoptobj)+'     '+str(setupfirst)+'     '+str(setupupN)+'     '+str(Reprodcost)+'     '+str(REholdcost)+'     '+str(REbackcost)+'     '+str(REtranscost)+'     '+str(RCcostopen))
zaman.write('\n')
zaman.close()
#%% simulation starts
dh=np.zeros((rr,jj,tt))
realized_demand=np.zeros((rr,jj,tt))
worsecase_demand=sum(sum(sum((D_W[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt))
for iteration in range(100):
    #calculate samplede demand for each iteration (realized_demand)
    for r in range(rr):
        for j in range(jj):
            for t in range(tt): 
                dh[r][j][t]=dz[r][j][t]*(np.random.randint(-zeta*1000,zeta*1000)/1000) 
                #dh[r][j][t]=dz[r][j][t]*0.2
                realized_demand[r][j][t]=dh[r][j][t]+dz[r][j][t]
    total_realized_demand=sum(sum(sum( realized_demand[r][j][t] for r in range (rr)) for j in range (jj)) for t in range (tt))     
    #deterministic feasiblity check
    def feasiblitycheck():
        mdl=Model('feasiblitycheck')
        Bsimdeter=mdl.continuous_var_dict([(r,j,t) for r in range(rr) for j in range(jj) for t in range(tt)], name='Bsimdeter')
        #Inventory
        mdl.minimize(mdl.sum(mdl.sum(mdl.sum(c_b[j]*Bsimdeter[(r,j,t)]for r in range(rr))for j in range(jj))for t in range(tt)))        
        mdl.add_constraints(mdl.sum(mdl.sum(Xdeter[(w,r,j,k,t-q[r])] for w in range(ww))for k in range(kk))>= realized_demand[r][j][t]-Bsimdeter[(r,j,t)]+Bsimdeter[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t!=0))
        mdl.add_constraints(mdl.sum(mdl.sum(Xdeter[(w,r,j,k,t-q[r])] for w in range(ww))for k in range(kk))>= realized_demand[r][j][t]-Bsimdeter[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t==0))
        mdl.add_constraints(0>= realized_demand[r][j][t]-Bsimdeter[(r,j,t)]+Bsimdeter[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t!=0))
        mdl.add_constraints(0>= realized_demand[r][j][t]-Bsimdeter[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t==0))
        
    #4 10latex     
        mdl.add_constraint(mdl.sum(mdl.sum(Bsimdeter[(r,j,tt-1)] for r in range(rr))for j in range(jj))== 0 )
        solution = mdl.solve(log_output= True)     
        unsatisfieddemand=np.zeros((rr,jj,tt))
        Counter=np.zeros((rr,jj,tt))
        try:
          feasible = 1
          simdetbackorder = solution.get_objective_value()
        except:
          feasible = 0
          obj=-125          
        if feasible == 1:            
            deterransported_amount=sum(sum(sum(sum(sum(Xdeter[w][r][j][k][t]for r in range (rr))for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))
            deter_penalty=deterransported_amount-total_realized_demand
            deter_objective=objdeter-deterbackorder+simdetbackorder
            total_un_realized_demand=0
            total_un_realized_demand_ratio=0
        else:
                total_un_realized_demand=total_realized_demand-deterransported_amount
                total_un_realized_demand_ratio=(total_realized_demand-deterransported_amount)/total_realized_demand
                deter_penalty=0
                deter_objective=0
                for r in range(rr):
                    for j in range(jj):
                        for t in range (tt):
                            if (t-q[r]>=0 and t!=0):
                                unsatisfieddemand[r][j][t]=sum(sum(Xdeter[w][r][j][k][t-q[r]] for w in range (ww))for k in range(kk))- realized_demand[r][j][t]+Bdeter[r][j][t]-Bdeter[r][j][t-1]
                            elif   (t-q[r]>=0 and t==0):
                                unsatisfieddemand[r][j][t]=sum(sum(Xdeter[w][r][j][k][t-q[r]] for w in range (ww))for k in range(kk))- realized_demand[r][j][t]+Bdeter[r][j][t]
                            elif  (t-q[r]<0 and t!=0):
                                unsatisfieddemand[r][j][t]= realized_demand[r][j][t]-Bdeter[r][j][t]+Bdeter[r][j][t-1] 
                            else:
                                unsatisfieddemand[r][j][t]= realized_demand[r][j][t]-Bdeter[r][j][t]                            
                            if unsatisfieddemand[r][j][t]<0:
                                Counter[r][j][t]=1
                            else:
                                Counter[r][j][t]=0
        return(deter_objective,deter_penalty,Counter,total_un_realized_demand,total_un_realized_demand_ratio)
    deter_objective,deter_penalty,Counter,total_un_realized_demand,total_un_realized_demand_ratio=feasiblitycheck()       
    #ARC
    def SimRC():
        mdl=Model('SimRC')
        B_SimRC=mdl.continuous_var_dict([(r,j,t) for r in range(rr) for j in range(jj) for t in range(tt)], name='B_SimRC')                                   
        #backorder    
        mdl.minimize(mdl.sum(mdl.sum(mdl.sum(c_b[j]*B_SimRC[(r,j,t)]for r in range(rr))for j in range(jj))for t in range(tt)))    
        # Constraints
        mdl.add_constraints(mdl.sum(mdl.sum(XRC[(w,r,j,k,t-q[r])] for k in range(kk))for w in range(ww))>= realized_demand[r][j][t]-B_SimRC[(r,j,t)]+B_SimRC[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t!=0))
        mdl.add_constraints(mdl.sum(mdl.sum(XRC[(w,r,j,k,t-q[r])] for k in range(kk))for w in range(ww))>= realized_demand[r][j][t]-B_SimRC[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t==0))
        mdl.add_constraints(0>= realized_demand[r][j][t]-B_SimRC[(r,j,t)]+B_SimRC[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t!=0))
        mdl.add_constraints(0>= realized_demand[r][j][t]-B_SimRC[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t==0))
        #4 10latex     
        mdl.add_constraint(mdl.sum(mdl.sum(B_SimRC[(r,j,tt-1)] for r in range(rr))for j in range(jj))== 0 )
        solution = mdl.solve(log_output= True)
        B_SimRCvalue=np.zeros((rr,jj,tt))
        simRCbackordercost = solution.get_objective_value()
        for r in range (rr):
            for j in range (jj):
                for t in range (tt):
                    B_SimRCvalue[r][j][t]=float(B_SimRC[(r,j,t)].solution_value)
        simRCbackorderamount=sum(sum(sum( B_SimRCvalue[r][j][t]for r in range (rr))for j in range (jj))for t in range (tt))    
        return(simRCbackordercost,simRCbackorderamount,B_SimRCvalue) 
    simRCbackordercost,simRCbackorderamount,B_SimRCvalue=SimRC()
    worsecase_penaly=worsecase_demand-total_realized_demand
    worsecase_objective=objRC+simRCbackordercost-RCrbackorder
    RCOBJIMPSIM=RCrbackorder-simRCbackordercost
    def SimARC():
        mdl=Model('SimARC')
        B_SimARC=mdl.continuous_var_dict([(r,j,t) for r in range(rr) for j in range(jj) for t in range(tt)], name='B_SimARC')                                     
        #backorder    
        mdl.minimize(mdl.sum(mdl.sum(mdl.sum(c_b[j]*B_SimARC[(r,j,t)]for r in range(rr))for j in range(jj))for t in range(tt)))    
        # Constraints
        #3 9latex     
        mdl.add_constraints(mdl.sum(mdl.sum(XARC[(w,r,j,k,t-q[r])]+mdl.sum(mdl.sum(mdl.sum((betaARC[(w,r,j,k,t-q[r],rp,jp,tp)]*realized_demand[rp][jp][tp])for rp in range (rr))for jp in range (jj))for tp in range (tt)) for w in range (ww))for k in range(kk))>= realized_demand[r][j][t]-B_SimARC[(r,j,t)]+B_SimARC[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t!=0))
        mdl.add_constraints(mdl.sum(mdl.sum(XARC[(w,r,j,k,t-q[r])]+mdl.sum(mdl.sum(mdl.sum((betaARC[(w,r,j,k,t-q[r],rp,jp,tp)]*realized_demand[rp][jp][tp])for rp in range (rr))for jp in range (jj))for tp in range (tt)) for w in range (ww))for k in range(kk))>= realized_demand[r][j][t]-B_SimARC[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t==0))
        mdl.add_constraints(0>= realized_demand[r][j][t]-B_SimARC[(r,j,t)]+B_SimARC[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t!=0))
        mdl.add_constraints(0>= realized_demand[r][j][t]-B_SimARC[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t==0))    
       #4 10latex   
        mdl.add_constraint(mdl.sum(mdl.sum(B_SimARC[(r,j,tt-1)] for r in range(rr))for j in range(jj))== 0 )    
        mdl.export('SimARC4')
        try:
          feasible = 1
          solution = mdl.solve(log_output= True)
          B_SimARC_value=np.zeros((rr,jj,tt))
          simARCbackordercost = solution.get_objective_value()
          for r in range (rr):
              for j in range (jj):
                  for t in range (tt):
                      B_SimARC_value[r][j][t]=float(B_SimARC[(r,j,t)].solution_value)
        except:
          feasible = 0
          ARCobj=-125
        simARCbackorderamount=sum(sum(sum( B_SimARC_value[r][j][t]for r in range (rr))for j in range (jj))for t in range (tt))
        return(simARCbackordercost) 
    simARCbackordercost=SimARC()   
    
    ARCproduction_cost=sum(sum(sum(sum(p[j]*(YARC[w][j][k][t]+sum(sum(sum(alphaARC[w][j][k][t][rp][jp][tp]*(dz[rp][jp][tp]+dh[rp][jp][tp]) for rp in range (rr)) for jp in range (jj)) for tp in range (tt)))for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))
    ARCproduction_amount=sum(sum(sum(sum((YARC[w][j][k][t]+sum(sum(sum(alphaARC[w][j][k][t][rp][jp][tp]*(dz[rp][jp][tp]+dh[rp][jp][tp]) for rp in range (rr)) for jp in range (jj)) for tp in range (tt)))for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))
    ARCinventory_cost=sum(sum(sum(h[k]*sum(sum((YARC[w][j][k][tp]+sum(sum(sum(alphaARC[w][j][k][tp][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt)))for tp in range (t+1))for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))-(sum(sum(sum(sum(h[k]*sum(sum((XARC[w][r][j][k][tp]+sum(sum(sum(betaARC[w][r][j][k][tp][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt)))for tp in range (t+1))for r in range (rr))for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt)))
    ARCinventory_amount=sum(sum(sum(sum(sum((YARC[w][j][k][tp]+sum(sum(sum(alphaARC[w][j][k][tp][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt)))for tp in range (t+1))for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))-sum((sum(sum(sum(sum(sum((XARC[w][r][j][k][tp]+sum(sum(sum(betaARC[w][r][j][k][tp][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt)))for tp in range (t+1))for w in range (ww))for r in range (rr))for j in range (jj))for k in range (kk))for t in range (tt)))
    ARCtransportcost=sum(sum(sum(sum(sum(c[k]*(XARC[w][r][j][k][t]+sum(sum(sum(betaARC[w][r][j][k][t][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt)))for w in range (ww))for r in range (rr))for j in range (jj))for k in range (kk))for t in range (tt))    
    ARCtransported_amount=sum(sum(sum(sum(sum(XARC[w][r][j][k][t]+sum(sum(sum(betaARC[w][r][j][k][t][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt))for w in range (ww))for r in range (rr))for j in range (jj))for k in range (kk))for t in range (tt))                
    ARCbackordercost=sum(sum(sum(c_b[j]*(BARC[r][j][t]+sum(sum(sum(gammaARC[r][j][t][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt)))for r in range (rr))for j in range (jj))for t in range (tt))  
    ARCbackorder_amount=sum(sum(sum((BARC[r][j][t]+sum(sum(sum(gammaARC[r][j][t][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt)))for r in range (rr))for j in range (jj))for t in range (tt))  
    ARC_penalty=ARCtransported_amount- total_realized_demand
    simARC_obj=ARCproduction_cost+setupupN+ARCinventory_cost+setupfirst+ARCtransportcost+simARCbackordercost
    ARCOBJIMPSIM=ARCbackordercost-simARCbackordercost
    def SimREARC():
        mdl=Model('SimREARC')
        B_SimREARC=mdl.continuous_var_dict([(r,j,t) for r in range(rr) for j in range(jj) for t in range(tt)], name='B_SimREARC')                                        
        #backorder    
        mdl.minimize(mdl.sum(mdl.sum(mdl.sum(c_b[j]*B_SimREARC[(r,j,t)]for r in range(rr))for j in range(jj))for t in range(tt)))    
        # Constraints
        #3 9latex     
        mdl.add_constraints(mdl.sum(mdl.sum(XB[(w,r,j,k,t-q[r])]+mdl.sum(mdl.sum(mdl.sum((betaN[(w,r,j,k,t-q[r],rp,jp,tp)]*realized_demand[rp][jp][tp])for rp in range (rr))for jp in range (jj))for tp in range (tt)) for k in range(kk))for w in range(ww))>= realized_demand[r][j][t]-B_SimREARC[(r,j,t)]+B_SimREARC[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t!=0))
        mdl.add_constraints(mdl.sum(mdl.sum(XB[(w,r,j,k,t-q[r])]+mdl.sum(mdl.sum(mdl.sum((betaN[(w,r,j,k,t-q[r],rp,jp,tp)]*realized_demand[rp][jp][tp])for rp in range (rr))for jp in range (jj))for tp in range (tt)) for k in range(kk))for w in range(ww))>= realized_demand[r][j][t]-B_SimREARC[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t==0))
        mdl.add_constraints(0>= realized_demand[r][j][t]-B_SimREARC[(r,j,t)]+B_SimREARC[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t!=0))
        mdl.add_constraints(0>= realized_demand[r][j][t]-B_SimREARC[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t==0))    
       #4 10latex   
        mdl.add_constraint(mdl.sum(mdl.sum(B_SimREARC[(r,j,tt-1)] for r in range(rr))for j in range(jj))== 0 )    
        try:
          feasible = 1
          solution = mdl.solve(log_output= True)
          simREARCbackordercost = solution.get_objective_value()
        except:
          feasible = 0
          ARCobj=-125
        return(simREARCbackordercost) 
    simREARCbackordercost=SimREARC()
    #REARC
    REproduction_cost=sum(sum(sum(sum(p[j]*(YB[w][j][k][t]+sum(sum(sum(alphaN[w][j][k][t][rp][jp][tp]*(dz[rp][jp][tp]+dh[rp][jp][tp]) for rp in range (rr)) for jp in range (jj)) for tp in range (tt)))for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))
    REinventory_cost=sum(sum(sum(h[k]*sum(sum((YB[w][j][k][tp]+sum(sum(sum(alphaN[w][j][k][tp][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt)))for tp in range (t+1))for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))-(sum(sum(sum(h[k]*sum(sum(sum((XB[w][r][j][k][tp]+sum(sum(sum(betaN[w][r][j][k][tp][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt)))for tp in range (t+1))for r in range (rr))for j in range (jj))for k in range (kk))for w in range (ww))for t in range (tt)))
    REtransport=sum(sum(sum(sum(sum(c[k]*(XB[w][r][j][k][t]+sum(sum(sum(betaN[w][r][j][k][t][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt)))for w in range (ww))for r in range (rr))for j in range (jj))for k in range (kk))for t in range (tt))    
    REtransported_amount=sum(sum(sum(sum(sum(XB[w][r][j][k][t]+sum(sum(sum(betaN[w][r][j][k][t][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt))for w in range (ww))for r in range (rr))for j in range (jj))for k in range (kk))for t in range (tt))
    REbackordercost=sum(sum(sum(c_b[j]*(BB[r][j][t]+sum(sum(sum(gammaN[r][j][t][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt)))for r in range (rr))for j in range (jj))for t in range (tt))  
    REbackorder_amount=sum(sum(sum((BB[r][j][t]+sum(sum(sum(gammaN[r][j][t][rp][jp][tz]*(dz[rp][jp][tz]+dh[rp][jp][tz]) for rp in range (rr)) for jp in range (jj)) for tz in range (tt)))for r in range (rr))for j in range (jj))for t in range (tt))  
    RE_penalty=REtransported_amount- total_realized_demand
    REOP_objective=setupupN+setupfirst+REproduction_cost+REinventory_cost+simREARCbackordercost+REtransport
    REOBJIMPSIM=REbackordercost-simREARCbackordercost
    benevis=open('Resulcase.txt','a')
    unsatisfieddemandcounter=sum(sum(sum(Counter[r][j][t]for r in range (rr))for j in range (jj))for t in range (tt))
    benevis.write(str(deter_objective)+'  '+str(deter_penalty)+'  '+str(unsatisfieddemandcounter)+'  '+str(total_un_realized_demand)+'  '+str(total_un_realized_demand_ratio)+'  '
                  +str(worsecase_objective)+'  '+str(worsecase_penaly)+'  '+str(RCOBJIMPSIM)+'  '+str(RCprodcost)+'  '+str(RCholding)+'  '+str(simRCbackordercost)+'  '+str(RCtransport)+'  '
                  +str(simARC_obj)+'  '+str(ARC_penalty)+'  '+str(ARCOBJIMPSIM)+'  '+str(ARCproduction_cost)+'  '+str(ARCinventory_cost)+'  '+str(simARCbackordercost)+'  '+str(ARCtransportcost)+'  '
                  +str(REOP_objective)+'  '+str(RE_penalty)+'  '+str(REOBJIMPSIM)+'  '+str(REproduction_cost)+'  '+str(REinventory_cost)+'  '+str(simREARCbackordercost)+'  '+str(REtransport))
    benevis.write('\n')
    benevis.close() 
    
end= time.time()
SIMtime=end-SIMstart
zaman=open('timingcase.txt','a')
zaman.write(str(SIMtime)+'     '+str(SIMtime)+'     '+str(rr)+'     '+str(jj)+'     '+str(kk)+'     '+str(ww)+'     '+str(tt))
zaman.write('\n')
zaman.close()
print(setupfirst,setupupN)

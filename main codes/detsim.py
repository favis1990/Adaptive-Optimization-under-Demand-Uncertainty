

import time
import numpy as np
#import matplotlib.pyplot as plt
from docplex.mp.model import Model
def detsim(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, D_W, o_w, q, sigma, bigM, C_O_W, Y_det, podcost_det, dis_per_unit, ini_inv, carpi):
    RCstart = time.time()
    ''' Model '''
    mdl=Model('DeterministicMilk')
    mdl.parameters.mip.strategy.file=3
    # Decision Variables
    #material flow
    X=mdl.continuous_var_dict([(r,j,t) for r in range(rr) for j in range(jj) for t in range(tt)], name='x')
    #exxtra baraye emhaye ezafe ha disposal
    exxtra=mdl.continuous_var_dict([(r,j,t)  for r in range(rr) for j in range(jj)  for t in range(tt)], name='extra')    
    #american komak be production
    american=mdl.continuous_var_dict([(j,t)  for j in range(jj)  for t in range(tt)], name='american')    
    #Inventory
    I=mdl.continuous_var_dict([(w,j,k,t)  for w in range(ww) for j in range(jj) for k in range(kk) for t in range(tt)], name='I')
    #backorder
    B=mdl.continuous_var_dict([(r,j,t) for r in range(rr) for j in range(jj) for t in range(tt)], name='B')
    # Auxiliary variable for linearization os setup cost
    AX=mdl.continuous_var_dict([(w,k,kprime,t) for w in range(ww) for k in range(kk) for kprime in range(kk) for t in range(tt)], name='AX')
    #openning Binary variable
    Z_P=mdl.binary_var_dict([(w,k,t) for w in range(ww) for k in range(kk) for t in range(tt)], name='Z_P')
    #cost variables is used to seprate objective function for the sake of simplicity
    americancost = mdl.continuous_var(name='americancost')
    setupfirst = mdl.continuous_var(name='setupfirst')
    holding = mdl.continuous_var(name='holding')
    setup = mdl.continuous_var(name='setup')
    openingw = mdl.continuous_var(name='openingw')
    transport = mdl.continuous_var(name='transport')
    backorder = mdl.continuous_var(name='backorder')
    disposecost = mdl.continuous_var(name='disposecost')
    # Objective Function
    mdl.minimize(podcost_det+ setupfirst+ holding+ setup+ transport+ backorder+ openingw+ americancost+ disposecost)
    mdl.add_constraint(americancost == mdl.sum(mdl.sum(p[j]*carpi*american[(j,t)] for j in range(jj))for t in range(tt)))
    mdl.add_constraint(setupfirst == mdl.sum(mdl.sum(sigma[0][kprime]*Z_P[(w,kprime,0)]for w in range(ww))for kprime in range(kk)))
    mdl.add_constraint(holding == mdl.sum(mdl.sum(mdl.sum(mdl.sum(h[k][j]*I[(w,j,k,t)]for w in range(ww))for j in range(jj))for k in range(kk))for t in range(tt)))
    mdl.add_constraint(setup == mdl.sum(mdl.sum(mdl.sum(mdl.sum(sigma[k][kprime]*AX[(w,k,kprime,t)]for w in range(ww))for k in range(kk))for kprime in range(kk))for t in range(tt-1)))
    mdl.add_constraint(transport == mdl.sum(mdl.sum(mdl.sum(transperunit[j][r]*X[(r,j,t)]for r in range(rr))for j in range(jj))for t in range(tt)))
    mdl.add_constraint(backorder == mdl.sum(mdl.sum(mdl.sum(c_b[j]*B[(r,j,t)]for r in range(rr))for j in range(jj))for t in range(tt)))
    mdl.add_constraint(openingw == mdl.sum(mdl.sum(mdl.sum(Z_P[(w,k,t)]*C_O_W[(k)]for w in range(ww))for k in range(kk))for t in range(tt)))
    mdl.add_constraint(disposecost == mdl.sum(mdl.sum(mdl.sum(dis_per_unit*exxtra[(r,j,t)]for r in range(rr))for j in range(jj))for t in range(tt)))
    #2 8latex
    mdl.add_constraints(mdl.sum(mdl.sum(I[(w,j,k,t)]for j in range(jj))for k in range(kk))<=q_w[w][t]  for w in range(ww) for t in range(tt))
    #3
    mdl.add_constraints(X[(r,j,t-q[r])] == D_W[r][j][t]-B[(r,j,t)]+B[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t!=0))
    mdl.add_constraints(X[(r,j,t-q[r])] == D_W[r][j][t]-B[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t==0))
    mdl.add_constraints(0== D_W[r][j][t]-B[(r,j,t)]+B[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t!=0))
    mdl.add_constraints(0== D_W[r][j][t]-B[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t==0))
    #lead time cosntraint
    mdl.add_constraints(mdl.sum(X[(r,j,tp)] for tp in range(t+1))>= mdl.sum(D_W[r][j][tp] for tp in range(t+1))for r in range (rr) for j in range(jj) for t in range(tt))
    #4     
    mdl.add_constraint(mdl.sum(mdl.sum(B[(r,j,tt-1)] for r in range(rr))for j in range(jj))== 0 )
    #5 
    mdl.add_constraints(Y_det[j][0]+american[(j,0)]== mdl.sum(mdl.sum(I[(w,j,k,0)] -ini_inv[(w,j,k)] for k in range (kk))for w in range (ww))+mdl.sum(X[(r,j,0)]+exxtra[(r,j,0)] for r in range (rr)) for j in range (jj))
    mdl.add_constraints(Y_det[j][t]+american[(j,t)]==  mdl.sum(mdl.sum(I[(w,j,k,t)] -I[(w,j,k,t-1)]for k in range (kk))for w in range (ww)) +mdl.sum(X[(r,j,t)]+exxtra[(r,j,t)]for r in range (rr)) for j in range (jj) for t in range (1,tt))
    #6 
    mdl.add_constraints(mdl.sum(mdl.sum (I[(w,j,k,t)] for k in range (kk)) for w in range (ww)) <=CAP[j] for j in range (jj)for t in range(tt))
    #7
    mdl.add_constraints(mdl.sum(I[(w,j,k,t)] for j in range(jj))<=bigM *Z_P[(w,k,t)] for w in range(ww) for k in range(kk) for t in range(tt))
    #8 
    mdl.add_constraints(mdl.sum(Z_P[(w,k,t)] for k in range(kk))==1 for w in range(ww) for t in range(tt))
    #11 15
    mdl.add_constraints(mdl.sum(mdl.sum(I[(w,j,k,t)]for w in range(ww))for t in range(tt))<= bigM *o_w[j][k] for j in range(jj)for k in range(kk))
        # to finish with some nventory 
    mdl.add_constraints(I[(w,j,k,tt-1)]==ini_inv[w][j][k] for w in range (ww) for j in range (jj) for k in range (kk))
    #linearization of setup cost
    mdl.add_constraints(Z_P[(w,k,t)]+Z_P[(w,kprime,t+1)]-1<=AX[(w,k,kprime,t)]  for w in range(ww) for k in range(kk)for kprime in range(kk) for t in range(tt-1))
    mdl.add_constraints(AX[(w,k,kprime,t)]<=Z_P[(w,kprime,t+1)]  for w in range(ww) for k in range(kk)for kprime in range(kk) for t in range(tt-1))
    mdl.add_constraints(AX[(w,k,kprime,t)]<=Z_P[(w,k,t)]  for w in range(ww) for k in range(kk)for kprime in range(kk) for t in range(tt-1))
    
    mdl.export('DeterministicMilk')
    modelingRC = time.time()
    Z_P_V=np.zeros((ww,kk,tt))
    X_value=np.zeros((rr,jj,tt))
    disposal_valu=np.zeros((rr,jj,tt))
    I_value=np.zeros((ww,jj,kk,tt))
    American_value=np.zeros((jj,tt))
    B_value=np.zeros((rr,jj,tt))
    chooseRC = time.time()
    DETmodeling=modelingRC-RCstart
    DETtotal=chooseRC-modelingRC
       
    try:
      solution = mdl.solve(log_output= True)
      feasible = 1
      objdeter = solution.get_objective_value()

                  
    except:
      feasible = 0
      obj=-125

    for j in range (jj):
        for t in range (tt):
            American_value[j][t]=float(american[(j,t)].solution_value)
                       
                      
    for w in range(ww):
        for k in range (kk):
            for t in range (tt):
               Z_P_V[w][k][t]=float(Z_P[(w,k,t)].solution_value)
               for j in range(jj):
                   I_value[w][j][k][t]=float(I[(w,j,k,t)].solution_value)
    for r in range (rr) :
        for j in range (jj):
            for t in range (tt):
                B_value[r][j][t]=float(B[(r,j,t)].solution_value)
                X_value[r][j][t]=float(X[(r,j,t)].solution_value)
                disposal_valu[r][j][t]=float(exxtra[(r,j,t)].solution_value)
    costsetupfirst=float(setupfirst)
    costsetupd=float(setup)
    costholding=float(holding)
    costopeningwarehouse=float(openingw)
    costbackorder=float(backorder)
    costransport=float(transport)
    costdispose=float(disposecost)
    costamrican=float(americancost)
    total_Y = mdl.sum(mdl.sum(Y_det[j][t]for j in range (jj))for t in range (tt))
    total_X = mdl.sum(mdl.sum(mdl.sum(X_value[r][j][t]for j in range (jj))for t in range (tt))for r in range (rr))
    total_I = mdl.sum(mdl.sum(mdl.sum(mdl.sum(I_value[w][j][k][t]for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))
    final_b = mdl.sum(mdl.sum(B_value[r][j][tt-1]for r in range (rr))for j in range (jj))
    real_demand = mdl.sum(mdl.sum(mdl.sum(D_W[r][j][t]for r in range (rr))for j in range (jj)) for t in range (tt))
    print(solution)
    return(American_value,disposal_valu,objdeter,costamrican,costdispose,costbackorder,costransport,costsetupfirst,costsetupd,costholding,costopeningwarehouse)



#sol=detRCsim()
# Deterministic (Robust) function 
#if in create data file ou set zeta = 0 it is detereministic

import time
import numpy as np
from docplex.mp.model import Model
def RC(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, D_W, o_w, q, sigma, bigM, C_O_W, ini_inv):
    RCstart = time.time()
    ''' Model '''
    mdl=Model('DeterministicMilk')
    mdl.parameters.mip.strategy.file=3
    # Decision Variables
    #material flow
    X=mdl.continuous_var_dict([(r,j,t) for r in range(rr) for j in range(jj) for t in range(tt)], name='x')
    #production amount
    Y=mdl.continuous_var_dict([(j,t) for j in range(jj) for t in range(tt)], name='Y')
    #Inventory
    I=mdl.continuous_var_dict([(w,j,k,t)  for w in range(ww) for j in range(jj) for k in range(kk) for t in range(tt)], name='I')
    #backorder
    B=mdl.continuous_var_dict([(r,j,t) for r in range(rr) for j in range(jj) for t in range(tt)], name='B')
    # Auxiliary variable for linearization os setup cost
    AX=mdl.continuous_var_dict([(w,k,kprime,t) for w in range(ww) for k in range(kk) for kprime in range(kk) for t in range(tt)], name='AX')
    #openning Binary variable
    Z_P=mdl.binary_var_dict([(w,k,t) for w in range(ww) for k in range(kk) for t in range(tt)], name='Z_P')
    #cost variables is used to seprate objective function for the sake of simplicity
    prodcost=mdl.continuous_var(name='prodcost')
    setupfirst=mdl.continuous_var(name='setupfirst')
    holding=mdl.continuous_var(name='holding')
    setup=mdl.continuous_var(name='setup')
    openingw=mdl.continuous_var(name='openingw')
    transport=mdl.continuous_var(name='transport')
    backorder=mdl.continuous_var(name='backorder')
    # Objective Function
    mdl.minimize(prodcost+setupfirst+holding+setup+transport+backorder+openingw)
    mdl.add_constraint(prodcost == mdl.sum(mdl.sum( p[j]*Y[(j,t)]  for j in range(jj))for t in range(tt)))
    mdl.add_constraint(setupfirst == mdl.sum(mdl.sum(sigma[0][kprime] * Z_P[(w,kprime,0)]for w in range(ww))for kprime in range(kk)))
    mdl.add_constraint(holding == mdl.sum(mdl.sum(mdl.sum(mdl.sum(h[k][j]*I[(w,j,k,t)]for w in range(ww))for j in range(jj))for k in range(kk))for t in range(tt)))
    mdl.add_constraint(setup == mdl.sum(mdl.sum(mdl.sum(mdl.sum(sigma[k][kprime]*AX[(w,k,kprime,t)]for w in range(ww))for k in range(kk))for kprime in range(kk))for t in range(tt-1)))
    mdl.add_constraint(transport == mdl.sum(mdl.sum(mdl.sum(transperunit[j][r]*X[(r,j,t)] for j in range(jj)) for r in range(rr)) for t in range(tt)))
    mdl.add_constraint(backorder==mdl.sum(mdl.sum(mdl.sum(c_b[j]*B[(r,j,t)]for r in range(rr))for j in range(jj))for t in range(tt)))
    mdl.add_constraint(openingw==mdl.sum(mdl.sum(mdl.sum(Z_P[(w,k,t)]*C_O_W[(k)]for w in range(ww))for k in range(kk))for t in range(tt)))
    # Constraints
    #warehouse capacity constraint 
    mdl.add_constraints(mdl.sum(mdl.sum(I[(w,j,k,t)]for j in range(jj))for k in range(kk))<=q_w[w][t]  for w in range(ww) for t in range(tt))
    # # #production capacity constraint
    mdl.add_constraints(Y[(j,t)] <= q_p[j][t]  for j in range(jj) for t in range(tt))
    #inventory balance constraint eliminated in ARC 
    mdl.add_constraints(Y[(j,0)] == mdl.sum(mdl.sum(I[(w,j,k,0)] - ini_inv[(w,j,k)] for k in  range (kk)) for w in range (ww)) + mdl.sum(X[(r,j,0)] for r in range (rr)) for j in range (jj))
    mdl.add_constraints(Y[(j,t)] == mdl.sum(mdl.sum(I[(w,j,k,t)] - I[(w,j,k,t-1)]  for k in  range (kk)) for w in range (ww)) + mdl.sum(X[(r,j,t)] for r in range (rr)) for j in range (jj)  for t in range (1,tt))
    #inventory balance constraint2 backorder 
    mdl.add_constraints( X[(r,j,t-q[r])] >= D_W[r][j][t]-B[(r,j,t)]+B[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t!=0))
    mdl.add_constraints( X[(r,j,t-q[r])] >= D_W[r][j][t]-B[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t==0))
    mdl.add_constraints(0>= D_W[r][j][t]-B[(r,j,t)]+B[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t!=0))
    mdl.add_constraints(0>= D_W[r][j][t]-B[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t==0))
    #last period backorder   
    mdl.add_constraint(mdl.sum(mdl.sum(B[(r,j,tt-1)] for r in range(rr))for j in range(jj))== 0 )
    #lead time cosntraint
    mdl.add_constraints(mdl.sum(X[(r,j,tp)] for tp in range(t+1))>= mdl.sum(D_W[r][j][tp] for tp in range(t+1))for r in range (rr) for j in range(jj) for t in range(tt))
    # inventory capacity 2
    mdl.add_constraints(mdl.sum(mdl.sum (I[(w,j,k,t)] for k in range (kk)) for w in range (ww)) <=CAP[j] for t in range(tt)for j in range(jj))
    # #Inventory openning variable
    mdl.add_constraints( mdl.sum(I[(w,j,k,t)] for j in range(jj))<=bigM *Z_P[(w,k,t)] for w in range(ww) for k in range(kk) for t in range(tt))
    #Inventory openning variable2
    mdl.add_constraints(mdl.sum(Z_P[(w,k,t)] for k in range(kk))==1 for w in range(ww) for t in range(tt))
    # cosntraint for k_1
    mdl.add_constraints(mdl.sum(mdl.sum(I[(w,j,k,t)]for w in range(ww))for t in range(tt))<= bigM *o_w[j][k] for j in range(jj)for k in range(kk))
    # to finish with same nventory 
    mdl.add_constraints(I[(w,j,k,tt-1)]==ini_inv[(w,j,k)] for w in range (ww)for k in range (kk) for j in range (jj))
    #linearization of setup cost
    mdl.add_constraints(Z_P[(w,k,t)]+Z_P[(w,kprime,t+1)]-1<=AX[(w,k,kprime,t)]  for w in range(ww) for k in range(kk)for kprime in range(kk) for t in range(tt-1))
    mdl.add_constraints(AX[(w,k,kprime,t)]<=Z_P[(w,kprime,t+1)]  for w in range(ww) for k in range(kk)for kprime in range(kk) for t in range(tt-1))
    mdl.add_constraints(AX[(w,k,kprime,t)]<=Z_P[(w,k,t)]  for w in range(ww) for k in range(kk)for kprime in range(kk) for t in range(tt-1))
    #mdl.export('DeterministicMilk1')
    modelingRC = time.time()
    Z_P_V= np.zeros((ww,kk,tt))
    X_value= np.zeros((rr,jj,tt))
    I_value= np.zeros((ww,jj,kk,tt))
    y_value= np.zeros((jj,tt))
    B_value= np.zeros((rr,jj,tt))
    chooseRC = time.time()
    DETmodeling= modelingRC-RCstart
    DETtotal= chooseRC-modelingRC
    with open('model898.txt', 'w') as file:
        mdl.export_as_lp(file.name)

       
    try:
      solution = mdl.solve(log_output= True)
      feasible = 1
      objdeter = solution.get_objective_value()

                  
    except:
      feasible = 0
      obj=-125

    for j in range (jj):
        for t in range (tt):
            y_value[j][t]=float(Y[(j,t)].solution_value)
            for w in range (ww):
               for k in range (kk):
                   I_value[w][j][k][t]=(float(I[(w,j,k,t)].solution_value))
    for w in range(ww):
        for k in range (kk):
            for t in range (tt):
               if float(Z_P[(w,k,t)].solution_value)>= 0.01:
                   Z_P_V[w][k][t]=1
    for r in range (rr) :
        for j in range (jj):
            for t in range (tt):
                B_value[r][j][t]=float(B[(r,j,t)].solution_value)
                X_value[r][j][t]=float(X[(r,j,t)].solution_value)

                
    total_Y = sum(sum( y_value[j][t] for j in range (jj))for t in range (tt))
    total_X = sum(sum(sum( X_value[r][j][t] for j in range (jj)) for r in range (rr)) for t in range (tt))
    total_B = sum(sum(sum( B_value[r][j][t] for j in range (jj)) for r in range (rr)) for t in range (tt))
    total_I = sum(sum(sum(sum( I_value[w][j][k][t] for w in range (ww)) for j in range (jj)) for k in range (kk)) for t in range (tt))
    totsl_initi= sum(sum(sum (ini_inv[w][j][k] for w in range (ww))for j in range (jj))for k in range (kk))
    costsetupfirst = float(setupfirst)
    costsetupd = float(setup)
    costbackorder = float(backorder)
    costprodcost = float(prodcost)
    costholding = float(holding)
    costransport = float(transport)
    costopeningwarehouse = float(openingw)
    number_of_DV_RC_model=(rr*jj*tt)+(jj*tt)+(ww*jj*kk*tt)+(rr*jj*tt)+(ww*kk*kk*tt)+(ww*kk*tt)


    
    return(objdeter, Z_P_V, y_value, X_value, B_value, I_value, costsetupfirst, costsetupd, DETmodeling, DETtotal, costbackorder, costprodcost, costholding, costransport, costopeningwarehouse, number_of_DV_RC_model)


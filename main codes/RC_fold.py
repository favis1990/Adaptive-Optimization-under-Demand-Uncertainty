# folding horizon 
import time
import numpy as np
#import matplotlib.pyplot as plt
from docplex.mp.model import Model
def RC_fold(carpi, jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, D_W, o_p , q, sigma, bigM, C_O_W, ini_inv, fold_X, fold_Y, fold_Z, fold_I, period, dis_per_unit):
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
    #initial inventory
    #ini_inv=mdl.continuous_var_dict([(w,j,k) for w in range(ww) for j in range(jj) for k in range(kk)], name='init')
     #exxtra baraye emhaye ezafe ha disposal
    exxtra_p=mdl.continuous_var_dict([(r,j,t)  for r in range(rr) for j in range(jj)  for t in range(tt)], name='extra_p')    
    #american komak be production
    exxtra_n=mdl.continuous_var_dict([(r,j,t)  for r in range(rr) for j in range(jj)  for t in range(tt)], name='extra_n')    

    #cost variables is used to seprate objective function for the sake of simplicity
    prodcost=mdl.continuous_var(name='prodcost')
    setupfirst=mdl.continuous_var(name='setupfirst')
    holding=mdl.continuous_var(name='holding')
    setup=mdl.continuous_var(name='setup')
    openingw=mdl.continuous_var(name='openingw')
    transport=mdl.continuous_var(name='transport')
    backorder=mdl.continuous_var(name='backorder')
    americancost=mdl.continuous_var(name='americancost')
    disposecost=mdl.continuous_var(name='disposecost')

    # Objective Function
    mdl.minimize(prodcost+setupfirst+holding+setup+transport+backorder+openingw+americancost+disposecost)
    mdl.add_constraint(americancost == mdl.sum(mdl.sum(mdl.sum(p[j]*carpi*exxtra_p[(r,j,t)]for r in range(rr))for j in range(jj))for t in range(tt)))
    mdl.add_constraint(disposecost == mdl.sum(mdl.sum(mdl.sum(dis_per_unit*exxtra_n[(r,j,t)]for r in range(rr))for j in range(jj))for t in range(tt)))
    mdl.add_constraint(prodcost == mdl.sum(mdl.sum( p[j]*Y[(j,t)]  for j in range(jj))for t in range(tt)))
    mdl.add_constraint(setupfirst == mdl.sum(mdl.sum(sigma[0][kprime] * Z_P[(w,kprime,0)]for w in range(ww))for kprime in range(kk)))
    mdl.add_constraint(holding == mdl.sum(mdl.sum(mdl.sum(mdl.sum(h[k][j]*I[(w,j,k,t)]for w in range(ww))for j in range(jj))for k in range(kk))for t in range(tt)))
    mdl.add_constraint(setup == mdl.sum(mdl.sum(mdl.sum(mdl.sum(sigma[k][kprime]*AX[(w,k,kprime,t)]for w in range(ww))for k in range(kk))for kprime in range(kk))for t in range(tt-1)))
    mdl.add_constraint(transport == mdl.sum(mdl.sum(mdl.sum(transperunit[j][r]*X[(r,j,t)] for j in range(jj)) for r in range(rr)) for t in range(tt)))
    mdl.add_constraint(backorder == mdl.sum(mdl.sum(mdl.sum(c_b[j]*B[(r,j,t)]for r in range(rr))for j in range(jj))for t in range(tt)))
    mdl.add_constraint(openingw == mdl.sum(mdl.sum(mdl.sum(Z_P[(w,k,t)]*C_O_W[(k)]for w in range(ww))for k in range(kk))for t in range(tt)))
    # Constraints
    #production capacity constraint
    mdl.add_constraints(Y[(j,t)] <= q_p[j][t]  for j in range(jj) for t in range(tt))
    #warehouse capacity constraint
    mdl.add_constraints(mdl.sum(mdl.sum(I[(w,j,k,t)]for j in range(jj))for k in range(kk))<=q_w[w][t]  for w in range(ww) for t in range(tt))
    #inventory balance constraint eliminated in ARC 
    mdl.add_constraints(Y[(j,0)] == mdl.sum(mdl.sum(I[(w,j,k,0)] - ini_inv[(w,j,k)] for k in  range (kk)) for w in range (ww)) + mdl.sum(X[(r,j,0)] for r in range (rr)) for j in range (jj))
    mdl.add_constraints(Y[(j,t)] == mdl.sum(mdl.sum(I[(w,j,k,t)] - I[(w,j,k,t-1)]  for k in  range (kk)) for w in range (ww)) + mdl.sum(X[(r,j,t)] for r in range (rr)) for j in range (jj)  for t in range (1,tt))
     #6 
    #inventory balance constraint2 backorder 
    mdl.add_constraints( X[(r,j,t-q[r])]+exxtra_p[(r,j,t-q[r])]-exxtra_n[(r,j,t-q[r])] == D_W[r][j][t]-B[(r,j,t)]+B[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t!=0))
    mdl.add_constraints( X[(r,j,t-q[r])]+exxtra_p[(r,j,t-q[r])]-exxtra_n[(r,j,t-q[r])]== D_W[r][j][t]-B[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t==0))
    mdl.add_constraints(0 == D_W[r][j][t]-B[(r,j,t)]+B[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t!=0))
    mdl.add_constraints(0 == D_W[r][j][t]-B[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t==0))
    #last period backorder   
    mdl.add_constraint(mdl.sum(mdl.sum(B[(r,j,tt-1)] for r in range(rr))for j in range(jj)) == 0 )
    #lead time cosntraint
    mdl.add_constraints(mdl.sum(X[(r,j,tp)]+exxtra_p[(r,j,tp)]-exxtra_n[(r,j,tp)] for tp in range(t+1))>= mdl.sum(D_W[r][j][tp] for tp in range(t+1))for r in range (rr) for j in range(jj) for t in range(tt))
    # inventory capacity 2
    mdl.add_constraints(mdl.sum(mdl.sum (I[(w,j,k,t)] for k in range (kk)) for w in range (ww)) <=CAP[j] for t in range(tt)for j in range(jj))
    #Inventory openning variable
    mdl.add_constraints( mdl.sum(I[(w,j,k,t)] for j in range(jj))<=bigM *Z_P[(w,k,t)] for w in range(ww) for k in range(kk) for t in range(tt))
    #Inventory openning variable2
    mdl.add_constraints(mdl.sum(Z_P[(w,k,t)] for k in range(kk))==1 for w in range(ww) for t in range(tt))
    #
    mdl.add_constraints(mdl.sum(mdl.sum(I[(w,j,k,t)]for w in range(ww))for t in range(tt))<= bigM *o_p[j][k] for j in range(jj)for k in range(kk))
    # to finish with some nventory 
    mdl.add_constraints(I[(w,j,k,tt-1)] == ini_inv[w][j][k] for w in range (ww) for j in range (jj)for k in range (kk))
    #folding constraint 
    mdl.add_constraints(Y[(j,t)]==fold_Y[(j,t)] for j in range (jj) for t in range (period))
    mdl.add_constraints(X[(r,j,t)]==fold_X[(r,j,t)] for r in range (rr) for j in range (jj) for t in range (period))
    mdl.add_constraints(I[(w,j,k,t)]==fold_I[(w,j,k,t)] for k in range (kk) for j in range (jj)for w in range (ww) for t in range (period))
    mdl.add_constraints(Z_P[(w,k,t)]==fold_Z[(w,k,t)] for w in range (ww) for k in range (kk) for t in range (period))
    #linearization of setup cost
    mdl.add_constraints(Z_P[(w,k,t)]+Z_P[(w,kprime,t+1)]-1<=AX[(w,k,kprime,t)]  for w in range(ww) for k in range(kk)for kprime in range(kk) for t in range(tt-1))
    mdl.add_constraints(AX[(w,k,kprime,t)]<=Z_P[(w,kprime,t+1)]  for w in range(ww) for k in range(kk)for kprime in range(kk) for t in range(tt-1))
    mdl.add_constraints(AX[(w,k,kprime,t)]<=Z_P[(w,k,t)]  for w in range(ww) for k in range(kk)for kprime in range(kk) for t in range(tt-1))
    #mdl.export('DeterministicMilk1')
    modelingRC = time.time()
    Z_P_V=np.zeros((ww,kk,tt))
    X_value=np.zeros((rr,jj,tt))
    I_value=np.zeros((ww,jj,kk,tt))
    y_value=np.zeros((jj,tt))
    B_value=np.zeros((rr,jj,tt))
    chooseRC = time.time()
    DETmodeling=modelingRC-RCstart
    DETtotal=chooseRC-modelingRC
    ini_inv_amount= np.zeros((ww,jj,kk))
    ex_valu_n=np.zeros((rr,jj,tt))
    ex_valu_p=np.zeros((rr,jj,tt))
     
       
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
               Z_P_V[w][k][t]=float(Z_P[(w,k,t)].solution_value)
    for r in range (rr) :
        for j in range (jj):
            for t in range (tt):
                B_value[r][j][t]=float(B[(r,j,t)].solution_value)
                X_value[r][j][t]=float(X[(r,j,t)].solution_value)
                ex_valu_n[r][j][t]=float(exxtra_n[(r,j,t)].solution_value)
                ex_valu_p[r][j][t]=float(exxtra_p[(r,j,t)].solution_value)
# =============================================================================
#     for w in range(ww):
#         for j in range (jj):
#             for k in range (kk):
#                ini_inv_amount[w][j][k]=float(ini_inv[(w,j,k)].solution_value)
# =============================================================================
                
# =============================================================================
#     total_Y = sum(sum( y_value[j][t] for j in range (jj))for t in range (tt))
#     total_X = sum(sum(sum( X_value[r][j][t] for j in range (jj)) for r in range (rr)) for t in range (tt))
#     total_B = sum(sum(sum( B_value[r][j][t] for j in range (jj)) for r in range (rr)) for t in range (tt))
#     total_I = sum(sum(sum(sum( I_value[w][j][k][t] for w in range (ww)) for j in range (jj)) for k in range (kk)) for t in range (tt))
#     totsl_initi= sum(sum(sum (ini_inv_amount[w][j][k] for w in range (ww))for j in range (jj))for k in range (kk))
# =============================================================================
    costsetupfirst = float(setupfirst)
    costsetupd = float(setup)
    costbackorder = float(backorder)
    costprodcost = float(prodcost)
    costholding = float(holding)
    costransport = float(transport)
    costopeningwarehouse = float(openingw)
    costdispose=float(disposecost)
    costamrican=float(americancost)
    number_of_DV_RC_model=(rr*jj*tt)+(jj*tt)+(ww*jj*kk*tt)+(rr*jj*tt)+(ww*kk*kk*tt)+(ww*kk*tt)
    print(number_of_DV_RC_model,number_of_DV_RC_model,number_of_DV_RC_model,number_of_DV_RC_model)
    #farzad=ateshia
    print(solution)    
# =============================================================================
# #first period
#     total_X_0_0 = X_value[0][0][0]
#     total_X_0_1 = X_value[0][1][0]
#     total_X_0_2 = X_value[0][2][0]
#     total_X_0_3 = X_value[0][3][0]
#     total_X_0_4 = X_value[0][4][0]
#     print('Hamadan',total_X_0_0, total_X_0_1, total_X_0_2, total_X_0_3, total_X_0_4)                    
#     total_X_1_0 = X_value[1][0][0]
#     total_X_1_1 = X_value[1][1][0]
#     total_X_1_2 = X_value[1][2][0]
#     total_X_1_3 = X_value[1][3][0]
#     total_X_1_4 = X_value[1][4][0]
#     print('Tehran', total_X_1_0, total_X_1_1, total_X_1_2, total_X_1_3, total_X_1_4) 
#     total_X_2_0 = X_value[2][0][0]
#     total_X_2_1 = X_value[2][1][0]
#     total_X_2_2 = X_value[2][2][0]
#     total_X_2_3 = X_value[2][3][0]
#     total_X_2_4 = X_value[2][4][0]
#     print('Khorasan',total_X_2_0, total_X_0_1, total_X_2_2, total_X_2_3, total_X_2_4) 
#     total_X_3_0 = X_value[3][0][0]
#     total_X_3_1 = X_value[3][1][0]
#     total_X_3_2 = X_value[3][2][0]
#     total_X_3_3 = X_value[3][3][0]
#     total_X_3_4 = X_value[3][4][0]
#     print('Birjand',total_X_3_0, total_X_3_1, total_X_3_2, total_X_3_3, total_X_3_4) 
#     total_X_4_0 = X_value[4][0][0]
#     total_X_4_1 = X_value[4][1][0]
#     total_X_4_2 = X_value[4][2][0]
#     total_X_4_3 = X_value[4][3][0]
#     total_X_4_4 = X_value[4][4][0]
#     print('Kordestan',total_X_4_0, total_X_4_1, total_X_4_2, total_X_4_3, total_X_4_4)                     
#     total_X_5_0 = X_value[0][0][0]
#     total_X_5_1 = X_value[5][1][0]
#     total_X_5_2 = X_value[5][2][0]
#     total_X_5_3 = X_value[5][3][0]
#     total_X_5_4 = X_value[5][4][0]
#     print('Isfahan',total_X_5_0, total_X_5_1, total_X_5_2, total_X_5_3, total_X_5_4) 
#     total_X_6_0 = X_value[0][0][0]
#     total_X_6_1 = X_value[6][1][0]
#     total_X_6_2 = X_value[6][2][0]
#     total_X_6_3 = X_value[6][3][0]
#     total_X_6_4 = X_value[6][4][0]
#     print('Qom',total_X_6_0, total_X_6_1, total_X_6_2, total_X_6_3, total_X_6_4) 
#     total_X_7_0 = X_value[0][0][0]
#     total_X_7_1 = X_value[7][1][0]
#     total_X_7_2 = X_value[7][2][0]
#     total_X_7_3 = X_value[7][3][0]
#     total_X_7_4 = X_value[7][4][0]
#     print('Khoramabad',total_X_7_0, total_X_7_1, total_X_7_2, total_X_7_3, total_X_7_4) 
# =============================================================================
    return(objdeter,ex_valu_p,ex_valu_n,Z_P_V,y_value,X_value,B_value,I_value,costdispose,costamrican,costsetupfirst,costsetupd,DETmodeling,DETtotal,costbackorder,costprodcost,costholding,costransport,costopeningwarehouse,number_of_DV_RC_model)

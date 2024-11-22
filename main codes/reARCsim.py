
import time
import numpy as np
#import matplotlib.pyplot as plt
from docplex.mp.model import Model
def reARCsim(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, q, D_W, setupf, setup, openning, I, alpha, X,beta, dis_per_unit, ini_inv):
    RCstart = time.time()
    ''' Model '''
    mdl=Model('DeterministicMilk')
    mdl.parameters.mip.strategy.file=3
    # Decision Variables
    #backorder
    B=mdl.continuous_var_dict([(r,j,t) for r in range(rr) for j in range(jj) for t in range(tt)], name='B')
    #openning Binary variable
    Z_P=mdl.binary_var_dict([(w,k,t) for w in range(ww) for k in range(kk) for t in range(tt)], name='Z_P')
    #cost variables is used to seprate objective function for the sake of simplicity
    # Objective Function
    prodcost = (sum(p[j]*((sum(sum(I[w][j][k][0]-ini_inv[w][j][k]+sum(sum(sum( alpha[w][j][k][0][rp][jp][tp]*D_W[rp][jp][tp]for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for w in range (ww))for k in range(kk)))+
    sum(X[r][j][0]+sum(sum(sum((beta[r][j][0][rp][jp][tp])*D_W[rp][jp][tp]for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for r in range (rr)))for j in range (jj)) + sum(sum(p[j]*((sum(sum(I[w][j][k][t]-I[w][j][k][t-1]+sum(sum(sum(((alpha[w][j][k][t][rp][jp][tp]-alpha[w][j][k][t-1][rp][jp][tp])*D_W[rp][jp][tp])for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for w in range (ww))for k in range(kk)))+
    sum(X[r][j][t]+sum(sum(sum((beta[r][j][t][rp][jp][tp]*D_W[rp][jp][tp])for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for r in range (rr)))for j in range (jj))for t in range (1,tt)) )
                                                                                                                                                                                                     
    holdingcost = sum(sum(sum(sum(h[k][j]*(I[w][j][k][t]+sum(sum(sum(alpha[w][j][k][t][rp][jp][tp]*D_W[rp][jp][tp] for rp in range (rr)) for jp in range (jj)) for tp in range (tt)))for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))
   
    distcost = sum(sum(sum(transperunit[j][r]*(X[r][j][t]+sum(sum(sum(beta[r][j][t][rp][jp][tp]*(D_W[rp][jp][tp]) for rp in range (rr)) for jp in range (jj)) for tp in range (tt)))for r in range (rr))for j in range (jj))for t in range (tt))    
   
    mdl.minimize(distcost+prodcost+setupf+setup+openning+mdl.sum(mdl.sum(mdl.sum(c_b[j]*B[(r,j,t)]for r in range (rr))for j in range (jj))for t in range (tt)))
    # Constraints

    #3 9latex 
    mdl.add_constraints(X[r][j][t-q[r]]+mdl.sum(mdl.sum(mdl.sum((beta[r][j][t-q[r]][rp][jp][tp]*D_W[rp][jp][tp]) for rp in range (rr))for jp in range (jj))for tp in range (tt)) >= D_W[r][j][t]-B[(r,j,t)]+B[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t!=0))
    mdl.add_constraints(X[r][j][t-q[r]]+mdl.sum(mdl.sum(mdl.sum((beta[r][j][t-q[r]][rp][jp][tp]*D_W[rp][jp][tp]) for rp in range (rr))for jp in range (jj))for tp in range (tt)) >= D_W[r][j][t]-B[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t==0))
    mdl.add_constraints(0>= D_W[r][j][t]-B[(r,j,t)]+B[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t!=0))
    mdl.add_constraints(0>= D_W[r][j][t]-B[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t==0))
    
   #4 10latex   
    mdl.add_constraint(mdl.sum(mdl.sum(B[(r,j,tt-1)] for r in range(rr))for j in range(jj))== 0 )
    #6 12latex mored dar

    mdl.export('DeterministicMilk')
    modelingRC = time.time()
    B_value=np.zeros((rr,jj,tt))
    chooseRC = time.time()
    DETmodeling=modelingRC-RCstart
    DETtotal=chooseRC-modelingRC
       
    try:
      solution = mdl.solve(log_output= True)
      feasible = 1
      OBJ = solution.get_objective_value()

                  
    except:
      feasible = 0
      obj=-125

# =============================================================================
#     total_Y=mdl.sum(mdl.sum(mdl.sum(mdl.sum( y_value[w][j][k][t]for w in range (ww))for j in range (jj))for t in range (tt))for k in range (kk))
#     total_X=mdl.sum(mdl.sum(mdl.sum(mdl.sum(mdl.sum( X_value[w][r][j][k][t]for w in range (ww))for j in range (jj))for t in range (tt))for k in range (kk))for r in range (rr))
#     total_I=mdl.sum(mdl.sum(mdl.sum(mdl.sum(I_value[w][j][k][t]for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))
#     final_b=mdl.sum(mdl.sum(B_value[r][j][tt-1]for r in range (rr))for j in range (jj))
#     real_demand=mdl.sum(mdl.sum(mdl.sum(D_W[r][j][t]for r in range (rr))for j in range (jj)) for t in range (tt))
# =============================================================================
    print(solution)    
    X_real=np.zeros((rr,jj,tt))
    for j in range (jj):
        for t in range (tt):
            for r in range (rr):
                X_real[r][j][t]=X[r][j][t]+sum(sum(sum(beta[r][j][t][rp][jp][tp]*D_W[rp][jp][tp]for rp in range (rr))for jp in range (jj))for tp in range (tt))    
    exxtra_value_2=np.sum(np.sum(np.sum(X_real)))-np.sum(np.sum(np.sum(D_W)))
    for r in range (rr):
        for j in range (jj):
            for t in range (tt):
                B_value[r][j][t]=float(B[(r,j,t)])
                
    costbackorder = mdl.sum(mdl.sum(mdl.sum(c_b[j]* B_value[r][j][t]for r in range (rr))for j in range (jj))for t in range (tt))
    disp_cost=exxtra_value_2*dis_per_unit

    return(exxtra_value_2,X_real,OBJ,disp_cost,prodcost,distcost,costbackorder,setupf,setup,holdingcost,openning)


#sol=reARCsim()
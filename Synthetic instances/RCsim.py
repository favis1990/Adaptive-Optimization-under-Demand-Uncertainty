
import time
import numpy as np
#import matplotlib.pyplot as plt
from docplex.mp.model import Model
def RCsim(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, D_W, q, obj, backordercost, X, dis_per_unit):
    RCstart = time.time()
    ''' Model '''
    mdl=Model('DeterministicMilk')
    mdl.parameters.mip.strategy.file=3
    # Decision Variables
    #material flow    
    B=mdl.continuous_var_dict([(r,j,t) for r in range(rr) for j in range(jj) for t in range(tt)], name='B')
    # Auxiliary variable for linearization os setup cost
    backorder=mdl.continuous_var(name='backorder')
    # Objective Function
    mdl.minimize(mdl.sum(mdl.sum(mdl.sum(c_b[j]*B[(r,j,t)]for r in range(rr))for j in range(jj))for t in range(tt)))
    # Constraints
    mdl.add_constraints(X[(r,j,t-q[r])] >= D_W[r][j][t]-B[(r,j,t)]+B[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t!=0))
    mdl.add_constraints(X[(r,j,t-q[r])] >= D_W[r][j][t]-B[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]>=0 and t==0))
    mdl.add_constraints(0>= D_W[r][j][t]-B[(r,j,t)]+B[(r,j,t-1)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t!=0))
    mdl.add_constraints(0>= D_W[r][j][t]-B[(r,j,t)] for r in range (rr) for j in range(jj) for t in range(tt) if (t-q[r]<0 and t==0))
    
    mdl.add_constraint(mdl.sum(mdl.sum(B[(r,j,tt-1)] for r in range(rr))for j in range(jj))== 0 )
    
    B_value=np.zeros((rr,jj,tt))
    chooseRC = time.time()
       
    try:
      solution = mdl.solve(log_output= True)
      feasible = 1
      newbackordercost = solution.get_objective_value()
                  
    except:
      feasible = 0
      obj=-125

    for r in range (rr) :
        for j in range (jj):
            for t in range (tt):
                B_value[r][j][t]=float(B[(r,j,t)].solution_value)
                
    costbackorder=float(backorder)
    
# =============================================================================
#     total_Y=mdl.sum(mdl.sum(mdl.sum(mdl.sum( y_value[w][j][k][t]for w in range (ww))for j in range (jj))for t in range (tt))for k in range (kk))
#     total_X=mdl.sum(mdl.sum(mdl.sum(mdl.sum(mdl.sum( X_value[w][r][j][k][t]for w in range (ww))for j in range (jj))for t in range (tt))for k in range (kk))for r in range (rr))
#     total_I=mdl.sum(mdl.sum(mdl.sum(mdl.sum(I_value[w][j][k][t]for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))
#     final_b=mdl.sum(mdl.sum(B_value[r][j][tt-1]for r in range (rr))for j in range (jj))
#     real_demand=mdl.sum(mdl.sum(mdl.sum(D_W[r][j][t]for r in range (rr))for j in range (jj)) for t in range (tt))
# =============================================================================
    print(solution)
    exxtra_value_rc=np.sum(np.sum(np.sum(X)))-np.sum(np.sum(np.sum(D_W)))
    costdispose_c_rc=dis_per_unit*exxtra_value_rc
    newobj=obj-backordercost+newbackordercost+costdispose_c_rc
    return(newobj,exxtra_value_rc,costdispose_c_rc,newbackordercost)

#sol=RC()
#ARC milk job functon using fixed binary Z_P_V 
#[ARPDP](RPDP)

import time 
import numpy as np
#import matplotlib.pyplot as plt
from docplex.mp.model import Model

def ARC(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, D_W, dz, dh, o_w, q, sigma, bigM, openingPN,  setupfirstN,     costopenN,   Z_P_V, ini_inv) :    
    ARCstart = time.time()
    ''' Model '''
    axi=np.zeros((rr,jj,tt,rr,jj,tt))
    axi2=np.zeros((rr,jj,tt,tt,rr,jj,tt))
    for r in range(rr):
        for j in range(jj):
             for t in range(tt):              
                 for rp in range(rr):
                     for jp in range(jj):
                         for tp in range(tt):
                            if rp==r and jp==j and tp==t:
                                axi[r][j][t][rp][jp][tp]=1                                
                            else:
                                 axi[r][j][t][rp][jp][tp]=0                                                                  
                            for tz in range(tt):  
                                if rp==r and jp==j and tp==t and tz==t:
                                     axi2[r][j][t][tp][rp][jp][tz]=1
                                else:
                                     axi2[r][j][t][tp][rp][jp][tz]=0
    mdl=Model('ARCMilk')  
    # Decision Variables
    #material flow
    X_bar=mdl.continuous_var_dict([(r,j,t) for r in range(rr) for j in range(jj) for t in range(tt)], name='x_bar')
    beta=mdl.continuous_var_dict([(r,j,t,rp,jp,tz) for r in range(rr) for j in range(jj) for t in range(tt) for rp in range(rr) for jp in range(jj) for tz in range(tt)],lb=-float('inf'),  name='beta')
    #production amount
    I_bar=mdl.continuous_var_dict([(w,j,k,t) for w in range(ww) for j in range(jj) for k in range(kk) for t in range(tt)],lb=-float('inf'), name='I_bar')
    alpha=mdl.continuous_var_dict([(w,j,k,t,r,jp,tp) for w in range(ww) for j in range(jj) for k in range(kk) for t in range(tt) for r in range(rr) for jp in range(jj) for tp in range(tt)],lb=-float('inf'), name='alpha')
    #backorder
    B_bar=mdl.continuous_var_dict([(r,j,t) for r in range(rr) for j in range(jj) for t in range(tt)],lb=-float('inf'), name='B_bar')
    gamma=mdl.continuous_var_dict([(r,j,t,rp,jp,tp) for r in range(rr) for j in range(jj) for t in range(tt) for rp in range(rr) for jp in range(jj) for tp in range(tt)],lb=-float('inf'),  name='gamma')
    #cost variables is used to seprate objective function for the sake of simplicity
    OTHR=mdl.continuous_var(name='OTHR')
    #auxiliries variables for roboust counterpart
    AXOBJ=mdl.continuous_var_dict([(rp,jp,tp) for rp in range(rr) for jp in range(jj) for tp in range(tt)], name='AXOBJ')
    ax1  = mdl.continuous_var_dict([(j,t,r,jp,tp) for j in range(jj) for t in range(tt) for r in range(rr) for jp in range(jj) for tp in range(tt)], name='ax1')
    ax2  = mdl.continuous_var_dict([(w,t,rp,jp,tz) for w in range(ww)for t in range(tt) for rp in range(rr) for jp in range(jj) for tz in range(tt)], name='ax2')
    ax3  = mdl.continuous_var_dict([(r,j,t,rp,jp,tp)for r in range(rr) for j in range(jj) for t in range(tt) for jp in range(jj) for rp in range(rr)for tp in range(tt)], name='ax3')    
    ax4  = mdl.continuous_var_dict([(r,j,t,rp,jp,tz,tp)for r in range(rr)for j in range(jj) for t in range(tt) for rp in range(rr) for jp in range(jj) for tz in range(tt)for tp in range(tt)], name='ax4')
    ax5  = mdl.continuous_var_dict([(t,j,rp,jp,tp)  for t in range(tt)for j in range(jj) for rp in range(rr) for jp in range(jj) for tp in range(tt)], name='ax5')
    ax6  = mdl.continuous_var_dict([(w,k,t,rp,jp,tp) for w in range (ww) for k in range (kk)for t in range(tt) for jp in range(jj) for rp in range(rr)for tp in range(tt)], name='ax6')
    ax7  = mdl.continuous_var_dict([(r,j,t,rp,jp,tp) for r in range (rr)for j in range (jj) for t in range(tt) for tp in range(tt) for jp in range(jj) for rp in range(rr)], name='ax7')
    ax8  = mdl.continuous_var_dict([(w,j,k,t,rp,jp,tp) for w in range (ww) for j in range (jj) for k in range (kk)for t in range(tt) for jp in range(jj) for rp in range(rr)for tp in range(tt)], name='ax8')
    ax9  = mdl.continuous_var_dict([(r,j,t,rp,jp,tp)for r in range(rr) for j in range(jj) for t in range(tt) for jp in range(jj) for rp in range(rr)for tp in range(tt)], name='ax9')
    ax10 = mdl.continuous_var_dict([(j,t,rp,jp,tp) for j in range(jj) for t in range(tt) for jp in range(jj) for rp in range(rr)for tp in range(tt)], name='ax10')
    ax11 = mdl.continuous_var_dict([(r,j,t,rp,jp,tp)for r in range(rr) for j in range(jj) for t in range(tt)for jp in range(jj) for rp in range(rr)for tp in range(tt)], name='ax11')
    ax12 = mdl.continuous_var_dict([(j,k,rp,jp,tp)for j in range(jj) for k in range(kk) for jp in range(jj) for rp in range(rr)for tp in range(tt)], name='ax12')
        
     # Objective Function
    mdl.minimize(openingPN+setupfirstN+costopenN+OTHR)
    mdl.add_constraint(OTHR == mdl.sum(p[j]*((mdl.sum(mdl.sum(I_bar[(w,j,k,0)]-ini_inv[w][j][k]+mdl.sum(mdl.sum(mdl.sum(alpha[(w,j,k,0,rp,jp,tp)]*dz[rp][jp][tp]for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for w in range (ww))for k in range (kk)))
    + mdl.sum(X_bar[(r,j,0)]+mdl.sum(mdl.sum(mdl.sum((beta[(r,j,0,rp,jp,tp)])*dz[rp][jp][tp]for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for r in range (rr)))for j in range (jj))
    + mdl.sum(mdl.sum(p[j]*((mdl.sum(mdl.sum(I_bar[(w,j,k,t)]-I_bar[(w,j,k,t-1)]+mdl.sum(mdl.sum(mdl.sum(((alpha[(w,j,k,t,rp,jp,tp)]-alpha[(w,j,k,t-1,rp,jp,tp)])*dz[rp][jp][tp])for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for w in range (ww))for k in range (kk)))
    + mdl.sum(X_bar[(r,j,t)]+mdl.sum(mdl.sum(mdl.sum((beta[(r,j,t,rp,jp,tp)]*dz[rp][jp][tp])for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for r in range (rr)))for j in range (jj))for t in range (1,tt))
    + mdl.sum(mdl.sum(mdl.sum(mdl.sum( h[k][j]*(I_bar[(w,j,k,t)]+mdl.sum(mdl.sum(mdl.sum(alpha[(w,j,k,t,rp,jp,tp)]*dz[rp][jp][tp]for rp in range (rr)) for jp in range (jj)) for tp in range (tt)))for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))
    + mdl.sum(mdl.sum(mdl.sum( transperunit[j][r] * (X_bar[(r,j,t)]+mdl.sum(mdl.sum(mdl.sum(beta[(r,j,t,rp,jp,tp)]*dz[rp][jp][tp]for rp in range (rr)) for jp in range (jj)) for tp in range (tt)))for r in range (rr)) for j in range (jj)) for t in range (tt))
    + mdl.sum(mdl.sum(mdl.sum( c_b[j] *(B_bar[(r,j,t)]+mdl.sum(mdl.sum(mdl.sum(gamma[(r,j,t,rp,jp,tp)]*dz[rp][jp][tp] for rp in range (rr)) for jp in range (jj)) for tp in range (tt)))for r in range (rr))for j in range (jj))for t in range (tt))
    + mdl.sum(mdl.sum(mdl.sum(AXOBJ[(rp,jp,tp)]*dh[rp][jp][tp]for rp in range (rr))for jp in range (jj))for tp in range (tt)))
    
    
    
    mdl.add_constraints(AXOBJ[(rp,jp,tz)] >= mdl.sum(p[j]*(mdl.sum(mdl.sum(alpha[(w,j,k,0,rp,jp,tz)]for w in range (ww))for k in range (kk))
    + mdl.sum(beta[(r,j,0,rp,jp,tz)] for r in range (rr)))for j in range (jj))
    + mdl.sum(mdl.sum(p[j]*(mdl.sum(mdl.sum((alpha[(w,j,k,t,rp,jp,tz)]-alpha[(w,j,k,t-1,rp,jp,tz)])for w in range (ww))for k in range (kk))
    + mdl.sum(beta[(r,j,t,rp,jp,tz)]for r in range (rr)))for j in range (jj))for t in range (1,tt))
    + mdl.sum(mdl.sum(mdl.sum(mdl.sum( h[k][j]*alpha[(w,j,k,t,rp,jp,tz)]for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))
    + mdl.sum(mdl.sum(mdl.sum( transperunit[j][r] * beta[(r,j,t,rp,jp,tz)]for r in range (rr)) for j in range (jj)) for t in range (tt))
    + mdl.sum(mdl.sum(mdl.sum(c_b[j]*gamma[(r,j,t,rp,jp,tz)] for r in range (rr)) for j in range (jj)) for t in range (tt)) for rp in range (rr)for jp in range (jj)for tz in range (tt))
    
    mdl.add_constraints(AXOBJ[(rp,jp,tz)] >= -(mdl.sum(p[j]*(mdl.sum(mdl.sum(alpha[(w,j,k,0,rp,jp,tz)]for w in range (ww))for k in range (kk))
    + mdl.sum(beta[(r,j,0,rp,jp,tz)] for r in range (rr)))for j in range (jj))
    + mdl.sum(mdl.sum(p[j]*(mdl.sum(mdl.sum((alpha[(w,j,k,t,rp,jp,tz)]-alpha[(w,j,k,t-1,rp,jp,tz)])for w in range (ww))for k in  range (kk))
    + mdl.sum(beta[(r,j,t,rp,jp,tz)]for r in range (rr)))for j in range (jj))for t in range (1,tt))
    + mdl.sum(mdl.sum(mdl.sum(mdl.sum( h[k][j]*alpha[(w,j,k,t,rp,jp,tz)]for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))
    + mdl.sum(mdl.sum(mdl.sum( transperunit[j][r] * beta[(r,j,t,rp,jp,tz)]for r in range (rr)) for j in range (jj)) for t in range (tt))
    + mdl.sum(mdl.sum(mdl.sum(c_b[j]*gamma[(r,j,t,rp,jp,tz)] for r in range (rr)) for j in range (jj)) for t in range (tt))) for rp in range (rr)for jp in range (jj)for tz in range (tt))
  
    # Constraints
   
    # #production capacity constraint #1
    mdl.add_constraints(mdl.sum(mdl.sum((I_bar[(w,j,k,0)]-ini_inv[w][j][k] +mdl.sum(mdl.sum(mdl.sum((alpha[(w,j,k,0,rp,jp,tp)]*dz[rp][jp][tp]) for rp in range (rr))for jp in range (jj))for tp in range (tt))) for k in range(kk)) for w in range (ww))
    + (mdl.sum(X_bar[(r,j,0)]+mdl.sum(mdl.sum(mdl.sum(beta[(r,j,0,rp,jp,tp)]*dz[rp][jp][tp] for rp in range (rr))for jp in range (jj))for tp in range (tt))for r in range (rr))) 
    + mdl.sum(mdl.sum(mdl.sum(ax1[(j,0,rp,jp,tp)]*dh[rp][jp][tp] for rp in range (rr))for jp in range (jj))for tp in range (tt))<= q_p[j][0] for j in range (jj) )    
    mdl.add_constraints(mdl.sum(mdl.sum((I_bar[(w,j,k,t)]-I_bar[(w,j,k,t-1)]+mdl.sum(mdl.sum(mdl.sum(((alpha[(w,j,k,t,rp,jp,tp)]-alpha[(w,j,k,t-1,rp,jp,tp)])*dz[rp][jp][tp]) for rp in range (rr))for jp in range (jj))for tp in range (tt))) for k in range(kk)) for w in range (ww))
    + (mdl.sum(X_bar[(r,j,t)]+mdl.sum(mdl.sum(mdl.sum(beta[(r,j,t,rp,jp,tp)]*dz[rp][jp][tp] for rp in range (rr))for jp in range (jj))for tp in range (tt))for r in range (rr)))
    + mdl.sum(mdl.sum(mdl.sum(ax1[(j,t,rp,jp,tp)]*dh[rp][jp][tp] for rp in range (rr))for jp in range (jj))for tp in range (tt))<= q_p[j][t] for j in range (jj) for t in range (1,tt))    
    
    mdl.add_constraints(ax1[(j,0,rp,jp,tp)]>=   mdl.sum(mdl.sum( alpha[(w,j,k,0,rp,jp,tp)] for w in range (ww)) for k in range(kk)) + mdl.sum(beta[(r,j,0,rp,jp,tp)] for r in range (rr)) for j in range (jj)  for jp in range (jj)for tp in range (tt)for rp in range (rr))
    mdl.add_constraints(ax1[(j,0,rp,jp,tp)]>= - (mdl.sum(mdl.sum( alpha[(w,j,k,0,rp,jp,tp)] for w in range (ww)) for k in range(kk)) + mdl.sum(beta[(r,j,0,rp,jp,tp)] for r in range (rr))) for j in range (jj)  for jp in range (jj)for tp in range (tt)for rp in range (rr))
    
    mdl.add_constraints(ax1[(j,t,rp,jp,tp)]>=   mdl.sum(mdl.sum( (alpha[(w,j,k,t,rp,jp,tp)]-alpha[(w,j,k,t-1,rp,jp,tp)]) for w in range (ww)) for k in range(kk)) + mdl.sum(beta[(r,j,t,rp,jp,tp)] for r in range (rr)) for j in range (jj) for t in range (1,tt) for jp in range (jj)for tp in range (tt)for rp in range (rr))
    mdl.add_constraints(ax1[(j,t,rp,jp,tp)]>= - (mdl.sum(mdl.sum( (alpha[(w,j,k,t,rp,jp,tp)]-alpha[(w,j,k,t-1,rp,jp,tp)]) for w in range (ww)) for k in range(kk)) + mdl.sum(beta[(r,j,t,rp,jp,tp)] for r in range (rr))) for j in range (jj) for t in range (1,tt) for jp in range (jj)for tp in range (tt)for rp in range (rr))
    

    #  #warehouse capacity constraint #2
    mdl.add_constraints(mdl.sum(mdl.sum(I_bar[(w,j,k,t)]+mdl.sum(mdl.sum(mdl.sum(alpha[(w,j,k,t,r,jp,tp)]*dz[rp][jp][tp]for tp in range (tt))for jp in range (jj))for rp in range (rr))for k in range(kk))for j in range(jj))
    +mdl.sum(mdl.sum(mdl.sum(ax2[(w,t,rp,jp,tp)]*dh[rp][jp][tp]for tp in range (tt))for jp in range (jj))for rp in range (rr))<=q_w[w][t]  for w in range(ww) for t in range(tt))
    mdl.add_constraints(ax2[(w,t,rp,jp,tp)]>= mdl.sum(mdl.sum(alpha[(w,j,k,t,rp,jp,tp)] for j in range (jj)) for k in range (kk))  for w in range (ww)for t in range (tt) for tp in range (tt) for jp in range (jj)for rp in range (rr))
    mdl.add_constraints(ax2[(w,t,rp,jp,tp)]>=-mdl.sum(mdl.sum(alpha[(w,j,k,t,rp,jp,tp)] for j in range (jj)) for k in range (kk))  for w in range (ww)for t in range (tt) for tp in range (tt) for jp in range (jj)for rp in range (rr))
    
   
    #inventory balance constraint2 backorder 
    mdl.add_constraints(X_bar[(r,j,t-q[r])]+B_bar[(r,j,t)]-B_bar[(r,j,t-1)]+ mdl.sum(mdl.sum(mdl.sum(beta[(r,j,t-q[r],rp,jp,tp)]*dz[rp][jp][tp]
    -gamma[(r,j,t-1,rp,jp,tp)]*dz[rp][jp][tp]+gamma[(r,j,t,rp,jp,tp)]*dz[rp][jp][tp]-ax3[(r,j,t,rp,jp,tp)]*dh[rp][jp][tp] for rp in range (rr)) for jp in range (jj)) for tp in range (tt))-dz[r][j][t] >=0 for r in range (rr)for j in range (jj)for t in range (tt) if (t-q[r]>=0 and t!=0))
    mdl.add_constraints(ax3[(r,j,t,rp,jp,tp)]>= beta[(r,j,t-q[r],rp,jp,tp)] +gamma[(r,j,t,rp,jp,tp)]-gamma[(r,j,t-1,rp,jp,tp)]-axi[r][j][t][rp][jp][tp] for r in range (rr)for j in range (jj)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt)if (t-q[r]>=0 and t!=0))
    mdl.add_constraints(ax3[(r,j,t,rp,jp,tp)]>=-(beta[(r,j,t-q[r],rp,jp,tp)] +gamma[(r,j,t,rp,jp,tp)]-gamma[(r,j,t-1,rp,jp,tp)]-axi[r][j][t][rp][jp][tp]) for r in range (rr)for j in range (jj)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt)if (t-q[r]>=0 and t!=0))

    mdl.add_constraints(X_bar[(r,j,t-q[r])]+B_bar[(r,j,t)]+ mdl.sum(mdl.sum(mdl.sum(beta[(r,j,t-q[r],rp,jp,tp)]*dz[rp][jp][tp]
    +gamma[(r,j,t,rp,jp,tp)]*dz[rp][jp][tp]-ax3[(r,j,t,rp,jp,tp)]*dh[rp][jp][tp] for rp in range (rr)) for jp in range (jj)) for tp in range (tt))-dz[r][j][t] >=0 for r in range (rr)for j in range (jj)for t in range (tt) if (t-q[r]>=0 and t==0))
    mdl.add_constraints(ax3[(r,j,t,rp,jp,tp)]>=  beta[(r,j,t-q[r],rp,jp,tp)]+gamma[(r,j,t,rp,jp,tp)]-axi[r][j][t][rp][jp][tp]for r in range (rr)for j in range (jj)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt)if (t-q[r]>=0 and t==0))
    mdl.add_constraints(ax3[(r,j,t,rp,jp,tp)]>= -(beta[(r,j,t-q[r],rp,jp,tp)]+gamma[(r,j,t,rp,jp,tp)]-axi[r][j][t][rp][jp][tp]) for r in range (rr)for j in range (jj)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt)if (t-q[r]>=0 and t==0))
    
    mdl.add_constraints(B_bar[(r,j,t)]-B_bar[(r,j,t-1)]+mdl.sum(mdl.sum(mdl.sum(gamma[(r,j,t,rp,jp,tp)]*dz[rp][jp][tp]-gamma[(r,j,t-1,rp,jp,tp)]*dz[rp][jp][tp]-ax3[(r,j,t,rp,jp,tp)]*dh[rp][jp][tp] for rp in range (rr)) for jp in range (jj)) for tp in range (tt))-dz[r][j][t] >=0 for r in range (rr)for j in range (jj)for t in range (tt) if (t-q[r]<0 and t!=0))
    mdl.add_constraints(ax3[(r,j,t,rp,jp,tp)]>=gamma[(r,j,t,rp,jp,tp)]-gamma[(r,j,t-1,rp,jp,tp)]-axi[r][j][t][rp][jp][tp] for r in range (rr)for j in range (jj)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt)if (t-q[r]<0 and t!=0))
    mdl.add_constraints(ax3[(r,j,t,rp,jp,tp)]>=-(gamma[(r,j,t,rp,jp,tp)]-gamma[(r,j,t-1,rp,jp,tp)]-axi[r][j][t][rp][jp][tp]) for r in range (rr)for j in range (jj)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt)if (t-q[r]<0 and t!=0))
    
    mdl.add_constraints(B_bar[(r,j,t)]+ mdl.sum(mdl.sum(mdl.sum(gamma[(r,j,t,rp,jp,tp)]*dz[rp][jp][tp]-ax3[(r,j,t,rp,jp,tp)]*dh[rp][jp][tp] for rp in range (rr)) for jp in range (jj)) for tp in range (tt))-dz[r][j][t] >=0 for r in range (rr)for j in range (jj)for t in range (tt) if (t-q[r]<0 and t==0))
    mdl.add_constraints(ax3[(r,j,t,rp,jp,tp)]>=gamma[(r,j,t,rp,jp,tp)]-axi[r][j][t][rp][jp][tp] for r in range (rr)for j in range (jj)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt)if (t-q[r]<0 and t==0))
    mdl.add_constraints(ax3[(r,j,t,rp,jp,tp)]>=-(gamma[(r,j,t,rp,jp,tp)]-axi[r][j][t][rp][jp][tp]) for r in range (rr)for j in range (jj)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt)if (t-q[r]<0 and t==0))
    
   #lead time cosntraint
    mdl.add_constraints(mdl.sum((X_bar[(r,j,tp)]+mdl.sum(mdl.sum(mdl.sum((beta[(r,j,tp,rp,jp,tz)]*dz[rp][jp][tz])for rp in range (rr))for jp in range (jj))for tz in range (tp+1)))for tp in range(t+1))
    -dz[r][j][t]-mdl.sum(mdl.sum(mdl.sum(mdl.sum(ax4[(r,j,t,rp,jp,tz,tp)]*dh[rp][jp][tz]for jp in range (jj))for rp in range (rr))for tz in range (tp+1))for tp in range (t+1))>= 0 for r in range (rr)for j in range (jj)for t in range (tt))    
    mdl.add_constraints(ax4[(r,j,t,rp,jp,tz,tp)]>= beta[(r,j,tp,rp,jp,tz)]-axi2[r][j][t][tp][rp][jp][tz]for r in range (rr)for j in range (jj)for t in range (tt)for tp in range (tt)for tz in range (tt)for jp in range (jj) for rp in range (rr))
    mdl.add_constraints(ax4[(r,j,t,rp,jp,tz,tp)]>= -(beta[(r,j,tp,rp,jp,tz)]-axi2[r][j][t][tp][rp][jp][tz])for r in range (rr)for j in range (jj)for t in range (tt)for tp in range (tt)for tz in range (tt)for jp in range (jj) for rp in range (rr))
    
    #last period backorder   
    mdl.add_constraints(B_bar[(r,j,tt-1)] == 0 for r in range (rr)for j in range (jj))
    mdl.add_constraints(gamma[(r,j,tt-1,rp,jp,tp)]==0 for rp in range (rr)for jp in range (jj)for tp in range (tt)for r in range (rr)for j in range (jj))
    
    # inventory capacity 2
    mdl.add_constraints(mdl.sum(mdl.sum(I_bar[(w,j,k,t)]+mdl.sum(mdl.sum(mdl.sum(alpha[(w,j,k,t,rp,jp,tp)]*dz[rp][jp][tp]for tp in range (tt))for rp in range (rr))for jp in range (jj))for k in range(kk))for w in range(ww))
    +mdl.sum(mdl.sum(mdl.sum(ax5[(t,j,rp,jp,tp)]*dh[rp][jp][tp]for tp in range (tt))for jp in range (jj))for rp in range (rr))
    <=CAP[j]  for t in range(tt) for j in range(jj))
    mdl.add_constraints(ax5[(t,j,rp,jp,tp)]>= mdl.sum(mdl.sum(alpha[(w,j,k,t,rp,jp,tp)]  for k in range (kk))  for w in range (ww)) for j in range(jj)for t in range (tt) for tp in range (tt) for jp in range (jj)for rp in range (rr))
    mdl.add_constraints(ax5[(t,j,rp,jp,tp)]>=-mdl.sum(mdl.sum(alpha[(w,j,k,t,rp,jp,tp)]  for k in range (kk))  for w in range (ww)) for j in range(jj)for t in range (tt) for tp in range (tt) for jp in range (jj)for rp in range (rr))

    # # #Inventory openning #8
    mdl.add_constraints(mdl.sum(I_bar[(w,j,k,t)]+mdl.sum(mdl.sum(mdl.sum(alpha[(w,j,k,t,r,jp,tp)]*dz[rp][jp][tp]for tp in range (tt))for jp in range (jj))for j in range(jj))for j in range (jj))
    +mdl.sum(mdl.sum(mdl.sum(ax6[(w,k,t,rp,jp,tp)]*dh[rp][jp][tp]for tp in range (tt))for jp in range (jj))for rp in range (rr))
    <=bigM *Z_P_V[(w,k,t)]  for w in range(ww)for k in range(kk)for t in range(tt))
    
    
    
    mdl.add_constraints(ax6[(w,k,t,rp,jp,tp)]>= mdl.sum(alpha[(w,j,k,t,rp,jp,tp)] for j in range (jj)) for k in range (kk)  for w in range (ww)for t in range (tt) for tp in range (tt) for jp in range (jj)for rp in range (rr))
    mdl.add_constraints(ax6[(w,k,t,rp,jp,tp)]>=-mdl.sum(alpha[(w,j,k,t,rp,jp,tp)] for j in range (jj)) for k in range (kk)  for w in range (ww)for t in range (tt) for tp in range (tt) for jp in range (jj)for rp in range (rr))

    # to finish with initial inventory #10 
    mdl.add_constraints(I_bar[(w,j,k,tt-1)] == ini_inv[w][j][k] for w in range (ww) for j in range (jj)for k in range (kk))
    mdl.add_constraints(alpha[(w,j,k,tt-1,rp,jp,tp)]==0 for rp in range (rr)for jp in range (jj)for tp in range (tt)for w in range (ww) for j in range (jj)for k in range (kk))
    

    #non negtivity X
    mdl.add_constraints(X_bar[(r,j,t)]+mdl.sum(mdl.sum(mdl.sum((beta[(r,j,t,rp,jp,tp)]*dz[rp][jp][tp]-ax7[(r,j,t,rp,jp,tp)]*dh[rp][jp][tp]) for rp in range (rr))for jp in range (jj))for tp in range (tt))>=0 for w in range (ww)for r in range (rr)for j in range (jj)for k in range (kk)for t in range (tt))
    mdl.add_constraints(ax7[(r,j,t,rp,jp,tp)] >= beta[(r,j,t,rp,jp,tp)]for r in range (rr)for j in range (jj)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt))
    mdl.add_constraints(ax7[(r,j,t,rp,jp,tp)] >= -beta[(r,j,t,rp,jp,tp)]for r in range (rr)for j in range (jj)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt))
    
    mdl.add_constraints(I_bar[(w,j,k,t)]+mdl.sum(mdl.sum(mdl.sum((alpha[(w,j,k,t,rp,jp,tp)]*dz[rp][jp][tp]-ax8[(w,j,k,t,rp,jp,tp)]*dh[rp][jp][tp]) for rp in range (rr))for jp in range (jj))for tp in range (tt))>=0 for w in range (ww)for j in range (jj)for k in range (kk)for t in range (tt))
    mdl.add_constraints(ax8[(w,j,k,t,rp,jp,tp)]>=alpha[(w,j,k,t,rp,jp,tp)]for w in range (ww)for j in range (jj)for k in range (kk)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt))
    mdl.add_constraints(ax8[(w,j,k,t,rp,jp,tp)]>=-alpha[(w,j,k,t,rp,jp,tp)]for w in range (ww)for j in range (jj)for k in range (kk)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt))
    #non negtivity I
    mdl.add_constraints(B_bar[(r,j,t)]+mdl.sum(mdl.sum(mdl.sum((gamma[(r,j,t,rp,jp,tp)]*dz[rp][jp][tp]-ax9[(r,j,t,rp,jp,tp)]*dh[rp][jp][tp]) for rp in range (rr))for jp in range (jj))for tp in range (tt))>=0 for r in range (rr)for j in range  (jj)for t in range (tt))
    mdl.add_constraints(ax9[(r,j,t,rp,jp,tp)]>=gamma[(r,j,t,rp,jp,tp)]for j in range (jj)for r in range (rr)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt))
    mdl.add_constraints(ax9[(r,j,t,rp,jp,tp)]>=-gamma[(r,j,t,rp,jp,tp)]for j in range (jj)for r in range (rr)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt))
    #non negtivity Y
    mdl.add_constraints(mdl.sum(mdl.sum((I_bar[(w,j,k,0)]-ini_inv[w][j][k] +mdl.sum(mdl.sum(mdl.sum((alpha[(w,j,k,0,rp,jp,tp)]*dz[rp][jp][tp]) for rp in range (rr))for jp in range (jj))for tp in range (tt))) for k in range(kk)) for w in range (ww))
    + (mdl.sum(X_bar[(r,j,0)]+mdl.sum(mdl.sum(mdl.sum(beta[(r,j,0,rp,jp,tp)]*dz[rp][jp][tp] for rp in range (rr))for jp in range (jj))for tp in range (tt))for r in range (rr))) 
    - mdl.sum(mdl.sum(mdl.sum(ax10[(j,0,rp,jp,tp)]*dh[rp][jp][tp] for rp in range (rr))for jp in range (jj))for tp in range (tt))>= 0 for j in range (jj) )    

    mdl.add_constraints(mdl.sum(mdl.sum((I_bar[(w,j,k,t)]-I_bar[(w,j,k,t-1)]+mdl.sum(mdl.sum(mdl.sum(((alpha[(w,j,k,t,rp,jp,tp)]-alpha[(w,j,k,t-1,rp,jp,tp)])*dz[rp][jp][tp]) for rp in range (rr))for jp in range (jj))for tp in range (tt))) for k in range(kk)) for w in range (ww))
    + (mdl.sum(X_bar[(r,j,t)]+mdl.sum(mdl.sum(mdl.sum(beta[(r,j,t,rp,jp,tp)]*dz[rp][jp][tp] for rp in range (rr))for jp in range (jj))for tp in range (tt))for r in range (rr)))
    - mdl.sum(mdl.sum(mdl.sum(ax10[(j,t,rp,jp,tp)]*dh[rp][jp][tp] for rp in range (rr))for jp in range (jj))for tp in range (tt))>= 0 for j in range (jj) for t in range (1,tt))    
    
    mdl.add_constraints(ax10[(j,0,rp,jp,tp)]>=   mdl.sum(mdl.sum( alpha[(w,j,k,0,rp,jp,tp)] for w in range (ww)) for k in range(kk)) + mdl.sum(beta[(r,j,0,rp,jp,tz)] for r in range (rr)) for j in range (jj) for jp in range (jj)for tp in range (tt)for rp in range (rr))
    mdl.add_constraints(ax10[(j,0,rp,jp,tp)]>= - (mdl.sum(mdl.sum( alpha[(w,j,k,0,rp,jp,tp)] for w in range (ww)) for k in range(kk)) + mdl.sum(beta[(r,j,0,rp,jp,tz)] for r in range (rr))) for j in range (jj)  for jp in range (jj)for tp in range (tt)for rp in range (rr))
    
    mdl.add_constraints(ax10[(j,t,rp,jp,tp)]>=   mdl.sum(mdl.sum( (alpha[(w,j,k,t,rp,jp,tp)]-alpha[(w,j,k,t-1,rp,jp,tp)]) for w in range (ww)) for k in range(kk)) + mdl.sum(beta[(r,j,t,rp,jp,tz)] for r in range (rr)) for j in range (jj) for t in range (1,tt) for jp in range (jj)for tp in range (tt)for rp in range (rr))
    mdl.add_constraints(ax10[(j,t,rp,jp,tp)]>= - (mdl.sum(mdl.sum( (alpha[(w,j,k,t,rp,jp,tp)]-alpha[(w,j,k,t-1,rp,jp,tp)]) for w in range (ww)) for k in range(kk)) + mdl.sum(beta[(r,j,t,rp,jp,tz)] for r in range (rr))) for j in range (jj) for t in range (1,tt) for jp in range (jj)for tp in range (tt)for rp in range (rr))
    #non negtivity B
    mdl.add_constraints(B_bar[(r,j,t)]+mdl.sum(mdl.sum(mdl.sum((gamma[(r,j,t,rp,jp,tp)]*dz[rp][jp][tp]-ax11[(r,j,t,rp,jp,tp)]*dh[rp][jp][tp]) for rp in range (rr))for jp in range (jj))for tp in range (tt))>=0 for r in range (rr)for j in range  (jj)for t in range (tt))
    mdl.add_constraints(ax11[(r,j,t,rp,jp,tp)]>=gamma[(r,j,t,rp,jp,tp)]for j in range (jj)for r in range (rr)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt))
    mdl.add_constraints(ax11[(r,j,t,rp,jp,tp)]>=-gamma[(r,j,t,rp,jp,tp)]for j in range (jj)for r in range (rr)for t in range (tt)for rp in range (rr)for jp in range (jj)for tp in range (tt))
    
    # greater than t
    mdl.add_constraints(alpha[(w,j,k,t,rp,jp,tp)]==0  for w in range (ww)for j in range (jj) for k in range (kk) for t in range (tt) for rp in range (rr)for jp in range (jj)for tp in range (t+1,tt))
    mdl.add_constraints(beta[(r,j,t,rp,jp,tp)]==0   for r in range (rr)for j in range (jj) for t in range (tt) for rp in range (rr)for jp in range (jj)for tp in range (t+1,tt))
    
    
    # cosntraint for k_1
    mdl.add_constraints(mdl.sum(mdl.sum(I_bar[(w,j,k,t)]+mdl.sum(mdl.sum(mdl.sum(alpha[(w,j,k,t,r,jp,tp)]*dz[rp][jp][tp]for tp in range (tt))for jp in range (jj))for rp in range (rr))for w in range(ww))for t in range(tt))
    +mdl.sum(mdl.sum(mdl.sum(ax12[(j,k,rp,jp,tp)]*dh[rp][jp][tp]for tp in range (tt))for jp in range (jj))for rp in range (rr)) <=bigM *o_w[j][k] for j in range(jj)for k in range(kk))
    mdl.add_constraints(ax12[(j,k,rp,jp,tp)]>= mdl.sum(mdl.sum(mdl.sum(alpha[(w,j,k,t,rp,jp,tp)] for w in range (ww)) for t in range (tt))  for j in range (jj))for k in range (kk) for tp in range (tt) for jp in range (jj)for rp in range (rr))
    mdl.add_constraints(ax12[(j,k,rp,jp,tp)]>=-mdl.sum(mdl.sum(mdl.sum(alpha[(w,j,k,t,rp,jp,tp)] for w in range (ww)) for t in range (tt))  for j in range (jj))for k in range (kk) for tp in range (tt) for jp in range (jj)for rp in range (rr))

   #coefficient zero
# =============================================================================
#     mdl.add_constraints(beta[(r,j,t,rp,jp,tp)]==0 for rp in range (rr)for jp in range (jj)for tp in range (tt) for r in range (rr)for j in range (jj) for t in range (tt))
#     mdl.add_constraints(alpha[(w,j,k,t,rp,jp,tp)]==0 for rp in range (rr)for jp in range (jj)for tp in range (tt) for w in range (ww)for j in range (jj) for k in range (kk) for t in range (tt))
#     mdl.add_constraints(gamma[(r,j,t,rp,jp,tp)]==0 for rp in range (rr)for jp in range (jj)for tp in range (tt) for r in range (rr)for j in range (jj) for t in range (tt))
# =======================
    
    #mdl.export('ARCMilk')
    XARC=np.zeros((rr,jj,tt))
    IARC=np.zeros((ww,jj,kk,tt))
    BARC=np.zeros((rr,jj,tt))
    betaARC=np.zeros((rr,jj,tt,rr,jj,tt))
    alphaARC=np.zeros((ww,jj,kk,tt,rr,jj,tt))
    gammaARC=np.zeros((rr,jj,tt,rr,jj,tt))
    ARCmodeling = time.time()
    try:
      feasible = 1
      solution = mdl.solve(log_output= True)
      ARCobj = solution.get_objective_value()
      print(solution)
    except:
      feasible = 0
      ARCobj=-125
    if feasible == 1:
      ARCchoose = time.time()
      ARCgenrate=ARCmodeling-ARCstart
      ARCsolve=ARCchoose-ARCmodeling
      for w in range(ww):
          for k in range (kk):
              for t in range (tt):                            
                  for j in range (jj):
                      IARC[w][j][k][t]=float(I_bar[(w,j,k,t)].solution_value)
                      for rp in range (rr):
                          for jp in range (jj):
                              for tp in range (tt):
                                  alphaARC[w][j][k][t][rp][jp][tp]=float(alpha[(w,j,k,t,rp,jp,tp)].solution_value)
                                      
      for r in range(rr):
          for j in range (jj):
              for t in range (tt):
                  BARC[r][j][t]=float(B_bar[(r,j,t)].solution_value)
                  XARC[r][j][t]=float(X_bar[(r,j,t)].solution_value)
                  for rp in range (rr):
                          for jp in range (jj):
                              for tp in range (tt):
                                  gammaARC[r][j][t][rp][jp][tp]=float(gamma[(r,j,t,rp,jp,tp)].solution_value)
                                  betaARC[r][j][t][rp][jp][tp]=float(beta[(r,j,t,rp,jp,tp)].solution_value)
    #production    
    total_I  = sum(sum(sum(sum( IARC[w][j][k][t]+mdl.sum(mdl.sum(mdl.sum(alphaARC[w][j][k][t][rp][jp][tp]*D_W[rp][jp][tp]for rp in range (rr))for jp in range (jj))for tp in range (tt))for w in range (ww))for j in range (jj))for t in range (tt))for k in range (kk))
    holdcost = sum(sum(sum(sum(h[k][j]*(IARC[w][j][k][t]+sum(sum(sum(alphaARC[w][j][k][t][rp][jp][tp]*D_W[rp][jp][tp] for rp in range (rr)) for jp in range (jj)) for tp in range (tt)))for w in range (ww))for j in range (jj))for k in range (kk))for t in range (tt))
    
    #transportation
    total_X =   sum(sum(sum( XARC[r][j][t]+sum(sum(sum(betaARC[r][j][t][rp][jp][tp]*D_W[rp][jp][tp]for rp in range (rr))for jp in range (jj))for tp in range (tt))for j in range (jj))for t in range (tt))for r in range (rr))
    transcost = sum(sum(sum(transperunit[j][r]*(XARC[r][j][t]+sum(sum(sum(betaARC[r][j][t][rp][jp][tp]*(dz[rp][jp][tp]+dh[rp][jp][tp]) for rp in range (rr)) for jp in range (jj)) for tp in range (tt)))for r in range (rr))for j in range (jj))for t in range (tt))    
    #holding
    
    total_Y = sum(sum(sum(IARC[w][j][k][0]-ini_inv[w][j][k]+sum(sum(sum( alphaARC[w][j][k][0][rp][jp][tp]*D_W[rp][jp][tp]for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for w in range (ww))for k in range(kk))+
    sum(XARC[r][j][0]+sum(sum(sum((betaARC[r][j][0][rp][jp][tp])*D_W[rp][jp][tp]for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for r in range (rr))for j in range (jj)) + sum(sum(sum(sum(IARC[w][j][k][t]-IARC[w][j][k][t-1]+sum(sum(sum(((alphaARC[w][j][k][t][rp][jp][tp]-alphaARC[w][j][k][t-1][rp][jp][tp])*D_W[rp][jp][tp])for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for w in range (ww))for k in range(kk))+
    sum(XARC[r][j][t]+sum(sum(sum((betaARC[r][j][t][rp][jp][tp]*D_W[rp][jp][tp])for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for r in range (rr))for j in range (jj))for t in range (1,tt)) 

    prodcost = sum(p[j]*((sum(sum(IARC[w][j][k][0]-ini_inv[w][j][k]+sum(sum(sum( alphaARC[w][j][k][0][rp][jp][tp]*D_W[rp][jp][tp]for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for w in range (ww))for k in range(kk)))+
    sum(XARC[r][j][0]+sum(sum(sum((betaARC[r][j][0][rp][jp][tp])*D_W[rp][jp][tp]for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for r in range (rr)))for j in range (jj)) + sum(sum(p[j]*((sum(sum(IARC[w][j][k][t]-IARC[w][j][k][t-1]+sum(sum(sum(((alphaARC[w][j][k][t][rp][jp][tp]-alphaARC[w][j][k][t-1][rp][jp][tp])*D_W[rp][jp][tp])for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for w in range (ww))for k in range(kk)))+
    sum(XARC[r][j][t]+sum(sum(sum((betaARC[r][j][t][rp][jp][tp]*D_W[rp][jp][tp])for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for r in range (rr)))for j in range (jj))for t in range (1,tt)) 
    #backorde
    total_B = sum(sum(sum(BARC[r][j][t]+sum(sum(sum(gammaARC[r][j][t][rp][jp][tp]*(dz[rp][jp][tp]+dh[rp][jp][tp]) for rp in range (rr)) for jp in range (jj)) for tp in range (tt))for r in range (rr))for j in range (jj))for t in range (tt))  
    backcost = sum(sum(sum(c_b[j]*(BARC[r][j][t]+sum(sum(sum(gammaARC[r][j][t][rp][jp][tp]*(dz[rp][jp][tp]+dh[rp][jp][tp]) for rp in range (rr)) for jp in range (jj)) for tp in range (tt)))for r in range (rr))for j in range (jj))for t in range (tt))  
    real_demand = (sum(sum(sum(D_W[r][j][t]for r in range (rr))for j in range (jj)) for t in range (tt)))
    #number of decision variables 
    print(ARCobj,feasible)
    return(ARCobj,IARC,alphaARC,prodcost,XARC,betaARC,transcost,BARC,gammaARC,backcost,holdcost,ARCgenrate,ARCsolve)

#sol=ARC()

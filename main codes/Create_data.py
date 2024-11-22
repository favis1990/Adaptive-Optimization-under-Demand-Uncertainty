#this code creates input data for synthetic instances

import numpy as np
import time



def data_create (ww, kk, jj, rr, tt, zeta, seed):
    rnd = np.random    
    rnd.seed(seed)
    
    CAP = np.zeros((jj))
    #COEF is the coeficent of the cost for sensitivity analysis
    #Demand
    dz=  demand_matrix = [[[np.random.randint(5, 26) for t in range(tt)] for j in range(jj)] for r in range(rr)]
    D_W= np.zeros((rr,jj,tt))
    dh= np.zeros((rr,jj,tt))
    
    CAP = [sum(sum(demand_matrix[r][j][t] for r in range(rr)) for t in (0, 1)) for j in range(jj)]
    D_W = [[[demand_matrix[r][j][t] * (1 + zeta) for t in range(tt)] for j in range(jj)] for r in range(rr)]
    dh = [[[demand_matrix[r][j][t] * zeta for t in range(tt)] for j in range(jj)] for r in range(rr)]

    fixed_cap= 1
   
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
    
    return CAP, zeta, q_p, p, q, transperunit, c_b, h, dz, D_W, dh, ini_inv, o_w, sigma, C_O_W, bigM, total_demand, q_w
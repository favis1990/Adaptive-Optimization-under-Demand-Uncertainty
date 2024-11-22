#we can use this script to run different models 

from Create_data import data_create


ww = 2
kk = 3
jj = 5
rr = 5
tt = 5
zeta = 0.2
seed = 1

#this is for generating input data
CAP, zeta, q_p, p, q, transperunit, c_b, h, dz, D_W, dh, ini_inv, o_w, sigma, C_O_W, bigM, total_demand, q_w = data_create (ww, kk, jj, rr, tt, zeta, seed)


#This is used to get the deterministic solution
from NewRCfunction import RC
objective_value, Z_P_Values, Ydeter, Xdeter, Bdeter, Ideter, setup_cost_first, setup_cost, DETmodeling, DETtotal, backorder_det, prod_det, holding_det, transport_det, opening_cost, number_of_DV_RC_model = RC(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, D_W, o_w, q, sigma, bigM, C_O_W, ini_inv)


from ARPDP1 import ARC
#total setup cost of holding opening_cost  (openning periods) , setup_cost_first (set up first period), setup_cost (setup transission)
ARCobj,IARC,alphaARC,prodcost,XARC,betaARC,transcost,BARC,gammaARC,backcost,holdcost,ARCgenrate,ARCsolve = ARC(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, D_W, dz, dh, o_w,q, sigma, bigM, opening_cost, setup_cost_first, setup_cost, Z_P_Values, ini_inv)


from ARPDP2 import ARC2
#total setup cost of holding opening_cost  (openning periods) , setup_cost_first (set up first period), setup_cost (setup transission)
ARCobj,IARC,alphaARC,prodcost,XARC,betaARC,transcost,BARC,gammaARC,backcost,holdcost,ARCgenrate,ARCsolve = ARC2(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, D_W, dz, dh, o_w,q, sigma, bigM, C_O_W, ini_inv)


from NewREARCfunction import REARC
RARCobj, IRARC, alphaRARC, Rprodcost, XRARC, betaRARC, Rtranscost, BRARC, gammaRARC, Rbackcost, Rholdcost, RARCgenrate, RARCsolve = REARC(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, D_W, dz, dh, o_w, q, sigma, bigM, opening_cost, setup_cost_first, setup_cost, Z_P_Values, objective_value, ini_inv)


from NewREARCfunction2 import REARC2
RARCobj, IRARC, alphaRARC, Rprodcost, XRARC, betaRARC, Rtranscost, BRARC, gammaRARC, Rbackcost, Rholdcost, RARCgenrate, RARCsolve = REARC2(jj, kk, tt, rr, ww, q_p, q_w, c_b, h, p, transperunit, CAP, D_W, dz, dh, o_w, q, sigma, bigM, C_O_W, objective_value, ini_inv)
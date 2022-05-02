#ReaxFF Calculation for Hydrogen System
import numpy as np
from prettytable import PrettyTable
import math
import matplotlib.pyplot as plt #For graphing
myTable = PrettyTable(["Radius", "E-Bond","E-VDW","E-COl","E-TOT"])
param_BOE={"p_bo1":-0.016,"p_bo2":5.98,"p_boc3":5.02,"p_boc4":18.32,"p_boc5":8.32,"diss_Energy":168.4,
"p_be1":-0.310,"p_be2":10.25,"r_sig":0.656}
param_VDW={"vd_E1":1.69,"lmd_w_C":1.41,"lmd_w_H":5.36,"lmd_w_CH":3.385,
"alpha_CC":10.71,"alpha_CH":10.385,"alpha_HH":10.06,"d_HC":0.0862,"d_HH":0.0528,"d_HH":0.0194,
"r_vdw_CC":3.912,"r_vdw_CH":3.7805,"r_vdw_HH":3.649}
param_CIE={"q_C":-0.409128,"q_H":0.146488,"del_C":0.69,"del_H":0.37,"del_CH":0.58,"shielded_Constant":1.279}
e_rad_graph=[]
e_tot_graph=[]
e_bond_graph=[]
e_vdw_graph=[]
e_C_graph=[]

r_ij1=float(input(("Enter the starting radius range for calculation: ")))
r_ij2=float(input(("Enter the ending radius range for calculation: ")))
rad_incre=float(input(("Enter the radius gap for calculation: ")))
print(str(r_ij1)+"---"+str(r_ij2))
r_ij=r_ij1
while(r_ij<=r_ij2):
    
    #Bond energy calculation
    bo_u=(math.exp(param_BOE["p_bo1"]*((r_ij/param_BOE["r_sig"])**param_BOE["p_bo2"])))#Bond Order Uncorrected
    d_H=-1+bo_u#Delta C Value
    f_2=2*math.exp(-50*d_H)
    f_3=d_H
    f_1=(1+f_2)/(1+f_2+f_3)
    f_4=1/(1+math.exp(((-1*param_BOE["p_boc3"])*((param_BOE["p_boc4"]*bo_u*bo_u)-d_H))+param_BOE["p_boc5"]))
    f_5=f_4
    bo_c=bo_u*f_1*f_4*f_5
    e_bond=(math.exp(param_BOE["p_be1"]*(1-(bo_c**param_BOE["p_be2"])))*(-1)*param_BOE["diss_Energy"]*bo_c)
    
    #vanderwaal energy calculation
    f_13_HH=(((r_ij)**param_VDW["vd_E1"])+((1/param_VDW["lmd_w_H"])**param_VDW["vd_E1"]))**(1/param_VDW["vd_E1"]) #correction term
    e_vdw=param_VDW["d_HH"]*((math.exp(param_VDW["alpha_HH"]*(1-(f_13_HH/param_VDW["r_vdw_HH"]))))-(2*(math.exp(0.5*param_VDW["alpha_HH"]*(1-(f_13_HH/param_VDW["r_vdw_HH"]))))))
    

    #coloumb Interaction Energy
    e_C=param_CIE["shielded_Constant"]*param_CIE["q_H"]*param_CIE["q_H"]/((((r_ij)**3)+((1/param_CIE["del_H"])**3))**(1/3))
    

    e_tot=100.01175828774689+e_bond+e_vdw+e_C
    #storing the values for plotting
    #e_C_graph.append(e_C)
    e_bond_graph.append(e_bond)
    e_vdw_graph.append(e_vdw)
    e_rad_graph.append(r_ij)
    e_tot_graph.append(e_tot)
    

    myTable.add_row([r_ij, e_bond,e_vdw,e_C,e_tot])
    r_ij=r_ij+rad_incre
print(myTable)
# plotting the Energy points
plt.plot(e_rad_graph, e_tot_graph, label = "Total energy graph")
#plt.plot(e_rad_graph, e_bond_graph, label = "Bond energy graph")
#plt.plot(e_rad_graph, e_vdw_graph, label = "VanDerWaal energy graph")
#clplt.plot(e_C_graph, e_vdw_graph, label = "Coulomb Energy Line")



#File_Object=open("Output.txt","w")
#File_Object.write(myTable)

# naming the x axis
plt.xlabel('H-H radius in angstroms')
# naming the y axis
plt.ylabel('energy in kcal/mol')
# giving a title to my graph
plt.title('Energy calculations for H-H System')
 
# show a legend on the plot
plt.legend()
 
# function to show the plot
plt.show()

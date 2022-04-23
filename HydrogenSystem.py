import math
import matplotlib.pyplot as plt
#CALCULATION FOR H-H bond energy SYSTEM
p_bo1=-0.016
r_sig=0.7853
p_bo2=5.98
p_boc3=5.02
p_boc4=18.32
p_boc5=8.32
diss_Energy=168.4
p_be1=-0.310
p_be2=10.25
#calculation for h-h vanderwaal energy calculation
alpha_hh=10.06
dbond_HH=0.0194
r_vdw=3.649
lmd_w=5.36#lambda w values
p_vdw=1.69#vanderwaal energy
e_bond_graph=[]
e_vdw_graph=[]
e_rad_graph=[]
e_tot_graph=[]

r_ij1=float(input(("Enter the starting radius range for calculation")))
r_ij2=float(input(("Enter the ending radius range for calculation")))
rad_incre=float(input(("Enter the radius gap for calculation")))
print(str(r_ij1)+"---"+str(r_ij2))
r_ij=r_ij1
while(r_ij<=r_ij2):
    bo_u=math.exp(p_bo1*((r_ij/r_sig)**p_bo2))#Bond Order Uncorrected
    d_H=-1+bo_u#Delta H Value
    f_2=2*math.exp(-50*d_H)
    f_3=d_H
    f_1=(1+f_2)/(1+f_2+f_3)
    f_4=1/(1+math.exp(((-1*p_boc3)*((p_boc4*bo_u*bo_u)-d_H))+p_boc5))
    f_5=f_4
    bo_c=bo_u*f_1*f_4*f_5
    #Bond energy calculation
    e_bond=math.exp(p_be1*(1-(bo_c**p_be2)))*(-1)*diss_Energy*bo_c

    #vanderwaal energy calculation
    f_13=(((r_ij)**p_vdw)+((1/lmd_w)**p_vdw))**(1/p_vdw) #correction term
    e_vdw=dbond_HH*((math.exp(alpha_hh*(1-(f_13/r_vdw))))-(2*(math.exp(0.5*alpha_hh*(1-(f_13/r_vdw))))))
    e_tot=e_vdw+e_bond
    #storing the values for plotting
    e_bond_graph.append(e_bond)
    e_vdw_graph.append(e_vdw)
    e_rad_graph.append(r_ij)
    e_tot_graph.append(e_tot)
    

    print("Radius:"+str(r_ij)+"Bond Energy:"+str(e_bond)+"van der waals energy"+str(e_vdw))
    r_ij=r_ij+rad_incre

# plotting the line 1 points
plt.plot(e_rad_graph, e_bond_graph, label = " bond energy line")

# plotting the line 2 points
plt.plot(e_rad_graph, e_vdw_graph, label = "van der waal graph")
# plotting the line 3 points

plt.plot(e_rad_graph, e_tot_graph, label = "total energy graph")

# naming the x axis
plt.xlabel('radius in angstroms')
# naming the y axis
plt.ylabel('energy in kcal/mol')
# giving a title to my graph
plt.title('Energy calculations for H-H System')
 
# show a legend on the plot
plt.legend()
 
# function to show the plot
plt.show()




#ReaxFF Calculation for Ethylene System
#gen_param={"l1":50.0,"l2":15.61,"l3":5.02,"l4":18.32,"l6":8.32,"l7":1.94,"l8":-3.47,"l9":5.79,"l10":12.38,"l11":1.49,"l12":1.28,"l13":6.30,"l14":2.72,"l15":33.87,"l16":6.70,"l17":1.06,"l":,"l5":,"l5":,"l5":,"l5":,"l5":,"l5":,"l5":,"l5":,"l5":,"l5":,"l5":}
import cmath
import math
import matplotlib.pyplot as plt #For graphing
param_BOE={"p_bo1":-0.097,"p_bo2":38,"p_bo3":-0.26,"p_bo4":9.37,
"p_boc3":5.02,"p_boc4":18.32,"p_boc5":8.32,"diss_Energy":145.2,
"p_be1":0.318,"p_be2":0.65,"r_pi":1.266,"bo_CH":0.8703022784}
param_UCE={"p_under":29.4,"uc_E1":1.94,"uc_E2":-3.47,"uc_E3":5.79,"uc_E4":12.38}
param_VAE={"va_E1":1.49,"va_E2":1.28,"va_E3":6.30,"va_E4":2.72,
"va_E5":33.87,"va_E6":6.70,"va_E7":1.06,"va_E8":2.04,
"thet_HCH":69.94,"thet_HCC":71.56,"ka_HCC":29.65,
"kb_HCC":5.29,"ka_HCH":69.94,"kb_HCH":1.00}
param_TAE={"ta_E1":3.17,"ta_E2":10.00,"ta_E3":0.90,"v2":26.5,"v3":0.37,"pt":-2.33}
tors_angle=[0,-180,-180,-0]
param_CSE={"cs_E1":-1.14,"cs_E2":2.17}
param_VDW={"vd_E1":1.69,"lmd_w_C":1.41,"lmd_w_H":5.36,"lmd_w_CH":3.385,
"alpha_CC":10.71,"alpha_CH":10.385,"alpha_HH":10.06,"d_CC":0.0862,"d_CH":0.0528,"d_HH":0.0194,
"r_vdw_CC":3.912,"r_vdw_CH":3.7805,"r_vdw_HH":3.649}
param_CIE={"q_C":-0.409128,"q_H":0.146488,"del_C":0.69,"del_H":0.37,"del_CH":0.58,"shielded_constant":1.48}
#calculation for h-h vanderwaal energy calculation
#alpha_hh=10.06
#dbond_HH=0.0194
#r_vdw=3.649
#lmd_w=5.36#lambda w values
#p_vdw=1.69#vanderwaal energy
#e_bond_graph=[]
#e_vdw_graph=[]
e_rad_graph=[]
e_tot_graph=[]


r_ij1=float(input(("Enter the starting radius range for calculation: ")))
r_ij2=float(input(("Enter the ending radius range for calculation: ")))
rad_incre=float(input(("Enter the radius gap for calculation: ")))
print(str(r_ij1)+"---"+str(r_ij2))
r_ij=r_ij1
while(r_ij<=r_ij2):
    #Bond energy calculation
    bo_u=(math.exp(param_BOE["p_bo1"]*((r_ij/param_BOE["r_pi"])**param_BOE["p_bo2"])))+(math.exp(param_BOE["p_bo3"]*((r_ij/param_BOE["r_pi"])**param_BOE["p_bo4"])))#Bond Order Uncorrected
    d_C=-4+bo_u+(2*param_BOE["bo_CH"])#Delta C Value
    f_2=2*math.exp(-50*d_C)
    f_3=d_C
    f_1=(1+f_2)/(1+f_2+f_3)
    f_4=1/(1+math.exp(((-1*param_BOE["p_boc3"])*((param_BOE["p_boc4"]*bo_u*bo_u)-d_C))+param_BOE["p_boc5"]))
    f_5=f_4
    bo_c=bo_u*f_1*f_4*f_5
    e_bond=(math.exp(param_BOE["p_be1"]*(1-(bo_c**param_BOE["p_be2"])))*(-1)*param_BOE["diss_Energy"]*bo_c)
    
    #UnderCoordinationEnergy Calculation
    f_6=1/(1+(param_UCE["uc_E3"]*math.exp(param_UCE["uc_E4"]*d_C*bo_c)))
    e_under=-1*param_UCE["p_under"]*(1-math.exp(param_UCE["uc_E1"]*d_C))*f_6*(1/(1+math.exp(-1*param_UCE["uc_E3"]*d_C)))
 
    #Valence Angle Energy Calculation
    #HCC angles and HCH angles energy
    #f_7_HC=1-math.exp(-1*param_VAE["va_E1"]*(param_BOE["bo_CH"]**param_VAE["va_E2"]))
    #f_7_CC=1-math.exp(-1*param_VAE["va_E1"]*(bo_c**param_VAE["va_E2"]))
    #f_8=((2+math.exp(-1*param_VAE["va_E3"]*d_C))/(1+math.exp(-1*param_VAE["va_E3"]*d_C)))*(param_VAE["va_E4"]-((param_VAE["va_E4"]-1)*((2+math.exp(param_VAE["va_E5"]*d_C))/(1+math.exp(param_VAE["va_E5"]*d_C)))))
    #sbo_HCC=d_C-(2*(1-(cmath.exp(-5*((0.5*d_C)**param_VAE["va_E6"])))))+bo_c
    #sbo_HCH=d_C-(2*(1-cmath.exp(-5*((0.5*d_C)**param_VAE["va_E6"]))))
    #if float(sbo_HCC)>0:
    #    sbo2_HCC=sbo_HCC**param_VAE["va_E7"]
    #else:
     #   sbo2_HCC=0
    #if sbo_HCH>0:
     #   sbo2_HCH=sbo_HCH**param_VAE["va_E7"]
    #else:
    #    sbo2_HCH=0
    #thet0_HCC=180-(param_VAE["thet_HCC"]*(1-math.exp(-1*param_VAE["va_E8"]*(2-sbo2_HCC))))
    #thet0_HCH=180-(param_VAE["thet_HCH"]*(1-math.exp(-1*param_VAE["va_E8"]*(2-sbo2_HCH))))
    #E_val_HCC=f_7_CC*f_7_HC*f_8*(param_VAE["ka_HCC"]-(param_VAE["ka_HCC"]*math.exp(-1*param_VAE["kb_HCC"]*((thet0_HCC-param_VAE["thet_HCC"])**2))))
    #_val_HCH=f_7_HC*f_7_HC*f_8*(param_VAE["ka_HCH"]-(param_VAE["ka_HCH"]*math.exp(-1*param_VAE["kb_HCH"]*((thet0_HCH-param_VAE["thet_HCH"])**2))))
    #E_val=2*(E_val_HCC+E_val_HCH)
    
    #torsion angle energy
    f_10=(1-math.exp(-1*param_TAE["ta_E1"]*param_BOE["bo_CH"]))*(1-math.exp(-1*param_TAE["ta_E1"]*bo_c))*(1-math.exp(-1*param_TAE["ta_E1"]*param_BOE["bo_CH"]))
    f_11=(2+math.exp(-1*param_TAE["ta_E2"]*2*d_C))/(1+(2+math.exp(-1*param_TAE["ta_E2"]*2*d_C))+math.exp(param_TAE["ta_E3"]*2*d_C))
    E_tors=0
    for i in range(0,len(tors_angle)):
        E_tors+=f_10*(math.sin(math.radians(param_VAE["thet_HCC"])))*(math.sin(math.radians(param_VAE["thet_HCC"])))*((0.5*param_TAE["v2"]*math.exp(param_TAE["pt"]*((bo_c-3+f_11)**2)))*(1-math.cos(2*math.radians(tors_angle[i])))+(0.5*param_TAE["v3"]*(1+math.cos(math.radians(tors_angle[i])))))
    
    #Conjugated System Energy
    f_12=math.exp(-1*param_CSE["cs_E2"]*(((param_BOE["bo_CH"])-1.5)**2))*math.exp(-1*param_CSE["cs_E2"]*(((param_BOE["bo_CH"])-1.5)**2))*math.exp(-1*param_CSE["cs_E2"]*((bo_c-1.5)**2))
    E_cse=0
    for i in range(0,len(tors_angle)):
        E_cse+=f_12*param_CSE["cs_E1"]*(1+(((math.cos(math.radians(tors_angle[i])))**2)-1)*param_VAE["thet_HCC"]*param_VAE["thet_HCC"])
    
    #vanderwaal energy calculation
    f_13_CC=(((r_ij)**param_VDW["vd_E1"])+((1/param_VDW["lmd_w_C"])**param_VDW["vd_E1"]))**(1/param_VDW["vd_E1"]) #correction term
    f_13_CH=(((r_ij)**param_VDW["vd_E1"])+((1/param_VDW["lmd_w_CH"])**param_VDW["vd_E1"]))**(1/param_VDW["vd_E1"]) #correction term
    f_13_HH=(((r_ij)**param_VDW["vd_E1"])+((1/param_VDW["lmd_w_H"])**param_VDW["vd_E1"]))**(1/param_VDW["vd_E1"]) #correction term
    e_vdw_CC=param_VDW["d_CC"]*((math.exp(param_VDW["alpha_CC"]*(1-(f_13_CC/param_VDW["r_vdw_CC"]))))-(2*(math.exp(0.5*param_VDW["alpha_CC"]*(1-(f_13_CC/param_VDW["r_vdw_CC"]))))))
    e_vdw_CH=param_VDW["d_CH"]*((math.exp(param_VDW["alpha_CH"]*(1-(f_13_CH/param_VDW["r_vdw_CH"]))))-(2*(math.exp(0.5*param_VDW["alpha_CH"]*(1-(f_13_CH/param_VDW["r_vdw_CH"]))))))
    e_vdw_CC=param_VDW["d_HH"]*((math.exp(param_VDW["alpha_HH"]*(1-(f_13_HH/param_VDW["r_vdw_HH"]))))-(2*(math.exp(0.5*param_VDW["alpha_HH"]*(1-(f_13_HH/param_VDW["r_vdw_HH"]))))))
    e_vdw=e_vdw_CC+e_vdw_CH+e_vdw_CC

    #coloumb Interaction Energy
    e_C_CH=param_CIE["shielded_constant"]*param_CIE["q_C"]*param_CIE["q_H"]/((((r_ij)**3)+((1/param_CIE["del_CH"])**3))**(1/3))
    e_C_CC=param_CIE["shielded_constant"]*param_CIE["q_C"]*param_CIE["q_C"]/((((r_ij)**3)+((1/param_CIE["del_C"])**3))**(1/3))
    e_C=(4*e_C_CH)+e_C_CC

    e_tot=e_bond+e_under+E_tors+E_cse+e_vdw+e_C
    #storing the values for plotting
    #e_bond_graph.append(e_bond)
    #e_vdw_graph.append(e_vdw)
    e_rad_graph.append(r_ij)
    e_tot_graph.append(e_tot)
    

    print("Radius:  "+str(r_ij)+"Bond Energy:   "+str(e_tot))
    r_ij=r_ij+rad_incre

# plotting the line 1 points
#plt.plot(e_rad_graph, e_bond_graph, label = " C=C bond energy")

# plotting the line 2 points
#plt.plot(e_rad_graph, e_vdw_graph, label = "van der waal graph")
# plotting the line 3 points
plt.plot(e_rad_graph, e_tot_graph, label = "total energy graph")

# naming the x axis
plt.xlabel('c=c radius in angstroms')
# naming the y axis
plt.ylabel('energy in kcal/mol')
# giving a title to my graph
plt.title('Energy calculations for C=C System')
 
# show a legend on the plot
plt.legend()
 
# function to show the plot
plt.show()
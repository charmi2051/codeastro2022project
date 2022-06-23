# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
from astropy import constants as const
import matplotlib.pyplot as plt



class RotationalSpectra:


    def __init__(self, T, A, B, Jmax = 30, Kmax = 5):


        self.T = T 
        self.A = A 
        self.B = B
        self.C = A ## Just note C = A!

        self.Jmax = Jmax 
        self.Kmax = Kmax 

        self.h = const.h.cgs.value
        self.c = const.c.to('cm/s').value
        self.k = const.k_B.cgs.value

        self.origin = 15120

        # cm^{-1}
        self.ground_B = 0.0111
        self.excited_B = 0.010767

        self.round_C = 0.005552
        self.excited_C = 0.00538544



        #selection_rules =

        pass 



#Input constants:
    
T = 2
origin = 15120

ground_B = 0.0111
excited_B = 0.010767
delta_B = excited_B - ground_B

ground_C = 0.005552
excited_C = 0.00538544
delta_C = excited_C - ground_C

h = const.h.cgs.value
c = const.c.to('cm/s').value
k = const.k_B.cgs.value

#creating empty arrays: 
    
HL_Factor_pP_Branch = []
HL_Factor_rP_Branch = []
HL_Factor_pQ_Branch = []
HL_Factor_rQ_Branch = []
HL_Factor_pR_Branch = []
HL_Factor_rR_Branch = []
            
BD_pP_Branch = []
BD_rP_Branch = []
BD_pQ_Branch = []
BD_rQ_Branch = []
BD_pR_Branch = []
BD_rR_Branch = []

Intensity_pP_Branch = []
Intensity_rP_Branch = []
Intensity_pQ_Branch = []
Intensity_rQ_Branch = []
Intensity_pR_Branch = []
Intensity_rR_Branch = []

deviation_pP = []
deviation_rP = []
deviation_pQ = []
deviation_rQ = []
deviation_pR = []
deviation_rR = []

ground_energy_pP = []
ground_energy_rP = []
ground_energy_pQ = []
ground_energy_rQ = []
ground_energy_pR = []
ground_energy_rR = []

excited_energy_pP = []
excited_energy_rP = []
excited_energy_pQ = []
excited_energy_rQ = []
excited_energy_pR = []
excited_energy_rR = []


ground_energies = (ground_energy_pP, ground_energy_rP, ground_energy_pQ, ground_energy_rQ, ground_energy_pR, ground_energy_rR)
excited_energies = (excited_energy_pP, excited_energy_rP, excited_energy_pQ, excited_energy_rQ, excited_energy_pR, excited_energy_rR)
HL_factor = [HL_Factor_pP_Branch, HL_Factor_rP_Branch, HL_Factor_pQ_Branch, HL_Factor_rQ_Branch, HL_Factor_pR_Branch, HL_Factor_rR_Branch]
BD_factor = (BD_pP_Branch, BD_rP_Branch, BD_pQ_Branch, BD_rQ_Branch, BD_pR_Branch, BD_rR_Branch)
Intensities = (Intensity_pP_Branch, Intensity_rP_Branch, Intensity_pQ_Branch, Intensity_rQ_Branch, Intensity_pR_Branch, Intensity_rR_Branch)
deviation = (deviation_pP, deviation_rP, deviation_pQ, deviation_rQ, deviation_pR, deviation_rR)
label = ['pP', 'rP', 'pQ', 'rQ', 'pR', 'rR']
number_of_branches = [0,1,2,3,4,5]




#reading in PGopher file

#kmax = 1
#coronene = pd.read_csv(r"C:\Users\Charmi Bhatt\Desktop\Pgopher practice plots\Coronene P,Q and R Branches\coronene line lists pgopher\coronene_strengthtest1.txt", delim_whitespace= True)

#kmax = 5
#coronene = pd.read_csv(r"C:\Users\Charmi Bhatt\Desktop\Pgopher practice plots\Coronene P,Q and R Branches\coronene line lists pgopher\coronene_strengthtest2_k5.txt", delim_whitespace= True)

#normalized intenisty at kmax = 5 and origin = 15120
coronene = pd.read_csv(r"C:\Users\Charmi Bhatt\Desktop\Pgopher practice plots\Coronene P,Q and R Branches\coronene line lists pgopher\norminten-atk5-coronene.txt", delim_whitespace = True)


#honl london at kmax = 5 and origin = 15120
#coronene = pd.read_csv(r"C:\Users\Charmi Bhatt\Desktop\Pgopher practice plots\Coronene P,Q and R Branches\coronene line lists pgopher\hl-atk5-coronene.txt", delim_whitespace=(True))




'''>>>>>>>>>>>>>>>>>> Defining P, Q and R Branches <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'''

#sorting by value of J in ascending order
coronene = coronene.sort_values(by=['J'], ascending=True)

#coronene = coronene[(coronene['K'] != 0)] #different scale factors found for k = 0 and K > 0
#coronene = coronene[(coronene['Jprime'] != 0)]  #to avoid 'division by zero' error

P_Branch = coronene[(coronene['label'].str[1] == "P")]
Q_Branch = coronene[(coronene['label'].str[1] == "Q")]
R_Branch = coronene[(coronene['label'].str[1] == "R")]


pP_Branch = coronene[(coronene['label'].str[1] == "P") & (coronene['label'].str[0] == "p")]
rP_Branch = coronene[(coronene['label'].str[1] == "P") & (coronene['label'].str[0] == "r")]

pQ_Branch = coronene[(coronene['label'].str[1] == "Q") & (coronene['label'].str[0] == "p")]
rQ_Branch = coronene[(coronene['label'].str[1] == "Q") & (coronene['label'].str[0] == "r")]

pR_Branch = coronene[(coronene['label'].str[1] == "R") & (coronene['label'].str[0] == "p")]
rR_Branch = coronene[(coronene['label'].str[1] == "R") & (coronene['label'].str[0] == "r")]
rR_Branch = rR_Branch.iloc[1: , :]


Branches = (pP_Branch, rP_Branch, pQ_Branch, rQ_Branch, pR_Branch, rR_Branch)

lines = [] #number of transitions
index_of_lines = []
for n in number_of_branches:
    length = len(Branches[n])
    lines.append(length)
    
    
'''wavenumber'''

calculated_position_pP = []
calculated_position_rP = []
calculated_position_pQ = []
calculated_position_rQ = []
calculated_position_pR = []
calculated_position_rR = []

#Calculating Positions

for J, K in zip(pP_Branch['J'], pP_Branch['K']):
    wavenumber_pP = origin -J*(2*ground_B + delta_B) + (J**2)*delta_B + (K**2)*(delta_C - delta_B) + K*(2*excited_B - 2*excited_C) + (excited_C - excited_B)
    calculated_position_pP = np.append(calculated_position_pP, wavenumber_pP)

for J, K in zip(rP_Branch['J'], rP_Branch['K']):
    wavenumber_rP = origin - J*(2*ground_B + delta_B) + (J**2)*delta_B + (K**2)*(delta_C - delta_B) + K*(2*excited_C - 2*excited_B) + (excited_C - excited_B)
    calculated_position_rP = np.append(calculated_position_rP, wavenumber_rP)
    
for J, K in zip(pQ_Branch['J'], pQ_Branch['K']):
    wavenumber_pQ = origin + J*(J + 1)*delta_B + (K**2)*(delta_C - delta_B) + K*(2*excited_B - 2*excited_C) + (excited_C - excited_B)   
    calculated_position_pQ = np.append(calculated_position_pQ, wavenumber_pQ)

for J, K in zip(rQ_Branch['J'], rQ_Branch['K']):
    wavenumber_rQ = origin + J*(J + 1)*delta_B + (K**2)*(delta_C - delta_B) + K*(2*excited_C - 2*excited_B) + (excited_C - excited_B)
    calculated_position_rQ = np.append(calculated_position_rQ, wavenumber_rQ)
  
for J, K in zip(pR_Branch['J'], pR_Branch['K']):
    wavenumber_pR = origin + (J**2)*delta_B + J*(3*delta_B + 2*ground_B) + 2*excited_B + (K**2)*(delta_C - delta_B) + K*(2*excited_B - 2*excited_C) + (excited_C - excited_B)
    calculated_position_pR = np.append(calculated_position_pR, wavenumber_pR)

#for J in range(30):
for J, K in zip(rR_Branch['J'], rR_Branch['K']):           
     wavenumber_rR = origin + (J**2)*delta_B + J*(3*delta_B + 2*ground_B) + 2*excited_B + (K**2)*(delta_C - delta_B) + K*(2*excited_C - 2*excited_B) + (excited_C - excited_B)
     calculated_position_rR = np.append(calculated_position_rR, wavenumber_rR)

Calculated_positions = ( calculated_position_pP, calculated_position_rP, calculated_position_pQ, calculated_position_rQ, calculated_position_pR, calculated_position_rR)


    

'''Honl-London Factor'''

for n in number_of_branches:
        for J, K in zip(Branches[n]['J'], Branches[n]['K']): 
            #if K != 0:
                HL_equation_for_pP = ((J - 1 + K)*(J + K))/(2*J*(J+1))
                HL_equation_for_rP = ((J - 1 - K)*(J - K))/(2*J*(J+1))
                HL_equation_for_pQ = (J+1-K)*(J+K)/ (2*J*(J + 1))
                HL_equation_for_rQ = (J+1+K)*(J-K) / (2*J*(J + 1))
                HL_equation_for_pR = ( J + 2 - K)*( J + 1 - K)/ (2*(J + 1)*(2*J +1))
                HL_equation_for_rR = ( J + 2 + K)*( J + 1 + K)/ (2*(J + 1)*(2*J +1))
                
                HL_equations = (HL_equation_for_pP, HL_equation_for_rP, HL_equation_for_pQ, HL_equation_for_rQ, HL_equation_for_pR, HL_equation_for_rR)
    
                HL_factor[n] = np.append(HL_factor[n], HL_equations[n])
                
    
"HL factor using excited_J and excited_K"
# for n in number_of_branches:
#         for Jprime, Kprime in zip(Branches[n]['Jprime'], Branches[n]['Kprime']): 
            
#                 HL_equation_for_pP = ((Jprime - 1 + Kprime)*(Jprime + Kprime))/(2*Jprime*(Jprime+1))
#                 HL_equation_for_rP = ((Jprime - 1 - Kprime)*(Jprime - Kprime))/(2*Jprime*(Jprime+1))
#                 HL_equation_for_pQ = (Jprime+1-Kprime)*(Jprime+Kprime)/ (2*Jprime*(Jprime + 1))
#                 HL_equation_for_rQ = (Jprime+1+Kprime)*(Jprime-Kprime) / (2*Jprime*(Jprime + 1))
#                 HL_equation_for_pR = ( Jprime + 2 - Kprime)*( Jprime + 1 - Kprime)/ (2*(Jprime + 1)*(2*Jprime +1))
#                 HL_equation_for_rR = ( Jprime + 2 + Kprime)*( Jprime + 1 + Kprime)/ (2*(Jprime + 1)*(2*Jprime +1))
                
#                 HL_equations = (HL_equation_for_pP, HL_equation_for_rP, HL_equation_for_pQ, HL_equation_for_rQ, HL_equation_for_pR, HL_equation_for_rR)
    
#                 HL_factor[n] = np.append(HL_factor[n], HL_equations[n])
            
            

'''BD_factor'''

for n in number_of_branches:
    for J, K in zip(Branches[n]['J'], Branches[n]['K']):
                  E = ground_B*J*(J + 1) + (ground_C - ground_B)*(K**2)
                  ground_energies[n].append(E)
                  boltzmann_equation = ((2*J) + 1)*(np.exp((-h * c * E) / (k*T)))
                  BD_factor[n].append(boltzmann_equation)


       
#BD_factor using excited J and K 

# for n in number_of_branches:
#     for Jprime, Kprime in zip(Branches[n]['Jprime'], Branches[n]['Kprime']):
#                 Eprime = ground_B*Jprime*(Jprime + 1) + (ground_C - ground_B)*(Kprime**2)
#                 excited_energies[n].append(Eprime)
#                 boltzmann_equation = ((2*Jprime) + 1)*(np.exp((-h * c * Eprime ) / (k*T)))
#                 BD_factor[n].append(boltzmann_equation)
                
                
# for n in number_of_branches:
#     for Jprime, Kprime in zip(Branches[n]['Jprime'], Branches[n]['Kprime']):
#                 Eprime = ground_B*Jprime*(Jprime + 1) + (ground_C - ground_B)*(Kprime**2)
#                 excited_energies[n].append(Eprime)

'''Intensity'''

for n,l in zip(number_of_branches, lines):   
    #if n == 5:
        index = list(range(l))
        for i in index:
            Intensity_equaltion =  HL_factor[n][i] * BD_factor[n][i]
            Intensities[n].append(Intensity_equaltion) 

# print(ground_energies[0])
# print(excited_energies[0])
        



'''Deviation'''
        
for n,l in zip(number_of_branches, lines):   
    #if n ==0:
        index = list(range(l))
        for i in index:
            #ratio
            ratio = (Intensities[n][i] / Branches[n]['strength'].iloc[i])
            
            #difference
            # Normalized_Intenisty = (Intensities[n][i] /(max(Intensities[n]))) 
            # Normalized_Pgopher_strength = (Branches[n]['strength'][i]/(max(Branches[n]['strength'])))
            # difference = (Normalized_Intenisty - Normalized_Pgopher_strength)
            
            # deviation[n].append(difference)
            deviation[n].append(ratio)

            
# '''Plotting Deviation'''

for n in number_of_branches:
    #if n == 0:
        plt.scatter(Branches[n]['J'], deviation[n], label = label[n])
        plt.xlabel('J')
        plt.ylabel('deviation')
        plt.legend()
        #plt.ylim(0,4)
        # for i, l in enumerate(Branches[n]['label']):
        #         plt.annotate(Branches[n]['K'][i], (Branches[n]['J'][i], deviation[n][i]+ 0.05))
        #         plt.legend() 
            
        
        
print('max deviation in pP is' + str(max(deviation[0])))
print('min deviation in pP is' + str(min(deviation[0])))
     
# ''' Stem Plots'''

for n in number_of_branches:
    if n == 0:    
        fig, ax = plt.subplots(figsize=(25,6))
        ax.set(title = "Coronene  " + label[n] + "  Branch",
                xlabel = "Wavenumber",
                ylabel = "Normalized Inetensity")
            
        ax.stem(Branches[n]['position'], Intensities[n]/(max(Intensities[n])),linefmt=('red'), markerfmt= 'ro', label = 'calculated')
        ax.stem(Branches[n]['position'], Branches[n]['strength']/(max(Branches[n]['strength'])),linefmt=('blue'), markerfmt= 'bo', label = 'Pgopher')
        #ax.annotate('Both BD factor and HL fcator are calculated at ground levels', xy = (-0.8+1.500e3,0.8), xytext = (-0.8+1.500e3,0.8), fontsize=15)
        plt.legend()
        plt.show()
        

    
# '''Linelist'''
# #all_branches_table = {}
# for n in number_of_branches:
#     #if n == 5:
#             table_contents = { 'Branch' : label[n],
#                               'J' : Branches[n]['J'],                        
#                               'K' : Branches[n]['K'],
#                               'ground_E' : ground_energies[n],
#                               'Jprime' : Branches[n]['Jprime'],
#                               'Kprime' : Branches[n]['Kprime'],
#                               'excited_E' : excited_energies[n],
#                               'wavenumber' : Calculated_positions[n],
#                               'HL_factor' : HL_factor[n],
#                               'BD_factor' : BD_factor[n],
#                               'Intenisty' : Intensities[n],
#                               "PGopher strength" : Branches[n]['strength']}
   
           
                                
#             #linelist = pd.DataFrame(data = table_contents )
            
#             # print('--------------------------------')
#             # print(label[n])
#             # print('--------------------------------')
#             #print(linelist.to_string())
#             # print('--------------------------------')
            
#             #exporting line list
#             # linelist = linelist.round(decimals=6) 
#             # l#inelist.style.format(precision=6)  
#             # #linelist.style.set_caption("hello").format('{:.2f}')
# # linelist.to_excel(r"C:\Users\Charmi Bhatt\Desktop\rR.xlsx",  index=None) #sheet_name = 'rP') #sep='\t', mode='a')
# # print(len(coronene))


'''wavenumber'''

waveno_pP = []
waveno_rP = []
waveno_pQ = []
waveno_rQ = []
waveno_pR = []
waveno_rR = []

wavenos = (waveno_pP, waveno_rP, waveno_pQ, waveno_rQ, waveno_pR, waveno_rR )


for n in number_of_branches:
    for Jprime, Kprime in zip(Branches[n]['Jprime'], Branches[n]['Kprime']):
                Eprime = excited_B*Jprime*(Jprime + 1) + (excited_C - excited_B)*(Kprime**2)
                excited_energies[n].append(Eprime)
                
for n in number_of_branches:
   for J, K in zip(Branches[n]['J'], Branches[n]['K']):
                  E = ground_B*J*(J + 1) + (ground_C - ground_B)*(K**2)
                  ground_energies[n].append(E)
                  
for n,l in zip(number_of_branches, lines):   
    #if n == 5:
        index = list(range(l))
        for i in index:
            waveno = origin + excited_energies[n][i] - ground_energies[n][i]
            wavenos[n].append(waveno)
            
# print(wavenos[0])
# print(Calculated_positions[0])
# print(Branches[0]['position'])

# for n in number_of_branches:
#     if n == 0:    
#         fig, ax = plt.subplots(figsize=(25,6))
#         ax.set(title = "Coronene  " + label[n] + "  Branch",
#                 xlabel = "Wavenumber",
#                 ylabel = "Normalized Inetensity")
            
#         ax.stem(wavenos[n], Branches[n]['strength']/(max(Branches[n]['strength'])),linefmt=('red'), markerfmt= 'ro', label = 'calculated')
#         ax.stem(Calculated_positions[n], Branches[n]['strength']/(max(Branches[n]['strength'])),linefmt=('blue'), markerfmt= 'bo', label = 'Pgopher')
#         plt.legend()
#         plt.show()
        
        
# for n,l in zip(number_of_branches, lines):    
#         index = list(range(l))
#         for i in index:
#             difference = ((Calculated_positions[n][i]) - wavenos[n][i])
#             deviation[n].append(difference)
# #%%



# #print(deviation[0])
          
# #Plotting Offset
# for n in number_of_branches:
#     #if n == 2:
#           plt.plot(Branches[n]['J'], deviation[n], label = label[n]) 
#           plt.title('k_max = 20')
#           plt.legend()
#           plt.xlabel('J')
#           plt.ylabel('Offset')
          
# for n in number_of_branches:
#     if n == 5:
#             table_contents = { 'Branch' : label[n],
#                               'J' : Branches[n]['J'],                        
#                               'K' : Branches[n]['K'],
#                               'ground_E' : ground_energies[n],
#                               'Jprime' : Branches[n]['Jprime'],
#                               'Kprime' : Branches[n]['Kprime'],
#                               'excited_E' : excited_energies[n],
#                               'wavenumber' : Calculated_positions[n],
#                               'waveno' : wavenos[n]}
   
           
                                
#             linelist = pd.DataFrame(data = table_contents )
            
# linelist.to_excel(r"C:\Users\Charmi Bhatt\Desktop\rR.xlsx",  index=None) #sheet_name = 'rP') #sep='\t', mode='a')
       
          

                                                   
                       
                           
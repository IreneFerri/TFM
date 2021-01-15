
# It calculates the states energy density associated to the alpha3states
# opinion model and prints it in the screen, and writes it to a file.
#
import numpy as np
import math as m
import time
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Energy function for the fully connected graph
#
def energy(num_1, num0, num1, alpha2):
 return -0.5*((num1 - num_1)**2 + alpha2*num0**2 - num1 - num_1 - alpha2*num0)
# ----------------------------------------------------------------------------
start=time.time()       # Start measuring time
N = 100 # Number of agents
alpha = 0.2 # Neutrality parameter
#
filename = 'dos_N'+str(N)+'_a'+str(alpha)+'.csv'
alpha2 = alpha*alpha
total_states = 3**N
# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------
energy_list = []  # Stores the energy corresponding to every possible combination of agents
deg_list = []     # Stores the number of states associated to each combination
#
#
for num0 in range (N+1):
  for num1 in range (N+1-num0):
    num_1 = N - num1 - num0
    ene = energy(num_1, num0, num1, alpha2) 
    deg = m.factorial(N)//(m.factorial(num1)*m.factorial(num0)*m.factorial(num_1))  # Combinations
    energy_list.append(ene)
    deg_list.append(deg)
#  ----------------------------------------------------------------------------
#  Eliminate degeneration, sum states corresponding to the same energy -----------------------------
#
unique_ene = []  # Stores every different energy value
unique_deg = []  # Stores the number of states associated to every different energy value
min_E_val = min(energy_list) 
max_E_val = max(energy_list)
while (min_E_val<max_E_val):
  deg = 0
  remove_index = []
  for i in range(len(energy_list)):
    if (energy_list[i] == min_E_val):
      deg = deg + deg_list[i]  
      remove_index.append(i)
  unique_deg.append(deg)
  unique_ene.append(min_E_val)         
#
# Removes all occurrencies of the current min_E_val from energy_list
  remove_index = np.array(remove_index)
  energy_list = np.array(energy_list)
  deg_list = np.array(deg_list)
  energy_list = np.delete(energy_list, remove_index)
  deg_list = np.delete(deg_list, remove_index)
  energy_list = energy_list.tolist() 
  min_E_val = min(energy_list)  # New minimum value
unique_ene.append(energy_list[0])  # Add the maximum value 
deg = 0
for i in range(len(deg_list)): # Add the maximum value degeneracy
  deg = deg + deg_list[i] 
unique_deg.append(deg)
#
#  Compute the normalized dos and the logarithm of dos ------------------------------------------------
#
final_deg = [x / total_states for x in unique_deg]  # Normalize the states density
ln_deg = [m.log(x) for x in unique_deg]    # Takes the logarithm od the states density
#
#  Output --------------------------------------------------------------------------------------------
#
file=open(filename, 'w')
for i in range (len(unique_ene)):
  print(unique_ene[i], ln_deg[i], unique_deg[i])
  file.write(str(unique_ene[i]) + ' ' + str(ln_deg[i]) +'\n')
file.close()
#
end=time.time()                       # End measuring time
time_elapsed=end-start
print ('time elapsed',time_elapsed,'seconds')



# SPEA2-Paralell
Implementation of parallel Multiobjective Algorithm SPEA2 with c++ and MPI.
This algorithm uses different approaches to modify the population:

#Cooperative
  -> Execute GA, PSO and DE and chooses good individuals of all three.
  
#Competitive
  -> Get the population of best hipervolume
  
#Fuzzy
  -> Use the difference between real hipervolum and the found by algorithm.

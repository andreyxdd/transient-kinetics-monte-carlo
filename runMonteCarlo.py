import numpy as np

def runMonteCarlo(speciesVector):
  
  while True:
    rnd1=np.random.random()
    if rnd1 <= P1:
        prot1 += 1 #protonation 1
        while True:
            rnd2=np.random.random()
            if rnd2 <= P2: # protonation 1 over reverse
                prot1 -= 1
                int1 += 1 #protonated specie leading enol 1-2 or enol 2-3
                while True:
                    rnd3=np.random.random()
                    if rnd3 <= P3:
                        int1 -= 1
                        enol12 += 1
                    else: # enolation 2-3 (deprotonation step for enolation)
                        int1 -= 1
                        enol23 += 1
                        break
                    break
                break
            else:
                break
            break
        break
    else:
        prot2 += 1 #protonation 2
        while True:
            rnd2=np.random.random()
            if rnd2 <= P4: # protonation 2 over reverse
                prot2 -= 1
                int2 += 1 #protonated specie leading enol 1-2 or enol 2-3
                while True:
                    rnd3=np.random.random()
                    if rnd3 <= P5:
                        int2 -= 1
                        enol45 += 1
                    else: # enolation 2-3 (deprotonation step for enolation)
                        int2 -= 1
                        enol56 += 1
                        break
                    break
                break
            else:
                break
            break
        break
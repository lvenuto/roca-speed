#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
from PyRoca import PyRoca
from ParameterFinder import ParameterFinder
import logging
import time

if __name__ == "__main__":
    
    #Set the debug level you want:
    #DEBUG is very verbose 
    #INFO prints only the most important information
    #WARNING disables logs on the console
    logging.basicConfig(format='%(levelname)s: PyRoca %(message)s',level=logging.DEBUG)
    
    #Let's assume we want to factor this RSA-512 Key (chosen to complete in very little time, 
    #for random keys, expect an execution time of around 20-30 minutes with a relatively modern computer)
    modulus=0x958f12d6d12af295d46c6094ba292f020e77befdef0ac12b0c9b8c624fcba59fa36bc3e95a850eb3824d8ef4dac4a8f74b6d824151262d43913b0b9313ee0945
 
    #Before we start factoring it, we need the important parameters, M prime, mm and tt. You could run it with the full_m, but it would
    #take alot more to factor, so why not optimize the parameters as much as possible? This is where the ParameterFinder class comes in
    
    PF = ParameterFinder(n=modulus, general_prime=False) #we are assuming this is a RoCa key, so the generator used to generate that key is 65537
    
    start=time.time()
    print PF.is_roca() # We can also check if the n we passed is vulnerable to the roca attack
    end = time.time()
    print end-start #in very little time! 0.01 seconds for rsa-512

    #Now we can calculate the m_prime for that keysize (it's not gonna be tied to that specific modulus, these parameters will work for ALL rsa-512)
    #Thats why it's also possible to calculate the parameters (m_prime,mm,tt) specifying only the keysize ex. PF = ParameterFinder(512)
    
    #Check under the function definition for additional information
    time_and_m, parameters = PF.calculate_m_prime(steps=2, times_to_test=5, max_mm = 20, advanced_brute=False) 
    
    #Now time_and_m is a list ordered by time (the lower the better), and parameters is a dictionary containing the parameters for any M contained in time_and_m
    
    best_m = (time_and_m[0])[1] #Picking the "best" m prime
    mm,tt = parameters[best_m]  
    
    #We can now factorize the key
    RC = PyRoca(modulus, m=best_m, mm=mm, tt=tt, general_prime=False) #use the parameters we found earlier
    #by default it will use all cores avaliable in the processor
    solution = RC.factorize() #factor the key
    print 'The factored key is {}'.format(solution)
    

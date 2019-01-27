from functools import reduce
from operator import mul
from math import log
from flint import *
from PyRoca import PyRoca
from gmpy2 import is_prime
from random import randint
import time
import logging
import operator

class ParameterFinder:
        def __init__(self, keylen=None, n=None, generator=None, general_prime = True):
            if keylen is None:
                if n is not None:
                    #keylen = (int(math.log(n, 256)) + 1)*8
                    self.keylen = n.bit_length()
                    self.n = n
                else:
                    raise ValueError('Please specify keylen or n')
            else:
                self.keylen = keylen
                self.n = None
                
            if generator is None:
                logging.debug('Generator not specified, using the default of 65537')
                self.generator = 65537 #default generator for RoCa vulnerable keys
            
            if general_prime: #determines the size of XX
                self.general_prime = True 
            else: self.general_prime = False
            
            # Attributes necessary during the calculation of the best possible m_prime,
            # They are created by the greedy heuristic procedure, and later used in the
            # calculate_m_prime method reward_cost are numbers that are not picked by the
            # greedy heuristic (maybe because their reward_cost was low), while goback are
            # the prime numbers which are picked by the greedy heuristic.
            self.reward_cost = None
            self.goback = None 
        
        def is_roca(self):
            #check if the specified n is roca vulnerable
            if not self.n:
                raise ValueError('Please specify an n..') 
            full_m = self.calculate_full_m(self.keylen)
            order = self.multOrder(self.generator, full_m, fmpz(full_m).factor())
            decomp_order = fmpz(order).factor() 
            d = PyRoca.pohlig_hellman(self.n, self.generator, order, decomp_order, full_m)
            if d == None:
                return False
            return True
        
        @property
        def n(self):
            return self.n
             
        @n.setter
        def keylen(self, n):
            self.n = n
            self.keylen = n.bit_length()
            
        @property
        def keylen(self):
            return self.keylen
             
        @keylen.setter
        def keylen(self, keylen):
            self.keylen = keylen
                
        def _rwh_primes(self, n_primes):
            # https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
            # Returns  a list of primes < n_primes 
            n = n_primes
            sieve = [True] * n
            for i in xrange(3,int(n**0.5)+1,2):
                if sieve[i]:
                    sieve[i*i::2*i]=[False]*((n-i*i-1)/(2*i)+1)
            return [2] + [i for i in xrange(3,n,2) if sieve[i]]

        def primorial(self, n):
            # The product of the first n primes
            if not hasattr(self, 'primelist'):
                self.primelist = list(self._rwh_primes(2000))
            return reduce(operator.mul, self.primelist[:n], 1)
        
        def _gcd(self, a, b):
            while b != 0:
                a, b = b, a % b
            return a
 
        def _lcm(self, a, b):
            return (a*b) / self._gcd(a, b)
   
        def multOrder(self, a, m, mfactored):
            # Faster multiplicative order modified from https://rosettacode.org/wiki/Multiplicative_order#Python
            # Initially it was using a very slow way of finding the factorization of a number (it was trying all the prime numbers!)
            # ,which was working well until you start working with big numbers... luckily we have the FLINT library at our disposal to
            # factorize a number. a is the generator, like 65537 for RoCa vulnerable keys
   
            assert self._gcd(a, m) == 1
            mofs = []
            for r in mfactored:
                p=int(r[0])
                e=int(r[1])
                m = p**e
                t = (p-1)*(p**(e-1)) #  = Phi(p**e) where p prime
                qs = [1,]
                factored = fmpz(t).factor()
                for f in factored:
                    qs = [ q * int(f[0])**j for j in range(1+int(f[1])) for q in qs ]
                qs.sort()
                for q in qs:
                    if pow( a, q, m )==1:
                        mofs.append(q)
                        break
            return reduce(self._lcm, mofs, 1)
            
        def algorithm2(self, full_m, generator, divisors_ord_m):   
            # Returns a dictionary with the divisors of ord_m as the key (ex. (83,1) (53,1)) and as the value
            # a list with the divisors of full_m to remove for that given key (ex[167] for (83,1))    
            mfactored = fmpz(full_m).factor() #fix this, do the factorization once in the whole class..
            divisors_m = []
            for prime in mfactored:
                divisors_m.append(int(prime[0]))
            divisors_m.sort(reverse=True) #not needed right?
            removal_dict = {}
            for div_ord_m in divisors_ord_m:
                    if div_ord_m[1] > 1: #if the power is greater than 1
                         for power in xrange(1,div_ord_m[1]+1):
                            #Initialize dict values
                            removal_dict[(div_ord_m[0], power)] = [] #ex removal_dict[3 ** 4]
                            for prime in divisors_m:
                                ord_pi = self.multOrder(generator, prime, fmpz(prime).factor()) #es 166, p.s. check if you really need the factor
                                if (ord_pi % div_ord_m[0] ** power) == 0:# and (ord_m % div_ord_m) != 0: #es 166 % 83 == 0?
                                    removal_dict[(div_ord_m[0], power)].append(prime)   
                    else: #normal divisor of ord_m with power 1... like 83!     
                        #initialize dict values
                        removal_dict[div_ord_m] = [] # here it's like removal_dict[(83,1)]       
                        for prime in divisors_m:
                            ord_pi = self.multOrder(generator, prime, fmpz(prime).factor()) #es 166, p.s. check if you really need the factor
                            if (ord_pi % div_ord_m[0]) == 0:# and (ord_m % div_ord_m) != 0: #es 166 % 83 == 0?
                                removal_dict[div_ord_m].append(prime)    
            return removal_dict

        def rewardcost_calculator(self, removal_dict, divisors_ord_m):
            # Returns a decreasigly ordered list of tuples with (reward_at_cost, div_ord_m) 
            # for every div_ord_m
            reward_cost = []   
            max_power_of_order = {} #dict to avoid searching alot in the list divisors_ord_m
            for key in removal_dict:
                if key[1] > 1: # exponent is greater than 1, reward cost is different..
                    if key[0] not in max_power_of_order: #this div_ord_m is not present in the precomputed dict
                        #we do a full search for it
                        tuple_order_power = next((x for x in divisors_ord_m if x[0] == key[0]), None) #return for example (5,2) or (3,4)
                        max_power_of_order[tuple_order_power[0]] = tuple_order_power[1] #save it in the dict for possible future use
                    current_order = max_power_of_order[key[0]] # now we can use the dict to extract the max power of the order
                    numerator = log(key[0] ** (current_order - (key[1] - 1)), 2) #for the reward_at_cost calculation the numerator is current order - (the order we remove - 1)
                else:
                    numerator = log(key[0],2) #The numerator is easy, since the power is 1
                denominator = 0 
                for prime in removal_dict[key]:
                    denominator += log(prime,2) 
                if(denominator == 0): 
                    result = 0
                else:
                    result = numerator/denominator
                reward_cost.append((result,key))
                
            reward_cost.sort(reverse=True) #order so the best one is always in the head of the list
            return reward_cost

        #ex. [(83, 1), (53, 1), (41, 1), (37, 1), (29, 1), (23, 1), (17, 1), (13, 1), (11, 1), (7, 1), (5, 2), (3, 4), (2, 4)]
        def get_div_ord_m(self, full_m, generator, mfactored):
            # Returns a list with all the divisors of the order of full_m
            # it's ordereded in decreasing order, but its not really needed
   
            ord_m = self.multOrder(generator, full_m, mfactored)
            divisors_ord_m = []
            for elem in fmpz(ord_m).factor():
                divisors_ord_m.append((int(elem[0]), int(elem[1])))
            divisors_ord_m.reverse() #is it really needed?
            return divisors_ord_m  
                  
        def greedy_heuristic_m_finder(self, generator, full_m):
            # With the greedy heuristic described in the paper, we calculate a solution based on the reward cost of 
            # the ord prime 
            
            primes_removed = [] #sometimes different ord primes, remove the same prime(divisor of m), this is to avoid that a prime gets removed twice
            goback = {} #we need this later in another function to undo some operations
            mfactored = fmpz(full_m).factor()
            max_bits = log((2**self.keylen), 2) / 4 #most conservative limit, you could use log(n,2) / 4 here, but it shouldn't matter
            logging.debug("Greedy heuristic limit set at {} bits".format(max_bits))
            result = []
            candidate = None
            divisors_ord_m = self.get_div_ord_m(full_m, generator, mfactored)
            removal_dict = self.algorithm2(full_m, generator, divisors_ord_m)
            self.removal_dict = removal_dict
            reward_cost = self.rewardcost_calculator(removal_dict, divisors_ord_m)
            elements = len(removal_dict) #we need the real number of elements, in divisors_ord_m there are less elements(powers are grouped)
                      
            for i in range(0,elements):
                if candidate == reward_cost[0][1]:
                    logging.debug("We are in a stall situation, ending the search")
                    break
                candidate = reward_cost[0][1] # the list is always ordered by the function, we always pick the best candidate, ex. (83,1) or (3,2)
                logging.debug("Candidate is: {}".format(candidate))
                logging.debug("Its reward cost is: {}".format(reward_cost[0][0]))
                temp = full_m
                skip = False
                temp_primes_removed = [] 
                for prime in removal_dict[candidate]: #all the primes that we have to be able to remove
                    if prime not in primes_removed:
                        if log(full_m/prime,2) < max_bits: #we are going under the required n_bits, cannot pick this one
                            skip = True
                            break
                        else:
                            logging.debug("Dividing M by this prime: {}".format(prime))
                            full_m = full_m/prime
                            temp_primes_removed.append(prime)
                if skip: #if we skipped we have to restore the previous m
                    full_m = temp #restoring full_m   
                    logging.debug("Restoring for: {}".format(candidate))
                else:
                    goback[candidate] = temp_primes_removed
                    primes_removed.extend(temp_primes_removed)
                    result.append(candidate)
                    if candidate[1] > 1: #if the prime we removed has a power>1, we have to remove the chosen power, and insert back the powers < than power removed
                        logging.debug("Modifying: {}".format(divisors_ord_m))
                        tuple_order_power = next((x for x in divisors_ord_m if x[0] == candidate[0]), None)
                        divisors_ord_m.remove(tuple_order_power) #remove the power we chosen
                        #Now we must re-introduce in the list the div ord m with power < power we removed
                        divisors_ord_m.append((candidate[0],candidate[1] - 1))
                        divisors_ord_m.sort(reverse=True)
                        logging.debug("Modified: {}".format(divisors_ord_m))
                    else:
                        divisors_ord_m.remove(candidate)
                #Recalculate stuff, we have to recalculate the reward at cost every time!
                mfactored = fmpz(full_m).factor()
                removal_dict = self.algorithm2(full_m, generator, divisors_ord_m)
                reward_cost = self.rewardcost_calculator(removal_dict, divisors_ord_m)
                logging.debug("Current M: {}".format(full_m))
                logging.debug("bitsisze of M: {}".format(log(full_m,2)))    
            #Checking the result if two prime powers are removed (ex. (7,2) and (7,3)), we pick (7,2) but the result is right anyway.
            div_extra = {} #dict with 7 as the key, as the values a list with all the (prime,power) contained in result
            for div_ord in result:
                if div_ord[1] > 1: #only if power > 1 
                    if div_ord[0] not in div_extra: 
                        div_extra[div_ord[0]] = [] #remember to create the list 
                    div_extra[div_ord[0]].append(div_ord)
            for key in div_extra: 
                if len(div_extra[key]) > 1:
                    div_extra[key].sort()
                    keep = div_extra[key][0] #it makes sense to keep only the smallest prime power, it will remove more primes than greater powers 
                    for div_ord_power in div_extra[key]:
                        if div_ord_power != keep:
                            result.remove(div_ord_power)
                            logging.debug("Final_ Removing from the result {}, duplicate".format(div_ord_power))
            self.reward_cost = reward_cost
            self.goback = goback
            return result, full_m #note full_m is now greedym, TODO: fix variable names
   
        def calculate_XX(self, m):
            length = self.keylen
            if self.general_prime:
                XX = int(2**(length / 2)) // m 
            else:
                XX = int(2**(length / 2 - 1) + 2**(length / 2 - 2) + 2**(length / 2 - 4)) // m 
            return XX
        
        def calculate_k0guess(self, m, p, q):
            # From an m and the primes p and q we can calculate the k0guesses 
            # that are required to find those primes. Remember that one of them will
            # end up outside of the search interval used by PyRoca, we don't care
            # which one it is here, because we just want to find the correct mm and tt
            # and for that purpose either will work
            coppersmith_root1 = p/m
            coppersmith_root2 = q/m 
            k0_guess1 = p-coppersmith_root1*m
            k0_guess2 = q-coppersmith_root2*m
            return k0_guess1, k0_guess2
    
        def calculate_guess(self, k0_guess, order, decomp_order, m):
            # Get the guess used to calculate that k0guess, could be useful 
            # for debugging purposes, if k0_guess = n, we calculate the 
            # starting guess used by PyRoca
            #order = self.multOrder(self.generator, m, fmpz(m).factor())
            #decomp_order = fmpz(order).factor()   
            return PyRoca.pohlig_hellman(k0_guess, self.generator, order, decomp_order, m)
        
        def test_mm_tt(self, m, mm, tt, times=1):
            # You know mm and tt, thanks to brute_mm_tt, but you are not sure
            # that they work 100% of the times. Test them here
            # Watch out that exec_time refers to the time it takes to make a correct guess,
            # you should't use that time to formulate hypotheses on the time it takes to crack
            # a key, because the time for a correct guess >> time for a wrong guess
            time_list = []
            for i in xrange(0,times):
                n,p,q = self.gen_vuln_key()
                k0_guess1, k0_guess2 = self.calculate_k0guess(m,p,q)
                logging.debug('K0_guesses are: {} and {}'.format(k0_guess1,k0_guess2))
                m_inv = PyRoca.EuclidExt(m,n)
                expr = fmpz_poly([(k0_guess1*m_inv)%n,1]) #use any of the two guesses, it's ok either way
                XX = self.calculate_XX(m)  
                start=time.time()   
                result = PyRoca.coppersmith_howgrave_univariate(expr, n, mm, tt, XX)
                end=time.time()
                time_list.append(end-start) #to calculate an average of time a correct attempt takes..
                if (result != []):
                     #if the parameters are good we should find p
                    for root in result:
                        test_p = k0_guess1 + abs(root)*m
                        test_q = k0_guess2 + abs(root)*m
                        if (not p==test_p) and (not p==test_q) and not q==test_p and not q==test_q:
                            logging.info('Testing of the parameters failed, prime p is different')
                            return False, None
                else:
                     logging.info('Testing of the parameters failed')
                     return False, None
            logging.info('Parameters are working correctly {} times'.format(times))
            meantime = sum(time_list) / len(time_list)
            return True, meantime
              
        def brute_mm_tt(self, n, m, p, q, times=0, max_mm=40):
            # brute forcing parameters mm,tt
            # As I said before It doesn't matter if you use p or q here, don't bother to
            # find if a guess is !before the initial guess (AKA you will never find it
            # during a normal run of the algo)

            tt=0     
            #Check only mm,tt maximum 8,9 
            for mm in xrange(1,max_mm):
                tt=mm+1
                success, meantime = self.test_mm_tt(m, mm, tt, 1) 
                if success:
                    if times > 0:
                        logging.debug('These parameters mm = {}, tt = {}, work in this execution, lets check them {} more times!'.format(mm,tt,times))
                        success_again, meantime = self.test_mm_tt(m, mm, tt, times) 
                        if success_again:
                            logging.debug('Ok these mm={} tt={} really seem to work!'.format(mm,tt))
                            return mm,tt,meantime
                        else:
                            logging.debug('See? I told you to check again, this is a false positive mm={} tt={}'.format(mm,tt))
                    else:
                        logging.debug('Consider checking mm = {}, tt = {} again, although they work in this execution they are not guaranteed to work always!'.format(mm,tt))
                        return mm,tt,meantime
                
            logging.debug('Really couldnt find any mm, tt!')
            return None,None,meantime 
        
        def binary_search_mm_tt(self, m, times=5, max_mm=40):
            # Attempt a binary search to find the parameters mm,tt
            mm_list = range(1, max_mm+1)
            first = 0
            last = len(mm_list)-1
            while True:
                midpoint = (first + last)//2
                mm = mm_list[midpoint]
                tt = mm_list[midpoint]+1
                logging.debug('Current mm = {}, tt = {}'.format(mm,tt))
                result, exec_time = self.test_mm_tt(m, mm, tt, times)
                if result and midpoint == 0: #best mm,tt is mm=1,tt=1
                    return mm, tt, exec_time
                if not result and midpoint == 0:
                    return None, None, None #couldn't find any mm,tt combination in the search interval    
                if not result and midpoint == len(mm_list)-1:
                    return None, None, None #couldn't find any mm,tt combination in the search interval    
                if result:
                    logging.debug('Checking is smaller values work mm-1, tt-1')
                    result, exec_time2 = self.test_mm_tt(m, mm-1, tt-1, times) #check if a smaller mm and tt is possible
                    if not result: #we are in the best possible spot
                        return mm, tt, exec_time 
                    else: # a better solution is possible
                        last = midpoint-1 #we search in the first half
                else: 
                    first = midpoint+1 #we search in the second half
                                    
        def calculate_full_m(self, keylen):
            # Returns the full_m to use with a specific key
            if keylen >= 3968:
                return self.primorial(225)  # real bound 3968
            if keylen >= 1984:
                return self.primorial(126)  # 701 is used from 1984 rather than 2048
            if keylen >= 992:
                return self.primorial(71)  # 353 is used from 992 rather than 1024
            if keylen >= 512:
                return self.primorial(39)
            #if keylen >= 128:
            #    return self.primorial(13)
            #if keylen >= 64:
            #    return self.primorial(7)
            return self.primorial(39)  # no data for <512
            
        
        def gen_vuln_key(self, keylen=None, full_m=None):
            # Generate a vulnerable key consisting of modulus n, primes p and q
            # the format of the primes p and q in binary in RSAlib always begins in 0b11XX.. 
            # in this implementation, we generate primes between 0b1100.. to 0b1101..
            # many rsa libraries want primes in the form 0b11.., to ensure the p,q are always of the required size
            
            if keylen < 512 and self.keylen < 512:
                print 'Cannot generate a key smaller than 512 bit, there are no real RoCa keys with that size!'\
                      'if you want you could pass the full_m if you know it uses that specific m'\
                      ' for example if you have a sample key, you could use easily find how it is generated.'
                return None
            if not keylen:
                keylen = self.keylen
            if not full_m:
                full_m = self.calculate_full_m(keylen)
            
            order = self.multOrder(self.generator, full_m, fmpz(full_m).factor())
              
            min_k = 1 << keylen // 2 - 1 | 1 << keylen // 2 - 2  #set msb and msb-1 to 1, the rest is 0
            max_k = min_k | 1 << keylen // 2 - 4 #max k in bin is 0b111000..0XXX.. keylen cannot exceed keylen size (remember that we have to multiply full_m and sum pow)
                   
            min_k = min_k // full_m
            max_k = max_k // full_m
                   
            p=1
            q=1
            start=time.time()
            while not is_prime(p): # or log(p,2) > 255.5: #or log(p,2) > 255:
                k = randint(min_k, max_k) 
                a = randint(0, order-1)
                p = k*full_m + pow(65537,a,full_m)
            end=time.time()
            logging.debug('p prime generation time is {}'.format(end-start))
            logging.debug('p size in bits is {}'.format(log(p,2)))
            
            if keylen % 2 != 0: #keylen is not even, one prime must be 1 bit bigger than the other
            
                min_k = 1 << keylen // 2 | 1 << keylen // 2 - 1
                max_k = min_k | 1 << keylen // 2 - 2        
                min_k = min_k // full_m
                max_k = max_k // full_m
            start=time.time()
            
            while not is_prime(q) or q==p:
                k = randint(min_k, max_k) 
                a = randint(0, order-1)
                q = k*full_m + pow(65537,a,full_m)
            end=time.time()
            logging.debug('q prime generation time is {}'.format(end-start))
            logging.debug('q size in bits is {}'.format(log(q,2)))
            n=p*q
            logging.debug('Key generated is: ')
            logging.debug('Prime p = {}, Prime q = {}, Modulus = {}'.format(p,q,n))
            return n,p,q
                    
        def calculate_m_prime(self, steps=2, times_to_test=5, max_mm = 20, advanced_brute=False):
            # Try to calculate the best m_prime (which has the least execution time)
            # Takes the full_m, steps which is the number of ord_prime to "pop" from
            # the greedym, if advanced brute force is True, these are also the number of
            # ord prime we bruteforce (test all possible ord_primes, which "fit" in the
            # required number of bit for that key). If advanced brute is false, we first
            # test greedym as is, the we "pop" steps=2 ord primes from greedym, and test
            # the resulting solution. This is certainly alot faster, but sometimes you
            # don't get the best solution. times_to_test is the number of times to test the 
            # paramters found, PLEASE don't use a number smaller than 2, because sometimes
            # there are false positive and the combination of m and mm tt don't work on all 
            # keys!
            # max_mm is the maximum mm to search for, the smaller it is, the lesser the 
            # time to find the parameters is, then again, if you don't specify an mm big enough,
            # even if we find a perfect m, we cannot find the parameters mm,tt that work for that m
            # TODO implement a timeout for the testing of parameters, when the time is > 
            # than the best solution in the list, abort the process. Plus you can make it
            # multiprocess (the execution time could be wrong tough)
            
            full_m = self.calculate_full_m(self.keylen)
            full_m_factored = fmpz(full_m).factor()
            greedysol, greedym = self.greedy_heuristic_m_finder(self.generator, full_m)
            goback = self.goback #self.goback is populated by greedy_heuristic
            reward_cost_remain = self.reward_cost #remaining ord prime which are not picked by greedy heuristic
            removal_dict = self.removal_dict #removal dict for all the ord primes

            greedysol.reverse() #putting in the head the last values picked, for convenience's sake
            ord_to_add = []
                       
            m_and_tot_time = [] #contains the tuple (m,total time for that m)
            the_parameters = {} #to store parameters found
            m_list = [] #list containing the m obtained by bruteforcing (or not)
                        
            # if we don't test here this solution will never be checked, if we do bruteforcing it comes up later
            # UNLESS greedysol is empty, which happens for example for rsa-864 bit, hence why the addition of
            # or not greedysol          
            if not advanced_brute or not greedysol: 
                logging.debug('We will test greedym {} as is.'.format(greedym))
                m_list.append(greedym)
                       
            ##### REMOVING "steps=2" ord_primes from the solution #####
            step = 0 # n of ord prime to go back, this increases bruteforce time         
            while greedysol and step < steps:
                removed = greedysol.pop(0) #remove first element in the solution
                logging.debug('Removing {}'.format(removed))
                ord_to_add.append(removed) #appending to the ord_primes to test, only used if we brute force later
                #This is needed because we need to re-add the powers > , they were removed in the greedy heuristic...
                if removed[1] > 1: #power is greater than 1
                    for ord_prime in removal_dict:
                        if ord_prime[0] == removed[0] and ord_prime[1] > removed[1]:
                            ord_to_add.append(ord_prime)
                primes_to_multiply = goback[removed] #primes to re-add to the solution, we are going up
                for prime in primes_to_multiply:
                    logging.debug('Multiplying by {}'.format(prime))
                    greedym *= prime
                logging.debug('We will test this m {}, we are up {} step of a total of {} steps'.format(greedym,step,steps))
                m_list.append(greedym)
                step = step+1
            for ord_prime in reward_cost_remain: #adding all the ord_primes not chosen, so now these are the ones we have to test the combination of
                ord_to_add.append(ord_prime[1])  #this is used only if we brute force later
            ###########################################################
                            
            ##### Advanced brute force of the last 2 ord primes #####
            if advanced_brute:
                primes_in_greedym = []
                for prime,power in fmpz(greedym).factor(): #we need this because we don't want to remove a prime that is NOT there!
                    primes_in_greedym.append(prime)
                                  
                # Now greedym is up 2 steps, ord_to_add contains all the ord_prime combinations we need to test, 
                # and removal_dict the primes corresponding to the ord_primes. Remove ord_primes which are too big 
                # to add to the solution
                max_bits = log((2**self.keylen), 2) / 4
                new_ord_to_add = []
                new_removal_dict = {}
                for elem in ord_to_add:
                    bits = log(reduce(operator.mul, removal_dict[elem], 1),2)
                    if log(greedym,2)-bits > max_bits: #we don't go under the required max_bits if we pick this
                        new_removal_dict[elem] = removal_dict[elem]
                        new_ord_to_add.append(elem)
                removal_dict = new_removal_dict
                ord_to_add = new_ord_to_add
            
                # Adding single orders, effectively bruteforcing 1 ord prime
                for elem in ord_to_add:
                    logging.debug('Adding {} to the solution'.format(elem))
                    tempm = greedym 
                    for prime in removal_dict[elem]:
                        if prime in primes_in_greedym:
                            logging.debug('Dividing by {}'.format(prime))
                            tempm /= prime
                        else:
                            logging.debug('This prime is already removed') #happens when that prime is already removed by something already in solution
                    logging.debug('M is now {}'.format(tempm))
                    if tempm not in m_list:
                        m_list.append(tempm)
                    else: logging.debug('M is already in the list to test') #happens when primes in that order are already removed
                                        
                #Try to add 2 ord prime at the same time, effectively brute forcing the last two ord primes
                for i in xrange(0,len(ord_to_add)):
                    for j in xrange(i+1, len(ord_to_add)):
                        if not (ord_to_add[i][0] == ord_to_add[j][0]): #we are NOT combining (2,1) with(2,2) for example, which is not a possible solution
                            tempm = greedym
                            primes_added = [] #to avoid multplying twice for the same prime
                            logging.debug('Checking if {} and {} are over the required bits'.format(ord_to_add[i],ord_to_add[j]))
                            bits_i = log(reduce(operator.mul, removal_dict[ord_to_add[i]], 1),2)
                            bits_j = log(reduce(operator.mul, removal_dict[ord_to_add[j]], 1),2)
                            if log(tempm,2)-(bits_i+bits_j) > max_bits: #we can test this solution
                                for prime in removal_dict[ord_to_add[i]]:
                                    if prime in primes_in_greedym:
                                        tempm /= prime
                                        logging.debug('divide by {}'.format(prime))
                                    else: logging.info('This prime {} is not in greedym! Lucky we have this check'.format(prime))
                                for prime in removal_dict[ord_to_add[j]]:
                                    if prime in primes_in_greedym:
                                        tempm /= prime
                                        logging.debug('divide by {}'.format(prime))
                                    else: logging.info('This prime {} is not in greedym! Lucky we have this check'.format(prime))
                            if tempm not in m_list:
                                m_list.append(tempm)

            time_and_m = []
            parameters = {}
            logging.info('There are a total of {} m to consider'.format(len(m_list)))
            for m in m_list:
                logging.debug('====Finding parameters for this m {}===='.format(m))
                mm, tt, exec_time = self.binary_search_mm_tt(m, times_to_test, max_mm)
                if mm and tt:
                    order = self.multOrder(self.generator, m, fmpz(m).factor())
                    max_attempts = (order + 1) // 2 + 1
                    time_and_m.append((exec_time * order,m))
                    parameters[m] = ((mm,tt))
            
            time_and_m.sort()
            return time_and_m, parameters
                                

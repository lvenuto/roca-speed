#!/usr/bin/env python2.7
# -*- coding: utf-8 -*

import sys
import math
from fpylll import *
import time
from flint import *
from multiprocessing import Process
from multiprocessing import Pool
from multiprocessing import Queue
import multiprocessing as mp
#import signal

class Roca:

    '''default mm and tt values depending on keylen'''
    mm_tt = { 512: (5, 6),
              128: (2, 3), 
              64 : (3, 4) 
    }

    '''default M values depending on keylen'''

    m_prime = { 512: 0x1b3e6c9433a7735fa5fc479ffe4027e13bea,
                128: 0x6bfc675e31a, #43 bits primorial(12)
                64 : 0x7ca2e  #19 bits primorial(7)
    }

    def __init__(self, n, m=None, mm=None, tt=None, generator=None, k0_guess=None, nprocess=None, batch_size=None):
        keylen = (int(math.log(n, 256)) + 1)*8
        if n is None:
          raise ValueError('n must be specified.')
        else:
          self.n=n #modulus of the public key
        if m is None: #if no m is specified, try using default values
          if keylen in self.m_prime:
              self.m = self.m_prime[keylen]
          else:
              raise ValueError('m not specified and keylen of {} bit is not in the default values'.format(keylen))
        else:
          self.m = m
        if generator is None: #for Roca the generator is 65537, but if it's not roca, you can specify it
          self.generator = 65537
        else:
          self.generator=generator    
        if mm is None: #default mm_tt are precomputed, but if you have a weird keylen, and you computed a different m, you can specify it
          if keylen in self.mm_tt:
              self.mm = self.mm_tt[keylen][0]
          else:
              raise ValueError('mm not specified and default mm to use for keylen of {} bit is not in the default values'.format(keylen))              
        if tt is None:
          if keylen in self.mm_tt:
              self.tt = self.mm_tt[keylen][1]
          else:
              raise ValueError('mm not specified and default tt to use for keylen of {} bit is not in the default values'.format(keylen))     
        if k0_guess is not None: #starting guess, used if you want to specify a starting point, usually for finding best m
          self.k0_guess = k0_guess
        if nprocess is None: #you can specify the number of threads to use, if you want to not stress the cpu choose nprocess < n of cores 
          try:
              self.nprocess = mp.cpu_count() 
          except NotImplementedError:
              'Cant get cpu_count, using the default of 2 processes.'
              self.nprocess = 2
        else:
          self.nprocess = nprocess
        #batch size can be optimized to increase performance, usually depends on the time taken by coppersmith alg
        #to prevent starvation of the processes
        if batch_size is None: 
          self.batch_size = 100
        else:
          self.batch_size = batch_size
    
              
    '''Faster implementation of EuclidExt algorithm?
    Used in the Rosetta Code example'''
    '''
    def mul_inv(a, b):
        b0 = b
        x0, x1 = 0, 1
        if b == 1: return 1
        while a > 1:
            q = a / b
            a, b = b, a%b
            x0, x1 = x1 - q * x0, x0
        if x1 < 0: x1 += b0
        return x1
    '''

    #Computes the multiplicative inverse c mod d
    '''https://secgroup.dais.unive.it/teaching/cryptography/the-rsa-cipher/'''
    def _EuclidExt(self, c, d):
        d0 = d
        e = 1
        f = 0
        while d != 0:
            q = c/d		# integer division
            tmp = c - q*d	# this is c % d
            c = d
            d = tmp
            tmp = e - q*f	# new computation for the inverse
            e = f
            f = tmp
        if c == 1:
	        return e % d0	# if gcd is 1 we have that e is the inverse

    '''Implementation from https://rosettacode.org/wiki/Chinese_remainder_theorem#Python
    To be possibly replaced by the FLINT implementation, not much speed to gain though'''

    def _chinese_remainder(self, n, a):
        sum = 0
        prod = reduce(lambda a, b: a*b, n)
     
        for n_i, a_i in zip(n, a):
            p = prod / n_i
            sum += a_i * self._EuclidExt(p, n_i) * p
        return sum % prod
        
    '''Coppersmith_howgrave implementation using FLINT polynomials, because Sympy polynomials are 100 times slower
       Doesn't work if you work with negative numbers at the moment, the polynomials are not in mod n
       For the Sage implementation goto: https://github.com/mimoo/RSA-and-LLL-attacks
       But it's 40% slower partly because of LLL implementation in Sage, partly because the factorization here is faster (why?)
    '''

    def coppersmith_howgrave_univariate(self, expr, n, mm, tt, XX):

        dd = expr.degree()
        nn = dd * mm + tt
             
        start = time.time()
        # compute polynomials
        gg = []
        for ii in range(mm):
            for jj in range(dd):
                gg.append((fmpz_poly([0,XX])**jj) * fmpz_poly([n**(mm-ii)]) * (expr(fmpz_poly([0,XX]))**ii) )
        for ii in range(tt):
            gg.append((fmpz_poly([0,XX])**ii) * (expr(fmpz_poly([0,XX]))**mm) )
       
        BB = IntegerMatrix(nn,nn)
        
        for ii in range(nn):
            for jj in range(ii+1):
                BB[ii, jj] = int(gg[ii][jj])
        # LLL
        W=LLL.Wrapper(BB,delta=0.5)
        W()
          
        # transform shortest vector in polynomial    
        new_pol = fmpz_poly([0])
       
        for ii in range(nn):
            new_pol += (fmpz_poly([0,1]))**ii * fmpz(BB[0, ii] / XX**ii)

        roots = []
        
        factorization = new_pol.factor()
        
        for elem in factorization:
            if isinstance(elem,type(fmpz(2))): 
                if elem != -1 and elem != 1:
                    roots.append(elem)
            if isinstance(elem,list):
                for elem2 in elem:
                    if(elem2[0].degree() == 1):
                        roots.append((-elem2[0][0]))               
        end = time.time()
        exec_time = end-start #exec time is useful if you are searching for a perfect m, normally not useful
        return roots, exec_time 
        
    def try_guess(self, n, m, k0_guess, k0_guess_times_m_inv, mm, tt, XX):
        
        k0_g_inv_mod = k0_guess_times_m_inv%n
        expr = fmpz_poly([k0_g_inv_mod,1])
        roots, exec_time = self.coppersmith_howgrave_univariate(expr, n, mm, tt, XX)
        for root in roots:
            factor1 = k0_guess + abs(root) * m
            if (n % factor1) == 0: #found our prime p
                print 'Correct guess is: ' + str(k0_guess)
                #print 'ABS root is: ' + str(abs(root)) 
                factor2 = n // factor1
                return [int(factor1), int(factor2)]
        return None
        
    '''TODO: Correct typing for use with FLINT library'''
    def generator_order(self, generator, m, phi_n, decomp_phi_n):
        order = int(phi_n)
        if generator == 1:
            return 1 # by definition
        if pow(generator, order, m) != 1:
            return None #not an element of the group
        #for factor,power in decomp_phi_n.items(): for sympy dictionary use .items()
        for factor,power in decomp_phi_n:
            for p in range (1, power + 1):
                next_order = order // int(factor)
                if pow(generator, next_order, m) == 1:
                    order = next_order           
                else:
                    break 
        #print 'order is : ' + str(order)   
        return order

    #Pohlig hellman algorithm for an effecient computation of the discrete logarithm
    def pohlig_hellman(self, modulus, generator, generator_order, generator_order_decomposition, m):
        
        if pow(modulus, generator_order, m) != 1:
            return None
        moduli = []
        remainders = []
        for my_tuple in generator_order_decomposition:
            prime = int(my_tuple[0])
            power = int(my_tuple[1])
            prime_to_power = prime ** power
            order_div_prime_power = generator_order // prime_to_power
            g_dash = pow(generator, order_div_prime_power, m)
            h_dash = pow(modulus, order_div_prime_power, m)
            found=False
            for i in xrange(0, prime_to_power):
                if pow(g_dash, i, m) == h_dash:
                    remainders.append(i)
                    moduli.append(prime_to_power)
                    found = True
                    break
            if not found:
                return None
        #start = time.time()
        #result = crt(moduli, remainders, symmetric=False)[0]
        result = self._chinese_remainder(moduli,remainders)
        #end = time.time()
        #print "crt time is: " + str(end-start) 
        return result

    '''This function is primarily used for the processes spawned by this class
    they get a batch of n elements, and perform operations. A multiprocessing
    queue is used to share data between them.

    TODO: do a checkpoint to resume (E.g. for amazon EC2 spot instance resume)
    and put the max_attempts in the queue, so we know how many have been done
    '''
    def get_batch(self, lock, queue, n, m, m_inv, generator, mm, tt, XX):
        batch = []
        #global max_attempts    
        lock.acquire()
        next_k0_to_assign = queue.get()["next_k0_to_assign"]
        attempts = queue.get()["attempts"] 
        if attempts > self.max_attempts:
            #Restore queue
            queue.put({"next_k0_to_assign" : next_k0_to_assign})
            queue.put({"attempts" : attempts})
            lock.release()
            return None #no more jobs for you
        k0_guess = next_k0_to_assign
        k0_guess_times = int(k0_guess * m_inv)
        batch.append((k0_guess, k0_guess_times))
        for i in xrange(0, self.batch_size-1): #we have already assigned one job
            k0_guess = ((k0_guess) * generator) % m 
            k0_guess_times = int(k0_guess * m_inv)
            batch.append((k0_guess,k0_guess_times))
        attempts = attempts + len(batch)
        next_k0_to_assign = ((k0_guess) * generator) % m
        queue.put({"next_k0_to_assign" : next_k0_to_assign})
        queue.put({"attempts" : attempts})
        lock.release()
        #print attempts
        return batch
       
    def __separate_process(self, proc_id, lock, stop_event, queue, q_solution, n, m, m_inv, generator, mm, tt, XX):
        batch = self.get_batch(lock, queue, n, m, m_inv, generator, mm, tt, XX)
        while batch: #while there is work to do
            for guess in batch:
                if(stop_event.is_set()):
                    print "Process n " + str(proc_id) + " has been requested to terminate, exiting now"
                    exit(0)
                factors = self.try_guess(n, m, guess[0], guess[1], mm, tt, XX) #guess[0] = k0_guess, guess[1]= k0_guess_times_m_inv
                if factors is not None: #we found the solution
                    print "Process n " + str(proc_id) + " found the solution!" 
                    stop_event.set()
                    q_solution.put({"solution" : factors})
                    exit(0)
            batch = self.get_batch(lock, queue, n, m, m_inv, generator, mm, tt, XX)
        print "Process n " + str(proc_id) + " found nothing, exiting now"

    def factorize(self):
        n = self.n
        m = self.m
        mm = self.mm
        tt = self.tt
        generator = self.generator
        general_prime = False #Should be true for rsalib keys, check this
        #generator is known: p=k * m + (generator^a mod m)
        phi_n = fmpz.euler_phi(fmpz(m)) #Euler totient in FLINT
        decomp_phi_n = phi_n.factor() # list with prime factors & multiplicities
        order = self.generator_order(generator, m, phi_n, decomp_phi_n)
        decomp_order = fmpz(order).factor() 
        d = self.pohlig_hellman(n, generator, order, decomp_order, m)
        guess = d // 2
        print 'guess is : ' + str(guess)
        self.max_attempts = ((order + 1) // 2 + 1)
        print 'Max_attemps are: ' + str(self.max_attempts)
        m_inv = int(self._EuclidExt(m, n))
        length = int(math.ceil(math.log(n, 2)))
        if general_prime:
            # any prime of |n|/2 bits
            XX = int(2**(length / 2)) // m
        else:
            # primes of the form 0b1100...
            XX = int(2**(length / 2 - 1) + 2**(length / 2 - 2) + 2**(length / 2 - 4)) // m
      
        k0_guess = int(pow(generator, guess, m))
        lock = mp.Lock()
        stop_event = mp.Event()
        self.stop_event = stop_event
        queue = mp.Queue()
        q_solution = mp.Queue()
        queue.put({"next_k0_to_assign" : k0_guess})
        queue.put({"attempts" : 0})
     
        processes = [mp.Process(target=self.__separate_process, args=(
        proc_id, lock, stop_event, queue, q_solution, n, m, m_inv, generator, mm, tt, XX)) 
        for proc_id in xrange(1,self.nprocess+1)]
        
        try:
            start = time.time()
            print("Main Thread starting the processes")
            for p in processes:
                p.start()
            for p in processes:
                p.join()
            print("All processes terminated.")
            end = time.time()
            print "Real time of execution: " + str(end-start)
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, asking processes to terminate")
            for p in processes:
                p.terminate()
            exit(0)
        next_k0_to_assign = queue.get()["next_k0_to_assign"]
        total_attempts = queue.get()["attempts"]
        print "Total attempts are: " + str(total_attempts)
        print "Attempts per second: " + str(int(total_attempts / (end-start))) 
        #Old stuff but important for the writeup, shit doesn't work correctly
        #results = [pool.apply(try_guess, args=(n, m, k0_guess_list[w], k0_guess_times_list[w], mm, tt, XX)) for w in xrange(0,n_process)]
        #results = [p.get() for p in results]   
        if not q_solution.empty():
            solution = q_solution.get()["solution"]
            return solution
        return None
        
        

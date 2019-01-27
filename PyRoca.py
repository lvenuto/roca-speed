#!/usr/bin/env python2.7
# -*- coding: utf-8 -*

import sys
import math
import time
from fpylll import *
from flint import *
from multiprocessing import Process
from multiprocessing import Pool
from multiprocessing import Queue
import multiprocessing as mp
import logging
#import signal

class PyRoca:

    '''default mm and tt values depending on keylen'''
    mm_tt = { 
              1024: (4, 5),
              512: (5, 6),
    }

    '''default M values depending on keylen'''

    m_prime = { 
                1024 : 0x24683144f41188c2b1d6a217f81f12888e4e6513c43f3f60e72af8bd9728807483425d1eL,
                512 : 0x1b3e6c9433a7735fa5fc479ffe4027e13bea,            
    }

    def __init__(self, n, m=None, mm=None, tt=None, generator=None, general_prime=True, guess=None, nprocess=None, batch_size=None):
       
        #keylen = (int(math.log(n, 256)) + 1)*8
        keylen = n.bit_length()
        logging.debug('keylen is {}'.format(keylen))
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
        if mm is None: #some default mm_tt are precomputed, but if you have a weird keylen, and you computed a different m, you can specify it
            if keylen in self.mm_tt:
                self.mm = self.mm_tt[keylen][0]
            else:
                raise ValueError('mm not specified and default mm to use for keylen of {} bit is not in the default values'.format(keylen))   
        else:
            self.mm=mm           
        if tt is None:
            if keylen in self.mm_tt:
                self.tt = self.mm_tt[keylen][1]
            else:
                raise ValueError('mm not specified and default tt to use for keylen of {} bit is not in the default values'.format(keylen)) 
        else:
            self.tt=tt    
        
        if general_prime: #determines the size of XX, check the factorize function for more info
            self.general_prime = True 
        else: self.general_prime = False
        
        if guess is not None: #starting guess, used if you want to specify a starting point, usually for finding best m
            self.guess = guess
        else: self.guess = None
        
        if nprocess is None: #you can specify the number of threads to use, if you want to not stress the cpu choose nprocess < n of cores 
            try:
                self.nprocess = mp.cpu_count() 
            except NotImplementedError:
                logging.warning('Cant get cpu_count, using the default of 1 process.')
                self.nprocess = 1
        else:
            self.nprocess = nprocess
        #batch size can be optimized to increase performance, usually depends on the time taken by coppersmith alg
        #to prevent starvation of the processes
        if batch_size is None: 
            self.batch_size = 100
        else:
            self.batch_size = batch_size
    
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
  
    @classmethod
    def EuclidExt(cls, c, d):
    # Computes the multiplicative inverse c mod d
    # reference: https://secgroup.dais.unive.it/teaching/cryptography/the-rsa-cipher/   
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
    
    @classmethod
    def chinese_remainder(cls, n, a):
    # Implementation from https://rosettacode.org/wiki/Chinese_remainder_theorem#Python
    # To be possibly replaced in the future by the FLINT implementation, not much speed to gain though
        sum = 0
        prod = reduce(lambda a, b: a*b, n)
     
        for n_i, a_i in zip(n, a):
            p = prod / n_i
            sum += a_i * cls.EuclidExt(p, n_i) * p
        return sum % prod
        
    @staticmethod
    def coppersmith_howgrave_univariate(expr, n, mm, tt, XX):
    # Coppersmith_howgrave implementation using FLINT polynomials, because Sympy polynomials are more than an order of magnitude slower
    # Doesn't work if you work with negative numbers at the moment, the polynomials are not modulo n
    # For the Sage implementation the reference is: https://github.com/mimoo/RSA-and-LLL-attacks
    # But it's 20% slower partly because of LLL implementation in Sage, partly because the factorization with FLINT in Python is faster
        dd = expr.degree()
        nn = dd * mm + tt

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
                  
        return roots 
        
    def try_guess(self, n, m, k0_guess, mm, tt, XX):
    # Given a guess and parameters m,mm,tt,XX run coppersmith attack
        k0_g_inv_mod = (k0_guess * self.m_inv) % n
        expr = fmpz_poly([k0_g_inv_mod,1])
        roots = self.coppersmith_howgrave_univariate(expr, n, mm, tt, XX)
        for root in roots:
            factor1 = k0_guess + abs(root) * m
            if (n % factor1) == 0: #found our prime p
                logging.debug('Correct guess is: {}'.format(k0_guess))
                factor2 = n // factor1
                return [int(factor1), int(factor2)]
        return None
        
    def generator_order(self, generator, m, phi_n, decomp_phi_n):
    # TODO: Correct typing for use with FLINT library
        order = int(phi_n)
        if generator == 1:
            return 1 # by definition
        if pow(generator, order, m) != 1:
            return None #not an element of the group
        for factor,power in decomp_phi_n:
            for p in range (1, power + 1):
                next_order = order // int(factor)
                if pow(generator, next_order, m) == 1:
                    order = next_order           
                else:
                    break 
        return order
    
    @staticmethod
    def pohlig_hellman(modulus, generator, generator_order, generator_order_decomposition, m):
    # Pohlig hellman algorithm for an effecient computation of the discrete logarithm
        
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
        result = PyRoca.chinese_remainder(moduli,remainders)
        return result

    def get_batch(self, lock, queue, n, m, generator, mm, tt, XX):
    # This function is primarily used for the processes spawned by this class
    # they get a batch of n elements, and perform operations. A multiprocessing
    # queue is used to share data between them.
    # TODO: do a checkpoint to resume (E.g. for amazon EC2 spot instance resume)
    # and put the max_attempts in the queue, so we know how many have been done
    
        batch = []
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
        #k0_guess_times = int(k0_guess * m_inv)
        batch.append(k0_guess)
        for i in xrange(0, self.batch_size-1): #we have already assigned one job
            k0_guess = ((k0_guess) * generator) % m 
            #k0_guess_times = int(k0_guess * m_inv)
            #batch.append((k0_guess,k0_guess_times))
            batch.append(k0_guess)
        attempts = attempts + len(batch)
        next_k0_to_assign = ((k0_guess) * generator) % m
        queue.put({"next_k0_to_assign" : next_k0_to_assign})
        queue.put({"attempts" : attempts})
        lock.release()
        return batch
       
    def __separate_process(self, proc_id, lock, stop_event, queue, q_solution, n, m, generator, mm, tt, XX):
    # Every process will run this function until finished or asked to terminate
        batch = self.get_batch(lock, queue, n, m, generator, mm, tt, XX)
        while batch: #while there is work to do
            for guess in batch:
                if(stop_event.is_set()):
                    lock.acquire()
                    next_k0_to_assign = queue.get()["next_k0_to_assign"]
                    total_attempts = queue.get()["attempts"]
                    queue.put({"next_k0_to_assign" : next_k0_to_assign})
                    queue.put({"attempts" : total_attempts-len(batch)})#use this if you want to know (more or less) where the correct guess is in the batches
                    #queue.put({"attempts" : total_attempts-(batch.index(guess)+1)}) if you want more accurate attempts/sec, uncomment this
                    lock.release()         
                    logging.debug('Process n {} has been requested to terminate, exiting now'.format(proc_id))
                    exit(0)
                factors = self.try_guess(n, m, guess, mm, tt, XX) 
                if factors is not None: #we found the solution
                    logging.debug('Process n {} found the solution!'.format(proc_id))
                    stop_event.set() #asking others to terminate
                    q_solution.put({"solution" : factors}) 
                    #Fixing total attemps made
                    lock.acquire()
                    next_k0_to_assign = queue.get()["next_k0_to_assign"]
                    total_attempts = queue.get()["attempts"]
                    total_attempts = total_attempts - (len(batch) - (batch.index(guess)+1)) #len(batch)-remaining items in batch
                    queue.put({"next_k0_to_assign" : next_k0_to_assign})
                    queue.put({"attempts" : total_attempts})
                    lock.release()
                    exit(0)
            batch = self.get_batch(lock, queue, n, m, generator, mm, tt, XX)
        logging.debug('Process n {} found nothing, exiting now'.format(proc_id))

    def factorize(self):
        n = self.n
        m = self.m
        mm = self.mm
        tt = self.tt
        generator = self.generator
        #generator is known: p=k * m + (generator^a mod m)
        phi_n = euler_phi(fmpz(m)) #Euler totient in FLINT
        decomp_phi_n = phi_n.factor() # list with prime factors & multiplicities
        order = self.generator_order(generator, m, phi_n, decomp_phi_n)
        decomp_order = fmpz(order).factor() 
        d = self.pohlig_hellman(n, generator, order, decomp_order, m)
        logging.debug('Order of the generator is: {}'.format(order))    
        if self.guess:
            guess = self.guess #if we have specified a starting guess, useful for resuming
        else:
            guess = d // 2
        logging.debug('Starting guess is: {}'.format(guess))
        self.max_attempts = ((order + 1) // 2 + 1)
        logging.debug('Max attemps are: {}'.format(self.max_attempts))
        self.m_inv = int(self.EuclidExt(m, n)) 
        length = int(math.ceil(math.log(n, 2)))
        if self.general_prime:
            # any prime of |n|/2 bits
            XX = int(2**(length / 2)) // m
        else:
            #used in vulnerable keys found in rsalib, as far as I know, it's a small optimization anyway
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
        proc_id, lock, stop_event, queue, q_solution, n, m, generator, mm, tt, XX)) 
        for proc_id in xrange(1,self.nprocess+1)]
        
        try:
            start = time.time()
            logging.info("Main Thread starting the processes")
            for p in processes:
                p.start()
            for p in processes:
                p.join()
            logging.info("All processes terminated.")
            end = time.time()
            logging.debug("Real time of execution: {}".format(end-start)) 
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating processes")
            for p in processes:
                p.terminate()
            exit(0)
        next_k0_to_assign = queue.get()["next_k0_to_assign"]
        total_attempts = queue.get()["attempts"]
        logging.info('Total attempts made are: {}'.format(total_attempts)) #not really true, it's more or less right depending on the scheduler...
        logging.info('Attempts per second: {}'.format(total_attempts / (end-start)))
        #Old stuff but important for the writeup, this doesn't work correctly
        #results = [pool.apply(try_guess, args=(n, m, k0_guess_list[w], k0_guess_times_list[w], mm, tt, XX)) for w in xrange(0,n_process)]
        #results = [p.get() for p in results]   
        if not q_solution.empty():
            solution = q_solution.get()["solution"]
            return solution
        return None
        
        

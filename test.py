#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
import PyRoca as Rc
import ParameterFinder

if __name__ == "__main__":
    
    n1024 = 0x9af4cf8a8f5aa9faa6c15439c864d984f29246f038df6853c1e3397d9dd5b693f5f12e21014f5ad6c2b09978ffa6a4f7ccd03a134346ed708316c2f889354ecfaff37e7cf50db9c9497e928b5795cbc9ac8a926adce6ef9b1835769674b9ca2e6e2e0a0cc7bc65090bade141a776c16551f8926b562dd0aa9d12c38bf151f041L
    n512 = 0x958f12d6d12af295d46c6094ba292f020e77befdef0ac12b0c9b8c624fcba59fa36bc3e95a850eb3824d8ef4dac4a8f74b6d824151262d43913b0b9313ee0945
    n128 = 0x91d4b818e149fe83c7b97e4afce1e213
    n64 = 0xa5d49535af970e69
    
    #print roca.n
    
    
    print 'Factoring RSA-64 bit key'
    roca = Rc.PyRoca(n64, debug=True)
    print 'The two primes p and q are: {}'.format(roca.factorize())
    print '\n'
    print 'Factoring RSA-128 bit key'
    roca = Rc.PyRoca(n128, debug=True)
    print 'The two primes p and q are: {}'.format(roca.factorize())
    
    print '\n'
    print 'Factoring RSA-512 bit key'
    roca = Rc.PyRoca(n512, batch_size=100, nprocess=4, debug=True)
    print 'The two primes p and q are: {}'.format(roca.factorize())
  
    
    #PF = ParameterFinder.ParameterFinder(4096, debug=True)
    #result, full_m = PF.greedy_heuristic_m_finder(65537,PF._primorial(225))
    #print "Result is: {}".format(result)
    #print "full_m is: {}".format(full_m)

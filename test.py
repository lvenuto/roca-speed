#!/usr/bin/env python2.7
# -*- coding: utf-8 -*
import roca_class

if __name__ == "__main__":
    n512 = 0x958f12d6d12af295d46c6094ba292f020e77befdef0ac12b0c9b8c624fcba59fa36bc3e95a850eb3824d8ef4dac4a8f74b6d824151262d43913b0b9313ee0945
    n128 = 0x91d4b818e149fe83c7b97e4afce1e213
    n64 = 0xa5d49535af970e69
    #print roca.n
    print 'Factoring RSA-64 bit key'
    roca = roca_class.Roca(n64)
    print 'The two primes p and q are: {}'.format(roca.factorize())
    print '\n'
    print 'Factoring RSA-128 bit key'
    roca = roca_class.Roca(n128)
    print 'The two primes p and q are: {}'.format(roca.factorize())
    print '\n'
    print 'Factoring RSA-512 bit key'
    roca = roca_class.Roca(n512)
    print 'The two primes p and q are: {}'.format(roca.factorize())
    
    
    

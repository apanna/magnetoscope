# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 14:40:14 2012

@author: schlammi
"""
from __future__ import division
import math
import sys




def CalcTwoPoleCoeff(f0,fs,n,feq='Butterworth',ftype='lp'):
    """
    f0: is the cut off frequency
    fs: The sample frequency
    n: How often is the filter applied consecutively
    ftype: Up to now only 'Butterworth' is implemented
    
    Note, this function uses internally a different 
    notation.they need to beconverted at the return
    """

    # general form:         H(s)=g/(s^2+ps+g)
    # Butterworth:          H(s) = 1/(s^2 + sqrt(2) s +1)
    # Critically damped"   H(s) = 1/(s^2 + 2 s +1)
    # Bessel:                H(s) = 3/(s^2 + 3 s +3)

    if feq.upper()=="BUTTERWORTH" :
        g=1.0
        p=math.sqrt(2.0)
    elif feq.upper()=="CRITICAL" :
        g=1.0
        p=2.0
    elif feq.upper()=="BESSEL" :
        g=3.0
        p=3.0
        
    

    c2 = 2.0/(2.0*g-p*p+((2.0*g-p*p)**2-4*g*g*(1-2**(1.0/n)))**0.5)
    c = math.sqrt(c2)

    filterStable = True
    
    if ftype.upper()=="LP":
        fstar = c * f0/fs
        if fstar<1.0/8.0:
                filterStable = True
    
    else:    
        fstar = 1.0/2.0 - c * f0/fs            
        if fstar<1.0/2.0 and fstar > 3.0/8.0:
                filterStable = True
        
    print fstar,filterStable
    
    if ftype.upper()=="LP":
        o0 = 1.0*math.tan(math.pi *c *f0/fs)
    else:     
        c=1.0/c
        o0 = 1.0*math.tan(math.pi *(1.0/2.0-c *f0/fs))
        
    a0 = g*o0*o0/(1.0+p*o0+g*o0*o0)
    a1 = 2.0*a0
    a2 = a0
    b1 = 2.0*a0*(1.0/(g*o0*o0)-1.0)
    b2 = 1.0- (a0+a1+a2+b1)
    if ftype.upper()=="HP":
        a1 = -a1
        b1 = -b1
        
    # Note what is here b, is in reality -a
    return ((-b1,-b2),(a0,a1,a2))


def convert_2nd_order(A,B,f0,fs):
    """
    converts a Laplace Transfer function in the form:
    F(s) = ( B0+B1*S+B2*S*S  ) / (1+A1*S+A2*S*S)
    into a z-transfer function:
        
                         -1     -2
                b0 + b1 z + b2 z  
        H(z) = -----------------------
                         -1        -2 
                1 + a1 z     + a2 z   
        
    """
    
    A1=A[0]
    A2=A[1]
    B0=B[0]
    B1=B[1]
    B2=B[2]
    
    T = 1.0/fs
    omega_c = 2.0* math.pi*f0
  
    # To calculate the a's and b's we use an algorithm described in
    # D. Schlichthaertle, Digital Filters, p. 166
    # converts H(S)=(B0+B1*S+B2*S*S)/(1+A1*P+A2*P*P) into z

    C1 = 2.0*A1 /(omega_c * T)
    C2 = 4.0*A2 /(omega_c * T)/(omega_c * T)
    D1 = 2.0*B1/(omega_c * T)
    D2 = 4.0*B2/(omega_c * T)/(omega_c * T)
    
    V = (B0+D1+D2)/(1.0+C1+C2)
    b0 = 1.0
    b1 = 2.0*(B0-D2)/(B0+D1+D2)
    b2= (B0-D1+D2)/(B0+D1+D2)
    a1 = 2.0*(1.0-C2)/(1.0+C1+C2)
    a2 = (1.0-C1+C2)/(1.0+C1+C2)
    
    return ((a1,a2),(b0*V,b1*V,b2*V))

        

def CalcChebyshev(eps,f0,fs):
    """
    Note, so far only implemented for N=2
    eps is the pass band riple, 0.1 is a good value
    """
    cof = gcheby(eps,f0,fs,2)
    c=cof[0]
    return (c[0],c[1])



def gcheby(eps,f0,fs,N):
    """
     Calculates the filter coefficients for a Chebychef low pass
     eps is the amount of passband riple
     f0 is the corner frequency
     fs is the sample frequency
     N is the filter order. Currently N has to be of even order
    """
    if N % 2!= 0:
        print "Program does not work for odd Ns"
        sys.exit()
        
    coeffs=[]
    nu0 = math.asinh(1.0/eps)/N
    V=1.0/math.sqrt(1.0+eps*eps)
    for j in range(N,N+N//2):
        k=j
        re_a =   math.sinh(nu0)*math.sin(math.pi*(2*k+1)/2/N)
        im_a =   math.cosh(nu0)*math.cos(math.pi*(2*k+1)/2/N)
        #print "first pole: ",k,re_a,im_a
        #print nu0
        #k=2*N-1-(j-N)
        #re_b =   math.sinh(nu0)*math.sin(math.pi*(2*k+1)/2/N)
        #im_b =   math.cosh(nu0)*math.cos(math.pi*(2*k+1)/2/N)
        #print "second pole: ",k,re_b,im_b
        re_b=re_a
        im_b=-im_a
        ss = re_a*re_b-im_a*im_b
        if j==N:
            B0 = 1.0#V/ss
        else:
            B0 = 1.0#1.0/ss
        B1 = 0
        B2 = 0
        A1 = -2*re_a/ss
        A2 = 1.0/ss
        print "B0 = ",B0,"  A1 = ",A1,"     A2 = ",A2        
        a,b = convert_2nd_order((A1,A2),(B0,B1,B2),f0,fs)
        coeffs.append((a,b))
    return coeffs

def gbutter(f0,fs,N):
    """
     Calculates the filter coefficients for a butterworth low
     pass filter
     f0 is the corner frequency
     fs is the sample frequency
     N is the filter order. Currently N has to be of even order
    """

    if N % 2!= 0:
        print "Program does not work for odd Ns"
        sys.exit()
        
    coeffs=[]
    for j in range(N//2):
        k=j
        re_a =   -math.sin(math.pi*(2*k+1)/2/N)
        im_a =   math.cos(math.pi*(2*k+1)/2/N)
        #print "first pole: ",k,re_a,im_a
        #print nu0
        #k=2*N-1-(j-N)
        #re_b =   math.sinh(nu0)*math.sin(math.pi*(2*k+1)/2/N)
        #im_b =   math.cosh(nu0)*math.cos(math.pi*(2*k+1)/2/N)
        #print "second pole: ",k,re_b,im_b
        re_b=re_a
        im_b=-im_a
        ss = re_a*re_b-im_a*im_b
        B0 = 1.0/ss
        B0 = 1.0/ss
        B1 = 0
        B2 = 0
        A1 = -2*re_a/ss
        A2 = 1.0/ss
        a,b = convert_2nd_order((A1,A2),(B0,B1,B2),f0,fs)
        coeffs.append((a,b))
    return coeffs
    

class IIR:
    """
    An infinite impulse respose filter
    y(n)=b0*x(n)+b1*x(n-1)+b2*x(n-2)+...
            -a1*y(n-1)-a2*y(n-2)-.....
            
    Response function of the filter is:
                         -1      -2
                b0 + b1 z + b2 z   + ....
        H(z) = ------------------------
                         -1        -2 
                1 + a1 z     + a2 z   + ...
        
    """
    def __init__(self,a,b):
        self.a = a
        self.b = b
        self.inData=[]
        self.outData=[]
        for i in self.b:
            self.inData.append(0.0)
        for i in self.a:
            self.outData.append(0.0)
    
    def filtered(self,new):
        self.inData.pop(-1)
        self.inData.insert(0,new)
        mysum=0
        for x,y in zip(self.inData,self.b):
            mysum += x*y
        for x,y in zip(self.outData,self.a):
            mysum -= x*y
        self.outData.pop(-1)
        self.outData.insert(0,mysum) 
        return mysum
    
    def filterArray(self,inA):
        outA=[]
        for i in inA:
            outA.append(self.filtered(i))
        return outA
        
        
        
class nIIR:
    
    """
    A class that applies an IIR filter repeadetly
    """
    def __init__(self,a,b,n):        
        self.a = a
        self.b = b
        self.myFilterBank=[]
        for i in range(n):
            self.myFilterBank.append(IIR(a,b))
    
    def filtered(self,new):
        myinput =new
        for i in self.myFilterBank:
            myinput = i.filtered(myinput)
        return myinput
        
    
    def filterArray(self,inA):
        outA=[]
        for i in inA:
            myinput = i
            for j in self.myFilterBank:
                myinput = j.filtered(myinput)
            outA.append(myinput)
        return outA
        
class gIIR:
    
    """
    A class that applie different IIR filters
    """
    def __init__(self,coeff):        
        self.myFilterBank=[]        
        for i in coeff:
            a=i[0]
            b=i[1]
            self.myFilterBank.append(IIR(a,b))
    
    def filtered(self,new):
        myinput =new
        for i in self.myFilterBank:
            myinput = i.filtered(myinput)
        return myinput
        
    
    def filterArray(self,inA):
        outA=[]
        for i in inA:
            myinput = i
            for j in self.myFilterBank:
                myinput = j.filtered(myinput)
            outA.append(myinput)
        return outA
        
        
        

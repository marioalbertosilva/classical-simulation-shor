import numpy as np
from fractions import Fraction
import random
import math
import cmath
import time


"DETERMINE NUMBER OF QUBITS"

def size_input(N):
    "Determines the dimensions t,n of the initial state."
    n = math.ceil(math.log(N, 2))
    t = math.ceil(2*math.log(N, 2))
    return n, t


def obtain_x(Z):
    "Outputs a random x between 1 and N such that gcd(x,N)=1"
    y = random.randint(2,Z-1)
    while math.gcd(y,Z) != 1:
        y = random.randint(2,Z-1)
    return y


"QUANTUM OPERATIONS AND MEASUREMENTS"

def mods_xn_and_measure(x,N,two_t):
    "Runs operator V_x until cycle is found, then a value is measured and all respective states are returned"
    xmodn=[]
    xmodn.append(1)
    temp=1
    while True:
        temp=(temp*x)%N
        if temp==1:
            break
        xmodn.append(temp)
    print("Measuring state")
    measure=random.randint(0,len(xmodn)-1)
    print("Measured:",xmodn[measure])
    collapse = range(measure,two_t,len(xmodn))
    result = np.zeros((two_t,1))
    for i in collapse:
        result[i]=1
    return result


def nonuniform_measure(prob):
    "Measures a value from a normalized probability distribution"
    r=random.random()
    for i in range(len(prob)):
        r-=prob[i]
        if r<=0:
            return i
    return 0




"CALCULATE CONTINUED FRACTIONS"

def list_convergents(x,y):
    "Returns the list of all convergents of the fraction x/y"
    z=find_continued_fraction(x,y)
    K=[]
    m=len(z)
    while m!=0:
        K+=[get_continued_fraction(z,m)]
        m-=1
    return K


def find_continued_fraction(x,y,Z=[]):
    "Returns the list of values in the continued fraction of x/y"
    q=y//x
    Z+=[q]
    if y-q*x!=1 and y-q*x!=0:
        return find_continued_fraction(y-q*x,x,Z)
    else:
        Z+=[x]
        return Z

def get_continued_fraction(Y,n):
    "Returns the fraction associated with a set Y of values in a continued fraction, using the first n values of Y or all values if n>len(Y)"
    Y=Y[0:n]
    if len(Y)==1:
        return Fraction(1,Y[0])
    else:
        return Fraction(1,Y[0]+get_continued_fraction(Y[1:],n))




"FAST FOURIER TRANSFORM"

def fft( states ):
    " Computes the Fast Fourier Transform, using the Cooley-Tukey method, "
    " that runs in O(N log N) time. "
    
    N = len(states)
    
    if N == 1 :
        return states
    
    else :
        W = 1
        phi = (2 * math.pi ) / N
        Wn = complex(math.cos(phi), math.sin(phi))
        Aeven = states[0::2]
        Aodd = states[1::2]
        Yeven = fft(Aeven)
        Yodd = fft(Aodd)
        Y=np.empty((N,1),dtype=complex)
        for j in range(0, N//2):
            Y[j] = Yeven[j] + W * Yodd[j]
            Y[j + N//2] = Yeven[j] - W * Yodd[j]
            W = W * Wn

    return Y


"PROBABILITIES COMPUTATION"

def allprobs(dftvals, twot):
    " Computes the PDF for the values of j, follows the DFT. "
    probs = np.empty((twot,1))
    total_prob=0 #account for approximation errors
    for i in range(0,twot):
        #print(i,dftvals[i])
        probs[i]=abs(dftvals[i])**2
        total_prob+=probs[i]
    return probs/total_prob



"INPUT TESTS"

def PrimeQ(n):
    """Primality Test: fast probabilistic test then deterministic test"""

    #Fermat's primality test for base 2
    a= 2
    if (a**(n-1))%n!=1:
        return False

    #Deterministic test
    if n == 2:return True
    if n == 3:return True
    if n % 2 == 0:return False
    if n % 3 == 0:return False
    i = 5
    w = 2
    while i * i <= n:
        if n % i == 0:return False
        i += w
        w = 6 - w
    return True

def PowerPrimeQ(n):
    "Returns True iff n is a power of a prime"
    power = int(math.log(n,2))
    if 2**power==n: # n is a power of 2
        print("%i = 2^%i" % (n,power))
        return True
    for i in range(1,power+1):
        #calculate root i of n until i is integer
        if i==1:m=n
        elif i==2:m=int(math.sqrt(n))
        else:m=int(math.exp(math.log(n)/i))
        if n==(m**i):
            if PrimeQ(m):
                if i!=1:print("%i = %i^%i" % (n,m,i))
                else:print("Prime.")
                return True
            elif i!=1:return False
    return False


"SHOR'S ALGORITHM"

def Shor(N):
    
    begin=time.time()

    n,t = size_input(N)
    two_t = 2**t
    
    print("Number of qubits =", n+t)

    done= False # is True when the factorization is done
    tries=0 # number of tries
    T=10 # maximum number of tries

    while not done and tries<T:
        print("===================================\n")
        tries += 1
        print("Attempt:",tries)

        x = obtain_x(N)

        print("Picked x =", x)

        print("Computing x^j mod N")

        measured_states= mods_xn_and_measure(x,N,two_t)

        print("Running Fast Fourier Transform")
        fftRes = fft(measured_states)
        
        print("Calculating proabilities")
        pdf = allprobs(fftRes, two_t)

        #Plot the probability distribution after Fourier Transform
        #plt.plot(pdf)
        #plt.show()

        print("Measuring state")

        measured= nonuniform_measure(pdf)

        while measured==0:
            print("Obtained zero, measuring again")
            measured = nonuniform_measure(pdf)
        
        print("Measured:",measured)
        fracs=list_convergents(measured,two_t)
        r=0

        for f in fracs:
            if f.denominator < N:
                r = f.denominator
                break
        
        print("r =",r)
            
        if r%2==1 and 2*r<N:
            r=2*r
        
        factor=-1
        if (x**r) % N!= N-1:
            if r%2==0:
                a=math.gcd(x**(r//2)+1,N)
                b=math.gcd(x**(r//2)+n-1,N)
                factor = max(a,b)
                                            
            if factor==1 or factor == N:
                print("trivial factor: %i, starting again" % factor)
                        
            elif factor!=-1:
                factor2 = N//factor
                print("%i = %i * %i" %(N,factor,factor2))
                print("----------------------------------")
                timetaken=int(time.time()-begin)
                if timetaken >= 60:
                    print("time taken:", timetaken//60, "min", timetaken%60, "sec" )
                else:
                    print("time taken:",timetaken,"sec")
                done=True

        else:
            print("r is odd, starting again")
    
    if tries==T:
        print("Tried %i times and failed!" % T)
        return
            

    


def main():
    print("SHOR'S ALGORITHM SIMULATION\n")

    N=int(input("Input number to be factorized\n"))
    
    
    print("Testing PowerPrimality...")
    if PowerPrimeQ(N):
        return

    Shor(N)


if __name__ == "__main__":
    main()
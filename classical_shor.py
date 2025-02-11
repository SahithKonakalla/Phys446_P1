import numpy as np

def isPrime(n):
    if n % 2 == 0 and n > 2: 
        return False
    return all(n % i for i in range(3, int(np.sqrt(n)) + 1, 2))

def isRoot(N):
    for a in range(2, int(np.ceil(np.log2(N)))+1):
        #print(np.float_power(N, 1/a))
        if np.float_power(N, 1/a).is_integer():
            return True
    return False

def classicalShor(N):
    if N % 2 == 0:
        return -1
    if isPrime(N):
        return -1
    if isRoot(N):
        return -1
    
    while True:
        x = np.random.randint(2, N)
        if np.gcd(x, N) != 1:
            return(np.gcd(x, N), N//np.gcd(x, N))
        
        r = 2
        while x**r % N != 1:
            r += 1

        if r % 2 == 1:
            continue

        A = np.gcd(int((x**(r//2)-1) % N),N)
        B = np.gcd(int((x**(r//2)+1) % N),N)
        if A != 1 and B != 1:
            return (A, B)
            #return (A, B, x, r)

for i in range(2,2**10):
    factors = classicalShor(i)
    if factors != -1:
        print(i, factors)

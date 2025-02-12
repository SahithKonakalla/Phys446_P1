import numpy as np
from fractions import Fraction
import sys
import phase_est_circuit
import quantum_simulators

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

def findPeriod(x, N):
    r = 2
    while x**r % N != 1:
        r += 1
    return r

def findPeriod2(x, N):
    eig = np.linalg.eigvals(buildUnitary(x,N))
    phases = (np.angle(eig) + 2*np.pi)/(2*np.pi)

    rand_eig = 0
    while rand_eig == 0:
        rand_eig = np.random.choice(phases)

    r = (Fraction(rand_eig).limit_denominator(39)).denominator
    return r

def findPeriodQ(x, N):
    u_wires = int(np.ceil(np.log2(N)))
    u_circuit = str(u_wires) + "\n xyModN 0 " + str(u_wires) + " " + str(x) + " " + str(N)
    top_wires = 2*u_wires + 1
    num_wires = top_wires + u_wires
    circuit = phase_est_circuit.makeArbitraryPhaseEstCircuit(top_wires,u_circuit)
    input_state = format(N, "0" + str(num_wires) + "b")

    quantum_simulators.writeCircuit(circuit, "Shor_Circuits/shor-"+ str(x) + "-" + str(N))
    out = quantum_simulators.Simulator_s(str(num_wires) + "\n" + "INITSTATE BASIS |" + input_state + "> \n" + circuit)
    #print(out)

    estim = 0
    amp = 0
    for i in out:
        temp_amp = np.conj(i[0])*i[0]
        if temp_amp > amp:
            amp = temp_amp
            estim = int(i[1][0:top_wires], 2)/(2**(top_wires))
    #print(estim)

    r = (Fraction(estim).limit_denominator(39)).denominator
    return r

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
        
        r = findPeriod(x,N)
        #print(x, N, r)

        if r % 2 == 1:
            continue

        A = np.gcd(int((x**(r//2)-1) % N),N)
        B = np.gcd(int((x**(r//2)+1) % N),N)
        if A != 1 and B != 1:
            #return (A, B)
            return (A, B, x, r)

def classicalShor2(N):
    if N % 2 == 0:
        return -1
    if isPrime(N):
        return -1
    if isRoot(N):
        return -1
    
    while True:
        x = np.random.randint(2, N)
        if np.gcd(x, N) != 1:
            continue
            #return(np.gcd(x, N), N//np.gcd(x, N))
        
        r = findPeriod2(x,N)
        #print(x, N, r)

        if r % 2 == 1:
            continue

        A = np.gcd(int((x**(r//2)-1) % N),N)
        B = np.gcd(int((x**(r//2)+1) % N),N)
        if A != 1 and B != 1:
            #return (A, B)
            return (A, B, x, r)

def QuantumShor(N):
    if N % 2 == 0:
        return -1
    if isPrime(N):
        return -1
    if isRoot(N):
        return -1
    
    while True:
        x = np.random.randint(2, N)
        if np.gcd(x, N) != 1:
            continue
            #return(np.gcd(x, N), N//np.gcd(x, N))
        
        print("Finding Period for: x =", x, ", N =", N)
        r = findPeriodQ(x,N)
        print(x, N, r)

        if r % 2 == 1:
            continue

        A = np.gcd(int((x**(r//2)-1) % N),N)
        B = np.gcd(int((x**(r//2)+1) % N),N)
        if A != 1 and B != 1:
            #return (A, B)
            return (A, B, x, r)

def buildUnitary(x, N):
    matrix_dim = 2**int(np.ceil(np.log2(N)))
    matrix = np.zeros((matrix_dim, matrix_dim))
    for i in range(matrix_dim):
        if i < N:
            matrix[i*x % N][i] = 1
        else:
            matrix[i][i] = 1
    return matrix

print(findPeriodQ(3,8))
print(findPeriodQ(3,5))
print(findPeriodQ(3,7))

# Quantum Shor

#QuantumShor(15)

# Classical with Unitary Matrix

'''for i in range(2,2**10):
    factors2 = classicalShor2(i)
    factors = classicalShor(i)
    if factors != -1:
        print(i, factors, " --- ", factors2)'''


# Eigenvalues

'''eig = np.linalg.eigvals(buildUnitary(2,5))
phases = (np.angle(eig) + 2*np.pi)/(2*np.pi)
print(phases)
r = 4 # 2^4 = 16 and 16 mod 5 = 1
print(phases*r)

eig = np.linalg.eigvals(buildUnitary(7,39))
phases = (np.angle(eig) + 2*np.pi)/(2*np.pi)
print(phases)
r = 12 # 2^4 = 16 and 16 mod 5 = 1
print(phases*r)'''

# Unitary

#print(buildUnitary(2,5))

# Classical

'''for i in range(2,2**10):
    factors = classicalShor(i)
    if factors != -1:
        print(i, factors)'''

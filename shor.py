import timeit
import numpy as np
from fractions import Fraction
import sys
import phase_est_circuit
import quantum_simulators
import os
import matplotlib.pyplot as plt

def buildUnitary(x, N):
    matrix_dim = 2**int(np.ceil(np.log2(N)))
    matrix = np.zeros((matrix_dim, matrix_dim))
    for i in range(matrix_dim):
        if i < N:
            matrix[i*x % N][i] = 1
        else:
            matrix[i][i] = 1
    return matrix

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

    r = (Fraction(rand_eig).limit_denominator(N)).denominator
    return r

def getTopWires(u_wires):
    return 2*u_wires - 1

def generateShorCircuit(x, N):
    u_wires = int(np.ceil(np.log2(N)))
    u_circuit = str(u_wires) + "\n xyModN 0 " + str(u_wires) + " " + str(x) + " " + str(N)
    top_wires = getTopWires(u_wires)
    num_wires = top_wires + u_wires
    circuit = phase_est_circuit.makeArbitraryPhaseEstCircuit(top_wires,u_circuit)
    input_state = format(1, "0" + str(num_wires) + "b")

    circuit = str(num_wires) + "\n" + "INITSTATE BASIS |" + input_state + "> \n" + circuit
    quantum_simulators.writeCircuit(circuit, "Shor_Circuits/shor-" + str(x) + "-" + str(N))
    return circuit

def generateShorFasterCircuit(x, N):
    u_wires = int(np.ceil(np.log2(N)))
    top_wires = getTopWires(u_wires)
    num_wires = top_wires + u_wires

    u_circuit = []
    
    #for i in range(top_wires):
    #    new_x = x
    #    for j in range(i):
    #        #print(i, j)
    #        u_circuit.append(str(u_wires) + "\n xyModN 0 " + str(u_wires) + " " + str(new_x) + " " + str(N))
    #        new_x = (new_x*new_x) % N
    new_x = x
    for i in range(top_wires):
        u_circuit.append(str(u_wires) + "\n xyModN 0 " + str(u_wires) + " " + str(new_x) + " " + str(N))
        new_x = (new_x*new_x) % N

    #print(u_circuit)

    circuit = phase_est_circuit.makeArbitraryShortPhaseEstCircuit(top_wires,u_circuit)
    input_state = format(1, "0" + str(num_wires) + "b")

    circuit = str(num_wires) + "\n" + "INITSTATE BASIS |" + input_state + "> \n" + circuit
    #print(circuit)
    quantum_simulators.writeCircuit(circuit, "Shor_Circuits/shorfast-" + str(x) + "-" + str(N))
    return circuit

#generateShorFasterCircuit(3, 5)

def findPeriodQ(x, N):
    u_wires = int(np.ceil(np.log2(N)))
    top_wires = getTopWires(u_wires)
    
    out = quantum_simulators.Simulator_s(generateShorCircuit(x, N))
    #print(out)

    estim = 0
    amp = 0
    for i in out:
        if int(i[1][0:top_wires], 2) == 0:
            continue
        temp_amp = np.conj(i[0])*i[0]
        if temp_amp > amp:
            amp = temp_amp
            estim = int(i[1][0:top_wires], 2)/(2**(top_wires))
    # print(estim)

    r = (Fraction(estim).limit_denominator(N)).denominator
    return r

def findPeriodQF(x, N):
    u_wires = int(np.ceil(np.log2(N)))
    top_wires = getTopWires(u_wires)
    
    out = quantum_simulators.Simulator_s(generateShorFasterCircuit(x, N) + "\n" + "MEASURE")
    #print(out)

    """ estim = 0
    amp = 0
    for i in out:
        #print(i[0], i[1], int(i[1][0:top_wires], 2)/(2**(top_wires)), int(i[1][0:top_wires], 2)/(2**(top_wires))*2, int(i[1][0:top_wires], 2)/(2**(top_wires))*4)
        if int(i[1][0:top_wires], 2) == 0:
            continue
        temp_amp = np.conj(i[0])*i[0]
        if temp_amp > amp:
            amp = temp_amp
            estim = int(i[1][0:top_wires], 2)/(2**(top_wires))
    print(estim) """
    estim = int(out[0:top_wires],2)/(2**(top_wires))
    print(estim)

    r = (Fraction(estim).limit_denominator(N)).denominator
    return r


def classicalShor(N):
    if disregardEasy(N):
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
        if A != 1 and B != 1 and A != N and B != N:
            return (A, B)
            #return (A, B, x, r)

def classicalShor2(N):
    if disregardEasy(N):
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
        if A != 1 and B != 1 and A != N and B != N:
            return (A, B)
            #return (A, B, x, r)

def QuantumShor(N):
    if disregardEasy(N):
        return -1
    
    while True:
        x = np.random.randint(2, N)
        if np.gcd(x, N) != 1:
            continue
        
        print("Finding Period for: x =", x, ", N =", N)
        r = findPeriodQ(x,N)
        print("Found Period of: r =", r)

        if r % 2 == 1:
            continue

        if x**r % N != 1:
            continue

        A = np.gcd(int((x**(r//2)-1) % N),N)
        B = np.gcd(int((x**(r//2)+1) % N),N)
        if A != 1 and B != 1 and A != N and B != N:
            return (A, B)
            #return (A, B, x, r)

def QuantumShorF(N):
    if disregardEasy(N):
        return -1
    
    while True:
        x = np.random.randint(2, N)
        if np.gcd(x, N) != 1:
            continue
        x = 7

        print("Finding Period for: x =", x, ", N =", N)
        r = findPeriodQF(x,N)
        print("Found Period of: r =", r)

        if r % 2 == 1:
            continue

        if x**r % N != 1:
            continue

        A = np.gcd(int((x**(r//2)-1) % N),N)
        B = np.gcd(int((x**(r//2)+1) % N),N)
        if A != 1 and B != 1 and A != N and B != N:
            return (A, B)
            #return (A, B, x, r)

def disregardEasy(N):
    if N % 2 == 0:
        return True
    if isPrime(N):
        return True
    if isRoot(N):
        return True


#print(findPeriodQ(3,8))
#print(findPeriodQ(3,5))
#print(findPeriodQ(3,10))
#print(findPeriodQ(14,15))
#print(findPeriodQ(5,21))

#print(QuantumShor(15))
#print(QuantumShorF(21))

# Gate Sizes

QuantumShorF(15)

'''for i in range(2,100):
    generateShorCircuit(1, i)

generated_circuits = os.listdir("Shor_Circuits")

n_list = []
gate_list = []
for circuit in generated_circuits:
    #print(circuit, circuit.split("-")[-1][:-8], len(open("Shor_Circuits/" + circuit).readlines()))
    n_list.append(int(circuit.split("-")[-1][:-8]))
    gate_list.append(len(open("Shor_Circuits/" + circuit).readlines()))

plt.figure(0)
plt.scatter(n_list,gate_list)
plt.xlabel("N")
plt.ylabel("Number of Gates")
plt.title("Gates needed for using Shor's Algorithm")
plt.savefig("Images/shorgates2")'''

# Speed Up times

""" start_time = timeit.default_timer()
print(QuantumShor(15))
elapsed = timeit.default_timer() - start_time
print("(SLOW) N=15 Time elapsed:", elapsed)

start_time = timeit.default_timer()
print(QuantumShorF(15))
elapsed = timeit.default_timer() - start_time
print("(FAST) N=15 Time elapsed:", elapsed) """

# Quantum Shor

#print(QuantumShor(21))

'''for i in range(7,25):
    print("Factoring:", i)
    factors3 = QuantumShor(i)
    factors2 = classicalShor2(i)
    factors = classicalShor(i)
    if factors != -1 or factors2 != -1 or factors3 != -1:
        print("Non-Trivial! Factors:")
        print(factors, "---", factors2, "---", factors3)'''

# CxyModN Gates

""" for i in range(0,7):
    for j in range(0,32,16):
        circuit = "5 \nINITSTATE BASIS |" + format(i+j, "0" + str(5) + "b") + "> \nCxyModN 0 1 4 2 8"

        print(circuit)
        print(quantum_simulators.Simulator_s(circuit))
 """
# xyModN Gates

""" for i in range(0,7):
    circuit = "4 \nINITSTATE BASIS |" + format(i, "0" + str(4) + "b") + "> \nxyModN 0 4 2 8"

    print(circuit)
    print(quantum_simulators.Simulator_s(circuit)) """


# Classical with Unitary Matrix

""" for i in range(2,2**10):
    factors2 = classicalShor2(i)
    #factors = classicalShor(i)
    if factors2 != -1:
        print(i, factors2) """

# Eigenvalues

""" for N in range(8,2**5):
    while True:
        x = np.random.randint(2, N)
        if np.gcd(x, N) != 1:
            continue
        
        r = findPeriod2(x,N)

        if x**r % N != 1:
            continue

        print("N =", N, "x =", x, "r =", r)
        break """

""" eig = np.linalg.eigvals(buildUnitary(2,5))
phases = (np.angle(eig) + 2*np.pi)/(2*np.pi)
print("x =",2,"N =",5)
print("e = ", phases)
r = 4 # 2^4 = 16 and 16 mod 5 = 1
print("r = ", r)
print("e*r = ", phases*r)

eig = np.linalg.eigvals(buildUnitary(7,39))
phases = (np.angle(eig) + 2*np.pi)/(2*np.pi)
print("x =",7,"N =",39)
print("e = ", phases)
r = 12 # 2^4 = 16 and 16 mod 5 = 1
print("r = ", r)
print("e*r = ", phases*r) """

# Unitary

#print(buildUnitary(2,5))

# Classical

""" for N in range(8,2**5):
    while True:
        x = np.random.randint(2, N)
        if np.gcd(x, N) != 1:
            continue
        
        r = findPeriod(x,N)

        if x**r % N != 1:
            continue

        print("N =", N, "x =", x, "r =", r)
        break """

#print(classicalShor(33))

""" for i in range(2,2**10):
    factors = classicalShor(i)
    if factors != -1:
        print(i, factors) """

print(QuantumShorF(15))
#print(findPeriodQ(13,15))
#print(findPeriodQF(14,15))
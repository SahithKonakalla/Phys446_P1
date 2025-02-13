import os
import numpy as np
import matplotlib.pyplot as plt
import timeit
import phase_est_circuit

def H(wire,inputState):
  newState = []

  for element in inputState:
    if int(element[1][wire]) == 0:
      first = element[1][:wire] + "0" + element[1][wire+1:]
      second = element[1][:wire] + "1" + element[1][wire+1:]
      newState.append(((1/np.sqrt(2))*element[0], first))
      newState.append(((1/np.sqrt(2))*element[0], second))
    if int(element[1][wire]) == 1:
      first = element[1][:wire] + "0" + element[1][wire+1:]
      second = element[1][:wire] + "1" + element[1][wire+1:]
      newState.append(((1/np.sqrt(2))*element[0], first))
      newState.append((-(1/np.sqrt(2))*element[0], second))
  return newState

def P(wire, theta, inputState):
    newState = []

    for element in inputState:
        amplitude, basis = element
        if int(basis[wire]) == 0:
            newState.append((amplitude, basis))
        elif int(basis[wire]) == 1:
            newAmplitude = amplitude * np.exp(1j * theta)
            newState.append((newAmplitude, basis))

    return newState

def CNOT(controlWire,notWire,inputState):
    newState = []
    for element in inputState:
        if int(element[1][controlWire]) == 0:
          newState.append((element[0], element[1]))
        elif int(element[1][controlWire]) == 1:
          changed = element[1][:notWire] + ("1" if element[1][notWire] == "0" else "0") + element[1][notWire+1:]
          newState.append((element[0], changed))
    return newState

def xyModN(firstWire, numberOfWires, x, N, inputState):
    newState = []

    for element in inputState:
        amplitude, basis = element
        int_basis = int(basis[firstWire:(firstWire+numberOfWires)], 2)
        if int_basis < N:
            #print(basis[:firstWire], format(int_basis*x % N, "0" + str(numberOfWires) + "b"), basis[(firstWire+numberOfWires):])
            newState.append((amplitude, basis[:firstWire] + format(int_basis*x % N, "0" + str(numberOfWires) + "b") + basis[(firstWire+numberOfWires):]))
        else:
            newState.append(element)
    return newState

def CxyModN(controlWire, firstWire, numberOfWires, x, N, inputState):
    newState = []

    for element in inputState:
        if int(element[1][controlWire]) == 0:
            newState.append(element)
            continue
        
        amplitude, basis = element
        int_basis = int(basis[firstWire:(firstWire+numberOfWires)], 2)
        if int_basis < N:
            newState.append((amplitude, basis[:firstWire] + format(int_basis*x % N, "0" + str(numberOfWires) + "b") + basis[(firstWire+numberOfWires):]))
        else:
            newState.append(element)
    return newState

def AddDuplicates(myState):
  num_bits = len(myState[0][1])
  vec = [0 for i in range(2**num_bits)]
  for state in myState:
    vec[int(state[1],2)] += np.round(state[0],4)
  newState = []
  for i in range(len(vec)):
    elem = np.round(vec[i],4)
    if elem != 0:
        newState.append((elem, format(i, "0" + str(num_bits) + "b")))
  return newState

def ReadInputString(myInput_lines):
  myInput=[]
  myInput_lines=myInput_lines.split('\n')
  myInput_lines = [ i for i in myInput_lines if i!='']
  numberOfWires=int(myInput_lines[0])
  for line in myInput_lines[1:]:
      myInput.append(line.split())
  return (numberOfWires,myInput)

# M-a
def HadamardArray(wire, num_wires):
    # Define the Hadamard gate matrix
    hadamard_matrix = np.array([[1/np.sqrt(2), 1/np.sqrt(2)],
                                  [1/np.sqrt(2), -1/np.sqrt(2)]])

    myMatrix = 1  # Start with a scalar identity
    for i in range(num_wires):
        if i == wire:
            myMatrix = np.kron(myMatrix, hadamard_matrix)  # Apply Hadamard on target wire
        else:
            myMatrix = np.kron(myMatrix, np.identity(2))  # Apply identity on other wires

    return myMatrix


def PhaseArray(wire, num_wires, theta):
    # Define the Hadamard gate matrix
    phase_matrix = np.array([[1, 0],
                            [0, np.exp(1j*theta)]])

    myMatrix = 1  # Start with a scalar identity
    for i in range(num_wires):
        if i == wire:
            myMatrix = np.kron(myMatrix, phase_matrix)  # Apply Hadamard on target wire
        else:
            myMatrix = np.kron(myMatrix, np.identity(2))  # Apply identity on other wires

    return myMatrix

def CNOTArray(control, output, num_wires):
    size = 2**num_wires
    cnot = []

    # Choose whether or not to use the flipped CNOT
    if control < output:
        cnot = np.matrix([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,0,1],
                        [0,0,1,0]])
    if control > output:
        cnot = np.matrix([[1,0,0,0],
                        [0,0,0,1],
                        [0,0,1,0],
                        [0,1,0,0]])
    
    myMatrix=1
    for i in range(num_wires):
        newMatrix = np.zeros((size, size))
        if control < output:
            if i == control:
                newMatrix = cnot
            elif i-1 == control:
                continue
            else:
                newMatrix = np.identity(2)
        else:
            if i == control:
                newMatrix = cnot
            elif i+1 == control:
                continue
            else:
                newMatrix = np.identity(2)
        
        myMatrix = np.kron(myMatrix, newMatrix)
    return myMatrix

def ReadInputString(myInput_lines):
    myInput=[]
    myInput_lines=myInput_lines.split('\n')
    myInput_lines = [ i for i in myInput_lines if i!='']
    numberOfWires=int(myInput_lines[0])
    measure = False
    myState = [1 if i == 0 else 0 for i in range(2**numberOfWires)]
    for line in myInput_lines[1:]:
        if len(line.split()) == 0:
            continue
        if line == "MEASURE":
            measure = True
        elif line.split()[0] == "INITSTATE":
            if line.split()[1] == "FILE":
                myState = ReadInputStateFromFile(line.split()[2])
            elif line.split()[1] == "BASIS":
                myState = ReadInputStateFromBasis(line.split()[2])
        else:
            myInput.append(line.split())
    return (numberOfWires,myInput,myState, measure)

def ReadInputStateFromFile(file):
    input_lines = open(file).read()
    input_lines = input_lines.split('\n')
    input_lines = [ i for i in input_lines if i!='']

    state = []
    for line in input_lines:
        real = float(line.split()[0])
        imag = float(line.split()[1])
        state.append(real + 1j*imag)

    return state

def MakeFileFromBasis(basis, filename):
    f = open(filename, "w")
    for elem in basis:
        f.write(str(np.real(elem)) + " " + str(np.imag(elem)) + "\n")

def ReadInputStateFromBasis(basis):
    index = int(basis[1:(len(basis)-1)], 2)

    myState = [1 if i == index else 0 for i in range(2**(len(basis)-2))]

    return myState

def measureState(myState):

    myStateConj = np.conjugate(myState)
    myState = np.real(np.multiply(myStateConj, myState))
    return format(np.where(np.random.multinomial(1, myState) == 1)[0][0], "0" + str(int(np.log2(len(myState)))) + "b")

def Simulator_s(circuit):
    numberOfWires,myInput,myState,measure=ReadInputString(circuit)

    myInput = precompile(myInput)

    myState = VecToDirac(myState)

    circuit_len = len(myInput)
    print(circuit_len, " gates... ")
    #print(sum(1 for n in myInput if n[0] != "CxyModN"))

    i = 0
    for gate in myInput:
        match gate[0]:
            case "H":
                myState = H(int(gate[1]), myState)
            case "P":
                myState = P(int(gate[1]), float(gate[2]), myState)
            case "CNOT":
                myState = CNOT(int(gate[1]), int(gate[2]), myState)
            case "xyModN":
                myState = xyModN(int(gate[1]), int(gate[2]), int(gate[3]), int(gate[4]), myState)
            case "CxyModN":
                myState = CxyModN(int(gate[1]), int(gate[2]), int(gate[3]), int(gate[4]), int(gate[5]), myState)
    
        i += 1
        #print(i, circuit_len, gate)
        printProgressBar(i, circuit_len, prefix = 'Progress:', suffix = 'Complete', length = 50)
        myState = AddDuplicates(myState)

    if measure:
        return measureState(DiracToVec(myState))
    return myState

def Simulator_a(circuit, print_matrix = False):
    numberOfWires,myInput,myState,measure=ReadInputString(circuit)

    myInput = precompile(myInput)
    
    myMatrix = np.identity(2**numberOfWires)

    circuit_len = len(myInput)
    i = 0
    for gate in myInput:
        match gate[0]:
            case "H":
                myMatrix = HadamardArray(int(gate[1]), numberOfWires).dot(myMatrix)
            case "P":
                myMatrix = PhaseArray(int(gate[1]), numberOfWires, float(gate[2])).dot(myMatrix)
            case "CNOT":
                myMatrix = CNOTArray(int(gate[1]), int(gate[2]), numberOfWires).dot(myMatrix)
        i += 1
        printProgressBar(i, circuit_len, prefix = 'Progress:', suffix = 'Complete', length = 50)
    if print_matrix:
        if measure:
            return measureState(np.ravel(np.asarray(myMatrix.dot(myState)))), np.round(myMatrix, 3)
        return VecToDirac(np.round(np.ravel(np.asarray(myMatrix.dot(myState))), 4)), np.round(myMatrix,3)
    else:
        if measure:
            return measureState(np.ravel(np.asarray(myMatrix.dot(myState))))
        return VecToDirac(np.round(np.ravel(np.asarray(myMatrix.dot(myState))), 4))

def Simulator_b(circuit):
    numberOfWires,myInput,myState,measure=ReadInputString(circuit)

    myInput = precompile(myInput)
    
    myState = np.transpose([myState])

    circuit_len = len(myInput)
    i = 0
    for gate in myInput:
        match gate[0]:
            case "H":
                myState = HadamardArray(int(gate[1]), numberOfWires).dot(myState)
            case "P":
                myState = PhaseArray(int(gate[1]), numberOfWires, float(gate[2])).dot(myState)
            case "CNOT":
                myState = CNOTArray(int(gate[1]), int(gate[2]), numberOfWires).dot(myState)
        i += 1
        printProgressBar(i, circuit_len, prefix = 'Progress:', suffix = 'Complete', length = 50)
        return measureState(np.ravel(myState))
    return VecToDirac(np.round(np.ravel(np.transpose(myState)), 4))

def precompile(circuit):
    for i in range(len(circuit)):
        for j in range(len(circuit[i])):
            if "pi" in circuit[i][j]:
                circuit[i][j] = eval(circuit[i][j])

    depth = 1
    j = 0
    while j < depth:
        i = 0
        while i < len(circuit):
            match circuit[i][0]:
                case "CNOT":
                    control = circuit[i][1]
                    output = circuit[i][2]

                    if abs(int(control) - int(output)) > 1: # Long Range
                        tempgate = "0"
                        if int(control) > int(output):
                            tempgate = str(int(output)+1)
                        else:
                            tempgate = str(int(control)+1)
                        
                        circuit[i] = ["SWAP", tempgate, control]
                        circuit.insert(i+1, ["CNOT", tempgate, output])
                        circuit.insert(i+2, ["SWAP", tempgate, control])
                        depth += 1

                        #circuit[i] = ["CNOT", tempgate, control] # Swap 1
                        #circuit.insert(i+1, ["CNOT", control, tempgate])
                        #circuit.insert(i+2, ["CNOT", tempgate, control])

                        #circuit.insert(i+3, ["CNOT", tempgate, output]) # CNOT

                        #circuit.insert(i+4, ["CNOT", tempgate, control]) #Swap 2
                        #circuit.insert(i+5, ["CNOT", control, tempgate])
                        #circuit.insert(i+6, ["CNOT", tempgate, control])
                case "NOT":
                    gate = circuit[i][1]
                    circuit[i] = ["H", gate]
                    circuit.insert(i+1, ["P", gate, str(np.pi)])
                    circuit.insert(i+2, ["H", gate])
                case "RZ":
                    gate = circuit[i][1]
                    phase = float(circuit[i][2])
                    circuit[i] = ["NOT", circuit[i][1]]
                    circuit.insert(i+1, ["P", gate, str(-phase/2)])
                    circuit.insert(i+2, ["NOT", gate])
                    circuit.insert(i+3, ["P", gate, str(phase/2)])
                    depth += 1
                case "CRZ":
                    control = circuit[i][1]
                    output = circuit[i][2]
                    phase = float(circuit[i][3])
                    circuit[i] = ["CNOT", control, output]
                    circuit.insert(i+1, ["P", output, str(-phase/2)])
                    circuit.insert(i+2, ["CNOT", control, output])
                    circuit.insert(i+3, ["P", output, str(phase/2)])
                    depth += 1
                case "CPHASE":
                    control = circuit[i][1]
                    output = circuit[i][2]
                    phase = float(circuit[i][3])
                    circuit[i] = ["CRZ", control, output, str(phase)]
                    circuit.insert(i+1, ["P", control, str(phase/2)])  
                    depth += 1
                case "SWAP":
                    gate1 = circuit[i][1]
                    gate2 = circuit[i][2]

                    count = 0
                    if abs(int(gate1) - int(gate2)) > 1:
                        for g in range(int(gate1), int(gate2)):
                            if count == 0:
                                circuit[i] = ["SWAP", str(g), str(g+1)]
                            else:
                                circuit.insert(i+count, ["SWAP", str(g), str(g+1)])
                            count += 1
                        for g in range(int(gate2)-1, int(gate1), -1):
                            circuit.insert(i+count, ["SWAP", str(g), str(g-1)])
                            count += 1
                        depth += 1
                    else:
                        circuit[i] = ["CNOT", gate1, gate2]
                        circuit.insert(i+1, ["CNOT", gate2, gate1])
                        circuit.insert(i+2, ["CNOT", gate1, gate2])
                        
            i += 1
        j += 1
    
    return circuit


def DiracToVec(myState):
    num_bits = len(myState[0][1])
    vec = [0 for i in range(2**num_bits)]
    for state in myState:
        index = int(state[1],2)
        vec[index] += np.round(state[0],3)
    return np.array(vec)

def VecToDirac(myVec):
    state = []
    num_bits = int(np.log2(len(myVec)))
    for i in range(len(myVec)):
        elem = myVec[i]
        if elem != 0:
            state.append((elem, format(i, "0" + str(num_bits) + "b")))
    return state

def isDirac(myState):
    if isinstance(myState[0], (int,float,complex)):
        return False
    return True

def writeCircuit(circuit, filename):
    f = open(filename + ".circuit", "w")
    f.write(circuit)
    f.close()

def writeCircuitArray(circuit, filename):
    f = open(filename + ".circuit", "w")

    output = ""
    for gate in circuit:
        output += " ".join(gate) + "\n"
    f.write(output)
    f.close()

def phaseEstBoxer(top_wires,num_wires,out):
    vals = [0 for i in range(2**top_wires)]

    for state in out:
        pos = int(state[1][:top_wires], 2)
        vals[pos] += np.conj(state[0])*state[0]
    
    return vals

def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

np.set_printoptions(formatter={'all': lambda x: "{:.4g}".format(x)})

circuit = '''
5
SWAP 2 5
'''

numberOfWires,myInput,myState,measure=ReadInputString(circuit)
writeCircuitArray(precompile(myInput), "Circuits/swap25")


#circuit = open('Circuits/example.circuit').read()

#print(Simulator_a(circuit,True)[0])

# Shor's

# Gates vs N

# CxyModN


'''x = 2
N = 5
mat = 3

for q in range(mat+1, 8):
    for w in range(1,q-mat+1):
        print("qubits: ", q, ", first wire: ", w)
        for i in range(0,2**q):
            circuit = str(q) + " \n INITSTATE BASIS |" + format(i, "0" + str(q) + "b") + "> \n CxyModN 0 " + str(w) + " " + str(mat) + " " + str(x) + " " + str(N) 
            amplitude, basis = Simulator_s(circuit)[0]
            int_basis = int(format(i, "0" + str(q) + "b")[w:(w+mat)],2)
            if int(basis[0]) == 1 and int(basis[w:(w+mat)], 2) == (int_basis*x % N if int_basis < N else int_basis):
                print(format(i, "0" + str(q) + "b"), "c", int(basis[w:(w+mat)], 2), "==", (int_basis*x % N if int_basis < N else int_basis))
            elif int(basis[0]) == 0 and int(basis[w:(w+mat)], 2) == int_basis:
                print(format(i, "0" + str(q) + "b"), "nc", int(basis[w:(w+mat)], 2) , "==", int_basis)'''


# xyModN

'''
x = 2
N = 5
mat = 3

for q in range(mat, 8):
    for w in range(q-mat+1):
        print("qubits: ", q, ", first wire: ", w)
        for i in range(0,2**q):
            circuit = str(q) + " \n INITSTATE BASIS |" + format(i, "0" + str(q) + "b") + "> \n xyModN " + str(w) + " " + str(mat) + " " + str(x) + " " + str(N) 
            amplitude, basis = Simulator_s(circuit)[0]
            int_basis = int(format(i, "0" + str(q) + "b")[w:(w+mat)],2)
            print(i, int(basis[w:(w+mat)], 2) == (int_basis*x % N if int_basis < N else int_basis), int(basis[w:(w+mat)], 2) , " == ", (int_basis*x % N if int_basis < N else int_basis))
'''

# Quantum Fourier Transform

# 4 Qubits

#print(phase_est_circuit.makeQFTCircuit(4))

# 5 Qubits

'''# qftcircuit = phase_est_circuit.makeQFTCircuit(5)
# writeCircuit("5 \n" + "INITSTATE FILE myInputState.txt \n"+ qftcircuit, "Circuits\qft5")
circuit = open('Circuits/qft5.circuit').read()
print(Simulator_s(circuit))'''

# 3 Qubits

'''circuit = open("Circuits/qft3.circuit").read()

A = Simulator_a(circuit, True)[1]

B = np.zeros((8,8), dtype=complex)
for i in range(len(B)):
    for j in range(len(B[0])):
        B[i][j] = np.round(np.power(np.sqrt(0+1j),(i*j)), 2)/np.sqrt(8)

print(np.round(A-B,2))'''

# Phase Estimation

# Arbitrary 

'''phase_estimation_circuit = "8 \n INITSTATE BASIS |00000011> \n" + phase_est_circuit.makeArbitraryPhaseEstCircuit(6, "NOT 0 \n P 0 0.3 \n NOT 1 \n")
out = Simulator_s(phase_estimation_circuit)
print(out)

steps = [i/(2**(len(out[0][1])-2)) for i in range(2**6+1)]
vals = phaseEstBoxer(6,8,out)

plt.figure(0)
plt.stairs(vals, steps, fill=True)
plt.title("Phase Estimation with New Circuit")
plt.xlabel("Estimated Phase")
plt.ylabel("Probability")
plt.savefig("Images/phaseesthistarbit.png")
plt.show()'''

# Speed

'''phase_estimation_circuit = "7 \n INITSTATE BASIS |0000001> \n" + phase_est_circuit.makePhaseEstCircuit(6,7, 0.5*2*np.pi)
out = Simulator_s(phase_estimation_circuit)
print(out)

steps = [i/(2**(len(out[0][1])-1)) for i in range(2**6+1)]
vals = phaseEstBoxer(6,7,out)

#print(vals)
#print(steps)

plt.figure(0)
plt.stairs(vals, steps, fill=True)
plt.axvline(x=0.5, color='r', linestyle='--')
plt.title("Phase Estimation of 0.5 with New Phase Trick")
plt.xlabel("Estimated Phase")
plt.ylabel("Probability")
plt.savefig("Images/phaseesthistspeed.png")
plt.show()'''

# Eigenstates

'''testvec = [0 for i in range(2**7)]
testvec[0] = np.sqrt(0.3)
testvec[1] = np.sqrt(0.7)
MakeFileFromBasis(testvec, "eigenstates_test.txt")
phase_estimation_circuit = "7 \n INITSTATE FILE eigenstates_test.txt \n" + phase_est_circuit.makePhaseEstCircuit(6,7, 0.5*2*np.pi)
out = Simulator_s(phase_estimation_circuit)
print(out)

steps = [i/(2**(len(out[0][1])-1)) for i in range(2**6+1)]
vals = phaseEstBoxer(6,7,out)

#print(vals)
#print(steps)

plt.figure(0)
plt.stairs(vals, steps, fill=True)
plt.axvline(x=0.5, color='r', linestyle='--')
plt.title("Phase Estimation of 0.5 with Different Eigenstate")
plt.xlabel("Estimated Phase")
plt.ylabel("Probability")
plt.savefig("Images/phaseesthisteigen.png")
plt.show()'''

# 6 Wires

'''phase = np.linspace(0, 2*np.pi, 9)
estim_list = []

for p in phase:
    print(p/(2*np.pi))
    phase_estimation_circuit = "7 \n INITSTATE BASIS |0000001> " + phase_est_circuit.makePhaseEstCircuit(6,7, p)

    out = Simulator_s(phase_estimation_circuit)
    estim = 0
    amp = 0

    for i in out:
        temp_amp = np.conj(i[0])*i[0]
        if temp_amp > amp:
            amp = temp_amp
            estim = int(i[1][0:len(i[1])-1], 2)/(2**(len(i[1])-1))

    estim_list.append(estim)

plt.figure(0)
plt.plot(phase/(2*np.pi), estim_list)
plt.title("Phase Estimation with Six Wires")
plt.xlabel("Real Phase")
plt.ylabel("Estimated Phase")
plt.savefig("Images/phaseest6.png")

p = 0.1432394487827058*2*np.pi
phase_estimation_circuit = "7 \n INITSTATE BASIS |0000001> \n" + phase_est_circuit.makePhaseEstCircuit(6,7, p)
out = Simulator_s(phase_estimation_circuit)

vals = phaseEstBoxer(6,7,out)
steps = [i/(2**(len(out[0][1])-1)) for i in range(len(out)+1)]

#print(vals)
#print(steps)

plt.figure(1)
plt.stairs(vals, steps, fill=True)
plt.axvline(x=0.1432394487827058, color='r', linestyle='--')
plt.title("Phase Estimation of 0.1432394487827058")
plt.xlabel("Estimated Phase")
plt.ylabel("Probability")
plt.savefig("Images/phaseesthist6.png")
plt.show()'''

# Two Wires

'''phase = np.linspace(0, 2*np.pi, 100)
estim_list = []

for p in phase:
    phase_estimation_circuit = "3 \n INITSTATE BASIS |001> \n H 0 \n H 1 \n CPHASE 1 2 " + str(p) + " \n CPHASE 0 2 " + str(p) + " \n CPHASE 0 2 " + str(p) + "\n H 0 \n CPHASE 0 1 " + str(-np.pi/2) + " \n H 1 \n SWAP 0 1"

    out = Simulator_s(phase_estimation_circuit)

    estim = 0
    amp = 0

    for i in out:
        temp_amp = np.conj(i[0])*i[0]
        if temp_amp > amp:
            amp = temp_amp
            estim = int(i[1][0:2], 2)/(2**(len(i[1])-1))

    estim_list.append(estim)

plt.figure(0)
plt.plot(phase/(2*np.pi), estim_list)
plt.title("Phase Estimation with Two Wires")
plt.xlabel("Real Phase")
plt.ylabel("Estimated Phase")
plt.savefig("Images/phaseest2.png")

p = 0.1432394487827058*2*np.pi
phase_estimation_circuit = "3 \n INITSTATE BASIS |001> \n H 0 \n H 1 \n CPHASE 1 2 " + str(p) + " \n CPHASE 0 2 " + str(p) + " \n CPHASE 0 2 " + str(p) + "\n H 0 \n CPHASE 0 1 " + str(-np.pi/2) + " \n H 1 \n SWAP 0 1"
out = Simulator_s(phase_estimation_circuit)

vals = phaseEstBoxer(2,3,out)
steps = [i/(2**(len(out[0][1])-1)) for i in range(len(out)+1)]


plt.figure(1)
plt.stairs(vals, steps, fill=True)
plt.axvline(x=0.1432394487827058, color='r', linestyle='--')
plt.title("Phase Estimation of 0.1432394487827058")
plt.xlabel("Estimated Phase")
plt.ylabel("Probability")
plt.savefig("Images/phaseesthist2.png")
plt.show()'''

# One Wire

'''phase = np.linspace(0, 2*np.pi, 100)
estim_list = []

for p in phase:
    phase_estimation_circuit = "2 \n INITSTATE BASIS |01> \n H 0 \n CPHASE 0 1 " + str(p) + "\n H 0"
    out = Simulator_s(phase_estimation_circuit)
    print(out)
    
    estim = 0
    amp = 0

    for i in out:
        temp_amp = np.conj(i[0])*i[0]
        if temp_amp > amp:
            amp = temp_amp
            estim = int(i[1][0])/2
    
    estim_list.append(estim)

plt.figure(0)
plt.plot(phase/(2*np.pi), estim_list)
plt.title("Phase Estimation with One Wire")
plt.xlabel("Real Phase")
plt.ylabel("Estimated Phase")
plt.savefig("phaseest1.png")'''

# Histogram

'''phase_estimation_circuit = "2 \n INITSTATE BASIS |01> \n H 0 \n CPHASE 0 1 " + str(2*np.pi*0.1432394487827058) + "\n H 0"
out = Simulator_s(phase_estimation_circuit)

vals = phaseEstBoxer(1,2,out)
steps = [i/len(out[0][1]) for i in range(len(out)+1)]

print(vals)
print(steps)

plt.figure(1)
plt.stairs(vals, steps, fill=True)
plt.axvline(x=0.1432394487827058, color='r', linestyle='--')
plt.title("Phase Estimation of 0.1432394487827058")
plt.xlabel("Estimated Phase")
plt.ylabel("Probability")
plt.savefig("phaseesthist1.png")
plt.show()'''

# Measure Time Complexity

'''skip_m = 12
s_times = []
a_times = []
b_times = []

max_num_bits = 20
num_gates = 10
for i in range(2,max_num_bits+1):
    circuit = open("Timed_Circuits/circuit-" + str(i) + "-" + str(num_gates) + ".circuit").read()

    print(i, "qubits")

    start_time = timeit.default_timer()
    out1 = Simulator_s(circuit)
    elapsed = timeit.default_timer() - start_time
    print("s-s:", elapsed)
    s_times.append(elapsed)

    if i >= skip_m:
        continue

    start_time = timeit.default_timer()
    out2 = Simulator_a(circuit)
    elapsed = timeit.default_timer() - start_time
    print("s-a:", elapsed)
    a_times.append(elapsed)

    start_time = timeit.default_timer()
    out3 = Simulator_b(circuit)
    elapsed = timeit.default_timer() - start_time
    print("s-b:", elapsed)
    b_times.append(elapsed)

plt.figure(0)
plt.plot([i for i in range(2, max_num_bits + 1)], s_times)
plt.plot([i for i in range(2, skip_m)], a_times)
plt.plot([i for i in range(2, skip_m)], b_times)

plt.legend(["Simulator S", "Simulator M-a", "Simulator M-b"])
plt.title("Time Complexity of Different Quantum Simulators")
plt.xlabel("Number of Qubits")
plt.ylabel("Time (s)")
plt.savefig("timecomp2.png")'''

# Measure Probability Distribution

'''outputState = DiracToVec(Simulator_a(measure_circuit2, myState))

outputStateConj = np.conjugate(outputState)
outputStateSquared = np.real(np.multiply(outputStateConj, outputState))

print(outputStateSquared)
#plt.stairs(outputStateSquared, [i for i in range(32)])
bins = [i for i in range(32)]

plt.figure(1)
plt.title("Expected Probability Distribution")
plt.ylabel("Probabilty")
plt.xlabel("State")
plt.hist(bins, bins, weights=outputStateSquared)
plt.savefig("ExpPD")

prob_dist = np.zeros(32)
for i in range(1000):
    prob_dist[Simulator_a(measure_circuit, myState)] += 1

plt.figure(2)
plt.title("Measured Probability Distribution")
plt.ylabel("Counts")
plt.xlabel("State")
plt.hist(bins, bins, weights=prob_dist)
plt.savefig("MeasPD")

#print(Simulator_b(measure_circuit, myState))'''

# Measure

'''circuit = open("Circuits/input.circuit").read()
print(Simulator_s(circuit))'''
import numpy as np
import matplotlib.pyplot as plt


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
        myInput.append(line.split())
        if line == "MEASURE":
            measure = True
        elif line.split()[0] == "INITSTATE":
            if line.split()[1] == "FILE":
                myState = ReadInputStateFromFile(line.split()[2])
            elif line.split()[1] == "BASIS":
                myState = ReadInputStateFromBasis(line.split()[2])
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

def ReadInputStateFromBasis(basis):
    index = int(basis[1:(len(basis)-1)], 2)

    myState = [1 if i == index else 0 for i in range(2**(len(basis)-2))]

    return myState

def measureState(myState):

    myStateConj = np.conjugate(myState)
    myState = np.real(np.multiply(myStateConj, myState))

    return np.random.choice(len(myState), 1, p=myState)

def Simulator_a(circuit):
    numberOfWires,myInput,myState,measure=ReadInputString(circuit)
    
    myMatrix = np.identity(2**numberOfWires)
    for gate in myInput:
        match gate[0]:
            case "H":
                myMatrix = HadamardArray(int(gate[1]), numberOfWires).dot(myMatrix)
            case "P":
                myMatrix = PhaseArray(int(gate[1]), numberOfWires, float(gate[2])).dot(myMatrix)
            case "CNOT":
                myMatrix = CNOTArray(int(gate[1]), int(gate[2]), numberOfWires).dot(myMatrix)
    if measure:
        return measureState(np.ravel(np.asarray(myMatrix.dot(myState))))
    return VecToDirac(np.ravel(np.asarray(myMatrix.dot(myState))))

def Simulator_b(circuit):
    numberOfWires,myInput,myState,measure=ReadInputString(circuit)
    
    myState = np.transpose([myState])
    for gate in myInput:
        match gate[0]:
            case "H":
                myState = HadamardArray(int(gate[1]), numberOfWires).dot(myState)
            case "P":
                myState = PhaseArray(int(gate[1]), numberOfWires, float(gate[2])).dot(myState)
            case "CNOT":
                myState = CNOTArray(int(gate[1]), int(gate[2]), numberOfWires).dot(myState)
    if measure:
        return measureState(np.ravel(myState))
    return VecToDirac(np.ravel(np.transpose(myState)))

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

circuit = '''
3
H 1
H 2
P 2 0.3
CNOT 2 1
H 1
H 2
MEASURE
'''

measure_circuit = open('measure.circuit').read()
measure_circuit2 = open('measure_full.circuit').read()
input_circuit = open('input.circuit').read()
input_circuit2 = open('input2.circuit').read()
not_test_circuit = open('not_test.circuit').read()


np.set_printoptions(formatter={'all': lambda x: "{:.4g}".format(x)})

print(Simulator_a(not_test_circuit))
#print(Simulator_b(circuit, myState))

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
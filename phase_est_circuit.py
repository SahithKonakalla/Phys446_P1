import numpy as np

def makeQFTCircuit(num_wires):
    output = ""

    for i in range(num_wires-1):
        for j in range(num_wires-i-1):
            output += "SWAP " + str(j) + " " + str(j+1) + "\n"

    
    for i in range(num_wires):
        output += "H " + str(num_wires-i-1) + "\n"
        for j in range(num_wires-i-1):
            output += "CPHASE " + str(num_wires-i-1) + " " + str(num_wires-j-i-2) + " " + str(-np.pi/(2**(j+1))) + "\n"
    return output

def invert(circuit):
    lines = circuit.split('\n')
    reversed_lines = lines[::-1]
    reversed_text = '\n'.join(reversed_lines)
    
    return reversed_text

def makePhaseEstCircuit(top_wires, num_wires, phase):
    #f = open("Phase_Circuits/phase_est_circuit-" + str(num_wires) + "-" + str(phase) +  ".circuit", "w")
    output = str(num_wires) + "\n"

    # Hadamard Section
    for i in range(top_wires):
        output += "H " + str(i) + "\n"

    # U Section
    for i in range(top_wires):
        for j in range(2**i):
            output += "CPHASE " + str(top_wires-i-1) + " " + str(top_wires) + " " + str(phase) + "\n"
    
    # Inverse QFT
    output += invert(makeQFTCircuit(top_wires))
    return output

print(makePhaseEstCircuit(2,3, np.pi))
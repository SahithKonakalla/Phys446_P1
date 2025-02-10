

def makeCircuit(top_wires, num_wires):
    f = open("Phase_Circuits/phase_est_circuit-" + str(num_wires) + ".circuit", "w")
    f.write(str(num_wires) + "\n")

    # Hadamard Section
    for i in range(top_wires):
        f.write("H " + str(i) + "\n")

    # U Section
    for i in range(top_wires):
        for j in range(2**i):
            f.write("CPHASE " + str(top_wires-i-1) +" "+ str(top_wires)+ "\n")
    
    # Inverse QFT
    f.write("Inverse QFT " + "\n")

makeCircuit(6,7)
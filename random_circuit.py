import numpy as np

gate_list = ["H", "P", "CNOT"] #, "NOT", "RZ", "CRZ", "CPHASE", "SWAP"

def makeCircuit(num_wires, num_gates):

    rand_list = np.random.randint(0, len(gate_list), num_gates)

    f = open("Timed_Circuits/circuit-" + str(num_wires) + "-" + str(num_gates) + ".circuit", "w")
    f.write(str(num_wires) + "\n")

    gates = []
    for rand in rand_list:
        gates.append(gate_list[rand])

    for gate in gates:
        args = [gate]
        match gate:
            case "H":
                args.append(np.random.randint(0, num_wires))
            case "P":
                args.append(np.random.randint(0, num_wires))
                args.append(np.random.rand()*2*np.pi)
            case "CNOT":
                arg1 = np.random.randint(0, num_wires)
                args.append(arg1)
                args.append(np.random.choice([i for i in range(0,num_wires) if i not in [arg1]]))

        for arg in args:
            f.write(str(arg) + " ")
        
        f.write("\n")

for i in range(2, 21):
    makeCircuit(i, 10)

'''case "NOT":
    args.append(np.random.randint(0, num_wires))
case "RZ":
    args.append(np.random.randint(0, num_wires))
    args.append(np.random.rand()*2*np.pi)
case "CRZ":
    arg1 = np.random.randint(0, num_wires)
    args.append(arg1)
    args.append(np.random.choice([i for i in range(0,num_wires) if i not in [arg1]]))
    args.append(np.random.rand()*2*np.pi)
case "CPHASE":
    arg1 = np.random.randint(0, num_wires)
    args.append(arg1)
    args.append(np.random.choice([i for i in range(0,num_wires) if i not in [arg1]]))
    args.append(np.random.rand()*2*np.pi)
case "SWAP":
    arg1 = np.random.randint(0, num_wires)
    args.append(arg1)
    args.append(np.random.choice([i for i in range(0,num_wires) if i not in [arg1]]))'''

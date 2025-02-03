import numpy as np


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

    # Iterate through each element in the input state
    for element in inputState:
        amplitude, basis = element  # Unpack amplitude and basis state
        if int(basis[wire]) == 0:
            # If the qubit is in state |0>, no phase change
            newState.append((amplitude, basis))
        elif int(basis[wire]) == 1:
            # If the qubit is in state |1>, apply phase shift e^(i * theta)
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

def AddDuplicates(myState):
  num_bits = len(myState[0][1])
  vec = [0 for i in range(2**num_bits)]
  for state in myState:
    vec[int(state[1],2)] += np.round(state[0],3)
  newState = []
  for i in range(len(vec)):
    elem = vec[i]
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

def Simulator(myInput, myState):
  for gate in myInput:
    match gate[0]:
      case "H":
        myState = H(int(gate[1]), myState)
      case "P":
        myState = P(int(gate[1]), float(gate[2]), myState)
      case "CNOT":
        myState = CNOT(int(gate[1]), int(gate[2]), myState)
        
  
  return AddDuplicates(myState)

#numberOfWires,myInput=ReadInputString(open('example.circuit').read())

numberOfWires,myInput=ReadInputString('''
3
H 1
H 2
P 2 0.3
CNOT 2 1
H 1
H 2

''')#CNOT 2 0

myState=[
  (1, '000')
]

np.set_printoptions(formatter={'all': lambda x: "{:.4g}".format(x)})

print(Simulator(myInput, myState))
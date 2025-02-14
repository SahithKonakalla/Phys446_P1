import numpy as np

def PrettyPrintBinary(myState):
    print("( ", end="")
    i = 0
    for state in myState:
        print(state[0], "|" + str(state[1]) + ">", end="")
        if i < len(myState) - 1:
            print(" + ", end="")
        i += 1
    print(" )")

def PrettyPrintInteger(myState):
    print("( ", end="")
    i = 0
    for state in myState:
        print(state[0], "|" + str(int(str(state[1]),2)) + ">", end="")
        if i < len(myState) - 1:
            print(" + ", end="")
        i += 1
    print(" )")


def DiracToVec(myState):
    num_bits = len(myState[0][1])
    vec = [0 for i in range(2**num_bits)]
    for state in myState:
        index = int(state[1],2)
        vec[index] += np.round(state[0],3)
    return vec

def VecToDirac(myVec):
    state = []
    num_bits = int(np.log2(len(myVec)))
    for i in range(len(myVec)):
        elem = myVec[i]
        if elem != 0:
            state.append((elem, format(i, "0" + str(num_bits) + "b")))
    return state

myState2=[
   (np.sqrt(0.1)*1.j, '101'),
   (np.sqrt(0.5), '000') ,
   (-np.sqrt(0.4), '010' )] 

#PrettyPrintBinary(myState2)
#PrettyPrintInteger(myState2)

#print(DiracToVec(myState2))
#print(VecToDirac(DiracToVec(myState2)))
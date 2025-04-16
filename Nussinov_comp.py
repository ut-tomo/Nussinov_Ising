import numpy as np
from datetime import timedelta
import amplify as amp

n = 8
bases = np.random.choice(['A', 'U', 'G', 'C'], size=n) #randomly generate RNA sequence
print("RNA sequence:", ''.join(bases))

#Possible Pairs
def IsValidPair(bi, bj):
    return (bi == 'A' and bj == 'U') or \
           (bi == 'U' and bj == 'A') or \
           (bi == 'G' and bj == 'C') or \
           (bi == 'C' and bj == 'G') or \
           (bi == 'G' and bj == 'U') or \
           (bi == 'U' and bj == 'G')
ValidPairs = [(i, j) for i in range(n) for j in range(n) if IsValidPair(bases[i], bases[j])]

#Generate binary variables
gen = amp.VariableGenerator()
x = gen.array("Binary", shape = (n, n))
for i in range(n):
    for j in range(n):
        if (IsValidPair(bases[i], bases[j]) == False):
            x[i, j] = 0
        if (i > j):
            x[i, j] = x[j, i]

            
score = sum(
    x[i, j] 
    for i in range(n) for j in range(i + 1, n)
) 
obj = - score

"""
CONSTRAINTS  
Constraint 1:  
A single base can pair with only one other base >> The sum of each row and column must be at most 1  

Constraint 2:  
No crossing pairs allowed  

Constraint 3:  
Bases that are too close (less than 3 bases apart) cannot form a pair  

"""

row_constraints = amp.sum(amp.less_equal(amp.sum(x[i, :]), 1) for i in range(n))
col_constraints = amp.sum(amp.less_equal(amp.sum(x[:, j]), 1) for j in range(n))
constraint1 = row_constraints + col_constraints

constraint2a = sum(
    x[i][j] * x[k][l]
    for i in range(n) for j in range(i+1, n)
    for k in range(0, i)
    for l in range(i+1, j) 
)
constraint2b = sum(
    x[i][j] * x[k][l]
    for i in range(n) for j in range(i+1, n)
    for k in range(i + 1, j) 
    for l in range(j + 1, n)
)   
constraint2 = constraint2a + constraint2b

constraint3 = sum(
    x[i][j]
    for i in range(n) for j in range(i + 1, min(n, i + 4))
)

weight1 = 1e10
weight2 = 1e10
weight3 = 1e10
model = obj + constraint1*weight1 + constraint2*weight2 + constraint3*weight3

def PredictedStructure():
    client = amp.FixstarsClient()
    client.parameters.timeout = timedelta(milliseconds=10000)
    client.token = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" # Replace with your token
    result = amp.solve(model, client=client)

    values = result.best.values
    selected_values = x.evaluate(values)


    structure = ['.'] * n
    paired = []
    for i, j in ValidPairs:
        if selected_values[i, j] == 1 and i < j:
            structure[i] = '('
            structure[j] = ')'
            paired.append((i, j))
            
    return structure, int(-result.best.objective)


"""
Nussinov Algorithm
This algorithm is a dynamic programming approach to predict the secondary structure of RNA sequences.
"""


DP = np.zeros((n, n), dtype=int)

def delta(i, j):
    if j - i - 1 < 3: #avoid sharp turns
        return 0
    else:
        if (bases[i] == 'A' and bases[j] == 'U') or \
           (bases[i] == 'U' and bases[j] == 'A') or \
           (bases[i] == 'G' and bases[j] == 'C') or \
           (bases[i] == 'C' and bases[j] == 'G') or \
           (bases[i] == 'G' and bases[j] == 'U') or \
           (bases[i] == 'U' and bases[j] == 'G'):
            return 1
        else:
            return 0

def recursion(i, j):
    if i >= j:
        return 0
    else:
        case1 = DP[i + 1][j]
        case2 = DP[i][j - 1]
        case3 = DP[i + 1][j - 1] + delta(i, j)
        case4 = 0
        for k in range(i, j):
            temp = DP[i][k] + DP[k + 1][j]
            if temp > case4:
                case4 = temp    
        
        return max(case1, case2, case3, case4)
    

def traceback(i, j, structure):
    if i >= j:
        return
    if DP[i][j] == DP[i + 1][j]:
        structure[i] = '.'
        traceback(i + 1, j, structure)
    elif DP[i][j] == DP[i][j - 1]:
        structure[j] = '.'
        traceback(i, j - 1, structure)
    elif DP[i + 1][j - 1] + delta(i, j) == DP[i][j]:
        structure[i] = '('
        structure[j] = ')'
        traceback(i + 1, j - 1, structure)
    else:
        for k in range(i + 1, j):
            if DP[i][j] == DP[i][k] + DP[k + 1][j]:
                traceback(i, k, structure)
                traceback(k + 1, j, structure)
                break

def nussinov():
    for i in range(n):
        I = 0
        for j in range(i+1, n):
            DP[I][j] = recursion(I, j)
            I += 1
    structure = ['.'] * n
    traceback(0, n - 1, structure)
    return structure, DP[0][n - 1]


def main():
    iter = 5 
    match = 0
    mismatch_low = 0
    mismatch_high = 0
    nu, goal_score = nussinov()
    print("Score by DP: ", goal_score)
    print("Structure: ", "".join(nu))
    for i in range(iter):
        amp, score = PredictedStructure()
        if amp == nu:
            match += 1
            print("MATCHED")
        else:
            if score < goal_score:
                mismatch_low += 1
                print("Score by Amplify: ", score)
            else:
                mismatch_high += 1
                print("Structure Amplify: ", "".join(amp))
    print("result: ", match/iter, (match + mismatch_high)/iter)
    
if __name__ == "__main__":
    main()
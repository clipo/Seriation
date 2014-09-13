
import random

N=20
edges=5
X = range(N)
random.seed() # uses system time to initialize random number generator, or you can pass in a deterministic seed as an argument if you want

# code to use to generate K pairs
A = random.sample(X,2*edges) # now you have a list of 2*K unique elements from 0 to N-1
pairs = zip(A[0:edges],A[edges:(2*edges)]) # now you have your pairs
print pairs

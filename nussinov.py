import argparse, os, sys
import numpy as np

def delta(xi, xj, seq):
    # xi and xj are characters from a sequence
    # return 0 or 1 indicating whether they are a complementary pair
    if (seq[xi] == 'A' and seq[xj] == 'U') or (seq[xi]=='U' and seq[xj]=='A'):
        return 1;
    elif (seq[xi] =='C' and seq[xj]=='G') or (seq[xi]=='G' and seq[xj]=='C'):
        return 1;
    else:
        return 0;

def bifurcation(i,j,seq):
    # Don't want to accidentally look at gamma[-1] (which is actually the last position)
    k = range(i+1,j)
    maximum=0
    for n in k:
        temp = setGamma(i,n,seq)+setGamma(n+1,j,seq)
        if maximum > temp:
            maximum=temp
            
    return maximum
    
## DEBUG: Not sure if this is working correctly, need to get someone to check ###
    # Apparently there are 1549056 calls to setGamma for example 1 with 12 bases
def setGamma(i,j,seq):
    L = len(seq)
    setGamma.counter += 1
    if (i > L) or (j > L):
        return 0
    elif (i==j) or (i-1==j and i>=1):
        return 0
    else:
        m = max(setGamma(i+1,j, seq), setGamma(i,j-1, seq), setGamma(i+1,j-1,seq)+delta(i,j,seq))
        b = 0
        if (j-i) > 1:
            b = bifurcation(i,j,seq)
        return max(m,b)
            

def main(args):
    seq_file_path = args.sequence
    output_file_path = args.out
    setGamma.counter = 0
    
    sys.setrecursionlimit(4000)
    # File import and split string into character list
    seq_file = open(seq_file_path,"r")
    seq = list(seq_file.read())
    L = len(seq)
  
    # Recursively populate matirx
    gamma = np.empty((L,L))
    gamma[:] = np.NAN
    for i in range(L):
        gamma[i,i] = 0
        if i>0:
            gamma[i,i-1] = 0
    
    for l in range(L):
        for i in range(L-l):
            j = i+l
            gamma[i][j] = setGamma(i,j, seq)

    #print(gamma)

    # Traceback
    stack = []
    trace = []
    stack.append((0,L-1))
    # stack will be false when it's empty
    while stack:
        t = stack.pop()
        i = t[0]
        j = t[1]
        if i >= j:
            continue
        elif gamma[i+1][j] == gamma[i][j]:
            stack.append((i+1,j))
        elif gamma[i][j-1]==gamma[i][j]:
            stack.append((i,j-1))
        # ASK WHEN DUE DATE IS CLOSER: How to tie-break between (5,10) and (6,10)?
        # Current method returns 6, 10, example returns 5,10
        elif gamma[i+1][j-1]+delta(i,j,seq) == gamma[i][j]:
            # Eliminate pairings with less than four bases between them
            if abs(i-j) >= 4:
                trace.append((i,j))
            stack.append((i+1,j-1))
        else:
            for k in range(i+1,j):
                if gamma[i][k]+gamma[k+1][j] == gamma[i][j]:
                    stack.append((k+1,j))
                    stack.append((i,k))
                    break
    #print(trace)
    
    # File output
    with open(output_file_path,'w') as fp:
        fp.write('\n'.join('{} {}'.format(x[0]+1,x[1]+1) for x in trace))


if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('sequence',
                      help='sequence file path.',
                      type=str)
  parser.add_argument('--out',
                      help='output file path.',
                      type=str,
                      default='out.txt')

  args = parser.parse_args()
  main(args)
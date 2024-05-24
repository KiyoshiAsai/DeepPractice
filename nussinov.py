# nussinov algorithm
# prediction of the rna secondary structure.
# input: rna sequence
# output: rna secondary structure

def nussinov(seq):
    n = len(seq)
    # initialize the matrix
    mat = [[0 for i in range(n)] for j in range(n)]
    
    # fill the matrix
    for length in range(1, n):
        for i in range(n - length):
            j = i + length
            if j - i < 4:
                mat[i][j] = 0
            else:
                pair = (seq[i] == 'A' and seq[j] == 'U') or (seq[i] == 'U' and seq[j] == 'A') or (seq[i] == 'C' and seq[j] == 'G') or (seq[i] == 'G' and seq[j] == 'C')
                mat[i][j] = max(mat[i + 1][j], mat[i][j - 1], mat[i + 1][j - 1] + pair)
                for k in range(i + 1, j):
                    mat[i][j] = max(mat[i][j], mat[i][k] + mat[k + 1][j])
    
    # traceback
    def traceback(i, j, pairs):
        if i < j:
            if mat[i][j] == mat[i + 1][j]:
                traceback(i + 1, j, pairs)
            elif mat[i][j] == mat[i][j - 1]:
                traceback(i, j - 1, pairs)
            elif mat[i][j] == mat[i + 1][j - 1] + ((seq[i] == 'A' and seq[j] == 'U') or (seq[i] == 'U' and seq[j] == 'A') or (seq[i] == 'C' and seq[j] == 'G') or (seq[i] == 'G' and seq[j] == 'C')):
                pairs.append((i, j))
                traceback(i + 1, j - 1, pairs)
            else:
                for k in range(i + 1, j):
                    if mat[i][j] == mat[i][k] + mat[k + 1][j]:
                        traceback(i, k, pairs)
                        traceback(k + 1, j, pairs)
                        break

    pairs = []
    traceback(0, n - 1, pairs)
    
    # express secondary structure by . and ()
    ss = ['.' for _ in range(n)]
    for i, j in pairs:
        ss[i] = '('
        ss[j] = ')'
    
    return ''.join(ss)

# create a test case
seq = 'GGGAAAUCC'
print(nussinov(seq))

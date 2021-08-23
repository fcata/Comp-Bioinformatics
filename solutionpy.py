def compute_matrices(seq1, seq2, substitution_matrix, gap_penalty, modality):
    '''Takes in input two sequences, a substitution matrix, a gap penalty and a modality (nw|sw).
    Returns a score matrix and a backtrack matrix following respectively the needleman-wunsh or the smith-waterman algorithm, depending on the modality'''

    a = len(seq1)
    b = len(seq2)

    score_matrix = [[0 for j in range(b + 1)] for i in range(a + 1)]
    bt_matrix = [[0 for j in range(b + 1)] for i in range(a + 1)]

    if modality == "nw":
        for i in range(1, a + 1):
            score_matrix[i][0] = score_matrix[i - 1][0] - gap_penalty
            bt_matrix[i][0] = 3
            score_matrix[0][i] = score_matrix[0][i - 1] - gap_penalty
            bt_matrix[0][i] = 1

    for i in range(1, a + 1):
        for j in range(1, b + 1):
            one = score_matrix[i][j - 1] - gap_penalty  # gap in first sequence
            two = score_matrix[i - 1][j - 1] + substitution_matrix[seq1[i - 1]][seq2[j - 1]]  # match
            three = score_matrix[i - 1][j] - gap_penalty  # gap in second sequence

            if modality == "nw":
                m = max(one, two, three)
            else:
                m = max(0, one, two, three)

            score_matrix[i][j] = m

            if m == one:
                bt_matrix[i][j] = 1
            elif m == two:
                bt_matrix[i][j] = 2
            elif m == three:
                bt_matrix[i][j] = 3
        # else if m == 0:
        # bt_matrix[i][j] remains equal to 0

    return score_matrix, bt_matrix


def sw_back_tracking(seq1, seq2, score_matrix, bt_matrix, threshold):
    '''Takes in input two sequences, the matrices given in output by the function compute_matrices in modality sw and a threshold.
    Returns every local alignment with a score higher or equal to the threshold. If threshold<0, returns only the best scoring alignment'''

    if threshold < 0:
        for row in score_matrix:
            for score in row:
                if score > threshold:
                    threshold = score

    return_list = []

    for a in range(len(score_matrix)):
        for b in range(len(score_matrix[a])):
            if score_matrix[a][b] >= threshold:  # Begin backtracking from this position
                A = ""
                B = ""
                i = a
                j = b
                while bt_matrix[i][j] != 0:  # Stops when we find a 0
                    if bt_matrix[i][j] == 1:  # Go left
                        A = '-' + A
                        B = seq2[j - 1] + B
                        j -= 1
                    elif bt_matrix[i][j] == 2:  # Go diagonal
                        A = seq1[i - 1] + A
                        B = seq2[j - 1] + B
                        j -= 1
                        i -= 1
                    elif bt_matrix[i][j] == 3:  # Go up
                        A = seq1[i - 1] + A
                        B = '-' + B
                        i -= 1
                return_list.append((A, B, score_matrix[a][b]))
    return return_list


def nw_back_tracking(seq1, seq2, score_matrix, bt_matrix):
    '''Takes in input two sequences and the matrices given in output by the function compute_matrices in modality nw.
    Returns the global alignment and its score'''

    return_list = []
    A = ""
    B = ""
    i = len(seq1)
    j = len(seq2)
    while i != 0 or j != 0:  # Stop when you arrive at the opposite corner
        if bt_matrix[i][j] == 1:  # Go left
            A = '-' + A
            B = seq2[j - 1] + B
            j -= 1
        elif bt_matrix[i][j] == 2:  # Go diagonal
            A = seq1[i - 1] + A
            B = seq2[j - 1] + B
            j -= 1
            i -= 1
        elif bt_matrix[i][j] == 3:  # Go right
            A = seq1[i - 1] + A
            B = '-' + B
            i -= 1
    return_list.append((A, B, score_matrix[len(seq1)][len(seq2)]))

    return return_list


def printer1(matrix):
    for row in matrix:
        print(row)
    print()


def printer2(alns):
    for aln in alns:
        print("Aln1:\t" + aln[0])
        print("Aln2:\t" + aln[1])
        print("Score:\t" + str(aln[2]))
    print()


from solution1 import compute_matrices, printer

score_matrix, bt_matrix = compute_matrices(seq1, seq1, substitution_matrix, gap_penalty, "sw")
printer(score_matrix)
print(bt_matrix)
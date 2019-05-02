
#!/usr/bin/env python
import numpy as np
import math
import itertools

# Dynamic Programming Iterative Implementation of Coin Change Problem
def num_of_ways(coins, num_coins, amount): 
  
    # Initialize all table values as 0. Table will store number of solutions
    info = [0 for k in range(amount+1)] 
  
    # Base case = if amount is 0
    info[0] = 1
  
    # Pick coins one by one and update the info[] values 
    for i in range(0,num_coins): 
        for j in range(coins[i],amount+1): 
            info[j] += info[j-coins[i]] 
  
    return info[amount] 

def combinations(sequence, length, NULL=object()):
    """ Find all combinations of the specified length from a sequence. """
    if length <= 0:
        combos = [NULL]
    else:
        combos = []
        for i, item in enumerate(sequence, 1):
            rem_items = sequence[i:]
            rem_combos = combinations(rem_items, length-1)
            combos.extend(item if combo is NULL else [item, combo]
                            for combo in rem_combos)
    return combos


def pairs(*lists):
    for t in itertools.combinations(lists, 2):
        for pair in itertools.product(*t):
            #Don't output pairs containing duplicated elements 
            if pair[0] != pair[1]:
                yield pair

def sum_combos(combos):
    G_win_count = 0
    for x in combos:
        G_sum = 0
        O_sum = 0
        for y in x[0]:
            G_sum = G_sum + y
        for z in x[1]:
            O_sum = O_sum + z
        if G_sum > O_sum:
            G_win_count += 1
    return G_win_count / len(combos)

def H_algorithm(matrix):
    mat_size = 15
    # Step 0.5: Subtract Largest Element
    max_val = np.max(matrix)
    new_matrix = max_val - matrix

    # Step 1: Subtract Local Minima
    row_minimums = np.amin(new_matrix, axis = 1)  
    count = 0
    subtracted_rows_matrix = np.empty_like(new_matrix)
    for x in new_matrix:
        x = x - row_minimums[count]
        temp = np.array([x])
        if count == 0:
            subtracted_rows_matrix = temp
        else:
            subtracted_rows_matrix = np.vstack((subtracted_rows_matrix, temp))
        count += 1

    # Step 2: Subtract Column Minima
    column_minimums = np.amin(subtracted_rows_matrix, axis = 0)
    count = 0
    subtracted_columns_matrix = np.empty_like(new_matrix)
    for y in subtracted_rows_matrix.T:
        if column_minimums[count] != 0:
            y = y - column_minimums[count]
        temp = np.array([y])
        if count == 0:
            subtracted_columns_matrix = temp
        else:
            subtracted_columns_matrix = np.concatenate((subtracted_columns_matrix, temp), axis = 0)
        count += 1
    matrix_to_cover = subtracted_columns_matrix.T

    # Step 3: Cover Zeros w/ Min Number of Lines
    cover_count = 0
    while cover_count < mat_size:
        cover_count = 0
        covered = np.zeros((mat_size,mat_size), dtype = int)
        # Make Horizontal Line If >= 2 Zeros In Row
        for i in range(0, mat_size):
            num_zeros_rows = np.count_nonzero(matrix_to_cover[i] == 0)
            if num_zeros_rows >= 2:
                covered[i] = covered[i] + 1
                cover_count += 1
        # Make Vertical Line If Zero Present
        for i in range(0, mat_size):
            index_zeros_col = np.nonzero(matrix_to_cover.T[i] == 0)[0]
            for h in range(0, np.size(index_zeros_col)):
                if covered.T[i][index_zeros_col[h]] == 0:
                    covered.T[i] = covered.T[i] + 1
                    cover_count += 1
        # If # of Lines Less Than Dimensions
        if cover_count < mat_size:
            matrix_to_cover = new_zeros(matrix_to_cover, covered, mat_size)
        else:
        # Number of Lines Equal to Dimensions
            return optimized_sum(matrix_to_cover, matrix, mat_size)

# Find Optimal Assignment from Orginal Matrix
def optimized_sum(matrix_to_cover, matrix, mat_size):
    sum = 0
    visited = [[], []]
    for x in range(0, mat_size):
        visited[0].append(False)
        visited[1].append(False)
    # Complete process until all rows and columns visited
    while False in visited[0] and False in visited[1]:
        # If 1 0 in row, add value in matrix to sum
        for x in range(0, mat_size):
            index_zeros_rows = np.nonzero(matrix_to_cover[x] == 0)[0]
            if np.size(index_zeros_rows) == 1:
                sum = sum + matrix[x][index_zeros_rows[0]]
                matrix_to_cover[x][index_zeros_rows[0]] += 1
                visited[0][x] = True
                index_zeros_col = np.nonzero(matrix_to_cover.T[index_zeros_rows[0]] == 0)[0]
                # Remove 0 in same column
                for y in index_zeros_col:
                    if y != index_zeros_rows[0]:
                        matrix_to_cover[x][y] = matrix_to_cover[x][y] + 1
        # If 1 0 in column, add value to sum
        for x in range(0, mat_size):
            index_zeros_col = np.nonzero(matrix_to_cover.T[x] == 0)[0]
            if np.size(index_zeros_col) == 1:
                sum = sum + matrix.T[x][index_zeros_col[0]]
                matrix_to_cover.T[x][index_zeros_col[0]] += 1
                visited[1][x] = True
                index_zeros_rows = np.nonzero(matrix_to_cover[index_zeros_col[0]] == 0)[0]
                # Remove 0 from same row
                for y in index_zeros_rows:
                    if y != index_zeros_rows[0]:
                        matrix_to_cover.T[x][y] = matrix_to_cover.T[x][y] + 1

    return sum

def new_zeros(matrix_to_cover, covered, mat_size):
    # Step 4: Create Additional Zeros
    uncovered_min = -1
    # Find Lowest Value Uncovered
    all_values = np.argsort(matrix_to_cover.flat)
    for x in all_values:
        if uncovered_min == -1:
            if covered.flat[x] == 0:
                uncovered_min = matrix_to_cover.flat[x]
    for row in range(0, mat_size): 
        for column in range(0, mat_size):
            # Subtract Min From Uncovered
            if covered[row, column] == 0:
                matrix_to_cover[row, column] = matrix_to_cover[row, column] - uncovered_min
            # Add Min To Covered Twice
            elif covered[row, column] == 2:
                matrix_to_cover[row, column] = matrix_to_cover[row, column] + uncovered_min
    return matrix_to_cover


def main():
    #1. Pandora Currency Problem
    Target_Value = 500
    Denominations = [1, 5, 10, 20, 50, 100]
    
    possible_ways = num_of_ways(Denominations, len(Denominations), Target_Value)
    print(possible_ways)

    #2. Different Number of Differing Sided Dice
    five_sided_dice = [1, 2, 3, 4, 5]
    ten_sided_dice = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    Gregor = [1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
              4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5]
    Oberyn = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7,
              8, 8, 8, 8, 9, 9, 9, 9, 10, 10, 10, 10]
    G_poss = set(itertools.combinations(Gregor, 8))
    O_poss = set(itertools.combinations(Oberyn, 4))

    poss_pairs = set(pairs(G_poss, O_poss))
    prob = sum_combos(poss_pairs)
    print(prob)
    

    #3. Hungarian Algorithm Sum of Matrix
    matrix = np.array([[7, 53, 183, 439, 863, 497, 383, 563, 79, 973, 287, 63, 343, 169, 583],
                      [627, 343, 773, 959, 943, 767, 473, 103, 699, 303, 957, 703, 583, 639, 913],
                      [447, 283, 463, 29, 23, 487, 463, 993, 119, 883, 327, 493, 423, 159, 743],
                      [217, 623, 3, 399, 853, 407, 103, 983, 89, 463, 290, 516, 212, 462, 350], 
                      [960, 376, 682, 962, 300, 780, 486, 502, 912, 800, 250, 346, 172, 812, 350],
                      [870, 456, 192, 162, 593, 473, 915, 45, 989, 873, 823, 965, 425, 329, 803],
                      [973, 965, 905, 919, 133, 673, 665, 235, 509, 613, 673, 815, 165, 992, 326],
                      [322, 148, 972, 962, 286, 255, 941, 541, 265, 323, 925, 281, 601, 95, 973],
                      [445, 721, 11, 525, 473, 65, 511, 164, 138, 672, 18, 428, 154, 448, 848],
                      [414, 456, 310, 312, 798, 104, 566, 520, 302, 248, 694, 976, 430, 392, 198],
                      [184, 829, 373, 181, 631, 101, 969, 613, 840, 740, 778, 458, 284, 760, 390],
                      [821, 461, 843, 513, 17, 901, 711, 993, 293, 157, 274, 94, 192, 156, 574],
                      [34, 124, 4, 878, 450, 476, 712, 914, 838, 669, 875, 299, 823, 329, 699],
                      [815, 559, 813, 459, 522, 788, 168, 586, 966, 232, 308, 833, 251, 631, 107],
                      [813, 883, 451, 509, 615, 77, 281, 613, 459, 205, 380, 274, 302, 35, 805]])

    max_sum = H_algorithm(matrix)
    print(max_sum)


main()  
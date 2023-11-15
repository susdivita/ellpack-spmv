"""
Generate a diagonal matrix with random values in (0,1], setting how many diagonals
it should have.
Also generates a random vector for that size, with values in [0,1].
"""
import argparse
import random


parser = argparse.ArgumentParser(prog="Diagonal matrix generator")
parser.add_argument('-s', '--size', required=True, type=int,
                    help='Size of matrix')
parser.add_argument('-d', '--diag', required=True, type=int,
                    help="Number of non-zero diagonals in the matrix")
parser.add_argument('--sym', action='store_true', help="Symetric matrix")
args = parser.parse_args()

assert 2 <= args.size, "Size must be at least 2"
assert 0 < args.diag <= (2 * args.size - 1), "Can't have more than 2N-1 diags"


nz_diags = [] # indices of non-zero diagonals
nz_elems = 0 # number of non-zero elements

if args.sym:
    if args.diag % 2 == 1:
        nz_diags.append(0)
    tmp = random.sample(list(range(1, args.size)), args.diag // 2)
    tmp_neg = [-i for i in tmp]
    nz_diags += tmp + tmp_neg
else:
    nz_diags = random.sample(list(range(-args.size + 1, args.size)), args.diag)
nz_elems = sum([args.size - abs(i) for i in nz_diags])



PATH = "../inputs/SquareDiagonal"
MTX_FILE = f"{PATH}/test_diag_nzdiag{args.diag}_{args.size}.mtx"
VCT_FILE = f"{PATH}/test_vec_{args.size}_random_in.txt"


with open(MTX_FILE, "wt") as fout:
    fout.write(f"{args.size} {args.size} {nz_elems}\n")
    for d in nz_diags:
        for i in range(args.size - abs(d)):
            row, col = i, i
            if d > 0:
                col += d
            else:
                row += (-d)
            fout.write(f"{col+1} {row+1} {random.randrange(0,100)/100}\n")

with open(VCT_FILE, "wt") as fout:
    fout.write(f"{args.size}\n")
    for _ in range(args.size):
        fout.write(f"{random.randrange(0,100)/100}\n")
import argparse

# ---------------------
# Parse input arguments
# ---------------------
parser = argparse.ArgumentParser(description="Convert sfs (Tanya file format) to dadi file format for sfs")
parser.add_argument("--num_bin",required=True,help="Number of bin, which is equal to the number of unfolded sfs bins plus "
                                                 "one bin for the monomorphic bin. Ex: if there are 7 individuals, input "
                                                 "15 here because there are 7*2=14 unfolded sfs bins plus one bin for the monomorphic bin.")
parser.add_argument("--folded_or_unfolded",required=True,help="Input either folded or unfolded.")
parser.add_argument("--population_name",required=True,help="Input the name of the population.")
parser.add_argument("--sfs_filename",required=True,help="Input the path to the sfs file with Tanya file format.")
parser.add_argument("--num_individuals",required=True,help="Input the number of individuals in the sample.")
parser.add_argument("--out_filename",required=True,help="Input the path to the output file.")


args = parser.parse_args()

outfile = open(args.out_filename, "w")

# Print first line
print (args.num_bin + " " + args.folded_or_unfolded + " " + '"' + args.population_name + '"', file=outfile)

# Print second line
sfs_bins = [str(0)] #initilize with 0 for the monomorphic bin, because we are going to mask it out anyway
mask = [str(1)]
with open(args.sfs_filename, "r") as f:
    for line in f:
        sfs_bins.append(str(int(float(line.rstrip("\n").split("\t")[1]))))
        mask.append(str(0))

# append the 0s to the rest of the bin
for i in range(int(args.num_individuals)):
    sfs_bins.append(str(0))
    mask.append(str(1))

print (" ".join(sfs_bins), file=outfile)
print (" ".join(mask), file=outfile)



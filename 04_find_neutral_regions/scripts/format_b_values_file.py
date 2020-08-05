# In this script, I want to format the b value file. In the file that was downloaded, it's just 1, 2, 3, etc and I want to convert it to chr1, chr2, chr3, etc...

outfile = open("/agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/bscores_chrfmt_high.tsv", "w")
with open("/agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/download/bscores.tsv", "r") as f:
    for line in f:
        items = line.rstrip("\n").split("\t")
        if float(items[3]) >= 900:
            if items[0] != "23":
                new_items = ["chr" + items[0], items[1], items[2], items[3]]
                print ("\t".join(new_items), file=outfile)
            else:
                new_items = ["chrX", items[1], items[2], items[3]]
                print("\t".join(new_items), file=outfile)
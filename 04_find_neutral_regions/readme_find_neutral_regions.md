## This file documents how I am obtaining the putatively neutral regions

1. B values
   1. Obtaining B values:
     - The file `bscores.tsv.gz` was obtained from Adriana Arneson who downloaded it from the CADD annotation.
     - The format of this file is chr, start, end, b values
2. Find regions from UCSC tracks
  1. UCSC tracks:
    1. https://genome.ucsc.edu/cgi-bin/hgTables
  2. **Genic regions**: clade: Mammal; genome: Human; assembly: Dec. 2013 (GRCh38/hg38); group: Genes and Gene Predictions; track: GENCODE v32; table: knownGene; region: genome; output format: BED - browser extensible data; output file: GRCh38_gencodev32_genome_knownGene; file type returned: gzip compressed; Create one BED record per: Whole Gene
  3. **Conserved regions**: clade: Mammal; genome: Human; assembly: Dec. 2013 (GRCh38/hg38); group: Comparative Genomics; track: Conservation; table: 100 Vert. El (phastConsElements100way); region: genome; output format: BED - browser extensible data; output file: GRCh38_phastConsElements100way; file type returned: gzip compressed
  4. **Repeats**: clade: Mammal; genome: Human; assembly: Dec. 2013 (GRCh38/hg38); group: Repeats; track: RepeatMasker; table: msk; region: genome; output format: BED - browser extensible data; output file: GRCh38_repeats; file type returned: gzip compressed
  5. **CpG islands**: clade: Mammal; genome: Human; assembly: Dec. 2013 (GRCh38/hg38); group: regulation; track: CpG islands; table: cpgIslandExt; region: genome; output format: BED - browser extensible data; output file: GRCh38_cpgislands; file type returned: gzip compressed
3. Unzip these files:
  ```
  pwd
  /home/tphung3/scratch/Turkana_Demography/04_find_neutral_regions/download
  ```

  ```
  ls -lrth
  total 185M
  -rw-rw-r-- 1 tphung3 tphung3 9.3M Nov 27 06:11 GRCh38_gencodev32_genome_knownGene.gz
  -rw-rw-r-- 1 tphung3 tphung3  82M Nov 27 06:40 GRCh38_phastConsElements100way.gz
  -rw-rw-r-- 1 tphung3 tphung3  61M Nov 27 06:45 GRCh38_repeats.gz
  -rw-rw-r-- 1 tphung3 tphung3 320K Nov 27 07:35 GRCh38_cpgislands.gz
  -rw-rw-r-- 1 tphung3 tphung3  34M Nov 27 07:43 bscores.tsv.gz
  ```

  ```
  gunzip *
  ```

  ```
  ls -lrth
  total 682M
  -rw-rw-r-- 1 tphung3 tphung3   30M Nov 27 06:11 GRCh38_gencodev32_genome_knownGene
  -rw-rw-r-- 1 tphung3 tphung3  350M Nov 27 06:40 GRCh38_phastConsElements100way
  -rw-rw-r-- 1 tphung3 tphung3  198M Nov 27 06:45 GRCh38_repeats
  -rw-rw-r-- 1 tphung3 tphung3 1005K Nov 27 07:35 GRCh38_cpgislands
  -rw-rw-r-- 1 tphung3 tphung3  105M Nov 27 07:43 bscores.tsv
  ```
4. Process the bscore file:
  - We want to convert the naming of chromosome from 1 to chr1.
  - Also, we just want to keep the sites where bscore is greater or equal to 0.9
  ```
  python ~/scratch/Turkana_Demography/04_find_neutral_regions/scripts/format_b_values_file.py
  ```
5. Separate by chromosome for parallelizing
  ```
  ./separate_by_chr.sh
  ```
6. Lift over the coordinates for bscore from hg19 to GRCh38
  ```
  ~/scratch/Turkana_Demography/04_find_neutral_regions/scripts/bscores_hg19tohg38_liftover.sh
  ```
7. Create bed files representing putatively neutral regions
  1. Date: 11/27/2019, the putatively neutral regions are defined as:
    1. Not overlapping with genic regions, conserved regions, CpG islands, and repeats
    2. Bscores are greater or equal to 0.9
  2. Subtract out the genic regions, conserved regions, CpG islands, and repeats from the genome
    1. Obtain the chrom size from UCSC: `wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes`
    ```
    /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/scripts/generate_bed_for_chrom_size.sh
    ```
    2. Run bedtools:
    ```
    /agavescratch/tphung3/Turkana_Demography/04_find_neutral_regions/scripts/obtain_putatively_neutral_regions.sh
    ```

# Supermatrix_Pipeline-NM123
A Pipeline for creating supermatrix of NM123 datasets.

**Required software:**
CD-HIT v4.8.1, BLAST v2.16, HMMER v3.4, MAFFT v7.525, TrimAl v1.5, FastTree v2.1.11, and PhyloPyPruner v1.2.6

**Test system:**
Ubuntu 22.04.5 LTS (windows subsystem)

**Installation guide:**
These software packages can be installed using conda, and there are no strict version requirements.

**Demo data:**
Demo_data.zip

**Run time for Demo:**
9 minutes

**Required dataset:**
Ref-datasets.zip

**Instruction for use:**
python Supermatrix_Pipeline-NM123.py   (You can change the pathway in script)

**Expected output:**
Demo_data-NM123-concatenation.fasta, Demo_data-NM123-concatenation.fasttree-ref.tre, Demo_data-NM123-concatenation_bygene.nex, Gene_Count_Matrix.csv, partition.txt, statistics.csv, and Supermatrix_Pipeline.log.

**Note:**
The original NM126 dataset contained three duplicated genes; therefore, we corrected it here to the NM123 dataset.

**Reference:**
Baker et al. Phylogenomic analyses indicate the archaeal superphylum DPANN originated from free-living euryarchaeal-like ancestors. Nature Microbiology 10, 1593-1604 (2025).

**More information:** https://doi.org/10.6084/m9.figshare.31878814



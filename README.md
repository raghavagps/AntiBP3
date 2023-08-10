# **AntiBP3**
A tool for predicting and designing Antibacterial peptides for three groups of bacteria using the sequence information.
## Introduction
Antibacterial peptides (ABPs) are specific cationic AMPs against various bacterial species that disturb the bacterial cell membrane with different mechanisms. The ABPs can also hamper intracellular processes of the bacterial pathogen in different ways, such as modulating the enzymatic activity, protein degradation & synthesis and nucleic acid synthesis, which gives an advantage over traditional antibiotics in developing resistance. Hence, develop an Antibacterial prediction tool for gram-positive, gram-negative and gram-variable bacteria.
AntiBP3 is also available as a web server at https://webs.iiitd.edu.in/raghava/antibp3. Please read/cite the content about the AntiBP3 for complete information, including the algorithm behind the approach.
## Reference
Bajiya N, Choudhury S, Dhall A, Raghava GPS. AntiBP3: A hybrid method for predicting antibacterial peptides against gram-positive/negative/variable bacteria. <a href="https://europepmc.org/article/PPR/PPR697637">bioRxiv; 2023. DOI: 10.1101/2023.07.25.550443</a>.

## Pip installation
The pip version of AntiBP3 is also available for easy installation and usage of the tool. The following command is required to install the package 
```
pip install antibp3
```
To know about the available option for the pip package, type the following command:
```
antibp3 -h
```

## Standalone
The Standalone version of AntiBP3 is written in python3, and the following libraries are necessary for the successful run:
- scikit-learn
- Pandas
- Numpy
- blastp

## Minimum USAGE
To know about the available option for the standalone, type the following command:
```
python antibp3.py -h
```
To run the example, type the following command:
```
python3 antibp3.py -i example_file.fasta
```
This will predict if the submitted sequences are Antibacterial or non-Antibacterial for the chosen class of bacteria. It will use other parameters by default. It will save the output in "outfile.csv" in CSV (comma-separated variables).

## Full Usage
```
usage: antibp3.py [-h] 
                       [-i INPUT 
                       [-o OUTPUT]
                       [-s {1,2,3}]
		       [-j {1,2,3,4,5}]
		       [-t THRESHOLD]
                       [-e EVAL]
		       [-d {1,2}]
		       [-wd WORKING]
```
```
Please provide the following arguments for the successful run

Optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input: protein or peptide sequence(s) in FASTA format or single sequence per line in single letter code
  -o OUTPUT, --output OUTPUT
                        Output: File for saving results by default outfile.csv
  -s {1,2,3}, --source {1,2,3}
                        Source: 1:GP ABPs, 2:GN ABPs, 3:GV ABPs by default 1
 -j {1,2,3,4,5}, --job {1,2,3,4,5}
                        Job Type: 1:Predict, 2:Design, 3:BLAST Search 4:Motif Scan, 5:Protein Scan ; by default 1
  -t THRESHOLD, --threshold THRESHOLD
                        Threshold: Value between 0 to 1 by default 0.5 for GP ABPs, 0.45 for GN ABPs and 0.51 for GV ABPs
  -e EVAL, --eval EVAL  E-value for Blast search (Blast Search only), by default 0.01 for GP ABPs, 0.01 for GN ABPs and 0.001 for GV ABPs
  -w {8,9,10,11,12,13,14,15,16,17,18,19,20}, --winleng {8,9,10,11,12,13,14,15,16,17,18,19,20}
                        Window Length: 8 to 20 (scan mode only), by default 8
  -d {1,2}, --display {1,2}
                        Display: 1:ABPs only, 2: All peptides, by default 1
  -wd WORKING, --working WORKING
                        Working Directory: Location for writing results
```

**Input File:** It allows users to provide input in the FASTA format.

**Output File:** Program will save the results in the CSV format; if the user does not provide the output file name, it will be stored in "outfile.csv".

**Source:** User should provide sources 1,2, and 3; by default, it's 1 for GP ABPs.

**Working Directory:** This option allows users to set the working directory in which they want to get the output files.

**Job:** User is allowed to choose between three different modules, such as 1 for prediction, 2 for Designing, 3 for Blast scan, 4 for Motif scan and 5 for Protein scan; by default, it's 1.

**Threshold:** User should provide a threshold between 0 and 1; by default, 0.5 for GP ABPs, 0.45 for GN ABPs and 0.51 for GV ABPs

**Window length**: User can choose any pattern length between 8 and 20 in long sequences. This option is available for only protein scan module.

**e-value:** User is allowed to choose an e-value for Blast search between  0.0001 to 1e-20; by default, 0.01 for GP ABPs, 0.01 for GN ABPs and 0.001 for GV ABPs

**Display type:** This option allow users to fetch either only Antibacterial peptides by choosing option 1 or prediction against all peptides by choosing option 2.


AntiBP3 Package Files
=======================
It contains the following files; a brief description of these files is given below

INSTALLATION                    : Installations instructions

LICENSE                         : License information

README.md                       : This file provides information about this package

envfile                                : This file comprises paths for the database and blastp executable

antibp3.py                        : Main Python program.    

model                            : This folder contains the models

motif                              : This folder contains the list of motifs

example_file.fasta                       : Test example file contains peptide sequences in FASTA format

example_predict_output.csv      : Example output file for predict module

example_blast_output.csv    : Example output file for Blast scan module

example_proteinscan_output.csv  : Example output file for Protein scan module

example_motifscan_output.csv.     : Example output file for motif scan module

example_design_output.csv       : Example output file for design module

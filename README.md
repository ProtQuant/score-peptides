# SCORE_PEPTIDES
**Stage one of the project : Generate scores for peptides**

* srounce code:  **src/ProduceLabel9.py**. 

* **output** folder contains output files of source code

  * **score.csv** is the needed file for Stage two

* **flydata_example3**  folder contains all the peptides information

* **uniprot-proteome_UP000000803.fasta** contains all the protein formulas

  


# How To Use ProduceLabel9

* Software dependencies

  ```python
  # python 3.8
  # need to install openpyxl
  import copy
  import heapq
  import os
  import numpy
  import pandas as pd
  import time
  ```

* Input values

  ```python
  """
  Two input values:
  	fastaFile: string, path of one .fasta file that contains protein formulas
  	folderPath: string, path of one folder that contains the many files of peptides' infomation
  """
  fastaFile = '../uniprot-proteome_UP000000803.fasta'
  folderPath = '../flydata_example3'
  ```

* Output file **score.csv**, needed for stage two of the project

  * in default output path **../saved data/`scource_file_name`/output** folder
  * columns = [(index), peptide, score]

* To run the code

  * go to the folder contains the source file 
  * type `python ProduceLabel9.py` at command line
  * Total processing time with default input value (above) is around 8 hours

(Uncomment line 598~614 to see the accumulative output for each 10 peptides files)



# Scoring Rules

( Protein splitting rule: After Argenine (R) or Lysine (K), as long as they aren’t followed by Proline (P) )

1. Record the different proteins from a .fasta file

2. For each peptides file (.csv)

   ​	a.	Find the different peptides in the file and record their largest `intensity` value

   ​	b.	For each protein, record its component peptide(s) in the file, and sum up the largest 3 intensity of the component peptides as `protein_total` value.

   ​	c.	For each peptides, calculate the `rate(S)` of `intensity/protein_total`, according to its parent protein(s)

   ​			\-	score = 20 if rate == 100%	(The proteins were once considered as ones that should be ignored.)

   ​				 score = 50 if rate <= 20%

   ​				 score = 100 if rate > 20%

   ​	d.	If some peptides in the file is exclusive to some proteins, these proteins are considered to be present. For these proteins, delete their component peptides appearing in the file, and score their rest  component peptides (split according to the splitting rule) as 0.

3. Collect the peptides and their scores from all the .csv files.

   * The peptides scored as 20 should not duplicate with peptides scored as 50 or 100.
   * The peptides scored as 0 should not duplicate with other scored peptides.



# Default File Structure

```
- uniprot-proteome_UP000000803.fasta (input)
- flydata_example3 (input, dowanload from: https://pgb.liv.ac.uk/~tony/flydata/)
- src
	- ProduceLabel9.py (source code)
- saved data (output folder)
	- ProduceLabel9 (folder name will be the same as the source code)
		- Output
		    - score.csv
		      (columns = [(index), peptide, score])
		    - ignore.csv
		      (columns = [(index), filename, #ignoredProtein, #valuedProtein], but actually we doesn't ignore the proteins having only one peptide in a .csv file)
		- proteinInfo (if set proteinInfo to True)
			- proteinInfo.xlsx
				- sheet1: protein index, columns = [(index), protein]
				  sheet2: simple peptides, columns = [(index), sPeptide, protIndex]
				- You could use this file to check the corresponding protein formula for each protIndex,
				  and all the simple peptides with indices of their parent protein(S)
		- peptideInfo (if set peptideInfo to True) 
			- fileIndex_Debug filename.xlsx (info for each input .csv file)
				- sheet1: protein infomation, columns = [(index), protein, pepIndex, count, unseenPeptide, list_intensity, prot_total]
				  sheet2: peptide score, columns = [(index), index, peptide, intensity, occurrence, protIndex, prot_total, score_0, score]
				  sheet3: unseen peptide, columns = [(index), score, protIndex] 
				- You could use this file to check the scoring results for single file
```

* Meanings of different columns will be explained later.



# Explanation of ProduceLabel9


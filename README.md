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

   * The peptides scored as 20 should not duplicate with peptides scored as 50 or 200.
   * The peptides scored as 0 should not duplicate with other scored peptides.



# Default File Structure


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

   ​	d.	If some peptides in the file is exclusive to some proteins, these proteins are considered to be present. For these proteins, delete their component peptides appearing in the file, and score their rest  unseen component peptides (split according to the splitting rule) as 0.

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
				  and see all the simple peptides with indices of their parent protein(S)
		- peptideInfo (if set peptideInfo to True) 
			- fileIndex_Debug filename.xlsx (info for each input .csv file)
				- sheet1: protein infomation, columns = [(index), protein, pepIndex, count, unseenPeptide, list_intensity, prot_total]
				  sheet2: peptide score, columns = [(index), index, peptide, intensity, occurrence, protIndex, prot_total, score_0, score]
				  sheet3: unseen peptide, columns = [(index), score, protIndex] 
				- You could use this file to check the scoring results for single file
```

* Meanings of different columns will be explained later.



# Main Data Structures In ProduceLabel9

* **`protList`**

  A list of all the proteins from the proteins file (.fasta) without duplicates

  For storing the index of proteins

* **`protDict`**

  A dictionary with all the different proteins as its keys

  ```
  protein_formula: {'pepIndex': a list of indices of the protein's component peptides that can be seen from the .csv file,
                    'count': int, how many peptides from the .csv file are in this protein
                    'leftProtein': empty, no use for now
                    'unseenPeptides': a list of unseen peptides formula (an empty value does not mean all its component peptides appear in the .csv file, it is not empty only when we know their parent proteins are present)
                    'list_intensity': a list of intensity values of the seen peptides
                    'prot_total':int, sum of the largest 3 intensity value, set to -1 if list_intensity is empty}
  ```

* **`allPossiblePepDict`**

  A dictionary storing all the **simple peptides** and lists of indices of their parent protein(s)

  `simplePep_formula: [protein_indices]`

* **`protpepDictList`**

  A list storing the simple peptides in each protein in the same order  as `ProtList`

  ```
  [{simplePep_formula: protein_index_0, simplePep_formula: protein_index_0, simplePep_formula: protein_index_0, ...}
   {simplePep_formula: protein_index_1, simplePep_formula: protein_index_1, simplePep_formula: protein_index_1, ...}
   {simplePep_formula: protein_index_2, simplePep_formula: protein_index_2, simplePep_formula: protein_index_2, ...}
   ...]
  ```

  When updating the occurrence of peptides in proteins for each .csv file, the (concatenation of) simple peptides that appear in the file will be deleted from `protpepDictList`. This data will then be used to generate `unseen peptides` later.

* **`pepList`**

  A list of all the peptides from one peptides file (.csv) without duplicates

  For storing the index of peptides in one file

* **`pepDict`**

  A dictionary with all peptides from one peptides file (.csv) as its keys

  ```
  peptide_formula: {'intensity': largest intensity of the peptide from one .csv file,
                    'occurrence': the number of its parent proteins,
                    'protIndex': a list of the indices of its parent proteins}
  ```

* **`fileInfoDf`**

  A data frame storing the ignored proteins in peptides files (.csv). The ignored proteins have only one of its component peptides appearing in a peptides file. Actually we still use these peptides, so its only for recording the number.

  `columns=['fileName', '#ignoredProtein', '#valuedProtein']`

* **`finalScoreDf`**

  A data frame storing the final scores for peptides.

  `columns=['peptide', 'score']`



# Dataflow in ProduceLabel9
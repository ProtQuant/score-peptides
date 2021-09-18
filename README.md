# score-peptides
**Stage one of the project : Generate scores for peptides**

* srounce code:  **src/ProduceLabel9.py**. 

* **output** folder contains output files of source code

  * **score.csv** is the needed file for Stage two

* **flydata_example3**  folder contains all the peptides information

* **uniprot-proteome_UP000000803.fasta** contains all the protein formulas

  


# ProduceLabel9.py

## How to use ProduceLabel9

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
  * `python ProduceLabel9.py`




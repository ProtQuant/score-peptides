# score-peptides
**stage one of the project**

* The latest version of code is src/ProduceLabel7.py. It can run directly and  will process files in `folderpath`.
* The code is not finished yet. 
  * more test cases are needed



other info:

* flydata_example3 contains all the peptides information from the new data

* uniprot-proteome_UP000000803.fasta contains protein formula

* ignore Test.py, Testcases, Testcases2, readme.pdf

* will delete some useless files in the repository later.

  

# Cases that may cause problems

* peptide score = [20, 20, 20, 100, 100, 50]

  final score is 20

* protein = 'ABBKBBK'

  pep = 'BBK'

  protein.replace(pep, '')  --> leftProtein = 'A'

  Any protein with duplicate sub-string may have this problem.

  

# Explanation of ProduceLabel7.py

## How to use the code

1. **ProduceLabel7.py should be able to run directly without passing parameters to it**

   

2. **Two required input**

  * `fastaFile`:	string, path of one .fasta file that contains protein formulas
  * `folderPath`:    string, path of one folder that contains the many (.csv/.tsv) files of peptides' infomation

  

3. **Default output file structure**

* `flyquant`:    The folder should contain all the .csv files on peptides.

* `output`:    The folder contains two output files.

  * score.csv:    columns = ['peptide', 'score']

  * ignore.csv:    columns = ['filename', '#ignoredProtein', '#valuedProtein']

    ​	#ignoredProtein:    the number of proteins that contain only one peptide from the .csv file ('filename')

    ​	#valuedProtein:    he number of proteins that contain more than one peptide from the .csv file ('filename')

* `proteinInfo`:    The folder will contain two files if`proteinInfo` is set to be true

  ​	Details of the files will be explained below.

* `peptideInfo`:    The folder will contain sub-folders whose name is the same as the .csv files if`peptideInfo` is set to be true. 3 files will be in each sub-folder.

  ​	Details of the files will be explained below.
  
  

4. **Two debugging options**:

  * `proteinInfo`:    boolean, will generate two files in protein folder if set as True

    * protein index.csv:    columns = ['index', 'protein']

      ​	You could use this file to check the corresponding protein formula for each protIndex.

    * simple peptides.csv:	columns = ['sPeptide', 'protIndex']

      ​	You could see in this file all the simple peptides broke up from proteins in protein index.csv and which protein has the simple peptide.

  * `peptideInfo`:    boolean, will generate 3 files in peptide/*filename* folder for each .csv file if set as True

    * protein information.csv

      | COLUMN         | DESCRIPTION                                                  |
      | -------------- | ------------------------------------------------------------ |
      | index          | protein index                                                |
      | protein        | protein formula                                              |
      | pepIndex       | a list of indices of the protein's peptides that can be seen from the .csv files<br/>the peptide formula corresponding to the indices can be checked in peptide score.csv |
      | count          | how many peptides from the .csv file are in this protein     |
      | leftProtein    | the protein formula without the seen Peptides (not empty only when we know there are unseenPeptides from the protein) |
      | unseenPeptides | a list of peptides formula broken up from leftProtein        |
      | list_intensity | a list of intensity value of the seen peptides               |
      | prot_total     | sum of the largest 3 intensity value , <br>if there are less than 3 peptides, set this number as the sum of all the intensity values, <br/>set to -1 if list_intensity is empty, otherwise, the value should be positive |

    * peptide score.csv

      | COLUMN                                | DESCRIPTION                                                  |
      | ------------------------------------- | ------------------------------------------------------------ |
      | index                                 | peptide index (may be 2 columns)                             |
      | peptide                               | peptide formula                                              |
      | intensity                             | largest intensity of the peptide from the .csv file          |
      | occurrence                            | the number of proteins that contain the peptide (should obey the splitting rule) |
      | protIndex                             | index of one protein that contain the peptide, <br/>if the peptide appears x times in different proteins, the peptide should appear in x rows in this file. |
      | prot_total                            | prot_total of the protein corresponding to the protIndex     |
      | score_0<br/>$\in \{NaN\} \cup(0,1]$   | = intensity divides prot_total, if prot_total > 0<br/>= NaN, if  if prot_total < 0 |
      | score<br/>$\in \{NaN, 20, 50, 100\} $ | = 20, if score_0 == 1<br/>= 100, if 0.2 < score_0 < 1<br/>= 50, if score_0 <= 0.2<br/>= NaN, if score_0 == NaN |

    * unseen peptide.csv

      | COLUMN    | DESCRIPTION                                              |
      | --------- | -------------------------------------------------------- |
      | peptide   | peptide formula                                          |
      | score     | set to 0                                                 |
      | protIndex | list of indices of the proteins that contain the peptide |

    Microsoft Excel may have problem with opening files that contain a list of protein indices, try WPS office.
    
    

## Data Structure

Three major data structures are used in the code: list, dictionary and data frame.



1. **Lists are used for maintaining the index for peptides and proteins.**

  * `protList` is for indexing proteins in .fasta file.

  * The many `pepList` are for indexing the peptides in the .csv files respectively.

    

2. **Dictionary are used for storing the information of peptides and proteins without duplication.**

  * `protDict`

    * protein_formula: {
      
      ​	'pepIndex':    [], 
      
      ​	'count':    int, 
      
      ​	'leftProtein':    string, 
      
      ​	'unseenPeptides':    [], 
      
      ​	'list_intensity':    [],
      
      ​	'prot_total'：   int,    }
      
      The meaning of each value is the same as the ones explained above in protein information.csv.
      
    * The keys and its order (protein_formula) won't be changed after initialization.
    
      The values will be updated from scratch when processing each .csv file.
    
  * `pepDict`
  
    * peptide_formula: {
      
      ​	'intensity':    int,
      
      ​	'occurrence':    int,
      
      ​	'protIndex':    [],    }
      
      The meaning of each value is the same as the ones explained above in unseen peptide.csv, except that protIndex is a list of  indices of all related proteins rather then one index.
      
    * Items in this dictionary will be generated according to each .csv file.
    
  * `unseenPepDict`:    record unseen peptides for each .csv file.
  
    * peptide_formula: {
    
      ​	'score':    0, 
    
      ​	'protIndex':    [], list of indices of the proteins that contain the peptide,    }
      
      

3. **Data frame is used for collect & gather the information generated according to each .csv file.**

  * `all_scoreDf`:    columns=['peptide', 'score']

    ​	For each .csv file, the seen peptides and their scores will be appended to the end of `all_scoreDf`.

    ​	One peptide may appear in more than one rows with different scores.

  * `all_unsPepDf`:    columns=['peptide', 'score']

    ​	For each .csv file, the unseen peptides and their scores (=0) will be appended to the end of `all_unsPepDf`. Peptides may duplicate in different files.

  * `finalScoreDf`:    columns=['peptide', 'score']
    
    ​	Record the scores for the seen and unseen peptides for all the files without duplications.

  * `fileInfoDf`:   columns=['fileName', '#ignoredProtein', '#valuedProtein']

    ​	For each .csv file, record the number of proteins that have only one peptide appearing in this file.



## Procedures

For the tables in this part, first column is the name of the data structure. Other columns in grey are the data that will be updated in that step. Data that hasn't been updated from the initial default value is not presented.

1. `init_proteins()`

* input:   the .fasta file

* output:

  * | `protList`   | `index`     | `protein`    |
    | ---- | ---- | ---- |

  * | `protDict`   | `protein_formula`    |
    | ---- | ---- |

2. `find_all_possible_simple_peptides()`

* input：`protList`

* output:
  * | `allPossiblePepDict`   | `simplePep_formula`     | `protein_index`     |
    | ---- | ---- | ---- |

3. For each .csv file:

   a.  `init_peptides_from_one_file()`

   * input:    one .csv file
   * output:
     * | `pepList`   | `index`     | `peptide`     |
       | ---- | ---- | ---- |
     * | `pepDict`   | `peptide_formula`     | `intensity`     |
       | ---- | ---- | ---- |

   b.   `update_occurrence_of_one_file()`

   * input:    `allPossiblePepDict`, `protList`, `pepList`

   * output:   

     * | `pepDict`   | peptide_formula     | intensity    | `occurrence`    | `protIndex`    |
       | ---- | ---- | ---- | ---- | ---- |
     * | `protDict`   | protein_formula     | `pepIndex`     | `count`    | `list_intensity`     | `prot_total`      |
       | ---- | ---- |---- |---- |---- |---- |

     * `ignored`, `valued`

   c. `update_unseen_peptides_of_one_file()`

   * input:    `protList`, `pepList`, `pepDict[pep]['occurrence']`, `pepDict[pep]['protIndex']`, `protDict[prot]['pepIndex']` 
   * output:
     * | `protDict` | protein_formula | pepIndex  | count | `leftProtein` | `unseenPeptides` | list_intensity | prot_total |
       | ---- | ---- |---- |---- |---- |---- |---- |---- |
     * | `unseenPepDict`   | `peptide_formula`     | `score`     | `protIndex` |
       | ---- | ---- | ---- |---- |

   d.  `calculate_score_of_one_file()`

   ​	will expand the list of `protIndex` and fill in `pepDf['prot_total']` according to the index

   * input:   `pepDict[pep]['intensity']`, `pepDict[pep]['protIndex']`, `protDict[prot]['prot_total']`
   * output:    
     * | `pepDf`   | peptide | intensity | occurrence| `protIndex` | `prot_total` | `score_0` | `score` |          
       | ---- | ---- | ---- | ---- | ---- | ---- | ---- |---- |

   c.   gather information

   * append `pepDf['peptide']` and `pepDf['score']` to `all_scoreDf`

   * append  `unseenPepDf` to `all_unsPepDf`
   * append `filename`, `ignored`, `valued` to `fileInfoDf` 

4. Process the collected data

   * choose the most common score of one peptide in `all_scoreDf` 
   * drop duplications in `all_unsPepDf`
   * In `all_unsPepDf`, delete the peptides that appear in `all_scoreDf`
   * combine `all_unsPepDf` and`all_unsPepDf` as `finalScoreDf`

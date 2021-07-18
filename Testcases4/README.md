# Check List

* read in proteins correctly
  * one protein may have many lines
  * many proteins
  * break up the proteins correctly
* read in peptides correctly (one file)
  * largest intensity
  * read in all pep
* record occurrence, prot_total, score from one file
  * occurrence should obey splitting rule ('K', 'R', 'P')
  * may miss 'M' at the beginning
  * choose the sum of largest 3 intensity as prot_total
  * one peptide may have different scores in different proteins
  * score rule (0, 20, 50, 100)
* initial protein info should not change while processing different pep files
  * empty values of protDict
  * protpepDcitList
* choose the most common score among dif files as the final score 



*  write a program and some code to automatically test the scoring process
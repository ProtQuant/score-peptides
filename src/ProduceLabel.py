import pandas as pd
import csv
import re


def find_peptides_from_a_protein(proteinFormula):
    # read in a string of protein
    # return a list of possible peptides
    # peptides may be duplicate

    peptides = []
    pep = ""
    for amino_acid in proteinFormula:
        pep += amino_acid
        if (amino_acid == 'K' or amino_acid == 'R'):
            peptides.append(pep)
            pep = ""
    if not (proteinFormula.endswith('K') or proteinFormula.endswith('R')):
        peptides.append(pep)

    # peptides = re.split(r"(K|R)", proteinFormula)
    # print(type(peptides))
    # print(peptides)
    return peptides


def find_peptidesName_from_fasta(filepath, sort=False):
    # read in .fasta file to get the possible peptides
    # return a (sorted) dictionary of dif pepName with value of 0

    # fl = open('..\\flydata_example\\uniprot-proteome_UP000000803.fasta')
    # fl = open('..\\flydata_example\\toy_example.fasta')
    fl = open(filepath)
    peptides = {}
    quantify = 0
    protein = ''
    line = fl.readline().strip()
    while line:
        if '>' in line:
            if len(protein)>0:
                pepName = find_peptides_from_a_protein(protein)  # break down the protein before '>'
                for pep in pepName:
                    peptides[pep] = quantify
                protein = ''  # '>' is the start of a new protein
        else:
            protein += line
        line = fl.readline().strip()
    # the last protein
    pepName = find_peptides_from_a_protein(protein)
    for pep in pepName:
        peptides[pep] = quantify

    if sort:
        # peptides = sorted(peptides) # sorted list of keys
        peptides = {k: peptides[k] for k in sorted(peptides)}

    # print(peptides)
    return peptides

def update_peptidesCount_we_can_see_from_tsv(tsvFiles, pepDict):
    # read in the tsv file path and a dict of all possible peptides
    # return a dict with updated count

    for filepath in tsvFiles:
        df = pd.read_csv(filepath, sep='\t')
        # print(df)
        for i in range(len(df)):
            pep = df['peptide'][i]
            count = df['tot_num_ions'][i]
            if (pep not in pepDict):
                print(pep + " is not found in .fasta file, add to the dict")
                pepDict[pep] = count
            elif (pepDict[pep] < count):
                pepDict[pep] = count


    # print(pepDict1)
    return pepDict

def generate_outputFile(pepCount, filename):
    dataframe = pd.DataFrame({'peptide': pepCount.keys(), 'tot_num_ions': pepCount.values()})
    dataframe.to_csv(filename, index=False, sep='\t')

def main():

    ## toy

    # fastaFile = '..\\flydata_example\\toy_example.fasta'
    fastaFile = '..\\flydata_example\\test.fasta'
    #
    # tsvFiles = ['..\\flydata_example\\toy_example_search1.tsv',
    #             '..\\flydata_example\\toy_example_search2.tsv']
    #
    # outputFile = '..\\output\\o1.tsv'
    #
    pepDict = find_peptidesName_from_fasta(fastaFile, False)
    print(pepDict)
    # pepCount = update_peptidesCount_we_can_see_from_tsv(tsvFiles, pepDict1)
    # # print(pepCount)
    # generate_outputFile(pepCount, outputFile)

    ## fruit-flies
    # fastaFile = '..\\flydata_example\\uniprot-proteome_UP000000803.fasta'
    # tsvFiles = ['..\\flydata_example\\QE752_MSQ988_20170315_BenHopkins_R1V8.tsv',
    #             '..\\flydata_example\\QE752_MSQ988_20170315_BenHopkins_R2M8.tsv',
    #             '..\\flydata_example\\QE752_MSQ988_20170315_BenHopkins_R4M2.tsv']
    # outputFile = '..\\output\\o2.tsv'
    # pepDict1 = find_peptidesName_from_fasta(fastaFile, False)
    # pepCount = update_peptidesCount_we_can_see_from_tsv(tsvFiles, pepDict1)
    # generate_outputFile(pepCount, outputFile)


if __name__ == '__main__':
    main()

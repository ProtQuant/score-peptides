# 5/7/2021 using the second version of data
#  [[]] * len(pepDf) vs [[] for _ in xrange(len(pepDf))]: the first express : change one will change all

import pandas as pd
import time

from joblib.numpy_pickle_utils import xrange


def init_proteins(filepath):
    """
    initialize proteins and their index from .fasta file
    :param
        filepath: path of one .fasta file that describes all of the proteins believed to be in the species
    :return
        protDf: pandas dataframe
                    column: index, protein
    """

    # store the proteins into a list
    protList = []
    fl = open(filepath)
    protein = ''
    line = fl.readline().strip()
    while line:
        if '>' in line:
            if len(protein) > 0:  # when read in the first line, protein is empty
                if protein not in protList:
                    protList.append(protein)
                protein = ''  # '>' is the start of a new protein
        else:
            protein += line
        line = fl.readline().strip()
    # the last protein
    if protein not in protList:
        protList.append(protein)
    print("There are " + str(len(protList)) + " proteins in total from .fasta files")  # 22117

    # translate the list into a dataframe
    protDf = pd.DataFrame({'protein': protList})
    # print(protDf)

    return protDf


def generate_proteinIndex_file(protDf, filepath):
    """
    generate a file that contains the proteins and their indices
    """
    protDf.to_csv(filepath, index=True, header=True, sep=',')


def read_in_stored_proteins(filepath, colname=["protein"]):
    """

    :param filepath:
    :param colname:
    :return:
    """
    # read in the file
    if (filepath.endswith(".csv")):
        protDf = pd.read_csv(filepath, sep=',', usecols=colname)
    elif (filepath.endswith(".tsv")):
        protDf = pd.read_csv(filepath, sep='\t', usecols=colname)
    print(print("Read in " + str(len(protDf)) + " proteins."))
    # print(protDf)
    return protDf


def split_protein(protein):
    """
    break up the protein after 'K' or 'R' except that it is followed by 'P'
    :param
        protein: a string of protein formula
    :return:
        peptides: a list of peptides formula that broken up from the input
    """
    peptides = []
    pep = ""
    for amino_acid in protein:
        pep += amino_acid
        if (amino_acid == 'K' or amino_acid == 'R'):
            peptides.append(pep)
            if (pep.startswith('P') and len(peptides) > 1):
                pep2 = peptides.pop()
                pep1 = peptides.pop()
                peptides.append(pep1 + pep2)
            pep = ""
    if not (protein.endswith('K') or protein.endswith('R')):
        peptides.append(pep)
        if (pep.startswith('P') and len(peptides) > 1):
            pep2 = peptides.pop()
            pep1 = peptides.pop()
            peptides.append(pep1 + pep2)

    return peptides


def find_all_possible_simple_peptides(protDf):
    """
    break up all the proteins from .fasta file to get all possible simple peptides
    :param
        protDf: pandas dataframe, should have a column "protein" that contains all the protein formula
    :return:
        allPossiblePepDict: a dictionary storing the formula of all the possible simple peptides (no duplicate)
                            and the indices of proteins that generates them
                                simplePep_formula: [protein_index]
    """
    protIndex = 0
    allPossiblePepDict = {}
    for i in range(len(protDf)):
        prot = protDf['protein'][i]
        simplePepList = split_protein(prot)
        for pep in simplePepList:
            if pep not in allPossiblePepDict:
                allPossiblePepDict[pep] = []
            allPossiblePepDict[pep].append(protIndex)
        protIndex += 1

    print("There are " + str(len(allPossiblePepDict)) + " simple peptides in total from all the proteins")  # 464696

    return allPossiblePepDict


def generate_simplePep_file(allPossiblePepDict, filepath):
    peptide = allPossiblePepDict.keys()
    protIndex = allPossiblePepDict.values()
    sPepDf = pd.DataFrame({'peptide': peptide,
                              'protIndex': protIndex})
    # print(sPepDf.head(n=10))
    # print(list(allPossiblePepDict.keys())[:5])
    # print(list(allPossiblePepDict.values())[:5])
    sPepDf.to_csv(filepath, index=True, sep=',')


def read_in_stored_simple_peptides(filepath, colname=["peptide", "protIndex"]):
    """
    problem exits:
    1. I cannot use MST Excel to open the file on my compuer (WPS can do)
    2. When casting dataframe to dictionary, list of numbers are in the format of strings
    :param filepath:
    :param colname: first peptide
    :return:
    """
    # read in the file
    if (filepath.endswith(".csv")):
        sPepDf = pd.read_csv(filepath, sep=',', usecols=colname)
    elif (filepath.endswith(".tsv")):
        sPepDf = pd.read_csv(filepath, sep='\t', usecols=colname)
    allPossiblePepDict = sPepDf.set_index(colname[0]).T.to_dict('list')
    print(print("Read in " + str(len(allPossiblePepDict)) + " simple peptides."))
    print("The method still have problems to solve")
    # print(sPepDf.head(n=10))
    # print(list(allPossiblePepDict.keys())[:5])
    # print(list(allPossiblePepDict.values())[:5])
    return allPossiblePepDict


def init_peptides_from_one_file(filepath, colname=['peptide', 'intensity']):
    """
    read in the peptides from one .csv file, drop duplicate and choose the largest intensity
    :param
        filepath: string,
        colname: list of strings, the names of columns that need to be read in
                    the first element should be peptide_column
                    the second element should be the column for scoring
    :return:
        pepDf: panda dataframe
                column: index, peptide, intensity
    """
    # read in the file
    if (filepath.endswith(".csv")):
        pepDf = pd.read_csv(filepath, sep=',', usecols=colname)
    elif (filepath.endswith(".tsv")):
        pepDf = pd.read_csv(filepath, sep='\t', usecols=colname)

    # drop the lines that contain NaN
    pepDf = pepDf.dropna(axis=0, how='any')
    # print(pepDf)

    # find the max intensity for each peptide
    max_intensity = pepDf.groupby([colname[0]])[colname[1]].transform(max)
    # print(max_intensity)

    # drop duplicate
    pepDf = pepDf.loc[pepDf[colname[1]] == max_intensity]  # duplication may still exist since two max intensity of the same value may exist
    pepDf = pepDf.drop_duplicates()
    # print(pepDf)

    # reset index
    pepDf = pepDf.reset_index(drop=True)
    # print(pepDf)

    return pepDf


def update_occurrence_of_one_file(pepDf, allPossiblePepDict, protDf):
    """
    old version:
    will update pepDict[pep]['occurrence'],  pepDict[pep]['proteinIndex'], protDict_temp[prot]['seenPeptides'] + count
    :param pepDf:
    :param allPossiblePepDict:
    :param protDf:
    :return:
        pepDf: (index), peptide, intensity, protIndex, occurrence
        protDf_temp: (index), protein, pepIndex, count
    """

    # protDf should not be changed when handle different csv files # need test
    protDf_temp = protDf.copy(deep=True)

    # add the needed columns
    # pepDf["protIndex"] = [[]] * len(pepDf)  # unchangeable, all the [] refer to the same address
    pepDf["protIndex"] = [[] for _ in xrange(len(pepDf))]  # list of protein index that contain the peptide
    pepDf["occurrence"] = [0] * len(pepDf)
    protDf_temp["pepIndex"] = [[]for _ in xrange(len(protDf_temp))]
    protDf_temp["count"] = [0] * len(protDf_temp)  # how many peptides from pepDf is in this protein

    # print(pepDf)
    # print(protDf_temp)

    # update occurrence
    for pepIndex in range(len(pepDf)):
        pep = pepDf['peptide'][pepIndex]

        # first search it in allPossiblePepDict
        if pep in allPossiblePepDict:
            protIndex = allPossiblePepDict[pep]
            for i in protIndex:
                pepDf.loc[pepIndex, "protIndex"].append(i)
            # pepDf.loc[pepIndex, "protIndex"] = protIndex  #  ValueError: Must have equal len keys and value when setting with an ndarray
            pepDf.loc[pepIndex, "occurrence"] = len(protIndex)  # loc change the original df, df[][] may only change a view
            for i in protIndex:
                protDf_temp.loc[i, 'pepIndex'].append(pepIndex)

        # or search it through all the proteins
        else:
            occurrence = 0
            for protIndex in range(len(protDf_temp)):
                prot = protDf_temp['protein'][protIndex]

                # the next letter should not be a P
                if prot.find(pep + 'P') == -1:  # ignore for now the case that a peptide may appear twice in a protein, and one can be counted
                    # may locate at the head of a protein
                    if prot.startswith(pep):
                        occurrence += 1
                        pepDf.loc[pepIndex, "protIndex"].append(protIndex)
                        protDf_temp.loc[protIndex, 'pepIndex'].append(pepIndex)
                    # may locate at the head of a protein missing an M
                    elif prot.startswith('M' + pep):
                        occurrence += 1
                        pepDf.loc[pepIndex, "protIndex"].append(protIndex)
                        protDf_temp.loc[protIndex, 'pepIndex'].append(pepIndex)
                    # may locate in the protein following a K
                    elif prot.find('K' + pep) != -1:
                        occurrence += 1
                        pepDf.loc[pepIndex, "protIndex"].append(protIndex)
                        protDf_temp.loc[protIndex, 'pepIndex'].append(pepIndex)
                    # may locate in the protein following an R
                    elif prot.find('R' + pep) != -1:
                        occurrence += 1
                        pepDf.loc[pepIndex, "protIndex"].append(protIndex)
                        protDf_temp.loc[protIndex, 'pepIndex'].append(pepIndex)
            if occurrence == 0:
                print('can not find the ' + str(pepIndex) + 'th peptides: ' + pep + ' in the proteins')
            pepDf.loc[pepIndex, "occurrence"] = occurrence

    # print(protDf_temp.head(n=100)['pepIndex'])

    # update count
    protDf_temp['count'] = protDf_temp.apply(lambda x: len(x['pepIndex']), axis=1)
    # print(protDf_temp.head(n=30))

    return pepDf, protDf_temp




def main():
    start = time.perf_counter()

    # # init proteins
    fastaFile = '..\\flydata_example\\uniprot-proteome_UP000000803.fasta'
    protDf = init_proteins(fastaFile)
    # # store proteins
    # protIndexFile = '.\\stored data\\protein index.csv'
    # generate_proteinIndex_file(protDf, protIndexFile)

    # # split proteins
    allPossiblePepDict = find_all_possible_simple_peptides(protDf)
    # # store simple peptides
    # simplePepFile = '.\\stored data\\simple peptides.csv'
    # generate_simplePep_file(allPossiblePepDict, simplePepFile)

    # # read in stored proteins
    # protIndexFile = '.\\stored data\\protein index.csv'
    # protDf = read_in_stored_proteins(protIndexFile)
    # # read in stored simple peptides, problem exists
    # simplePepFile = '.\\stored data\\simple peptides.csv'
    # read_in_stored_simple_peptides(simplePepFile, colname=["peptide", "protIndex"])

    # # init pep from one file
    filepath = '.\\Testcases2\\Test1.csv'
    pepDf = init_peptides_from_one_file(filepath)

    # # update occurrence
    pepDf, protDf = update_occurrence_of_one_file(pepDf, allPossiblePepDict, protDf)
    pepDf.to_csv('.\\Testcases2\\Test1_pepDf.csv', index=True, sep=',')
    protDf.to_csv('.\\Testcases2\\Test1_protDf.csv', index=True, sep=',')

    end = time.perf_counter()
    print('finish running fruit-fly example, total processing time: \t' + str(end - start))





if __name__ == '__main__':
    main()

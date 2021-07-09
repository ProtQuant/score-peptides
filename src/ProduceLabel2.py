# 2021/6/21

# !!wrong, if using 'value', all pepDict[pep] will point to one address, change one will change all
#  pepDict[pep] = value

# problems: aaaakp'seenPep'aaak

import pandas as pd


def init_proteins(filepath):
    # index the proteins from .fasta file via a list
    # return a list and a dictionary of proteins

    protList = []

    protDict = {}
    # seenPeptides = []  # index of peptides that can be seen in the .tsv file
    # leftProtein = ''  # the protein formula that generated by remonving the seenPeptides from original protein
    # unseenPeptides = []  # formula of unseen peptides
    # value = {'seenPeptides': [], 'leftProtein': '', 'unseenPeptides': []}

    fl = open(filepath)
    protein = ''
    line = fl.readline().strip()
    while line:
        if '>' in line:
            if len(protein) > 0:  # when read in the first line, protein is empty
                protList.append(protein)
                # !!wrong, if using 'value', all protDict_temp[protein] will point to one address, change one will change all
                # protDict_temp[protein] = value
                protDict[protein] = {'seenPeptides': [], 'leftProtein': '', 'unseenPeptides': []}
                protein = ''  # '>' is the start of a new protein
        else:
            protein += line
        line = fl.readline().strip()
    # the last protein
    protList.append(protein)
    protDict[protein] = {'seenPeptides': [], 'leftProtein': '', 'unseenPeptides': []}

    # if need sorting, just sort protList
    return protList, protDict


def init_peptides(filepaths):
    # index the peptides from the .tsv files
    # return a list and dictionary of peptides

    pepList = []  # store peptides and its index

    pepDict = {}
    # value = {'total_num_ions': 0, 'occurrence': 0, 'proteinIndex': []}

    for filepath in filepaths:
        df = pd.read_csv(filepath, sep='\t')
        for i in range(len(df)):
            pep = df['peptide'][i]
            count = df['tot_num_ions'][i]
            if (pep not in pepDict):
                # !!wrong, if using 'value', all pepDict[pep] will point to one address, change one will change all
                # pepDict[pep] = value
                pepDict[pep] = {'total_num_ions': 0, 'occurrence': 0, 'proteinIndex': []}
                pepDict[pep]['total_num_ions'] = count
                pepList.append(pep)
            elif (pepDict[pep]['total_num_ions'] < count):
                pepDict[pep]['total_num_ions'] = count

    # if sorting is needed, just sort pepList
    return pepList, pepDict


def update_occurrence(protList, protDict, pepList, pepDict):
    # find the occurance of peptides in proteins
    # update pepDict1[pep]['occurrence'],  pepDict1[pep]['proteinIndex'], protDict1[prot]['seenPeptides']
    # return updated protDict1 and pepDict1

    pepIndex = 0
    for pep in pepList:
        occurrence = 0
        protIndex = 0
        for prot in protList:
            if prot.find(pep) != -1:
                occurrence += 1
                pepDict[pep]['proteinIndex'].append(protIndex)
                protDict[prot]['seenPeptides'].append(pepIndex)
            protIndex += 1
        pepDict[pep]['occurrence'] = occurrence
        pepIndex += 1

    return protDict, pepDict


def split_protein(protein):
    peptides = []
    pep = ""
    for amino_acid in protein:
        pep += amino_acid
        if (amino_acid == 'K' or amino_acid == 'R'):
            peptides.append(pep)
            pep = ""
            if (pep.startswith('P') and len(peptides) > 1):
                pep2 = peptides.pop()
                pep1 = peptides.pop()
                peptides.append(pep1 + pep2)
    if not (protein.endswith('K') or protein.endswith('R')):
        peptides.append(pep)
        if (pep.startswith('P') and len(peptides) > 1):
            pep2 = peptides.pop()
            pep1 = peptides.pop()
            peptides.append(pep1 + pep2)
    return peptides


def find_unseen_peptides_from_one_protein(protein, protDict, pepList):
    # update protDict_temp
    # return unseen peptides

    leftPortein = protein
    for pepIndex in protDict[protein]['seenPeptides']:
        leftPortein = leftPortein.replace(pepList[pepIndex], '')
    protDict[protein]['leftPortein'] = leftPortein
    unseenPeptides = split_protein(leftPortein)
    protDict[protein]['unseenPeptides'] = unseenPeptides
    return unseenPeptides


def update_all_unseen_peptides(protList, protDict, pepList, pepDict):
    for pep in pepList:
        if pepDict[pep]['occurrence'] == 1:
            protIndex = pepDict[pep]['proteinIndex'][0]
            protein = protList[protIndex]
            unseenPeptides = find_unseen_peptides_from_one_protein(protein, protDict, pepList)
            for pep in unseenPeptides:
                if pep not in pepDict:
                    pepDict[pep] = {'total_num_ions': 0, 'occurrence': 0, 'proteinIndex': []}
                pepDict[pep]['proteinIndex'].append(protIndex)
    return pepDict

def generate_outputFile(pepDict, filename):
    dataframe = pd.DataFrame({'peptide': pepDict.keys(),
                              'tot_num_ions': [pepDict[k]['total_num_ions'] for k in pepDict]})
    dataframe.to_csv(filename, index=False, sep='\t')


def main():
    protList, protDict = init_proteins('..\\flydata_example\\toy_example.fasta')
    print(protList)
    print(protDict)
    pepList, pepDict = init_peptides(['..\\flydata_example\\toy_example_search1.tsv',
                                      '..\\flydata_example\\toy_example_search2.tsv'])
    print(pepList)
    print(pepDict)
    protDict, pepDict = update_occurrence(protList, protDict, pepList, pepDict)
    print(protDict)
    print(pepDict)
    pepDict = update_all_unseen_peptides(protList, protDict, pepList, pepDict)
    print(pepDict)
    generate_outputFile(pepDict, '..\\output\\o1.tsv')

if __name__ == '__main__':
    main()

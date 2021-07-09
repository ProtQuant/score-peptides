# 30/06/2021 try to figure out errors in ProduceLabel3.py, and try speeding up the process

# no error in PL3.... I made sth wrong when testing...
# speed up: break up the protein before updating occurrences
# changes: one newly added method: find_all_possible_simple_peptides(protList)
#          one altered method: update_occurrence(allPossiblePepDict, protList, protDict_temp, pepList, pepDict)
#                                   also change the way to define "occur"


import pandas as pd
import time


def init_proteins(filepath):
    """
    initialize proteins and their index from .fasta file
    :param
        filepath: path of one .fasta file that describes all of the proteins believed to be in the species
    :return
        protList: a list of all proteins from the input file
                    for storing the indices of the proteins
        protDict_temp: a dictionary of all proteins from the input file
                    protein_formula: {'seenPeptides': a list of indices of the protein's peptides that can be seen
                                                      from the .tsv files,
                                      'leftProtein': the protein formula without seenPeptides (not empty only
                                                     when we know there are unseenPeptides from the protein),
                                      'unseenPeptides': a list of peptides formula broken up from leftProtein }
    """
    protList = []
    protDict = {}

    fl = open(filepath)
    protein = ''
    line = fl.readline().strip()
    while line:
        if '>' in line:
            if len(protein) > 0:  # when read in the first line, protein is empty
                protList.append(protein)
                protDict[protein] = {'seenPeptides': [], 'leftProtein': '',
                                     'unseenPeptides': []}  # will be updated later
                protein = ''  # '>' is the start of a new protein
        else:
            protein += line
        line = fl.readline().strip()
    # the last protein
    protList.append(protein)
    protDict[protein] = {'seenPeptides': [], 'leftProtein': '', 'unseenPeptides': []}

    print("There are " + str(len(protList)) + " proteins in total from .fasta files")

    # if sorting is needed, just sort protList
    return protList, protDict


def init_peptides(filepaths):
    """
    initialize peptides and their index from the .tsv files
    :param
        filepaths: a list of paths of the .tsv files of one specie
    :return:
        pepList: a list of all peptides from the input files without duplicates
                    for storing the indices of the peptides
        pepDict: a dictionary of all peptides from the input files
                    peptide_formula: {'total_num_ions': largest total_num_ions of the peptide from .tsv files,
                                      'occurrence': the number of the proteins that contain the peptide
                                      'proteinIndex': the list of the indices of the proteins that contain the peptide }
    """
    pepList = []
    pepDict = {}

    for filepath in filepaths:
        df = pd.read_csv(filepath, sep='\t')
        for i in range(len(df)):
            pep = df['peptide'][i]
            count = df['tot_num_ions'][i]
            if (pep not in pepDict):
                pepDict[pep] = {'total_num_ions': 0, 'occurrence': 0, 'proteinIndex': []}  # will be updated later
                pepDict[pep]['total_num_ions'] = count
                pepList.append(pep)
            elif (pepDict[pep]['total_num_ions'] < count):
                pepDict[pep]['total_num_ions'] = count

    print("There are " + str(len(pepList)) + " peptides in total from .tsv files")

    # if sorting is needed, just sort pepList
    return pepList, pepDict


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


def find_all_possible_simple_peptides(protList):
    """
    a newly added method
    break up all the proteins from .fasta file to get all possible simple peptides
    in order to speed up update_occurrence()
    :param
        protList: the first return value of init_proteins()
    :return:
        allPossiblePepDict: a dictionary storing the formula of all the possible simple peptides and the index of
                            a protein that generates them
                                simplePep_formula: [protein_index]
    """
    protIndex = 0
    allPossiblePepDict = {}
    for prot in protList:
        simplePepList = split_protein(prot)
        for pep in simplePepList:
            if pep not in allPossiblePepDict:
                allPossiblePepDict[pep] = []
            allPossiblePepDict[pep].append(protIndex)
        protIndex += 1

    print("There are " + str(len(allPossiblePepDict)) + "simple peptides in total from all the proteins")

    return allPossiblePepDict


def update_occurrence(allPossiblePepDict, protList, protDict, pepList, pepDict):
    """
    record the occurrence of peptides in proteins
    for each peptides in .tsv file (pepList), first search it in allPossiblePepDict, or search it through protDict_temp
    will update pepDict[pep]['occurrence'],  pepDict[pep]['proteinIndex'], protDict_temp[prot]['seenPeptides']
    :param
        returns of find_all_possible_simple_peptides(), init_proteins() and init_peptides()
    :return:
        updated protDict_temp and pepDict
    """
    pepIndex = 0
    for pep in pepList:
        # first search it in allPossiblePepDict
        if pep in allPossiblePepDict:
            protIndex = allPossiblePepDict[pep]
            pepDict[pep]['proteinIndex'] = protIndex
            pepDict[pep]['occurrence'] = len(protIndex)
            for i in protIndex:
                prot = protList[i]
                protDict[prot]['seenPeptides'].append(pepIndex)
        else:  # or search it through all the proteins
            occurrence = 0
            protIndex = 0
            for prot in protList:
                # if prot.find(pep) != -1: # two cases are correct
                flag = prot.find('K' + pep) != -1 or prot.find('R' + pep) != -1 \
                        or prot.startswith('M' + pep) or prot.startswith(pep)
                if flag:  # only the second case are correct
                    occurrence += 1
                    pepDict[pep]['proteinIndex'].append(protIndex)
                    protDict[prot]['seenPeptides'].append(pepIndex)
                protIndex += 1
            if occurrence == 0:
                print('can not find the ' + str(pepIndex) + 'th peptides: ' + pep + ' in the proteins')
            pepDict[pep]['occurrence'] = occurrence
        pepIndex += 1

    return protDict, pepDict


def find_unseen_peptides_from_one_protein(protein, protDict, pepList):
    """
    find unseen peptides of a present protein:
        use pepList and the indices in protDict_temp[prot]['seenPeptides'] to find the formula of the seen peptides
        update protDict_temp[prot]['leftProtein'] by removing seen peptides from original protein
        update protDict_temp[prot]['unseenPeptides'] by breaking up the leftProtein via split_protein()
    :param
        protein: a String formula of one original present protein
        protDict_temp: first return value of update_occurrence()
        pepList: second return value of update_occurrence()
    :return
        same as the return of split_protein()
    """
    leftProtein = protein

    seenPepFormula = []
    for pepIndex in protDict[protein]['seenPeptides']:
        seenPepFormula.append(pepList[pepIndex])
    seenPepFormula.sort(key=lambda x: len(x),
                        reverse=True)  # delete the longer peptides first  # when AAK and AAAK are in the same protein

    for pep in seenPepFormula:
        leftProtein = leftProtein.replace(pep, '')  # delete the seen peptides
    protDict[protein]['leftProtein'] = leftProtein
    if len(leftProtein) == 0:
        return []
    unseenPeptides = split_protein(leftProtein)
    protDict[protein]['unseenPeptides'] = unseenPeptides
    return unseenPeptides


def update_all_unseen_peptides(protList, protDict, pepList, pepDict):
    """
    add all unseen peptides into pepDict
    pepDict[pep]['proteinIndex'] of these unseen peptides, whose occurrence == 0, only records the present
        (not all) proteins, which are the proteins that produce the unseen peptides
    :param
        protList: first return value of init_proteins()
        protDict_temp: first return value of update_occurrence()
        pepList: first return value of init_peptides()
        pepDict: second return value of update_occurrence()
    :return:
        updated pepDict
    """
    for pep in pepList:
        if pepDict[pep]['occurrence'] == 1:
            protIndex = pepDict[pep]['proteinIndex'][0]
            protein = protList[protIndex]  # the protein that can be proved to be present
            if len(protDict[protein]['leftProtein']) != 0:  # unseen peptides in this protein has been found
                continue
            unseenPeptides = find_unseen_peptides_from_one_protein(protein, protDict, pepList)
            for pep in unseenPeptides:
                if pep not in pepDict:
                    pepDict[pep] = {'total_num_ions': 0, 'occurrence': 0, 'proteinIndex': []}
                pepDict[pep]['proteinIndex'].append(protIndex)
    return pepDict


def generate_outputFile(pepDict, filename):
    """
    generate the needed output file
    call generate_proteinIndex_file() if need to know the protein formulas
    """
    dataframe = pd.DataFrame({'peptide': pepDict.keys(),
                              'tot_num_ions': [pepDict[k]['total_num_ions'] for k in pepDict],
                              'occurrence': [pepDict[k]['occurrence'] for k in pepDict],
                              'proteinIndex': [pepDict[k]['proteinIndex'] for k in pepDict]})
    dataframe.to_csv(filename, index=False, sep='\t')


def generate_proteinIndex_file(protList, filename):
    """
    generate a flie that contains the proteins and their indices
    """
    dataframe = pd.DataFrame({'index': list(range(0, len(protList))),
                              'protein': protList})
    dataframe.to_csv(filename, index=False, sep='\t')


def main():
    """
    whole process (fruit-fly example)
        takes about 680 seconds, most of the time are spent at updating the occurrence of peptides in proteins
    """
    start = time.perf_counter()

    '''
    input files
    '''
    fastaFile = '..\\flydata_example\\uniprot-proteome_UP000000803.fasta'
    tsvFiles = ['..\\flydata_example\\QE752_MSQ988_20170315_BenHopkins_R1V8.tsv',
                '..\\flydata_example\\QE752_MSQ988_20170315_BenHopkins_R2M8.tsv',
                '..\\flydata_example\\QE752_MSQ988_20170315_BenHopkins_R4M2.tsv']
    '''
    generated files
    '''
    protIndexFile = '..\\flydata_example\\LP4_try_fruit_fly_protIndex.tsv'  # proteins and their index
    outputFile = '..\\flydata_example\\LP4_try_fruit_fly_output.tsv'  # the needed output

    protList, protDict = init_proteins(fastaFile)
    pepList, pepDict = init_peptides(tsvFiles)
    t1 = time.perf_counter()
    print('finish initializing....\t' + str(t1 - start))  # record time needed for each part
    allPossiblePepDict = find_all_possible_simple_peptides(protList)
    t2 = time.perf_counter()
    print('finish breaking up all proteins....\t' + str(t2 - t1))
    protDict, pepDict = update_occurrence(allPossiblePepDict, protList, protDict, pepList, pepDict)
    t3 = time.perf_counter()
    print('finish updating occurrence....\t' + str(t3 - t2))  # 99% of total time
    pepDict = update_all_unseen_peptides(protList, protDict, pepList, pepDict)
    t4 = time.perf_counter()
    print('finish adding unseen peptides....\t' + str(t4 - t3))

    generate_proteinIndex_file(protList, protIndexFile)
    t5 = time.perf_counter()
    print('finish generating protein index file....\t' + str(t5 - t4))
    generate_outputFile(pepDict, outputFile)
    t6 = time.perf_counter()
    print('finish generating output file....\t' + str(t6 - t5))

    end = time.perf_counter()
    print("")
    print('finish running fruit-fly example, total processing time: \t' + str(end - start))

    print("")


if __name__ == '__main__':
    main()

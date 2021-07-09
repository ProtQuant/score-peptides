# 2021/6/25 reformat code from ProduceLabel2.py and some changes

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

    # if sorting is needed, just sort pepList
    return pepList, pepDict


def update_occurrence(protList, protDict, pepList, pepDict):
    """
    record the occurrence of peptides in proteins
    will update pepDict[pep]['occurrence'],  pepDict[pep]['proteinIndex'], protDict_temp[prot]['seenPeptides']
    :param
        returns of init_proteins() and init_peptides()
    :return:
        updated protDict_temp and pepDict
    """
    pepIndex = 0
    for pep in pepList:
        occurrence = 0
        protIndex = 0
        for prot in protList:
            # if prot.find(pep) != -1:
            if prot.find('K' + pep) != -1 or prot.find('R' + pep) != -1 or \
                    prot.startswith('M' + pep) or prot.startswith(pep):
                occurrence += 1
                pepDict[pep]['proteinIndex'].append(protIndex)
                protDict[prot]['seenPeptides'].append(pepIndex)
            protIndex += 1
        if occurrence == 0:
            print('can not find the ' + str(pepIndex) + 'th peptides: ' + pep + ' in the proteins')
        pepDict[pep]['occurrence'] = occurrence
        pepIndex += 1

    return protDict, pepDict


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
            protein = protList[protIndex]  # the protein that can be prove to be present
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
        takes about 400 seconds, most of the time are spent at updating the occurrence of peptides in proteins
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
    protIndexFile = '..\\flydata_example\\LP3_try_fruit_fly_protIndex.tsv'  # proteins and their index
    outputFile = '..\\flydata_example\\LP3_try_fruit_fly_output.tsv'  # the needed output

    protList, protDict = init_proteins(fastaFile)
    pepList, pepDict = init_peptides(tsvFiles)
    t1 = time.perf_counter()
    print('finish initializing....\t' + str(t1 - start))  # record time needed for each part
    protDict, pepDict = update_occurrence(protList, protDict, pepList, pepDict)
    t2 = time.perf_counter()
    print('finish updating occurrence....\t' + str(t2 - t1))  # 99% of total time
    pepDict = update_all_unseen_peptides(protList, protDict, pepList, pepDict)
    t3 = time.perf_counter()
    print('finish adding unseen peptides....\t' + str(t3 - t2))

    generate_proteinIndex_file(protList, protIndexFile)
    t4 = time.perf_counter()
    print('finish generating protein index file....\t' + str(t4 - t3))
    generate_outputFile(pepDict, outputFile)
    t5 = time.perf_counter()
    print('finish generating output file....\t' + str(t5 - t4))

    end = time.perf_counter()
    print("")
    print('finish running fruit-fly example, total processing time: \t' + str(end - start))

    print("")
    """
    Test cases
    """
    print("************************** Testing **************************")
    '''
    init_proteins(filepath)
        1. read in all proteins 
        2. read in one protein of many lines 
        3. ignore the content in the line that contains '>'
    '''
    print('================== Test init_proteins ==================')
    fastaFile = '.\\Testcases\\test.fasta'
    protList, protDict = init_proteins(fastaFile)
    print('protList: \t' + str(protList))
    print('protDict_temp: \t' + str(protDict))
    '''
    init_peptides(filepaths)
        1. read in all peptides from different files 
        2. no duplicate 
        3. record the largest 'total_num_ions'
    '''
    print('================== Test init_peptides ==================')
    tsvFiles = ['.\\Testcases\\test_search1.tsv', '.\\Testcases\\test_search2.tsv']
    pepList, pepDict = init_peptides(tsvFiles)
    print('pepList: \t' + str(pepList))
    print('pepDict: \t' + str(pepDict))
    '''
    update_occurrence(protList, protDict_temp, pepList, pepDict)
        1. the number of occurrence of one protein in the proteins: pepDict[pep]['occurrence']
        2. index of the proteins that contain one particular peptides: pepDict[pep]['proteinIndex']
        3. index of peptides in one chosen protein: protDict_temp[prot]['seenPeptides']
        4. print error when a peptide can not be found in the proteins
    '''
    print('================== Test update_occurrence ==================')
    protDict, pepDict = update_occurrence(protList, protDict, pepList, pepDict)
    print('protDict_temp: \t' + str(protDict))
    print('pepDict: \t' + str(pepDict))
    '''
    split_protein(protein)
        1. break after K or R without P following it
        2. not missing the last peptides
    '''
    print('================== Test split_protein ==================')
    print(split_protein('AAKAAR'))
    print(split_protein('AAKAARAAA'))
    print(split_protein('AAKPAARAAARPAAARAAAKKKPA'))
    '''
    find_unseen_peptides_from_one_protein(protein, protDict_temp, pepList)
        1. correctly get the leftProtein by removing the seen peptides from original protein
        2. correctly get the unseen peptides
        3. proteins that do not have unseen peptides
        4. when one seen peptide contains another
        
        will change pepDict, this part need to be commented before testing the next part
    '''
    # print('================== Test find_unseen_peptides_from_one_protein ==================')
    # print(find_unseen_peptides_from_one_protein('AAMRAADRAAUKAAKAAWRAAU', protDict_temp, pepList))
    # print(protDict_temp['AAMRAADRAAUKAAKAAWRAAU']['leftProtein'])
    # print(find_unseen_peptides_from_one_protein('AARAAK', protDict_temp, pepList))
    # print(protDict_temp['AARAAK']['leftProtein'])
    # print(find_unseen_peptides_from_one_protein('AAAKAAU', protDict_temp, pepList))
    # print(protDict_temp['AAAKAAU']['leftProtein'])
    '''
    update_all_unseen_peptides(protList, protDict_temp, pepList, pepDict)
        1. do not find unseen peptides from one protein more than once
    '''
    print('================== Test update_all_unseen_peptides ==================')
    print(update_all_unseen_peptides(protList, protDict, pepList, pepDict))

    protIndexFile = '.\\Testcases\\protIndex.tsv'
    generate_proteinIndex_file(protList, protIndexFile)

    outputFile = '.\\Testcases\\output.tsv'
    generate_outputFile(pepDict, outputFile)


if __name__ == '__main__':
    main()

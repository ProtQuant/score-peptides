#  produce score each ten files

# 2021/7/21
# need to install openpyxl
# BBK is in MBBKAKKCCKDDK, but if MBBK is also find,  how to calculate prot_total. Problem persists.
# scoring: [20, 20] -- > 20; [20, 20, 50] --> 50. Solved

import copy
import heapq
import os
import numpy
import pandas as pd
import time

base = '../saved data/'+ os.path.basename(__file__).split('.')[0]

def init_proteins(in_filepath, generateFile=False, out_filePath= base + '/init_proteins.csv'):
    """
    * initialize proteins and their index from .fasta file, note that duplicate protein may exist in the file
    * will create protList and protDict
    :param
        in_filepath: string, path of one .fasta file that describes all of the proteins believed to be in the species
        generateFile: boolean, will output a file containing the items in protDict if true
        out_filePath: string, path of the above file, needs to be located in an existing directory
    :return
        protList: a list of all proteins from the input file
                    for storing the indices of the proteins
        protDict: a dictionary of all proteins from the input file
                    protein_formula: {'pepIndex': a list of indices of the protein's peptides that can be seen
                                                      from the .csv files,
                                      'count': int, how many peptides from one .csv file are in this protein
                                      'leftProtein': no use for now, was:
                                                     string, the protein formula without seenPeptides (not empty only
                                                     when we know there are unseenPeptides from the protein),
                                      'unseenPeptides': a list of peptides formula broken up from leftProtein,
                                      'list_intensity': a list of intensity value of the seen peptides
                                      'prot_total':int, sum of the largest 3 intensity value
                                                        set to -1 if list_intensity is empty}
    """
    protList = []
    protDict = {}

    fl = open(in_filepath)
    protein = ''
    line = fl.readline().strip()
    while line:
        if '>' in line:
            if len(protein) > 0:  # when read in the first line, protein is empty
                if protein not in protDict:
                    protList.append(protein)
                    protDict[protein] = {'pepIndex': [], 'count': 0, 'leftProtein': '', 'unseenPeptides': [],
                                         'list_intensity': [], 'prot_total': -1}
                protein = ''  # '>' is the start of a new protein
        else:
            protein += line
        line = fl.readline().strip()
    # the last protein
    if protein not in protDict:
        protList.append(protein)
        protDict[protein] = {'pepIndex': [], 'count': 0, 'leftProtein': '', 'unseenPeptides': [],
                             'list_intensity': [], 'prot_total': -1}

    print("There are " + str(len(protList)) + " proteins in total from .fasta files")  # 21180

    if generateFile:
        protDf = pd.DataFrame({'protein': protDict.keys(),
                               'pepIndex': [protDict[k]['pepIndex'] for k in protDict],
                               'count': [protDict[k]['count'] for k in protDict],
                               'leftProtein': [protDict[k]['leftProtein'] for k in protDict],
                               'unseenPeptides': [protDict[k]['unseenPeptides'] for k in protDict],
                               'list_intensity': [protDict[k]['list_intensity'] for k in protDict],
                               'prot_total': [protDict[k]['prot_total'] for k in protDict]})
        protDf.to_csv(out_filePath, index=True, sep=',')

    return protList, protDict


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
    for amino_acid in protein:  # for each letter in the string
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


def find_all_possible_simple_peptides(protList, generateFile=False,
                                      out_filePath= base + '/simple_peptides.csv'):
    """
    break up all the proteins from .fasta file to get all possible simple peptides
    will create allPossiblePepDict
    :param
        protList: a list containing all the protein formula, should be the first return value of init_proteins()
        generateFile: boolean, will output a file containing the items in allPossiblePepDict if true
        out_filePath: string, path of the above file, need to be located in an existing directory
    :return:
        allPossiblePepDict: a dictionary storing the formula of all the possible simple peptides and the indices of
                            proteins that contain them
                                simplePep_formula: [protein_indices]
        protpepDictList: a list storing the simple peptides in each protein from protList
                                [simplePepDcit] = [{simplePep_formula: protein_index}]
    """
    protIndex = 0
    allPossiblePepDict = {}  # simple pep dict for all proteins
    protpepDictList = []

    for prot in protList:
        simplePepDcit = {}  # simple pep dict for one protein
        simplePepList = split_protein(prot)
        for pep in simplePepList:
            simplePepDcit[pep] = protIndex
            if pep not in allPossiblePepDict:
                allPossiblePepDict[pep] = []
            allPossiblePepDict[pep].append(protIndex)
        protIndex += 1
        protpepDictList.append(simplePepDcit)

    print("There are " + str(len(allPossiblePepDict)) + " simple peptides from all the proteins")
    # print(len(protpepDictList))

    if generateFile:
        sPepDf = pd.DataFrame({'sPeptide': allPossiblePepDict.keys(),
                               'protIndex': allPossiblePepDict.values()})
        sPepDf.to_csv(out_filePath, index=True, sep=',')

    return allPossiblePepDict, protpepDictList


def init_peptides_from_one_file(filepath, colname=['peptide', 'intensity'],
                                generateFile=False, out_filePath=base+'/init_peptides_from_one_file.csv'):
    """
    initialize peptides and their index from the one of the .csv files
    will create pepList and pepDict
    will update pepDict['intensity']
    :param
        filepath: string, path of a .csv file storing information of peptides
        colname: list of string, the needed column names,
                    the first element should be the name of the column of peptides' formula
                    the second element should be the name of the column of values for scoring
        generateFile: boolean, will output a file containing items in pepDict if true
        out_filePath: string, path of the above file, needs to be located in an existing directory
    :return:
        pepList: a list of all peptides from the input files without duplicates
                    for storing the indices of the peptides
        pepDict: a dictionary of all peptides from the input files
                    peptide_formula: {'intensity': largest intensity of the peptide from one .csv file,
                                      'occurrence': the number of the proteins that contain the peptide,
                                      'protIndex': a list of the indices of the proteins that contain the peptide, }
    """
    pepList = []
    pepDict = {}

    # read in the file
    if filepath.endswith(".csv"):
        df = pd.read_csv(filepath, sep=',', usecols=colname)  # dtype?
    elif filepath.endswith(".tsv"):
        df = pd.read_csv(filepath, sep='\t', usecols=colname)

    for i in range(len(df)):
        if numpy.isnan(df[colname[1]][i]):
            continue
        pep = df[colname[0]][i]
        intensity = df[colname[1]][i]
        if pep not in pepDict:
            pepDict[pep] = {'intensity': 0, 'occurrence': 0, 'protIndex': []}
            pepDict[pep]['intensity'] = intensity
            pepList.append(pep)
        elif pepDict[pep]['intensity'] < intensity:
            pepDict[pep]['intensity'] = intensity

    print("\t There are " + str(len(pepList)) + " peptides in this .csv file")

    if generateFile:
        pepDf = pd.DataFrame({'peptide': pepDict.keys(),
                              'intensity': [pepDict[k]['intensity'] for k in pepDict],
                              'occurrence': [pepDict[k]['occurrence'] for k in pepDict],
                              'protIndex': [pepDict[k]['protIndex'] for k in pepDict]})
        pepDf.to_csv(out_filePath, index=True, sep=',')

    return pepList, pepDict


def update_occurrence_of_one_file(allPossiblePepDict, protpepDictList, protList, protDict, pepList, pepDict,
                                  generateFile=False,
                                  out_filePath_pep=base+'/pep_update_occurrence_of_one_file.csv',
                                  out_filePath_prot=base+'/prot_update_occurrence_of_one_file.csv'):
    """
    record the occurrence of peptides in proteins, find prot_total for each peptides
    delete the seen peptides in protains in protpepDictList_temp
    for each peptides in .csv file (pepList), first search it among simple peptides, or search it in each protein
    will update pepDict[pep]['occurrence'],  pepDict[pep]['protIndex'],
                protDict_temp[prot]['pepIndex'], protDict_temp[prot]['count'], protDict_temp[prot]['list_intensity'],
                    protDict[prot]['prot_total']
                protpepDictList_temp
    :param
        return values of find_all_possible_simple_peptides(), init_proteins() and init_peptides_from_one_file()
    :return:
        updated protpepDictList_temp: deleted peptides from pepList
        updated protDict_temp: has same column with protDict,
        updated pepDict,
        ignored: int, the number of proteins that has only one peptide appears in this .csv file
        valued: int, the number of proteins that has more than one peptides appears in this .csv file
    """
    # the original protDict should not change when processing different files
    protDict_temp = copy.deepcopy(protDict)
    protpepDictList_temp = copy.deepcopy(protpepDictList)

    # update occurrence
    pepIndex = 0
    for pep in pepList:  # pepList contains all the peptides from .tsv/.csv files

        isSimple = False
        # first search it in allPossiblePepDict
        if pep in allPossiblePepDict:  # allPossiblePepDict contains all the simple peptides
            isSimple = True
            protIndex = allPossiblePepDict[pep]
            pepDict[pep]['protIndex'] = protIndex
            pepDict[pep]['occurrence'] = len(protIndex)
            for i in protIndex:
                prot = protList[i]
                protDict_temp[prot]['pepIndex'].append(pepIndex)
                protDict_temp[prot]['list_intensity'].append(pepDict[pep]['intensity'])
                try:
                    del protpepDictList_temp[i][pep]  # try..except..
                except Exception:
                    pass

        # or search it through all the proteins
        else:
            occurrence = 0
            protIndex = 0
            for prot in protList:
                # the next letter should not be a P
                if prot.find(
                        pep + 'P') == -1:  # ignore for now the case that a peptide may appear twice in a protein, and one can be counted
                    # if prot.startswith(pep) or prot.startswith('M' + pep) or prot.find('K' + pep) or prot.find(
                    #         'R' + pep):  # will raise memoryError??
                    #     occurrence += 1
                    #     pepDict[pep]['protIndex'].append(protIndex)
                    #     protDict_temp[prot]['pepIndex'].append(pepIndex)
                    #     protDict_temp[prot]['list_intensity'].append(pepDict[pep]['intensity'])

                    # may locate at the head of a protein
                    if prot.startswith(pep):
                        occurrence += 1
                        pepDict[pep]['protIndex'].append(protIndex)
                        protDict_temp[prot]['pepIndex'].append(pepIndex)
                        protDict_temp[prot]['list_intensity'].append(pepDict[pep]['intensity'])
                        peps = split_protein(pep)
                        for peptide in peps:
                            try:
                                del protpepDictList_temp[protIndex][pep]
                            except Exception:
                                pass
                    # may locate at the head of a protein missing an M
                    # problems appear when AAK, MAAK both exist
                    elif prot.startswith('M' + pep):
                        # print("M problem")
                        occurrence += 1
                        pepDict[pep]['protIndex'].append(protIndex)
                        protDict_temp[prot]['pepIndex'].append(pepIndex)
                        protDict_temp[prot]['list_intensity'].append(pepDict[pep]['intensity'])
                        peps = split_protein('M' + pep)
                        for peptide in peps:
                            try:
                                del protpepDictList_temp[protIndex][pep]
                            except Exception:
                                pass
                    # may locate in the protein following a K
                    elif prot.find('K' + pep) != -1:
                        occurrence += 1
                        pepDict[pep]['protIndex'].append(protIndex)
                        protDict_temp[prot]['pepIndex'].append(pepIndex)
                        protDict_temp[prot]['list_intensity'].append(pepDict[pep]['intensity'])
                        peps = split_protein(pep)
                        for peptide in peps:
                            try:
                                del protpepDictList_temp[protIndex][pep]
                            except Exception:
                                pass
                    # may locate in the protein following an R
                    elif prot.find('R' + pep) != -1:
                        occurrence += 1
                        pepDict[pep]['protIndex'].append(protIndex)
                        protDict_temp[prot]['pepIndex'].append(pepIndex)
                        protDict_temp[prot]['list_intensity'].append(pepDict[pep]['intensity'])
                        peps = split_protein(pep)
                        for peptide in peps:
                            try:
                                del protpepDictList_temp[protIndex][pep]
                            except Exception:
                                pass
                protIndex += 1
            if occurrence == 0:
                print('\t can not find the ' + str(pepIndex) + 'th peptides: ' + pep + ' in the proteins')
            pepDict[pep]['occurrence'] = occurrence
        pepIndex += 1

    # update count and prot_total
    ignored = 0
    valued = 0
    for prot in protDict_temp:
        protDict_temp[prot]['count'] = len(protDict_temp[prot]['pepIndex'])
        count = protDict_temp[prot]['count']
        if count == 1:
            ignored += 1
            protDict_temp[prot]['prot_total'] = sum(protDict_temp[prot]['list_intensity'])
        elif count == 2:
            protDict_temp[prot]['prot_total'] = sum(protDict_temp[prot]['list_intensity'])
            valued += 1
        elif count >= 3:
            protDict_temp[prot]['prot_total'] = sum(heapq.nlargest(3, protDict_temp[prot]['list_intensity']))
            valued += 1

    if generateFile:
        pepDf = pd.DataFrame({'peptide': pepDict.keys(),
                              'intensity': [pepDict[k]['intensity'] for k in pepDict],
                              'occurrence': [pepDict[k]['occurrence'] for k in pepDict],
                              'protIndex': [pepDict[k]['protIndex'] for k in pepDict]})
        pepDf.to_csv(out_filePath_pep, index=True, sep=',')
        protDf = pd.DataFrame({'protein': protDict_temp.keys(),
                               'pepIndex': [protDict_temp[k]['pepIndex'] for k in protDict_temp],
                               'count': [protDict_temp[k]['count'] for k in protDict_temp],
                               'leftProtein': [protDict_temp[k]['leftProtein'] for k in protDict_temp],
                               'unseenPeptides': [protDict_temp[k]['unseenPeptides'] for k in protDict_temp],
                               'list_intensity': [protDict_temp[k]['list_intensity'] for k in protDict_temp],
                               'prot_total': [protDict_temp[k]['prot_total'] for k in protDict_temp]})
        protDf.to_csv(out_filePath_prot, index=True, sep=',')

    return protpepDictList_temp, protDict_temp, pepDict, ignored, valued


def update_unseen_peptides_of_one_file(protList, protDict_temp, pepList, pepDict, allPossiblePepDict, protpepDictList_temp,
                                       generateFile=False,
                                       out_filePath_prot=base+'/prot_update_unseen_peptides_of_one_file.csv',
                                       out_filePath_unsPep=base+'/unsPep_update_unseen_peptides_of_one_file.csv'
                                       ):
    """
    will create unseenPepDict for unseen peptides
    will update unseenPepDict [pep]['protIndex'] , which are the indices of proteins that produce the unseen peptides
    will update protDict_temp[prot]['unseenPeptides'] when calling find_unseen_peptides_from_one_protein()
    :param
        protList: first return value of init_proteins()
        protDict_temp: second return value of update_occurrence_of_one_file()
        pepList: first return value of init_peptides()
        pepDict: second return value of update_occurrence()
        allPossiblePepDict:
        protpepDictList_temp: first return value of update_occurrence_of_one_file()
        generateFile: boolean, will output two files containing items in protDict_temp and unseenPepDict respectively if true
        out_filePath_prot: string, path for the above first file, needs to be located in an existing directory
        out_filePath_unsPep: string, path for the above first file, needs to be located in an existing directory
    :return:
        updated protDict_temp
        unseenPepDict: a dictionary for the unseen peptides in one .csv file
            peptide_formula: {'score': 0, 'protIndex': list of indices of the proteins that contain the peptide}
        unsPepDf: pandas dataframe, containing the unseen peptides according to this .csv file
                    columns = ['peptide', 'score']
    """
    unseenPepDict = {}
    for pep in pepList:
        if pepDict[pep]['occurrence'] == 1:
            protIndex = pepDict[pep]['protIndex'][0]
            protein = protList[protIndex]  # the protein that can be proved to be present
            if len(protDict_temp[protein]['unseenPeptides']) != 0:
                continue
            unseenPeptides = protpepDictList_temp[protIndex].keys()
            protDict_temp[protein]['unseenPeptides'] = unseenPeptides
            for pep in unseenPeptides:
                if pep not in unseenPepDict:
                    unseenPepDict[pep] = {'score': 0, 'protIndex': []}
                unseenPepDict[pep]['protIndex'].append(protIndex)
                # if pep not in allPossiblePepDict:
                #     print('\t invalid unseen peptides: ' + pep)

    if generateFile:
        protDf = pd.DataFrame({'protein': protDict_temp.keys(),
                               'pepIndex': [protDict_temp[k]['pepIndex'] for k in protDict_temp],
                               'count': [protDict_temp[k]['count'] for k in protDict_temp],
                               'leftProtein': [protDict_temp[k]['leftProtein'] for k in protDict_temp],
                               'unseenPeptides': [protDict_temp[k]['unseenPeptides'] for k in protDict_temp],
                               'list_intensity': [protDict_temp[k]['list_intensity'] for k in protDict_temp],
                               'prot_total': [protDict_temp[k]['prot_total'] for k in protDict_temp]})
        protDf.to_csv(out_filePath_prot, index=True, sep=',')
        unsPepDf = pd.DataFrame({'peptide': unseenPepDict.keys(),
                                 'score': [unseenPepDict[k]['score'] for k in unseenPepDict],
                                 'protIndex': [unseenPepDict[k]['protIndex'] for k in
                                               unseenPepDict]})  # may duplicate due to dif file
        unsPepDf.to_csv(out_filePath_unsPep, index=True, sep=',')

    unsPepDf = pd.DataFrame({'peptide': unseenPepDict.keys(),
                             'score': [unseenPepDict[k]['score'] for k in unseenPepDict]})

    return protDict_temp, unseenPepDict, unsPepDf


def calculate_score_of_one_file(pepDict, protDict_temp, generateFile=False,
                                out_filePath=base+'/7_calculate_score_of_one_file.csv'):
    """
        One peptide may have many scores according to pepDict[pep]['protIndex']
        Score the peptide as 20 if it has an intensity value of less than 20% of the protein total
            (The corresponding protein has only this one peptide appear in this .csv file)
        Score the peptide as 100 if it has an intensity value of more than 20% of the protein total
        Score the peptide as 50 if it has an intensity value of less than 20% of the protein total
    :param
        pepDict: second return value of update_occurrence_of_one_file()
        protDict_temp: first return value of update_occurrence_of_one_file() or update_unseen_peptides_of_file()
        generateFile: boolean, will generate file containing all the relevant infomation of peptides
        out_filePath: string, path for the above file, needs to be located in an existing directory
    :return:
        pepDf: pandas dataframe, contains all the peptides from the .csv file
                column = ['peptide', 'intensity', 'occurrence', 'protIndex', 'score_0', 'score']
        pepDf_score: pandas dataframe, only contains peptides scored 50 or 100
                column = ['peptide', 'score']
        pepDf_ignore: pandas dataframe, only contains peptides scored 20
                column = ['peptide', 'score']
    """
    pepDf = pd.DataFrame({'peptide': pepDict.keys(),
                          'intensity': [pepDict[k]['intensity'] for k in pepDict],
                          'occurrence': [pepDict[k]['occurrence'] for k in pepDict],
                          'protIndex': [pepDict[k]['protIndex'] for k in pepDict]})

    pepDf = pepDf.explode('protIndex').reset_index()  # NaN if protIndex is []

    protDf = pd.DataFrame({'prot_total': [protDict_temp[k]['prot_total'] for k in protDict_temp]},
                          index=list(range(0, len(protDict_temp))))

    # add 'prot_total' to pepDf according to 'protIndex'
    pepDf = pepDf.join(protDf, on="protIndex")

    def calculate_score(x):
        if x['prot_total'] > 0:
            s = x['intensity'] / x['prot_total']
            return s

    def define_level(x):
        s = x['score_0']
        if s == 1:
            return 20
        elif s > 0.2:
            return 100
        elif s <= 0.2:
            return 50  # NaN is still NaN

    pepDf['score_0'] = pepDf.apply(lambda x: calculate_score(x), axis=1)
    pepDf['score'] = pepDf.apply(lambda x: define_level(x), axis=1)

    if generateFile:
        pepDf.to_csv(out_filePath, index=True, sep=',')

    pepDf_have_score = pd.DataFrame(pepDf, columns=['peptide', 'score'])
    pepDf_have_score = pepDf_have_score.dropna().reset_index(drop=True)
    pepDf_score = pepDf_have_score.drop(pepDf_have_score[pepDf_have_score['score'] == 20].index)
    pepDf_ignore = pepDf_have_score.drop(pepDf_have_score[pepDf_have_score['score'] != 20].index)

    return pepDf, pepDf_score, pepDf_ignore


def main():
    start = time.perf_counter()

    """
    Two input values:
        fastaFile: string, path of one .fasta file that contains protein formulas
        folderPath: string, path of one folder that contains the many files of peptides' infomation
    """
    fastaFile = '../uniprot-proteome_UP000000803.fasta'
    folderPath = '../flydata_example3'

    # fastaFile = '../Testcases4/T5/test.fasta'
    # folderPath = '../Testcases4/T5/flyquant'
    """
    Two debugging options:
        proteinInfo: boolean, will generate one excel file with two sheets in proteinInfo folder if set as true
                        one sheet for protein index and protein formula
                        one sheet for all simple peptides and their parent protein indices 
        peptideInfo: boolean, will generate one excel file with 3 sheets in peptideInfo folder for each input .csv file 
                     if set as true
                        one sheet for protein ['',]
                        one sheet for scored peptides []
                        one sheet for unseen peptides []
    Remember to close the generated files before starting the programme again.
    """
    proteinInfo = False
    peptideInfo = False

    # proteinInfo = True
    # peptideInfo = True

    """
    Main Process
    """
    # create default output folder
    base = '../saved data/' + os.path.basename(__file__).split('.')[0]
    if not os.path.exists(base):
        os.makedirs(base)

    # initialize proteins and simple peptides
    generateFile = False
    protList, protDict = init_proteins(fastaFile, generateFile)
    allPossiblePepDict, protpepDictList = find_all_possible_simple_peptides(protList)

    # output debugging info
    if proteinInfo:
        proteinInfoFolder = base + "/proteinInfo"
        if not os.path.exists(proteinInfoFolder):
            os.makedirs(proteinInfoFolder)
        protDf = pd.DataFrame({'protein': protDict.keys()})
        sPepDf = pd.DataFrame({'sPeptide': allPossiblePepDict.keys(),
                               'protIndex': allPossiblePepDict.values()})

        writer = pd.ExcelWriter(proteinInfoFolder + '/proteinInfo.xlsx')
        protDf.to_excel(writer, sheet_name='protein index', index=True)
        sPepDf.to_excel(writer, sheet_name='simple peptides', index=True)
        writer.save()

    protTime = time.perf_counter()
    print('Finish initializing proteins, time: ' + str(protTime))

    # 3 dataframe to store the info from each .csv file
    fileInfoDf = pd.DataFrame(columns=['fileName', '#ignoredProtein', '#valuedProtein']) # actually we didn't ignore them
    all_scoreDf = pd.DataFrame(columns=['peptide', 'score'])
    all_ignoreDf = pd.DataFrame(columns=['peptide', 'score'])
    all_unsPepDf = pd.DataFrame(columns=['peptide', 'score'])

    # process the .csv files one by one
    filelist = os.listdir(folderPath)
    fileIndex = 0
    for file in filelist:  # empty file may lead to problem
        if file.endswith('.csv') or file.endswith('.tsv'):
            pepTime1 = time.perf_counter()
            fileIndex += 1
            print('processing file ' +str(fileIndex)+' : '+ file)
            filePath = folderPath + '/' + file
            pepList, pepDict = init_peptides_from_one_file(filePath)
            protpepDictList_temp, protDict_temp, pepDict, ignored, valued = update_occurrence_of_one_file(
                allPossiblePepDict, protpepDictList, protList, protDict, pepList, pepDict)
            protDict_temp, unseenPepDict, unsPepDf = update_unseen_peptides_of_one_file(protList, protDict_temp,
                                                                                        pepList, pepDict,
                                                                                        allPossiblePepDict,
                                                                                        protpepDictList_temp)
            pepDf, pepDf_score, pepDf_ignore = calculate_score_of_one_file(pepDict, protDict_temp)

            all_scoreDf = all_scoreDf.append(pepDf_score, ignore_index=True)
            all_ignoreDf = all_ignoreDf.append(pepDf_ignore, ignore_index=True)
            all_unsPepDf = all_unsPepDf.append(unsPepDf, ignore_index=True)
            fileInfoDf = fileInfoDf.append([{'fileName': file, '#ignoredProtein': ignored, '#valuedProtein': valued}],
                                           ignore_index=True)

            # output debugging info
            if peptideInfo:
                protDf = pd.DataFrame({'protein': protDict_temp.keys(),
                                       'pepIndex': [protDict_temp[k]['pepIndex'] for k in protDict_temp],
                                       'count': [protDict_temp[k]['count'] for k in protDict_temp],
                                       # 'leftProtein': [protDict_temp[k]['leftProtein'] for k in protDict_temp],
                                       'unseenPeptides': [protDict_temp[k]['unseenPeptides'] for k in protDict_temp],
                                       'list_intensity': [protDict_temp[k]['list_intensity'] for k in protDict_temp],
                                       'prot_total': [protDict_temp[k]['prot_total'] for k in protDict_temp]})
                unseenPepDf = pd.DataFrame({'peptide': unseenPepDict.keys(),
                                            'score': [unseenPepDict[k]['score'] for k in unseenPepDict],
                                            'protIndex': [unseenPepDict[k]['protIndex'] for k in
                                                          unseenPepDict]})
                peptideInfoFolder = base + "/peptideInfo"
                if not os.path.exists(peptideInfoFolder):
                    os.makedirs(peptideInfoFolder)
                peptideInfoFile = base + "/peptideInfo/" + str(fileIndex) + "_Debug " + file.replace('.csv',
                                                                                                              '')
                writer = pd.ExcelWriter(peptideInfoFile + '.xlsx')
                protDf.to_excel(writer, sheet_name='protein information', index=True)
                pepDf.to_excel(writer, sheet_name='peptide score', index=True)
                unseenPepDf.to_excel(writer, sheet_name='unseen peptide', index=True)
                writer.save()

            # uncomment below if you want to output accumulative data after processing every 10 .csv files
            # if fileIndex % 10 == 0:
            #     print('* create output data for '+ str(fileIndex) + 'files')
            #     all_scoreDf1 = all_scoreDf.groupby(['peptide'])['score'].agg(
            #         lambda x: x.value_counts().index[0]).reset_index(
            #         drop=False)
            #     all_ignoreDf1 = all_ignoreDf.drop_duplicates().reset_index(drop=True)
            #     all_ignoreDf1 = all_ignoreDf1[~ all_ignoreDf1['peptide'].isin(all_scoreDf1['peptide'])]
            #     all_unsPepDf1 = all_unsPepDf.drop_duplicates().reset_index(drop=True)
            #     all_unsPepDf1 = all_unsPepDf1[~ all_unsPepDf1['peptide'].isin(all_ignoreDf1['peptide'])]
            #     all_unsPepDf1 = all_unsPepDf1[~ all_unsPepDf1['peptide'].isin(all_scoreDf1['peptide'])]
            #     finalScoreDf = all_scoreDf1.append(all_ignoreDf1, ignore_index=True).append(all_unsPepDf1,
            #                                                                                 ignore_index=True)
            #     outputFolder = folderPath + '/../output'
            #     if not os.path.exists(outputFolder):
            #         os.makedirs(outputFolder)
            #     finalScoreDf.to_csv(outputFolder + '/' + str(fileIndex) + '_score.csv', index=True, sep=',')
            #     fileInfoDf.to_csv(outputFolder + '/' + str(fileIndex) + '_ignore.csv', index=True, sep=',')

            pepTime2 = time.perf_counter()
            print('\t Time needed: ' + str(pepTime2-pepTime1))

    # print(all_scoreDf)
    # print(all_scoreDf.value_counts())

    # choose the most common score, will choose the larger score if the quantities of each value are the same
    all_scoreDf = all_scoreDf.groupby(['peptide'])['score'].agg(
        lambda x: x.value_counts().index[0]).reset_index(
        drop=False)  # !! wrong if add x.value_counts(sort=False)  # reset index to make it still a dataframe rather than a series
    # print(all_scoreDf)

    # delete duplicate unseen peptides and delete the ones appearing in all_scoreDf
    all_ignoreDf = all_ignoreDf.drop_duplicates().reset_index(drop=True)
    all_ignoreDf = all_ignoreDf[~ all_ignoreDf['peptide'].isin(all_scoreDf['peptide'])]

    # delete duplicate unseen peptides and delete the ones appearing in all_scoreDf and all_ignoreDf
    all_unsPepDf = all_unsPepDf.drop_duplicates().reset_index(drop=True)
    all_unsPepDf = all_unsPepDf[~ all_unsPepDf['peptide'].isin(all_ignoreDf['peptide'])]
    all_unsPepDf = all_unsPepDf[~ all_unsPepDf['peptide'].isin(all_scoreDf['peptide'])]

    # concatenating scored peptides and unseen peptides
    finalScoreDf = all_scoreDf.append(all_ignoreDf, ignore_index=True).append(all_unsPepDf, ignore_index=True)

    outputFolder = base + '/output'
    if not os.path.exists(outputFolder):
        os.makedirs(outputFolder)
    finalScoreDf.to_csv(outputFolder + '/score.csv', index=True, sep=',')
    fileInfoDf.to_csv(outputFolder + '/ignore.csv', index=True, sep=',')

    end = time.perf_counter()
    print('Finish all, total processing time: \t' + str(end - start))


if __name__ == '__main__':
    main()

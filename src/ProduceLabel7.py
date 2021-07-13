# 2021/7/7 using dict is quicker than dataframe
# try 3rd version of data
import copy
import heapq
import os

import numpy
import pandas as pd
import time


def init_proteins(in_filepath, generateFile=False, out_filePath='.\\Testcases2\\7_init_proteins.csv'):
    """
    initialize proteins and their index from .fasta file, note that duplicate protein may exist in the file
    will create protList and protDict
    :param
        filepath: string, path of one .fasta file that describes all of the proteins believed to be in the species
        generateFile: boolean, will output a file containing the items in protDict if true
        out_filePath: string, path of the above file, needs to be located in an existing directory
    :return
        protList: a list of all proteins from the input file
                    for storing the indices of the proteins
        protDict: a dictionary of all proteins from the input file
                    protein_formula: {'pepIndex': a list of indices of the protein's peptides that can be seen
                                                      from the .csv files,
                                      'count': int, how many peptides from one .csv file are in this protein
                                      'leftProtein': string, the protein formula without seenPeptides (not empty only
                                                     when we know there are unseenPeptides from the protein),
                                      'unseenPeptides': a list of peptides formula broken up from leftProtein,
                                      'list_intensity': a list of intensity value of the seen peptides
                                      'prot_total'：int, sum of the largest 3 intensity value
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

    # if sorting is needed, just sort protList
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


def find_all_possible_simple_peptides(protList, generateFile=False,
                                      out_filePath='.\\Testcases2\\7_simple_peptides.csv'):
    """
    break up all the proteins from .fasta file to get all possible simple peptides
    will create allPossiblePepDict
    :param
        protList: a list containing all the protein formula, should be the first return value of init_proteins()
        generateFile: boolean, will output a file containing the items in allPossiblePepDict if true
        out_filePath: string, path of the obove file, need to be located in an existing directory
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

    print("There are " + str(len(allPossiblePepDict)) + " simple peptides from all the proteins")

    if generateFile:
        sPepDf = pd.DataFrame({'sPeptide': allPossiblePepDict.keys(),
                               'protIndex': allPossiblePepDict.values()})
        sPepDf.to_csv(out_filePath, index=True, sep=',')

    return allPossiblePepDict


def init_peptides_from_one_file(filepath, colname=['peptide', 'intensity'],
                                generateFile=False, out_filePath='.\\Testcases2\\7_init_peptides_from_one_file.csv'):
    """
    initialize peptides and their index from the one of the .csv files
    will create pepList and pepDict
    will update pepDict['intensity']
    :param
        filepath: string, path of a .csv file storing information of peptides
        colname: list of string, the needed colume names,
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
        if numpy.isnan(df['intensity'][i]):
            continue
        pep = df['peptide'][i]
        intensity = df['intensity'][i]
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

    # if sorting is needed, just sort pepList
    return pepList, pepDict


def update_occurrence_of_one_file(allPossiblePepDict, protList, protDict, pepList, pepDict,
                                  generateFile=False,
                                  out_filePath_pep='.\\Testcases2\\7_pep_update_occurrence_of_one_file.csv',
                                  out_filePath_prot='.\\Testcases2\\7_prot_update_occurrence_of_one_file.csv'):
    """
    record the occurrence of peptides in proteins, find prot_total for each peptides
    for each peptides in .csv file (pepList), first search it among simple peptides, or search it in each protein
    will update pepDict[pep]['occurrence'],  pepDict[pep]['protIndex'],
                protDict_temp[prot]['pepIndex'], protDict_temp[prot]['count'], protDict_temp[prot]['list_intensity'],
                    protDict[prot]['prot_total']
    :param
        return values of find_all_possible_simple_peptides(), init_proteins() and init_peptides_from_one_file()
    :return:
        updated protDict_temp: has same column with protDict,
        updated pepDict,
        ignored: int, the number of proteins that has only one peptide appears in this .csv file
        valued: int, the number of proteins that has more than one peptides appears in this .csv file
    """
    # the original protDict should not change when processing different files
    protDict_temp = copy.deepcopy(protDict)

    # update occurrence
    pepIndex = 0
    for pep in pepList:  # pepList contains all the peptides from .tsv files

        # first search it in allPossiblePepDict
        if pep in allPossiblePepDict:  # allPossiblePepDict contains all the simple peptides
            protIndex = allPossiblePepDict[pep]
            pepDict[pep]['protIndex'] = protIndex
            pepDict[pep]['occurrence'] = len(protIndex)
            for i in protIndex:
                prot = protList[i]
                protDict_temp[prot]['pepIndex'].append(pepIndex)
                protDict_temp[prot]['list_intensity'].append(pepDict[pep]['intensity'])

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
                    # may locate at the head of a protein missing an M
                    elif prot.startswith('M' + pep):
                        occurrence += 1
                        pepDict[pep]['protIndex'].append(protIndex)
                        protDict_temp[prot]['pepIndex'].append(pepIndex)
                        protDict_temp[prot]['list_intensity'].append(pepDict[pep]['intensity'])
                    # may locate in the protein following a K
                    elif prot.find('K' + pep) != -1:
                        occurrence += 1
                        pepDict[pep]['protIndex'].append(protIndex)
                        protDict_temp[prot]['pepIndex'].append(pepIndex)
                        protDict_temp[prot]['list_intensity'].append(pepDict[pep]['intensity'])
                    # may locate in the protein following an R
                    elif prot.find('R' + pep) != -1:
                        occurrence += 1
                        pepDict[pep]['protIndex'].append(protIndex)
                        protDict_temp[prot]['pepIndex'].append(pepIndex)
                        protDict_temp[prot]['list_intensity'].append(pepDict[pep]['intensity'])
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

    return protDict_temp, pepDict, ignored, valued


def find_unseen_peptides_from_one_protein(protein, protDict_temp, pepList):
    """
    find unseen peptides of a present protein:
        use pepList and the indices in protDict_temp[prot]['pepIndex'] to find the formula of the seen peptides
        will update protDict_temp[prot]['leftProtein'] by removing seen peptides from original protein
        will update protDict_temp[prot]['unseenPeptides'] by breaking up the leftProtein via split_protein()
    :param
        protein: a String formula of one original present protein
        protDict_temp: first return value of update_occurrence_of_one_file()
        pepList: first return value of init_peptides_from_one_file()
    :return
        unseenPeptides: same as the return of split_protein(), a list of peptides broke up from leftProtein
    """
    leftProtein = protein

    seenPepFormula = []
    for pepIndex in protDict_temp[protein]['pepIndex']:
        seenPepFormula.append(pepList[pepIndex])
    seenPepFormula.sort(key=lambda x: len(x),
                        reverse=True)  # delete the longer peptides first  # when AAK and AAAK are in the same protein

    for pep in seenPepFormula:
        leftProtein = leftProtein.replace(pep, '')  # delete the seen peptides
    protDict_temp[protein]['leftProtein'] = leftProtein
    if len(leftProtein) == 0:
        return []
    unseenPeptides = split_protein(leftProtein)
    protDict_temp[protein]['unseenPeptides'] = unseenPeptides
    return unseenPeptides


def update_unseen_peptides_of_one_file(protList, protDict_temp, pepList, pepDict,
                                       generateFile=False,
                                       out_filePath_prot='.\\Testcases2\\7_prot_update_unseen_peptides_of_one_file.csv',
                                       out_filePath_unsPep='.\\Testcases2\\7_unsPep_update_unseen_peptides_of_one_file.csv'
                                       ):
    """
    will create unseenPepDict for unseen peptides
    will update unseenPepDict [pep]['protIndex'] , which are the proteins that produce the unseen peptides
    will update protDict_temp[prot]['leftProtein'] and protDict_temp[prot]['unseenPeptides'] when calling find_unseen_peptides_from_one_protein()
    :param
        protList: first return value of init_proteins()
        protDict_temp: first return value of update_occurrence()
        pepList: first return value of init_peptides()
        pepDict: second return value of update_occurrence()
        generateFile: boolean, will output two files containing items in protDict_temp and unseenPepDict respectively if true
        out_filePath_prot: string, path for the above first file, needs to be located in an existing directory
        out_filePath_unsPep: string, path for the above first file, needs to be located in an existing directory
    :return:
        updated protDict_temp
        unseenPepDict: a dictionary for the unseen peptides in one .csv file
            peptide_formula: {'score': 0, 'protIndex': list of indices of the proteins that contain the peptide}
        unsPepDf: pandas dataframe, containing the unseen peptides of this .csv file
                    columns = ['peptide', 'score']
    """
    unseenPepDict = {}
    for pep in pepList:
        if pepDict[pep]['occurrence'] == 1:
            protIndex = pepDict[pep]['protIndex'][0]
            protein = protList[protIndex]  # the protein that can be proved to be present
            if len(protDict_temp[protein]['leftProtein']) != 0:  # unseen peptides in this protein has been found
                continue
            unseenPeptides = find_unseen_peptides_from_one_protein(protein, protDict_temp, pepList)
            for pep in unseenPeptides:
                if pep not in unseenPepDict:
                    unseenPepDict[pep] = {'score': 0, 'protIndex': []}
                unseenPepDict[pep]['protIndex'].append(protIndex)

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
                                out_filePath='.\\Testcases2\\7_calculate_score_of_one_file.csv'):
    """
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
        pepDf: pandas dataframe
                column = ['peptide', 'intensity', 'occurrence', 'protIndex', 'score_0', 'score']
        pepDf_score: pandas dataframe
                column = ['peptide', 'score']
    """
    pepDf = pd.DataFrame({'peptide': pepDict.keys(),
                          'intensity': [pepDict[k]['intensity'] for k in pepDict],
                          'occurrence': [pepDict[k]['occurrence'] for k in pepDict],
                          'protIndex': [pepDict[k]['protIndex'] for k in pepDict]})
    pepDf = pepDf.explode('protIndex').reset_index()  # NaN if protIndex is []

    protDf = pd.DataFrame({'prot_total': [protDict_temp[k]['prot_total'] for k in protDict_temp]},
                          index=list(range(0, len(protDict_temp))))

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

    pepDf_score = pd.DataFrame(pepDf, columns=['peptide', 'score'])
    pepDf_score = pepDf_score.dropna().reset_index(drop=True)

    return pepDf, pepDf_score


# def gather_all_files(folderPath, protList, protDict, allPossiblePepDict):
#   in main method


def main():
    start = time.perf_counter()

    """
    Two input values:
        fastaFile: string, path of one .fasta file that contains protein formulas
        folderPath: string, path of one folder that contains the many files of peptides' infomation
    """

    fastaFile = '../uniprot-proteome_UP000000803.fasta'
    folderPath = '../Testcases3/flyquant'

    """
    Two debugging options:
        proteinInfo: boolean, will generate two files in protein folder if true
                        one for protein index and protein formula
                        one for simple peptides broke up from the proteins
        peptideInfo: boolean, will generate 3 files in peptide/filename folder for each file if true
                        one for protein ['',]
                        one for scored peptides []
                        one for unseen peptides []
    Remember to close the generated files before starting the programme again.
    """
    proteinInfo = True
    peptideInfo = True

    # initialize proteins and simple peptides
    generateFile = False
    protList, protDict = init_proteins(fastaFile, generateFile)
    allPossiblePepDict = find_all_possible_simple_peptides(protList)

    if proteinInfo:
        proteinInfoFolder = folderPath + "/../proteinInfo"
        if not os.path.exists(proteinInfoFolder):
            os.makedirs(proteinInfoFolder)
        protDf = pd.DataFrame({'protein': protDict.keys()})
        sPepDf = pd.DataFrame({'sPeptide': allPossiblePepDict.keys(),
                               'protIndex': allPossiblePepDict.values()})
        sPepDf.to_csv(proteinInfoFolder + '/simple peptides.csv', index=True, sep=',')
        protDf.to_csv(proteinInfoFolder + '/protein index.csv', index=True, sep=',')

    protTime = time.perf_counter()
    print('Finish initializing proteins, time: ' + str(protTime))

    # 3 dataframe to store the info from each .csv file
    fileInfoDf = pd.DataFrame(columns=['fileName', '#ignoredProtein', '#valuedProtein'])
    all_scoreDf = pd.DataFrame(columns=['peptide', 'score'])
    all_unsPepDf = pd.DataFrame(columns=['peptide', 'score'])

    # process the .csv files one by one
    filelist = os.listdir(folderPath)
    for file in filelist:
        if file.endswith('.csv') or file.endswith('.tsv'):
            print('processing file: ' + file)
            filePath = folderPath + '/' + file
            pepList, pepDict = init_peptides_from_one_file(filePath)
            protDict_temp, pepDict, ignored, valued = update_occurrence_of_one_file(allPossiblePepDict, protList,
                                                                                   protDict,
                                                                                   pepList, pepDict)
            protDict_temp, unseenPepDict, unsPepDf = update_unseen_peptides_of_one_file(protList, protDict_temp, pepList, pepDict)
            pepDf, pepDf_score = calculate_score_of_one_file(pepDict, protDict_temp)

            all_scoreDf = all_scoreDf.append(pepDf_score, ignore_index=True)
            all_unsPepDf = all_unsPepDf.append(unsPepDf,
                                               ignore_index=True)  # need drop duplicate, and delete the ones in all_scoreDf
            fileInfoDf = fileInfoDf.append([{'fileName': file, '#ignoredProtein': ignored, '#valuedProtein': valued}],
                                           ignore_index=True)

            if peptideInfo:
                peptideInfoFolder = folderPath + "/../peptideInfo/" + file
                if not os.path.exists(peptideInfoFolder):
                    os.makedirs(peptideInfoFolder)
                pepDf.to_csv(peptideInfoFolder + "/peptide score.csv", index=True, sep=',')
                protDf = pd.DataFrame({'protein': protDict_temp.keys(),
                                       'pepIndex': [protDict_temp[k]['pepIndex'] for k in protDict_temp],
                                       'count': [protDict_temp[k]['count'] for k in protDict_temp],
                                       'leftProtein': [protDict_temp[k]['leftProtein'] for k in protDict_temp],
                                       'unseenPeptides': [protDict_temp[k]['unseenPeptides'] for k in protDict_temp],
                                       'list_intensity': [protDict_temp[k]['list_intensity'] for k in protDict_temp],
                                       'prot_total': [protDict_temp[k]['prot_total'] for k in protDict_temp]})
                protDf.to_csv(peptideInfoFolder + "/protein information.csv", index=True, sep=',')
                unseenPepDf = pd.DataFrame({'peptide': unseenPepDict.keys(),
                                         'score': [unseenPepDict[k]['score'] for k in unseenPepDict],
                                         'protIndex': [unseenPepDict[k]['protIndex'] for k in
                                                       unseenPepDict]})
                unseenPepDf.to_csv(peptideInfoFolder + "/unseen peptide.csv", index=True, sep=',')

            pepTime = time.perf_counter()
            print('\t Finish time: ' + str(pepTime))

    # choose the most common score, will choose 100 if the quantities of 100 and 50 are the same
    all_scoreDf = all_scoreDf.groupby(['peptide'])['score'].agg(
        lambda x: x.value_counts(sort=False).index[0]).reset_index(
        drop=False)  # x.value_counts(sort=False)  # reset index to make it still a dataframe rather than a series

    # delete duplicate unseen peptides and delete the ones appearing in all_scoreDf
    all_unsPepDf = all_unsPepDf.drop_duplicates().reset_index(drop=True)
    all_unsPepDf_filter = all_unsPepDf[~ all_unsPepDf['peptide'].isin(all_scoreDf['peptide'])]

    # concatenating scored peptides and unseen peptides
    finalScoreDf = all_scoreDf.append(all_unsPepDf_filter, ignore_index=True)

    outputFolder = folderPath + '/../output'
    if not os.path.exists(outputFolder):
        os.makedirs(outputFolder)
    finalScoreDf.to_csv(outputFolder + '/score.csv', index=True, sep=',')
    fileInfoDf.to_csv(outputFolder + '/ignore.csv', index=True, sep=',')

    end = time.perf_counter()
    print('total processing time: \t' + str(end - start))


if __name__ == '__main__':
    main()

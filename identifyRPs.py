#!/usr/bin/env python3

import re, sys, csv, argparse
from os import mkdir, path, walk, listdir
from Bio import SeqIO

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('hmm_output', help='Folder storing tblout outputs from hmmsearch.')
    parser.add_argument('-o','--output_path', help='Output path for results. Default ($hmm_output/RPs.csv)')
    parser.add_argument('-fasta','--output_fasta', help='Provide location of original sequence files, this will write out folder with fasta files of extracted sequences.')
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def readFile(filePath):
    dataList = list()
    with open(filePath,'r') as fh:
        for curline in fh:
            if curline.startswith("#"):
                continue
            dataList.append(curline)
    return dataList

def constructFileArch(outPath):
    fileArchitecture = dict()
    rpList = list()
    for r, d, f in walk(outPath):
        if outPath == r:
            continue
        for file in f:
            if r.split(outPath+'/')[1] not in fileArchitecture.keys():
                fileArchitecture[r.split(outPath+'/')[1]] = dict()
            if file.endswith(".out"):
                rpList.append(file.replace('.out',''))
                fileArchitecture[r.split(outPath+'/')[1]][file.replace('.out','')] = readFile(path.join(r, file))
    return fileArchitecture, sorted(set(rpList))

def main (commandline_arguments):
    '''Main entry point'''
    comm_args = create_and_parse_argument_options(commandline_arguments)
    if comm_args.output_path is None:
        comm_args.output_path = f'{comm_args.hmm_output}/RPs.csv'
    fileArch, rpList = constructFileArch(comm_args.hmm_output)
    csvList, headerList = list(), ['Sequence file/rProtein topHit;numHits']
    headerList.extend(rpList)
    for seqFile, rProtOutputs in fileArch.items():
        seqFileOutput = list()
        seqFileOutput.append(seqFile)
        for rp in rpList:
            if len(rProtOutputs[rp]) > 0:
                topHit = rProtOutputs[rp][0].split()[0].split('|')[0]
            else:
                topHit = ''
            seqFileOutput.append(f'{topHit};{len(rProtOutputs[rp])}')
        csvList.append(seqFileOutput)
    with open(comm_args.output_path, 'w') as f:
        write = csv.writer(f)
        write.writerow(headerList)
        write.writerows(csvList)

    if comm_args.output_fasta:
        if not path.isdir(f'{comm_args.output_fasta}/RPhmmFASTA'):
            mkdir(f'{comm_args.output_fasta}/RPhmmFASTA')
        for rpNum, rpName in enumerate(headerList[1:], start=1):
            outSeqs = list()
            fileAndSeq = [(csvList[seqNum][0], sub[rpNum]) for seqNum,sub  in enumerate(csvList)]
            for origFasFile, seqName in fileAndSeq:
                for origFileName in listdir(comm_args.output_fasta):
                    if re.match(origFasFile, origFileName):
                        filteredSeqs = [seq for seq in SeqIO.parse(f'{comm_args.output_fasta}/{origFileName}', "fasta") if re.match(seqName.split(';')[0], seq.id.split("|")[0])]
                        break
                    else:
                        filteredSeqs = list()
                if len(filteredSeqs) == 1:
                    outSeqs.append(filteredSeqs[0])
            if len(outSeqs) > 0:
                with open(f'{comm_args.output_fasta}/RPhmmFASTA/{rpName}.faa', "w") as output_handle:
                    SeqIO.write(outSeqs, output_handle, "fasta")

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
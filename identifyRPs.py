#!/usr/bin/env python3

import sys, csv, json, argparse
from subprocess import Popen, PIPE
from os import remove, path, mkdir, walk

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('hmm_output', help='Folder storing tblout outputs from hmmsearch.')
    parser.add_argument('-o','--output_path', help='Output path for results. Default ($hmm_output/RPs.csv)')
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



if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
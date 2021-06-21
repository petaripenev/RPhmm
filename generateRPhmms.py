#!/usr/bin/env python3

import sys, json, argparse, urllib.request
from subprocess import Popen, PIPE
from os import remove, path, mkdir

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-o','--output_path', help='Output path to folder where hmms will be stored.')
    parser.add_argument('-t','--taxonomic_ids', help='List taxonomic ids for groups available on proteovision (Default: 2,2157,2759).', default='2,2157,2759')
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def get_avail_alns(taxIDs):
    '''Given a list of taxIDs, retrieves available alignments and produces a list of names to alnIDs.
    Skips aln names starting with ae or ab or eb.
    '''
    taxIDtoNames = dict()
    for taxID in taxIDs:
        try:
            with urllib.request.urlopen(f'https://proteovision.chemistry.gatech.edu/desire-api/taxonomic-groups/?format=json&taxgroup_id__in={taxID}') as url:
                data = json.loads(url.read().decode())
        except:
            sys.exit(f'Failed to fetch available alignments for taxID url {taxID}.')
        nameToids = list()
        for alnID in data['results'][0]['alignment_ids']:
            if alnID[1][:2] == 'ae' or alnID[1][:2] == 'ab' or alnID[1][:2] == 'be' or alnID[1][-2:] == 'ab':
                continue
            if alnID[2] != 'PROMALS3D':
                continue
            nameToids.append((alnID[1], alnID[0]))
        taxIDtoNames[taxID] = nameToids
    return taxIDtoNames

def get_alns (taxID, nameToids):
    nameToAlns = list()
    for name, alnID in nameToids:
        try:
            with urllib.request.urlopen(f"https://proteovision.chemistry.gatech.edu/ortholog-aln-api/{alnID}/{taxID}") as url:
                data = json.loads(url.read().decode())
        except:
            sys.exit(f'Failed to fetch alignment for taxID {taxID} and alnID {alnID}.')
        
        nameToAlns.append((name, data["Alignment"]))
    return nameToAlns

def makeHMMfromAln (taxID, alnString, hmmName, outpath):
    
    fh = open(f'{hmmName}_temp.faa', "w")
    fh.write(alnString)
    fh.close()

    if not path.isdir(f"{outpath}/{taxID}/"):
        mkdir(f"{outpath}/{taxID}/")

    pipe = Popen(f"trimal -in {hmmName}_temp.faa -out {hmmName}_temp_nogap.faa -nogaps", stdout=PIPE, shell=True)
    output = pipe.communicate()[0]

    pipe = Popen(f"hmmbuild {outpath}/{taxID}/{hmmName}.hmm {hmmName}_temp_nogap.faa", stdout=PIPE, shell=True)
    output = pipe.communicate()[0]
    remove(f'{hmmName}_temp.faa')
    remove(f'{hmmName}_temp_nogap.faa')

def main(commandline_arguments):
    '''Main entry point'''
    comm_args = create_and_parse_argument_options(commandline_arguments)
    taxIDs = comm_args.taxonomic_ids.split(',')
    taxIDtoNames = get_avail_alns(taxIDs)
    for taxID,nameToids  in taxIDtoNames.items():
        nameToAlns = get_alns (taxID, nameToids)
        for name, alnString in nameToAlns:
            makeHMMfromAln (taxID, alnString, name, comm_args.output_path)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
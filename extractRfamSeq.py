#!/usr/bin/env python3
'''
This script will take an Rfam output (tblout format 2), an Rfam name, and the path to the original contigs file.
It will then extract the matching Rfam hits to STDOUT.'''
import re, sys, argparse, gzip
from Bio import SeqIO
from copy import deepcopy

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('rfam_results', help='Path to Rfam results in tblout format 2.')
    parser.add_argument('rfam_name', help='Filter results by this string (e.g. LSU_rRNA_bacteria). Can do regex (e.g. LSU_rRNA_*).')
    parser.add_argument('input_contigs', help='Path to fasta file used to generate rfam_results. Can be gzipped.')
    parser.add_argument('-jp','--join_percentage', help='Join rfam results as long as their combined length times this value is larger than the gap between them (Default: 0.333).', default=0.333, type=float)
    parser.add_argument('-io','--include_overlaps', help='Include rfam results which are not marked with ^ (Default: False).', action='store_true')
    parser.add_argument('-ib','--include_between', help='When merging ranges add the region between them, instead of concatenating (Default: False).', action='store_true')
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def is_gz_file(filepath):
    'Check if file is gzipped.'
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'

def filterHitsByName(rfam_results, rfam_name, contigs, io):
    '''Collects the hits in a list of seq records'''
    hits = list()
    with open(rfam_results, "r") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            line = line.strip().split()
            if re.match(rfam_name,line[1]) and (line[19] == "^" or io):
                contig_name, strand = line[3], line[11]
                seqStart, seqEnd = int(line[9]), int(line[10])
                modelStart, modelEnd = int(line[7]), int(line[8])
                contigRecord = extractSequenceFromContigRecord(contigs, contig_name, strand, seqStart, seqEnd)
                contigRecord.id+=f':{seqStart}-{seqEnd}{strand}'
                contigRecord.description+=f' rfamName={line[1]}:{modelStart}-{modelEnd}'
                hits.append(contigRecord)
    return hits

def organizeHitsByContigName(hits):
    '''Organize hits in a dictionary with key their contig name and value a list of hits'''
    hitsByContig = dict()
    for hit in hits:
        name = hit.name
        if name not in hitsByContig.keys():
            hitsByContig[name] = list()
        hitsByContig[name].append(hit)
    return hitsByContig

def extractSequenceFromContigRecord(contigs, contig_name, strand, seqStart, seqEnd):
    filteredContigs = [x for x in contigs if x.name == contig_name]
    if len(filteredContigs) > 1:
        raise ValueError("More than one contig with the same name found in the contigs file.")
    if len(filteredContigs) == 0:
        raise ValueError("No contigs matched the the names from the rfam output! Did you input the correct rfam output file?")
    contigRecord = deepcopy(filteredContigs[0])
    if strand == "-":
        contigRecord.seq = contigRecord.seq[seqEnd:seqStart].reverse_complement()
    else:
        contigRecord.seq = contigRecord.seq[seqStart:seqEnd]
    return contigRecord

def organizeHitsByRange(contigHits):
    rangeToHit = dict()
    for hit in contigHits:
        if len(hit.id.split(';')) > 1:
            concatRanges = re.split(':|;', hit.id[:-1])[1:]
            sortedRange = sorted([int(max(x.split('-'))) for x in concatRanges])
        else:
            sortedRange = sorted([int(x) for x in hit.id.split(':')[1][:-1].split('-')])
        rangeToHit[tuple(sortedRange)] = hit
    return rangeToHit

def findInBetweenRanges(sortedRanges):
    inbetweenRanges = list()
    for i,x in enumerate(sortedRanges):
        if i+1 < len(sortedRanges):
            betweenRange = [x[1],sortedRanges[i+1][0]]
            if betweenRange[0] < betweenRange[1]:
                inbetweenRanges.append(betweenRange)
    inbetweenRanges.append([betweenRange[1]+1,betweenRange[1]+2])
    return inbetweenRanges

def generateModelRange(rangeToHit, sortedRanges, i, rfamString):
    if len(rfamString.split(';')) > 1:
        mRanges = [x for k,x in enumerate(re.split(':|;', rfamString)) if k%2 == 1]
        return [min([int(x) for x in '-'.join(mRanges).split('-')]), max([int(x) for x in '-'.join(mRanges).split('-')])]
    return [int(x) for x in rangeToHit[sortedRanges[i]].description.split(" rfamName=")[1].split(':')[1].split('-')]

def generateCombinedEntry(includeBetween, rangeToHit, sortedRanges, range1, range2, firstIndex, secondIndex):
    newOne = deepcopy(rangeToHit[sortedRanges[firstIndex]])
    if includeBetween:
        strand = newOne.id[-1]
        newOne.id = newOne.id.split(':')[0] + ':' + str(range1[0]) + '-' + str(range2[1]) + strand
        if strand == '-':
            newOne.seq = extractSequenceFromContigRecord(includeBetween, newOne.id.split(':')[0], strand, range2[1], range1[0]).seq
        else:
            newOne.seq = extractSequenceFromContigRecord(includeBetween, newOne.id.split(':')[0], strand, range1[0], range2[1]).seq
    else:
        newOne.id = newOne.id[:-1]+';'+rangeToHit[sortedRanges[secondIndex]].id.split(':')[1]
        newOne.seq = newOne.seq+rangeToHit[sortedRanges[secondIndex]].seq
    newOne.description = newOne.description+';'+rangeToHit[sortedRanges[secondIndex]].description.split(' rfamName=')[1]
    return newOne

def combineHits(contigHits, lengthPerc=0.333, includeBetween=False):
    '''Recursively combine hits which are at distance shorter than 0.333 length of the pair of hits.'''
    combinedHits = list()
    rangeToHit = organizeHitsByRange(contigHits)
    sortedRanges = sorted(rangeToHit.keys())
    
    inbetweenRanges = findInBetweenRanges(sortedRanges)

    for i,inbetweenRange in enumerate(inbetweenRanges):
        #This handles a case where there is only one hit in the range
        if i+1 == len(sortedRanges):
            combinedHits.append(rangeToHit[sortedRanges[i]])
            break
        range1, range2 = [int(x) for x in sortedRanges[i]], [int(x) for x in sortedRanges[i+1]]
        rfam1 = rangeToHit[sortedRanges[i]].description.split(" rfamName=")[1]
        rfam2 = rangeToHit[sortedRanges[i+1]].description.split(" rfamName=")[1]
        modelRange1 = generateModelRange(rangeToHit, sortedRanges, i, rfam1)
        modelRange2 = generateModelRange(rangeToHit, sortedRanges, i+1, rfam2)
        r1Set = set(range(modelRange1[0], modelRange1[1]))
        modelIntersect = r1Set.intersection(range(modelRange2[0], modelRange2[1]))
        
        #Order the hits by their model ranges
        if modelRange1 < modelRange2:
            firstIndex, secondIndex = i, i+1
        else:
            firstIndex, secondIndex = i+1, i
        
        #Check if combination should happen
        if len(modelIntersect) == 0 and rangeToHit[sortedRanges[i]].id[-1] == rangeToHit[sortedRanges[i+1]].id[-1] and \
        (inbetweenRange[1]-inbetweenRange[0] < lengthPerc*(range1[1]-range1[0]+range2[1]-range2[0])):     
            newCombinedEntry = generateCombinedEntry(includeBetween, rangeToHit, sortedRanges, range1, range2, firstIndex, secondIndex)
            combinedHits.append(newCombinedEntry)
            
            if i+2 < len(sortedRanges): #Make sure all other uncombined hits are added before checking for recursion.
                for leftOver in range(i+2, len(sortedRanges)):
                    combinedHits.append(rangeToHit[sortedRanges[leftOver]])
            if len(combinedHits) > 1:   #If there are more than one hits after combination, recurse. 
                return combineHits(combinedHits, lengthPerc=lengthPerc, includeBetween=includeBetween)
            break                       #If there is only one hit after combination, nothing to do => break
        else:                           #If combination should not happen, add the hit to the list
            combinedHits.append(rangeToHit[sortedRanges[i]])
    return combinedHits

def main(commandline_arguments):
    comm_args = create_and_parse_argument_options(commandline_arguments)
    rfam_results = comm_args.rfam_results
    rfam_name = comm_args.rfam_name
    input_contigs = comm_args.input_contigs

    if is_gz_file(input_contigs):
        input_contigs = gzip.open(input_contigs, mode='rt')
    
    contigs = list(SeqIO.parse(input_contigs, 'fasta'))
    hits = filterHitsByName(rfam_results, rfam_name, contigs, comm_args.include_overlaps)
    hitsByContig = organizeHitsByContigName(hits)

    combinedHitsByContig = dict()
    for contig in hitsByContig.keys():
        if len(hitsByContig[contig]) > 1:
            if comm_args.include_between:
                combinedHits = combineHits(hitsByContig[contig], lengthPerc=comm_args.join_percentage, includeBetween=contigs)
            else:
                combinedHits = combineHits(hitsByContig[contig], lengthPerc=comm_args.join_percentage)
            combinedHitsByContig[contig] = combinedHits

    for final_seq in combinedHitsByContig:
        SeqIO.write(combinedHitsByContig[final_seq], sys.stdout,'fasta')  

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
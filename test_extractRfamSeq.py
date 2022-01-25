import unittest
from extractRfamSeq import filterHitsByName, combineHits
from Bio import SeqIO, SeqRecord

TEST_CONTIGS = 'test_data/test_contigs.fasta'
TEST_RFAM_RESULTS = 'test_data/test_rfam_results.tblout'
TEST_RFAM_NAME = 'LSU_rRNA.*'

class TestExtraction(unittest.TestCase):

    def test_filterHitsByName(self):
        contigs = list(SeqIO.parse(TEST_CONTIGS, 'fasta'))
        hits = filterHitsByName(TEST_RFAM_RESULTS, TEST_RFAM_NAME, contigs, True)
        self.assertEqual(len(hits), 4)
        self.assertEqual(hits[0].id, '13_1_20cm_full_scaffold_1:3168-2079-')
        self.assertEqual(hits[1].id, '13_1_20cm_full_scaffold_1:53461-53302-')
        self.assertEqual(hits[2].id, '13_1_20cm_full_scaffold_1:53267-53074-')
        self.assertEqual(hits[3].id, '13_1_20cm_full_scaffold_1:3603-3521-')

    def test_combineHits(self):
        contigs = list(SeqIO.parse(TEST_CONTIGS, 'fasta'))
        hits = filterHitsByName(TEST_RFAM_RESULTS, TEST_RFAM_NAME, contigs, True)
        # Think about using fake data here
        # hits = [SeqRecord(Seq('ACGT'), id='test_contigs.fasta:1-100-', description='test_contigs.fasta:1-100-'),
        #         SeqRecord(Seq('ACGT'), id='test_contigs.fasta:1-100-', description='test_contigs.fasta:1-100-'),
        #         SeqRecord(Seq('ACGT'), id='test_contigs.fasta:1-100-', description='test_contigs.fasta:1-100-')]
        combinedHits = combineHits(hits)
        self.assertEqual(len(combinedHits), 2)
        self.assertEqual(combinedHits[0].id, '13_1_20cm_full_scaffold_1:3603-3521;3168-2079-')
        self.assertEqual(combinedHits[1].id, '13_1_20cm_full_scaffold_1:53461-53302;53267-53074-')
        self.assertEqual(len(combinedHits[0].seq), 1171)
        
        combinedHits = combineHits(hits, lengthPerc=0.1)
        self.assertEqual(len(combinedHits), 3)

        combinedHits = combineHits(hits, includeBetween=contigs)
        self.assertEqual(len(combinedHits), 2)
        self.assertEqual(combinedHits[0].id, '13_1_20cm_full_scaffold_1:2079-3603-')
        self.assertEqual(combinedHits[1].id, '13_1_20cm_full_scaffold_1:53074-53461-')
        self.assertEqual(len(combinedHits[0].seq), 1524)

    def test_sorting_combineHits(self):
        contigs = list(SeqIO.parse('test_data/test_contigs_sorting.fasta', 'fasta'))
        hits = filterHitsByName('test_data/test_sorting_rfam_results.tblout', TEST_RFAM_NAME, contigs, True)
        self.assertEqual(len(hits), 2)
        combinedHits = combineHits(hits)
        self.assertEqual(len(combinedHits), 2)
        combinedHits = combineHits(hits, includeBetween=contigs)
        self.assertEqual(len(combinedHits), 2)

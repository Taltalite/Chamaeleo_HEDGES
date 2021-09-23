import random
import unittest

from Chamaeleo.methods.flowed import DNAFountain


class TestEncodeDecode(unittest.TestCase):

    def setUp(self):
        random.seed(30)
        self.tool = DNAFountain(need_logs=False, redundancy=0.5)
        self.tool.segment_length = 200
        self.test_list = [[random.randint(0, 1) for _ in range(200)] for _ in range(10)]

    def test_list_to_sequence(self):
        dna_sequences = self.tool.silicon_to_carbon(self.test_list, 200 * 10).get("dna")
        self.assertEqual(dna_sequences, [
            ['A', 'G', 'G', 'A', 'A', 'G', 'G', 'G', 'A', 'A', 'A', 'G', 'A', 'A', 'A', 'A', 'T', 'A', 'T', 'A', 'A',
             'G', 'G', 'T', 'C', 'A', 'T', 'G', 'C', 'A', 'T', 'G', 'C', 'T', 'T', 'T', 'T', 'G', 'T', 'C', 'A', 'T',
             'A', 'T', 'T', 'G', 'T', 'G', 'A', 'A', 'G', 'A', 'G', 'C', 'A', 'A', 'T', 'A', 'C', 'A', 'A', 'A', 'C',
             'A', 'C', 'G', 'T', 'G', 'G', 'T', 'C', 'G', 'C', 'T', 'C', 'G', 'G', 'C', 'G', 'T', 'C', 'C', 'C', 'A',
             'A', 'T', 'A', 'A', 'A', 'A', 'G', 'A', 'G', 'A', 'A', 'C', 'C', 'T', 'T', 'T', 'G', 'G', 'A', 'T', 'C',
             'G', 'C', 'G', 'T', 'C', 'A', 'C', 'A', 'C', 'G', 'A'],
            ['A', 'C', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A', 'A', 'G', 'C', 'C', 'G', 'C',
             'A', 'A', 'C', 'G', 'G', 'C', 'G', 'G', 'T', 'C', 'T', 'A', 'C', 'A', 'T', 'A', 'T', 'T', 'G', 'A', 'T',
             'T', 'T', 'A', 'C', 'G', 'C', 'T', 'T', 'T', 'T', 'G', 'C', 'G', 'A', 'G', 'T', 'C', 'A', 'T', 'T', 'A',
             'G', 'G', 'C', 'G', 'C', 'C', 'C', 'G', 'T', 'C', 'G', 'T', 'G', 'A', 'T', 'G', 'A', 'G', 'C', 'A', 'C',
             'A', 'T', 'C', 'A', 'A', 'C', 'A', 'G', 'A', 'T', 'C', 'C', 'G', 'T', 'C', 'G', 'T', 'C', 'A', 'C', 'A',
             'A', 'T', 'A', 'G', 'T', 'T', 'G', 'A', 'A', 'A', 'C'],
            ['A', 'A', 'G', 'G', 'A', 'A', 'G', 'G', 'G', 'A', 'A', 'A', 'G', 'A', 'A', 'A', 'C', 'C', 'G', 'C', 'G',
             'C', 'G', 'T', 'C', 'T', 'A', 'G', 'C', 'C', 'T', 'G', 'G', 'T', 'C', 'C', 'G', 'G', 'A', 'A', 'G', 'G',
             'G', 'G', 'A', 'A', 'A', 'A', 'T', 'G', 'G', 'T', 'G', 'C', 'A', 'A', 'C', 'G', 'A', 'T', 'A', 'C', 'C',
             'A', 'G', 'G', 'T', 'G', 'C', 'C', 'G', 'A', 'C', 'A', 'A', 'G', 'A', 'C', 'G', 'G', 'T', 'G', 'G', 'C',
             'C', 'A', 'C', 'G', 'C', 'A', 'T', 'T', 'G', 'C', 'G', 'T', 'A', 'G', 'C', 'A', 'G', 'G', 'G', 'G', 'T',
             'A', 'G', 'C', 'A', 'T', 'G', 'A', 'G', 'T', 'C', 'G'],
            ['A', 'A', 'C', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'T', 'C', 'C', 'A',
             'A', 'A', 'T', 'G', 'C', 'C', 'T', 'G', 'C', 'G', 'C', 'A', 'G', 'T', 'G', 'G', 'G', 'G', 'T', 'G', 'A',
             'C', 'C', 'A', 'G', 'T', 'G', 'A', 'G', 'C', 'G', 'A', 'T', 'G', 'T', 'A', 'G', 'C', 'C', 'G', 'C', 'C',
             'G', 'G', 'A', 'A', 'G', 'A', 'T', 'A', 'C', 'C', 'G', 'A', 'G', 'A', 'C', 'G', 'G', 'G', 'T', 'T', 'A',
             'G', 'A', 'C', 'G', 'A', 'T', 'G', 'T', 'G', 'A', 'G', 'T', 'C', 'A', 'G', 'C', 'G', 'G', 'A', 'C', 'T',
             'A', 'C', 'A', 'T', 'G', 'G', 'A', 'C', 'A', 'G', 'C'],
            ['A', 'A', 'A', 'C', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'T', 'G', 'C', 'C', 'T',
             'C', 'T', 'G', 'C', 'T', 'A', 'C', 'T', 'G', 'C', 'G', 'A', 'C', 'T', 'G', 'C', 'A', 'C', 'C', 'A', 'A',
             'C', 'G', 'A', 'A', 'C', 'A', 'C', 'C', 'T', 'T', 'C', 'G', 'T', 'A', 'A', 'A', 'G', 'T', 'A', 'A', 'G',
             'T', 'G', 'C', 'T', 'C', 'G', 'A', 'A', 'C', 'T', 'A', 'T', 'T', 'A', 'G', 'C', 'G', 'G', 'A', 'A', 'A',
             'C', 'T', 'C', 'T', 'T', 'A', 'G', 'T', 'T', 'C', 'G', 'C', 'A', 'A', 'G', 'G', 'A', 'A', 'A', 'G', 'A',
             'G', 'T', 'T', 'A', 'T', 'T', 'G', 'A', 'G', 'T', 'G'],
            ['A', 'A', 'A', 'A', 'C', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'A', 'A', 'C', 'A', 'C', 'T', 'C', 'A', 'T',
             'G', 'A', 'C', 'T', 'C', 'C', 'C', 'G', 'C', 'C', 'G', 'A', 'C', 'C', 'T', 'T', 'T', 'A', 'A', 'C', 'G',
             'G', 'C', 'C', 'C', 'G', 'T', 'T', 'C', 'C', 'A', 'T', 'G', 'C', 'A', 'A', 'T', 'G', 'C', 'C', 'C', 'G',
             'G', 'T', 'C', 'C', 'T', 'G', 'A', 'G', 'C', 'T', 'C', 'T', 'A', 'A', 'G', 'G', 'C', 'A', 'C', 'A', 'G',
             'A', 'A', 'A', 'C', 'T', 'G', 'A', 'C', 'A', 'A', 'A', 'G', 'A', 'G', 'A', 'A', 'G', 'C', 'T', 'G', 'A',
             'G', 'T', 'G', 'T', 'A', 'G', 'G', 'C', 'G', 'A', 'C'],
            ['C', 'C', 'A', 'C', 'G', 'A', 'C', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'A', 'A', 'T', 'A', 'A', 'C', 'G',
             'C', 'C', 'T', 'G', 'C', 'A', 'A', 'T', 'A', 'T', 'A', 'T', 'A', 'A', 'T', 'T', 'A', 'T', 'G', 'A', 'G',
             'G', 'A', 'G', 'A', 'C', 'G', 'G', 'C', 'G', 'G', 'G', 'T', 'T', 'G', 'T', 'C', 'G', 'A', 'T', 'G', 'C',
             'A', 'G', 'T', 'C', 'T', 'A', 'T', 'C', 'C', 'A', 'A', 'T', 'A', 'A', 'C', 'T', 'C', 'G', 'T', 'T', 'C',
             'G', 'G', 'C', 'T', 'T', 'A', 'G', 'C', 'T', 'T', 'C', 'A', 'G', 'T', 'G', 'G', 'A', 'T', 'A', 'G', 'A',
             'C', 'A', 'C', 'T', 'G', 'C', 'G', 'A', 'T', 'C', 'C'],
            ['A', 'A', 'C', 'C', 'A', 'C', 'G', 'A', 'C', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'G', 'C', 'G', 'A', 'G',
             'A', 'A', 'G', 'C', 'A', 'C', 'T', 'G', 'T', 'T', 'C', 'T', 'A', 'T', 'C', 'T', 'T', 'C', 'T', 'C', 'G',
             'G', 'G', 'C', 'A', 'C', 'T', 'G', 'G', 'G', 'A', 'G', 'T', 'A', 'A', 'T', 'A', 'C', 'G', 'T', 'G', 'T',
             'T', 'G', 'T', 'A', 'T', 'G', 'C', 'A', 'G', 'T', 'T', 'A', 'G', 'T', 'G', 'A', 'G', 'A', 'G', 'T', 'G',
             'G', 'A', 'G', 'C', 'A', 'G', 'A', 'C', 'C', 'C', 'T', 'G', 'A', 'C', 'G', 'C', 'T', 'G', 'T', 'G', 'A',
             'G', 'C', 'A', 'G', 'G', 'T', 'G', 'G', 'T', 'G', 'T'],
            ['A', 'A', 'A', 'G', 'G', 'A', 'T', 'A', 'A', 'G', 'G', 'A', 'A', 'G', 'G', 'G', 'T', 'G', 'C', 'C', 'T',
             'C', 'T', 'G', 'C', 'T', 'A', 'C', 'T', 'G', 'C', 'G', 'A', 'C', 'T', 'G', 'C', 'A', 'C', 'C', 'A', 'A',
             'C', 'G', 'A', 'A', 'C', 'A', 'C', 'C', 'T', 'T', 'C', 'G', 'T', 'A', 'A', 'A', 'G', 'T', 'A', 'A', 'G',
             'T', 'G', 'C', 'T', 'C', 'G', 'A', 'A', 'C', 'T', 'A', 'T', 'T', 'A', 'G', 'C', 'G', 'G', 'A', 'A', 'A',
             'C', 'T', 'C', 'T', 'T', 'A', 'G', 'T', 'T', 'C', 'G', 'C', 'A', 'A', 'G', 'G', 'A', 'A', 'A', 'G', 'A',
             'G', 'T', 'T', 'A', 'T', 'T', 'G', 'A', 'G', 'T', 'G'],
            ['G', 'G', 'A', 'T', 'G', 'G', 'A', 'T', 'A', 'A', 'G', 'G', 'A', 'A', 'G', 'G', 'T', 'A', 'C', 'T', 'T',
             'C', 'C', 'G', 'C', 'G', 'G', 'C', 'T', 'A', 'G', 'G', 'T', 'G', 'C', 'G', 'C', 'A', 'C', 'C', 'G', 'G',
             'T', 'T', 'T', 'C', 'T', 'A', 'C', 'C', 'C', 'G', 'C', 'A', 'T', 'C', 'T', 'C', 'C', 'A', 'G', 'C', 'C',
             'A', 'G', 'G', 'G', 'G', 'T', 'A', 'T', 'G', 'C', 'G', 'T', 'T', 'G', 'G', 'A', 'C', 'T', 'T', 'T', 'A',
             'T', 'T', 'T', 'A', 'G', 'T', 'C', 'G', 'C', 'A', 'C', 'T', 'G', 'T', 'T', 'G', 'G', 'C', 'T', 'T', 'A',
             'C', 'C', 'A', 'A', 'C', 'T', 'G', 'G', 'T', 'T', 'T'],
            ['G', 'A', 'G', 'T', 'T', 'G', 'G', 'A', 'T', 'A', 'A', 'G', 'G', 'A', 'A', 'G', 'A', 'A', 'G', 'G', 'G',
             'A', 'A', 'A', 'C', 'A', 'G', 'G', 'T', 'C', 'A', 'C', 'G', 'T', 'C', 'T', 'G', 'T', 'A', 'C', 'A', 'C',
             'A', 'G', 'A', 'G', 'C', 'T', 'T', 'G', 'C', 'G', 'G', 'G', 'C', 'T', 'G', 'G', 'G', 'A', 'A', 'G', 'A',
             'C', 'A', 'T', 'G', 'G', 'A', 'G', 'T', 'G', 'C', 'G', 'A', 'G', 'C', 'G', 'T', 'G', 'T', 'A', 'C', 'T',
             'G', 'T', 'A', 'A', 'G', 'C', 'G', 'T', 'T', 'A', 'A', 'G', 'A', 'C', 'G', 'C', 'G', 'C', 'G', 'A', 'G',
             'C', 'G', 'G', 'C', 'T', 'A', 'C', 'T', 'A', 'C', 'A'],
            ['C', 'A', 'C', 'C', 'T', 'T', 'C', 'A', 'C', 'G', 'A', 'C', 'C', 'A', 'A', 'C', 'T', 'G', 'C', 'C', 'T',
             'C', 'T', 'G', 'C', 'T', 'A', 'C', 'T', 'G', 'C', 'G', 'A', 'C', 'T', 'G', 'C', 'A', 'C', 'C', 'A', 'A',
             'C', 'G', 'A', 'A', 'C', 'A', 'C', 'C', 'T', 'T', 'C', 'G', 'T', 'A', 'A', 'A', 'G', 'T', 'A', 'A', 'G',
             'T', 'G', 'C', 'T', 'C', 'G', 'A', 'A', 'C', 'T', 'A', 'T', 'T', 'A', 'G', 'C', 'G', 'G', 'A', 'A', 'A',
             'C', 'T', 'C', 'T', 'T', 'A', 'G', 'T', 'T', 'C', 'G', 'C', 'A', 'A', 'G', 'G', 'A', 'A', 'A', 'G', 'A',
             'G', 'T', 'T', 'A', 'T', 'T', 'G', 'A', 'G', 'T', 'G'],
            ['G', 'A', 'A', 'C', 'T', 'T', 'G', 'G', 'A', 'T', 'A', 'A', 'G', 'G', 'A', 'A', 'T', 'G', 'C', 'A', 'A',
             'G', 'T', 'C', 'G', 'T', 'G', 'C', 'G', 'C', 'T', 'A', 'G', 'A', 'G', 'T', 'T', 'A', 'C', 'C', 'C', 'G',
             'A', 'A', 'T', 'C', 'A', 'A', 'A', 'T', 'T', 'G', 'G', 'A', 'G', 'A', 'T', 'A', 'T', 'G', 'G', 'G', 'C',
             'C', 'T', 'A', 'G', 'C', 'A', 'G', 'G', 'A', 'A', 'A', 'G', 'A', 'G', 'C', 'C', 'T', 'G', 'C', 'C', 'C',
             'C', 'G', 'A', 'T', 'C', 'A', 'G', 'C', 'T', 'G', 'T', 'A', 'C', 'C', 'C', 'C', 'T', 'C', 'A', 'A', 'A',
             'T', 'G', 'T', 'T', 'T', 'A', 'A', 'C', 'T', 'C', 'C'],
            ['C', 'A', 'A', 'A', 'T', 'T', 'T', 'C', 'A', 'C', 'G', 'A', 'C', 'C', 'A', 'A', 'T', 'G', 'A', 'C', 'T',
             'C', 'T', 'T', 'C', 'A', 'T', 'A', 'G', 'G', 'G', 'C', 'C', 'G', 'C', 'G', 'A', 'G', 'C', 'A', 'A', 'C',
             'G', 'A', 'C', 'T', 'C', 'A', 'T', 'G', 'T', 'G', 'C', 'T', 'G', 'G', 'T', 'G', 'A', 'T', 'C', 'C', 'C',
             'C', 'G', 'T', 'T', 'C', 'G', 'C', 'G', 'A', 'C', 'C', 'C', 'A', 'C', 'A', 'T', 'C', 'T', 'C', 'C', 'C',
             'A', 'C', 'T', 'T', 'A', 'A', 'A', 'C', 'A', 'G', 'G', 'C', 'T', 'T', 'T', 'A', 'T', 'G', 'A', 'C', 'C',
             'T', 'T', 'C', 'C', 'C', 'A', 'T', 'T', 'T', 'C', 'T'],
            ['A', 'G', 'A', 'A', 'C', 'T', 'T', 'G', 'G', 'A', 'T', 'A', 'A', 'G', 'G', 'A', 'A', 'G', 'T', 'G', 'G',
             'C', 'G', 'G', 'G', 'T', 'C', 'T', 'A', 'C', 'G', 'G', 'T', 'C', 'C', 'T', 'T', 'T', 'G', 'A', 'T', 'G',
             'T', 'T', 'A', 'A', 'T', 'C', 'G', 'T', 'A', 'A', 'T', 'C', 'C', 'T', 'G', 'C', 'A', 'A', 'T', 'A', 'T',
             'G', 'A', 'G', 'A', 'C', 'C', 'T', 'C', 'C', 'G', 'T', 'T', 'C', 'A', 'A', 'G', 'A', 'C', 'G', 'T', 'A',
             'G', 'C', 'C', 'G', 'G', 'G', 'A', 'T', 'A', 'A', 'G', 'G', 'T', 'C', 'T', 'T', 'A', 'A', 'A', 'C', 'G',
             'T', 'G', 'G', 'C', 'A', 'G', 'C', 'A', 'A', 'G', 'T']
        ])

    def test_sequence_to_list(self):
        self.tool.decode_packets = 10
        bit_segments = self.tool.carbon_to_silicon([
            ['A', 'G', 'G', 'A', 'A', 'G', 'G', 'G', 'A', 'A', 'A', 'G', 'A', 'A', 'A', 'A', 'T', 'A', 'T', 'A', 'A',
             'G', 'G', 'T', 'C', 'A', 'T', 'G', 'C', 'A', 'T', 'G', 'C', 'T', 'T', 'T', 'T', 'G', 'T', 'C', 'A', 'T',
             'A', 'T', 'T', 'G', 'T', 'G', 'A', 'A', 'G', 'A', 'G', 'C', 'A', 'A', 'T', 'A', 'C', 'A', 'A', 'A', 'C',
             'A', 'C', 'G', 'T', 'G', 'G', 'T', 'C', 'G', 'C', 'T', 'C', 'G', 'G', 'C', 'G', 'T', 'C', 'C', 'C', 'A',
             'A', 'T', 'A', 'A', 'A', 'A', 'G', 'A', 'G', 'A', 'A', 'C', 'C', 'T', 'T', 'T', 'G', 'G', 'A', 'T', 'C',
             'G', 'C', 'G', 'T', 'C', 'A', 'C', 'A', 'C', 'G', 'A'],
            ['A', 'C', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A', 'A', 'G', 'C', 'C', 'G', 'C',
             'A', 'A', 'C', 'G', 'G', 'C', 'G', 'G', 'T', 'C', 'T', 'A', 'C', 'A', 'T', 'A', 'T', 'T', 'G', 'A', 'T',
             'T', 'T', 'A', 'C', 'G', 'C', 'T', 'T', 'T', 'T', 'G', 'C', 'G', 'A', 'G', 'T', 'C', 'A', 'T', 'T', 'A',
             'G', 'G', 'C', 'G', 'C', 'C', 'C', 'G', 'T', 'C', 'G', 'T', 'G', 'A', 'T', 'G', 'A', 'G', 'C', 'A', 'C',
             'A', 'T', 'C', 'A', 'A', 'C', 'A', 'G', 'A', 'T', 'C', 'C', 'G', 'T', 'C', 'G', 'T', 'C', 'A', 'C', 'A',
             'A', 'T', 'A', 'G', 'T', 'T', 'G', 'A', 'A', 'A', 'C'],
            ['A', 'A', 'G', 'G', 'A', 'A', 'G', 'G', 'G', 'A', 'A', 'A', 'G', 'A', 'A', 'A', 'C', 'C', 'G', 'C', 'G',
             'C', 'G', 'T', 'C', 'T', 'A', 'G', 'C', 'C', 'T', 'G', 'G', 'T', 'C', 'C', 'G', 'G', 'A', 'A', 'G', 'G',
             'G', 'G', 'A', 'A', 'A', 'A', 'T', 'G', 'G', 'T', 'G', 'C', 'A', 'A', 'C', 'G', 'A', 'T', 'A', 'C', 'C',
             'A', 'G', 'G', 'T', 'G', 'C', 'C', 'G', 'A', 'C', 'A', 'A', 'G', 'A', 'C', 'G', 'G', 'T', 'G', 'G', 'C',
             'C', 'A', 'C', 'G', 'C', 'A', 'T', 'T', 'G', 'C', 'G', 'T', 'A', 'G', 'C', 'A', 'G', 'G', 'G', 'G', 'T',
             'A', 'G', 'C', 'A', 'T', 'G', 'A', 'G', 'T', 'C', 'G'],
            ['A', 'A', 'C', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'T', 'C', 'C', 'A',
             'A', 'A', 'T', 'G', 'C', 'C', 'T', 'G', 'C', 'G', 'C', 'A', 'G', 'T', 'G', 'G', 'G', 'G', 'T', 'G', 'A',
             'C', 'C', 'A', 'G', 'T', 'G', 'A', 'G', 'C', 'G', 'A', 'T', 'G', 'T', 'A', 'G', 'C', 'C', 'G', 'C', 'C',
             'G', 'G', 'A', 'A', 'G', 'A', 'T', 'A', 'C', 'C', 'G', 'A', 'G', 'A', 'C', 'G', 'G', 'G', 'T', 'T', 'A',
             'G', 'A', 'C', 'G', 'A', 'T', 'G', 'T', 'G', 'A', 'G', 'T', 'C', 'A', 'G', 'C', 'G', 'G', 'A', 'C', 'T',
             'A', 'C', 'A', 'T', 'G', 'G', 'A', 'C', 'A', 'G', 'C'],
            ['A', 'A', 'A', 'C', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'T', 'G', 'C', 'C', 'T',
             'C', 'T', 'G', 'C', 'T', 'A', 'C', 'T', 'G', 'C', 'G', 'A', 'C', 'T', 'G', 'C', 'A', 'C', 'C', 'A', 'A',
             'C', 'G', 'A', 'A', 'C', 'A', 'C', 'C', 'T', 'T', 'C', 'G', 'T', 'A', 'A', 'A', 'G', 'T', 'A', 'A', 'G',
             'T', 'G', 'C', 'T', 'C', 'G', 'A', 'A', 'C', 'T', 'A', 'T', 'T', 'A', 'G', 'C', 'G', 'G', 'A', 'A', 'A',
             'C', 'T', 'C', 'T', 'T', 'A', 'G', 'T', 'T', 'C', 'G', 'C', 'A', 'A', 'G', 'G', 'A', 'A', 'A', 'G', 'A',
             'G', 'T', 'T', 'A', 'T', 'T', 'G', 'A', 'G', 'T', 'G'],
            ['A', 'A', 'A', 'A', 'C', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'A', 'A', 'C', 'A', 'C', 'T', 'C', 'A', 'T',
             'G', 'A', 'C', 'T', 'C', 'C', 'C', 'G', 'C', 'C', 'G', 'A', 'C', 'C', 'T', 'T', 'T', 'A', 'A', 'C', 'G',
             'G', 'C', 'C', 'C', 'G', 'T', 'T', 'C', 'C', 'A', 'T', 'G', 'C', 'A', 'A', 'T', 'G', 'C', 'C', 'C', 'G',
             'G', 'T', 'C', 'C', 'T', 'G', 'A', 'G', 'C', 'T', 'C', 'T', 'A', 'A', 'G', 'G', 'C', 'A', 'C', 'A', 'G',
             'A', 'A', 'A', 'C', 'T', 'G', 'A', 'C', 'A', 'A', 'A', 'G', 'A', 'G', 'A', 'A', 'G', 'C', 'T', 'G', 'A',
             'G', 'T', 'G', 'T', 'A', 'G', 'G', 'C', 'G', 'A', 'C'],
            ['C', 'C', 'A', 'C', 'G', 'A', 'C', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'A', 'A', 'T', 'A', 'A', 'C', 'G',
             'C', 'C', 'T', 'G', 'C', 'A', 'A', 'T', 'A', 'T', 'A', 'T', 'A', 'A', 'T', 'T', 'A', 'T', 'G', 'A', 'G',
             'G', 'A', 'G', 'A', 'C', 'G', 'G', 'C', 'G', 'G', 'G', 'T', 'T', 'G', 'T', 'C', 'G', 'A', 'T', 'G', 'C',
             'A', 'G', 'T', 'C', 'T', 'A', 'T', 'C', 'C', 'A', 'A', 'T', 'A', 'A', 'C', 'T', 'C', 'G', 'T', 'T', 'C',
             'G', 'G', 'C', 'T', 'T', 'A', 'G', 'C', 'T', 'T', 'C', 'A', 'G', 'T', 'G', 'G', 'A', 'T', 'A', 'G', 'A',
             'C', 'A', 'C', 'T', 'G', 'C', 'G', 'A', 'T', 'C', 'C'],
            ['A', 'A', 'C', 'C', 'A', 'C', 'G', 'A', 'C', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'G', 'C', 'G', 'A', 'G',
             'A', 'A', 'G', 'C', 'A', 'C', 'T', 'G', 'T', 'T', 'C', 'T', 'A', 'T', 'C', 'T', 'T', 'C', 'T', 'C', 'G',
             'G', 'G', 'C', 'A', 'C', 'T', 'G', 'G', 'G', 'A', 'G', 'T', 'A', 'A', 'T', 'A', 'C', 'G', 'T', 'G', 'T',
             'T', 'G', 'T', 'A', 'T', 'G', 'C', 'A', 'G', 'T', 'T', 'A', 'G', 'T', 'G', 'A', 'G', 'A', 'G', 'T', 'G',
             'G', 'A', 'G', 'C', 'A', 'G', 'A', 'C', 'C', 'C', 'T', 'G', 'A', 'C', 'G', 'C', 'T', 'G', 'T', 'G', 'A',
             'G', 'C', 'A', 'G', 'G', 'T', 'G', 'G', 'T', 'G', 'T'],
            ['A', 'A', 'A', 'G', 'G', 'A', 'T', 'A', 'A', 'G', 'G', 'A', 'A', 'G', 'G', 'G', 'T', 'G', 'C', 'C', 'T',
             'C', 'T', 'G', 'C', 'T', 'A', 'C', 'T', 'G', 'C', 'G', 'A', 'C', 'T', 'G', 'C', 'A', 'C', 'C', 'A', 'A',
             'C', 'G', 'A', 'A', 'C', 'A', 'C', 'C', 'T', 'T', 'C', 'G', 'T', 'A', 'A', 'A', 'G', 'T', 'A', 'A', 'G',
             'T', 'G', 'C', 'T', 'C', 'G', 'A', 'A', 'C', 'T', 'A', 'T', 'T', 'A', 'G', 'C', 'G', 'G', 'A', 'A', 'A',
             'C', 'T', 'C', 'T', 'T', 'A', 'G', 'T', 'T', 'C', 'G', 'C', 'A', 'A', 'G', 'G', 'A', 'A', 'A', 'G', 'A',
             'G', 'T', 'T', 'A', 'T', 'T', 'G', 'A', 'G', 'T', 'G'],
            ['G', 'G', 'A', 'T', 'G', 'G', 'A', 'T', 'A', 'A', 'G', 'G', 'A', 'A', 'G', 'G', 'T', 'A', 'C', 'T', 'T',
             'C', 'C', 'G', 'C', 'G', 'G', 'C', 'T', 'A', 'G', 'G', 'T', 'G', 'C', 'G', 'C', 'A', 'C', 'C', 'G', 'G',
             'T', 'T', 'T', 'C', 'T', 'A', 'C', 'C', 'C', 'G', 'C', 'A', 'T', 'C', 'T', 'C', 'C', 'A', 'G', 'C', 'C',
             'A', 'G', 'G', 'G', 'G', 'T', 'A', 'T', 'G', 'C', 'G', 'T', 'T', 'G', 'G', 'A', 'C', 'T', 'T', 'T', 'A',
             'T', 'T', 'T', 'A', 'G', 'T', 'C', 'G', 'C', 'A', 'C', 'T', 'G', 'T', 'T', 'G', 'G', 'C', 'T', 'T', 'A',
             'C', 'C', 'A', 'A', 'C', 'T', 'G', 'G', 'T', 'T', 'T'],
            ['G', 'A', 'G', 'T', 'T', 'G', 'G', 'A', 'T', 'A', 'A', 'G', 'G', 'A', 'A', 'G', 'A', 'A', 'G', 'G', 'G',
             'A', 'A', 'A', 'C', 'A', 'G', 'G', 'T', 'C', 'A', 'C', 'G', 'T', 'C', 'T', 'G', 'T', 'A', 'C', 'A', 'C',
             'A', 'G', 'A', 'G', 'C', 'T', 'T', 'G', 'C', 'G', 'G', 'G', 'C', 'T', 'G', 'G', 'G', 'A', 'A', 'G', 'A',
             'C', 'A', 'T', 'G', 'G', 'A', 'G', 'T', 'G', 'C', 'G', 'A', 'G', 'C', 'G', 'T', 'G', 'T', 'A', 'C', 'T',
             'G', 'T', 'A', 'A', 'G', 'C', 'G', 'T', 'T', 'A', 'A', 'G', 'A', 'C', 'G', 'C', 'G', 'C', 'G', 'A', 'G',
             'C', 'G', 'G', 'C', 'T', 'A', 'C', 'T', 'A', 'C', 'A'],
            ['C', 'A', 'C', 'C', 'T', 'T', 'C', 'A', 'C', 'G', 'A', 'C', 'C', 'A', 'A', 'C', 'T', 'G', 'C', 'C', 'T',
             'C', 'T', 'G', 'C', 'T', 'A', 'C', 'T', 'G', 'C', 'G', 'A', 'C', 'T', 'G', 'C', 'A', 'C', 'C', 'A', 'A',
             'C', 'G', 'A', 'A', 'C', 'A', 'C', 'C', 'T', 'T', 'C', 'G', 'T', 'A', 'A', 'A', 'G', 'T', 'A', 'A', 'G',
             'T', 'G', 'C', 'T', 'C', 'G', 'A', 'A', 'C', 'T', 'A', 'T', 'T', 'A', 'G', 'C', 'G', 'G', 'A', 'A', 'A',
             'C', 'T', 'C', 'T', 'T', 'A', 'G', 'T', 'T', 'C', 'G', 'C', 'A', 'A', 'G', 'G', 'A', 'A', 'A', 'G', 'A',
             'G', 'T', 'T', 'A', 'T', 'T', 'G', 'A', 'G', 'T', 'G'],
            ['G', 'A', 'A', 'C', 'T', 'T', 'G', 'G', 'A', 'T', 'A', 'A', 'G', 'G', 'A', 'A', 'T', 'G', 'C', 'A', 'A',
             'G', 'T', 'C', 'G', 'T', 'G', 'C', 'G', 'C', 'T', 'A', 'G', 'A', 'G', 'T', 'T', 'A', 'C', 'C', 'C', 'G',
             'A', 'A', 'T', 'C', 'A', 'A', 'A', 'T', 'T', 'G', 'G', 'A', 'G', 'A', 'T', 'A', 'T', 'G', 'G', 'G', 'C',
             'C', 'T', 'A', 'G', 'C', 'A', 'G', 'G', 'A', 'A', 'A', 'G', 'A', 'G', 'C', 'C', 'T', 'G', 'C', 'C', 'C',
             'C', 'G', 'A', 'T', 'C', 'A', 'G', 'C', 'T', 'G', 'T', 'A', 'C', 'C', 'C', 'C', 'T', 'C', 'A', 'A', 'A',
             'T', 'G', 'T', 'T', 'T', 'A', 'A', 'C', 'T', 'C', 'C'],
            ['C', 'A', 'A', 'A', 'T', 'T', 'T', 'C', 'A', 'C', 'G', 'A', 'C', 'C', 'A', 'A', 'T', 'G', 'A', 'C', 'T',
             'C', 'T', 'T', 'C', 'A', 'T', 'A', 'G', 'G', 'G', 'C', 'C', 'G', 'C', 'G', 'A', 'G', 'C', 'A', 'A', 'C',
             'G', 'A', 'C', 'T', 'C', 'A', 'T', 'G', 'T', 'G', 'C', 'T', 'G', 'G', 'T', 'G', 'A', 'T', 'C', 'C', 'C',
             'C', 'G', 'T', 'T', 'C', 'G', 'C', 'G', 'A', 'C', 'C', 'C', 'A', 'C', 'A', 'T', 'C', 'T', 'C', 'C', 'C',
             'A', 'C', 'T', 'T', 'A', 'A', 'A', 'C', 'A', 'G', 'G', 'C', 'T', 'T', 'T', 'A', 'T', 'G', 'A', 'C', 'C',
             'T', 'T', 'C', 'C', 'C', 'A', 'T', 'T', 'T', 'C', 'T'],
            ['A', 'G', 'A', 'A', 'C', 'T', 'T', 'G', 'G', 'A', 'T', 'A', 'A', 'G', 'G', 'A', 'A', 'G', 'T', 'G', 'G',
             'C', 'G', 'G', 'G', 'T', 'C', 'T', 'A', 'C', 'G', 'G', 'T', 'C', 'C', 'T', 'T', 'T', 'G', 'A', 'T', 'G',
             'T', 'T', 'A', 'A', 'T', 'C', 'G', 'T', 'A', 'A', 'T', 'C', 'C', 'T', 'G', 'C', 'A', 'A', 'T', 'A', 'T',
             'G', 'A', 'G', 'A', 'C', 'C', 'T', 'C', 'C', 'G', 'T', 'T', 'C', 'A', 'A', 'G', 'A', 'C', 'G', 'T', 'A',
             'G', 'C', 'C', 'G', 'G', 'G', 'A', 'T', 'A', 'A', 'G', 'G', 'T', 'C', 'T', 'T', 'A', 'A', 'A', 'C', 'G',
             'T', 'G', 'G', 'C', 'A', 'G', 'C', 'A', 'A', 'G', 'T']
        ]).get("bit")

        self.assertEqual(bit_segments, self.test_list)
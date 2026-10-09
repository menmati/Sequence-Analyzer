"""
Suite de pruebas unitarias para el analizador de secuencias FASTA.
"""

import unittest
from analizador.fasta import parse_fasta, reverse_complement, FastaRecord
from analizador.codons import calculate_gc_content, find_codons, find_orfs
from analizador.translation import translate_dna, translate_all_frames
from analizador.mutations import identify_point_mutations, classify_nucleotide_substitution, classify_protein_effect


class TestFastaParser(unittest.TestCase):

    def test_parse_fasta_single(self):
        fasta_content = ">seq1 Secuencia de prueba\nATGCATGC\nGGCC\n"
        records = parse_fasta(fasta_content)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].id, "seq1")
        self.assertEqual(records[0].description, "Secuencia de prueba")
        self.assertEqual(records[0].sequence, "ATGCATGCGGCC")

    def test_parse_fasta_multiple(self):
        fasta_content = ">seq1\nATGC\n>seq2 Segunda\nCGTA\n"
        records = parse_fasta(fasta_content)
        self.assertEqual(len(records), 2)
        self.assertEqual(records[0].id, "seq1")
        self.assertEqual(records[1].id, "seq2")
        self.assertEqual(records[1].sequence, "CGTA")

    def test_reverse_complement(self):
        seq = "ATGC"
        # Complemento: TACG -> Reverso: GCAT
        self.assertEqual(reverse_complement(seq), "GCAT")


class TestGCContent(unittest.TestCase):

    def test_gc_percentage(self):
        res1 = calculate_gc_content("GCGC")
        self.assertEqual(res1["gc_percent"], 100.0)

        res2 = calculate_gc_content("ATAT")
        self.assertEqual(res2["gc_percent"], 0.0)

        res3 = calculate_gc_content("ATGC")
        self.assertEqual(res3["gc_percent"], 50.0)
        self.assertEqual(res3["counts"]["A"], 1)
        self.assertEqual(res3["counts"]["G"], 1)
        self.assertEqual(res3["counts"]["C"], 1)
        self.assertEqual(res3["counts"]["T"], 1)


class TestCodonsAndORFs(unittest.TestCase):

    def test_find_start_and_stop_codons(self):
        # ATG en pos 1 (marco 1), TGA en pos 2 (marco 2), TAA en pos 10 (marco 1)
        # ATG AAA CCC TAA
        seq = "ATGAAACCCTAA"
        # Búsqueda sólo en marco 1
        starts_f1, stops_f1 = find_codons(seq, frame=1)
        self.assertEqual(len(starts_f1), 1)
        self.assertEqual(starts_f1[0].position, 1)
        self.assertEqual(starts_f1[0].frame, 1)
        self.assertEqual(len(stops_f1), 1)
        self.assertEqual(stops_f1[0].position, 10)
        self.assertEqual(stops_f1[0].frame, 1)

        # Búsqueda en todos los marcos (encuentra TGA en marco 2 además de TAA en marco 1)
        starts_all, stops_all = find_codons(seq, frame=0)
        self.assertEqual(len(starts_all), 1)
        self.assertEqual(len(stops_all), 2)

    def test_find_orfs(self):
        seq = "ATGAAACCCTAA"
        orfs = find_orfs(seq, min_length_nt=12)
        self.assertEqual(len(orfs), 1)
        self.assertEqual(orfs[0].start_pos, 1)
        self.assertEqual(orfs[0].stop_pos, 12)
        self.assertEqual(orfs[0].length_nt, 12)
        self.assertEqual(orfs[0].length_aa, 4)


class TestTranslation(unittest.TestCase):

    def test_standard_translation(self):
        # ATG (M), GAA (E), GTC (V), TAA (*)
        seq = "ATGGAAGTCTAA"
        prot = translate_dna(seq, frame=1)
        self.assertEqual(prot, "MEV*")

    def test_translation_frames(self):
        seq = "CATGGAAGTCTAA"
        # Marco 1: CAT (H), GGA (G), AGT (S), CTA (L) -> HGSL
        # Marco 2: ATG (M), GAA (E), GTC (V), TAA (*) -> MEV*
        all_frames = translate_all_frames(seq)
        self.assertTrue(all_frames["+2"].startswith("MEV*"))


class TestMutations(unittest.TestCase):

    def test_transition_transversion(self):
        self.assertEqual(classify_nucleotide_substitution("A", "G"), "Transición")
        self.assertEqual(classify_nucleotide_substitution("C", "T"), "Transición")
        self.assertEqual(classify_nucleotide_substitution("A", "C"), "Transversión")
        self.assertEqual(classify_nucleotide_substitution("G", "T"), "Transversión")

    def test_synonymous_mutation(self):
        # GAA (E) -> GAG (E)
        ref = "ATGGAA"
        alt = "ATGGAG"
        muts = identify_point_mutations(ref, alt, frame=1)
        self.assertEqual(len(muts), 1)
        self.assertEqual(muts[0].position, 6)
        self.assertEqual(muts[0].ref_base, "A")
        self.assertEqual(muts[0].alt_base, "G")
        self.assertEqual(muts[0].type, "Transición")
        self.assertIn("Sinónima", muts[0].effect)

    def test_missense_mutation(self):
        # GAG (E) -> GTG (V) (Mutación drepanocitosis)
        ref = "ATGGAG"
        alt = "ATGGTG"
        muts = identify_point_mutations(ref, alt, frame=1)
        self.assertEqual(len(muts), 1)
        self.assertEqual(muts[0].position, 5)
        self.assertEqual(muts[0].ref_base, "A")
        self.assertEqual(muts[0].alt_base, "T")
        self.assertEqual(muts[0].type, "Transversión")
        self.assertIn("Missense", muts[0].effect)
        self.assertEqual(muts[0].ref_aa, "E")
        self.assertEqual(muts[0].alt_aa, "V")

    def test_nonsense_mutation(self):
        # AAG (K) -> TAG (*)
        ref = "ATGAAG"
        alt = "ATGTAG"
        muts = identify_point_mutations(ref, alt, frame=1)
        self.assertEqual(len(muts), 1)
        self.assertIn("Nonsense", muts[0].effect)


if __name__ == "__main__":
    unittest.main()

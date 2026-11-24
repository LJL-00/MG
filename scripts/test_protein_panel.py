"""Synthetic tests for protein-to-gene mapping and explicit panel selection."""
from pathlib import Path
from importlib import util
import tempfile
import unittest

spec = util.spec_from_file_location('protein_panel', Path(__file__).with_name('09_protein_panel.py'))
panel = util.module_from_spec(spec); spec.loader.exec_module(panel)


def hit(query, subject, score, identity=50, length=30, evalue=1e-5):
    return dict(qseqid=query, sseqid=subject, pident=str(identity), length=str(length), evalue=str(evalue), bitscore=str(score))


class PanelTests(unittest.TestCase):
    def test_distinct_query_best_and_union(self):
        rows=[hit('query_a','protein_a',50),hit('query_a','protein_b',20),hit('query_b','protein_a',40)]
        retained,best,genes=panel.build_panel(rows,{'protein_a':'gene_a','protein_b':'gene_b'},{'manual_gene'},{'shared_gene'})
        self.assertEqual(len(retained),3)
        self.assertEqual([x['gene_id'] for x in best],['gene_a','gene_a'])
        self.assertEqual({x['gene_id'] for x in genes},{'gene_a','gene_b','manual_gene','shared_gene'})
        self.assertEqual(sum(x['focal'] for x in genes),2)
        self.assertEqual(sum(x['extended'] for x in genes),3)

    def test_identity_length_then_input_ties(self):
        rows=[hit('q','a',10,40,20),hit('q','b',10,50,20),hit('q','c',10,50,30),hit('q','d',10,50,30)]
        _,best,_=panel.build_panel(rows,{x:x for x in 'abcd'})
        self.assertEqual(best[0]['sseqid'],'c')

    def test_explicit_thresholds(self):
        rows=[hit('q','a',100,evalue=1),hit('q','b',10,evalue=1e-8)]
        self.assertEqual(panel.build_panel(rows,{'a':'a','b':'b'})[1][0]['sseqid'],'a')
        self.assertEqual(panel.build_panel(rows,{'a':'a','b':'b'},max_evalue=1e-5)[1][0]['sseqid'],'b')

    def test_mapping_is_required(self):
        with self.assertRaisesRegex(ValueError,'No gene mapping'):
            panel.build_panel([hit('q','missing',10)],{})

    def test_parent_graph_and_multiple_exons(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'reference.gff'
            p.write_text('##gff-version 3\n'
                'ctg\tx\tgene\t1\t90\t.\t+\t.\tID=gene-alpha;locus_tag=alpha\n'
                'ctg\tx\tmRNA\t1\t90\t.\t+\t.\tID=rna-alpha;Parent=gene-alpha\n'
                'ctg\tx\tCDS\t1\t30\t.\t+\t0\tID=cds-a;Parent=rna-alpha;protein_id=protein_a\n'
                'ctg\tx\tCDS\t60\t90\t.\t+\t0\tID=cds-a;Parent=rna-alpha;protein_id=protein_a\n')
            self.assertEqual(panel.subject_mapping(gff=p),{'protein_a':'alpha'})

    def test_mapping_conflict(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'map.tsv';p.write_text('sseqid\tgene_id\na\tone\na\ttwo\n')
            with self.assertRaisesRegex(ValueError,'multiple genes'):
                panel.subject_mapping(mapping=p)

    def test_gene_list_header_and_prefix(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'genes.tsv';p.write_text('gene_id\nfake_a\ngene-fake_b\n')
            self.assertEqual(panel.gene_list(p),{'fake_a','fake_b'})


if __name__=='__main__': unittest.main()

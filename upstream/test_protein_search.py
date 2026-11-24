"""Synthetic BLAST command planning and FASTA input checks."""
from pathlib import Path
from unittest import mock
import json
import sys
import tempfile
import unittest
import run_protein_search as search


class ProteinSearchTests(unittest.TestCase):
    def test_empty_and_duplicate_records(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'protein.faa'
            for content in ('>\nABC\n','>a\n>b\nABC\n','>a\nABC\n>a\nDEF\n',''):
                p.write_text(content)
                with self.assertRaises(ValueError): search.fasta_ids(p)

    def test_complete_dry_run(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);queries=folder/'queries.faa';proteins=folder/'proteins.faa'
            queries.write_text('>query_example\nMSAA\n')
            proteins.write_text('>protein_one\nMSAA\n>protein_two\nMCCC\n')
            output=folder/'out'
            args=['run_protein_search.py','--queries',str(queries),'--proteins',str(proteins),
                  '--evalue','0.0001','--seg','no','--max-hsps','1','--output',str(output),'--dry-run']
            with mock.patch.object(sys,'argv',args),mock.patch('subprocess.run',side_effect=AssertionError('dry-run executed a tool')):
                search.main()
            plan=json.loads((output/'search_plan.json').read_text())
            self.assertEqual(plan['max_target_seqs'],2)
            self.assertIn('qcovhsp',plan['blastp_command'][plan['blastp_command'].index('-outfmt')+1])
            self.assertEqual(plan['blastp_command'][-2:],['-max_hsps','1'])
            self.assertEqual(plan['evalue'],0.0001)
            self.assertEqual(plan['seg'],'no')

    def test_default_search_settings(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);queries=folder/'queries.faa';proteins=folder/'proteins.faa'
            queries.write_text('>query_example\nMSAA\n')
            proteins.write_text('>protein_example\nMSAA\n')
            output=folder/'out'
            args=['run_protein_search.py','--queries',str(queries),'--proteins',str(proteins),
                  '--output',str(output),'--dry-run']
            with mock.patch.object(sys,'argv',args),mock.patch('subprocess.run',side_effect=AssertionError('dry-run executed a tool')):
                search.main()
            plan=json.loads((output/'search_plan.json').read_text())
            self.assertEqual(plan['evalue'],1e-5)
            self.assertEqual(plan['seg'],'yes')

    def test_server_version_configuration(self):
        config=search.read_config(search.BASE/'config/upstream.json')
        for name in ('BLASTP','Makeblastdb'):
            self.assertEqual(config['tools'][name]['required_version'],'2.17.0')


if __name__=='__main__': unittest.main()

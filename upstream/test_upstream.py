#!/usr/bin/env python3
"""Small synthetic checks for new entry points, separate from study-data replays."""
from pathlib import Path
import copy,csv,json,subprocess,sys,tempfile,unittest
import numpy as np
from common import BASE,Runner,check_tool,load_manifest
from strict_core import build_core,inspect_tree_alignment,qualified,read_fasta,sample_consensus

def vcf_row(chrom,pos,ref,alt,qual=100,dp=10,ao=9,gt='1/1'):
    return [chrom,str(pos),'.',ref,alt,str(qual),'PASS','.','GT:DP:AO',f'{gt}:{dp}:{ao}']

def write_vcf(path,records):
    with path.open('w',encoding='utf-8') as handle:
        handle.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n')
        for row in records:handle.write('\t'.join(row)+'\n')

def write_fasta(path,records):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_bytes(b''.join(b'>'+name.encode()+b'\n'+seq+b'\n' for name,seq in records.items()))

class UpstreamTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.reference={'chrA':b'A'*80,'chrB':b'C'*30}
        self.ref=self.root/'reference.fa';write_fasta(self.ref,self.reference)
        self.ids=['s1','s2','s3'];self.manifest=self.root/'reads.tsv'
        with self.manifest.open('w',encoding='utf-8',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=['sample','r1','r2','analysis_set','assemble','enabled'],delimiter='\t');writer.writeheader()
            for sid in self.ids:
                writer.writerow({'sample':sid,'r1':sid+'.1.fq','r2':sid+'.2.fq','analysis_set':'Main','assemble':'1','enabled':'1'})
        self.substitutions={
            's1':[vcf_row('chrA',5,'A','G'),vcf_row('chrA',7,'A','G')],
            's2':[vcf_row('chrA',12,'AA','GT'),vcf_row('chrB',25,'C','T')],
            's3':[vcf_row('chrA',20,'A','C',qual=99),vcf_row('chrA',21,'A','C',dp=9,ao=9),vcf_row('chrA',22,'A','C',ao=8),vcf_row('chrA',23,'A','C',gt='0/1')]}
        for sid in self.ids:
            folder=self.root/'snippy'/sid;folder.mkdir(parents=True)
            write_fasta(folder/'reference/ref.fa',self.reference)
            aligned={k:bytearray(v) for k,v in self.reference.items()}
            if sid=='s1':aligned['chrA'][31]=ord('n')
            if sid=='s2':aligned['chrA'][29]=ord('N')
            if sid=='s3':aligned['chrA'][30]=ord('-');aligned['chrA'][6]=ord('N')
            write_fasta(folder/'snps.aligned.fa',{k:bytes(v) for k,v in aligned.items()})
            write_vcf(folder/'snps.subs.vcf',self.substitutions[sid])
            indels=[vcf_row('chrA',50,'AA','A'),vcf_row('chrB',1,'CC','C')] if sid=='s1' else []
            write_vcf(folder/'snps.vcf',self.substitutions[sid]+indels)
        self.out=self.root/'core'

    def tearDown(self):self.tmp.cleanup()

    def test_quality_boundaries(self):
        row=vcf_row('chrA',1,'A','G')
        self.assertTrue(qualified(row,{'GT':'1/1','DP':'10','AO':'9'}))
        for fmt in [{'GT':'1/1','DP':'9','AO':'9'},{'GT':'1/1','DP':'10','AO':'8'},{'GT':'0/1','DP':'10','AO':'9'},{'GT':'1/1','DP':'10'}]:self.assertFalse(qualified(row,fmt))
        row[5]='99';self.assertFalse(qualified(row,{'GT':'1/1','DP':'10','AO':'9'}))

    def test_core_MNS_missing_and_indel_edges(self):
        report=build_core(self.manifest,self.ref,self.root/'snippy',self.out,expected_samples=3)
        self.assertEqual(report['core_bp'],68);self.assertEqual(report['variable_sites'],4);self.assertEqual(report['constant_sites'],64)
        seq=read_fasta(self.out/'full_core.fasta.gz')
        self.assertTrue(all(set(x)<=set(b'ACGT') for x in seq.values()))
        # The missing s3 call at chrA:7 removes this column in every isolate.
        self.assertEqual(sum(a!=b for a,b in zip(seq['s1'],seq['s3'])),1)
        self.assertEqual(sum(a!=b for a,b in zip(seq['s1'],seq['s2'])),4)
        self.assertEqual(sum(a!=b for a,b in zip(seq['s2'],seq['s3'])),3)
        offsets={'chrA':0,'chrB':80};_,mask,_=sample_consensus(self.root/'snippy/s1',self.reference,list(self.reference),offsets)
        self.assertTrue(np.all(mask[39:61]));self.assertFalse(mask[38]);self.assertFalse(mask[61])
        self.assertTrue(np.all(mask[80:92]));self.assertFalse(mask[79]);self.assertFalse(mask[92])

    def test_tree_full_sites_and_SNP_only_guard(self):
        build_core(self.manifest,self.ref,self.root/'snippy',self.out,expected_samples=3)
        info=inspect_tree_alignment(self.out/'full_core.fasta.gz',self.ids)
        self.assertEqual(info['constant_sites'],64)
        with self.assertRaisesRegex(ValueError,'SNP-only'):inspect_tree_alignment(self.out/'SNPs_only.fasta.gz',self.ids)
        with self.assertRaisesRegex(ValueError,'IDs'):inspect_tree_alignment(self.out/'full_core.fasta.gz',['s1','s2'])
        invalid=self.root/'ambiguous.fa';write_fasta(invalid,{'s1':b'ACN','s2':b'ACC','s3':b'ACT'})
        with self.assertRaisesRegex(ValueError,'ACGT'):inspect_tree_alignment(invalid,self.ids)

    def test_version_match_and_mismatch(self):
        cfg={'tools':{'FastQC':{'command':[sys.executable,'-c',"print('FastQC v0.12.1')"],'version_args':[],'version_regex':r'FastQC\s+v?([\d.]+)','required_version':'0.12.1','version_source':'server inventory','manuscript_version':'0.11.9'}}}
        self.assertEqual(check_tool(cfg,'FastQC')['status'],'MATCH')
        cfg['tools']['FastQC']['required_version']='0.11.9'
        with self.assertRaisesRegex(RuntimeError,'configuration requires'):check_tool(cfg,'FastQC')

    def test_server_versions_reject_previous_executables(self):
        cfg=json.loads((BASE/'config/upstream.json').read_text(encoding='utf-8'))
        previous={'Trimmomatic':('0.39','0.39'),'SPAdes':('SPAdes v3.15.5','3.15.5'),'BUSCO':('BUSCO 5.2.2','5.2.2'),'Bowtie2':('bowtie2-align-s version 2.5.2','2.5.2'),'SAMtools':('samtools 1.17','1.17')}
        for name,(output,_) in previous.items():
            spec=cfg['tools'][name];spec['command']=[sys.executable,'-c',f'print({output!r})'];spec['version_args']=[]
            with self.subTest(software=name),self.assertRaisesRegex(RuntimeError,'configuration requires'):check_tool(cfg,name)

    def test_snippy_dependency_failure_stops_before_analysis(self):
        cfg={'tools':{'Snippy':{'command':[sys.executable,'-c',"import sys; print('snippy 4.6.0'); sys.exit(1 if '--check' in sys.argv else 0)"],'version_args':[],'version_regex':r'snippy\s+([\d.]+)','required_version':'4.6.0','version_source':'server inventory','preflight_args':['--check']}}}
        with self.assertRaisesRegex(RuntimeError,'dependency preflight failed'):check_tool(cfg,'Snippy')

    def test_changed_version_does_not_reuse_previous_outputs(self):
        version=self.root/'version.txt';version.write_text('1.0')
        tool=self.root/'changing_tool.py'
        tool.write_text("from pathlib import Path\nimport sys\nv=Path(sys.argv[1]).read_text()\nif '--version' in sys.argv:print('Mock '+v)\nelse:Path(sys.argv[-1]).write_text(v)\n",encoding='utf-8')
        cfg={'tools':{'Mock':{'command':[sys.executable,str(tool),str(version)],'version_args':['--version'],'version_regex':r'Mock\s+([\d.]+)','required_version':'1.0','version_source':'test'}}}
        out=self.root/'new_outputs';product=out/'product.txt'
        first=Runner(out,cfg);first.run('mock',first.tool('Mock')+[str(product)],[product])
        self.assertEqual(product.read_text(),'1.0')
        again=Runner(out,cfg);again.run('mock',again.tool('Mock')+[str(product)],[product])
        version.write_text('2.0');cfg['tools']['Mock']['required_version']='2.0'
        changed=Runner(out,cfg)
        with self.assertRaisesRegex(RuntimeError,'completion record differs'):
            changed.run('mock',changed.tool('Mock')+[str(product)],[product])
        self.assertEqual(product.read_text(),'1.0')

    def test_manifest_duplicate_guard(self):
        self.assertEqual(len(load_manifest(self.manifest)),3)
        with self.manifest.open('a',encoding='utf-8') as handle:handle.write('s1\t1.fq\t2.fq\tMain\t1\t1\n')
        with self.assertRaisesRegex(ValueError,'duplicate'):load_manifest(self.manifest)

    def test_cli_plan_and_explicit_tree_choice(self):
        cfg=json.loads((BASE/'config/upstream.json').read_text(encoding='utf-8'))
        cfg['expected_main_isolates']=3;cfg['reference_fasta']=str(self.ref);cfg['adapters_fasta']=str(self.ref)
        config=self.root/'config.json';config.write_text(json.dumps(cfg),encoding='utf-8')
        command=[sys.executable,str(BASE/'upstream/workflow.py'),'all','--config',str(config),'--manifest',str(self.manifest),'--output',str(self.root/'planned'),'--dry-run']
        result=subprocess.run(command,capture_output=True,text=True);self.assertEqual(result.returncode,2);self.assertIn('explicitly',result.stderr)
        result=subprocess.run(command+['--new-tree-from-core'],capture_output=True,text=True);self.assertEqual(result.returncode,0,result.stderr)
        plan=json.loads((self.root/'planned/command_plan.json').read_text(encoding='utf-8'))
        self.assertEqual(plan['mode'],'PLAN_ONLY_NOT_EXECUTED')
        commands=[r['command'] for r in plan['commands'] if 'command' in r]
        snippy=[x for x in commands if '--minfrac' in x]
        self.assertEqual(len(snippy),3);self.assertTrue(all(x[x.index('--minfrac')+1]=='0.90' for x in snippy))
        self.assertTrue(any('--isolate' in x for x in commands));self.assertTrue(any('GTR+F+R4' in x and '-B' in x and '-alrt' in x for x in commands))
        self.assertFalse(any('--ctgs' in x for x in commands))
        multiqc=next(r for r in plan['commands'] if r.get('name')=='multiqc')
        self.assertEqual(len(multiqc['inputs']),12)
        self.assertTrue(all(x.endswith('_fastqc.zip') for x in multiqc['inputs']))
        self.assertNotIn(str(self.root/'planned/fastqc'),multiqc['command'])

    def test_missing_input_is_actionable(self):
        with self.assertRaisesRegex(FileNotFoundError,'needs nonempty input files'):
            build_core(self.manifest,self.ref,self.root/'missing',self.out,expected_samples=3)

    def test_mock_tool_processes_do_not_precreate_Snippy_or_BUSCO_dirs(self):
        mock=self.root/'mock_tools.py'
        mock.write_text('''from pathlib import Path
import sys
tool=sys.argv[1];args=sys.argv[2:]
if '--version' in args:
 print({'Snippy':'snippy 4.6.0','QUAST':'QUAST v5.3.0','BUSCO':'BUSCO 6.0.0'}[tool]);raise SystemExit(0)
if '--check' in args:
 print('synthetic dependency interface only');raise SystemExit(0)
if tool=='Snippy':
 folder=Path(args[args.index('--outdir')+1]);assert not folder.exists(), 'Snippy output directory was precreated'
 for name in ['snps.aligned.fa','snps.vcf','snps.subs.vcf','reference/ref.fa']:
  out=folder/name;out.parent.mkdir(parents=True,exist_ok=True);out.write_text('synthetic interface test only\\n')
elif tool=='QUAST':
 folder=Path(args[args.index('-o')+1]);folder.mkdir(parents=True,exist_ok=True);(folder/'report.tsv').write_text('synthetic interface test only\\n')
elif tool=='BUSCO':
 name=args[args.index('-o')+1];folder=Path(args[args.index('--out_path')+1])/name
 assert not folder.exists(), 'BUSCO output directory was precreated'
 folder.mkdir(parents=True);(folder/f'short_summary.specific.fungi_odb10.{name}.txt').write_text('synthetic interface test only\\n')
''',encoding='utf-8')
        cfg=json.loads((BASE/'config/upstream.json').read_text(encoding='utf-8'))
        cfg['expected_main_isolates']=3;cfg['reference_fasta']=str(self.ref)
        for tool in ['Snippy','QUAST','BUSCO']:cfg['tools'][tool]['command']=[sys.executable,str(mock),tool]
        config=self.root/'mock_config.json';config.write_text(json.dumps(cfg),encoding='utf-8')
        output=self.root/'mock_run'
        for sid in self.ids:
            write_fasta(output/'assembly'/sid/'contigs.fasta',self.reference)
            for mate in [1,2]:
                path=output/'trim'/sid/f'{sid}.R{mate}.paired.fastq.gz';path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(b'test interface only')
        for stage in ['snippy','assess']:
            command=[sys.executable,str(BASE/'upstream/workflow.py'),stage,'--config',str(config),'--manifest',str(self.manifest),'--output',str(output)]
            result=subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
        for sid in self.ids:
            self.assertTrue((output/'logs'/f'{sid}_snippy.completed.json').exists())
            self.assertTrue((output/'logs'/f'{sid}_busco.completed.json').exists())

if __name__=='__main__':unittest.main(verbosity=2)

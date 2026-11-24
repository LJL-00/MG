#!/usr/bin/env python3
"""Run new read QC, assembly, assessment, Snippy filtering and tree analyses.

The supplied server inventory controls new execution versions. Analysis
parameters follow the study methods.
"""
from pathlib import Path
import argparse,json,re,sys
from common import BASE,Runner,check_tool,ensure_files,load_manifest,read_config,save_json,sha256
from strict_core import build_core,inspect_tree_alignment

def paired(runner,sid):
    folder=runner.output/'trim'/sid
    return folder/(sid+'.R1.paired.fastq.gz'),folder/(sid+'.R2.paired.fastq.gz')

def fastqc_stem(path):
    name=Path(path).name
    if name.endswith('.gz'):name=name[:-3]
    return re.sub(r'\.(fastq|fq)$','',name)+'_fastqc.zip'

def trim(runner,rows):
    adapters=Path(runner.config['adapters_fasta']).expanduser()
    qc_directories=[];qc_inputs=[]
    if not runner.dry_run:
        ensure_files([adapters],'Trimmomatic adapters')
        ensure_files([r[k] for r in rows for k in ['r1','r2']],'FASTQ manifest')
    for row in rows:
        sid=row['sample'];raw=[Path(row['r1']),Path(row['r2'])]
        if any(not str(x) or str(x)=='.' for x in raw):raise ValueError('Fill both r1/r2 paths for '+sid)
        rawqc=runner.output/'fastqc/raw'/sid
        rawout=[rawqc/fastqc_stem(p) for p in raw]
        runner.run(sid+'_fastqc_raw',runner.tool('FastQC')+['--threads',str(runner.config['threads']),'--outdir',str(rawqc)]+list(map(str,raw)),rawout,raw)
        qc_directories.append(rawqc);qc_inputs.extend(rawout)
        p1,p2=paired(runner,sid);u1=p1.with_name(sid+'.R1.unpaired.fastq.gz');u2=p2.with_name(sid+'.R2.unpaired.fastq.gz')
        command=runner.tool('Trimmomatic')+['PE','-threads',str(runner.config['threads']),'-phred33',str(raw[0]),str(raw[1]),str(p1),str(u1),str(p2),str(u2),f'ILLUMINACLIP:{adapters}:2:30:10','LEADING:20','TRAILING:20','SLIDINGWINDOW:4:20','MINLEN:50']
        runner.run(sid+'_trim',command,[p1,p2,u1,u2],raw+[adapters])
        trimmedqc=runner.output/'fastqc/trimmed'/sid
        trimmedout=[trimmedqc/fastqc_stem(p1),trimmedqc/fastqc_stem(p2)]
        runner.run(sid+'_fastqc_trimmed',runner.tool('FastQC')+['--threads',str(runner.config['threads']),'--outdir',str(trimmedqc),str(p1),str(p2)],trimmedout,[p1,p2])
        qc_directories.append(trimmedqc);qc_inputs.extend(trimmedout)
    report=runner.output/'multiqc/multiqc_report.html'
    runner.run('multiqc',runner.tool('MultiQC')+list(map(str,qc_directories))+['-o',str(report.parent)], [report],qc_inputs)

def assemble(runner,rows):
    for row in rows:
        if row['assemble']!='1':continue
        sid=row['sample'];r1,r2=paired(runner,sid);folder=runner.output/'assembly'/sid
        runner.run(sid+'_spades',runner.tool('SPAdes')+['--isolate','-1',str(r1),'-2',str(r2),'-o',str(folder),'-t',str(runner.config['threads']),'-m',str(runner.config['memory_gb']),'-k',','.join(map(str,runner.config['spades_kmers']))],[folder/'contigs.fasta'],[r1,r2],create_output_parents=False)

def assess(runner,rows):
    assemblies=[(r['sample'],runner.output/'assembly'/r['sample']/'contigs.fasta') for r in rows if r['assemble']=='1']
    if not assemblies:raise ValueError('No rows have assemble=1')
    for sid,contigs in assemblies:
        folder=runner.output/'quast'/sid
        runner.run(sid+'_quast',runner.tool('QUAST')+['--fungus','--min-contig','500','--threads',str(runner.config['threads']),'-o',str(folder),str(contigs)],[folder/'report.tsv'],[contigs])
        parent=runner.output/'busco';name=sid
        command=runner.tool('BUSCO')+['-i',str(contigs),'-m','genome','-l','fungi_odb10','-c',str(runner.config['threads']),'-o',name,'--out_path',str(parent)]
        dataset=runner.config.get('busco_dataset_path','')
        if dataset:
            # An exact local fungi_odb10 directory can replace the dataset name.
            command[command.index('-l')+1]=dataset
        else:
            cache=runner.config.get('busco_download_path','')
            if cache:command+=['--download_path',cache]
        expected=parent/name/f'short_summary.specific.fungi_odb10.{name}.txt'
        runner.run(sid+'_busco',command,[expected],[contigs],create_output_parents=False)

def snippy(runner,rows):
    reference=runner.config.get('reference_for_snippy') or runner.config['reference_fasta']
    for row in rows:
        sid=row['sample'];r1,r2=paired(runner,sid);folder=runner.output/'snippy'/sid
        command=runner.tool('Snippy')+['--cpus',str(runner.config['threads']),'--outdir',str(folder),'--prefix','snps','--ref',reference,'--R1',str(r1),'--R2',str(r2),'--minqual','100','--mincov','10','--minfrac','0.90']
        runner.run(sid+'_snippy',command,[folder/n for n in ['snps.aligned.fa','snps.vcf','snps.subs.vcf','reference/ref.fa']],[r1,r2,reference],create_output_parents=False)

def core(runner,manifest):
    if runner.dry_run:
        runner.plan.append({'stage':'core','manifest':str(manifest),'QUAL':100,'DP':10,'AF':.9,'indel_flank_bp':10,'outputs':['core/full_core.fasta.gz','core/SNPs_only.fasta.gz','core/pairwise_SNP_distances.tsv','core/core_report.json']});return
    if (runner.output/'core/core_report.json').exists():
        raise RuntimeError('Core output already exists. Choose a fresh --output directory, or call tree separately with the existing alignment.')
    ensure_files([runner.config['reference_fasta']],'Core reference')
    build_core(manifest,Path(runner.config['reference_fasta']),runner.output/'snippy',runner.output/'core',runner.config['expected_main_isolates'])

def tree(runner,rows,args):
    if args.new_tree_from_core:
        alignment=runner.output/'core/full_core.fasta.gz';origin='new_strict_core_all_sites'
    else:
        if not args.alignment:raise ValueError('tree requires --alignment PATH or explicit --new-tree-from-core. Choose a full-sites alignment with invariant columns retained.')
        alignment=Path(args.alignment).resolve();origin='user_supplied_complete_case_alignment'
    expected=[r['sample'] for r in rows if r['analysis_set']=='Main']
    if len(expected)!=runner.config['expected_main_isolates']:raise ValueError('Tree requires exactly the expected number of enabled Main isolates')
    if runner.dry_run:
        runner.plan.append({'stage':'tree_input_validation','input':str(alignment),'origin':origin,'requirements':'Expected Main sample IDs, same length, complete ACGT, invariant columns retained; SHA256 will be recorded'})
    else:
        ensure_files([alignment],'Tree input')
        report=inspect_tree_alignment(alignment,expected)
        save_json(runner.output/'tree/input_audit.json',report|{'origin':origin,'invariant_sites_required':True})
    iqtree=runner.tool('IQ-TREE')
    if args.check_modelfinder:
        prefix=runner.output/'tree/modelfinder'
        runner.run('tree_modelfinder',iqtree+['-s',str(alignment),'-st','DNA','-m','MF','-merit','BIC','-T',str(runner.config['threads']),'--prefix',str(prefix)],[prefix.with_suffix('.iqtree')],[alignment])
        if not runner.dry_run:
            text=prefix.with_suffix('.iqtree').read_text(encoding='utf-8',errors='replace')
            selected=re.search(r'Best-fit model according to BIC:\s*(\S+)',text)
            if not selected or selected.group(1)!='GTR+F+R4':
                raise RuntimeError('ModelFinder did not confirm the manuscript GTR+F+R4 model for this input. Read '+str(prefix.with_suffix('.iqtree'))+'; select the model for this input before final inference.')
    prefix=runner.output/'tree/iqtree'
    command=iqtree+['-s',str(alignment),'-st','DNA','-m','GTR+F+R4','-B','1000','-alrt','1000','-T',str(runner.config['threads']),'--prefix',str(prefix)]
    runner.run('tree_inference',command,[prefix.with_suffix('.treefile'),prefix.with_suffix('.iqtree'),prefix.with_suffix('.log')],[alignment])
    if not runner.dry_run:
        save_json(runner.output/'tree/run_settings.json',{'scope':'new inference using manuscript software and model settings','software':'IQ-TREE 3.1.1','model':'GTR+F+R4','ultrafast_bootstrap':1000,'SH_aLRT':1000,'alignment_origin':origin,'alignment_sha256':sha256(alignment),'modelfinder_BIC_checked':args.check_modelfinder,'rooting':'no outgroup supplied; interpret topology as unrooted','seed':'IQ-TREE default, recorded by IQ-TREE in its log','alignment_type':'full-sites'})

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('stage',choices=['versions','trim','assemble','assess','snippy','core','tree','all'])
    ap.add_argument('--config',type=Path,default=BASE/'config/upstream.json')
    ap.add_argument('--manifest',type=Path,help='Sample TSV; required except for the versions stage')
    ap.add_argument('--output',type=Path,default=BASE/'results/upstream')
    ap.add_argument('--alignment',help='Confirmed full-sites tree alignment, with all Main sample IDs and invariant columns')
    ap.add_argument('--new-tree-from-core',action='store_true',help='Explicitly generate a new inference from the strict core including invariant sites; includes invariant columns')
    ap.add_argument('--check-modelfinder',action='store_true',help='Run ModelFinder/BIC first and require GTR+F+R4 before final tree inference')
    ap.add_argument('--dry-run',action='store_true',help='Write commands without executing tools or claiming that versions/results were verified')
    ap.add_argument('--allow-missing-tools',action='store_true',help='For versions inventory only; report unavailable tools without failing')
    args=ap.parse_args();config=read_config(args.config)
    if args.stage=='versions':
        reports={}
        for name in config['tools']:
            try:reports[name]=check_tool(config,name)
            except RuntimeError as exc:
                reports[name]={'status':'UNAVAILABLE_OR_MISMATCH','message':str(exc)}
                if not args.allow_missing_tools:raise
        save_json(args.output/'logs/software_versions.json',reports);print(json.dumps(reports,indent=2));return
    if args.manifest is None:raise ValueError('Set --manifest to a sample TSV')
    rows=load_manifest(args.manifest)
    if args.stage=='all' and not args.alignment and not args.new_tree_from_core:
        raise ValueError('all requires a confirmed --alignment PATH, or --new-tree-from-core to select a new full-sites alignment explicitly')
    runner=Runner(args.output,config,args.dry_run)
    if args.stage in ['trim','all']:trim(runner,rows)
    if args.stage in ['assemble','all']:assemble(runner,rows)
    if args.stage in ['assess','all']:assess(runner,rows)
    if args.stage in ['snippy','all']:snippy(runner,rows)
    if args.stage in ['core','all']:core(runner,args.manifest)
    if args.stage in ['tree','all']:tree(runner,rows,args)
    save_json(args.output/('command_plan.json' if args.dry_run else 'command_record.json'),{'stage':args.stage,'mode':'PLAN_ONLY_NOT_EXECUTED' if args.dry_run else 'EXECUTED','configuration':str(args.config.resolve()),'manifest':str(args.manifest.resolve()),'commands':runner.plan})
    print(json.dumps({'stage':args.stage,'dry_run':args.dry_run,'commands':len(runner.plan),'output':str(runner.output)}))

if __name__=='__main__':
    try:main()
    except (OSError,ValueError,RuntimeError) as exc:
        print('ERROR: '+str(exc),file=sys.stderr);raise SystemExit(2)


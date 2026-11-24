#!/usr/bin/env python3
"""Launch new STRUCTURE 2.3.4 chains with the manuscript 10k/20k settings."""
from pathlib import Path
import argparse,json,re,subprocess,sys,shutil
from common import BASE,check_tool,executable,read_config,save_json,sha256

def parameters(text,values):
    for key,value in values.items():
        pattern=rf'(?m)^(#define\s+{re.escape(key)}\s+)\S+'
        if not re.search(pattern,text):raise ValueError('Missing template setting '+key)
        text=re.sub(pattern,lambda match:match[1]+str(value),text)
    return text

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,default=BASE/'config/upstream.json')
    ap.add_argument('--input',type=Path,required=True,help='Labelled, haploid-encoded STRUCTURE matrix')
    ap.add_argument('--templates',type=Path,default=BASE/'config/structure')
    ap.add_argument('--output',type=Path,default=BASE/'results/new_structure')
    ap.add_argument('--dry-run',action='store_true')
    args=ap.parse_args();config=read_config(args.config)
    source=args.input.resolve();lines=[x.split() for x in source.read_text(encoding='utf-8').splitlines() if x.strip()]
    if not lines or len({r[0] for r in lines})!=len(lines) or len({len(r) for r in lines})!=1:raise ValueError('Expected unique labelled rows with equal marker counts')
    n_individuals=len(lines)
    n_loci=len(lines[0])-1
    if not n_loci or not all(x in ['1','2','3','4','-9'] for r in lines for x in r[1:]):raise ValueError('Expected one encoded allele per marker, with -9 for missing')
    if args.output.exists():raise ValueError('Choose a new output directory for new chains')
    tool_version=None if args.dry_run else check_tool(config,'STRUCTURE')
    main_text=(args.templates/'mainparams').read_text(encoding='utf-8')
    extra_text=(args.templates/'extraparams').read_text(encoding='utf-8');jobs=[]
    args.output.mkdir(parents=True)
    shutil.copy2(source,args.output/'input.structure')
    for k in range(1,9):
        for rep in range(1,6):
            folder=args.output.resolve()/f'K{k:02d}_rep{rep:02d}';folder.mkdir(parents=True)
            seed=20261017+1000*k+rep
            main=parameters(main_text,{'MAXPOPS':k,'BURNIN':10000,'NUMREPS':20000,'NUMINDS':n_individuals,'NUMLOCI':n_loci,'PLOIDY':1,'MISSING':-9,'ONEROWPERIND':1,'LABEL':1,'POPDATA':0,'POPFLAG':0,'LOCDATA':0,'PHENOTYPE':0,'EXTRACOLS':0,'MARKERNAMES':0,'MAPDISTANCES':0,'INFILE':'../input.structure','OUTFILE':'result'})
            extra=parameters(extra_text,{'NOADMIX':0,'LINKAGE':0,'USEPOPINFO':0,'LOCPRIOR':0,'LOCISPOP':0,'FREQSCORR':1,'INFERALPHA':1,'POPALPHAS':0,'ALPHA':1,'INFERLAMBDA':0,'LAMBDA':1,'RANDOMIZE':0,'SEED':seed,'PRINTNET':0,'PRINTLAMBDA':0,'PRINTQSUM':0,'PRINTLIKES':1,'UPDATEFREQ':500,'ECHODATA':0,'ANCESTDIST':0,'COMPUTEPROB':1})
            (folder/'mainparams').write_text(main,encoding='utf-8');(folder/'extraparams').write_text(extra,encoding='utf-8')
            command=executable(config,'STRUCTURE')+['-m','mainparams','-e','extraparams','-D',str(seed)]
            record={'scope':'new chain','K':k,'replicate':rep,'burnin':10000,'sampling':20000,'seed':seed,'command':command,'directory':str(folder),'dry_run':args.dry_run}
            if not args.dry_run:
                with (folder/'console.log').open('w',encoding='utf-8') as handle:
                    result=subprocess.run(command,cwd=folder,stdout=handle,stderr=subprocess.STDOUT,input='\n',text=True)
                record['exit_code']=result.returncode
                if result.returncode or not (folder/'result_f').is_file():
                    save_json(folder/'execution.json',record);raise RuntimeError('Chain failed; read '+str(folder/'console.log'))
            save_json(folder/'execution.json',record);jobs.append(record)
    save_json(args.output/'run_plan.json',{'scope':'STRUCTURE execution','input_sha256':sha256(source),'software':tool_version,'samples':n_individuals,'markers':n_loci,'burnin':10000,'sampling':20000,'K':list(range(1,9)),'runs_per_K':5,'jobs':jobs})
    print(json.dumps({'dry_run':args.dry_run,'chains':len(jobs),'burnin':10000,'sampling':20000}))

if __name__=='__main__':main()

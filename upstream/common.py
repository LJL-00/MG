"""Configuration, input validation and recorded command execution."""
from pathlib import Path
import csv, datetime, hashlib, json, re, subprocess

BASE = Path(__file__).resolve().parents[1]

def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def save_json(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')

def load_manifest(path):
    path=Path(path).resolve()
    with path.open(encoding='utf-8-sig',newline='') as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    required={'sample','r1','r2','analysis_set','assemble','enabled'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError('Manifest needs columns: '+', '.join(sorted(required)))
    seen=set();selected=[]
    for row in rows:
        sid=row['sample'].strip()
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',sid) or sid in seen:
            raise ValueError('Unsafe or duplicate sample identifier: '+sid)
        seen.add(sid)
        for field in ['assemble','enabled']:
            if row[field] not in ['0','1']:raise ValueError(field+' must be 0 or 1 for '+sid)
        if row['analysis_set'] not in ['Main','ANI assessment']:
            raise ValueError('analysis_set must be Main or ANI assessment for '+sid)
        for field in ['r1','r2']:
            raw=row[field].strip()
            row[field]=str((path.parent/raw).resolve()) if raw and not Path(raw).is_absolute() else raw
        if row['enabled']=='1':selected.append(row)
    if not selected:raise ValueError('No enabled samples in manifest')
    return selected

def read_config(path):
    value=json.loads(Path(path).read_text(encoding='utf-8'))
    for name,spec in value['tools'].items():
        if not isinstance(spec['command'],list) or not all(isinstance(x,str) for x in spec['command']):
            raise ValueError('tools.'+name+'.command must be a JSON string array')
    return value

def executable(config,name):return config['tools'][name]['command'][:]

def check_tool(config,name):
    spec=config['tools'][name]
    command=spec['command']+spec['version_args']
    try:
        result=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,errors='replace',timeout=60)
    except (OSError,subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f'{name}: cannot run {command!r}. Install its environment or edit config/upstream.json tools.{name}.command. {exc}') from exc
    output=result.stdout
    found=re.search(spec['version_regex'],output,re.I)
    if not found:raise RuntimeError(f'{name}: version was not found in output: {output[:1000]}')
    version=found.group(1).rstrip('.')
    required=spec.get('required_version')
    if required and version!=required:
        raise RuntimeError(f'{name}: configuration requires {required}, executable reports {version}. Select the matching executable in config/upstream.json.')
    if result.returncode and not spec.get('version_nonzero_ok',False):
        raise RuntimeError(f'{name}: version command failed with exit {result.returncode}: {output[:1000]}')
    report={'software':name,'required_version':required,'detected_version':version,'command':command,'output':output,'status':'MATCH' if required else 'RECORDED_WITHOUT_VERSION_TARGET'}
    if spec.get('preflight_args'):
        preflight=spec['command']+spec['preflight_args']
        try:
            result=subprocess.run(preflight,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,errors='replace',timeout=60)
        except (OSError,subprocess.TimeoutExpired) as exc:
            raise RuntimeError(f'{name}: dependency preflight could not run: {exc}') from exc
        report['dependency_check']={'command':preflight,'exit_code':result.returncode,'output':result.stdout}
        if result.returncode:
            raise RuntimeError(f'{name}: dependency preflight failed. Snippy requires BCFtools at least 1.7. Activate its compatible environment. Output: {result.stdout[:2000]}')
    return report

def ensure_files(paths,context):
    missing=[str(p) for p in paths if not Path(p).is_file() or Path(p).stat().st_size==0]
    if missing:raise FileNotFoundError(context+' needs nonempty input files. Set their paths in the manifest/configuration:\n'+'\n'.join(missing[:20]))

class Runner:
    def __init__(self, output, config, dry_run=False):
        self.output=Path(output).resolve();self.config=config;self.dry_run=dry_run;self.plan=[];self.versions={}
    def tool(self,name):
        if not self.dry_run and name not in self.versions:
            self.versions[name]=check_tool(self.config,name)
            save_json(self.output/'logs/software_versions.json',self.versions)
        return executable(self.config,name)
    def run(self,name,command,outputs,inputs=(),create_output_parents=True):
        outputs=[Path(p) for p in outputs];inputs=[Path(p) for p in inputs]
        info={'name':name,'command':[str(x) for x in command],'outputs':[str(x) for x in outputs],'inputs':[str(x) for x in inputs]}
        self.plan.append(info)
        if self.dry_run:return
        ensure_files(inputs,name)
        matching=[(tool,record) for tool,record in self.versions.items() if info['command'][:len(self.config['tools'][tool]['command'])]==self.config['tools'][tool]['command']]
        if not matching:raise RuntimeError(f'{name}: external command has no checked software version. Obtain it through Runner.tool before execution.')
        tool,record=max(matching,key=lambda item:len(self.config['tools'][item[0]]['command']))
        software={'name':tool,'required_version':record['required_version'],'detected_version':record['detected_version'],'version_output_sha256':hashlib.sha256(record['output'].encode()).hexdigest()}
        if 'dependency_check' in record:software['dependency_output_sha256']=hashlib.sha256(record['dependency_check']['output'].encode()).hexdigest()
        info['software']=software
        signature=hashlib.sha256(json.dumps({'command':info['command'],'software':software,'inputs':[(str(p),p.stat().st_size,p.stat().st_mtime_ns) for p in inputs]},sort_keys=True).encode()).hexdigest()
        marker=self.output/'logs'/(name+'.completed.json')
        if marker.exists():
            previous=json.loads(marker.read_text(encoding='utf-8'))
            if previous['signature']==signature and all(p.is_file() and p.stat().st_size for p in outputs):return
            raise RuntimeError(f'{name}: existing completion record differs from this command or inputs. Choose a new --output directory.')
        if any(p.exists() for p in outputs):raise RuntimeError(f'{name}: outputs exist without matching completion record. Choose a new --output directory.')
        if create_output_parents:
            for p in outputs:p.parent.mkdir(parents=True,exist_ok=True)
        log=self.output/'logs'/(name+'.log');log.parent.mkdir(parents=True,exist_ok=True)
        started=datetime.datetime.now(datetime.timezone.utc).isoformat()
        with log.open('w',encoding='utf-8') as handle:
            handle.write(json.dumps(info['command'])+'\n');handle.flush()
            result=subprocess.run(info['command'],stdout=handle,stderr=subprocess.STDOUT)
        if result.returncode:raise RuntimeError(f'{name} exited {result.returncode}. Read {log}')
        ensure_files(outputs,name+' output validation')
        save_json(marker,info|{'signature':signature,'started_UTC':started,'exit_code':0})


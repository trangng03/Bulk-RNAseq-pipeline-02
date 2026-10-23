#!/usr/bin/env python3
"""Stage-based paired-end RNA-seq runner. Standard library only; no installs."""
import argparse, csv, gzip, hashlib, json, math, os, re, shlex, shutil, subprocess, sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parent.parent

def die(message):
    raise SystemExit(message)

def open_text(path):
    return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)

def attrs(text):
    return dict(re.findall(r'(?:^|;\s*)(\w+)\s+"([^"]*)"',text))

def prepare_gtf(src, dest):
    """Make BED12 transcript models, tx2gene and symbols from one annotation."""
    dest.mkdir(parents=True,exist_ok=True)
    transcripts={}; txmap={}; symbols={}; cds=defaultdict(list)
    plain=dest/'annotation.gtf'
    with open_text(src) as fi, plain.open('w') as fo:
        for line in fi:
            fo.write(line)
            if line.startswith('#'): continue
            fields=line.rstrip('\n').split('\t')
            if len(fields)!=9: die('Malformed GTF row: expected nine tab-separated columns')
            a=attrs(fields[8]); gid=a.get('gene_id'); tid=a.get('transcript_id')
            if gid and (gid not in symbols or not symbols[gid]): symbols[gid]=a.get('gene_name','')
            if tid and gid:
                if tid in txmap and txmap[tid]!=gid: die('Transcript maps to multiple genes: '+tid)
                txmap[tid]=gid
            if fields[2]=='CDS' and tid:
                cds[tid].append((int(fields[3])-1,int(fields[4])))
            if fields[2]=='exon' and tid:
                key=(fields[0],fields[6]); entry=transcripts.setdefault(tid,[key,[]])
                if entry[0]!=key: die('Transcript has inconsistent chromosome/strand: '+tid)
                entry[1].append((int(fields[3])-1,int(fields[4])))
    if not txmap or not transcripts: die('No transcript/exon models found in GTF')
    bed=dest/'transcripts.bed'
    with bed.open('w') as fo:
        for tid,(key,exons) in sorted(transcripts.items()):
            merged=[]
            for start,end in sorted(set(exons)):
                if merged and start<=merged[-1][1]: merged[-1]=(merged[-1][0],max(end,merged[-1][1])) #overlapped intervals
                else: merged.append((start,end))
            start,end=merged[0][0],merged[-1][1]
            thick_start=min(s for s,e in cds[tid]) if cds[tid] else start
            thick_end=max(e for s,e in cds[tid]) if cds[tid] else start
            row=[key[0],str(start),str(end),tid,'0',key[1],str(thick_start),str(thick_end),'0',str(len(merged)),','.join(str(e-s) for s,e in merged)+',',','.join(str(s-start) for s,e in merged)+',']
            fo.write('\t'.join(row)+'\n')
    for filename,header,rows in [('tx2gene.tsv',['transcript_id','gene_id'],sorted(txmap.items())),('gene_symbols.tsv',['gene_id','gene_name'],sorted(symbols.items()))]:
        with (dest/filename).open('w') as fo:
            writer=csv.writer(fo,delimiter='\t',lineterminator='\n'); writer.writerow(header); writer.writerows(rows)
    return plain,bed

def sampled_lengths(path,n=1000):
    lengths=set()
    with open_text(path) as f:
        for i in range(n):
            header=f.readline()
            if not header: break
            seq=f.readline().rstrip('\r\n'); plus=f.readline(); qual=f.readline().rstrip('\r\n')
            if not header.startswith('@') or not plus.startswith('+') or len(seq)!=len(qual):
                die('Malformed or multiline FASTQ in '+str(path))
            lengths.add(len(seq))
    if not lengths: die('Empty FASTQ: '+str(path))
    return sorted(lengths)

def infer(path):
    text=Path(path).read_text()
    def fraction(label):
        m=re.search(r'Fraction of reads explained by "'+re.escape(label)+r'":\s*([0-9.]+)',text)
        return float(m.group(1)) if m else None
    f=fraction('1++,1--,2+-,2-+'); r=fraction('1+-,1-+,2++,2--')
    if f is None or r is None: return (None,None,'unclear')
    usable=f+r
    if usable<0.5: return(f,r,'unclear')
    ratio=f/usable
    if ratio>=0.8: strand='forward'
    elif ratio<=0.2: strand='reverse'
    elif 0.4<=ratio<=0.6: strand='unstranded'
    else: strand='unclear'
    return f,r,strand

class Pipeline:
    def __init__(self,path,force=False):
        self.config_path=Path(path).resolve(); self.base=self.config_path.parent
        self.c=json.loads(self.config_path.read_text()); self.force=force
        c=self.c
        # Run only RSeQC in its own Conda environment (no activation required).
        # Optional config key: "rseqc_env": "rseqc"; RSEQC_ENV is a fallback.
        c.setdefault('rseqc_env', os.environ.get('RSEQC_ENV', 'rseqc'))
        if not isinstance(c['rseqc_env'], str) or not c['rseqc_env'].strip():
            die('rseqc_env must be a nonempty Conda environment name')
        def resolve(p): return str((self.base/p).resolve()) if not Path(p).is_absolute() else p
        self.resolve=resolve
        for k in ['samples','results_dir','star_index','kallisto_index','gtf','gmt']:
            c[k]=resolve(c[k])
        if c.get('existing_star_dir'): c['existing_star_dir']=resolve(c['existing_star_dir'])
        self.out=Path(c['results_dir'])
        for d in ['00_reference','01_fastqc','02_star','03_rseqc','04_kallisto','05_deseq2','06_gsea','07_rmats','08_multiqc','logs']:
            (self.out/d).mkdir(parents=True,exist_ok=True)
        if c['case']==c['control']: die('case and control must differ')
        if not isinstance(c['threads'],int) or c['threads']<1: die('threads must be a positive integer')
        if c['strand'] not in ['auto','unstranded','forward','reverse']: die('Invalid strand setting')
        with open(c['samples']) as f:
            reader=csv.DictReader(f,delimiter='\t')
            if reader.fieldnames!=['Sample','Condition','Read1','Read2']: die('Sample sheet header must be Sample, Condition, Read1, Read2')
            self.samples=list(reader)
        seen=set()
        for s in self.samples:
            if None in s or any(not v for v in s.values()): die('Invalid sample row')
            if not re.fullmatch(r'[A-Za-z0-9_.-]+',s['Sample']) or s['Sample'] in seen: die('Invalid/duplicate sample ID')
            seen.add(s['Sample'])
            for key in ['Read1','Read2']:
                s[key]=resolve(s[key]); self.file(s[key])
        groups=Counter(s['Condition'] for s in self.samples)
        if set(groups)!={c['case'],c['control']} or min(groups.values())<2: die('Need exactly case/control and at least 2 biological replicates each')
        self.file(c['gtf'])
        with (self.out/'input_manifest.tsv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=['Sample','Condition','Read1','Read2'],delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(self.samples)
        # Completion records are invalidated when config, sample sheet or input metadata changes.
        metadata=[]
        for p in [c['gtf'],c['samples']]+[s[k] for s in self.samples for k in ['Read1','Read2']]:
            st=Path(p).stat(); metadata.append((p,st.st_size,st.st_mtime_ns))
        code_hashes={str(p.name):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'scripts/pipeline.py',ROOT/'scripts/analyze.R']}
        for extra in [c['gmt'],c['kallisto_index'],str(Path(c['star_index'])/'genomeParameters.txt')]:
            if Path(extra).is_file():
                st=Path(extra).stat(); metadata.append((extra,st.st_size,st.st_mtime_ns))
        if c.get('existing_star_dir'):
            for sample in self.samples:
                extra=str(Path(c['existing_star_dir'])/sample['Sample']/'Aligned.sortedByCoord.out.bam')
                if Path(extra).is_file():
                    st=Path(extra).stat(); metadata.append((extra,st.st_size,st.st_mtime_ns))
        self.fingerprint=hashlib.sha256(json.dumps([c,self.samples,metadata,code_hashes],sort_keys=True).encode()).hexdigest()
        (self.out/'resolved_config.json').write_text(json.dumps(c,indent=2)+'\n')
        self.gtf=self.out/'00_reference/annotation.gtf'; self.bed=self.out/'00_reference/transcripts.bed'
    def file(self,p):
        if not Path(p).is_file(): die('File missing: '+str(p))
    def need(self,*commands):
        missing=[x for x in commands if not shutil.which(x)]
        with (self.out/'logs/executables.tsv').open('a') as f:
            for command in commands: f.write(command+'\t'+str(shutil.which(command))+'\n')
        if missing: die('Missing commands in active environment: '+', '.join(missing)+'. Nothing was installed.')
    def run(self,cmd,name,output=None):
        self.need(cmd[0]); cmd=list(map(str,cmd))
        logfile=self.out/'logs'/name
        print(shlex.join(cmd),flush=True)
        with logfile.open('w') as log:
            log.write(shlex.join(cmd)+'\n'); log.flush()
            if output:
                with open(output,'w') as stdout: result=subprocess.run(cmd,stdout=stdout,stderr=log)
            else: result=subprocess.run(cmd,stdout=log,stderr=log)
        if result.returncode: die('Command failed; see '+str(logfile))
    def run_rseqc(self,cmd,name,output=None):
        # Argument lists preserve spaces in BAM/BED/output paths.
        wrapped=['conda','run','--no-capture-output','-n',self.c['rseqc_env']]+list(cmd)
        self.run(wrapped,name,output)
    def check_rseqc(self):
        tools=['infer_experiment.py']
        check=(
            'import shutil, sys; '
            'missing=[x for x in sys.argv[1:] if shutil.which(x) is None]; '
            'print("Missing RSeQC environment commands: " + ", ".join(missing) '
            'if missing else "RSeQC environment commands available", flush=True); '
            'sys.exit(1 if missing else 0)'
        )
        print('Checking RSeQC in Conda environment: '+self.c['rseqc_env'],flush=True)
        self.run_rseqc(['python','-c',check]+tools,'rseqc.dependencies.log')
    def bam(self,s):
        base=Path(self.c.get('existing_star_dir') or self.out/'02_star')
        return base/s['Sample']/'Aligned.sortedByCoord.out.bam'
    def strand(self):
        explicit=self.c['strand']
        if explicit!='auto': return explicit
        values=[infer(self.out/'03_rseqc'/s['Sample']/'infer_experiment.txt')[2] for s in self.samples]
        if 'unclear' in values or len(set(values))!=1:
            die('Strandedness ambiguous or inconsistent. Review 03_rseqc/strandedness.tsv and library protocol; set strand explicitly only after review.')
        return values[0]
    def qc(self):
        self.need('fastqc','samtools','multiqc','conda')
        self.check_rseqc()
        if not self.c.get('existing_star_dir'):
            self.need('STAR'); self.file(Path(self.c['star_index'])/'genomeParameters.txt')
        prepare_gtf(self.c['gtf'],self.out/'00_reference')
        summary=[]; read_lengths=[]
        for s in self.samples:
            name=s['Sample']
            for mate in ['Read1','Read2']:
                read_lengths.append([name,mate,','.join(map(str,sampled_lengths(s[mate])))])
            fqdir=self.out/'01_fastqc'/name; fqdir.mkdir(exist_ok=True)
            self.run(['fastqc','-t',str(min(2,self.c['threads'])),'-o',fqdir,s['Read1'],s['Read2']],name+'.fastqc.log')
            bam=self.bam(s)
            if not self.c.get('existing_star_dir'):
                d=bam.parent; d.mkdir(exist_ok=True)
                cmd=['STAR','--runThreadN',str(self.c['threads']),'--genomeDir',self.c['star_index'],'--sjdbGTFfile',self.gtf,'--readFilesIn',s['Read1'],s['Read2'],'--outSAMtype','BAM','SortedByCoordinate','--quantMode','GeneCounts','--outFileNamePrefix',str(d)+'/']
                compressed=[p.endswith('.gz') for p in (s['Read1'],s['Read2'])]
                if any(compressed):
                    if not all(compressed): die('Both mates must have the same compression type')
                    self.need('zcat'); cmd+=['--readFilesCommand','zcat']
                self.run(cmd,name+'.star.log')
            self.file(bam)
            self.run(['samtools','quickcheck','-v',bam],name+'.bam_check.log')
            self.run(['samtools','index','-@',str(self.c['threads']),bam],name+'.index.log')
            d=self.out/'03_rseqc'/name; d.mkdir(exist_ok=True)
            self.run(['samtools','flagstat',bam],name+'.flagstat.log',d/'flagstat.txt')
            self.run_rseqc(['infer_experiment.py','-i',bam,'-r',self.bed,'-s','200000'],name+'.infer.log',d/'infer_experiment.txt')
            f,r,decision=infer(d/'infer_experiment.txt'); summary.append([name,f,r,decision])
        with (self.out/'03_rseqc/strandedness.tsv').open('w') as f:
            w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(['Sample','forward_fraction','reverse_fraction','suggestion']); w.writerows(summary)
        with (self.out/'03_rseqc/read_lengths.tsv').open('w') as f:
            w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(['Sample','mate','lengths_in_first_1000_reads']); w.writerows(read_lengths)
        self.multiqc()
    def quant(self):
        self.need('kallisto'); self.file(self.c['kallisto_index']); strand=self.strand()
        for s in self.samples:
            d=self.out/'04_kallisto'/s['Sample']; d.mkdir(exist_ok=True)
            cmd=['kallisto','quant','-i',self.c['kallisto_index'],'-o',d,'-t',str(self.c['threads'])]
            if strand=='forward': cmd+=['--fr-stranded']
            elif strand=='reverse': cmd+=['--rf-stranded']
            self.run(cmd+[s['Read1'],s['Read2']],s['Sample']+'.kallisto.log')
        (self.out/'04_kallisto/strand_used.txt').write_text(strand+'\n')
    def analysis(self):
        self.need('Rscript'); self.file(self.c['gmt'])
        if not (self.out/'00_reference/tx2gene.tsv').exists(): prepare_gtf(self.c['gtf'],self.out/'00_reference')
        self.run(['Rscript',ROOT/'scripts/analyze.R',self.out/'input_manifest.tsv',self.out,self.c['gmt'],self.c['case'],self.c['control']], 'analysis.log')
    def rmats(self):
        c=self.c; strand=self.strand(); command=c['rmats_command']
        if not isinstance(command,list) or not command: die('rmats_command must be a nonempty JSON array')
        if not isinstance(c['read_length'],int) or c['read_length']<=0: die('Set read_length to your actual read length')
        if not c['variable_read_length']:
            for sample in self.samples:
                for mate in ['Read1','Read2']:
                    lengths=sampled_lengths(sample[mate])
                    if lengths!=[c['read_length']]: die('FASTQ read length differs from read_length for '+sample['Sample']+' '+mate+': '+str(lengths)+'. Verify configuration or enable variable_read_length for genuinely variable reads.')
        if not (0<c['rmats_fdr']<1 and 0<=c['rmats_delta_psi']<=1): die('Invalid rMATS filtering thresholds')
        if not self.gtf.exists(): prepare_gtf(c['gtf'],self.out/'00_reference')
        d=self.out/'07_rmats'; result=d/'raw'; tmp=d/'tmp'
        if result.exists() and any(result.iterdir()): die('rMATS raw directory already contains results. Use a new results_dir or move 07_rmats/raw and tmp before rerunning.')
        result.mkdir(exist_ok=True); tmp.mkdir(exist_ok=True)
        for group,key in [(c['case'],'case'),(c['control'],'control')]:
            bams=[str(self.bam(s)) for s in self.samples if s['Condition']==group]
            for b in bams:
                self.file(b)
                if ',' in b: die('rMATS BAM paths cannot contain commas')
            (d/(key+'_bams.txt')).write_text(','.join(bams)+'\n')
        lib={'unstranded':'fr-unstranded','forward':'fr-secondstrand','reverse':'fr-firststrand'}[strand]
        cmd=command+['--b1',str(d/'case_bams.txt'),'--b2',str(d/'control_bams.txt'),'--gtf',str(self.gtf),'-t','paired','--readLength',str(c['read_length']),'--libType',lib,'--nthread',str(c['threads']),'--od',str(result),'--tmp',str(tmp)]
        if c['variable_read_length']: cmd+=['--variable-read-length']
        if c['allow_clipping']: cmd+=['--allow-clipping']
        self.run(cmd,'rmats.log')
        summary=[]
        for p in sorted(result.glob('*.MATS.J*.txt')):
            with p.open() as f:
                reader=csv.DictReader(f,delimiter='\t'); fields=reader.fieldnames; rows=list(reader)
            significant=[]
            for row in rows:
                try: fdr=float(row['FDR']); delta=float(row['IncLevelDifference'])
                except (ValueError,KeyError): continue
                if math.isfinite(fdr) and math.isfinite(delta) and fdr<c['rmats_fdr'] and abs(delta)>=c['rmats_delta_psi']: significant.append(row)
            with (d/('significant_'+p.name)).open('w') as f:
                w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(significant)
            summary.append([p.name,len(rows),len(significant)])
        if not summary: die('No rMATS event tables found; inspect logs/rmats.log')
        with (d/'event_summary.tsv').open('w') as f:
            w=csv.writer(f,delimiter='\t',lineterminator='\n'); w.writerow(['table','events','significant_events']); w.writerows(summary)
        (d/'comparison.txt').write_text('Group 1 = '+c['case']+'\nGroup 2 = '+c['control']+'\nIncLevelDifference = inclusion(case) - inclusion(control). Interpret per event type, not as gene fold change.\n')
    def multiqc(self):
        self.run(['multiqc',self.out,'-o',self.out/'08_multiqc','--force','--ignore','08_multiqc'],'multiqc.log')
    def stage(self,name):
        marker=self.out/'logs'/('.'+name+'.done.json')
        if marker.exists() and not self.force:
            previous=json.loads(marker.read_text())
            if previous.get('fingerprint')!=self.fingerprint: die('Inputs/config changed since '+name+'. Use a fresh results_dir or rerun with --force (rMATS requires fresh raw/tmp folders).')
            print('Already completed: '+name+'; use --force to rerun',flush=True); return
        getattr(self,name)()
        marker.write_text(json.dumps({'fingerprint':self.fingerprint})+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['qc','quant','analysis','rmats','report','all'])
    p.add_argument('config'); p.add_argument('--force',action='store_true')
    a=p.parse_args()
    try:
        pipeline=Pipeline(a.config,a.force)
        if a.stage=='report': pipeline.multiqc()
        else:
            for stage in (['qc','quant','analysis','rmats'] if a.stage=='all' else [a.stage]): pipeline.stage(stage)
            if a.stage=='all': pipeline.multiqc()
        print('Done. Results: '+str(pipeline.out))
    except (OSError,KeyError,ValueError) as e:
        die('Pipeline error: '+str(e))

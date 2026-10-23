import csv, importlib.util, tempfile, unittest
from pathlib import Path
spec=importlib.util.spec_from_file_location('pipeline',Path(__file__).with_name('pipeline.py'))
p=importlib.util.module_from_spec(spec); spec.loader.exec_module(p)

class Helpers(unittest.TestCase):
 def test_annotation(self):
  with tempfile.TemporaryDirectory() as tmp:
   d=Path(tmp); gtf=d/'source.gtf'
   gtf.write_text('chr1\tx\tgene\t101\t400\t.\t-\t.\tgene_id "G1"; gene_name "NAME1";\n'
    'chr1\tx\texon\t101\t150\t.\t-\t.\tgene_id "G1"; transcript_id "T1.2";\n'
    'chr1\tx\texon\t301\t400\t.\t-\t.\tgene_id "G1"; transcript_id "T1.2";\n')
   plain,bed=p.prepare_gtf(gtf,d/'ref')
   fields=bed.read_text().strip().split('\t')
   self.assertEqual(fields,['chr1','100','400','T1.2','0','-','100','100','0','2','50,100,','0,200,'])
   self.assertIn('T1.2\tG1', (d/'ref/tx2gene.tsv').read_text())
   self.assertIn('G1\tNAME1', (d/'ref/gene_symbols.tsv').read_text())
 def test_infer(self):
  with tempfile.TemporaryDirectory() as tmp:
   f=Path(tmp)/'infer.txt'
   for forward,reverse,expected in [(0.04,0.9,'reverse'),(0.9,0.04,'forward'),(0.45,0.45,'unstranded'),(0.3,0.6,'unclear'),(0.05,0.1,'unclear')]:
    f.write_text(f'Fraction of reads explained by "1++,1--,2+-,2-+": {forward}\nFraction of reads explained by "1+-,1-+,2++,2--": {reverse}\n')
    self.assertEqual(p.infer(f)[2],expected)
   f.write_text('unexpected output'); self.assertEqual(p.infer(f)[2],'unclear')
 def test_stage_commands(self):
  import json
  with tempfile.TemporaryDirectory() as tmp:
   d=Path(tmp); (d/'index').mkdir(); (d/'index/genomeParameters.txt').write_text('index')
   (d/'kallisto.idx').write_text('index'); (d/'pathways.gmt').write_text('P\tdescription\tNAME1\n')
   (d/'genes.gtf').write_text('chr1\tx\texon\t101\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "NAME1";\n')
   rows=[]
   for condition in ['treated','control']:
    for i in range(2):
     sample=condition+str(i)
     mates=[]
     for mate in [1,2]:
      f=d/(sample+f'_{mate}.fastq'); f.write_text('@read\nACGT\n+\nIIII\n'); mates.append(str(f))
     rows.append([sample,condition]+mates)
   with (d/'samples.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['Sample','Condition','Read1','Read2']); w.writerows(rows)
   config={'samples':'samples.tsv','results_dir':'results','case':'treated','control':'control','threads':2,'star_index':'index','kallisto_index':'kallisto.idx','gtf':'genes.gtf','gmt':'pathways.gmt','strand':'auto','existing_star_dir':'','read_length':4,'variable_read_length':False,'allow_clipping':False,'rmats_command':['rmats.py'],'rmats_fdr':0.05,'rmats_delta_psi':0.1}
   cp=d/'config.json'; cp.write_text(json.dumps(config)); runner=p.Pipeline(cp)
   commands=[]
   runner.need=lambda *args: None
   def fake_run(cmd,name,output=None):
    cmd=list(map(str,cmd)); commands.append(cmd)
    if cmd[0]=='STAR':
     prefix=Path(cmd[cmd.index('--outFileNamePrefix')+1]); (prefix/'Aligned.sortedByCoord.out.bam').write_text('fake')
    if output:
     Path(output).write_text('Fraction of reads explained by "1++,1--,2+-,2-+": 0.04\nFraction of reads explained by "1+-,1-+,2++,2--": 0.90\n' if cmd[0]=='infer_experiment.py' else 'QC')
    if cmd[0]=='rmats.py':
     od=Path(cmd[cmd.index('--od')+1]); (od/'SE.MATS.JC.txt').write_text('ID\tFDR\tIncLevelDifference\n1\t0.001\t-0.2\n2\t0.9\t0.4\n')
   runner.run=fake_run
   runner.qc(); runner.quant(); runner.analysis(); runner.rmats()
   kallisto=[x for x in commands if x[0]=='kallisto']; self.assertEqual(len(kallisto),4)
   self.assertTrue(all('--rf-stranded' in x for x in kallisto))
   rmats=[x for x in commands if x[0]=='rmats.py'][0]
   self.assertEqual(rmats[rmats.index('--libType')+1],'fr-firststrand')
   case_bams=(runner.out/'07_rmats/case_bams.txt').read_text()
   self.assertIn('treated0',case_bams); self.assertNotIn('control0',case_bams)
   result=(runner.out/'07_rmats/significant_SE.MATS.JC.txt').read_text()
   self.assertIn('1\t0.001\t-0.2',result); self.assertNotIn('2\t0.9',result)
if __name__=='__main__': unittest.main()

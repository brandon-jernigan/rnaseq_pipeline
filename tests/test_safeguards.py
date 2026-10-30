"""Software-only regression tests. No downloads, installations, or real analyses.
Run: python3 -m unittest discover -s tests -v
"""
import gzip
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ['fastqc', 'trimmomatic', 'hisat2', 'hisat2-build', 'samtools',
         'featureCounts', 'java', 'hisat2_extract_snps_haplotypes_UCSC.py',
         'hisat2_extract_splice_sites.py', 'hisat2_extract_exons.py', 'curl']
# These executables only log arguments and create tiny placeholder files.
STUB = r'''#!/usr/bin/env python3
import gzip, json, os, pathlib, subprocess, sys
p = pathlib.Path
name = p(sys.argv[0]).name
a = sys.argv[1:]
versions = {'fastqc':'FastQC v0.11.9', 'trimmomatic':'0.39',
 'hisat2':'hisat2-align-s version 2.2.1', 'hisat2-build':'hisat2-build-s version 2.2.1',
 'samtools':'samtools 1.10\nUsing htslib 1.10',
 'featureCounts':'featureCounts v' + os.environ.get('TEST_SUBREAD', '2.0.1'),
 'java':'openjdk version "11.0.22" 2024-01-16'}
if a and a[0] in ['--version','-version','-v']:
    print(versions[name], file=sys.stderr if name in ['featureCounts','java'] else sys.stdout)
    sys.exit(0)
fd=os.open(os.environ['TEST_CALLS'],os.O_WRONLY|os.O_CREAT|os.O_APPEND,0o600)
os.write(fd,(json.dumps([name]+a)+'\n').encode()); os.close(fd)
def write(path, text='synthetic\n'):
    p(path).parent.mkdir(parents=True,exist_ok=True); p(path).write_text(text)
if name == 'curl':
    out=a[a.index('-o')+1]; url=a[-1]; basename=url.rsplit('/',1)[-1]
    downloads=p(os.environ['TEST_DOWNLOADS'])
    if basename=='CHECKSUMS':
        choices=list(downloads.glob('*.fa.gz' if '/fasta/' in url else '*.gtf.gz'))
        text=''
        for f in choices:
            fields=subprocess.check_output(['sum','-r',str(f)],text=True).split()
            text+=' '.join(fields[:2]+[f.name])+'\n'
        write(out,text)
    else:
        data=(downloads/basename).read_bytes()
        if os.environ.get('TEST_BAD_DOWNLOAD')=='1': data=data[:-8]
        p(out).write_bytes(data)
        if os.environ.get('TEST_FAIL_DOWNLOAD')=='1': sys.exit(22)
elif name == 'fastqc':
    for f in a[a.index('-o')+2:]:
        stem=p(f).name.removesuffix('.gz').removesuffix('.fastq').removesuffix('.fq')
        for ext in ['html','zip']: write(p(a[a.index('-o')+1])/(stem+'_fastqc.'+ext))
elif name == 'trimmomatic':
    for f in a[6:10]: write(f)
elif name == 'hisat2':
    write(a[a.index('--summary-file')+1]); print('synthetic SAM')
elif name == 'samtools':
    if a[0]=='sort': sys.stdin.read(); write(a[a.index('-o')+1])
    elif a[0]=='index': write(a[1]+'.bai')
elif name == 'featureCounts':
    out=a[a.index('-o')+1]; write(out); write(out+'.summary')
elif name == 'hisat2-build':
    for i in range(1,9):
        write(a[-1]+'.'+str(i)+'.ht2')
        if os.environ.get('TEST_FAIL_BUILD')=='1': sys.exit(9)
elif name == 'hisat2_extract_snps_haplotypes_UCSC.py':
    write(a[-1]+'.snp'); write(a[-1]+'.haplotype')
elif name in ['hisat2_extract_splice_sites.py','hisat2_extract_exons.py']:
    print('synthetic annotation')
else: raise SystemExit('Unexpected tool: '+name)
'''


class Safeguards(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='rnaseq-tests-')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)/'project with spaces'
        shutil.copytree(ROOT/'scripts', self.root/'scripts')
        self.bin = self.root/'mock-bin'; self.bin.mkdir()
        for tool in TOOLS:
            f=self.bin/tool; f.write_text(STUB); f.chmod(0o755)
        envdir=self.root/'.micromamba/envs/rnaseq'
        adapter=envdir/'share/trimmomatic/adapters/TruSeq3-PE.fa'
        adapter.parent.mkdir(parents=True); adapter.write_text('>adapter\nACGT\n')
        mamba=self.root/'bin/micromamba'; mamba.parent.mkdir()
        hook='micromamba() { export CONDA_PREFIX='+shlex.quote(str(envdir))+'; export PATH='+shlex.quote(str(self.bin))+':"$PATH"; }'
        mamba.write_text('#!/usr/bin/env bash\nprintf "%s\\n" '+shlex.quote(hook)+'\n'); mamba.chmod(0o755)
        self.fastq=self.root/'fastq'; self.fastq.mkdir()
        self.ref=self.root/'ref/ensembl_104'
        self.out=self.root/'results'
        self.calls=self.root/'calls.jsonl'
        self.downloads=self.root/'downloads'; self.downloads.mkdir()
        for name, data in [
            ('Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz',b'>1\nACGTACGT\n'),
            ('Homo_sapiens.GRCh38.104.gtf.gz',b'1\tsource\texon\t1\t8\t.\t-\t.\tgene_id "g";\n'),
            ('snp151Common.txt.gz',b'0\tchr1\t1\t2\trs1\t0\t+\tC\tC\tC/T\tgenomic\tsingle\n')]:
            (self.downloads/name).write_bytes(gzip.compress(data))
        self.env={k:v for k,v in os.environ.items() if k in ['PATH','HOME','LANG','LC_ALL','TMPDIR']}
        self.env.update(THREADS='2', TEST_CALLS=str(self.calls), TEST_DOWNLOADS=str(self.downloads), RNASEQ_CLEAN_BASH='1')
        self.pair('sample_R1.fastq.gz','sample_R2.fastq.gz')

    def pair(self, r1, r2):
        for name in [r1,r2]: (self.fastq/name).write_bytes(gzip.compress(b'@read\nACGT\n+\nIIII\n'))

    def run_script(self, script, *args, ok=True, **env):
        r=subprocess.run(['bash',str(self.root/'scripts'/script),*args],cwd=self.root,
                         env={**self.env,**env},capture_output=True,text=True)
        if ok: self.assertEqual(r.returncode,0,r.stdout+r.stderr)
        else: self.assertNotEqual(r.returncode,0,r.stdout+r.stderr)
        return r

    def commands(self):
        return [json.loads(l) for l in self.calls.read_text().splitlines()] if self.calls.exists() else []

    def prepare(self, mode='snp_tran'):
        self.run_script('01_download_reference.sh')
        self.run_script('02_build_hisat2_index.sh', BUILD_INDEXES=mode)

    def test_bash_syntax(self):
        for p in ROOT.glob('scripts/**/*.sh'):
            r=subprocess.run(['bash','-n',str(p)],capture_output=True,text=True)
            self.assertEqual(r.returncode,0,str(p)+r.stderr)

    def test_duplicate_names_fail_before_tools(self):
        self.pair('sample_1.fastq.gz','sample_2.fastq.gz')
        r=self.run_script('run_all.sh',ok=False)
        self.assertIn('Duplicate sample name',r.stderr)
        self.assertEqual(self.commands(),[])

    def test_missing_mate_and_empty_read(self):
        (self.fastq/'sample_R2.fastq.gz').unlink()
        self.assertIn('Missing or empty R2',self.run_script('run_all.sh',ok=False).stderr)
        self.pair('sample_R1.fastq.gz','sample_R2.fastq.gz')
        (self.fastq/'sample_R1.fastq.gz').write_bytes(b'')
        self.assertIn('R1 FASTQ is empty',self.run_script('run_all.sh',ok=False).stderr)
        self.assertEqual(self.commands(),[])

    def test_wrong_subread_rejected_before_tools(self):
        r=self.run_script('run_all.sh',ok=False,TEST_SUBREAD='2.0.2')
        self.assertIn('featureCounts must be 2.0.1',r.stderr)
        self.assertEqual(self.commands(),[])

    def test_missing_environment_never_installs(self):
        (self.root/'bin/micromamba').unlink()
        r=self.run_script('run_all.sh',ok=False)
        self.assertIn('00_install_tools.sh first',r.stderr)
        self.assertEqual(self.commands(),[])

    def test_all_intended_output_types_protected(self):
        paths=['align/sample.sorted.bam','align/sample.sorted.bam.bai',
               'trimmed/sample_R1.trimmed.fastq.gz','trimmed/sample_R2.trimmed.fastq.gz',
               'trimmed/sample_R1.unpaired.fastq.gz','trimmed/sample_R2.unpaired.fastq.gz',
               'qc_raw/sample_R1_fastqc.html','qc_raw/sample_R2_fastqc.zip',
               'qc_trimmed/sample_R1.trimmed_fastqc.html','qc_trimmed/sample_R2.trimmed_fastqc.zip',
               'logs/hisat2_sample.summary.txt','counts/gene_counts.featureCounts.txt',
               'counts/gene_counts.featureCounts.txt.summary']
        for path in paths:
            with self.subTest(path=path):
                f=self.out/path; f.parent.mkdir(parents=True,exist_ok=True); f.write_text('keep me')
                r=self.run_script('run_all.sh',ok=False)
                self.assertIn('Refusing to overwrite',r.stderr)
                self.assertEqual(f.read_text(),'keep me'); f.unlink()
        self.assertEqual(self.commands(),[])

    def test_dangling_output_symlink_protected(self):
        f=self.out/'align/sample.sorted.bam'; f.parent.mkdir(parents=True); f.symlink_to('missing')
        self.assertIn('Refusing to overwrite',self.run_script('run_all.sh',ok=False).stderr)

    def test_default_plain_and_both_build_selection(self):
        for mode, build, expected in [('snp_tran',None,1),('plain',None,1),('snp_tran','both',2)]:
            with self.subTest(mode=mode,build=build):
                if self.ref.exists(): shutil.rmtree(self.ref)
                if self.out.exists(): shutil.rmtree(self.out)
                self.calls.unlink(missing_ok=True)
                extra={'INDEX_TYPE':mode}
                if build: extra['BUILD_INDEXES']=build
                self.run_script('run_all.sh',**extra)
                commands=self.commands(); builds=[c for c in commands if c[0]=='hisat2-build']
                self.assertEqual(len(builds),expected)
                if not build: self.assertEqual('--snp' in builds[0],mode=='snp_tran')
                align=next(c for c in commands if c[0]=='hisat2')
                self.assertTrue(align[align.index('-x')+1].endswith('_'+mode))

    def test_incompatible_selection_fails_early(self):
        self.assertIn('BUILD_INDEXES must match',self.run_script('run_all.sh',ok=False,INDEX_TYPE='plain',BUILD_INDEXES='snp_tran').stderr)
        self.assertEqual(self.commands(),[])

    def test_truncated_and_failed_download_never_published(self):
        for key in ['TEST_BAD_DOWNLOAD','TEST_FAIL_DOWNLOAD']:
            with self.subTest(key=key):
                if self.ref.exists(): shutil.rmtree(self.ref)
                self.run_script('01_download_reference.sh',ok=False,**{key:'1'})
                self.assertFalse((self.ref/'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz').exists())
        self.run_script('01_download_reference.sh')
        self.assertTrue((self.ref/'genome.fa.verified.sha256').exists())

    def test_partial_legacy_extraction_rejected_not_adopted(self):
        self.ref.mkdir(parents=True)
        f=self.ref/'Homo_sapiens.GRCh38.dna.primary_assembly.fa'; f.write_text('>1\nAC')
        r=self.run_script('01_download_reference.sh',ok=False)
        self.assertIn('differs from freshly verified',r.stdout+r.stderr)
        self.assertEqual(f.read_text(),'>1\nAC')
        self.assertFalse(Path(str(f)+'.verified.sha256').exists())

    def test_valid_legacy_reference_adoption(self):
        self.ref.mkdir(parents=True)
        for f in self.downloads.glob('*.gz'): shutil.copyfile(f,self.ref/f.name)
        self.run_script('01_download_reference.sh')
        self.assertTrue((self.ref/'genome.fa.verified.sha256').exists())

    def test_reference_corruption_rejected_on_reuse(self):
        self.run_script('01_download_reference.sh')
        f=self.ref/'genome.fa'; f.write_text('changed\n')
        self.assertIn('checksum mismatch',self.run_script('01_download_reference.sh',ok=False).stdout)
        self.assertEqual(f.read_text(),'changed\n')

    def test_failed_index_is_not_published(self):
        self.run_script('01_download_reference.sh')
        self.run_script('02_build_hisat2_index.sh',ok=False,TEST_FAIL_BUILD='1')
        self.assertFalse(list(self.ref.glob('hisat2_index/*.ht2')))
        self.assertFalse(list(self.ref.glob('hisat2_index/.build.*')))
        self.run_script('02_build_hisat2_index.sh')

    def test_index_without_provenance_rejected(self):
        self.run_script('01_download_reference.sh')
        folder=self.ref/'hisat2_index'; folder.mkdir()
        f=folder/'grch38_ensembl104_snp_tran.1.ht2'; f.write_text('keep')
        self.assertIn('no verified provenance',self.run_script('02_build_hisat2_index.sh',ok=False).stdout)
        self.assertEqual(f.read_text(),'keep')

    def test_index_reused_only_while_unchanged(self):
        self.prepare(); self.calls.unlink()
        self.run_script('02_build_hisat2_index.sh')
        self.assertEqual(self.commands(),[])
        f=next(self.ref.glob('hisat2_index/*.ht2')); f.write_text('corrupt')
        self.assertIn('checksum mismatch',self.run_script('03_run_pipeline.sh',ok=False).stderr)
        self.assertEqual(self.commands(),[])

    def test_scientific_command_contract(self):
        self.prepare(); self.calls.unlink()
        self.run_script('03_run_pipeline.sh')
        c=self.commands()
        raw=[str(self.fastq/'sample_R1.fastq.gz'),str(self.fastq/'sample_R2.fastq.gz')]
        trims=[str(self.out/'trimmed'/f'sample_{mate}.{kind}.fastq.gz')
               for mate,kind in [('R1','trimmed'),('R1','unpaired'),('R2','trimmed'),('R2','unpaired')]]
        self.assertIn(['fastqc','-t','2','-o',str(self.out/'qc_raw'),*raw],c)
        self.assertIn(['trimmomatic','PE','-threads','2','-phred33',*raw,*trims,
                       'ILLUMINACLIP:'+str(self.root/'.micromamba/envs/rnaseq/share/trimmomatic/adapters/TruSeq3-PE.fa')+':2:30:10',
                       'TRAILING:15','MINLEN:20'],c)
        self.assertIn(['fastqc','-t','2','-o',str(self.out/'qc_trimmed'),trims[0],trims[2]],c)
        self.assertIn(['hisat2','-p','2','-x',str(self.ref/'hisat2_index/grch38_ensembl104_snp_tran'),
                       '-1',trims[0],'-2',trims[2],'--qc-filter','--rna-strandness','RF',
                       '--avoid-pseudogene','--dta-cufflinks','--summary-file',str(self.out/'logs/hisat2_sample.summary.txt')],c)
        bam=str(self.out/'align/sample.sorted.bam')
        self.assertIn(['samtools','sort','-@','2','-o',bam,'-'],c)
        self.assertIn(['samtools','index',bam],c)
        self.assertIn(['featureCounts','-T','2','-a',str(self.ref/'Homo_sapiens.GRCh38.104.gtf'),
                       '-o',str(self.out/'counts/gene_counts.featureCounts.txt'),'-C','-p','-s','2','-t','exon',bam],c)
        self.assertEqual(len(c),7)
        before=self.commands()
        self.assertIn('Refusing to overwrite',self.run_script('03_run_pipeline.sh',ok=False).stderr)
        self.assertEqual(self.commands(),before)


if __name__ == '__main__':
    unittest.main()

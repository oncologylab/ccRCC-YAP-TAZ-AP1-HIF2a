from pathlib import Path
import argparse
import json
import re
import shutil
import subprocess


def command(argv,stdout=None):
    return {'argv':[str(v) for v in argv], 'stdout':str(stdout) if stdout else None}


def build_plan(a):
    out=a.output.resolve();p=[];inputs=[]
    if a.stage=='align':
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',a.sample):raise ValueError('Invalid sample ID')
        if a.r1.resolve()==a.r2.resolve():raise ValueError('R1 and R2 must differ')
        inputs=[a.r1,a.r2,a.blacklist,a.picard]
        trimmed=out/'trimmed';sam=out/'aligned.sam';sorted_bam=out/'sorted.bam';rg=out/'readgroups.bam'
        marked=out/'marked.bam';filtered=out/'quality_nodup.bam';final=out/'filtered.bam'
        p=[command(['trim_galore','--paired','--gzip','--fastqc','--cores',a.threads,'--basename',a.sample,'-o',trimmed,a.r1.resolve(),a.r2.resolve()]),
           command(['bowtie2','-p',a.threads,'-x',a.index.resolve(),'-1',trimmed/(a.sample+'_val_1.fq.gz'),'-2',trimmed/(a.sample+'_val_2.fq.gz'),
                    '--very-sensitive-local','--no-mixed','--dovetail','--phred33','-X','1000','-S',sam]),
           command(['samtools','sort','-@',a.threads,'-o',sorted_bam,sam]),
           command(['samtools','addreplacerg','-r','ID:'+a.sample,'-r','SM:'+a.sample,'-r','PL:ILLUMINA','-o',rg,sorted_bam]),
           command(['java','-jar',a.picard.resolve(),'MarkDuplicates','I='+str(rg),'O='+str(marked),'M='+str(out/'duplicate_metrics.txt')]),
           command(['samtools','view','-b','-q',a.mapq,'-F','1028','-o',filtered,marked]),
           command(['bedtools','intersect','-v','-ubam','-a',filtered,'-b',a.blacklist.resolve()],final),
           command(['samtools','index',final])]
    elif a.stage=='peaks':
        inputs=[a.bam,a.chrom_sizes];tag=out/'TagDirectory';peaks=out/'peaks.txt'
        options=['-style','factor','-L','15','-localSize','150000'] if a.style=='factor' else [
            '-style','histone','-localSize','500000','-F','0','-L','0','-C','0','-size','150','-minDist','2500']
        p=[command(['makeTagDirectory',tag,a.bam.resolve(),'-sspe']),
           command(['findPeaks',tag,*options,'-fdr',a.fdr,'-o',peaks]),
           command(['pos2bed.pl',peaks],out/'peaks.bed'),
           command(['makeUCSCfile',tag,'-bigWig',a.chrom_sizes.resolve(),'-norm',a.norm_target,'-o',out/'signal.bw'])]
    elif a.stage=='profiles':
        inputs=a.signals+a.regions;matrix=out/'matrix.gz'
        p=[command(['computeMatrix','reference-point','--referencePoint','center','-S',*[x.resolve() for x in a.signals],
                    '-R',*[x.resolve() for x in a.regions],'-b',a.flank,'-a',a.flank,'--binSize',a.bin_size,'-p',a.threads,'-o',matrix]),
           command(['plotHeatmap','-m',matrix,'--sortRegions','descend','--sortUsing','mean','-out',out/'heatmap.pdf']),
           command(['plotProfile','-m',matrix,'-out',out/'profile.pdf','--outFileNameData',out/'profile.tsv'])]
    elif a.stage=='footprints':
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',a.condition):raise ValueError('Invalid condition ID')
        inputs=[a.signal,a.genome,a.peaks,a.motifs];score=out/'footprints.bw'
        p=[command(['TOBIAS','FootprintScores','--signal',a.signal.resolve(),'--regions',a.peaks.resolve(),'--output',score,'--cores',a.threads]),
           command(['TOBIAS','BINDetect','--motifs',a.motifs.resolve(),'--signals',score,'--genome',a.genome.resolve(),
                    '--peaks',a.peaks.resolve(),'--outdir',out/'binding','--cond_names',a.condition,'--cores',a.threads])]
    return p,[Path(v).resolve() for v in inputs]


def main():
    parser=argparse.ArgumentParser(description=None)
    parser.add_argument('--execute',action='store_true',help='Otherwise print the plan without creating outputs or running tools')
    sub=parser.add_subparsers(dest='stage',required=True)
    align=sub.add_parser('align',help='Paired-end trim, host alignment, duplicate/MAPQ/blacklist filtering')
    for option in ['r1','r2','index','blacklist','picard']:align.add_argument('--'+option,type=Path,required=True)
    align.add_argument('--sample',required=True);align.add_argument('--mapq',type=int,default=20)
    peaks=sub.add_parser('peaks',help='HOMER tag directory, peaks and normalized signal track')
    peaks.add_argument('--bam',type=Path,required=True);peaks.add_argument('--chrom-sizes',type=Path,required=True)
    peaks.add_argument('--style',choices=['factor','histone'],required=True)
    peaks.add_argument('--fdr',type=float,default=1e-5);peaks.add_argument('--norm-target',type=float,default=1e7)
    profiles=sub.add_parser('profiles',help='deepTools centred signal matrix, heatmap and profile')
    profiles.add_argument('--signals',type=Path,nargs='+',required=True);profiles.add_argument('--regions',type=Path,nargs='+',required=True)
    profiles.add_argument('--flank',type=int,default=1000);profiles.add_argument('--bin-size',type=int,default=10)
    footprints=sub.add_parser('footprints',help='TOBIAS scores/binding from an already bias-corrected signal')
    for option in ['signal','genome','peaks','motifs']:footprints.add_argument('--'+option,type=Path,required=True)
    footprints.add_argument('--condition',required=True)
    for stage in [align,peaks,profiles,footprints]:
        stage.add_argument('--output',type=Path,required=True);stage.add_argument('--threads',type=int,default=8)
    args=parser.parse_args()
    if args.threads<1:parser.error('threads must be positive')
    if args.stage=='align' and not 0<=args.mapq<=255:parser.error('invalid MAPQ')
    if args.stage=='peaks' and not (0<args.fdr<1 and args.norm_target>0):parser.error('invalid FDR/normalization target')
    if args.stage=='profiles' and min(args.flank,args.bin_size)<1:parser.error('flank and bin size must be positive')
    plan,inputs=build_plan(args)
    if not args.execute:
        print(json.dumps({'stage':args.stage,'execute':False,'commands':plan},indent=2));return
    for path in inputs:
        if not path.is_file():raise FileNotFoundError(path)
    if args.stage=='align' and not any(Path(str(args.index.resolve())+suffix).is_file() for suffix in ['.1.bt2','.1.bt2l']):
        raise FileNotFoundError('Bowtie2 index prefix: '+str(args.index))
    for name in {c['argv'][0] for c in plan}:
        if shutil.which(name) is None:raise FileNotFoundError('Required executable: '+name)
    if args.output.exists():raise FileExistsError('Choose a new output directory: '+str(args.output))
    args.output.mkdir(parents=True)
    if args.stage=='align':(args.output/'trimmed').mkdir()
    (args.output/'commands.json').write_text(json.dumps(plan,indent=2)+'\n')
    for i,step in enumerate(plan,1):
        stdout=Path(step['stdout']) if step['stdout'] else args.output/f'{i:02d}.stdout.log'
        with stdout.open('xb') as output,(args.output/f'{i:02d}.stderr.log').open('xb') as error:
            subprocess.run(step['argv'],check=True,stdout=output,stderr=error)


if __name__=='__main__':main()

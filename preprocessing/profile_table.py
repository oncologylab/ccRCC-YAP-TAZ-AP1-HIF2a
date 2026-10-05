from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def convert(path,before,after):
    rows=[];group=None
    for line in Path(path).read_text().splitlines():
        fields=line.split('\t')
        if not line.strip() or line.startswith('#') or fields[0] in {'bin labels','bins'}:continue
        if len(fields)==1:group=fields[0].strip();continue
        if len(fields)<3:raise ValueError('Expected sample, region and numeric bins')
        label=(group+' '+fields[0]).strip() if group else fields[0].strip()
        numbers=np.array([float(x) for x in fields[2:] if x!=''])
        if not len(numbers) or not np.isfinite(numbers).all():raise ValueError('Missing/non-finite bins')
        width=(before+after)/len(numbers)
        positions=-before+(np.arange(len(numbers))+.5)*width
        rows.extend({'panel':fields[1].strip(),'group':label,'x':x,'value':v} for x,v in zip(positions,numbers))
    out=pd.DataFrame(rows)
    if out.empty or out.duplicated(['panel','group','x']).any():raise ValueError('Empty or duplicated profiles')
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=None)
    parser.add_argument('--input',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--upstream',type=float,required=True);parser.add_argument('--downstream',type=float,required=True)
    args=parser.parse_args()
    if min(args.upstream,args.downstream)<0 or args.upstream+args.downstream<=0:parser.error('Invalid flanks')
    if args.output.exists():raise FileExistsError(args.output)
    result=convert(args.input,args.upstream,args.downstream)
    args.output.parent.mkdir(parents=True,exist_ok=True);result.to_csv(args.output,sep='\t',index=False)

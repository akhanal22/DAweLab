from pathlib import Path
import gzip,re,collections,bisect,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent.parent;out=root/'variant_enrichment';files=list((root/'upload').glob('*.gz'))
sources={'SNVs':next(f for f in files if f.name.startswith('private_biallelic_snps')),'Indels (1–50 bp)':next(f for f in files if f.name.startswith('private_biallelic_indels')),'SVs (>50 bp)':next(f for f in files if f.name.startswith('syri.final'))}
lengths={}
with gzip.open(sources['SNVs'],'rt') as h:
 for l in h:
  m=re.match(r'##contig=<ID=(chr(?:[1-9]|10)_RagTag),length=(\d+)',l)
  if m:lengths[m[1]]=int(m[2])
  if l.startswith('#CHROM'):break
assert len(lengths)==10
# Sweep annotation boundaries into nonoverlapping intervals. Priority CDS > UTR > other exon > intron > intergenic.
events={c:collections.defaultdict(lambda:np.zeros(4,dtype=int)) for c in lengths}
with open(root/'upload'/'J702.liftoff.gff3') as h:
 for l in h:
  if l.startswith('#'):continue
  f=l.rstrip().split('\t')
  if len(f)!=9 or f[0] not in lengths:continue
  typ=f[2];j={'gene':0,'exon':1,'CDS':2,'five_prime_UTR':3,'three_prime_UTR':3}.get(typ)
  if j is None:continue
  a=int(f[3])-1;b=int(f[4]);assert 0<=a<b<=lengths[f[0]]
  events[f[0]][a][j]+=1;events[f[0]][b][j]-=1
segments={};bp=collections.Counter()
for c,L in lengths.items():
 ev=events[c];ev[0];ev[L];active=np.zeros(4,dtype=int);starts=[];cats=[];last=0
 for p in sorted(ev):
  if p>last:
   cat='CDS' if active[2] else 'UTR' if active[3] else 'Other exon' if active[1] else 'Intron' if active[0] else 'Intergenic'
   starts.append(last);cats.append(cat);bp[cat]+=p-last
  active+=ev[p];assert np.all(active>=0);last=p
 segments[c]=(starts,cats)
assert sum(bp.values())==sum(lengths.values())
catorder=[c for c in ['CDS','UTR','Intron','Intergenic','Other exon'] if bp[c]]
print('Genome category bp:',dict(bp))
rows=[];audit=collections.Counter();samples=['B73D1','B73D2','B73D4','B73D7c','B73D10']
for cls,path in sources.items():
 with gzip.open(path,'rt') as h:
  for l in h:
   if l.startswith('#CHROM'):names=l.rstrip().split('\t')[9:];continue
   if l.startswith('#'):continue
   f=l.rstrip().split('\t');c=f[0]
   if c not in lengths:continue
   ref,alt=f[3:5]
   if ',' in alt or alt.startswith('<') or '[' in alt or ']' in alt:audit[(cls,'nonsequence_or_multiallelic')]+=1;continue
   if 'NNN' in ref.upper() or 'NNN' in alt.upper():audit[(cls,'NNN_excluded')]+=1;continue
   size=abs(len(alt)-len(ref))
   if cls=='SNVs' and not (len(ref)==len(alt)==1):continue
   if cls.startswith('Indels') and not (1<=size<=50):continue
   if cls.startswith('SVs') and size<=50:continue
   fmt=f[8].split(':');g=fmt.index('GT');carriers=[]
   for name,field in zip(names,f[9:]):
    gt=field.split(':')[g];alleles=re.split(r'[/|]',gt)
    if any(a not in ('.','0','') for a in alleles):carriers.append(name.replace('B73D7c1','B73D7c'))
   if len(carriers)!=1:audit[(cls,'not_exactly_one_ALT_carrier')]+=1;continue
   # Assign one locus per event: SNV base; first deleted/changed base; insertion's left junction base.
   pos=int(f[1])-1;prefix=0
   if cls!='SNVs':
    while prefix<min(len(ref),len(alt)) and ref[prefix]==alt[prefix]:prefix+=1
    pos+=prefix if len(ref)>len(alt) else max(0,prefix-1)
   assert 0<=pos<lengths[c]
   starts,cats=segments[c];cat=cats[bisect.bisect_right(starts,pos)-1]
   rows.append({'sample':carriers[0],'variant_class':cls,'chromosome':c,'VCF_POS':f[1],'id':f[2],'length_bp':size if cls!='SNVs' else 1,'assigned_locus_1based':pos+1,'category':cat})
df=pd.DataFrame(rows);df.to_csv(out/'variant_category_assignments.tsv',sep='\t',index=False)
genome=pd.DataFrame([{'category':c,'bp':bp[c],'genome_percent':100*bp[c]/sum(bp.values())} for c in catorder]);genome.to_csv(out/'genome_category_proportions.tsv',sep='\t',index=False)
stats=[]
for sample in samples+['All lines']:
 sub=df if sample=='All lines' else df[df['sample']==sample]
 for cls in sources:
  x=sub[sub.variant_class==cls];n=len(x)
  for c in catorder:
   k=int((x.category==c).sum());fraction=bp[c]/sum(bp.values());stats.append(dict(sample=sample,variant_class=cls,category=c,count=k,total=n,variant_percent=100*k/n if n else np.nan,genome_percent=100*fraction,enrichment=(k/n)/fraction if n else np.nan))
st=pd.DataFrame(stats);st.to_csv(out/'normalized_variant_distribution.tsv',sep='\t',index=False)
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
colors={'CDS':'#2878b5','UTR':'#D55E00','Intron':'#009e73','Intergenic':'#8797a6','Other exon':'#a66bbe'}
def plot(sample,ax):
 ss=st[st['sample']==sample];xs=np.arange(3);w=.8/len(catorder)
 for j,c in enumerate(catorder):
  vals=[ss[(ss.variant_class==cl)&(ss.category==c)].enrichment.iloc[0] for cl in sources]
  bars=ax.bar(xs-.4+w/2+j*w,vals,width=w*.93,color=colors[c],label=c)
  for bar,v in zip(bars,vals):ax.text(bar.get_x()+bar.get_width()/2,v+.025,f'{v:.2f}',ha='center',va='bottom',fontsize=8)
 ax.axhline(1,color='#333333',ls='--',lw=1);ax.set_xticks(xs,[f'{cl}\nn = {len(df[(df.variant_class==cl)&((df["sample"]==sample) if sample!="All lines" else True)]):,}' for cl in sources]);ax.set_ylabel('Fold enrichment');ax.set_title(sample if sample!='All lines' else 'Private variants normalized by genomic category size',fontweight='bold');ax.set_ylim(0,max(1.2,ss.enrichment.max()*1.22));ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
fig,ax=plt.subplots(figsize=(12,6));plot('All lines',ax);ax.legend(ncol=len(catorder),loc='upper center',bbox_to_anchor=(.5,1.0),frameon=False);fig.text(.5,.015,'Variant proportion ÷ genomic bp proportion • 1× = expected • chromosomes 1–10 • pooled across five lines',ha='center',fontsize=10);fig.tight_layout(rect=[0,.055,1,1]);fig.savefig(out/'normalized_variant_distribution.png',dpi=250);fig.savefig(out/'normalized_variant_distribution.pdf');plt.close(fig)
fig,axes=plt.subplots(3,2,figsize=(14,13));
for ax,s in zip(axes.flat,samples):plot(s,ax)
axes.flat[-1].axis('off');handles,labels=axes.flat[0].get_legend_handles_labels();axes.flat[-1].legend(handles,labels,loc='center',frameon=False);fig.suptitle('Genomic enrichment by mutation-accumulation line',fontsize=18);fig.tight_layout(rect=[0,.02,1,.96]);fig.savefig(out/'normalized_variant_distribution_by_line.png',dpi=220);fig.savefig(out/'normalized_variant_distribution_by_line.pdf');plt.close(fig)
(out/'METHODS.txt').write_text('Denominator: full reference chromosome 1–10 lengths from SNP VCF header, including any assembly gaps; not callable-normalized. Categories are nonoverlapping unions across all annotated isoforms, ignoring strand, priority CDS > annotated UTR > other exon > intron > intergenic. Intron is gene body outside all annotated exons, CDS and UTR. Other exon preserves exonic bases not annotated CDS/UTR. Counts are events assigned once by SNV base, first changed/deleted base, or insertion left junction. This measures locus enrichment, not full SV span overlap or affected bp. SNPs and private biallelic indels use their supplied files; only 1–50 bp indels retained. SyRI uses sequence-defined indels >50 bp. NNN alleles excluded. Exactly one ALT-carrying sample required. No additional PASS/DP/GQ filtering. All-lines pools event occurrences across lines. Fold enrichment = within-class variant fraction / genomic bp fraction. Descriptive; no significance test.\nAudit exclusions: '+str(dict(audit)))
print(df.groupby(['variant_class','sample']).size().to_string());print(genome.to_string(index=False));print('Audit',audit)

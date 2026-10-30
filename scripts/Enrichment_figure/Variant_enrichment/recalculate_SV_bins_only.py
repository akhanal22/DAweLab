from pathlib import Path
import runpy,collections,bisect
import pandas as pd,numpy as np
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parent.parent
v=runpy.run_path(str(root/'variant_enrichment'/'normalize_variants.py'))
out=root/'variant_enrichment_SV_bins_only';out.mkdir(exist_ok=True)
df=v['df'];samples=v['samples'];cats=v['catorder'];lengths=v['lengths'];segments=v['segments'];bp=v['bp'];sources=v['sources'];colors=v['colors']

validation=[]
for typ,pattern,depth,flags in [
 ('INS','insertion_validation*tsv','own_depth_ratio',['own_start_pass','own_end_pass','j702_junction_pass']),
 ('DEL','deletion_validation*tsv','j702_depth_ratio',['own_junction_pass','j702_start_pass','j702_end_pass'])]:
 x=pd.read_csv(next((root/'upload').glob(pattern)),sep='\t');x['sample']=x['sample'].replace({'B73D7c1':'B73D7c'});x['original_verdict']=x.verdict
 x['required_depth_ratio']=.70;x['updated_depth_pass']=np.where(x[depth]>=.70,'PASS','FAIL')
 # Preserve non-depth failures, including junction edge errors.
 edge=x.notes.fillna('').str.contains('junction_too_close_to_contig_edge_for_anchor')
 ok=x[flags].eq('PASS').all(axis=1)&(x[depth]>=.70)&~edge
 x['verdict']=np.where(ok,'PASS','FAIL');validation.append(x)
 x.to_csv(out/f'{typ}_validation_depth070.tsv',sep='\t',index=False)
 print(typ,'old PASS',sum(x.original_verdict=='PASS'),'new PASS',sum(ok),'rescued',sum(ok&(x.original_verdict!='PASS')))
val=pd.concat(validation,ignore_index=True);passed=val[val.verdict=='PASS'][['sample','chrom','pos','id','length']]
df['VCF_POS']=df['VCF_POS'].astype(int)
sv=df[df.variant_class=='SVs (>50 bp)'];matched=sv.merge(passed,left_on=['sample','chromosome','VCF_POS','id','length_bp'],right_on=['sample','chrom','pos','id','length'],validate='one_to_one')
assert len(matched)==len(passed)
df=pd.concat([df[df.variant_class!='SVs (>50 bp)'],matched[df.columns]],ignore_index=True)

bad={};audit=[];genome=[];fractions={}
for s in samples:
 path=root/'upload'/f'{s.replace("B73D7c","B73D7c1")}.bad_500kb_bins.bed';raw=collections.defaultdict(list)
 for l in path.read_text().splitlines():
  if not l or l.startswith(('#','track','browser')):continue
  c,a,b=l.split()[:3];a=int(a);b=int(b)
  if c in lengths:
   a=max(0,a);b=min(lengths[c],b)
   if a<b:raw[c].append((a,b))
 merged={}
 for c,items in raw.items():
  arr=[]
  for a,b in sorted(items):
   if arr and a<=arr[-1][1]:arr[-1]=(arr[-1][0],max(arr[-1][1],b))
   else:arr.append((a,b))
  merged[c]=arr
 bad[s]=merged;remaining=dict(bp);excluded=0
 for c,items in merged.items():
  starts,labels=segments[c]
  for a,b in items:
   excluded+=b-a;i=bisect.bisect_right(starts,a)-1;p=a
   while p<b:
    end=min(b,starts[i+1] if i+1<len(starts) else lengths[c]);remaining[labels[i]]-=end-p;p=end;i+=1
 total=sum(remaining.values());assert total==sum(lengths.values())-excluded
 fractions[s]={c:remaining[c]/total for c in cats}
 for c in cats:genome.append(dict(sample=s,category=c,retained_bp=remaining[c],retained_genome_bp=total,excluded_bp=excluded,retained_genome_percent=100*fractions[s][c]))
def inbad(r):
 intervals=bad[r['sample']].get(r.chromosome,[]);p=int(r.assigned_locus_1based)-1;i=bisect.bisect_right(intervals,(p,float('inf')))-1
 return i>=0 and intervals[i][0]<=p<intervals[i][1]
mask=df.apply(inbad,axis=1) & df.variant_class.eq('SVs (>50 bp)');removed=df[mask];df=df[~mask];removed.to_csv(out/'variants_in_bad_bins_excluded.tsv',sep='\t',index=False);df.to_csv(out/'retained_variant_assignments.tsv',sep='\t',index=False)
stat=[]
for s in samples:
 for cl in sources:
  d=df[(df['sample']==s)&(df.variant_class==cl)];n=len(d)
  audit.append(dict(sample=s,variant_class=cl,retained_variants=n,excluded_variants=len(removed[(removed['sample']==s)&(removed.variant_class==cl)])))
  for c in cats:
   k=int((d.category==c).sum());f=fractions[s][c] if cl=='SVs (>50 bp)' else bp[c]/sum(bp.values());stat.append(dict(sample=s,variant_class=cl,category=c,count=k,total=n,expected_count=n*f,background_percent=100*f,variant_percent=100*k/n if n else np.nan,enrichment=k/(n*f) if n and f else np.nan))
st=pd.DataFrame(stat)
for cl in sources:
 for c in cats:
  d=st[(st.variant_class==cl)&(st.category==c)];k=d['count'].sum();n=d.total.sum();expected=d.expected_count.sum();stat.append(dict(sample='All lines',variant_class=cl,category=c,count=k,total=n,expected_count=expected,background_percent=100*expected/n,variant_percent=100*k/n,enrichment=k/expected))
st=pd.DataFrame(stat);st.to_csv(out/'normalized_variant_distribution.tsv',sep='\t',index=False);pd.DataFrame(genome).to_csv(out/'retained_genome_proportions.tsv',sep='\t',index=False);pd.DataFrame(audit).to_csv(out/'bad_bin_variant_audit.tsv',sep='\t',index=False)
def plot(s,ax):
 d=st[st['sample']==s];xs=np.arange(3);w=.19
 for j,c in enumerate(cats):
  vals=[d[(d.variant_class==cl)&(d.category==c)].enrichment.iloc[0] for cl in sources];bars=ax.bar(xs-.285+j*w,vals,width=w*.93,color=colors[c],label=c)
  for bar,val in zip(bars,vals):ax.text(bar.get_x()+bar.get_width()/2,val+.025,f'{val:.2f}',ha='center',va='bottom',fontsize=8)
 ax.axhline(1,color='#333',ls='--',lw=1);ax.set_xticks(xs,[f'{cl}\nn = {len(df[(df.variant_class==cl)&((df["sample"]==s) if s!="All lines" else True)]):,}' for cl in sources]);ax.set_ylabel('Fold enrichment');ax.set_title(s if s!='All lines' else 'Private variant enrichment: bad bins excluded for SVs only',fontweight='bold');ax.set_ylim(0,max(1.2,d.enrichment.max()*1.22));ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
fig,ax=plt.subplots(figsize=(12,6));plot('All lines',ax);ax.legend(ncol=4,loc='upper center',bbox_to_anchor=(.5,1),frameon=False);fig.text(.5,.015,'DeepVariant: full-genome background • validated SyRI SVs: bad bins excluded • 1× = expected',ha='center',fontsize=10);fig.tight_layout(rect=[0,.055,1,1]);fig.savefig(out/'normalized_variant_distribution_validated.png',dpi=250);fig.savefig(out/'normalized_variant_distribution_validated.pdf');plt.close(fig)
fig,axes=plt.subplots(3,2,figsize=(14,13))
for ax,s in zip(axes.flat,samples):plot(s,ax)
axes.flat[-1].axis('off');h,l=axes.flat[0].get_legend_handles_labels();axes.flat[-1].legend(h,l,loc='center',frameon=False);fig.suptitle('Variant enrichment: bad bins excluded; validated SVs only',fontsize=18);fig.tight_layout(rect=[0,.02,1,.96]);fig.savefig(out/'normalized_variant_distribution_by_line_validated.png',dpi=220);fig.savefig(out/'normalized_variant_distribution_by_line_validated.pdf');plt.close(fig)
(out/'METHODS.txt').write_text('UPDATE: bad bins applied ONLY to SVs; DeepVariant SNVs and small indels use all supplied events and full chr1–10 backgrounds. The below general bad-bin notes apply to SVs only. '+'SV validation: depth ratio >=0.70 for INS own interval and DEL J702 interval (sequence-bearing genome); all existing junction pass requirements and edge exclusions preserved. SVs matched by sample, chromosome, position, ID and length. SNVs and small indels unchanged. '+'Bad bins are BED 0-based half-open, clipped to chr1–10 and merged within each line. Each line uses its own exclusions, not the union across lines. Category bp are recalculated outside those bins. Variant loci inside bad bins are removed. Remaining reference is not a full callable mask. Categories and point assignment match prior analysis. Pooled enrichment = sum(observed category counts) / sum(N_line_class × retained category bp_line / retained genome bp_line). This weights each background by its line/class event count. SVs assigned by insertion left junction or first deleted base; full SV overlap with bad bins is not tested.\n')
print(pd.DataFrame(audit).to_string(index=False));print(pd.DataFrame(genome).groupby('sample').first()[['excluded_bp','retained_genome_bp']].to_string());print(st[st['sample']=='All lines'][['variant_class','category','enrichment']].to_string(index=False))

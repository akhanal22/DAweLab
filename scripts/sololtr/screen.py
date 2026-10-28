import csv,collections,json
from pathlib import Path
R=Path('/workspace/scratch/48118ffd80f5');O=R/'solo_LTR_screen'
parents={};ltrs=collections.defaultdict(list)
for line in open(R/'upload/J702_flip.fa.mod.EDTA.intact(3).gff3'):
 if line.startswith('#'):continue
 a=line.rstrip().split('\t')
 if len(a)!=9:continue
 d=dict(x.split('=',1) for x in a[8].split(';') if '=' in x)
 if a[2]=='repeat_region':parents[d['ID']]={'chrom':a[0],**d}
 if a[2]=='long_terminal_repeat':ltrs[d['Parent']].append((int(a[3])-1,int(a[4])))
bychr=collections.defaultdict(list)
for pid,ls in ltrs.items():
 assert len(ls)==2,(pid,ls)
 ls.sort();bychr[parents[pid]['chrom']].append((pid,*ls))
rows=list(csv.DictReader(open(R/'upload/deletion_validation_withcigar_for500kbfilteredvcfagainstj702 (1)(2).tsv'),delimiter='\t'))
valid=[r for r in rows if int(r['length'])>50 and float(r['j702_depth_ratio'])>=.7 and all(r[k]=='PASS' for k in ['own_junction_pass','j702_start_pass','j702_end_pass'])]
result=[];overlap=set();counts=collections.Counter()
for r in valid:
 # normalized VCF simple deletion removes POS+1 through END inclusive;
 # reference cut sites in BED coordinates are POS and END.
 s=int(r['pos']);e=int(r['end']);assert e-s==int(r['length'])
 for pid,(a,b),(c,d) in bychr[r['chrom']]:
  if not (s<d and e>a):continue
  overlap.add((r['sample'],r['chrom'],r['pos']))
  dl=max(a-s,0,s-b);dr=max(c-e,0,e-d)
  if dl>1000 or dr>1000:continue
  strict=a<=s<=b and c<=e<=d
  tier='paired_LTR_breakpoints' if strict else ('within_100bp' if max(dl,dr)<=100 else 'within_1000bp')
  counts[tier]+=1
  result.append({'sample':r['sample'],'chrom':r['chrom'],'vcf_pos_1based':s,'deleted_start_1based':s+1,'deleted_end_1based':e,'deletion_length':r['length'],'id':r['id'],'element_id':pid,'classification':parents[pid].get('classification'),'left_LTR_start_1based':a+1,'left_LTR_end_1based':b,'right_LTR_start_1based':c+1,'right_LTR_end_1based':d,'breakpoint_distance_left_bp':dl,'breakpoint_distance_right_bp':dr,'tier':tier,'LTR_identity':parents[pid].get('ltr_identity'),'remaining_element_bp':d-a-(e-s),'left_LTR_length':b-a,'right_LTR_length':d-c,'depth_ratio':r['j702_depth_ratio']})
fields=list(result[0]) if result else ['sample','chrom','tier']
for name,rs in [('paired_LTR_candidates.tsv',[r for r in result if r['tier']=='paired_LTR_breakpoints']),('all_candidates_within1000bp.tsv',result)]:
 with open(O/name,'w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rs)
summary={'intact_LTR_elements':len(ltrs),'validated_deletions':len(valid),'deletions_overlapping_intact_LTR_elements':len(overlap),'candidate_pairs_by_tier':dict(counts),'strict_unique_deletions':len({(r['sample'],r['chrom'],r['vcf_pos_1based']) for r in result if r['tier']=='paired_LTR_breakpoints'})}
(O/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
for r in result:print(r)

import isotools
print(f'This is isotools version {isotools.__version__}')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logging, os, time, pickle

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
logger = logging.getLogger('isotools')

# ---- Paths ----
sample = 'alzheimer2019_isoseq'
bam = f'alzheimer/aligned/{sample}_aligned.sorted.bam'
genome = 'reference/GRCh38.p13.genome.fa'
anno = 'reference/gencode.v36.chr_patch_hapl_scaff.annotation_sorted'
pickle_file = f'alzheimer/pickle/{sample}_isotools.pkl'
ref_pickle = anno + '.isotools.pkl'
gff_file = anno + '.gff3.gz'

# ---- Load IsoTools object ----
try:
    print("Trying to load processed pickle...")
    isoseq = isotools.Transcriptome.load(pickle_file)
    print("Loaded existing IsoTools object.")
except FileNotFoundError:
    try:
        print("Trying to load reference pickle...")
        isoseq = isotools.Transcriptome.from_reference(ref_pickle)
    except FileNotFoundError:
        print("Building reference from GFF3...")
        isoseq = isotools.Transcriptome.from_reference(gff_file)
        isoseq.save_reference(ref_pickle)

    print("Adding sample from BAM (30–40 min)...")
    isoseq.add_sample_from_bam(bam, sample_name='alzheimer_1', group='alzheimer', min_align_fraction=0)

    print("Computing QC metrics (15–20 min)...")
    isoseq.add_qc_metrics(genome)

    isoseq.make_index()
    isoseq.save(pickle_file)
    print(f"Saved processed object to {pickle_file}")

# ---- Plot saturation & rarefaction ----
from isotools.plots import plot_saturation, plot_rarefaction, plot_bar, plot_distr

plt.rcParams["figure.figsize"] = (14,7)
fig1, axs = plt.subplots(1,2)

plot_saturation(isoseq, cov_th=2, x_range=(1e4,5e6,1e4), ax=axs[0])
rarefaction, total = isoseq.rarefaction(min_coverage=2, tr_filter={'query':'FSM'})
plot_rarefaction(rarefaction, total=total, ax=axs[1])

fig1.tight_layout()
plt.show(block=False)
plt.pause(3)   # show first plot briefly, then continue

# ---- Load pickle directly for timing test ----
t0 = time.time()
with open(pickle_file, "rb") as f:
    obj = pickle.load(f)
print(f"Loaded in {time.time()-t0:.2f} seconds")
print(type(obj))

# ---- Add filters ----
isoseq.add_filter("HIGH_SUPPORT", 'transcript_support_level=="1"', context='reference')
isoseq.add_filter("PROTEIN_CODING", 'transcript_type=="protein_coding"', context='reference')

# Print filters
for context in isoseq.filter:
    print(f'\n{context}\n{"="*len(context)}')
    for tag, expression in isoseq.filter[context].items():
        print(f'- {tag}:\t{expression}')

# ---- Compute stats ----
tr_stats = [
    isoseq.transcript_length_hist(groups=isoseq.groups(), add_reference=True,
                                  min_coverage=2, tr_filter=dict(query='not NOVEL_GENE'),
                                  ref_filter=dict(query='HIGH_SUPPORT')),
    isoseq.downstream_a_hist(groups=isoseq.groups(), transcript_filter=dict(query='FSM'),
                             ref_filter=dict(query='not REF_UNSPLICED')),
    isoseq.downstream_a_hist(groups=isoseq.groups(), transcript_filter=dict(query='NOVEL_GENE and UNSPLICED')),
    isoseq.direct_repeat_hist(groups=isoseq.groups(), bins=np.linspace(-.5,10.5,12))
]
tr_stats.append((pd.concat([tr_stats[1][0].add_suffix(' known'),
                            tr_stats[2][0].add_suffix(' novel')]), tr_stats[2][1]))

# Filter stats
f_stats = isoseq.filter_stats(weight_by_coverage=True, min_coverage=1,
                              tags=['INTERNAL_PRIMING', 'RTTS', 'FRAGMENT'])
f_stats[0].index = f_stats[0].index.str.replace('_', '\n')

# ---- QC plots ----
plt.rcParams["figure.figsize"] = (15,15)
plt.rcParams.update({'font.size': 14})

fig2, axs = plt.subplots(2,2)

# A) transcript length
plot_distr(tr_stats[0][0], smooth=3, ax=axs[0,0], **tr_stats[0][1])
# B) internal priming
plot_distr(tr_stats[4][0], smooth=3, ax=axs[0,1], density=True, fill=True, **tr_stats[4][1])
# C) RTTS
plot_distr(tr_stats[3][0], ax=axs[1,0], density=True, **tr_stats[3][1])
# D) frequency of artifacts
plot_bar(f_stats[0], ax=axs[1,1], drop_categories=['PASS'], colors=['blue'], **f_stats[1])

fig2.tight_layout()
plt.show()   # stays open until you close

ref=[[[12,20],[30,40], [50,60],[70,81]],
     [[11,20],[35,40],         [75,79]],
     [[10,20],[30,40], [50,60],[70,80]]]
novel={'FSM':         [[10,20],[30,40], [50,60],[70,80]],
       "5' fragment": [[33,40], [50,60],[70,80]],
       "3' fragment": [[10,20],[30,40], [50,55]],
       "mono exon"  : [[22,35]],
       "exon skipping"     :  [[10,20], [50,60],[70,80]],
       "intron retention"  :  [[10,40], [50,60],[70,80]],
       "novel combination" :  [[10,20],[35,40], [50,60],[70,80]],
       "novel junction"  :   [[10,20],[30,40], [50,60], [75,80]],
       "novel exonic TSS"  :  [[26,40], [50,60],[70,80]],
       "novel exonic PAS"  :  [[10,20],[30,40], [50,66]],
       "novel 5' splice site":[[10,24],[30,40], [50,60],[70,80]],
       "novel 3' splice site":[[10,20],[26,40], [50,60],[70,80]],
       "novel exon"  :        [[10,20],[30,40],[43,47], [50,60],[70,80]],
       "novel intronic TSS" : [[43,47],[50,60],[70,80]],
       "novel intronic PAS" : [[10,20],[30,40], [82,90]]}
ref={'transcripts':[{'exons':e, 'transcript_name':f'reference {i+1}'} for i,e in enumerate(ref)]}
transcripts=[{'exons':e, 'transcript_name':n} for n,e in novel.items()]
example=isotools.Gene(10,80,{'strand':'+','ID':'example','reference':ref, 'transcripts':transcripts},None)
f,axs=plt.subplots(2,figsize=(10,7), gridspec_kw={'height_ratios': [1, 4]})
cat=['FSM','ISM','NIC','NNC','novel gene']
sg=example.ref_segment_graph
for novel in example.transcripts:
    alt_splice=sg.get_alternative_splicing(novel['exons'])
    print(f"{novel['transcript_name']}: {alt_splice[1]}")
    novel['transcript_name']=f"{','.join(alt_splice[1])} ({cat[alt_splice[0]]}) "

example.gene_track(ax=axs[0], x_range=[10,90], title='')
example.gene_track(reference=False,ax=axs[1], x_range=[10,90], title='', color='green')
for ax in axs:
    ax.get_xaxis().set_visible(False)
f.tight_layout()
plt.show()

cnr={}
for gene, transcript_id, transcript in isoseq.iter_transcripts():
    for anno in transcript['annotation'][1]:
        cnr[anno]=min(cnr.get(anno,5),transcript['annotation'][0])
del cnr['FSM']
altsplice=[ isoseq.altsplice_stats(groups=isoseq.groups(), weight_by_coverage=True, min_coverage=1, tr_filter=dict( query="not( RTTS or FRAGMENT or INTERNAL_PRIMING)")),
            isoseq.altsplice_stats(groups=isoseq.groups(), weight_by_coverage=True, min_coverage=2, tr_filter=dict( query="not( RTTS or FRAGMENT or INTERNAL_PRIMING)")),
            isoseq.altsplice_stats(groups=isoseq.groups(), weight_by_coverage=False, min_coverage=20, tr_filter=dict( query="not( RTTS or FRAGMENT or INTERNAL_PRIMING)"))]
for i in range(3):
    altsplice[i][0].index=altsplice[i][0].index+[f'\n({cat[cnr[subcat]]})' if subcat in cnr else '' for subcat in altsplice[i][0].index]
    altsplice[i][0].index=altsplice[i][0].index.str.replace('splice ','\nsplice ')

from isotools.plots import plot_bar, plot_distr

plt.rcParams["figure.figsize"] = (20,10)
_=plot_bar(altsplice[0][0], bar_width=.9, ylabel='fraction of reads [%]',
           colors=['blue'], legend=False, rot=90, drop_categories=['FSM'])

plt.tight_layout()
plt.show()

#access gene of interest by Name
gene=isoseq['MAPT']
print(gene)
#this reveals that the MAPT has 267 different variants, according to isoseq.
total_cov=gene.coverage.sum()
#However, most are supported by few reads only
print(f'{sum([cov>total_cov *.01 for cov in gene.coverage][0])} transcripts contribute at least 1% to that gene')
#lets look at the primary transcript
max_i=np.argmax(gene.coverage)
print(f'The primary transcript is number {max_i} and contributes {gene.coverage[0,max_i]/total_cov:.2%}  ({gene.coverage[0,max_i]}/{total_cov})')
#all the information for this transcript are stored in this dict:
primary=gene.transcripts[max_i]
print(f'\nThese are the infos for this transcript:')
for k,v in primary.items():
    print(f'{k}: {str(v)[:100]}{"..." if len(str(v))>100 else ""}')
# this reveals that it is a mono exon in the 3'UTR of that gene.
second_i=np.argsort(gene.coverage[0])[-2]
second=gene.transcripts[second_i]
print(f'\nThese are the infos for the second transcript:')
for k,v in second.items():
    print(f'{k}: {str(v)[:100]}{"..." if len(str(v))>100 else ""}')

# this reveals that it is a FSM with reference transcript nr 7.

print(f'\nThe corresponding reference transcript: ')
for k,v in gene.ref_transcripts[second["annotation"][1]["FSM"][0]].items():
    print(f'{k}: {str(v)[:100]}{"..." if len(str(v))>100 else ""}')



#we can iterate over transcripts and filter with a query:
i=0
for gene, transcript_number, transcript in isoseq.iter_transcripts(query='INTERNAL_PRIMING'):
    print(f'transcript nr {transcript_number} of {gene} with a coverage of {gene.coverage[0,transcript_number]}')
    i+=1
    if i>10:
        break

#we can also iterate over transcripts and filter by region (chr[:start-end]), novelty category, or request a minimum coverage:
for gene, transcript_number, transcript in isoseq.iter_transcripts(region='chr1', query='NOVEL_EXON', min_coverage=100):
    print(f'Gene {gene.name} has a highly covered transcript with novel exon:')
    print(f'Transcript nr {transcript_number}  with a coverage of {gene.coverage[0,transcript_number]}/{gene.coverage.sum()}')
    print(f'The novel exon is at position {transcript["annotation"][1]["novel exon"]}')

plt.rcParams["figure.figsize"] = (20,10)
fig,axs=isoseq['MAPT'].sashimi_figure(x_range=[45960000,46015000])
fig.tight_layout()
plt.show()

#
import glob
import pandas as pd
from pathlib import Path
import ipdb

def is_nonempty(file):
  assert Path(file).stat.st_size
def is_over_size(file,n):
  assert Path(file).stat.st_size > n
def newfolder(file,newroot):
  file = Path(file)
  assert snakedir in file.parents , str(snakedir) + " doesn't seem to be in parents of " + str(file)
  return str(Path(*(newroot,)+file.relative_to(snakedir).parts[1:]))

shell.executable("bash")
shell.prefix("set -e  pipefail;")
# user set parameter

configfile: "../src/config.yaml"
config['root'] = Path(config['root'])


snakedir = Path(config['root']).resolve() / 'pipeline'
TMPDIR = Path('../tmp')


seqfilesdf = pd.read_csv(config['sample_files']).set_index("sample_id", drop=False)
sampledf = pd.read_csv(config['sample_parameter']).set_index("sample_id", drop=False)

assert sampledf.sample_id.is_unique
assert seqfilesdf.file.is_unique
assert set(seqfilesdf.sample_id) == set(seqfilesdf.sample_id)
for sample in sampledf.sample_id:
  for f in seqfilesdf.loc[[sample],'file']:
    assert snakedir/'input'/sample in Path(f).parents, 'This path should be input/sample/file :'+f


samples = list(sampledf['sample_id'].unique())
fastqs = list(seqfilesdf['file'].unique())
ribosamples = sampledf.sample_id[sampledf['assay']=='ribo']


#local copies of the annotation
REF = snakedir / Path(config['REF_orig']).with_suffix('.fa').name
GTF = snakedir / Path(config['GTF_orig']).with_suffix('.gtf').name
GFF = snakedir / Path(config['GFF_orig']).with_suffix('.gff3').name
CDSGTF = GTF.with_suffix('.cds.gtf')
RNAFASTA = GTF.with_suffix('.fa')
CODINGFASTA=GTF.with_suffix('.coding.fa')
PROTEINFASTA=GTF.with_suffix('.protein.fa')
CDSFASTA=GTF.with_suffix('.cds.fa')
BED=GTF.with_suffix('.bed')


rule all:
  input:
    seqfilesdf.file.unique(),
    expand("processed_reads/{sample}/.done", sample = sampledf.sample_id.unique()),
    # expand("fastqc/data/{sample}/.done", sample = samples),
    # expand("star/data/{sample}/.done", sample = samples),
    # expand("qc/data/{sample}/.done", sample = samples),
    # ("multiqc/multiqc_report.html"),
    # expand("feature_counts/data/{sample}/feature_counts", sample = samples),
    # expand("feature_counts/all_feature_counts"),
    # # expand("bigwigs/{group}/{strand}/{istrans}.done",group = samples,strand=strands,istrans=istransvals),
    # # expand("mergedbigwigs/{group}/{strand}/{istrans}.done",group = GROUPS,strand=strands,istrans=istransvals),
    # expand('riboqc/reports/{sample}/riboqcreport.html', sample = ribosamples+groupnames),
    # expand('groupedsatan/{group}.fasta', group = groupnames),


rule link_in_ref:
  input: config['REF_orig']
  output: REF
  shell:r"""
      ln -fs {config['REF_orig']} {REF}
      """

rule link_in_files:
  input: 'input/{sample}/{fastq}'
  output: 'preprocessed_reads/{sample}/{fastq}'
  run:  
    sample = wildcards['sample']
    fastq = wildcards['fastq']
    shell(r"""
      mkdir -p $(dirname {output})
      ln -sf $(readlink -f input/{sample}/{fastq}) {output}
    """)


rule cutadapt_reads:
  input: 'preprocessed_reads/{sample}/{fastq}'
  output: 'cutadapt_reads/{sample}/{fastq}'
  conda: '../envs/cutadapt.yml'
  params: MINREADLENGTH=config['MINREADLENGTH'],MAXREADLENGTH=config['MAXREADLENGTH'],QUALLIM=config['QUALLIM']
  shell: r"""    #   set -evx
      set -e
       mkdir -p cutadapt_reads/{wildcards.sample}/
        zcat {input} \
           | cutadapt \
             -a TGGAATTCTCGGGTGCCAAGG \
            --minimum-length {params.MINREADLENGTH} \
            --maximum-length {params.MAXREADLENGTH} \
            -q {params.QUALLIM} - \
        2> cutadapt_reads/{wildcards.sample}/{wildcards.fastq}.cutadaptstats.txt \
        | gzip  > {output}
"""

rule collapse_reads:
    input: 'cutadapt_reads/{sample}/{fastq}'
    output: 'collapse_reads/{sample}/{fastq}'
    run:
        sample = wildcards['sample']
        shell(r"""
       set -evx
     
       mkdir -p collapse_reads/{sample}/
     
       zcat {input}  \
         | ~/work/bin/collapse_reads.pl {wildcards.sample} \
         2> collapse_reads/{wildcards.sample}/{wildcards.fastq}.collreadstats.txt \
         | cat > {output}
     """)
        is_over_size(output[0],100)

#this finds small filse
 # find collapse_reads/ -name "*.fastq.gz" -size -100M
# find trim_reads/ -name "*.fastq.gz" -size -10M  | xargs ls -latr
# #this finds everything in a certain rule that's less than 10M and then quits
# for i in $(find trim_reads/ -name "*.fastq.gz" -size -10M);do   find . -name $(dirname $i | xargs basename) | grep -v input | grep -v cutadapt; done
rule trim_reads:
    input: 'collapse_reads/{sample}/{fastq}'
    output: 'trim_reads/{sample}/{fastq}'
    run:
        sample = wildcards['sample']
        shell(r"""
       set -evx
     
       OUTDIR=$(dirname {output})
       mkdir -p  $OUTDIR
     
       {REMOVE8NBIN} {input} {output}

       gzip -f {output}
       mv {output}.gz {output}

     """)



rule make_trna_rrna_indices:
  input: GTF,contaminants="../contaminants/contaminants.fa"
  output: touch('tRNA_rRNA_index/tRNA_rRNA_index.done')
  run:
    outprefix = output[0].replace('.done','')
    fafile =outprefix+'.fa'
    
    shell(r"""
      source activate tophat
      #get rRNAs and tRNAs, rename the tRNAs to exons for GenePRed,
      #get the gene type out and stick it front of the transcript id
      #for better names in the fasta
      cp {input.contaminants} {fafile}

       bowtie2-build {fafile} {outprefix} -p {threads}
       
      """)

rule make_bowtie_indices:
  input: REF
  output: touch('bowtie_index/.done')
  threads: 8
  run:
    outprefix = output[0].replace('.done','')
    fafile =outprefix+'.fa'
    shell(r"""
      #get rRNAs and tRNAs, rename the tRNAs to exons for GenePRed,
      #get the gene type out and stick it front of the transcript id
      #for better names in the fasta
       bowtie2-build {REF} bowtie_index/ref 
       cp {REF} bowtie_index/ref.fa
       bowtie2-build {RNAFASTA} bowtie_index/transcriptome
       cp {RNAFASTA} bowtie_index/transcriptome.fa 
       
      """)

collate_idxscript = "../exploration/collate_idx.R"

rule filter_tRNA_rRNA:
    input: 
      'trim_reads/{sample}/{fastq}',
      'tRNA_rRNA_index/tRNA_rRNA_index.done'  
      # filter_index    
    output: 'filter_reads/{sample}/{fastq}',
    threads: 8
    run:
      sample = wildcards['sample']
      indexname = input[1].replace('.done','')
      outdir = os.path.dirname(output[0])
      shell(r"""
       set -evx

       [ -f {outdir} ] && rm -rf {outdir}
     
       mkdir -p  {outdir}

      bowtie2 \
        -x {indexname} \
        -L 20  \
        -p {threads}  \
        -N 0 \
        -U  {input[0]} \
        --un-gz {output[0]} \
        --no-unal \
        2> {output[0]}.alignreport.log > {output[0]}.filtered_reads.sam


        samtools view -bh  {output[0]}.filtered_reads.sam \
        | samtools sort -@ {threads}  > {output[0]}.filtered_reads.bam

      samtools index {output[0]}.filtered_reads.bam
      
      #those which mismatch twice should not be included
      samtools view -hb {outdir}/filtered_reads.bam \
      | bamtools filter -tag XM:2-10 -in - -out /dev/stdout \
      | samtools view -H > {output[0]}.mm.sam
      #>> {outdir}/unmapped.sam
     
      #group the idx columns stuff is from 
      samtools idxstats {output[0]}.filtered_reads.bam \
      | perl -lanpe 's/^(\S+)_[^_\s]+\t/$1\t/' > {output[0]}.idxtmp

      #Rscript --vanilla {collate_idxscript} {output[0]}.idxtmp {indexname}.fa

      samtools stats {output[0]}.filtered_reads.bam > samtools stats {output[0]}.filtered_reads.bam.stats

    """)



def get_processed_files(wc): 
  if wc['sample'] in ribosamples:
    return [ newfolder(fq,'filter_reads') for fq in seqfilesdf.loc[ [wc['sample']],'file' ] ]
  else:
    return list(seqfilesdf['file'][wc['sample']])
#this rule is the 'signal spliter where we go from sample to indiv fastqs
rule link_processed_reads:
  input: get_processed_files
  output: touch('processed_reads/{sample}/.done')
  run:
    shell(r"""
        mkdir -p processed_reads/{wildcards.sample}
        ln -rifs $(readlink -f {input}) processed_reads/{wildcards.sample}/
    """)

rule fastqc:
     input: 'processed_reads/{sample}/.done'
     output: touch('fastqc/data/{sample}/.done')
     threads: 4
     log:'fastqc/reports/{sample}/fastqc.log'
     params:
      reads = lambda wc: [fq.replace('input/','processed_reads/') for fq in seqfilesdf['file'][wc['sample']]],
      outdir = lambda wc: 'fastqc/data/'+wc.sample+'/'
     shell: '''
          OUTDIR=$(dirname {output[0]})
          mkdir -p {params.outdir}
          wait $(for i in {params.reads}; do $( fastqc -o {params.outdir} $i ) & done) 
        '''

##################
#process the annotation
##################

rule gffread:
  input: REF=config['REF_orig'],GTF=config['GTF_orig'],GFF=config['GFF_orig']
  output: REF,GTF,CDSGTF,RNAFASTA,CDSFASTA,BED,GFF
  conda: '../envs/gffread.yml'
  shell: r""" 
      ln -s {input.REF} {REF}
      # set -x
      #with filtering output all sequences
      cat {input.GTF} \
      | sed -r  's/((transcript_id|gene_id|protein_id|ID|Parent|exon_id|havana_gene|havana_transcript)\W+\w+)\.[0-9]+/\1/g' \
      > {GTF}

      cat {input.GFF} \
      | sed -r  's/((transcript_id|gene_id|protein_id|ID|Parent|exon_id|havana_gene|havana_transcript)\W+\w+)\.[0-9]+/\1/g' \
      > {GFF}


      #needs gff - output exon sequences
      cat {input.GFF} |  grep -P -e'\texon\t|^\#' | gffread - -F -E -g {REF} -W -w {RNAFASTA} -o /dev/null

      #Note we are now minus the transcript and exon entries for these
      #now make GTF

      #| grep -P -e'\tCDS\t|^\#' 
     #with filtering, output the coding sequences filteirng out the ones that aren't in frame, have a stop codon, are pseudogenes etc.
      
      cat {input.GFF}  \
        | sed -r  's/((transcript_id|gene_id|protein_id|ID=\w+|Parent)\W+\w+)\.[0-9]+/\1/g' \
        | gffread - -C -V -J --no-pseudo  -F -E -g {REF} \
        -W -w {CODINGFASTA} -x {CDSFASTA} -y {PROTEINFASTA} -T \
        -o /dev/stdout \
        | awk -v FS="\t" -v OFS="\t" '{{if($3=="CDS"){{$3="exon";print $0}}}}' \
         > {CDSGTF}

      #now make bed
      cat {input.GTF} | awk '{{print $1,$4,$5,"name",$6,$7}}' > {BED}
    """

rule star_index:
  input: REF=REF,GTF=GTF
  output: touch('starindex/.done')
  threads: 8
  run:
    shell(r"""
      STAR \
      --runThreadN {threads} \
      --runMode genomeGenerate \
      --genomeDir $(dirname {output}) \
      --gtf {input.GTF} \
      --genomeFastaFiles {input.REF}
      """)  

def get_fastqops(inputdir,read_pattern,lstring='<( zcat ',rstring=')'):
  import glob as glob
  #get our input files, in either paired end or single end form
  assert '-f' in read_pattern
  fastqops = read_pattern.split('-')[1].replace('f ','',1).strip()
  fastqops = glob.glob(inputdir+'/'+fastqops)
  fastqops.sort()
  assert fastqops
  assert all([os.stat(fastq).st_size!=0 for fastq in fastqops])
  fastqops = ' '.join(fastqops)
  
  fastqops = lstring+fastqops+')'

  if '-q' in read_pattern:
    fastqops2 = read_pattern.split('-')[2].replace('q ','',1).strip()
    fastqops2 = glob.glob(inputdir+'/'+fastqops2)
    assert all([os.stat(fastq).st_size!=0 for fastq in fastqops2])
    fastqops2.sort()
    fastqops2 = ' '.join(fastqops2)
    fastqops2 = lstring+fastqops2+')'
    fastqops += ' ' + fastqops2
  return(fastqops)


    


rule star:
     input:
          fastqs='processed_reads/{sample}/.done',
          STARINDEX='starindex/.done',
          bowtie_index='bowtie_index/.done',
     output:
          done = touch('star/data/{sample,[^/]+}/.done')
     threads: 8
     run:
          input.STARINDEX=input.STARINDEX.replace('.done','')
          markdup = '' if ASSAY_DICT[wildcards['sample']] == 'ribo' else '-m'
          platform = 'NotSpecified'
          inputdir = os.path.dirname(input['fastqs'])
          outputdir = os.path.dirname(output[0])
          read_pattern = READ_PATTERN_DICT[wildcards['sample']]
          fastqops = get_fastqops(inputdir,read_pattern,lstring='<( zcat ',rstring=')')
          repdir = outputdir.replace('data','reports')
          tophatindex =input['bowtie_index'].replace('.done','')
          
          halfthreads = threads/2
          sortmem = str(int(5000/halfthreads))+'M'

          remap = '1' if ASSAY_DICT[wildcards['sample']] == 'ribo' else ''

          sample = wildcards['sample']
          shell(r"""
            set -x
         MY_TMP_DIR=$(mktemp -d)
        trap "set -x; rm -rf ${{MY_TMP_DIR}}" EXIT KILL TERM INT HUP

         mkdir -p $MY_TMP_DIR
        mkdir -p $MY_TMP_DIR/star
        mkdir -p $MY_TMP_DIR/tophat2

        #--outSAMmultNmax 20 --winAnchorMultimapNmax 50 --outFilterMultimapNmax 20 \

        STAR \
              --genomeDir {input.STARINDEX} \
              --runThreadN {threads} \
              --outSAMunmapped Within \
              --outFilterType BySJout \
              --outMultimapperOrder Random \
              --alignSJoverhangMin 8 \
              --alignSJDBoverhangMin 1 \
              --outFilterMismatchNmax 999 \
              --outFilterMismatchNoverLmax 0.04 \
              --alignIntronMin 20 \
              --alignIntronMax 1000000 \
              --alignMatesGapMax 1000000 \
              --genomeLoad NoSharedMemory \
              --quantMode GeneCounts \
              --outSAMattributes NH HI AS NM MD \
              --outSAMtype BAM  Unsorted\
              --outSAMattrRGline \"ID:{sample}\" \"SM:{sample}\" \"PL:{platform}\" \
              --outFileNamePrefix ${{MY_TMP_DIR}}/star/ \
              --outReadsUnmapped Fastx \
              --readFilesIn {fastqops}
          

          iftophat=
          if [ {remap} -eq 1 ] && [ ${{MY_TMP_DIR}}/star/Unmapped.out.mate* ]; then
            iftophat=true

            tophat2 \
                -p {threads} \
                -z0 \
                -g 100 \
                --output-dir ${{MY_TMP_DIR}}/tophat2 \
                --library-type fr-unstranded \
                --no-coverage-search \
                --transcriptome-index {tophatindex}/transcriptome \
                {tophatindex}/ref \
                ${{MY_TMP_DIR}}/star/Unmapped.out.mate*

            umapped=${{MY_TMP_DIR}}/tophat2/unmapped.bam
            tmapped=
            
            samtools merge \
              -@ {threads}  -f \
             ${{MY_TMP_DIR}}/all.bam \
             ${{MY_TMP_DIR}}/star/Aligned.out.bam \
             ${{MY_TMP_DIR}}/tophat2/*.bam
          else
            cp ${{MY_TMP_DIR}}/star/Aligned.out.bam ${{MY_TMP_DIR}}/all.bam
          fi
          
         samtools sort \
          -@ {halfthreads}\
          -m {sortmem} \
          -T ${{MY_TMP_DIR}} \
          -o {outputdir}/{sample}.bam \
          ${{MY_TMP_DIR}}/all.bam
      
        samtools index {outputdir}/{sample}.bam 

        mkdir -p {repdir}
        samtools stats {outputdir}/{sample}.bam > {repdir}/{sample}.bamstats.txt
        samtools flagstat {outputdir}/{sample}.bam > {repdir}/{sample}.flagstat.log
        samtools idxstats {outputdir}/{sample}.bam > {repdir}/{sample}.idxstats.log
        
        cp  ${{MY_TMP_DIR}}/star/ReadsPerGene.out.tab {outputdir}/ReadsPerGene.out.tab
        cp  ${{MY_TMP_DIR}}/star/SJ.out.tab {outputdir}/
        cp  ${{MY_TMP_DIR}}/star/{{Log.final.out,Log.out}} {repdir}/
        if [ $iftophat ] ;then cp ${{MY_TMP_DIR}}/tophat2/align_summary.txt {repdir} ;fi

          """)

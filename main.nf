#!/usr/bin/env nextflow

log.info """\
    ================================================================
    NF-varmit  –  Iterative WGS variant calling and consensus
    https://github.com/MozesBlom/nf-variant
    Authors: Mozes P.K. Blom ; Jon J. Hoffman ; Dylan DeBaun
    ================================================================
    indivs file   : ${params.indivs_file}
    chromos file  : ${params.chromos_file}
    reference     : ${params.ref_file}
    reads dir     : ${params.readsdir}
    output dir    : ${params.outputdir}

    Sex chromos   : ${params.sex_chromos}
    Auto min len  : ${params.auto_chromos_min_length} bp

    Rounds        : ${params.n_rounds}
    Start round   : ${params.start_round}

    Coverage stats       : ${params.calc_coverage}
    Variant calling      : ${params.indiv_var_call}
    Consensus calling    : ${params.call_consensus}
    Missing-data stats   : ${params.calc_missing_data}

    Variant filtering
    -----------------
    QUAL min      : ${params.qual_min}
    Allele balance: ${params.min_ab} – ${params.max_ab}
    Filter indels : ${params.filt_indels}
    Mask hets     : ${params.mask_hets}
    Mask coverage : ${params.mask_cov}

    Depth thresholds (auto-calculated from round-1 coverage):
      min = max(5, floor(mean × ${params.min_depth_multiplier}))
      max = ceil(mean × ${params.max_depth_multiplier})
    ================================================================
    """.stripIndent()

/*
 * ─────────────────────────────────────────────────────────────────────────────
 *  Input channels and validation
 * ─────────────────────────────────────────────────────────────────────────────
 */

params.n_rounds              = params.n_rounds              ? params.n_rounds.toInteger()    : 4
params.auto_chromos_min_length = params.auto_chromos_min_length ?: 0
params.min_depth_multiplier  = params.min_depth_multiplier  ?: 0.5
params.max_depth_multiplier  = params.max_depth_multiplier  ?: 2.0
params.qual_min              = params.qual_min              ?: 20
params.min_ab                = params.min_ab                ?: 0.2
params.max_ab                = params.max_ab                ?: 0.8

if (params.n_rounds < 1) error "params.n_rounds must be >= 1 (got: ${params.n_rounds})"

def ref_file_obj = file(params.ref_file)
if (!ref_file_obj.exists()) error "Reference not found: ${params.ref_file}"

indivs    = file(params.indivs_file).readLines()
indivs_ch = Channel.fromList(indivs)

chromos    = file(params.chromos_file).readLines()
chromos_ch = Channel.fromList(chromos)

ref_ch = Channel.value(ref_file_obj)

params.reads_suffix1 = params.reads_suffix1 ?: '_R1.fastq.gz'
params.reads_suffix2 = params.reads_suffix2 ?: '_R2.fastq.gz'

reads_ch = indivs_ch.map { indiv ->
    def r1 = file("${params.readsdir}/${indiv}${params.reads_suffix1}")
    def r2 = file("${params.readsdir}/${indiv}${params.reads_suffix2}")
    if (!r1.exists()) error "Missing R1 reads for ${indiv}: ${r1}"
    if (!r2.exists()) error "Missing R2 reads for ${indiv}: ${r2}"
    tuple(indiv, r1, r2)
}

/*
 * ─────────────────────────────────────────────────────────────────────────────
 *  Output directory layout
 *
 *  01.bams/            – sorted, indexed BAM for every individual × round
 *  02.variants/        – per-chromosome VCFs (and mask TSVs) for every round
 *  03.consensus/       – final consensus FASTAs and missing-data stats
 *  04.consensus_refs/  – per-individual FASTA references built after each
 *                        intermediate round (input for the next round)
 *  05.iter_consensus/  – per-chromosome dominant-allele FASTAs used to build
 *                        intermediate references (diagnostic; one file per
 *                        individual × round × chromosome)
 *  06.coverage/        – per-individual depth-threshold TSVs and summary plots
 *  07.iter_stats/      – per-round reference-bias TSVs and combined summary
 * ─────────────────────────────────────────────────────────────────────────────
 */

/*
 * ─────────────────────────────────────────────────────────────────────────────
 *  Processes
 * ─────────────────────────────────────────────────────────────────────────────
 */

process index_ref {

    label 'Mapping'
    tag "${reference.simpleName}"

    input:
    path reference

    output:
    tuple path(reference),
          path("${reference}.bwt.2bit.64"),
          path("${reference}.0123"),
          path("${reference}.pac"),
          path("${reference}.ann"),
          path("${reference}.amb"),
          path("${reference}.fai")

    script:
    """
    samtools faidx ${reference}
    bwa-mem2 index ${reference}
    """
}

// ── BWA-MEM2 mapping ──────────────────────────────────────────────────────────

process bwa_map {

    label 'Mapping'
    tag "${individual} round${round}"
    publishDir "${params.outputdir}/01.bams", mode: 'copy',
               saveAs: { fn -> "${individual}_round${round}.${fn.tokenize('.')[-1]}" }

    input:
    tuple val(individual),
          path(read1), path(read2),
          path(reference),
          path(bwt), path(zero123), path(pac), path(ann), path(amb), path(fai),
          val(round)

    output:
    tuple val(individual), val(round), path("${individual}.bam"), path("${individual}.bam.bai")

    script:
    """
    bwa-mem2 mem -t ${task.cpus} ${reference} ${read1} ${read2} \
        | samtools sort -@ ${task.cpus} -m 2G -o ${individual}.bam

    samtools index -@ ${task.cpus} ${individual}.bam
    """
}

// ── Coverage estimation ───────────────────────────────────────────────────────

process cov_estimate {

    tag "${indiv} ${chromo}"

    input:
    tuple val(indiv), path(indiv_bam), path(indiv_bam_bai), val(chromo), path(reference)

    output:
    tuple val(indiv), path("${chromo}.tsv")

    script:
    """
    samtools coverage --region ${chromo} --output ${chromo}.tsv ${indiv_bam}
    """
}

process cov_summary_INDIV {

    tag "${indiv}"
    publishDir "${params.outputdir}/06.coverage", mode: 'copy'

    input:
    tuple val(indiv), path(chromo_cov_tsv_list)

    output:
    path("${indiv}_coverage.tsv")

    script:
    """
    python \$NF_VARMIT_BIN/00_cov_stats_INDIV.py \
      -c ${chromo_cov_tsv_list} \
      -i ${indiv}
    """
}

process cov_summary_ALL {

    tag "all individuals"
    publishDir "${params.outputdir}/06.coverage", mode: 'copy'

    input:
    path(indivs_cov_tsv_list)

    output:
    file('*')

    script:
    """
    python \$NF_VARMIT_BIN/01_cov_stats_SUMMARY.py \
      -c ${indivs_cov_tsv_list} \
      -s ${params.sex_chromos} \
      -p ${params.plots_per_row}
    """
}

// ── Depth thresholds ──────────────────────────────────────────────────────────

process calc_depth_thresholds {

    tag "${individual}"
    publishDir "${params.outputdir}/06.coverage", mode: 'copy'

    input:
    tuple val(individual), path(chromo_cov_tsv_list)

    output:
    tuple val(individual), path("${individual}_depth_thresholds.tsv")

    script:
    """
    python \$NF_VARMIT_BIN/00_calc_depth_thresholds.py \
      -c ${chromo_cov_tsv_list} \
      -i ${individual} \
      -s ${params.sex_chromos} \
      -l ${params.auto_chromos_min_length} \
      --min_depth_multiplier ${params.min_depth_multiplier} \
      --max_depth_multiplier ${params.max_depth_multiplier}
    """
}

// ── Variant calling ───────────────────────────────────────────────────────────

process call_variants_CHROMO {

    label 'Endurance'
    tag "${individual} round${round} ${chromo}"
    publishDir "${params.outputdir}/02.variants", mode: 'copy',
               saveAs: { fn -> "${individual}/round${round}/${fn}" }

    input:
    tuple val(individual), val(round), path(indiv_bam), path(indiv_bam_bai), val(chromo), path(reference)

    output:
    tuple val(individual), val(round), val(chromo), path("${chromo}_vars_filt.vcf.gz")

    script:
    """
    freebayes -f ${reference} --region ${chromo} -m 10 -p 2 ${indiv_bam} | \
      vcffilter -f "QUAL < ${params.qual_min}" \
                -f "( AB > 0 ) & ( AB < ${params.min_ab} )" --invert --or | \
      vcfallelicprimitives -k -g | \
      bgzip -c > ${chromo}_vars_filt.vcf.gz
    """
}

process remove_indels {

    tag "${individual} round${round} ${chromo}"
    publishDir "${params.outputdir}/02.variants", mode: 'copy',
               saveAs: { fn -> "${individual}/round${round}/${fn}" }

    input:
    tuple val(individual), val(round), val(chromo), path(var_vcf)

    output:
    tuple val(individual), val(round), val(chromo), path("${chromo}_vars_filt_noindels.vcf.gz")

    script:
    """
    bgzip -d -c ${var_vcf} | \
      vcffilter -f 'TYPE = ins' -f 'TYPE = del' -f 'TYPE = complex' --invert --or | \
      bgzip -c > ${chromo}_vars_filt_noindels.vcf.gz
    """
}

// ── Masking ───────────────────────────────────────────────────────────────────

process mask_cov {

    tag "${individual} round${round} ${chromo}"
    publishDir "${params.outputdir}/02.variants", mode: 'copy',
               saveAs: { fn -> "${individual}/round${round}/${fn}" }

    input:
    tuple val(individual), val(round), path(indiv_bam), path(indiv_bam_bai),
          val(chromo), path(reference), path(depth_thresholds)

    output:
    tuple val(individual), val(round), val(chromo), path("${chromo}_cov.tsv")

    script:
    """
    is_sex=0
    for sc in \$(echo "${params.sex_chromos}" | tr ',' ' '); do
        [ "${chromo}" = "\$sc" ] && is_sex=1 && break
    done

    if [ "\$is_sex" -eq 1 ]; then
        min_cov=\$(awk -F'\\t' '\$1=="sex"  {print \$2}' ${depth_thresholds})
        max_cov=\$(awk -F'\\t' '\$1=="sex"  {print \$3}' ${depth_thresholds})
    else
        min_cov=\$(awk -F'\\t' '\$1=="auto" {print \$2}' ${depth_thresholds})
        max_cov=\$(awk -F'\\t' '\$1=="auto" {print \$3}' ${depth_thresholds})
    fi

    samtools depth -aa -Q 10 -r ${chromo} ${indiv_bam} | \
      awk -v min="\$min_cov" -v max="\$max_cov" \
        '(\$3 < min || \$3 > max) {print \$1,\$2}' > ${chromo}_cov.tsv
    """
}

process mask_hets {

    tag "${individual} round${round} ${chromo}"
    publishDir "${params.outputdir}/02.variants", mode: 'copy',
               saveAs: { fn -> "${individual}/round${round}/${fn}" }

    input:
    tuple val(individual), val(round), val(chromo), path(var_vcf)

    output:
    tuple val(individual), val(round), val(chromo), path("${chromo}_hets.tsv")

    script:
    """
    bgzip -d -c ${var_vcf} | \
      vcffilter -f "( AF < 1 ) & ( AB < ${params.max_ab} )" | \
      cut -f 1,2 > ${chromo}_hets.tsv
    """
}

process mask_merge {

    tag "${individual} round${round} ${chromo}"
    publishDir "${params.outputdir}/02.variants", mode: 'copy',
               saveAs: { fn -> "${individual}/round${round}/${fn}" }

    input:
    tuple val(individual), val(round), val(chromo), path(het_bed), path(cov_bed)

    output:
    tuple val(individual), val(round), val(chromo), path("${chromo}_cov_hets.tsv")

    script:
    """
    sed -i '/#CHROM/d' ${cov_bed}
    sed -i '/##/d' ${het_bed}
    sed -i 's/ /\\t/g' ${cov_bed}
    sed -i 's/ /\\t/g' ${het_bed}

    cat ${cov_bed} ${het_bed} | sort -Vk1 -Vk2 | uniq > ${chromo}_cov_hets.tsv

    sed -i '1i #CHROM\\tPOS' ${chromo}_cov_hets.tsv
    """
}

// ── Iterative consensus (intermediate rounds) ─────────────────────────────────

process build_iter_consensus {

    tag "${individual} round${round} ${chromo}"
    publishDir "${params.outputdir}/05.iter_consensus", mode: 'copy',
               saveAs: { fn -> "${individual}/round${round}/${fn}" }

    input:
    tuple val(individual), val(round), val(chromo), path(vcf_fn), path(cov_mask), path(reference)

    output:
    tuple val(individual), val(round), val(chromo), path("${individual}_${chromo}_iter_cons.fa")

    script:
    if (round == 1) {
        """
        sed 's/ /\\t/g' ${cov_mask} > bad_positions.tsv

        bgzip -d -c ${vcf_fn} | \
          awk 'NR==FNR {bad[\$1"\\t"\$2]=1; next}
               /^#/   {print; next}
               !(\$1"\\t"\$2 in bad)' \
          bad_positions.tsv - | \
          vcffilter -f 'AF >= 0.5' | \
          bgzip -c > ${chromo}_dominant.vcf.gz

        tabix -p vcf ${chromo}_dominant.vcf.gz

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${chromo}_dominant.vcf.gz \
            -o ${individual}_${chromo}_iter_cons.fa
        """
    } else {
        """
        tabix -p vcf ${vcf_fn}
        bcftools +fixref ${vcf_fn} -Oz -o ${chromo}_fixed.vcf.gz -- \
          -f ${reference} -m flip -d
        tabix -p vcf ${chromo}_fixed.vcf.gz

        sed 's/ /\\t/g' ${cov_mask} > bad_positions.tsv

        bgzip -d -c ${chromo}_fixed.vcf.gz | \
          awk 'NR==FNR {bad[\$1"\\t"\$2]=1; next}
               /^#/   {print; next}
               !(\$1"\\t"\$2 in bad)' \
          bad_positions.tsv - | \
          vcffilter -f 'AF >= 0.5' | \
          bgzip -c > ${chromo}_dominant.vcf.gz

        tabix -p vcf ${chromo}_dominant.vcf.gz

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${chromo}_dominant.vcf.gz \
            -o ${individual}_${chromo}_iter_cons.fa
        """
    }
}

process build_ref_from_consensus {

    tag "${individual} round${round}"
    publishDir "${params.outputdir}/04.consensus_refs", mode: 'copy',
               saveAs: { fn -> "${individual}/${fn}" }

    input:
    tuple val(individual), val(round), path(cons_list)

    output:
    tuple val(individual), val(round), path("${individual}_round${round}_ref.fa")

    script:
    """
    cat ${cons_list} > tmp.fa

    awk '/^>/ { print; next } { gsub(/[^ACGTacgt]/, "N"); print }' \
      tmp.fa > ${individual}_round${round}_ref.fa

    rm tmp.fa
    """
}

// ── Final consensus ───────────────────────────────────────────────────────────

process call_consensus {

    tag "${individual} round${round} ${chromo}"
    publishDir "${params.outputdir}/03.consensus", mode: 'copy',
               saveAs: { fn -> "${individual}/${fn}" }

    input:
    tuple val(individual), val(round), val(chromo), path(vcf_fn), path(reference)

    output:
    tuple val(individual), val(round), val(chromo), path("${individual}_${chromo}_cons.fa")

    script:
    if (round == 1) {
        """
        tabix -p vcf ${vcf_fn}

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${vcf_fn} -o ${individual}_${chromo}_cons.fa
        """
    } else {
        """
        tabix -p vcf ${vcf_fn}

        bcftools +fixref ${vcf_fn} -Oz -o ${chromo}_fixed.vcf.gz -- \
          -f ${reference} -m flip -d
        tabix -p vcf ${chromo}_fixed.vcf.gz

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${chromo}_fixed.vcf.gz \
            -o ${individual}_${chromo}_cons.fa
        """
    }
}

process call_consensus_MASK {

    tag "${individual} round${round} ${chromo}"
    publishDir "${params.outputdir}/03.consensus", mode: 'copy',
               saveAs: { fn -> "${individual}/${fn}" }

    input:
    tuple val(individual), val(round), val(chromo), path(vcf_fn), path(mask_fn), path(reference)

    output:
    tuple val(individual), val(round), val(chromo), path("${individual}_${chromo}_cons.fa")

    script:
    if (round == 1) {
        """
        tabix -p vcf ${vcf_fn}

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${vcf_fn} -m ${mask_fn} \
            -o ${individual}_${chromo}_cons.fa
        """
    } else {
        """
        tabix -p vcf ${vcf_fn}

        bcftools +fixref ${vcf_fn} -Oz -o ${chromo}_fixed.vcf.gz -- \
          -f ${reference} -m flip -d
        tabix -p vcf ${chromo}_fixed.vcf.gz

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${chromo}_fixed.vcf.gz -m ${mask_fn} \
            -o ${individual}_${chromo}_cons.fa
        """
    }
}

process call_consensus_dominant {

    tag "${individual} round${round} ${chromo} (dominant)"
    publishDir "${params.outputdir}/03.consensus", mode: 'copy',
               saveAs: { fn -> "${individual}/${fn}" }

    input:
    tuple val(individual), val(round), val(chromo), path(vcf_fn), path(reference)

    output:
    tuple val(individual), val(round), val(chromo), path("${individual}_${chromo}_dominant_cons.fa")

    script:
    if (round == 1) {
        """
        tabix -p vcf ${vcf_fn}

        bgzip -d -c ${vcf_fn} | vcffilter -f 'AF >= 0.5' | \
          bgzip -c > ${chromo}_dominant.vcf.gz
        tabix -p vcf ${chromo}_dominant.vcf.gz

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${chromo}_dominant.vcf.gz \
            -o ${individual}_${chromo}_dominant_cons.fa
        """
    } else {
        """
        tabix -p vcf ${vcf_fn}

        bcftools +fixref ${vcf_fn} -Oz -o ${chromo}_fixed.vcf.gz -- \
          -f ${reference} -m flip -d
        tabix -p vcf ${chromo}_fixed.vcf.gz

        bgzip -d -c ${chromo}_fixed.vcf.gz | vcffilter -f 'AF >= 0.5' | \
          bgzip -c > ${chromo}_dominant.vcf.gz
        tabix -p vcf ${chromo}_dominant.vcf.gz

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${chromo}_dominant.vcf.gz \
            -o ${individual}_${chromo}_dominant_cons.fa
        """
    }
}

process call_consensus_dominant_MASK {

    tag "${individual} round${round} ${chromo} (dominant+mask)"
    publishDir "${params.outputdir}/03.consensus", mode: 'copy',
               saveAs: { fn -> "${individual}/${fn}" }

    input:
    tuple val(individual), val(round), val(chromo), path(vcf_fn), path(mask_fn), path(reference)

    output:
    tuple val(individual), val(round), val(chromo), path("${individual}_${chromo}_dominant_cons.fa")

    script:
    if (round == 1) {
        """
        tabix -p vcf ${vcf_fn}

        bgzip -d -c ${vcf_fn} | vcffilter -f 'AF >= 0.5' | \
          bgzip -c > ${chromo}_dominant.vcf.gz
        tabix -p vcf ${chromo}_dominant.vcf.gz

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${chromo}_dominant.vcf.gz \
            -m ${mask_fn} \
            -o ${individual}_${chromo}_dominant_cons.fa
        """
    } else {
        """
        tabix -p vcf ${vcf_fn}

        bcftools +fixref ${vcf_fn} -Oz -o ${chromo}_fixed.vcf.gz -- \
          -f ${reference} -m flip -d
        tabix -p vcf ${chromo}_fixed.vcf.gz

        bgzip -d -c ${chromo}_fixed.vcf.gz | vcffilter -f 'AF >= 0.5' | \
          bgzip -c > ${chromo}_dominant.vcf.gz
        tabix -p vcf ${chromo}_dominant.vcf.gz

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${chromo}_dominant.vcf.gz \
            -m ${mask_fn} \
            -o ${individual}_${chromo}_dominant_cons.fa
        """
    }
}

process call_consensus_iupac {

    tag "${individual} round${round} ${chromo} (IUPAC)"
    publishDir "${params.outputdir}/03.consensus", mode: 'copy',
               saveAs: { fn -> "${individual}/${fn}" }

    input:
    tuple val(individual), val(round), val(chromo), path(vcf_fn), path(reference)

    output:
    tuple val(individual), val(round), val(chromo), path("${individual}_${chromo}_iupac_cons.fa")

    script:
    if (round == 1) {
        """
        tabix -p vcf ${vcf_fn}

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${vcf_fn} --iupac-codes \
            -o ${individual}_${chromo}_iupac_cons.fa
        """
    } else {
        """
        tabix -p vcf ${vcf_fn}

        bcftools +fixref ${vcf_fn} -Oz -o ${chromo}_fixed.vcf.gz -- \
          -f ${reference} -m flip -d
        tabix -p vcf ${chromo}_fixed.vcf.gz

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${chromo}_fixed.vcf.gz --iupac-codes \
            -o ${individual}_${chromo}_iupac_cons.fa
        """
    }
}

process call_consensus_iupac_MASK {

    tag "${individual} round${round} ${chromo} (IUPAC+mask)"
    publishDir "${params.outputdir}/03.consensus", mode: 'copy',
               saveAs: { fn -> "${individual}/${fn}" }

    input:
    tuple val(individual), val(round), val(chromo), path(vcf_fn), path(mask_fn), path(reference)

    output:
    tuple val(individual), val(round), val(chromo), path("${individual}_${chromo}_iupac_cons.fa")

    script:
    if (round == 1) {
        """
        tabix -p vcf ${vcf_fn}

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${vcf_fn} --iupac-codes \
            -m ${mask_fn} \
            -o ${individual}_${chromo}_iupac_cons.fa
        """
    } else {
        """
        tabix -p vcf ${vcf_fn}

        bcftools +fixref ${vcf_fn} -Oz -o ${chromo}_fixed.vcf.gz -- \
          -f ${reference} -m flip -d
        tabix -p vcf ${chromo}_fixed.vcf.gz

        samtools faidx ${reference} ${chromo} | \
          bcftools consensus ${chromo}_fixed.vcf.gz --iupac-codes \
            -m ${mask_fn} \
            -o ${individual}_${chromo}_iupac_cons.fa
        """
    }
}

// ── Missing data ──────────────────────────────────────────────────────────────

process calc_missing_data_INDIV {

    tag "${individual}"
    publishDir "${params.outputdir}/03.consensus", mode: 'copy',
               saveAs: { fn -> "${individual}/${fn}" }

    input:
    tuple val(individual), path(cons_fn_list)

    output:
    file("${individual}_missing_data.tsv")

    script:
    """
    python \$NF_VARMIT_BIN/02_cons_stats_INDIV.py \
      -c ${cons_fn_list} \
      -i ${individual}
    """
}

process calc_missing_data_SUMMARY {

    tag "all individuals"
    publishDir "${params.outputdir}/03.consensus", mode: 'copy'

    input:
    path(cons_fn_list)

    output:
    file("all_indivs_missing_data.tsv")
    file("all_indivs_missing_data.pdf")

    script:
    """
    unset DISPLAY

    python \$NF_VARMIT_BIN/03_cons_stats_SUMMARY.py \
      -c ${cons_fn_list}
    """
}

// ── Reference-bias iteration stats ───────────────────────────────────────────

process calc_iter_stats {

    tag "${individual} round${round}"
    publishDir "${params.outputdir}/07.iter_stats", mode: 'copy',
               saveAs: { fn -> fn }

    input:
    tuple val(individual), val(round), path(indiv_bam), path(indiv_bam_bai), path(vcf_list)

    output:
    tuple val(individual), val(round), path("${individual}_round${round}_iter_stats.tsv")

    script:
    """
    python \$NF_VARMIT_BIN/01_iter_stats.py \
      -v ${vcf_list} \
      -b ${indiv_bam} \
      -i ${individual} \
      -r ${round} \
      -s ${params.sex_chromos} \
      -l ${params.auto_chromos_min_length}
    """
}

process iter_stats_summary {

    tag "all rounds"
    publishDir "${params.outputdir}/07.iter_stats", mode: 'copy'

    input:
    path(stats_tsv_list)

    output:
    path("all_rounds_iter_stats.tsv")

    script:
    """
    set -- ${stats_tsv_list}
    head -n 1 \$1 > all_rounds_iter_stats.tsv

    for f in ${stats_tsv_list}; do
        tail -n +2 \$f
    done | sort -k1,1 -k2,2n -k3,3 >> all_rounds_iter_stats.tsv
    """
}

/*
 * ─────────────────────────────────────────────────────────────────────────────
 *  Sub-workflows
 * ─────────────────────────────────────────────────────────────────────────────
 */

/*
 * INTERMEDIATE ROUNDS (1 … n_rounds-1)
 *
 * Each round:
 *   1. Map reads to the current per-individual reference.
 *   2. Round 1 only: estimate coverage → calculate depth thresholds
 *      (mean × multipliers, separately for autosomes and sex chromosomes).
 *   3. Call SNPs; filter indels unconditionally (clean references need no indels).
 *   4. Mask bad-depth positions using auto-calculated thresholds.
 *   5. Apply dominant SNPs (AF ≥ 0.5) at well-covered positions → zero missing data.
 *   6. Concatenate per-chromosome sequences → new per-individual reference.
 */

workflow iteration {
    take:
    reads_ch             // [indiv, r1, r2]
    current_indiv_ref_ch // [indiv, ref]
    chromos_ch           // [chromo]
    round                // integer
    depth_thresholds_ch  // [indiv, thresholds_tsv]

    main:
    // Round 1: all individuals share the same starting reference → deduplicate before indexing.
    // Round > 1: each individual has their own per-individual consensus reference.
    def ref_index_input = round == 1
        ? current_indiv_ref_ch.map{ it[1] }.distinct()
        : current_indiv_ref_ch.map{ it[1] }

    indexed_ref_ch = ref_index_input | index_ref

    // Join current_indiv_ref_ch with indexed_ref_ch on filename to pair each individual
    // with their own reference.  Carry ref_indexed (the work-dir staged copy) forward —
    // NOT ref_orig — so that .fai / .bwt index files are co-located with the FASTA in
    // every downstream process work directory.
    def indiv_indexed_ref_ch = current_indiv_ref_ch
        .combine(indexed_ref_ch)
        .filter { indiv, ref_orig, ref_indexed, bwt, zero123, pac, ann, amb, fai ->
                  ref_orig.getName() == ref_indexed.getName() }
        .map    { indiv, ref_orig, ref_indexed, bwt, zero123, pac, ann, amb, fai ->
                  tuple(indiv, ref_indexed, bwt, zero123, pac, ann, amb, fai) }

    // BUG FIX: derive a keyed (indiv, ref) channel from the INDEXED copy of the ref,
    // not from current_indiv_ref_ch.  current_indiv_ref_ch carries the original file
    // path (e.g. params.ref_file in round 1, or the publishDir copy in round > 1).
    // Passing that path to bcftools/samtools faidx in downstream processes means those
    // tools look for the .fai index next to the original file — which either doesn't
    // exist (round 1, original ref may be unindexed) or is the wrong round's index.
    // Using the staged work-dir copy guarantees the correct index is always present.
    def indiv_ref_ch = indiv_indexed_ref_ch
        .map { indiv, ref, bwt, zero123, pac, ann, amb, fai -> tuple(indiv, ref) }

    map_input_ch = reads_ch
        .combine(indiv_indexed_ref_ch, by: 0)
        .map { indiv, r1, r2, ref, bwt, zero123, pac, ann, amb, fai ->
            tuple(indiv, r1, r2, ref, bwt, zero123, pac, ann, amb, fai, round)
        }

    bam_ch = bwa_map(map_input_ch)

    def bam_chromo_ref_ch = bam_ch
        .combine(chromos_ch)
        .combine(indiv_ref_ch, by: 0)
        .map { indiv, rnd, bam, bai, chromo, ref ->
            tuple(indiv, rnd, bam, bai, chromo, ref)
        }

    // Coverage + thresholds (round 1 only)
    def updated_thresholds_ch = depth_thresholds_ch
    if (round == 1) {
        cov_input_ch = bam_chromo_ref_ch.map { indiv, rnd, bam, bai, chromo, ref ->
            tuple(indiv, bam, bai, chromo, ref)
        }
        cov_estimate(cov_input_ch)
        updated_thresholds_ch = calc_depth_thresholds(cov_estimate.out.groupTuple())

        if (params.calc_coverage) {
            cov_estimate.out.groupTuple() | cov_summary_INDIV
            cov_summary_INDIV.out.collect() | cov_summary_ALL
        }
    }

    // Variant calling (always remove indels in intermediate rounds)
    vars_ch = call_variants_CHROMO(bam_chromo_ref_ch) | remove_indels

    // Coverage masking
    def bam_chromo_ref_thresh_ch = bam_chromo_ref_ch
        .combine(updated_thresholds_ch, by: 0)
        .map { indiv, rnd, bam, bai, chromo, ref, thresh ->
            tuple(indiv, rnd, bam, bai, chromo, ref, thresh)
        }

    cov_mask_ch = mask_cov(bam_chromo_ref_thresh_ch)

    // Build dominant-allele consensus for use as next-round reference.
    // Use indiv_ref_ch (staged, indexed copy) so samtools faidx finds the .fai.
    iter_cons_input_ch = vars_ch
        .combine(cov_mask_ch, by: [0, 1, 2])
        .combine(indiv_ref_ch, by: 0)
        .map { indiv, rnd, chromo, vcf, cov_mask, ref ->
            tuple(indiv, rnd, chromo, vcf, cov_mask, ref)
        }

    new_ref = build_iter_consensus(iter_cons_input_ch)
        .groupTuple(by: [0, 1])
        .map { indiv, rnd, chromo_list, fasta_list -> tuple(indiv, rnd, fasta_list.flatten()) }
        | build_ref_from_consensus
        | map { indiv, rnd, ref -> tuple(indiv, ref) }

    // Reference-bias stats
    vcfs_by_indiv_ch = vars_ch.groupTuple(by: [0, 1]).map { it -> [it[0], it[1], it[3].flatten()] }
    stats_ch = calc_iter_stats(bam_ch.combine(vcfs_by_indiv_ch, by: [0, 1])).map { it[2] }

    emit:
    new_ref            = new_ref
    updated_thresholds = updated_thresholds_ch
    stats              = stats_ch
}

workflow RUN_ITERATIONS {
    take:
    reads_ch
    ref_ch
    chromos_ch
    current_round
    max_rounds
    thresholds_ch
    stats_acc

    main:
    res = iteration(reads_ch, ref_ch, chromos_ch, current_round, thresholds_ch)
    def next_round      = current_round + 1
    def current_stats   = stats_acc.mix(res.stats)

    if (next_round < max_rounds) {
        RUN_ITERATIONS(
            reads_ch, res.new_ref, chromos_ch,
            next_round, max_rounds,
            res.updated_thresholds, current_stats
        )
        out_ref    = RUN_ITERATIONS.out.final_ref
        out_thresh = RUN_ITERATIONS.out.final_thresholds
        out_stats  = RUN_ITERATIONS.out.final_stats
    } else {
        out_ref    = res.new_ref
        out_thresh = res.updated_thresholds
        out_stats  = current_stats
    }

    emit:
    final_ref        = out_ref
    final_thresholds = out_thresh
    final_stats      = out_stats
}

/*
 * ─────────────────────────────────────────────────────────────────────────────
 *  Main workflow
 * ─────────────────────────────────────────────────────────────────────────────
 */

workflow {

    def rounds          = params.n_rounds as int
    def initial_ref_ch  = reads_ch.map { indiv, r1, r2 -> [indiv, file(params.ref_file)] }
    def all_iter_stats_ch = Channel.empty()

    // ── Intermediate rounds (1 … n_rounds-1) ─────────────────────────────────
    def final_inter_ref, final_inter_thresholds, all_stats

    if (rounds > 1) {
        RUN_ITERATIONS(
            reads_ch, initial_ref_ch, chromos_ch,
            1, rounds,
            Channel.empty(), Channel.empty()
        )
        final_inter_ref        = RUN_ITERATIONS.out.final_ref
        final_inter_thresholds = RUN_ITERATIONS.out.final_thresholds
        all_stats              = RUN_ITERATIONS.out.final_stats
    } else {
        final_inter_ref        = initial_ref_ch
        final_inter_thresholds = Channel.empty()
        all_stats              = Channel.empty()
    }

    // ── Final round (n_rounds) ────────────────────────────────────────────────
    /*
     * Map reads to the best per-individual reference (output of the last
     * intermediate round, or the original reference when n_rounds == 1).
     *
     * When n_rounds == 1 there were no intermediate rounds, so coverage
     * estimation and threshold calculation happen here.
     *
     * The final consensus encodes true variant calls:
     *   - Passing sites (REF or ALT with sufficient depth) → included
     *   - Sites without a passing call → N / missing data
     *
     * Masking and indel handling are user-configurable for this round only.
     */
    {
        def round = rounds

        log.info "===================================="
        log.info "  FINAL ROUND (round ${rounds})"
        log.info "===================================="

        // Index per-individual references for the final round.
        // Each individual has a uniquely named consensus ref from the last intermediate
        // round, so no .distinct() is needed or wanted here.
        def final_indexed_ref_ch = final_inter_ref
            .map { indiv, ref -> ref }
            | index_ref

        // BUG FIX (same as iteration workflow): carry the staged, indexed ref path forward —
        // NOT final_inter_ref — so .fai / .bwt files are co-located in every process work dir.
        def final_indiv_indexed_ref_ch = final_inter_ref
            .combine(final_indexed_ref_ch)
            .filter { indiv, ref_orig, ref_indexed, bwt, zero123, pac, ann, amb, fai ->
                      ref_orig.getName() == ref_indexed.getName() }
            .map    { indiv, ref_orig, ref_indexed, bwt, zero123, pac, ann, amb, fai ->
                      tuple(indiv, ref_indexed, bwt, zero123, pac, ann, amb, fai) }

        // Keyed (indiv, ref) channel using the staged copy — used for all downstream
        // processes that need the reference FASTA (variant calling, consensus, masking).
        def final_indiv_ref_ch = final_indiv_indexed_ref_ch
            .map { indiv, ref, bwt, zero123, pac, ann, amb, fai -> tuple(indiv, ref) }

        def map_input_ch = reads_ch
            .combine(final_indiv_indexed_ref_ch, by: 0)
            .map { indiv, r1_reads, r2_reads, ref, bwt, zero123, pac, ann, amb, fai ->
                tuple(indiv, r1_reads, r2_reads, ref, bwt, zero123, pac, ann, amb, fai, round)
            }

        def bam_ch = bwa_map(map_input_ch)

        def bam_chromo_ref_ch = bam_ch
            .combine(chromos_ch)
            .combine(final_indiv_ref_ch, by: 0)
            .map { indiv, rnd, bam, bai, chromo, ref ->
                tuple(indiv, rnd, bam, bai, chromo, ref)
            }

        // Coverage + thresholds (only when n_rounds == 1; otherwise from round 1)
        if (rounds == 1) {
            def cov_input_ch = bam_chromo_ref_ch.map { indiv, rnd, bam, bai, chromo, ref ->
                tuple(indiv, bam, bai, chromo, ref)
            }
            cov_estimate(cov_input_ch)
            final_inter_thresholds = cov_estimate.out.groupTuple() | calc_depth_thresholds

            if (params.calc_coverage) {
                cov_estimate.out.groupTuple() | cov_summary_INDIV
                cov_summary_INDIV.out.collect() | cov_summary_ALL
            }
        }

        // Variant calling
        def vars_filt_ch
        if (params.indiv_var_call && params.filt_indels) {
            vars_filt_ch = call_variants_CHROMO(bam_chromo_ref_ch) | remove_indels
        } else if (params.indiv_var_call && !params.filt_indels) {
            vars_filt_ch = call_variants_CHROMO(bam_chromo_ref_ch)
        } else {
            vars_filt_ch = Channel.empty()
        }

        // Reference-bias stats for the final round
        if (params.indiv_var_call) {
            def vcfs_by_indiv_final_ch = vars_filt_ch
                .groupTuple(by: [0, 1])
                .map { indiv, rnd, chromo_list, vcf_list ->
                    tuple(indiv, rnd, vcf_list.flatten())
                }

            def iter_stats_input_final_ch = bam_ch
                .combine(vcfs_by_indiv_final_ch, by: [0, 1])
                .map { indiv, rnd, bam, bai, vcf_list ->
                    tuple(indiv, rnd, bam, bai, vcf_list)
                }
            all_iter_stats_ch = all_iter_stats_ch
                .mix(calc_iter_stats(iter_stats_input_final_ch).map { indiv, rnd, f -> f })
        }

        // Masking
        def bam_chromo_ref_thresh_ch = bam_chromo_ref_ch
            .combine(final_inter_thresholds, by: 0)
            .map { indiv, rnd, bam, bai, chromo, ref, thresh ->
                tuple(indiv, rnd, bam, bai, chromo, ref, thresh)
            }

        def mask_fn_ch = null
        if (params.indiv_var_call && params.mask_hets && params.mask_cov) {
            def mask_het_ch  = vars_filt_ch | mask_hets
            def mask_cov_ch  = mask_cov(bam_chromo_ref_thresh_ch)
            def mask_comb_ch = mask_het_ch.combine(mask_cov_ch, by: [0, 1, 2])
            mask_fn_ch = mask_merge(mask_comb_ch)
        } else if (params.indiv_var_call && params.mask_hets && !params.mask_cov) {
            mask_fn_ch = vars_filt_ch | mask_hets
        } else if (params.indiv_var_call && !params.mask_hets && params.mask_cov) {
            mask_fn_ch = mask_cov(bam_chromo_ref_thresh_ch)
        }

        // Final consensus
        def consensus_ch = null
        if (params.indiv_var_call && params.call_consensus) {
            if (params.mask_hets || params.mask_cov) {
                def vars_mask_ref_ch = vars_filt_ch
                    .combine(mask_fn_ch, by: [0, 1, 2])
                    .combine(final_indiv_ref_ch, by: 0)
                    .map { indiv, rnd, chromo, vcf, mask, ref ->
                        tuple(indiv, rnd, chromo, vcf, mask, ref)
                    }
                consensus_ch = call_consensus_MASK(vars_mask_ref_ch)
            } else {
                def vcf_ref_ch = vars_filt_ch
                    .combine(final_indiv_ref_ch, by: 0)
                    .map { indiv, rnd, chromo, vcf, ref ->
                        tuple(indiv, rnd, chromo, vcf, ref)
                    }
                consensus_ch = call_consensus(vcf_ref_ch)
            }
        }

        // Dominant-allele and IUPAC consensus outputs (when hets are not masked)
        if (params.indiv_var_call && params.call_consensus && !params.mask_hets) {
            if (params.mask_cov) {
                def vcf_covmask_ref_ch = vars_filt_ch
                    .combine(mask_fn_ch, by: [0, 1, 2])
                    .combine(final_indiv_ref_ch, by: 0)
                    .map { indiv, rnd, chromo, vcf, mask, ref ->
                        tuple(indiv, rnd, chromo, vcf, mask, ref)
                    }
                call_consensus_dominant_MASK(vcf_covmask_ref_ch)
                call_consensus_iupac_MASK(vcf_covmask_ref_ch)
            } else {
                def vcf_ref_extra_ch = vars_filt_ch
                    .combine(final_indiv_ref_ch, by: 0)
                    .map { indiv, rnd, chromo, vcf, ref ->
                        tuple(indiv, rnd, chromo, vcf, ref)
                    }
                call_consensus_dominant(vcf_ref_extra_ch)
                call_consensus_iupac(vcf_ref_extra_ch)
            }
        }

        // Missing-data stats on final consensus
        if (params.calc_missing_data && params.indiv_var_call && params.call_consensus) {
            def md_input_ch = consensus_ch
                .groupTuple(by: 0)
                .map { indiv, rnd_list, chromo_list, cons_list ->
                    tuple(indiv, cons_list.flatten())
                }
            calc_missing_data_INDIV(md_input_ch)
            calc_missing_data_INDIV.out.collect() | calc_missing_data_SUMMARY
        }
    }

    // ── Cross-round reference-bias summary ────────────────────────────────────
    all_iter_stats_ch.collect() | iter_stats_summary
}

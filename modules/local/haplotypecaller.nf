
process HAPLOTYPECALLER {
    tag "${sample_id} - ${interval.getName()}"
    
    conda "bioconda::gatk4=4.6.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0':
        'biocontainers/gatk4:4.6.1.0--py310hdfd78af_0' }"

    
    input:
    tuple val(sample_id), path(bam), path(bam_bai)
    path reference_genome
    path reference_genome_fai
    path reference_genome_dict
    each path (interval) 
    
    output:
    tuple val(sample_id), path("${sample_id}*.g.vcf.gz"), path("${sample_id}*.vcf.gz.tbi"), emit: gvcf
    
    script:
    gvcf_basename = "${sample_id}.${interval.getName()}.g.vcf.gz"
    
    java_options = ""
    if (task.ext.java_options){
        java_options = "--java-options \"${task.ext.java_options}\""
    }
    
    """
    gatk ${java_options} \
        HaplotypeCaller \
        -R ${reference_genome} \
        -I ${bam} \
        -L ${interval} \
        -O ${gvcf_basename} \
        ${task.ext.args}

    """
}

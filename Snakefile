rule all:
    input:
        'success/make-aggregated-object_complete.txt'

rule makeFASTQ:
    input:
        'Snakefile'
    output:
        temp('success/make-FASTQ_complete.txt')
    conda:
        'COMUNEQAID_python.yml'
    threads:
        128
    script:
        'code/01_make-FASTQ.py'

rule mapquant:
    input:
        rules.makeFASTQ.output
    output:
        temp('success/map-quant_complete.txt')
    conda:
        'COMUNEQAID_python.yml'
    threads:
        128
    script:
        'code/02_map-quant.py'

rule makeFASTQtables:
    input:
        rules.mapquant.output
    output:
        temp('success/make-FASTQ-tables_complete.txt')
    conda:
        'COMUNEQAID_R.yml'
    script:
        'code/03_make-FASTQ-tables.R'

rule makereactionobject:
    input:
        rules.makeFASTQtables.output
    output:
        temp('success/make-reaction-object_complete.txt')
    conda:
        'COMUNEQAID_R.yml'
    script:
        'code/04_make-reaction-object.R'

rule makeaggregatedobject:
    input:
        rules.makereactionobject.output
    output:
        temp('success/make-aggregated-object_complete.txt')
    conda:
        'COMUNEQAID_R.yml'
    script:
        'code/05_make-aggregated-object.R'
        
rule makeqcoutput:
    input:
        rules.makeaggregatedobject.output
    output:
        temp('success/make-qc-output_complete.txt')
    conda:
        'COMUNEQAID_R.yml'
    script:
        'code/06_make-QC-output.R'


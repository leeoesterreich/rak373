configfile: "/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/ichorCNA/config38.yaml"
configfile: "/ix1/alee/LO_LAB/Personal/Rahul/cfDNA_project/AURORA/ichorCNA/samples.yaml"

rule all:
  input: 
      expand("results/ichorCNA/{tumor}/{tumor}.cna.seg", tumor=config["samples"]),
      expand("results/readDepth/{samples}.bin{binSize}.wig", samples=config["samples"], binSize=str(config["binSize"])),
      expand("results/genomeWideplot/{tumor}_genomeWide.pdf", tumor=config["samples"])

rule read_counter:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    output:
        "results/readDepth/{samples}.bin{binSize}.wig"        
    params:
        readCounter=config["readCounterScript"],
        binSize=config["binSize"],
        qual="20",
        chrs=config["chrs"]
    resources:
        mem=4
    log:
        "logs/readDepth/{samples}.bin{binSize}.log"
    shell:
        "{params.readCounter} {input} -c {params.chrs} -w {params.binSize} -q {params.qual} > {output} 2> {log}"

rule ichorCNA:
    retries: 2
    input:
        tum="results/readDepth/{tumor}.bin" + str(config["binSize"]) + ".wig",
        #norm=lambda wildcards: "results/readDepth/" + config["pairings"][wildcards.tumor] + ".bin" + str(config["binSize"]) + ".wig"
    output:
        #corrDepth="results/ichorCNA/{tumor}/{tumor}.correctedDepth.txt",
        #param="results/ichorCNA/{tumor}/{tumor}.params.txt",
        cna="results/ichorCNA/{tumor}/{tumor}.cna.seg",
        genomeWidePdf="results/ichorCNA/{tumor}/{tumor}/{tumor}_genomeWide.pdf"
        #segTxt="results/ichorCNA/{tumor}/{tumor}.seg.txt",
        #seg="results/ichorCNA/{tumor}/{tumor}.seg",
        #rdata="results/ichorCNA/{tumor}/{tumor}.RData"
    params:
        outDir="results/ichorCNA/{tumor}/",
        rscript=config["ichorCNA_rscript"],
        id="{tumor}",
        ploidy=config["ichorCNA_ploidy"],
        normal=config["ichorCNA_normal"],
        gcwig=config["ichorCNA_gcWig"],
        mapwig=config["ichorCNA_mapWig"],
        normalpanel=config["ichorCNA_normalPanel"],
        estimateNormal=config["ichorCNA_estimateNormal"],
        estimatePloidy=config["ichorCNA_estimatePloidy"],
        estimateClonality=config["ichorCNA_estimateClonality"],
        scStates=config["ichorCNA_scStates"],
        maxCN=config["ichorCNA_maxCN"],
        includeHOMD=config["ichorCNA_includeHOMD"],
        chrs=config["ichorCNA_chrs"],
        chrTrain=config["ichorCNA_chrTrain"],
        genomeBuild=config["ichorCNA_genomeBuild"],
        genomeStyle=config["ichorCNA_genomeStyle"],
        centromere=config["ichorCNA_centromere"],
        fracReadsChrYMale=config["ichorCNA_fracReadsInChrYForMale"],
        minMapScore=config["ichorCNA_minMapScore"],
        maxFracGenomeSubclone=config["ichorCNA_maxFracGenomeSubclone"],
        maxFracCNASubclone=config["ichorCNA_maxFracCNASubclone"],
        exons=config["ichorCNA_exons"],
        txnE=config["ichorCNA_txnE"],
        txnStrength=config["ichorCNA_txnStrength"],
        plotFileType=config["ichorCNA_plotFileType"],
        plotYlim=config["ichorCNA_plotYlim"],
        libdir=config["ichorCNA_libdir"]
    resources:
        mem=4
    log:
        "logs/ichorCNA/{tumor}.log"    
    shell:
        "mkdir -p {params.outDir} && " \
        "LD_LIBRARY_PATH=\"\" TMPDIR=\"${{TMPDIR:-/tmp}}\" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 " \
        "bash Scripts/ichorCNA_retry_rscript.sh {params.rscript} --id {params.id} --libdir {params.libdir} --WIG {input.tum} --gcWig {params.gcwig} --mapWig {params.mapwig} --normalPanel {params.normalpanel} --ploidy \"{params.ploidy}\" --normal \"{params.normal}\" --maxCN {params.maxCN} --includeHOMD {params.includeHOMD} --chrs \"{params.chrs}\" --chrTrain \"{params.chrTrain}\" --genomeStyle {params.genomeStyle} --genomeBuild {params.genomeBuild} --estimateNormal {params.estimateNormal} --estimatePloidy {params.estimatePloidy} --estimateScPrevalence {params.estimateClonality} --scStates \"{params.scStates}\" --centromere {params.centromere} --exons.bed {params.exons} --txnE {params.txnE} --txnStrength {params.txnStrength} --minMapScore {params.minMapScore} --fracReadsInChrYForMale {params.fracReadsChrYMale} --maxFracGenomeSubclone {params.maxFracGenomeSubclone} --maxFracCNASubclone {params.maxFracCNASubclone} --plotFileType {params.plotFileType} --plotYLim \"{params.plotYlim}\" --outDir {params.outDir} > {log} 2> {log}"

rule export_genomewide_plot:
    input:
        cna="results/ichorCNA/{tumor}/{tumor}.cna.seg",
        pdf_src="results/ichorCNA/{tumor}/{tumor}/{tumor}_genomeWide.pdf"
    output:
        "results/genomeWideplot/{tumor}_genomeWide.pdf"
    shell:
        "mkdir -p results/genomeWideplot && cp {input.pdf_src} {output}"
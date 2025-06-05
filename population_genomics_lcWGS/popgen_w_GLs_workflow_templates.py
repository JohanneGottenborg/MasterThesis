#!/bin/env python3
from gwf import AnonymousTarget
import os, glob

########################## Functions ##########################

def GLs_for_population_structure(bam_paths: str, GL_model: int, output_directory: str, snp_pval: int, minq: int, minmapq: int):
    """
    # SUMMARY
    # Reads BAM files from stru200ind_bams.txt
    # Computes genotype likelihoods using GATK model (-GL 2) --because they use this in Pecnerova and Garcia-Erill
    # Outputs the log genotype likelihoods to a file - 2: beagle genotype likelihood format (-doGlf 2)
    # Major and minor allele are inferred from the data (-doMaf 1)
    # Filters SNPs with a p-value threshold of 1e-6 (-SNP_pval 1e-6)
    # Uses 16 threads for faster processing (-P 16)

    # Filters on minInd, min depth and max depth can be added here
    # minInd 180 only includes sites that have information from at least 180 of the individuals
    # setMinDepth and setMaxDepth filters out sites, where the total depth (combined across all individuals) are below or above the chosen thresholds
    ## Thresholds are based on calculation:
    ## Mean samples coverage = 2.82. Mean total coverage = 202*2.82 ≈ 570
    ## setMaxDepth to 2 * mean total coverage = 2 * 570 = 1140
    ## setMinDepth to 1/3 * mean total coverage = 1/3 * 665 = 190

    """
    inputs = {'bam_paths': bam_paths}
    outputs = {'GLs': f'{output_directory}/202inds_GLs.beagle.gz'}
    options = {
        'cores': 16,
		'memory': '64g',
		'walltime': '72:00:00'
    }
    spec = f"""
	# Sources environment
	if [ "$USER" == "johannegottenborg" ]; then
		source /home/"$USER"/.bashrc
		source activate angsd_env
	fi

    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

    [ -d {output_directory} ] || mkdir -p {output_directory}
    
    angsd -bam {inputs['bam_paths']} -GL {GL_model} \\
        -out '{output_directory}/202inds_GLs' -nThreads 8 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doCounts 1 \\
        -SNP_pval {snp_pval} -minInd 180 -setMinDepth 190 -setMaxDepth 1140 \\
        -minMapQ {minmapq} -minQ {minq}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def PCAngsd_population_structure(GL_file: str, output_directory: str):
    """
    Running PCAngsd on pre-calculated genotype likelihoods.

    Output is covariance matrices and admixture proportions.
    """
    inputs = {'GLs': GL_file}
    outputs = {'default_PCs': f'{output_directory}/PCAngsd_default_PCs.cov',
               '2PCs': f'{output_directory}/PCAngsd_2PCs.cov',
               '1PC': f'{output_directory}/PCAngsd_1PC.cov'}
    options = {
        'cores': 16,
		'memory': '64g',
		'walltime': '3:00:00'
    }
    spec = f"""
	# Sources environment
	if [ "$USER" == "johannegottenborg" ]; then
		source /home/"$USER"/.bashrc
		source activate pcangsd
	fi

    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

    [ -d {output_directory} ] || mkdir -p {output_directory}

    # Default number of PCs
    pcangsd --beagle {inputs['GLs']} --threads 16 --out '{output_directory}/PCAngsd_default_PCs' --admix

    # 2 PCs
    pcangsd --beagle {inputs['GLs']} --eig 2 --threads 16 --out '{output_directory}/PCAngsd_2PCs' --admix
    
    # 1 PC
    pcangsd --beagle {inputs['GLs']} --eig 1 --threads 16 --out '{output_directory}/PCAngsd_1PC' --admix
    
    echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def nucleotide_diversity_pi(bam_paths: str, reference: str, gl_model: int, output_directory: str, minq: int, minmapq: int, minind: int, maxtotaldepth: int, mintotaldepth: int, minhwepval: float):
    """
    SUMMARY
    Reads BAM files from all_black_inds.txt
    Computes genotype likelihoods using GATK model (-GL 2) --because they use this in Pecnerova and Garcia-Erill
    -minHWEpval is chosen based on visual inspection of the HWE p-value distribution and results in a loss of approx. 0.15% of variants.

    MinInd is based on analysis of how many variant sites are called, when including different number of individuals. Out of 202, 180 was chosen as threshold 
    # since this still retained more than 88% of sites. There are 128 blacks, therefore (180/202)*128 ≈ 114

    setMinDepth and setMaxDepth filters out sites, where the total depth (combined across all individuals) are below or above the chosen thresholds
    # Thresholds are based on calculation:
    # Mean samples coverage = 3.6. Mean total coverage = 128*3.6 ≈ 461
    # setMaxDepth to 2 * mean total coverage = 2 * 461 = 922
    # setMinDepth to 1/3 * mean total coverage = 1/3 * 513 = 154
    # Mean samples coverage are from depth analysis

    minMaf is set to 0.005 (≈1/128), corresponding to the minor allele occuring at least twice in the dataset, as recommended in Hemstrøm et al, 2024

    No filtering on SNP_Pval or minMaf since all sites are informative (Lou et al)
    """
    inputs = {'bam_paths': bam_paths}
    outputs = {#'SAFs': f'{output_directory}/angsd_safs.saf.idx',
               'SFS': f'{output_directory}/realSFS.sfs',
               'thetas': f'{output_directory}/thetas_w_fold.thetas.idx',
               'global_pi': f'{output_directory}/global_pi_w_fold.txt'
               }
    options = {
        'cores': 8,
		'memory': '128g',
		'walltime': '48:00:00'
    }
    spec = f"""
	# Sources environment
	if [ "$USER" == "johannegottenborg" ]; then
		source /home/"$USER"/.bashrc
		source activate angsd_env
	fi

    echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

    [ -d {output_directory} ] || mkdir -p {output_directory}

    ## -doSaf 1
    #angsd -bam {inputs['bam_paths']} -doSaf 1 -anc {reference} -GL {gl_model} -P 8 \\
    #    -out '{output_directory}/angsd_safs' -doCounts 1 -doMajorMinor 1 -doHWE 1 -minHWEpval {minhwepval} \\
    #    -minMapQ {minmapq} -minQ {minq} -minInd {minind} -setMaxDepth {maxtotaldepth} -setMinDepth {mintotaldepth}
        

    ## Construct the site frequency spectrum
    #realSFS '{output_directory}/angsd_safs.saf.idx' -P 15 -fold 1 > '{output_directory}/realSFS.sfs'

    ## Calculate theta neutrality statistics for each site
    #realSFS saf2theta '{output_directory}/angsd_safs.saf.idx' -sfs {outputs['SFS']} -fold 1 -outname '{output_directory}/thetas_w_fold'

    ## Calculating pi
    #touch '{output_directory}/thetas_w_fold.txt'
    #touch '{output_directory}/global_pi_w_fold.txt'

    thetaStat print '{output_directory}/thetas_w_fold.thetas.idx' 2>/dev/null > '{output_directory}/thetas_w_fold.txt'

    awk '{{if(NR>1) print $4}}' '{output_directory}/thetas_w_fold.txt' | awk '{{sum+=exp($1)}} END {{print sum/NR}}' >> {outputs['global_pi']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def GLs_filtering_and_calling(bam_paths: str, output_directory: str, gl_model: int, minhwepval: float, minind: int, mintotaldepth: int, maxtotaldepth: int, minq: int, minmapq: int, minmaf: int, snp_pval: float):
    """
    SUMMARY
    Reads BAM files from stru200ind_bams.txt
    Computes genotype likelihoods using GATK model (-GL 2) --because they use this in Pecnerova and Garcia-Erill
    Outputs the log genotype likelihoods to a file - 2: beagle genotype likelihood format (-doGlf 2)
    Major and minor allele are inferred from the data (-doMaf 1)
    Filters SNPs with a p-value threshold of 1e-6 (-SNP_pval 1e-6), strict since we only want the variants we are the most certain about
    Uses 16 threads for faster processing (-P 16)

    Filters on minInd, min depth and max depth can be added here
    -minHWEpval is chosen based on visual inspection of the HWE p-value distribution and results in a loss of approx. 0.15% of variants.

    MinInd is based on analysis of how many variant sites are called, when including different number of individuals. Out of 202, 180 was chosen as threshold 
    # since this still retained more than 88% of sites. There are 128 blacks, therefore (180/202)*128 ≈ 114

    setMinDepth and setMaxDepth filters out sites, where the total depth (combined across all individuals) are below or above the chosen thresholds
    # Thresholds are based on calculation:
    # Mean samples coverage = 3.6. Mean total coverage = 128*3.6 ≈ 461
    # setMaxDepth to 2 * mean total coverage = 2 * 461 = 922
    # setMinDepth to 1/3 * mean total coverage = 1/3 * 513 = 154
    # Mean samples coverage are from depth analysis

    minMaf is set to 0.005 (≈1/128), corresponding to the minor allele occuring at least twice in the dataset, as recommended in Hemstrøm et al, 2024
    """
    inputs = {'bam_paths': bam_paths}
    outputs = {'beagle': f'{output_directory}/angsd_filtered_sites.beagle.gz'}
    options = {
        'cores': 8,
        'memory': '64g',
        'walltime': '24:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "johannegottenborg" ]; then
        source /home/"$USER"/.bashrc
        source activate angsd_env
    fi

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    [ -d {output_directory} ] || mkdir -p {output_directory}

    angsd -bam {inputs['bam_paths']} -GL {gl_model} -out '{output_directory}/angsd_filtered_sites' \\
        -doGlf 2 -doMajorMinor 1 -doMaf 1 -doHWE 1 -doCounts 1 -nThreads 8 \\
        -minHWEpval {minhwepval} -minInd {minind} -setMinDepth {mintotaldepth} -setMaxDepth {maxtotaldepth} \\
        -minMapq {minmapq} -minQ {minq} -minMaf {minmaf} -SNP_pval {snp_pval}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def calculating_linkage_ngsLD(beagle_file: str, bam_paths: str, output_directory: str):
    """
    Calculating LD between each site
    """
    inputs = {'beagle_GLs': beagle_file,
              'individuals': bam_paths}
    outputs = {'pre_filtering_position_file': f'{output_directory}/pre_filtering_positions.txt',
               'beagle_no_header': f'{output_directory}/GLs_no_header.beagle.gz',
               'linkage_file': f'{output_directory}/LD_estimation'}
    options = {
        'cores': 16,
        'memory': '128g',
        'walltime': '24:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "johannegottenborg" ]; then
        source /home/"$USER"/.bashrc
        source activate ngstools_env
    fi

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    [ -d {output_directory} ] || mkdir -p {output_directory}

    zcat {inputs['beagle_GLs']} | awk 'NR>1 {{split($1,a,"_"); print a[1]"_"a[2]"\t"a[3]}}' > {outputs['pre_filtering_position_file']}
    zcat {inputs['beagle_GLs']} | tail -n +2 | cut -f4- | gzip > {outputs['beagle_no_header']}

    NO_INDS=$(wc -l < {inputs['individuals']})
    NO_SITES=$(wc -l < {outputs['pre_filtering_position_file']})

    ngsLD --geno {outputs['beagle_no_header']} --probs true --n_ind $NO_INDS --n_sites $NO_SITES \\
        --out '{output_directory}/LD_estimation' --pos {outputs['pre_filtering_position_file']} --n_threads 15

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def LD_pruning_w_prunegraph(linkage_estimation: str, output_directory: str):
    """
    Filtering sites, such that only unlinked sites remain.
    """
    inputs = {'pre_calculated_linkage': linkage_estimation}
    outputs = {'unlinked_positions': f'{output_directory}/unlinked_positions.pos'}
    options = {
        'cores': 16,
        'memory': '128g',
        'walltime': '5-00:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "johannegottenborg" ]; then
        source /home/"$USER"/.bashrc
        source activate prunegraph_env
    fi

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    [ -d {output_directory} ] || mkdir -p {output_directory}

    prune_graph --header --in {inputs['pre_calculated_linkage']} \\
        --weight-field "r2" --weight-filter "dist <= 50000 && r2 >= 0.5" \\
        --out {outputs['unlinked_positions']} --n-threads 15

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def calling_pruned_SAFs_and_GLs_w_angsd(bam_paths: str, output_directory: str, pruned_positions: str, reference: str, gl_model: int):
    """
    Calculating Site Allele Frequencies (SAFs) and Genotype Likelihoods (GLs) for only sites not in linkage disequilibrium.
    """
    inputs = {'bam_list': bam_paths,
              'pos_file': pruned_positions}
    outputs = {'indexed_pos_file': f'{output_directory}/unlinked_positions_two_cols.pos',
               'beagle_file': f'{output_directory}/pruned_sites_for_selection.beagle.gz',
               'saf_file': f'{output_directory}/pruned_sites_for_selection.saf.idx'}
    options = {
        'cores': 9,
        'memory': '64g',
        'walltime': '24:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "johannegottenborg" ]; then
        source /home/"$USER"/.bashrc
        source activate angsd_env
    fi

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    [ -d {output_directory} ] || mkdir -p {output_directory}

    awk -F: '{{print $1 "\t" $2}}' {inputs['pos_file']} > {output_directory}/unlinked_positions_two_cols.pos

    angsd sites index {output_directory}/unlinked_positions_two_cols.pos

    angsd -bam {inputs['bam_list']} -GL {gl_model} -out '{output_directory}/pruned_sites_for_selection' \\
        -doGlf 2 -doMajorMinor 1 -doSaf 1 -nThreads 8 \\
        -anc {reference} -sites {outputs['indexed_pos_file']}

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def selection_scan_w_PCAngsd(beagle_file: str, output_directory: str):
    """
    Doing selection scan with PCAngsd
    """
    inputs = {'GLs_beagle_file': beagle_file}
    outputs = {'selection_statistics': f'{output_directory}/pcangsd_selection_output.selection'}
    options = {
        'cores': 32,
        'memory': '64g',
        'walltime': '1:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "johannegottenborg" ]; then
        source /home/"$USER"/.bashrc
        source activate pcangsd
    fi

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    [ -d {output_directory} ] || mkdir -p {output_directory}

    pcangsd --beagle {inputs['GLs_beagle_file']} --threads 32 \\
        --out '{output_directory}/pcangsd_selection_output' --selection --maf 0

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def selection_scan_based_on_TajimasD(SAF_file: str, output_directory: str):
    """
    Calculating TajimasD based on a sliding-window algorithm implemented in angsd.
    """
    inputs = {'safs': SAF_file}
    outputs = {'sfs_file': f'{output_directory}/realSFS.sfs',
               'theta_file': f'{output_directory}/realSFS_thetas.thetas.idx',
               'TajimasD_file': f'{output_directory}/tajimasD_windows.pestPG'}
    options = {
        'cores': 32,
        'memory': '64g',
        'walltime': '1:00:00'
    }
    spec = f"""
    # Sources environment
    if [ "$USER" == "johannegottenborg" ]; then
        source /home/"$USER"/.bashrc
        source activate angsd_env
    fi

    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"

    [ -d {output_directory} ] || mkdir -p {output_directory}

    realSFS {inputs['safs']} -P 15 -fold 1 > '{output_directory}/realSFS.sfs'

    realSFS saf2theta {inputs['safs']} -sfs {outputs['sfs_file']} -outname '{output_directory}/realSFS_thetas' -fold 1

    thetaStat do_stat {outputs['theta_file']} -win 50000 -step 10000 -outnames '{output_directory}/tajimasD_windows'

    echo "END: $(date)"
    echo "$(jobinfo "$SLURM_JOBID")"
    """
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)
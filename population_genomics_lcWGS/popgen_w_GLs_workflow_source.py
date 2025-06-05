#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from popgen_w_GLs_workflow_templates import *

def running_GL_population_genomics_analyses(config_file: str = glob.glob('*config.y*ml')[0]):

	"""
	Workflow: Analyze population structure, genetic diversity and selection with 
	tools designed to work with Genotype Likelihoods (GL)

	:param str config_file:
		Configuration file containing pre-defined set of variables
	"""

	# --------------------------------------------------
	#                  Configuration                                    # Alt det nedenstående er defineret i configuration filen
	# --------------------------------------------------

	CONFIG = yaml.safe_load(open(config_file))
	ACCOUNT: str = CONFIG['account']
	TAXONOMY: str = CONFIG['taxonomic_group'].lower()
	SPECIES_NAME: str = CONFIG['species_name']
	REFERENCE_GENOME: str = CONFIG['reference_genome_path']
	WORK_DIR: str = CONFIG['working_directory_path']
	ALL_BAMS: str = CONFIG['all_202_inds']
	SUBSPECIES: str = CONFIG['blacks']
	SUBSPECIES_NAME ='blacks'
	# SUBSPECIES: str = CONFIG['blues']
	# SUBSPECIES_NAME ='blues'
	# BLUES: str = CONFIG['blues']
	GL_MODEL: list = CONFIG['GL_model']
	BASE_QUAL: int = CONFIG['minBaseQ']
	MAP_QUAL: int = CONFIG['minMapQ']
	SNP_THRESHOLD: float = CONFIG['SNP_Pval']
	MIN_DEPTH: int = CONFIG['minTotalDepth']
	MAX_DEPTH: int = CONFIG['maxTotalDepth']
	MAF: float = CONFIG['MAF']
	MIN_IND: int = CONFIG['minInd']
	HWE_PVAL: float = CONFIG['HWE_pval']


	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------

	gwf = Workflow(defaults={'account': ACCOUNT})

	for GL in GL_MODEL:
		MODEL_NAME = f'GL{GL}'

		top_dir = f'{WORK_DIR}/{SUBSPECIES_NAME}/{MODEL_NAME}' ### Remember to change this

	call_GLs_for_population_structure = gwf.target_from_template(
		name = "GLs_for_pop_structure",
		template=GLs_for_population_structure(
			bam_paths=ALL_BAMS,
			GL_model=GL_MODEL,
			output_directory=f'{top_dir}/all_inds_population_structure',
			snp_pval=SNP_THRESHOLD,
			minq=BASE_QUAL,
			minmapq=MAP_QUAL
		)
	)

	population_structure_w_PCAngsd = gwf.target_from_template(
		name = "PCAngsd_pop_structure",
		template=PCAngsd_population_structure(
			GL_file=call_GLs_for_population_structure.outputs['GLs'],
			output_directory=f'{top_dir}/all_inds_population_structure'
		)
	)

		nucleotide_diversity = gwf.target_from_template(
			name = f'calculating_pi_{SUBSPECIES_NAME}_{MODEL_NAME}',
			template=nucleotide_diversity_pi(
				bam_paths=SUBSPECIES,
				output_directory=f'{top_dir}/SAFs_and_pi',
				reference=REFERENCE_GENOME,
				gl_model=GL,
				minq=BASE_QUAL,
				minmapq=MAP_QUAL,
				minind=MIN_IND,
				maxtotaldepth=MAX_DEPTH,
				mintotaldepth=MIN_DEPTH,
				minhwepval=HWE_PVAL
			)
		)

		calling_and_filtering_GLs_for_LD_pruning = gwf.target_from_template(
			name = f'calling_and_filtering_GLs_{SUBSPECIES_NAME}_{MODEL_NAME}',
			template=GLs_filtering_and_calling(
				bam_paths=SUBSPECIES,
				output_directory=f'{top_dir}/angsd_filtering',
				gl_model=GL, 
				minq=BASE_QUAL,
				minmapq=MAP_QUAL,
				minind=MIN_IND,
				maxtotaldepth=MAX_DEPTH,
				mintotaldepth=MIN_DEPTH,
				minhwepval=HWE_PVAL,
				minmaf=MAF,
				snp_pval=SNP_THRESHOLD
			)
		)

		calculating_linkage_w_ngsLD = gwf.target_from_template(
			name = f'calculating_linkage_{SUBSPECIES_NAME}_{MODEL_NAME}',
			template=calculating_linkage_ngsLD(
				beagle_file=calling_and_filtering_GLs_for_LD_pruning.outputs['beagle'],
				bam_paths=SUBSPECIES,
				output_directory=f'{top_dir}/LD_pruning'
			)
		)

		LD_pruning = gwf.target_from_template(
			name = f'pruning_sites_of_linkage_{SUBSPECIES_NAME}_{MODEL_NAME}',
			template=LD_pruning_w_prunegraph(
				linkage_estimation=calculating_linkage_w_ngsLD.outputs['linkage_file'],
				output_directory=f'{top_dir}/LD_pruning',
			)
		)

		calling_LD_pruned_SAFs_and_GLs = gwf.target_from_template(
			name = f'calling_GLs_and_SAFs_after_pruning_{SUBSPECIES_NAME}_{MODEL_NAME}',
			template=calling_pruned_SAFs_and_GLs_w_angsd(
				bam_paths=SUBSPECIES,
				pruned_positions=LD_pruning.outputs['unlinked_positions'],
				output_directory=f'{top_dir}/LD_pruned_SAFs_and_GLs',
				reference=REFERENCE_GENOME,
				gl_model=GL
			)
		)

		PC_based_selection_scan_w_PCAngsd = gwf.target_from_template(
			name=f'PCAngsd_selection_scan_{SUBSPECIES_NAME}_{MODEL_NAME}',
			template=selection_scan_w_PCAngsd(
				beagle_file=calling_LD_pruned_SAFs_and_GLs.outputs['beagle_file'],
				output_directory=f'{top_dir}/PCAngsd_selection'
			)
		)

		calculating_TajimasD_w_angsd = gwf.target_from_template(
			name=f'TajimasD_{SUBSPECIES_NAME}_{MODEL_NAME}',
			template=selection_scan_based_on_TajimasD(
				SAF_file=calling_LD_pruned_SAFs_and_GLs.outputs['saf_file'],
				output_directory=f'{top_dir}/angsd_TajimasD'
			)
		)


	return gwf

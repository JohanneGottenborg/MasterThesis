#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from mapping_filtering_workflow_templates import *

#----------- Requirements ------------
# A conda environment containing the following software: samtools, bwa, fastp, adapterremoval


def filtering_raw_reads_mapping_filtering_bams(config_file: str = glob.glob('*config.y*ml')[0]): # Her findes configuration filen
	"""
	Workflow: Align resequencing data to reference genome and do basic filtering.
	
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
	WORK_DIR: str = CONFIG['working_directory_path'][:len(CONFIG['working_directory_path']) - 1] if CONFIG['working_directory_path'].endswith('/') else CONFIG['working_directory_path']
	SAMPLE_LISTS: list = CONFIG['sample_lists']
	ADAPTERREMOVAL_SETTINGS: dict = CONFIG['adapterremoval_settings']
	AR_MIN_QUAL: int | None = ADAPTERREMOVAL_SETTINGS['min_quality'] if ADAPTERREMOVAL_SETTINGS['min_quality'] else 25
	AR_MIN_LENGTH: int | None = ADAPTERREMOVAL_SETTINGS['min_length'] if ADAPTERREMOVAL_SETTINGS['min_length'] else 20
	AR_SEQUENCE1: str | None = ADAPTERREMOVAL_SETTINGS['adaptersequence1'] if ADAPTERREMOVAL_SETTINGS['adaptersequence1'] else 'AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA'
	AR_SEQUENCE2: str | None = ADAPTERREMOVAL_SETTINGS['adaptersequence2'] if ADAPTERREMOVAL_SETTINGS['adaptersequence2'] else 'AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG'
	SAMTOOLS_FILTER_SETTINGS: dict = CONFIG['filter_settings']
	STF_EXCLUDE: int | None = SAMTOOLS_FILTER_SETTINGS['flags_excluded'] if SAMTOOLS_FILTER_SETTINGS['flags_excluded'] else 3844
	STF_REQUIRED: int | None = SAMTOOLS_FILTER_SETTINGS['flags_required'] if SAMTOOLS_FILTER_SETTINGS['flags_required'] else None
	STF_MIN_MQ: int | None = SAMTOOLS_FILTER_SETTINGS['min_mq'] if SAMTOOLS_FILTER_SETTINGS['min_mq'] else 20
	REPEAT_REGIONS: list = CONFIG['files_containing_repeat_regions']


	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(defaults={'account': ACCOUNT})
	
	top_dir = f'{WORK_DIR}/{SPECIES_NAME.replace(" ", "_")}'

	for GROUP in SAMPLE_LISTS:
		if not GROUP['sample_folder_list']:
			continue
		GROUP_NAME: str = GROUP['group_name'].lower()
		SAMPLE_FOLDER_LIST: list = GROUP['sample_folder_list']
		sample_list = [get_sample_data(path, data_type=2) for path in SAMPLE_FOLDER_LIST] 
        ### get_sample_data is a template from workflow_templates - takes the raw reads the way they are organized in folders
        # containing sequencing data in both directions and rearranges them into a dictionary.
		# within_group_qualimap = [] # a list used for QC later on??
		for sample in sample_list:	
            ### Adapterremoval template tager output dictionary fra "get_sample_data"-template
			fastp_adapterremoval_and_trimming = gwf.target_from_template(
				name=f'T_{sample["sample_name"].replace("-", "_")}_fastp_trimming',
				template=fastp_filtering_and_adapterremoval_pairedend(
					sample_name=sample['sample_name'],
					read1_files=sample['read1_files'],
					read2_files=sample['read2_files'],
					output_directory=f'{top_dir}/{GROUP_NAME}',
					min_quality=AR_MIN_QUAL,
					min_length=AR_MIN_LENGTH,
					adapter1=AR_SEQUENCE1,
					adapter2=AR_SEQUENCE2
				)
			)

			filtering_w_adapterremoval = gwf.target_from_template(
				name=f'T_{sample["sample_name"].replace("-", "_")}_AR_filtering',
				template=adapterremoval_filtering(
					sample_name=sample['sample_name'],
					read1_file=fastp_adapterremoval_and_trimming.outputs['pair1'],
					read2_file=fastp_adapterremoval_and_trimming.outputs['pair2'],
					output_directory=f'{top_dir}/{GROUP_NAME}'
				)
			)

			# fastqc_filtered_reads = gwf.target_from_template(
			# 	name=f'T_{sample["sample_name"].replace("-", "_")}_fastqc',
			# 	template=quality_control_w_fastqc(
			# 		sample_name=sample['sample_name'],
			# 		read1_file=filtering_w_adapterremoval.outputs['pair1'],
			# 		read2_file=filtering_w_adapterremoval.outputs['pair2'],
			# 		collapsed_file=filtering_w_adapterremoval.outputs['collapsed'][0],
			# 		collapsed_file_truncated=filtering_w_adapterremoval.outputs['collapsed'][1],
			# 		output_directory=f'{top_dir}/{GROUP_NAME}'
			# 	)
			# )

			align_paired = gwf.target_from_template(
				name=f'T_{sample["sample_name"].replace("-", "_")}_paired_alignment',
				template=alignment_pairedend(
					read1_file=filtering_w_adapterremoval.outputs['pair1'],
					read2_file=filtering_w_adapterremoval.outputs['pair2'],
					reference_genome_file=REFERENCE_GENOME,
					sample_name=sample['sample_name'],
					output_directory=f'{top_dir}/{GROUP_NAME}'
				)
			)

			align_collapsed = gwf.target_from_template(
				name=f'T_{sample["sample_name"].replace("-", "_")}_collapsed_alignment',
				template=alignment_collapsed(
					collapsed_read_files=filtering_w_adapterremoval.outputs['collapsed'],
					reference_genome_file=REFERENCE_GENOME,
					sample_name=sample['sample_name'],
					output_directory=f'{top_dir}/{GROUP_NAME}'
				)
			)

			merge = gwf.target_from_template(
				name=f'T_{sample["sample_name"].replace("-", "_")}_merge_alignments',
				template=merge_alignments(
					alignment_files=[align_paired.outputs['alignment'],
									 align_collapsed.outputs['alignment']],
					sample_name=sample['sample_name'],
					output_directory=f'{top_dir}/{GROUP_NAME}'
				)
			)


			mark_duplicates = gwf.target_from_template(
				name=f'T_{sample["sample_name"].replace("-", "_")}_markdup',
				template=mark_duplicates_samtools(
					alignment_file=merge.outputs['merged'],
					sample_name=sample['sample_name'],
					output_directory=f'{top_dir}/{GROUP_NAME}'
				)
			)

			unmapped = gwf.target_from_template(
				name=f'T_{sample["sample_name"].replace("-", "_")}_extract_unmapped_reads',
				template=extract_unmapped_reads(
					alignment_file=mark_duplicates.outputs['markdup'],
					sample_name=sample['sample_name'],
					output_directory=f'{top_dir}/{GROUP_NAME}/unmapped'
				)
			)

			# qualimap_bam_before_filtering = gwf.target_from_template(
			# 	name=f'T_{sample["sample_name"].replace("-", "_")}_qualimap_before_filtering',
			# 	template=qc_qualimap(
			# 		alignment_file=mark_duplicates.outputs['markdup'],
			# 		output_directory=f'{top_dir}/{GROUP_NAME}/qualimap_bam_before_filtering/{sample["sample_name"]}'
			# 	)
			# )

			filter_alignment = gwf.target_from_template(
				name=f'T_{sample["sample_name"].replace("-", "_")}_filtering',
				template=samtools_filter(
					alignment_file=mark_duplicates.outputs['markdup'],
					sample_name=sample['sample_name'],
					output_directory=f'{top_dir}/{GROUP_NAME}/filtered_alignment',
					#output_directory=f'{top_dir}/{GROUP_NAME}/filtered_alignment_mq50',
					flags_excluded=STF_EXCLUDE,
					flags_required=STF_REQUIRED,
					min_mq=STF_MIN_MQ
				)
			)

			repeat_filtering_alignment = gwf.target_from_template(
				name=f'T_{sample["sample_name"].replace("-", "_")}_repeat_filtering',
				template=repeat_filtering(
					alignment_file=filter_alignment.outputs['filtered'],
					sample_name=sample['sample_name'],
					output_directory=f'{top_dir}/{GROUP_NAME}/repeat_filtered_alignment',
					repbase_repeats=REPEAT_REGIONS[0],
					repeatmodeler_repeats=REPEAT_REGIONS[1],
					genmap_repetitive_regions=REPEAT_REGIONS[2]
				)
			)

			qualimap_after_filtering = gwf.target_from_template(
				name=f'T_{sample["sample_name"].replace("-", "_")}_qualimap_after_filtering',
				template=qc_qualimap(
					alignment_file=repeat_filtering_alignment.outputs['filtered'],
					output_directory=f'{top_dir}/{GROUP_NAME}/qualimap_after_all_filtering/{sample["sample_name"]}'
					#output_directory=f'{top_dir}/{GROUP_NAME}/qualimap_bam_after_filtering_mq50/{sample["sample_name"]}'
				)
			)

	return gwf
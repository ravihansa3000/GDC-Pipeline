import logging
import json
import os

from swag.core.readgroups import make_readgroup_dict

import apps

LOGGER = logging.getLogger(__name__)


class GDCPatientDNASeq:
    """A patient record with normal and tumor DNA sequence.

    Class Attributes:
        gdc_output_dir: Output directory where resulting files of the pipeline will be created.
        gdc_executables: Dictionary of executables (name and path) required to execute pipeline.
        gdc_reference: Human reference genome data
        gdc_known_sites: Human variations data
        gdc_params: GDC Workflow parameters

    Instance Attributes:
        patient: Patient ID in GDC database
        patient_workdir: Directory used to create patient specific data
        bams: Dictionary of tumor and normal BAM files for a specific patient

    """
    gdc_output_dir = None
    gdc_executables = {}
    gdc_data_files = {}
    gdc_params = {}

    def __init__(self, patient, bams, label=None):
        LOGGER.debug(f"GDCPatientDNASeq patient: {patient} | bams: {str(bams)}")
        if label is None or label == '':
            self.patient = patient
            self.patient_workdir = os.path.join(GDCPatientDNASeq.gdc_output_dir, patient)
        else:
            LOGGER.debug(f"GDCPatientDNASeq label: {label}")
            self.patient = f"patient_{label}"
            self.patient_workdir = os.path.join(GDCPatientDNASeq.gdc_output_dir, patient, label)

        self.bams = bams

    def process_patient_seq_data(self):
        if not os.path.exists(self.patient_workdir):
            os.makedirs(self.patient_workdir)
        LOGGER.debug(f"patient_workdir created: {self.patient_workdir}")

        # List of BAM files that belong to the same read group that will be merged into a single BAM file
        # and then cleaned (indel-realined with BQSR engine) via GDC co-cleaning workflow
        cleaned_bams = {}

        for tissue in self.bams:
            tissue_bam = self.bams[tissue]
            _, tail = os.path.split(tissue_bam)
            bam_workdir = os.path.join(self.patient_workdir, tail.replace('.bam', ''))
            if not os.path.exists(bam_workdir):
                os.makedirs(bam_workdir)
            LOGGER.debug(f"bam_workdir created | bam_workdir: {bam_workdir}, tissue_bam: {tissue_bam}, \
                         tail: {tail}")

            # Prior to alignment, BAM files that were submitted to the GDC are split
            # by read groups and converted to FASTQ format.
            bamtofastq_output = self.do_pre_alignment(tissue_bam, bam_workdir)

            # During sort and alignment BAM file is split by read group then
            # each read group is aligned to the reference genome separately
            aligned_bams = self.align_and_sort(tissue_bam, bam_workdir, bamtofastq_output)

            # All read group alignments that belong to a single aliquot are merged
            merge_output_dir = os.path.join(bam_workdir, 'merged')
            merged_bam_file = os.path.join(merge_output_dir, tail.replace('.bam', '.merged.bam'))
            LOGGER.info(
                f"Running merge_and_mark_duplicates: merged_bam_file: {merged_bam_file}, tissue_bam: {tissue_bam}")

            merge_mark_output = apps.merge_and_mark_duplicates(
                GDCPatientDNASeq.gdc_executables,
                os.path.join(merge_output_dir, 'tmp'),
                inputs=aligned_bams,
                outputs=[merged_bam_file],
                label=self.patient
            ).outputs[0]
            LOGGER.debug(f"merged_bam: patient: {self.patient} | {merge_mark_output}")

            co_clean_output_dir = os.path.join(bam_workdir, 'co_clean')
            cleaned_bam_file = os.path.join(co_clean_output_dir, tail.replace('.bam', '.cleaned.bam'))
            cleaned_bams[tissue] = apps.co_cleaning_pipeline(
                GDCPatientDNASeq.gdc_executables,
                GDCPatientDNASeq.gdc_data_files['gdc_reference_seq_fa'],
                GDCPatientDNASeq.gdc_data_files['dbsnp_known_snp_sites'],
                GDCPatientDNASeq.gdc_data_files['known_indels'],
                input_bam=merge_mark_output,
                output_dir=co_clean_output_dir,
                outputs=[cleaned_bam_file],
                label=self.patient
            ).outputs[0]
            LOGGER.debug(f"co_cleaning_pipeline: patient: {self.patient} | {cleaned_bams[tissue]}")

        return cleaned_bams

    def do_pre_alignment(self, tissue_bam, bam_workdir):
        _, tail = os.path.split(tissue_bam)
        sorted_bam_output_dir = os.path.join(bam_workdir, 'sorted')
        if not os.path.exists(sorted_bam_output_dir):
            os.makedirs(sorted_bam_output_dir)
        sorted_bam_filepath = os.path.join(sorted_bam_output_dir, tail.replace('.bam', 'sorted.bam'))
        LOGGER.info(f"Running sort_bam_by_queryname: tissue_bam: {tissue_bam}, \
                    sorted_bam_filepath: {sorted_bam_filepath}, bam_workdir: {bam_workdir}")

        sorted_bam_output = apps.sort_bam_by_queryname(
            GDCPatientDNASeq.gdc_executables,
            os.path.join(sorted_bam_output_dir, 'tmp'),
            bam_filepath=tissue_bam,
            outputs=[sorted_bam_filepath],
            label=self.patient
        ).outputs[0]

        fastq_output_dir = os.path.join(bam_workdir, 'fastq')
        LOGGER.info(f"Running bamtofastq: sorted_bam_filepath: {sorted_bam_filepath}, \
                    fastq_output_dir: {fastq_output_dir}")

        bamtofastq_output = apps.bamtofastq(
            GDCPatientDNASeq.gdc_executables,
            sorted_bam_output,
            outputs=[fastq_output_dir],
            label=self.patient
        ).outputs[0]

        LOGGER.debug(f"bamtofastq_output: tissue_bam: {tissue_bam} | {str(bamtofastq_output)}")
        return bamtofastq_output

    def align_and_sort(self, tissue_bam, bam_workdir, bamtofastq_output):
        _, tail = os.path.split(tissue_bam)
        sort_and_align_output_dir = os.path.join(bam_workdir, 'align_and_sort')
        if not os.path.exists(sort_and_align_output_dir):
            os.makedirs(sort_and_align_output_dir)

        # list of align_and_sort App Futures (BAM files by read group) that will be merged back together
        aligned_bams = []

        read_groups = make_readgroup_dict(tissue_bam, GDCPatientDNASeq.gdc_executables['samtools'])
        LOGGER.debug(
            f"patient: {self.patient}, sort_and_align_output_dir: {sort_and_align_output_dir}, \
            make_readgroup_dict: {json.dumps(read_groups, indent=4)}")

        for rg_id, rg_line in read_groups.items():
            align_and_sort_output_filepath = os.path.join(sort_and_align_output_dir, tail.replace(
                '.bam', '{}.aligned.sorted.bam'.format(rg_id)))

            LOGGER.info(
                f"Running align_and_sort: align_and_sort_output_filepath: {align_and_sort_output_filepath}, \
                patient: {self.patient}, rg_id: {rg_id}, rg_line: {rg_line}")

            align_and_sort_output = apps.align_and_sort(
                GDCPatientDNASeq.gdc_executables,
                os.path.join(sort_and_align_output_dir, 'tmp'),
                GDCPatientDNASeq.gdc_data_files['gdc_reference_seq_fa'],
                bamtofastq_output,
                rg_id,
                rg_line.replace('\t', '\\t'),
                outputs=[align_and_sort_output_filepath],
                label=self.patient,
            ).outputs[0]

            aligned_bams.append(align_and_sort_output)
            LOGGER.debug(f"align_and_sort_output: {align_and_sort_output}")

        LOGGER.debug(
            f"align_and_sort: patient: {self.patient} | aligned_bams: {str(aligned_bams)}")
        return aligned_bams

    def validate_analysis_input(self, cleaned_bams, has_tumor=True):
        if not 'normal' in cleaned_bams:
            raise RuntimeError("Invalid input format | patient {self.patient}: missing normal BAM file")

        if has_tumor and not 'tumor' in cleaned_bams:
            return False

        return True

    def run_annotation_workflow(self, input_vcf_future, input_vcf_list):
        annotation_out_dir = os.path.join(self.patient_workdir, 'annotations')
        if not os.path.exists(annotation_out_dir):
            os.makedirs(annotation_out_dir)

        for input_vcf_gz in input_vcf_list:
            _, input_vcf_gz_filename = os.path.split(input_vcf_gz)
            vcf_filename = input_vcf_gz_filename.replace('.gz', '')
            output_file = f"{annotation_out_dir}/{vcf_filename}"
            apps.annotatate_vep(
                GDCPatientDNASeq.gdc_executables,
                GDCPatientDNASeq.gdc_data_files['gdc_vep_cache_dir'],
                input_vcf_future,
                vcf_filename,
                annotation_out_dir,
                outputs=[output_file],
                label=f"{self.patient}-{vcf_filename}-annotations"
            )

    def run_variant_callers(self, cleaned_bams):
        if (GDCPatientDNASeq.gdc_params.get('gdc_somaticsniper_enabled', False) and
                self.validate_analysis_input(cleaned_bams)):
            LOGGER.info(f"Running somaticsniper: patient: {self.patient}")
            apps.somaticsniper(
                GDCPatientDNASeq.gdc_executables,
                GDCPatientDNASeq.gdc_data_files['gdc_reference_seq_fa'],
                cleaned_bams['normal'],
                cleaned_bams['tumor'],
                self.patient_workdir,
                label='{}-somaticsniper'.format(self.patient)
            )

        if (GDCPatientDNASeq.gdc_params.get('gdc_muse_enabled', False) and
                self.validate_analysis_input(cleaned_bams)):
            LOGGER.info(f"Running muse: patient: {self.patient}")
            apps.muse(
                GDCPatientDNASeq.gdc_executables,
                GDCPatientDNASeq.gdc_data_files['gdc_reference_seq_fa'],
                cleaned_bams['normal'],
                cleaned_bams['tumor'],
                GDCPatientDNASeq.gdc_data_files['dbsnp_known_snp_sites'],
                self.patient_workdir,
                label='{}-muse'.format(self.patient)
            )

        if (GDCPatientDNASeq.gdc_params.get('gdc_varscan_enabled', False) and
                self.validate_analysis_input(cleaned_bams)):
            LOGGER.info(f"Running varscan: patient: {self.patient}")
            apps.varscan(
                GDCPatientDNASeq.gdc_executables,
                os.path.join(self.patient_workdir, 'tmp'),
                GDCPatientDNASeq.gdc_data_files['gdc_reference_seq_fa'],
                cleaned_bams['normal'],
                cleaned_bams['tumor'],
                self.patient_workdir,
                label='{}-varscan'.format(self.patient)
            )

        if (GDCPatientDNASeq.gdc_params.get('gdc_strelka2_somatic_enabled', False) and
                self.validate_analysis_input(cleaned_bams)):
            somatic_analysis_path = os.path.join(self.patient_workdir, 'strelka2-analysis-somatic')
            if not os.path.exists(somatic_analysis_path):
                os.makedirs(somatic_analysis_path)

            indels_output = os.path.join(somatic_analysis_path, 'results', 'variants', 'somatic.indels.vcf.gz')
            snvs_output = os.path.join(somatic_analysis_path, 'results', 'variants', 'somatic.snvs.vcf.gz')
            somatic_output_list = [indels_output, snvs_output]
            LOGGER.info(
                f"Running strelka2 somatic analysis: patient: {self.patient}, \
                analysis_output: {somatic_analysis_path}")

            strelka2_somatic_vcf_future = apps.strelka2_somatic(
                GDCPatientDNASeq.gdc_executables,
                GDCPatientDNASeq.gdc_data_files['gdc_reference_seq_fa'],
                cleaned_bams['normal'],
                cleaned_bams['tumor'],
                somatic_analysis_path,
                outputs=somatic_output_list,
                label='{}-strelka2-somatic'.format(self.patient)
            ).outputs[0]
            LOGGER.debug(
                f"""strelka2_somatic: patient: {self.patient} |
                strelka2_somatic_vcf_future: {str(strelka2_somatic_vcf_future)}""")
            self.run_annotation_workflow(strelka2_somatic_vcf_future, somatic_output_list)

        if (GDCPatientDNASeq.gdc_params.get('gdc_strelka2_germline_enabled', False) and
                self.validate_analysis_input(cleaned_bams, False)):
            germline_analysis_path = os.path.join(self.patient_workdir, 'strelka2-analysis-germline')
            if not os.path.exists(germline_analysis_path):
                os.makedirs(germline_analysis_path)
            variant_output = os.path.join(germline_analysis_path, 'results', 'variants', 'variants.vcf.gz')
            germline_output_list = [variant_output]

            LOGGER.info(
                f"""Running strelka2 germline analysis: patient: {self.patient},
                analysis_output: {germline_analysis_path}""")

            strelka2_germline_vcf_future = apps.strelka2_germline(
                GDCPatientDNASeq.gdc_executables,
                GDCPatientDNASeq.gdc_data_files['gdc_reference_seq_fa'],
                cleaned_bams['normal'],
                germline_analysis_path,
                outputs=germline_output_list,
                label='{}-strelka2-germline'.format(self.patient)
            ).outputs[0]

            LOGGER.debug(
                f"""strelka2_germline: patient: {self.patient} |
                strelka2_germline_vcf_future: {str(strelka2_germline_vcf_future)}""")
            self.run_annotation_workflow(strelka2_germline_vcf_future, germline_output_list)

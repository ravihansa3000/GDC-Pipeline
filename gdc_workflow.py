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
        gdc_tmp_dir: Temporary directory used to store temporary files created
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
    gdc_tmp_dir = None
    gdc_executables = {}
    gdc_reference = None
    gdc_known_sites = None
    gdc_params = {}

    def __init__(self, patient, bams):
        self.patient = patient
        self.patient_workdir = os.path.join(GDCPatientDNASeq.gdc_output_dir, patient)
        self.bams = bams

    def process_patient_seq_data(self):
        if not os.path.exists(self.patient_workdir):
            os.makedirs(self.patient_workdir)
        LOGGER.debug(f"patient_workdir created: {self.patient_workdir}")

        # List of BAM files that belong to the same read group that will be merged into a single BAM file
        merged_bams = {}

        for tissue in ['tumor', 'normal']:
            tissue_bam = self.bams[tissue]
            head, tail = os.path.split(tissue_bam)
            bam_workdir = os.path.join(self.patient_workdir, tail.replace('.bam', ''))
            if not os.path.exists(bam_workdir):
                os.makedirs(bam_workdir)
            LOGGER.debug(f"bam_workdir created | bam_workdir: {bam_workdir}, tissue_bam: {tissue_bam},"
                         " head: {head}, tail: {tail}")

            # Prior to alignment, BAM files that were submitted to the GDC are split
            # by read groups and converted to FASTQ format.
            sorted_bam_filepath = os.path.join(bam_workdir, tail.replace('.bam', '.sorted.bam'))
            bamtofastq_output = self.do_pre_alignment(tissue_bam, sorted_bam_filepath, bam_workdir)

            # During sort and alignment BAM file is split by read group then
            # each read group is aligned to the reference genome separately
            aligned_bams = self.sort_and_align(tissue_bam, bamtofastq_output)

            # All read group alignments that belong to a single aliquot are merged
            merged_bam_file = os.path.join(bam_workdir, tail.replace('.bam', '.merged.bam'))
            LOGGER.info(
                f"Running merge_and_mark_duplicates: merged_bam_file: {merged_bam_file}, tissue_bam: {tissue_bam}")

            merged_bams[tissue] = apps.merge_and_mark_duplicates(
                GDCPatientDNASeq.gdc_executables,
                inputs=aligned_bams,
                outputs=[merged_bam_file],
                label=self.patient
            ).outputs[0]

        LOGGER.debug(f"merged_bams: patient: {self.patient} | {str(merged_bams)}")
        return merged_bams

    def do_pre_alignment(self, tissue_bam, sorted_bam_filepath, bam_workdir):
        LOGGER.info(f"Running sort_bam_by_queryname: tissue_bam: {tissue_bam}, "
                    "sorted_bam_filepath: {sorted_bam_filepath}, bam_workdir: {bam_workdir}")

        sorted_bam_output = apps.sort_bam_by_queryname(
            GDCPatientDNASeq.gdc_executables,
            GDCPatientDNASeq.gdc_tmp_dir,
            bam_filepath=tissue_bam,
            outputs=[sorted_bam_filepath],
            label=self.patient
        ).outputs[0]

        fastq_output_dir = os.path.join(bam_workdir, 'fastq')
        LOGGER.info(f"Running bamtofastq: sorted_bam_filepath: {sorted_bam_filepath},"
                    "fastq_output_dir: {fastq_output_dir}")

        bamtofastq_output = apps.bamtofastq(
            GDCPatientDNASeq.gdc_executables,
            sorted_bam_output,
            outputs=[fastq_output_dir],
            label=self.patient
        ).outputs[0]

        LOGGER.debug(f"bamtofastq_output: tissue_bam: {tissue_bam} | {str(bamtofastq_output)}")
        return bamtofastq_output

    def sort_and_align(self, tissue_bam, bamtofastq_output):
        head, tail = os.path.split(tissue_bam)
        bam_workdir = os.path.join(self.patient_workdir, tail.replace('.bam', ''))
        LOGGER.debug("sort_and_align created bam_workdir: {}")

        # list of align_and_sort App Futures (BAM files by read group) that will be merged back together
        aligned_bams = []

        read_groups = make_readgroup_dict(tissue_bam, GDCPatientDNASeq.gdc_executables['samtools'])
        LOGGER.debug(f"patient: {self.patient}, make_readgroup_dict: {json.dumps(read_groups, indent=4)}")

        for rg_id, rg_line in read_groups.items():
            align_and_sort_output_dir = os.path.join(bam_workdir, tail.replace(
                '.bam', '{}.sorted.aligned.bam'.format(rg_id)))

            LOGGER.info(
                f"Running align_and_sort: align_and_sort_output_dir: {align_and_sort_output_dir}, "
                "patient: {self.patient}, rg_id: {rg_id}, rg_line: {rg_line}")

            align_and_sort_output = apps.align_and_sort(
                GDCPatientDNASeq.gdc_executables,
                GDCPatientDNASeq.gdc_reference,
                bamtofastq_output,
                rg_id,
                rg_line.replace('\t', '\\t'),
                outputs=[align_and_sort_output_dir],
                label=self.patient,
            ).outputs[0]

            aligned_bams.append(align_and_sort_output)
            LOGGER.debug(f"align_and_sort_output: {align_and_sort_output}")

        LOGGER.debug(
            f"sort_and_align: patient: {self.patient} | aligned_bams: {str(aligned_bams)}")
        return aligned_bams

    def run_somatic_variant_callers(self, merged_bams):
        if GDCPatientDNASeq.gdc_params.get('gdc_somaticsniper_enabled', False):
            LOGGER.info(f"Running somaticsniper: patient: {self.patient}")
            apps.somaticsniper(
                GDCPatientDNASeq.gdc_executables,
                GDCPatientDNASeq.gdc_reference,
                merged_bams['normal'],
                merged_bams['tumor'],
                self.patient_workdir,
                label='{}-somaticsniper'.format(self.patient)
            )

        if GDCPatientDNASeq.gdc_params.get('gdc_muse_enabled', False):
            LOGGER.info(f"Running muse: patient: {self.patient}")
            apps.muse(
                GDCPatientDNASeq.gdc_executables,
                GDCPatientDNASeq.gdc_reference,
                merged_bams['normal'],
                merged_bams['tumor'],
                GDCPatientDNASeq.gdc_known_sites,
                self.patient_workdir,
                label='{}-muse'.format(self.patient)
            )

        if GDCPatientDNASeq.gdc_params.get('gdc_muse_enabled', False):
            LOGGER.info(f"Running varscan: patient: {self.patient}")
            apps.varscan(
                GDCPatientDNASeq.gdc_executables,
                GDCPatientDNASeq.gdc_reference,
                merged_bams['normal'],
                merged_bams['tumor'],
                self.patient_workdir,
                label='{}-varscan'.format(self.patient)
            )

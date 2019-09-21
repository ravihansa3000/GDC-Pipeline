import json
import os
import logging
import sys
import parsl
from parsl.utils import get_all_checkpoints

from gdc_workflow import GDCPatientDNASeq
import gdc_workflow

from config.parsl_config import config

LOGGER = logging.getLogger(__name__)


def setup_gdc_pipeline(params):
    gdc_output_dir = params['gdc_output_dir']
    gdc_run_dir = params['gdc_run_dir']
    gdc_tmp_dir = params['gdc_tmp_dir']

    if not os.path.exists(gdc_output_dir):
        os.makedirs(gdc_output_dir)

    if not os.path.exists(gdc_run_dir):
        os.makedirs(gdc_run_dir)
    config.run_dir = gdc_run_dir

    if not os.path.exists(gdc_tmp_dir):
        os.makedirs(gdc_tmp_dir)

    # Parsl checkpoint is created each time an app completes or fails
    config.checkpoint_files = get_all_checkpoints()
    config.checkpoint_mode = 'task_exit'


def run_gdc_pipeline(params):
    gdc_bam_files = params['gdc_bam_files']
    parsl.set_stream_logger()
    parsl.load(config)
    LOGGER.info("GDC Pipeline started!")

    GDCPatientDNASeq.gdc_output_dir = params['gdc_output_dir']
    GDCPatientDNASeq.gdc_tmp_dir = params['gdc_tmp_dir']
    GDCPatientDNASeq.gdc_executables = params['gdc_executables']
    GDCPatientDNASeq.gdc_data_files = params['gdc_data_files']
    GDCPatientDNASeq.gdc_params = params
    gdc_workflow.LOGGER = LOGGER

    for patient, bams in gdc_bam_files.items():
        gdc_patient = GDCPatientDNASeq(patient, bams)
        merged_bams = gdc_patient.process_patient_seq_data()
        gdc_patient.run_variant_callers(merged_bams)

    LOGGER.info("Waiting for GDC Pipeline tasks to complete...")
    parsl.wait_for_current_tasks()
    LOGGER.info("GDC Pipeline tasks completed!")


def validate_config(params):
    missing_bam_files = []
    gdc_bam_files = params['gdc_bam_files']
    for patient, bams in gdc_bam_files.items():
        for tissue in ['tumor', 'normal']:
            bam_file = bams[tissue]
            if not os.path.exists(bam_file):
                LOGGER.error(f"BAM file: {bam_file} not found for patient: {patient}")
                missing_bam_files.append(bam_file)
    if missing_bam_files:
        LOGGER.error(f"Missing BAM files: {str(missing_bam_files)}")
        raise RuntimeError("GDC Pipeline validation failed.")

    missing_executables = []
    gdc_executables = params['gdc_executables']
    for key, val in gdc_executables.items():
        if not os.path.exists(val):
            LOGGER.error(f"Executable: {val} not found for key: {key}")
            missing_executables.append(val)
    if missing_executables:
        LOGGER.error(f"Missing dependencies: {str(missing_executables)}")
        raise RuntimeError("GDC Pipeline validation failed.")

    missing_data_files = []

    gdc_data_files = params['gdc_data_files']
    for key, val in gdc_data_files.items():
        if not os.path.exists(val):
            LOGGER.error(f"Data File: {val} not found for key: {key}")
            missing_data_files.append(val)
    if missing_data_files:
        LOGGER.error(f"Missing data files: {str(missing_executables)}")
        raise RuntimeError("GDC Pipeline validation failed.")


def load_defaults(params, key, val):
    if (key not in params or params[key] is None or params[key] == ''):
        params[key] = val


def load_configs():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(dir_path, 'config', 'gdc_config.json')) as f:
        gdc_config = json.load(f)

    print(f" ======== Running GDC Pipeline ========")
    print(json.dumps(gdc_config, indent=4))

    load_defaults(gdc_config, 'gdc_output_dir',
                  os.path.join(dir_path, 'output'))

    gdc_output_dir = gdc_config['gdc_output_dir']
    load_defaults(gdc_config, 'gdc_run_dir',
                  os.path.join(gdc_output_dir, 'runinfo'))
    load_defaults(gdc_config, 'gdc_bam_files.json',
                  os.path.join(dir_path, 'documents', 'data.json'))
    load_defaults(gdc_config, 'gdc_executables.json',
                  os.path.join(dir_path, 'documents', 'executables.json'))
    load_defaults(gdc_config, 'gdc_tmp_dir',
                  os.path.join(gdc_output_dir, 'tmp'))

    with open(gdc_config['gdc_bam_files.json']) as f:
        gdc_config['gdc_bam_files'] = json.load(f)

    with open(gdc_config['gdc_executables.json']) as f:
        gdc_executables = json.load(f)

    for key, val in gdc_executables.items():
        gdc_executables[key] = os.path.join(gdc_config['gdc_executables_dir'], val)

    if gdc_config.get('gdc_strelka2_somatic_enabled', False) or gdc_config.get('gdc_strelka2_germline_enabled', False):
        strelka2_install_path = gdc_config.get('STRELKA_INSTALL_PATH', gdc_config['gdc_executables_dir'])
        gdc_executables['strelka2_somatic_configure'] = os.path.join(
            strelka2_install_path, 'bin', 'configureStrelkaSomaticWorkflow.py')
        gdc_executables['strelka2_germline_configure'] = os.path.join(
            strelka2_install_path, 'bin', 'configureStrelkaGermlineWorkflow.py')

    gdc_config['gdc_executables'] = gdc_executables

    gdc_data_files = {
        'gdc_reference_seq_fa': gdc_config['gdc_reference_seq_fa'],
        'dbsnp_known_snp_sites': gdc_config['dbsnp_known_snp_sites'],
        'dbsnp_known_snp_sites_index': "{}.tbi".format(gdc_config['dbsnp_known_snp_sites'])
    }
    gdc_config['gdc_data_files'] = gdc_data_files

    return gdc_config


def setup_logger(params, logger_obj):
    log_level = os.getenv('log_level', params.get('log_level', 'INFO')).upper()
    formatter = logging.Formatter(
        '%(asctime)s %(name)s [%(levelname)s]  %(message)s')
    logger_obj.setLevel(log_level)

    fh = logging.FileHandler(os.path.join(params['gdc_output_dir'], 'GDC_Pipeline.log'))
    fh.setLevel(log_level)
    fh.setFormatter(formatter)
    logger_obj.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level)
    ch.setFormatter(formatter)
    logger_obj.addHandler(ch)


def main(argv):
    gdc_config = load_configs()
    setup_logger(gdc_config, LOGGER)
    validate_config(gdc_config)
    setup_gdc_pipeline(gdc_config)
    run_gdc_pipeline(gdc_config)


if __name__ == "__main__":
    main(sys.argv[1:])

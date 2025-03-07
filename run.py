import json
import os
import logging
import sys

import parsl
from parsl.utils import get_all_checkpoints

import gdc_workflow
from gdc_workflow import GDCPatientDNASeq
from config.parsl_config import (
    get_parsl_config_nscc,
    get_parsl_config_local,
    get_parsl_config_csi
)

LOGGER = logging.getLogger(__name__)


def setup_gdc_pipeline(params):
    gdc_output_dir = params['gdc_output_dir']
    gdc_run_dir = params['gdc_run_dir']

    if not os.path.exists(gdc_run_dir):
        os.makedirs(gdc_run_dir)

    if params['parsl_config_env'] == 'nscc':
        parsl_config = get_parsl_config_nscc()
    elif params['parsl_config_env'] == 'local':
        parsl_config = get_parsl_config_local()
    elif params['parsl_config_env'] == 'csi':
        parsl_config = get_parsl_config_csi()

    parsl_config.run_dir = gdc_run_dir

    # Parsl checkpointing: resume using from all available checkpoints
    parsl_config.checkpoint_files = get_all_checkpoints(gdc_run_dir)

    # Setup monitoring
    if parsl_config.monitoring is not None:
        parsl_config.monitoring.logging_endpoint = "sqlite:///{}/monitoring.db".format(gdc_output_dir)

    params['parsl_config'] = parsl_config


def run_gdc_pipeline(params):
    gdc_bam_files = params['gdc_bam_files']
    parsl.set_stream_logger()
    parsl.load(params['parsl_config'])
    LOGGER.info("GDC Pipeline started!")

    GDCPatientDNASeq.gdc_output_dir = params['gdc_output_dir']
    GDCPatientDNASeq.gdc_executables = params['gdc_executables']
    GDCPatientDNASeq.gdc_data_files = params['gdc_data_files']
    GDCPatientDNASeq.gdc_params = params
    gdc_workflow.LOGGER = LOGGER

    def process_bam_pair(patient, bam_pair, label=None):
        if isinstance(bam_pair, list):
            bam_pair = bam_pair[0]

        gdc_patient = GDCPatientDNASeq(patient, bam_pair, label)
        cleaned_bam_pair = {}
        if ('cleaned' in bam_pair) and (bam_pair['cleaned']):
            cleaned_bam_pair = bam_pair
        else:
            cleaned_bam_pair = gdc_patient.process_patient_seq_data()

        gdc_patient.run_variant_callers(cleaned_bam_pair)

    for patient, bam_pair_list in gdc_bam_files.items():
        if isinstance(bam_pair_list, dict) or len(bam_pair_list) == 1:
            process_bam_pair(patient, bam_pair_list)

        else:
            count = 1
            for bam_pair in bam_pair_list:
                process_bam_pair(patient, bam_pair, str(count))
                count += 1

    LOGGER.info("Waiting for GDC Pipeline tasks to complete...")
    parsl.wait_for_current_tasks()
    LOGGER.info("GDC Pipeline tasks completed!")


def validate_config(params):
    missing_bam_files = []
    gdc_bam_files = params['gdc_bam_files']

    def validate_bam_pair(patient, bam_pair):
        if isinstance(bam_pair, list):
            bam_pair = bam_pair[0]

        if not 'normal' in bam_pair:
            raise RuntimeError(f"Patient {patient} record does not contain a normal BAM file record")

        for tissue in bam_pair:
            bam_file = bam_pair[tissue]
            if not os.path.exists(bam_file):
                LOGGER.error(f"BAM file: {bam_file} not found for patient: {patient}")
                missing_bam_files.append(bam_file)

    for patient, bam_pair_list in gdc_bam_files.items():
        if isinstance(bam_pair_list, dict) or len(bam_pair_list) == 1:
            validate_bam_pair(patient, bam_pair_list)
        else:
            for bam_pair in bam_pair_list:
                validate_bam_pair(patient, bam_pair)

    if missing_bam_files:
        LOGGER.error(f"Missing BAM files: {str(missing_bam_files)}")
        raise RuntimeError(f"GDC Pipeline validation failed | missing BAM files: {str(missing_bam_files)}")

    missing_executables = []
    gdc_executables = params['gdc_executables']
    for key, val in gdc_executables.items():
        if not os.path.exists(val):
            LOGGER.error(f"Executable: {val} not found for key: {key}")
            missing_executables.append(val)

    if missing_executables:
        LOGGER.error(f"Missing dependencies: {str(missing_executables)}")
        raise RuntimeError(f"GDC Pipeline validation failed | missing dependencies: {str(missing_executables)}")

    missing_data_files = []

    gdc_data_files = params['gdc_data_files']
    for key, val in gdc_data_files.items():
        if not os.path.exists(val):
            LOGGER.error(f"Data File: {val} not found for key: {key}")
            missing_data_files.append(val)

    if missing_data_files:
        LOGGER.error(f"Missing data files: {str(missing_data_files)}")
        raise RuntimeError(f"GDC Pipeline validation failed | missing data files: {str(missing_data_files)}")


def load_defaults(params, key, val):
    if (val is None or val == ''):
        return

    if (key not in params or params[key] is None or params[key] == ''):
        params[key] = val


def load_defaults_dir(params, key, val):
    load_defaults(params, key, val)
    params[key] = os.path.expanduser(params[key])


def load_configs():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    default_gdc_config_file = os.path.join(dir_path, 'config', 'gdc_config.json')
    gdc_config_file = os.environ.get('GDC_CONFIG_FILE', default_gdc_config_file)

    with open(gdc_config_file) as f:
        gdc_config = json.load(f)

    print(f" ======== Running GDC Pipeline ========")
    print(json.dumps(gdc_config, indent=4))

    gdc_bam_input = os.getenv('GDC_BAM_INPUT')
    load_defaults_dir(gdc_config, 'gdc_bam_files.json', gdc_bam_input)
    if gdc_config['gdc_bam_files.json'] is None or gdc_config['gdc_bam_files.json'] == '':
        raise Exception('GDC-Pipeline failed | input bam file is missing')

    load_defaults(gdc_config, 'parsl_config_env', 'local')
    _, tail = os.path.split(gdc_config['gdc_bam_files.json'])
    campaign_name = tail.replace('.json', '')
    load_defaults(gdc_config, 'campaign_name', campaign_name)

    load_defaults_dir(gdc_config, 'gdc_output_dir', os.path.join(dir_path, 'output'))
    gdc_config['gdc_output_dir'] = os.path.join(gdc_config['gdc_output_dir'], gdc_config['campaign_name'])

    load_defaults_dir(gdc_config, 'gdc_run_dir', os.path.join(gdc_config['gdc_output_dir'], 'runinfo'))
    load_defaults_dir(gdc_config, 'gdc_executables.json', os.path.join(dir_path, 'documents', 'executables.json'))
    load_defaults_dir(gdc_config, 'gdc_executables_dir', os.path.join('~', 'anaconda3', 'envs', 'gdc', 'bin'))
    load_defaults_dir(gdc_config, 'STRELKA_INSTALL_PATH', gdc_config['gdc_executables_dir'])

    with open(gdc_config['gdc_bam_files.json']) as f:
        gdc_config['gdc_bam_files'] = json.load(f)

    with open(gdc_config['gdc_executables.json']) as f:
        gdc_executables = json.load(f)

    for key, val in gdc_executables.items():
        gdc_executables[key] = os.path.join(gdc_config['gdc_executables_dir'], val)

    if gdc_config.get('gdc_strelka2_somatic_enabled', False) or gdc_config.get('gdc_strelka2_germline_enabled', False):
        strelka2_install_path = gdc_config['STRELKA_INSTALL_PATH']
        gdc_executables['strelka2_somatic_configure'] = os.path.join(
            strelka2_install_path, 'bin', 'configureStrelkaSomaticWorkflow.py')
        gdc_executables['strelka2_germline_configure'] = os.path.join(
            strelka2_install_path, 'bin', 'configureStrelkaGermlineWorkflow.py')

    gdc_config['gdc_executables'] = gdc_executables

    gdc_data_files = {
        'gdc_reference_seq_fa': gdc_config['gdc_reference_seq_fa'],
        'gdc_reference_seq_fa_fai': "{}.fai".format(gdc_config['gdc_reference_seq_fa']),
        'gdc_reference_seq_fa_dict': "{}.dict".format(gdc_config['gdc_reference_seq_fa']),
        'dbsnp_known_snp_sites': gdc_config['dbsnp_known_snp_sites'],
        'dbsnp_known_snp_sites_index': "{}.tbi".format(gdc_config['dbsnp_known_snp_sites']),
        'known_indels': gdc_config['known_indels'],
        'known_indels_index': "{}.tbi".format(gdc_config['known_indels']),
        'gdc_vep_cache_dir': gdc_config['gdc_vep_cache_dir']
    }
    gdc_config['gdc_data_files'] = gdc_data_files

    return gdc_config


def setup_logger(params, logger_obj):
    print("****** Setting up GDC Pipeline ******")
    print(json.dumps(params, indent=4))
    log_level = os.getenv('log_level', params.get('log_level', 'INFO')).upper()
    formatter = logging.Formatter(
        '%(asctime)s %(name)s [%(levelname)s]  %(message)s')
    logger_obj.setLevel(log_level)

    gdc_output_dir = params['gdc_output_dir']
    if not os.path.exists(gdc_output_dir):
        os.makedirs(gdc_output_dir)

    fh = logging.FileHandler(os.path.join(gdc_output_dir, 'GDC_Pipeline.log'))
    fh.setLevel(log_level)
    fh.setFormatter(formatter)
    logger_obj.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(log_level)
    ch.setFormatter(formatter)
    logger_obj.addHandler(ch)

    logger_obj.info("****** setup complete ******")


def main(argv):
    gdc_config = load_configs()
    setup_logger(gdc_config, LOGGER)
    validate_config(gdc_config)
    setup_gdc_pipeline(gdc_config)
    run_gdc_pipeline(gdc_config)


if __name__ == "__main__":
    main(sys.argv[1:])

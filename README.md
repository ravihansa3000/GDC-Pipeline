# GDC Pipeline Workflow based on Parsl

## QuickStart

### Create Parsl config file

Create a Python configuration file `parsl_config.py` in workspace `config` dir using one of the available templates (`parsl_config.py.template`)

### Create GDC Pipeline config file

Create a json configuration file `gdc_config.json` in workspace `config` dir using the template `gdc_config.json.template`

### Create run-time environment

- Download GDC-Pipeline resource files

    `./scripts/gdc-installer.sh`

- Generate indexes for genome reference

    `./scripts/gdc-generate-ref-index.sh`

- Install dependencies and setup conda environment

    `./scripts/gdc-installer.sh`


### Run GDC-Pipeline Workflow

- Use configuration defined in `gdc_config.json`
`python run.py`

- Override parameters defined in `gdc_config.json` using environment variables
`GDC_BAM_INPUT=<input_campaign.json> GDC_CONFIG_FILE=<gdc_config_env.json> python ./run.py`




# GDC Pipeline Workflow based on Parsl

## QuickStart

### Create Parsl config file

Create a Python configuration file `parsl_config.py` in workspace `config` dir using one of the available templates (`aspire_config.py.template`, `local_config.py.template`)

### Create GDC Pipeline config file

Create a json configuration file `gdc_config.json` in workspace `config` dir using the template `gdc_config.json.template`

### Create conda environment

```
conda create --name gdc python=3.7
source activate gdc
pip install --no-binary pyzmq pyzmq
pip install swag
```

### Install dependencies

Run `scripts/gdc-btools.sh` to install dependencies

### RUN command

`python run.py`

`GDC_BAM_INPUT=/home/akila/Documents/NUS_CSI/Parsl/workspace/TCGA_OV_1.json GDC_CONFIG_FILE=/home/akila/Documents/NUS_CSI/Parsl/GDC-Pipeline/config/gdc_config_local.json python ./run.py`




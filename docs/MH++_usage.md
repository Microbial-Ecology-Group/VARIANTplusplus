Usage
-----

### Display Help Message

The `help` parameter displays the available options and commands.
```
nextflow run main_AMR++.nf --help
```

# Parameter selection
AMR++ comes with a default selection of parameters to perform a demonstration using example data provided in the "data/" directory. The example command below uses two types of parameters. 

```
nextflow run main_AMR++.nf -profile singularity --pipeline demo
```

The ```-profile``` parameter only has a single dash (```-```), meaning it corresponds to a nextflow-speific parameter and the  ```--pipeline```. Examples of parameters with one dash include ```-profile```, ```-resume```, and ```-config```. Further details on using the profile parameter, which determines how AMR++ runs on a computing cluster, can be found in the configuration document. Below, we'll see how to change the pipeline-specific parameters which are denoted using two dashes (```--```), such as ```--pipeline``` and ```--reads```.

The AMR++ pipeline pulls information from various sources to determine the correct parameters for running the pipeline. AMR++ is written in nextflow and this allows for us to change how the pipeline runs in a variety of ways. This is the order in which nextflow will prioritize parameters it receives.

1. Parameters specified on the command line (--something value)

2. Parameters provided using the -params-file option (params.config by default)

3. Config file specified using the -c my_config option (e.g. config/local.config)

4. The config file named nextflow.config in the current directory

5. The config file named nextflow.config in the workflow project directory

6. The config file $HOME/.nextflow/config

7. Values defined within the pipeline script itself (e.g. main_AMR++.nf)


## Run Assembly


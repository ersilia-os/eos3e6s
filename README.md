# ChEMBL decoy sampler

This model samples decoy molecules from ChEMBL35 to build negative datasets. We look for molecules with low tanimoto similarity (<= 0.25) but similar physicochemical characteristics (MW, LogP, HDB, HAB, RB). For each input molecule, the model will return 100 decoys. The design has been informed by other methods such as DUD-E and LUDe.



## Information
### Identifiers
- **Ersilia Identifier:** `eos3e6s`
- **Slug:** `chembl-decoys`

### Domain
- **Task:** `Sampling`
- **Subtask:** `Similarity search`
- **Biomedical Area:** `Any`
- **Target Organism:** `Any`
- **Tags:** `ChEMBL`, `Similarity`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `1`
- **Output Consistency:** `Variable`
- **Interpretation:** 100 decoy molecules with similar physicochemical properties and low tanimoto similarity.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| smi_00 | string |  | Sampled molecule index 0 as a decoy |
| smi_01 | string |  | Sampled molecule index 1 as a decoy |
| smi_02 | string |  | Sampled molecule index 2 as a decoy |
| smi_03 | string |  | Sampled molecule index 3 as a decoy |
| smi_04 | string |  | Sampled molecule index 4 as a decoy |
| smi_05 | string |  | Sampled molecule index 5 as a decoy |
| smi_06 | string |  | Sampled molecule index 6 as a decoy |
| smi_07 | string |  | Sampled molecule index 7 as a decoy |
| smi_08 | string |  | Sampled molecule index 8 as a decoy |
| smi_09 | string |  | Sampled molecule index 9 as a decoy |

_10 of 100 columns are shown_
### Source and Deployment
- **Source:** `Local`
- **Source Type:** `Internal`

### Resource Consumption


### References
- **Source Code**: [https://ersilia.io](https://ersilia.io)
- **Publication**: [https://ersilia.io](https://ersilia.io)
- **Publication Type:** `Other`
- **Publication Year:** `2025`
- **Ersilia Contributor:** [GemmaTuron](https://github.com/GemmaTuron)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [GPL-3.0-or-later](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos3e6s
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos3e6s
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!

example data were obtained with the following script

```python
import pandas as pd
from seq2science.util import parse_samples

# load full dataframe, filter for some example columns
for input_file in ["GRCz11-counts.tsv", "GRCz11-TPM.tsv"]:
    df = pd.read_csv(input_file, sep="\t", index_col=0, comment='#')
    samples = pd.read_csv("example/GRCz11_rna_samples.tsv", sep='\t', dtype='str', comment='#')
    samples = parse_samples(samples, dict())
    
    # preserve column order
    cols = [col for col in df.columns if col in set(samples["technical_replicates"])]
    df = df[cols]
    df.to_csv(f"example/{input_file}", sep="\t")
```
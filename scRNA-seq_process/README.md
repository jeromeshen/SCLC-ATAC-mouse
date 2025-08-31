
# CellRanger
## prepare commands
```python
import pandas as pd
import os
import glob
pathin = './scRNA/rawdata'
pathout = './scRNA/cellranger'
cellranger = '/data/p/CellRanger/cellranger-7.1.0/bin/cellranger'
transcriptome = '/data/pub/cellranger/GRCm38_gencode25/'
ls_sample = os.listdir(pathin)
ls_cmds = [f'cd {pathout} && {cellranger} count --id={i} --transcriptome={transcriptome} --sample={i} --fastqs={pathin}{i}' for i in ls_sample]
open('./scRNA/rna_cellranger_count.cmds', 'w').write('\n'.join(ls_cmds) + '\n')
```
## run commands
```shell
bash /scRNA/rna_cellranger_count.cmds
```

# remove ambient RNA molecules
## prepare commands and run
```python
import pandas as pd
import os
import glob
pathin = './scRNA/cellranger'
output = './scRNA/cellbender'
ls_sample = os.listdir(pathin)

files = []
for i inls_sample:
    t = './scRNA/cellranger/'+i+'/outs/raw_feature_bc_matrix.h5'
    files.append(t)

matric = []
for i in ls_sample:
    t = './scRNA/cellranger/'+i+'/outs/metrics_summary.csv'
    matric.append(t)

output_folders = [f'{output}{id}' for id in ls_sample]
for i in output_folders:
    if not os.path.exists(i):
        os.makedirs(i)

df_use = pd.DataFrame()
df_use['sample'] = ls_sample
df_use['raw_matrix_h5'] = files
df_use['output_folders'] = output_folders
df_use['matric'] = matric

ls_cmds = [f'cd {output}{file.split("/")[-3]} && cellbender remove-background --input {file} --output {file.split("/")[-3]}.cellbender.h5 --cuda --epochs 150 --expected-cells {df_use[df_use.raw_matrix_h5 == file].expected_cell_count.iloc[0]}' for file in files]

open('./scRNA/cellbender_remove_background.sh','w').write('\n'.join(ls_cmds))

os.system("bash ./scRNA/cellbender_remove_background.sh")

```

## prepare for Seurat
```python
import pandas as pd
import os
import glob
folder = './scRNA/cellbender'
files = glob.glob(folder+ '/**/*.cellbender_filtered.h5',recursive=True)

ids = [os.path.basename(i).replace('.cellbender_filtered.h5','') for i in files]

ls_cmds = [f'cd {os.path.dirname(i)} && ptrepack --complevel 5 {os.path.basename(i)}:/matrix {os.path.basename(i).replace("_filtered","_filtered_seurat")}:/matrix' for i in files]

open('./scRNA/cellbender_ptrepack.cmds','w').write('\n'.join(ls_cmds))

os.system("bash ./scRNA/cellbender_ptrepack.cmds")

```

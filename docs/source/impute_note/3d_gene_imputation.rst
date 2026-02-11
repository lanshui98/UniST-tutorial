Gene Imputation (3D)
----

.. note::

   The 3D Gene Imputation is computationally intensive and may 
   require several hours to complete, depending on dataset size and 
   hardware configuration.

#### Installation

```
!git clone https://github.com/lanshui98/UniST.git
%cd UniST
!pip install -r requirements.txt
```

Note the the spatial information should be at adata.obsm['spatial'];

Put the adata under external/SUICA_pro/data/.

```
%cd external/SUICA_pro
```

#### Read the data

```
!pip install scanpy
import scanpy as sc

adata = sc.read('path to 3d adata')
genes = adata.var_names
```

Example: get the index of gene "Myl2"

```
gene_idx = genes.get_loc("Myl2")
print(gene_idx)
```

##### Visualization with `Spateo`

```
! pip install spateo-release
import spateo as st
```

```
adata.obs["Myl2"] = adata.X[:,gene_idx].toarray().copy()
```

```
pc, cmap = st.tdr.construct_pc(
    adata=adata,
    spatial_key="spatial",
    groupby="gene1", 
    colormap="hot_r", 
   )
```

```
st.pl.three_d_plot(
    model=pc,
    key="Myl2",  
    colormap="hot_r", 
    model_style="points",
    model_size=4.0,
    show_legend=True,
    jupyter="trame",
    opacity=0.5,
    )
```
<img src="https://github.com/lanshui98/UniST-tutorial/docs/source/impute_note/figs/3D_gene1.png">

"static" for static image, "trame" for interactive window (need to install `nest_asyncio2`).

For more 3D visualization/animation details, please go to [Animation]('https://unist-tutorial.readthedocs.io/en/latest/tutorial.html').

#### Step1: Train GAE

```
! python train.py --mode embedder --conf ./configs/ST/embedder_gae_3d_sparse.yaml
```

#### Step2: Train INR + fine-tune GAE

```
! python train.py --mode inr --conf ./configs/ST/inr_embd_3d_sparse.yaml
```
                                                                      
#### Step3: Prediction/Imputation

Prepare normalized custom coords

```
! python prepare_custom_coords.py --mode 3d --reference data/3D_data.h5ad --coords your_coords.xyz --output data/preprocessed_data/custom_coords_3d_norm.npy --keep_ratio True --preserve_z_scale True --z_scale_factor 1.5
```

Run prediction

```
! python predict.py --mode inr --conf ./configs/ST/inr_pred_3d_sparse.yaml       
```

Map reconstructed coords back to original space

```
! python map_coords_back.py --reconstructed reconstructed-custom-3d.h5ad --reference data/3D_data.h5ad --output reconstructed-original-3d.h5ad --mode 3d --keep_ratio True --preserve_z_scale True --z_scale_factor 1.5
```

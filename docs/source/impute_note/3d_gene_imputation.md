# Gene Imputation (3D)

```{note}
The 3D Gene Imputation is computationally intensive and may require several hours to complete, depending on dataset size and hardware configuration.
```

## Installation

```
!git clone https://github.com/lanshui98/UniST.git
%cd UniST
!pip install -r requirements.txt
```

Note that the spatial information should be at `adata.obsm['spatial']`.

Put the adata under `external/SUICA_pro/data/`.

```
%cd external/SUICA_pro
```

## Read the Data

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

### Visualization with `Spateo`

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
    groupby="Myl2",
    colormap="hot_r"
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
    opacity=0.5
)
```

```{figure} figs/3D_gene1.png
:width: 70%
:align: center
```

- `"static"` for static image
- `"trame"` for interactive window (need to install `nest_asyncio2`)

For more 3D visualization/animation details, please go to [Animation](https://unist-tutorial.readthedocs.io/en/latest/tutorial.html).

---

# Step1: Train GAE

```
! python train.py --mode embedder --conf ./configs/ST/embedder_gae_3d_sparse.yaml
```

These parameters are used to handle **sparse z-direction in 3D spatial transcriptomics data** (e.g., when slice spacing is large).

##### `use_anisotropic_knn: True`
- **Meaning**: Whether to use anisotropic KNN graph construction
- **Effect**: 
  - `True`: Use anisotropic method, accounting for differences between z-direction and xy-directions
  - `False`: Use standard KNN with equal weights for all directions

##### `z_weight: 2.0`
- **Meaning**: Weight factor for z-direction
- **Effect**: 
  - `> 1`: **Reduces** the influence of z-direction distance
  - Calculation: `weighted_z = z / z_weight`
  - For example, `z_weight=2.0` means z-direction distance is **halved**, making z-direction points more likely to become neighbors
- **Principle**: 
  ```
  Original distance: d = sqrt((x1-x2)² + (y1-y2)² + (z1-z2)²)
  Weighted:          d = sqrt((x1-x2)² + (y1-y2)² + (z1-z2)²/z_weight²)
  ```
  
##### `z_threshold: null`
- **Meaning**: Maximum connection distance threshold in z-direction
- **Effect**: 
  - `null/None`: Automatically set to 30% of z-direction range
  - Numeric value: Manually set maximum z-distance (in original coordinate units)
- **Principle**: 
  - Even after weighting, if two points are too far apart in z-direction, they should not be connected
  - For example: points from slice 1 and slice 10 should not be directly connected, even if xy-distance is close
- **Example**:
  ```
  # If z_range = 1000 (from z=0 to z=1000)
  # z_threshold = null → automatically set to 1000 * 0.3 = 300
  # This means points with z-direction distance > 300 will not be connected
  ```

##### `preserve_z_scale: True`
- **Meaning**: Whether to preserve original z-direction scale
- **Effect**:
  - `True`: z-direction is **not compressed**, maintaining a relatively larger range
  - `False`: z-direction is compressed together with xy-directions to the same range
- **Principle**:
  ```
  # preserve_z_scale = False (default)
  # All directions compressed to [-1, 1], maintaining aspect ratio
  scale_x = x_range / max(x_range, y_range, z_range)
  scale_y = y_range / max(x_range, y_range, z_range)
  scale_z = z_range / max(x_range, y_range, z_range)  # z is compressed
  
  # preserve_z_scale = True
  # z-direction maintains larger range, not compressed
  scale_x = x_range / max(x_range, y_range)
  scale_y = y_range / max(x_range, y_range)
  scale_z = z_scale_factor  # z uses amplification factor
  ```
- **Use Case**: When z-direction is sparse, preserve z-direction importance

##### `z_scale_factor: 1.5`
- **Meaning**: Scaling factor for z-direction (only effective when `preserve_z_scale=True`)
- **Effect**:
  - `= 1.0`: z-direction maintains original relative scale
  - `> 1.0`: **Amplifies** z-direction importance (recommended for sparse z-direction)
  - `< 1.0`: Reduces z-direction importance (not recommended)
- **Principle**:
  ```
  normalized_z = (z - z_min) / z_range  # Normalize to [0,1]
  normalized_z = (normalized_z - 0.5) * 2.0  # Transform to [-1,1]
  normalized_z = normalized_z * z_scale_factor  # Apply scaling factor
  ```

# Step2: Train INR + fine-tune GAE

```
! python train.py --mode inr --conf ./configs/ST/inr_embd_3d_sparse.yaml
```
                                                                      
# Step3: Prediction/Imputation

#### Prepare normalized custom coords

```
! python prepare_custom_coords.py --mode 3d --reference data/3D_data.h5ad --coords your_coords.xyz --output data/preprocessed_data/custom_coords_3d_norm.npy --keep_ratio True --preserve_z_scale True --z_scale_factor 1.5
```

#### Run prediction

```
! python predict.py --mode inr --conf ./configs/ST/inr_pred_3d_sparse.yaml       
```

#### Map reconstructed coords back to original space

```
! python map_coords_back.py --reconstructed reconstructed-custom-3d.h5ad --reference data/3D_data.h5ad --output reconstructed-original-3d.h5ad --mode 3d --keep_ratio True --preserve_z_scale True --z_scale_factor 1.5
```

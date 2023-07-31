import scanpy as sc
import scispy as scis

import warnings
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")

# export OPENBLAS_NUM_THREADS=4
# export OMP_NUM_THREADS=4

scale_percent = 10
#dirs = ['FR2302_C5_003','FR2302_C7_003','FR2302_P11_003','FR2302_P14r0_003','FR2302_P14r1_003','FR2302_P16b_003']
#libs = ['C5','C7','P11','P14r0','P14r1','P16b']
#groups = ['CTRL','CTRL','PAH','PAH','PAH','PAH']

dirs = ['FR2302_C5_003','FR2302_P11_003']
libs = ['C5','P11']
groups = ['CTRL','PAH']

size = len(libs)
for i in range(size):

    folder = '/data/analysis/data_lebrigand/000-vizgen/000-DATA/' + dirs[i] + "/"
    library_id = libs[i]
    group = groups[i]

    adata = scis.io.load_merscope(folder,library_id,scale_percent)
    adata = scis.pp.filter_and_run_scanpy(adata, min_counts = 20)

    # htap
    label_ref='celltype'
    scref = sc.read_h5ad('htap_230727.h5ad')
    scref.layers['counts'] = scref.raw.X.copy()
    sc.pp.filter_cells(scref, min_counts=10)
    #scref = scref[scref.obs["group"].str.contains(group)]
    #del scref.var['features']
    scref.shape
    
    # or discovair
    #label_ref='celltype_lv2_V4'
    #scref = sc.read_h5ad('/data/analysis/data_collin/Discovair/data/integrated_V9/integrated_V9.h5ad')
    #scref.layers['counts'] = scref.raw.X.copy()
    #scref.shape
    
    # or HLCA
    #label_ref='ann_finest_level'
    #scref = sc.read_h5ad('/data/analysis/data_collin/Discovair/data/HLCA_v1.h5ad')
    #scref.layers['counts'] = scref.X.copy()
    #scref.shape
    
    # annotation
    adata = scis.pp.annotate(adata, scref, label_ref=label_ref, label_key='celltype', metaref2add='population')

    adata.obs['library_id'] = library_id
    adata.obs['group'] = group

    #sc.pl.embedding(adata, 'umap', color = 'celltype', palette=scis.io.get_palette('celltype'), show = False)
    #sc.pl.embedding(adata, 'umap', color = 'population', palette=scis.io.get_palette('population'), show = False)
    #sc.pl.embedding(adata, 'umap', color = 'group', palette=scis.io.get_palette('group'), show = False)

    adata.write("./outs/ad"+library_id+".htap.h5ad")


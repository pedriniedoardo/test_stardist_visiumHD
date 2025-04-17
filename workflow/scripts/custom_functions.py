import numpy as np
import pandas as pd
import scanpy as sc
import geopandas as gpd
from stardist.models import StarDist2D
from shapely.geometry import Polygon, Point

def custom_normalize(x, pmin=3, pmax=99.8, axis=None, clip=False, eps=1e-20, dtype=np.float32):
    """Percentile-based image normalization."""

    mi = np.array([[[134.]]])#np.percentile(x,pmin,axis=axis,keepdims=True)
    ma = np.percentile(x,pmax,axis=axis,keepdims=True)
    return normalize_mi_ma(x, mi, ma, clip=clip, eps=eps, dtype=dtype)


def normalize_mi_ma(x, mi, ma, clip=False, eps=1e-20, dtype=np.float32):
    if dtype is not None:
        x   = x.astype(dtype,copy=False)
        mi  = dtype(mi) if np.isscalar(mi) else mi.astype(dtype,copy=False)
        ma  = dtype(ma) if np.isscalar(ma) else ma.astype(dtype,copy=False)
        eps = dtype(eps)

    try:
        import numexpr
        x = numexpr.evaluate("(x - mi) / ( ma - mi + eps )")
    except ImportError:
        x =                   (x - mi) / ( ma - mi + eps )

    if clip:
        x = np.clip(x,0,1)

    return x

def remove_destripe_artifacts(adata):
    idx = adata.X.sum(axis=1).A1 == float('inf')
    print('removing '+str(idx.sum())+' dots')
    return adata[idx==False]

def nuclei_detection(polys, tissue_position_file, adata):
    # Creating a list to store Polygon geometries
    geometries = []

    # Iterating through each nuclei in the 'polys' DataFrame
    for nuclei in range(len(polys['coord'])):

        # Extracting coordinates for the current nuclei and converting them to (y, x) format
        coords = [(y, x) for x, y in zip(polys['coord'][nuclei][0], polys['coord'][nuclei][1])]

        # Creating a Polygon geometry from the coordinates
        geometries.append(Polygon(coords))

    # Creating a GeoDataFrame using the Polygon geometries
    gdf = gpd.GeoDataFrame(geometry=geometries)
    gdf['id'] = [f"ID_{i+1}" for i, _ in enumerate(gdf.index)]

    ########################################################################

    # Load the Spatial Coordinates
    df_tissue_positions=pd.read_parquet(tissue_position_file)

    #Set the index of the dataframe to the barcodes
    df_tissue_positions = df_tissue_positions.set_index('barcode')

    # Create an index in the dataframe to check joins
    df_tissue_positions['index']=df_tissue_positions.index

    # Adding the tissue positions to the meta data
    adata.obs =  pd.merge(adata.obs, df_tissue_positions[['pxl_row_in_fullres', 'pxl_col_in_fullres']], left_index=True, right_index=True)

    # Create a GeoDataFrame from the DataFrame of coordinates
    geometry = [Point(xy) for xy in zip(df_tissue_positions['pxl_col_in_fullres'], df_tissue_positions['pxl_row_in_fullres'])]
    gdf_coordinates = gpd.GeoDataFrame(df_tissue_positions, geometry=geometry)

    #########################################################################

    # Perform a spatial join to check which coordinates are in a cell nucleus
    result_spatial_join = gpd.sjoin(gdf_coordinates, gdf, how='left', predicate='within')

    result_spatial_join["in_tissue"] = np.array(result_spatial_join["in_tissue"]).astype(bool)

    # Identify nuclei associated barcodes and find barcodes that are in more than one nucleus
    result_spatial_join['is_within_polygon'] = ~result_spatial_join['index_right'].isna()
    barcodes_in_overlaping_polygons = pd.unique(result_spatial_join[result_spatial_join.duplicated(subset=['index'])]['index'])
    result_spatial_join['is_not_in_an_polygon_overlap'] = ~result_spatial_join['index'].isin(barcodes_in_overlaping_polygons)

    # Remove barcodes in overlapping nuclei
    barcodes_in_one_polygon =  result_spatial_join[result_spatial_join['is_not_in_an_polygon_overlap'] & result_spatial_join["in_tissue"]] # result_spatial_join['is_within_polygon'] &
    # tengo anche quelli unassigned to any nucleus perchè potrebbero essere in seguito aggregati come citoplasma
    # tanto in bin2cell per l'aggregazione dei nuclei gli 0 non verranno aggregati: Integers, with 0 being unassigned to an object.

    # The AnnData object is filtered to only contain the barcodes that are in non-overlapping polygon regions
    filtered_obs_mask = adata.obs_names.isin(barcodes_in_one_polygon['index'])
    filtered_adata = adata[filtered_obs_mask,:]

    # Add the results of the point spatial join to the Anndata object
    filtered_adata.obs =  pd.merge(filtered_adata.obs, barcodes_in_one_polygon[['index','geometry','id','is_within_polygon','is_not_in_an_polygon_overlap']], left_index=True, right_index=True)

    filtered_adata.obs.id = filtered_adata.obs.id.str.replace('ID_', '', regex=False)
    # filtered_adata.obs = filtered_adata.obs.fillna(0)
    filtered_adata.obs.id = filtered_adata.obs.id.fillna(0)
    filtered_adata.obs.id = filtered_adata.obs.id.astype(str)
    print(f"Number of nuclei detected of hires image: {gdf.shape[0]}\nNumber of nuclei detected on VisiumHD slide: {len(np.unique(filtered_adata.obs.id))}")

    # filter and returns also the geometry dataframe associated to nuclei
    gdf['id'] = gdf['id'].str.replace('ID_', '', regex=False)
    gdf = gdf.loc[gdf['id'].isin(filtered_adata.obs.id)]
    gdf.set_index('id', inplace=True)

    return filtered_adata, gdf


def add_obs_variables(adata, unit, species='Hs'):
    x = adata.copy()
    x.obs['counts_per_' + unit] = np.sum(x.X, axis = 1)
    x.obs['features_per_' + unit] = np.sum(x.X >0, axis = 1)

    if "bin_count" in adata.obs.columns:
        x.obs["bin_count_log"] = np.log10(x.obs["bin_count"])
    x.obs["counts_per_" + unit + "_log"] = np.log10(x.obs["counts_per_" + unit])
    x.obs["features_per_" + unit + "_log"] = np.log10(x.obs["features_per_" + unit])

    if species=='Hs':
        x.var["mt"] = x.var_names.str.startswith("MT-")
    if species=='Mm':
        x.var["mt"] = x.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(
        x, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    x.obs.loc[np.isnan(x.obs.pct_counts_mt),'pct_counts_mt'] = 0
    return x

def max_diameter(polygon):
    from scipy.spatial.distance import pdist
    # Extract coordinates of the exterior boundary
    coords = np.array(polygon.exterior.coords)
    # Calculate all pairwise distances
    distances = pdist(coords)
    return np.max(distances)

def cell_geometry(point_geometries):
    
    from scipy.spatial import ConvexHull
    # Transform from point geometry to coordinates
    points = [x.coords[0] for x in point_geometries]
    # Compute the convex hull
    try:
        hull = ConvexHull(points)
    except Exception:
        return Point()
    # Extract the vertices that form the convex hull
    boundary_points = [points[i] for i in hull.vertices]
    
    # Step 1: Calculate the centroid of the boundary_points
    centroid_x = sum(x for x, y in boundary_points) / len(boundary_points)
    centroid_y = sum(y for x, y in boundary_points) / len(boundary_points)
    
    # Step 2: Define a function to calculate the angle from the centroid
    def angle_from_centroid(point):
        x, y = point
        return np.arctan2(y - centroid_y, x - centroid_x)
    
    # Step 3: Sort boundary_points by angle in descending order for clockwise
    boundary_points_sorted_clockwise = sorted(boundary_points, key=angle_from_centroid, reverse=True)
    
    # Ensure the polygon is closed by repeating the first point at the end
    if boundary_points_sorted_clockwise[0] != boundary_points_sorted_clockwise[-1]:
        boundary_points_sorted_clockwise.append(boundary_points_sorted_clockwise[0])
    
    # Create the polygon
    polygon = Polygon(boundary_points_sorted_clockwise)
    
    return polygon

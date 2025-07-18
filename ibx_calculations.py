#####################################################################################
# IBx functions
#####################################################################################
def run_ibx_calc(recursive_df, genotypes, group_json, output, individual_inf_ibx=None, intervals=None):
    df = pd.read_csv(recursive_df)
    if 'subset_replicate' in df.columns:
        df = df.drop(['subset_replicate', 'frac_poly'], axis=1, errors='ignore').drop_duplicates()

    gt = np.load(genotypes)
    # confirm genome files match
    if 'original_nid' in df.columns: 
        original_genomes = max(df['original_nid'].tolist())
        subset_genomes = max(df['recursive_nid'].tolist())
        if original_genomes == gt.shape[0]:
            print("Number of genotypes matches the original genotype array.\
            \nUpdating indices to match the original genome-df.csv index.")
            df['recursive_nid'] = df['original_nid']
        if subset_genomes > gt.shape[0]:
            sys.exit("The number of genotypes specified in the recursive dataframe is greater than the number of genotypes provided.\
            \nThis could be caused by different user specified subset modes upon re-run.\
            \nUpdate the snakemake command with -R generate genotype, root_interval_genome to update.")

    if df.recursive_nid.nunique() < gt.shape[0]:
        print("Genotype array contains", gt.shape[0], "rows when", df.recursive_nid.nunique(), "nodes provided. \
        \nSubsetting the genotype array by recursive node IDs.")
        gt = gt[df['recursive_nid'].tolist()]     


    if intervals is not None:
        intervals = np.load(intervals)

    with open(group_json) as reader:
        group_dict = json.load( reader )    

    # align the gt object to work with C hash functions 
    gt = idm.align_data(gt)
    ibx_index, pop_ibx = ibx_calc(df, gt, group_dict, individual_inf_ibx, intervals)

    ibx_index.to_csv(output[0], index=False)
    save_json(output[1], pop_ibx, indent = 4)
    
    return


def run_ibx_summary(ibx_files, timescale, output):

    ibx_dict = split_ibx_dict(ibx_files, timescale)

    summary_dict = {}
    summary_dict['aggregate'] = get_ibx_summary(ibx_dict['aggregate'])

    # add column names for each summary group
    gt = list(ibx_dict['aggregate'].keys())
    groups = next(iter(ibx_dict['aggregate'].values())).keys()
    for subset in next(iter(ibx_dict['aggregate'].values())).keys():
        summary_dict['aggregate'][subset]['columns']=ibx_dict['aggregate'][list(ibx_dict['aggregate'].keys())[0]][subset]['columns']

    # add the individuals infections, only changes the order  
    if len(ibx_dict['per_infIndex']) > 0:
        summary_dict['per_infIndex'] = get_inf_summary(ibx_dict['per_infIndex'])

    save_json(output, summary_dict, indent=4)
    
    return


def split_ibx_dict(ibx_files, timescale):
    """
    Reverses the order of groupings in ibx dicts for merged IBX json.
    """
    # group IBx files from multiple years
    if '-year' in ibx_files[0]:
        ibx_copy = ibx_files
        gt_sets = set([basename(x).split("-year")[0] for x in ibx_copy])
        ibx_files = [sorted([x for x in ibx_copy if y in x]) for y in gt_sets]

    agg_d, inf_d, merge_d = {}, {}, {}
    for gt_files in ibx_files:
        if isinstance(gt_files, list):
            name = basename(gt_files[0]).split("-year")[0]
            ibx_d = {}
            for i in gt_files:
                with open(i) as reader:
                    json_data = json.load( reader )
                ibx_d = merge(ibx_d, json_data)   
        else:
            name = basename(gt_files).split("-ibx")[0]
            with open(gt_files) as reader:
                ibx_d = json.load(reader)
        
        agg_d[name] = ibx_d['aggregate']
 
        if 'per_infIndex' in ibx_d.keys():
            inf_d[name] = ibx_d['per_infIndex']
        
    merge_d['aggregate'], merge_d['per_infIndex'] = agg_d, inf_d

    return merge_d    
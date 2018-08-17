import glob
import os
import pandas as pd
import numpy as np

def extract_genus_name(df):
    ''' Get the text string associated with the genus name by parsing 
        the "name" field from IDSeq output 
    
    Args: 
        df: the input dataframe from IDSeq
    Returns:
        a dataframe with additional 'genus_name' column
    '''
    
    names = df['name']
    genera = []
    for i in names:
        try:
            if 'Non-species-specific' in i:
                genera.append(i.split(' ')[-1])
            else:
                genera.append(i.split(' ')[0])
        except:
            genera.append(None)
    df['genus_name'] = genera
    return(df)
    
def remove_non_genus_specific_virus(df):
    ''' Viruses are evaluated only at the genus level - so remove all 
        "non-genus-specific" rows from IDSeq output 
    
    Args: 
        df: the input dataframe from IDSeq
    Returns: 
        a dataframe with rows where "name" contains string "Non-genus-specific", removed
    '''
    
    names = df['name']
    remove = []
    for i in names:
        if 'Non-genus-specific' in i:
            remove.append(True)
        else:
            remove.append(False)
            
    df['remove'] = remove
    df = df[df['remove'] == False]
    return(df)

def mark_possibly_pathogenic(df, reference_pathogens):
    ''' Match names against list of reference pathogens to mark the potentially 
        pathogenic species (in new column "possible_pathogen")
        
    Args: 
        df: the input dataframe from IDSeq
    Returns:
        a dataframe with additional 'possible_pathogen' column
    '''
    possibly_pathogen = []    
    for i in df.index:
        try:
            name = df.loc[i]['name'].split(' ')
            if(len(name) > 1):
                if ' '.join(name[0:2]) in reference_pathogens:
                    possibly_pathogen.append(True)
                elif df.loc[i]['genus_name'] in reference_pathogens:
                    possibly_pathogen.append(True)
                else:
                    possibly_pathogen.append(False)
            else:
                if name[0] in reference_pathogens:
                    possibly_pathogen.append(True)
                else:
                    possibly_pathogen.append(False)        
        except:
            possibly_pathogen.append(False)
    df['possible_pathogen'] = possibly_pathogen
    return(df)


def disperse_genus_rpm(df):
    
    final_df = []
    
    genus_taxids = list(set(df['genus_taxid']))
    for gtid in genus_taxids:
        
        sub = df[df['genus_taxid'] == gtid]
        to_be_dispersed = sub[sub['tax_id'] < 0]
        to_be_dispersed = to_be_dispersed[to_be_dispersed['genus_taxid'] > 0 ]
        
        keep = sub[sub['tax_id'] > 0]
        
        if to_be_dispersed.shape[0] > 0:

            disperse = to_be_dispersed.iloc[0]['NT_rpm']
            taxid = to_be_dispersed[['genus_taxid']]

            denom_keep = sum(keep['NT_rpm'])  # total sum of species rpm
            
            if denom_keep > 0:
                new_rpm = []
                for i in keep.index:
                    new_rpm.append(float(keep.loc[i]['NT_rpm']) + ((float(keep.loc[i]['NT_rpm'])/denom_keep)*disperse))

                keep['NT_rpm'] = new_rpm

        if(keep.shape[0] > 0):
            final_df.append(keep)
        
    return(pd.concat(final_df))


def get_microbe_stats(data_directory, reference_pathogens):
    stats_dict_overall = {}
    
    for file in glob.glob(data_directory): 

        #print(file)
        stats_dict = {}

        df = pd.read_csv(file)#, index_col=0)

        # separate out viruses up front
        df_vir = df[df['category_name'] == 'Viruses'] 
        df_vir = df_vir[df_vir['tax_level'] == 2] # use only genus-level virus counts
        df_vir = remove_non_genus_specific_virus(df_vir) 

        df_vir = mark_possibly_pathogenic(df_vir, reference_pathogens)
        stats_dict['virus_total'] = sum(df_vir['NT_rpm'])
        stats_dict['virus_pathogenic_total'] = sum(df_vir[df_vir['possible_pathogen'] == True]['NT_rpm'])

        #print(df_vir[['name','possible_pathogen','NT_rpm']])

        # process to remove other categories
        df = df[df['is_phage'] == 0]  # remove phage
        df = df[df['tax_level'] == 1]  # select only the species-level hits (for bact and eukaryotes)

        # separate out bacteria
        df_bact = df[df['category_name'] == 'Bacteria']
        df_bact = extract_genus_name(df_bact)
        df_bact = mark_possibly_pathogenic(df_bact, reference_pathogens)
        df_bact = disperse_genus_rpm(df_bact)

        stats_dict['bact_total'] = sum(df_bact['NT_rpm'])
        stats_dict['bact_pathogenic_total'] = sum(df_bact[df_bact['possible_pathogen'] == True]['NT_rpm'])

        # separate out all other microbes (eukaryotes)
        df_other = df[df['category_name'] != 'Bacteria']
        df_other = df_other[df_other['category_name'] != 'Viruses']
        df_other = extract_genus_name(df_other)
        df_other = mark_possibly_pathogenic(df_other, reference_pathogens)
        df_other = disperse_genus_rpm(df_other)

        stats_dict['other_total'] = sum(df_other['NT_rpm'])
        stats_dict['other_pathogenic_total'] = sum(df_other[df_other['possible_pathogen'] == True]['NT_rpm'])

        stats_dict_overall[file] = stats_dict
        
    return(stats_dict_overall)

def summarize_results(stats_dict_overall, studyID):
    summary_df = pd.DataFrame(stats_dict_overall)

    if studyID == 'BANG':
        summary_df.columns = [i.split('/')[1].split('.')[0].upper() for i in summary_df.columns]
    elif studyID == 'MBAL':
        summary_df.columns = [i.split('/')[-1].split('.')[0].upper() for i in summary_df.columns]
    else:
        return("Did not recognize studyID")
    summary_df = summary_df.transpose()
    summary_df['B'] = summary_df['bact_pathogenic_total']/summary_df['bact_total']
    summary_df['V'] = summary_df['virus_pathogenic_total']/summary_df['virus_total']
    summary_df['O'] = summary_df['other_pathogenic_total']/summary_df['other_total']

    total_microbial = (summary_df['bact_total'] + summary_df['virus_total'] + summary_df['other_total'])
    summary_df['B_tot'] = summary_df['bact_pathogenic_total'] / total_microbial
    summary_df['V_tot'] = summary_df['virus_pathogenic_total'] / total_microbial
    summary_df['O_tot'] = summary_df['other_pathogenic_total'] / total_microbial

    summary_df = summary_df.fillna(0)

    summary_df['log_pathogenic_bact_total'] = np.log10(summary_df['bact_pathogenic_total']+1)
    summary_df['log_pathogenic_virus_total'] = np.log10(summary_df['virus_pathogenic_total']+1)
    summary_df['log_pathogenic_other_total'] = np.log10(summary_df['other_pathogenic_total']+1)

    return(summary_df)
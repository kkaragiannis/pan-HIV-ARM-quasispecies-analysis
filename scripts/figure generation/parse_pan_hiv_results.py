#!/home/karagiannisk/bin/anaconda_ete/bin/python
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import seaborn as sns
import time
import re
import os
import csv
import random
import sys
import pdb
import getopt
import json
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance
from ete2 import ClusterTree, TreeStyle
from itertools import combinations

# from PyQt5.QtWidgets import QApplication, QWidget


sns.set(style="whitegrid", context="paper", font_scale=1.5)
output_path = ""
plt.switch_backend('TkAgg')


def is_float(my_input):
    try:
        num = float(my_input)
    except ValueError:
        return False
    return True


def is_int(my_input):
    try:
        num = int(my_input)
    except ValueError:
        return False
    return True


def split_data_frame_list(df, target_column, separator):
    """" df = dataframe to split,
    target_column = the column containing the values to split
    separator = the symbol used to perform the split
    returns: a dataframe with each entry for the target column separated, with each element moved into a new row. 
    The values in the other columns are duplicated across the newly divided rows.
    """
    def split_list_to_rows(row, row_accumulator, m_target_column, m_separator):
        if pd.isna(row[m_target_column]):
            pass
        else:
            split_row = row[m_target_column].split(m_separator)
            for s in split_row:
                new_row = row.to_dict()
                new_row[m_target_column] = s
                row_accumulator.append(new_row)
    new_rows = []
    df.apply(split_list_to_rows, axis=1, args=(new_rows, target_column, separator))
    new_df = pd.DataFrame(new_rows)
    return new_df


def fix_consensus_names(my_df, column_name):
    for i, row in my_df.iterrows():
        names = row[column_name]
        new_names = []
        pures = []
        crfs = []
        sbtps = 0
        if not isinstance(names, str):
            continue
        for nm in names.split(','):
            if sbtps >= 2:
                break
            sbtps += 1
            nm = nm.replace('CONSENSUS_', '')
            if "_" in nm:
                nm = "CRF"+nm
                crfs.append(nm)
            else:
                pures.append(nm[0])
                new_names.append(nm)
        
        for crf in crfs:
            do_break = False
            for pure in pures:
                if pure in crf[5:]:
                    do_break = True
                    break
            if not do_break:
                new_names.append(crf)

        my_df.set_value(i, column_name, ','.join(sorted(new_names)))
    return my_df


# bar chart of number of runs for each technology
def draw_runs_per_technology(my_df):
    global output_path
    func_name = sys._getframe().f_code.co_name
    plot_name = func_name.replace("draw_", "")

    my_df = my_df.drop_duplicates(subset=['srrID'], keep='first')

    sns.set(style="whitegrid", context="paper", font_scale=1.8)

    g = sns.countplot(x='Platform', data=my_df)
    # plt.title('Samples per technology')

    for v in g.get_yticklabels():
        v.set_rotation(45)

    g.set_ylabel("# of samples")
    g.set_xlabel("Technology")

    for p in g.patches:
        y = p.get_bbox().get_points()[1, 1]
        x = p.get_bbox().get_points()[:, 0]
        g.annotate(str(int(y)), (x.mean(), y), ha='center', va='bottom')

    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + func_name
    plt.clf()
    return True


# pie chart. number of runs for each genotype
def draw_runs_per_subtype(my_df):
    global output_path

    my_df = my_df.drop_duplicates(subset=['srrID'], keep='first')

    cm_df = pd.concat([pd.Series(row['Final Reference names'].split(',')) for _, row in my_df.iterrows() if
                       not pd.isna(row['Final Reference names'])]).reset_index()

    cm_df.columns = ['index', 'genotype']

    func_name = sys._getframe().f_code.co_name
    plot_name = func_name.replace("draw_", "")

    sns.set(style="whitegrid", context="paper", font_scale=1.5)

    g = sns.countplot(x='genotype', data=cm_df, order=cm_df['genotype'].value_counts().index)

    # plt.title('Distribution of Forms')

    for v in g.get_yticklabels():
        v.set_rotation(90)

    for v in g.get_xticklabels():
        v.set_rotation(90)

    g.set_ylabel("# of samples")
    g.set_xlabel("Forms")

    for p in g.patches:
        y = p.get_bbox().get_points()[1, 1]
        x = p.get_bbox().get_points()[:, 0]
        g.annotate(str(int(y)), (x.mean(), y), ha='center', va='bottom', rotation=45)

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + func_name
    plt.clf()
    return True


# categorical. Count of each genotype per technology
def draw_subtypes_per_technology(my_df):
    global output_path

    my_df = my_df.drop_duplicates(subset=['srrID'], keep='first')

    cm_df = pd.concat(
        [pd.Series(row['Platform'], row['Final Reference names'].split(',')) for _, row in my_df.iterrows() if
         not pd.isna(row['Final Reference names'])]).reset_index()

    cm_df.columns = ['genotype', 'Platform']

    func_name = sys._getframe().f_code.co_name
    plot_name = func_name.replace("draw_", "")

    my_order = ['B', 'CRF02_AG', 'C', 'F1', 'CRF01_AE', 'CRF07_BC', 'A1', 'D', 'CPZ', 'A2', 'CRF14_BG',
                'CRF11_CPX', 'G', 'CRF08_BC', 'CRF10_CD', 'CRF12_BF', 'CRF06_CPX', 'O', 'H', 'CRF04_CPX']

    sns.set(style="whitegrid", context="paper", font_scale=1.5)

    g = sns.factorplot(x='genotype', data=cm_df,
                       col='Platform',
                       sharey=False,
                       kind='count',
                       order=my_order,
                       col_wrap=2)

    # plt.title('Distribution of subtypes')

    g.set_xticklabels(rotation='90')
    # for v in g.get_yticklabels():
    #     v.set_rotation(90)
    #
    # for v in g.get_xticklabels():
    #     v.set_rotation(90)

    g.set_axis_labels("", "# of samples")


    # plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + func_name
    plt.clf()
    return True


# categorical (violin). Unaligned Percentage technology
def draw_unaligned_per_technology(my_df):
    global output_path

    func_name = sys._getframe().f_code.co_name
    plot_name = func_name.replace("draw_", "")

    my_df = my_df.drop_duplicates(subset=['srrID'], keep='first')

    def conv_to_perc(df, my_index, my_row):
        if is_float(my_row['Quality']) and is_float(my_row['Total Reads']):
            df.at[my_index, 'Quality'] = 100 * float(my_row['Quality']) / float(my_row['Total Reads'])
        if is_float(my_row['Final']) and is_float(my_row['Total Reads']):
            df.at[my_index, 'Final'] = 100 * float(my_row['Final']) / float(my_row['Total Reads'])
            if df.at[my_index, 'Final'] < 0 or df.at[my_index, 'Quality'] < 0:
                print " boom "

    sns.set(style="whitegrid", context="paper", font_scale=1.5)
    c_df = my_df[['Platform', 'Total Reads', 'Unaligned (qual)', 'Unaligned (final)']]
    c_df.columns = ['Platform', 'Total Reads', 'Quality', 'Final']
    [conv_to_perc(c_df, m_index, m_row) for m_index, m_row in c_df.iterrows()]

    c_df = c_df.melt(['Platform', 'Total Reads'])
    c_df.columns = ['Platform', 'Total Reads', 'Step', 'Unaligned']
    g = sns.violinplot(x='Platform',
                       y='Unaligned',
                       data=c_df,
                       hue='Step',
                       cut=0,
                       scale='count', scale_hue=True,
                       split=True)

    g.set_xlabel("Platform")
    g.set_ylabel("Unaligned reads (%)")
    g.legend(loc="upper center")

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + func_name
    plt.clf()
    return True


# categorical . Div stats per technology
def draw_nucleotide_diversity_per_technology(my_df):
    global output_path

    func_name = sys._getframe().f_code.co_name
    plot_name = func_name.replace("draw_", "")

    my_df = my_df.drop_duplicates(subset=['srrID', 'assembly'], keep='first')

    sns.set(style="whitegrid", context="paper", font_scale=1.5)
    c_df = my_df[['Platform', 'sample diversity', 'assembly']]
    c_df.columns = ['Platform', 'Sample diversity', 'Assembly']

    g = sns.violinplot(x='Platform',
                       y='Sample diversity',
                       data=c_df,
                       hue='Assembly',
                       cut=0,
                       split=True)

    g.set_xlabel("Platform")
    g.set_ylabel("Sample diversity")

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + func_name
    plt.clf()
    return True


# categorical . Div stats per genotype
def draw_diversity_per_subtype(my_df, metric, title, plot_name, y_group="Final Reference names",
                               y_group_title="Forms", which_assembly="both"):
    global output_path

    # func_name = sys._getframe().f_code.co_name
    # plot_name = func_name.replace("draw_", "")

    my_df = my_df.drop_duplicates(subset=['srrID', 'assembly'], keep='first')

    sns.set(style="whitegrid", context="paper", font_scale=1.5)
    c_df = my_df[['srrID', y_group, metric, 'assembly']]
    if y_group == 'Final Reference names':
        c_df = split_data_frame_list(c_df, y_group, ',')

    c_df = c_df[['srrID', "assembly", metric, y_group]]
    # c_df = c_df.set_index(['srrID', metric, 'assembly'])['Final Reference names'].apply(pd.Series).stack()
    # c_df = c_df.reset_index()
    c_df.columns = ['srrID', 'Assembly', metric, y_group_title]
    c_df = c_df[c_df[metric].notnull()]
    if which_assembly == "both":
        g = sns.swarmplot(x=y_group_title, y=metric, data=c_df,
                          hue='Assembly',
                           # cut=0,
                           split=True)
    else:
        if which_assembly == "global" or which_assembly == "contig":
            c_df = c_df[c_df['Assembly'] == which_assembly]

        my_order = ['B','CRF02_AG','C','F1','CRF01_AE','CRF07_BC','A1','D','CPZ','A2','CRF14_BG',
                    'CRF11_CPX','G','CRF08_BC','CRF10_CD','CRF12_BF','CRF06_CPX','O','H','CRF04_CPX']

        if y_group == 'Platform':
            my_order = ["LS454", "ILLUMINA", "PACBIO_SMRT"]

        g = sns.swarmplot(x=y_group_title, y=metric,
                          order=my_order, data=c_df)
    g.set_xlabel(y_group_title)
    g.set_ylabel(title)
    for v in g.get_xticklabels():
        v.set_rotation(90)

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True


def draw_subtype_coocurrence(my_df):
    global output_path

    func_name = sys._getframe().f_code.co_name
    plot_name = func_name.replace("draw_", "")

    # This is a Neighbour-joining tree without distance corrections.
    tr = ClusterTree(tree_path)
    idx_dict = {
        'O': 0,
        'CPZ': 1,
        'C': 2,
        'CRF07_BC': 3,
        'CRF08_BC': 4,
        'F1': 5,
        'CRF12_BF': 6,
        'B': 7,
        'D': 8,
        'CRF10_CD': 9,
        'CRF11_CPX': 10,
        'A2': 11,
        'CRF02_AG': 12,
        'A1': 13,
        'CRF01_AE': 14,
        'CRF06_CPX': 15,
        'G': 16,
        'CRF14_BG': 17,
        'CRF04_CPX': 18,
        'H': 19
    }

    sbtp_zeroed = dict(idx_dict)
    for ii in sbtp_zeroed:
        sbtp_zeroed[ii] = 0

    sbtp_dict = {}
    for m_sb in idx_dict:
        sbtp_dict[m_sb] = dict(sbtp_zeroed)

    def get_linkage_distance(m_tr, m_idx_d):
        leaves = m_tr.get_leaf_names()

        dmat = np.zeros((len(m_idx_d),len(m_idx_d)))

        for l1, l2 in combinations(leaves, 2):
            d = m_tr.get_distance(l1, l2)
            dmat[m_idx_d[l1], m_idx_d[l2]] = dmat[m_idx_d[l2], m_idx_d[l1]] = d

        return dmat

    lnk_matrix = get_linkage_distance(tr, idx_dict)
    schlink = sch.linkage(scipy.spatial.distance.squareform(lnk_matrix), method='average', metric='euclidean')

    def get_coocurrences(my_ref_str):
        if pd.isna(my_ref_str) is not True:
            my_ref_array = my_ref_str.split(',')
            if len(my_ref_array) <= 1 :
                return
            for c_ref in my_ref_array:
                for ic_ref in my_ref_array:
                    if ic_ref == c_ref:
                        continue
                    sbtp_dict[c_ref][ic_ref] += 1

    [get_coocurrences(m_row['Final Reference names']) for m_index, m_row in my_df.iterrows()]

    sbtp_tree_idx = [idx_dict.keys()[idx_dict.values().index(i)] for i in range(0, len(idx_dict))]
    mcl_df = pd.DataFrame(data=sbtp_dict)
    mcl_df = mcl_df.reindex(sbtp_tree_idx, columns=sbtp_tree_idx)
    g = sns.clustermap(mcl_df, row_linkage=schlink, col_linkage=schlink, linewidths=.2, annot=True, fmt='d',
                       cmap="Reds", figsize=(12, 12 ))

    plt.savefig(os.path.join(output_path, plot_name + ".png"))
    print "Done with " + func_name
    plt.clf()
    return True


# categorical . Counts of drug rest mutations frequency percentage bins
def draw_drm_frequency_bins(my_df, plot_name, my_bins=range(0, 110, 10), percentage=False):
    global output_path

    my_df = my_df.drop_duplicates(subset=["srrID", "assemly", "sequence"])

    def get_freq_bins(m_df, bins, m_perc):
        res = m_df[["assembly"]]
        res = pd.concat([res,pd.get_dummies(pd.cut(m_df["freq"], bins))],axis=1).\
            groupby(res["assembly"]).sum().reset_index()

        if m_perc:
            res = pd.concat([res["assembly"], res[res.columns[1:]].apply(lambda x : (100.*x)/x.sum(),axis=1)], axis=1)

        res = res.melt("assembly")
        res.columns = ["Assembly", "Range", "# of ARM"]
        return res

    c_df = get_freq_bins(my_df, my_bins, percentage)

    sns.set(style="whitegrid", context="paper", font_scale=1.5)
    # c_df = my_df[['Platform', 'sample diversity', 'assembly']]
    # c_df.columns = ['Platform', 'Sample diversity', 'Assembly']

    g = sns.barplot(x='Range',
                    y='# of ARM',
                    data=c_df,
                    hue='Assembly')

    if not percentage:
        g.set_yscale("log")
        g.set_ylabel("% of ARM")
    else:
        g.set_ylabel("Percentage of ARM (%)")
    g.set_xlabel("Abundance (%)")

    for v in g.get_xticklabels():
        v.set_rotation(45)

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True


# categorical . Counts of drug rest mutations frequency percentage bins per drug class and drug
def draw_drm_frequency_bins_per_class(my_df, plot_name, my_bins=range(0, 110, 10), percentage=False):
    global output_path

    my_df = my_df.drop_duplicates(subset=["srrID", "assembly", "sequence", "class"])

    def get_freq_bins(m_df, bins, m_perc):
        key_cols_array = ["assembly", "class"]
        t_res = m_df[key_cols_array]
        res = pd.concat([t_res, pd.get_dummies(pd.cut(m_df["freq"], bins))], axis=1).groupby(key_cols_array).sum().reset_index()

        if m_perc:
            res = pd.concat([res[key_cols_array],
                             res[res.columns[len(key_cols_array):]]
                            .apply(lambda x: (100.*x)/x.sum(),axis=1)], axis=1)

        res = res.melt(key_cols_array)
        res.columns = ["Assembly", "Class", "Range", "# of ARM"]
        return res

    c_df = get_freq_bins(my_df, my_bins, percentage)

    sns.set(style="whitegrid", context="paper", font_scale=1.5)
    g = sns.catplot(data=c_df, kind='bar', x='Range', y='# of ARM', hue='Assembly',
                    col='Class', col_wrap=2)

    if percentage:
        g.set_axis_labels("Abundance (%)", "# of ARM (%)")
    else:
        g.set_axis_labels("Abundance (%)", "# of ARM")
        g.set(yscale="log")
    g.set_xticklabels(rotation='90')

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True


# categorical . Counts of drug rest mutations frequency percentage bins per platform
def draw_drm_frequency_bins_per_platform(my_df, plot_name, my_bins=range(0, 110, 10), percentage=False):
    global output_path

    my_df = my_df.drop_duplicates(subset=["srrID", "assembly", "sequence", "Platform"])

    def get_freq_bins(m_df, bins, m_perc):
        key_cols_array = ["assembly", "Platform"]
        t_res = m_df[key_cols_array]
        res = pd.concat([t_res, pd.get_dummies(pd.cut(m_df["freq"], bins))], axis=1).groupby(key_cols_array).sum().reset_index()

        if m_perc:
            res = pd.concat([res[key_cols_array],
                             res[res.columns[len(key_cols_array):]]
                            .apply(lambda x: (100.*x)/x.sum(),axis=1)], axis=1)

        res = res.melt(key_cols_array)
        res.columns = ["Assembly", "Platform", "Range", "# of ARM"]
        return res

    c_df = get_freq_bins(my_df, my_bins, percentage)

    sns.set(style="whitegrid", context="paper", font_scale=1.5)
    g = sns.catplot(data=c_df, kind='bar', x='Range', y='# of ARM', hue='Assembly',
                    col='Platform', col_wrap=2, legend=False)
    # g.fig.legend(loc='lower right')
    handles,lines = g.axes[0].get_legend_handles_labels()
    g.fig.legend(handles,lines, title="Assembly", bbox_to_anchor=(0.75, 0.25), loc='center')
    # g.fig.legend(loc='lower right', bbox_to_anchor=(1, 1.1), ncol=2,
    #           borderaxespad=0, frameon=False)


    if percentage:
        g.set_axis_labels("Abundance (%)", "# of ARM (%)")
    else:
        g.set_axis_labels("Abundance (%)", "# of ARM")
        g.set(yscale="log")
    g.set_xticklabels(rotation='60')

    plt.tight_layout(pad=1)
    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True


# categorical . Counts of drug rest mutations frequency percentage bins per platform
def draw_drm_frequency_abundance_scatter_on_level_and_class(my_df, plot_name, how='min'):
    global output_path

    key_cols_array = ["srrID", "assembly", "sequence", "class"]
    my_df = my_df.dropna(subset=["level","freq"])

    # my_df = my_df.drop_duplicates(subset=["srrID", "assembly", "sequence", "Platform"])
    if how == 'max':
        my_df = my_df.groupby(key_cols_array, group_keys=False).\
            apply(lambda x: x.loc[x['level'].idxmax()]).reset_index(drop=True)
    elif how == "mean":
        my_df = my_df.groupby(key_cols_array, group_keys=False, as_index=False).mean()
    else:
        my_df = my_df.groupby(key_cols_array, group_keys=False).\
            apply(lambda x: x.loc[x['level'].idxmin()]).reset_index(drop=True)

    sns.set(style="ticks", font_scale=1.5)
    if how == "mean":
        g = sns.FacetGrid(my_df, hue='assembly', hue_kws={"marker": ["s", "+"]},  col='class', col_wrap=2,
                          palette='Set2', legend_out=True)
        g.map(plt.scatter, 'level', 'freq')
        g.add_legend()
    else:
        g = sns.catplot(data=my_df, y='freq', x='level', hue='assembly', kind='violin',
                        inner='stick', col='class', col_wrap=2, palette='Set2')

    g._legend.set_title("Assembly")
    if how == "mean":
        g.set_axis_labels("Severity level", "Abundance (%)")
    else:
        g.set_axis_labels("Severity level", "Abundance (%)")

    plt.ylim(0,110)
    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True


# categorical . Counts of drug rest mutations frequency percentage bins per platform
def draw_drm_frequency_abundance_scatter_on_level_and_subtype(my_df, plot_name, how='min'):
    global output_path
    my_df.rename(columns={'Final Reference names': 'Form'}, inplace=True)
    my_df = split_data_frame_list(my_df, 'Form', ',')
    key_cols_array = ["srrID", "assembly", "sequence", "Form"]
    my_df = my_df.dropna(subset=["level","freq"])

    # my_df = my_df.drop_duplicates(subset=["srrID", "assembly", "sequence", "Platform"])
    if how == 'max':
        my_df = my_df.groupby(key_cols_array, group_keys=False).\
            apply(lambda x: x.loc[x['level'].idxmax()]).reset_index(drop=True)
    elif how == "mean":
        my_df = my_df.groupby(key_cols_array, group_keys=False, as_index=False).mean()
    else:
        my_df = my_df.groupby(key_cols_array, group_keys=False).\
            apply(lambda x: x.loc[x['level'].idxmin()]).reset_index(drop=True)

    sns.set(style="ticks", font_scale=1.5)
    if how == "mean":
        g = sns.FacetGrid(my_df, hue='assembly', hue_kws={"marker": ["s", "+"]},  col='Form', col_wrap=4,
                          palette='Set2', legend_out=True)
        g.map(plt.scatter, 'level', 'freq')
        g.add_legend()
    else:
        g = sns.catplot(data=my_df, y='freq', x='level', hue='assembly', kind='violin',
                    inner='stick', col='Form', col_wrap=4, palette='Set2')

    g._legend.set_title("Assembly")
    if how == "mean":
        g.set_axis_labels("Severity level", "Abundance (%)")
    else:
        g.set_axis_labels("Severity level", "Abundance (%)")

    plt.ylim(0,110)
    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True


# categorical . Counts of drug rest mutations vs number of mutations lengths and severity
def draw_drm_severity_and_number_of_mutations_scatter_on_synergistic(my_df, plot_name, how='min'):
    global output_path
    # my_df.rename(columns={'Final Reference names': 'Form'}, inplace=True)
    # # my_df = split_data_frame_list(my_df, 'Form', ',')
    # key_cols_array = ["srrID", "assembly", "sequence", "mutations", "length of mutations"]
    # my_df = my_df.dropna(subset=["level", "freq"])
    # my_df = my_df.groupby(key_cols_array, group_keys=False)
    #
    # # my_df = my_df.drop_duplicates(subset=["srrID", "assembly", "sequence", "Platform"])
    # if how == 'max':
    #     t_df = my_df.apply(lambda x: x.loc[x['level'].idxmax()]).reset_index(drop=True)
    # elif how == "mean":
    #     t_df = my_df.mean()
    # else:
    #     t_df = my_df.apply(lambda x: x.loc[x['level'].idxmin()]).reset_index(drop=True)
    #
    # my_df = pd.concat([t_df, my_df.size().reset_index(name='count')['count']], axis=1)

    my_df = pd.read_csv("E:\\app\\HexaHedron_scripts\\pan_HIV_analysis\\results\\lanl_runs_subtypes_lengths_and_mutations_count_fixed.csv")
    my_df.rename(columns={'Subtype': 'Form'}, inplace=True)
    sns.set(style="whitegrid", context="paper", font_scale=2.5)
    if how == "mean":
        g = sns.FacetGrid(my_df, hue='assembly', hue_kws={"marker": ["s", "+"]}, col='Form', col_wrap=4,
                          palette='Set2', legend_out=True)
        g.map(plt.scatter, 'level', 'freq')
        g.add_legend()
    else:
        g = sns.catplot(data=my_df, kind='violin', y='count', x='length of mutations', hue='assembly',
                        col='level', col_wrap=3, palette='Set2', sharex=False, sharey=False)

        # g.map(plt.scatter, 'count', 'length of mutations')
        # g = sns.scatterplot(data=my_df, x='count', y='length of mutations', hue='assembly', size='level', palette='Set2')

    # g._legend.set_title("Assembly")
    # if how == "mean":
    # else:
    #     g.set_axis_labels("Severity level", "Abundance (%)")

    g.set_axis_labels("Cardinality of mutation sets", "# of mutation sets")
    for ax in g.axes:
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=14)

    # plt.ylim(0, 110)
    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True


# categorical . Counts of drug rest mutations vs number of mutations lengths and severity
def draw_drm_severity_and_cardinatlity_of_mutations(my_df, plot_name, how='min'):
    global output_path
    # my_df.rename(columns={'Final Reference names': 'Form'}, inplace=True)
    # # my_df = split_data_frame_list(my_df, 'Form', ',')
    # key_cols_array = ["srrID", "assembly", "sequence", "mutations", "length of mutations"]
    # my_df = my_df.dropna(subset=["level", "freq"])
    # my_df = my_df.groupby(key_cols_array, group_keys=False)
    #
    # # my_df = my_df.drop_duplicates(subset=["srrID", "assembly", "sequence", "Platform"])
    # if how == 'max':
    #     t_df = my_df.apply(lambda x: x.loc[x['level'].idxmax()]).reset_index(drop=True)
    # elif how == "mean":
    #     t_df = my_df.mean()
    # else:
    #     t_df = my_df.apply(lambda x: x.loc[x['level'].idxmin()]).reset_index(drop=True)
    #
    # my_df = pd.concat([t_df, my_df.size().reset_index(name='count')['count']], axis=1)
    my_df = pd.read_csv("E:\\app\\HexaHedron_scripts\\pan_HIV_analysis\\results\\lanl_runs_subtypes_lengths_and_mutations_count_fixed.csv")
    my_df = my_df.groupby(['level', 'length of mutations','assembly'], group_keys=False)['count'].\
        sum().reset_index(name='count')
    my_df.rename(columns={'level':'Severity level','Subtype':'Form'}, inplace=True)
    sns.set(style="whitegrid", context="paper", font_scale=1.5)
    if how == "mean":
        g = sns.FacetGrid(my_df, hue='assembly', hue_kws={"marker": ["s", "+"]}, col='Form', col_wrap=4,
                          palette='Set2', legend_out=True)
        g.map(plt.scatter, 'level', 'freq')
        g.add_legend()
    else:
        g = sns.FacetGrid(data=my_df, col='assembly', hue='Severity level',
                          legend_out=True, palette='Set2', sharex=False, sharey=False)
        g.map(plt.plot, 'length of mutations', 'count', marker=".").add_legend()

        # g = sns.catplot(data=my_df, kind='line', y='count', x='length of mutations', hue='assembly',
        #                 col='level', col_wrap=3, palette='Set2', sharex=False, sharey=False)
        # g = sns.relplot(data=my_df, kind='line', y='count', x='length of mutations', hue='level',
        #                 col='assembly', palette='Set2')
        # g.map(plt.scatter, 'count', 'length of mutations')
        # g = sns.scatterplot(data=my_df, x='count', y='length of mutations', hue='assembly', size='level', palette='Set2')

    # g._legend.set_title("Severity")
    # if how == "mean":
    # else:
    #     g.set_axis_labels("Severity level", "Abundance (%)")

    g.set_axis_labels("Cardinality of mutation set", "# of mutation sets")

    # plt.ylim(0, 110)
    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True


def draw_drm_severity_per_sample(my_df, plot_name, how='min'):
    global output_path

    cached_fl = "E:\\app\\HexaHedron_scripts\\pan_HIV_analysis\\results\\lanl_runs_min_levels_and_frequencies_per_sample.csv"
    t_df=""
    if os.path.exists(cached_fl):
        t_df = pd.read_csv(cached_fl)
    else:
        def add_missing_freqs(x):
            no_mutations_freq = 100 - x.groupby('sequence')['freq'].first().fillna(0).sum()
            if no_mutations_freq == 100:
                x.iat[0, 11] = 0
                x.iat[0, 89] = 100
            elif no_mutations_freq > 0:
                last_row = x.iloc[-1]
                last_row['freq'] = no_mutations_freq
                last_row['level'] = 0
                last_row['sequence'] = '-'
                x = x.append(last_row)
            return x

        my_df = my_df.dropna(subset=["sequence"])
        my_df = my_df[~my_df["sequence"].str.startswith("missing")]

        key_cols_array = ["srrID", "assembly"]
        t_df = my_df.groupby(key_cols_array, group_keys=False).apply(add_missing_freqs).reset_index(drop=True)

        def get_levels_and_freq_per_group(x):
            min_level = x['level'].min()
            x = x.loc[x['level'] == min_level]
            my_argmax = x['freq'].fillna(0).idxmax()
            return x.loc[my_argmax:my_argmax]

        key_cols_array_seq = ["srrID", "assembly", "sequence"]
        t_df = t_df.groupby(key_cols_array_seq, group_keys=False).apply(get_levels_and_freq_per_group)
        t_df = t_df.reset_index(drop=True)

        def get_levels_per_sample(x):
            min_level = x['level'].min()
            x = x.loc[x['level'] == min_level]
            my_argmax = x['freq'].fillna(0).idxmax()
            return x.loc[my_argmax:my_argmax]

        t_df = t_df.groupby(key_cols_array, group_keys=False).apply(get_levels_per_sample).reset_index(drop=True)
        t_df.to_csv(cached_fl)

    sns.set(style="whitegrid", context="paper", font_scale=1.5)
    g = sns.catplot(data=t_df, y='freq', x='level', kind='violin', cut=0,
                    inner='stick', col='assembly', col_wrap=2, palette='Set2')

    g.set_axis_labels("Severity level", "Abundance (%)")

    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True


def draw_drm_max_severity_per_sample(my_df, plot_name, how='min'):
    global output_path

    cached_fl = "E:\\app\\HexaHedron_scripts\\pan_HIV_analysis\\results\\lanl_runs_max_levels_and_frequencies_per_sample.csv"
    t_df=""
    if os.path.exists(cached_fl):
        t_df = pd.read_csv(cached_fl)
    else:
        def add_missing_freqs(x):
            no_mutations_freq = 100 - x.groupby('sequence')['freq'].first().fillna(0).sum()
            if no_mutations_freq == 100:
                x.iat[0, 11] = 0
                x.iat[0, 89] = 100
            elif no_mutations_freq > 0:
                last_row = x.iloc[-1]
                last_row['freq'] = no_mutations_freq
                last_row['level'] = 0
                last_row['sequence'] = '-'
                x = x.append(last_row)
            return x

        my_df = my_df.dropna(subset=["sequence"])
        my_df = my_df[~my_df["sequence"].str.startswith("missing")]

        key_cols_array = ["srrID", "assembly"]
        t_df = my_df.groupby(key_cols_array, group_keys=False).apply(add_missing_freqs).reset_index(drop=True)

        def get_levels_and_freq_per_group(x):
            max_level = x['level'].max()
            x = x.loc[x['level'] == max_level]
            my_argmax = x['freq'].fillna(0).idxmax()
            return x.loc[my_argmax:my_argmax]

        key_cols_array_seq = ["srrID", "assembly", "sequence"]
        t_df = t_df.groupby(key_cols_array_seq, group_keys=False).apply(get_levels_and_freq_per_group)
        t_df = t_df.reset_index(drop=True)

        def get_levels_per_sample(x):
            max_level = x['level'].max()
            x = x.loc[x['level'] == max_level]
            my_argmax = x['freq'].fillna(0).idxmax()
            return x.loc[my_argmax:my_argmax]

        t_df = t_df.groupby(key_cols_array, group_keys=False).apply(get_levels_per_sample).reset_index(drop=True)
        t_df.to_csv(cached_fl)

    sns.set(style="whitegrid", context="paper", font_scale=1.5)
    g = sns.catplot(data=t_df, y='freq', x='level', kind='violin', cut=0,
                    inner='stick', col='assembly', col_wrap=2, palette='Set2')

    g.set_axis_labels("Severity level", "Abundance (%)")

    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True

def draw_drm_severity_and_diversity_per_sample(my_df, plot_name):
    global output_path

    cached_fl = "E:\\app\\HexaHedron_scripts\\pan_HIV_analysis\\results\\lanl_runs_max_levels_and_diverisity_per_sample.csv"
    t_df=""
    if os.path.exists(cached_fl):
        t_df = pd.read_csv(cached_fl)
    else:
        def add_missing_freqs(x):
            no_mutations_freq = 100 - x.groupby('sequence')['freq'].first().fillna(0).sum()
            if no_mutations_freq == 100:
                x.iat[0, 11] = 0
                x.iat[0, 32] = 0
                x.iat[0, 89] = 100
            elif no_mutations_freq > 0:
                last_row = x.iloc[-1]
                last_row['freq'] = no_mutations_freq
                last_row['level'] = 0
                last_row['sequence'] = '-'
                x = x.append(last_row)
            return x

        my_df = my_df.dropna(subset=["sequence"])
        my_df = my_df[~my_df["sequence"].str.startswith("missing")]

        key_cols_array = ["srrID", "assembly"]
        t_df = my_df.groupby(key_cols_array, group_keys=False).apply(add_missing_freqs).reset_index(drop=True)

        def get_levels_and_freq_per_group(x):
            max_level = x['level'].max()
            x = x.loc[x['level'] == max_level]
            my_argmax = x['freq'].fillna(0).idxmax()
            return x.loc[my_argmax:my_argmax]

        key_cols_array_seq = ["srrID", "assembly", "sequence"]
        t_df = t_df.groupby(key_cols_array_seq, group_keys=False).apply(get_levels_and_freq_per_group)
        t_df = t_df.reset_index(drop=True)

        def get_levels_per_sample(x):
            max_level = x['level'].max()
            x = x.loc[x['level'] == max_level]
            my_argmax = x['freq'].fillna(0).idxmax()
            return x.loc[my_argmax:my_argmax]

        t_df = t_df.groupby(key_cols_array, group_keys=False).apply(get_levels_per_sample).reset_index(drop=True)
        t_df.to_csv(cached_fl)

    t_df.rename(columns={'assembly': 'Assembly', 'level':'Severity level'}, inplace=True)
    sns.set(style="whitegrid", context="paper", font_scale=1.5)
    g = sns.relplot(x="sample diversity", y="freq", col="Severity level", s=50,col_wrap=3,
                    hue="Assembly", data=t_df, legend="full")
    plt.ylim(0, None)
    plt.xlim(0, None)

    # g = sns.relplot(x="sample diversity", y="freq", size="Severity level", sizes=(10, 100),
    #                 col="Assembly", data=t_df, legend="full")
    # g = sns.catplot(data=t_df, y='freq', x='level', kind='violin', cut=0,
    #                 inner='stick', col='assembly', col_wrap=2, palette='Set2')

    g.set_axis_labels("Sample diversity", "Abundance (%)")

    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True


def draw_drm_min_severity_and_diversity_per_sample(my_df, plot_name):
    global output_path

    cached_fl = "E:\\app\\HexaHedron_scripts\\pan_HIV_analysis\\results\\lanl_runs_min_levels_and_diverisity_per_sample.csv"
    t_df=""
    if os.path.exists(cached_fl):
        t_df = pd.read_csv(cached_fl)
    else:
        def add_missing_freqs(x):
            no_mutations_freq = 100 - x.groupby('sequence')['freq'].first().fillna(0).sum()
            if no_mutations_freq == 100:
                x.iat[0, 11] = 0
                x.iat[0, 32] = 0
                x.iat[0, 89] = 100
            elif no_mutations_freq > 0:
                last_row = x.iloc[-1]
                last_row['freq'] = no_mutations_freq
                last_row['level'] = 0
                last_row['sequence'] = '-'
                x = x.append(last_row)
            return x

        my_df = my_df.dropna(subset=["sequence"])
        my_df = my_df[~my_df["sequence"].str.startswith("missing")]

        key_cols_array = ["srrID", "assembly"]
        t_df = my_df.groupby(key_cols_array, group_keys=False).apply(add_missing_freqs).reset_index(drop=True)

        def get_levels_and_freq_per_group(x):
            min_level = x['level'].min()
            x = x.loc[x['level'] == min_level]
            my_argmax = x['freq'].fillna(0).idxmax()
            return x.loc[my_argmax:my_argmax]

        key_cols_array_seq = ["srrID", "assembly", "sequence"]
        t_df = t_df.groupby(key_cols_array_seq, group_keys=False).apply(get_levels_and_freq_per_group)
        t_df = t_df.reset_index(drop=True)

        def get_levels_per_sample(x):
            min_level = x['level'].min()
            x = x.loc[x['level'] == min_level]
            my_argmax = x['freq'].fillna(0).idxmax()
            return x.loc[my_argmax:my_argmax]

        t_df = t_df.groupby(key_cols_array, group_keys=False).apply(get_levels_per_sample).reset_index(drop=True)
        t_df.to_csv(cached_fl)

    t_df.rename(columns={'assembly': 'Assembly', 'level': 'Severity level'}, inplace=True)
    sns.set(style="whitegrid", context="paper", font_scale=1.5)
    g = sns.relplot(x="sample diversity", y="freq", col="Severity level", col_wrap=3, hue="Assembly",
                    data=t_df, legend="full")
    # g = sns.catplot(data=t_df, y='freq', x='level', kind='violin', cut=0,
    #                 inner='stick', col='assembly', col_wrap=2, palette='Set2')

    g.set_axis_labels("Sample diversity", "Abundance (%)")

    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True



def draw_drm_who_severity_per_sample_as_individual_sequences(my_df, plot_name):
    global output_path

    cached_fl = "E:\\app\\HexaHedron_scripts\\pan_HIV_analysis\\results\\lanl_runs_who_levels_and_frequencies_per_sample.csv"
    t_df=""
    if os.path.exists(cached_fl):
        t_df = pd.read_csv(cached_fl)
    else:
        def add_missing_freqs(x):
            no_mutations_freq = 100 - x.groupby('sequence')['freq'].first().fillna(0).sum()
            if no_mutations_freq == 100:
                x.iat[0, 11] = 0
                x.iat[0, 89] = 100
            elif no_mutations_freq > 0:
                last_row = x.iloc[-1]
                last_row['freq'] = no_mutations_freq
                last_row['level'] = 0
                last_row['sequence'] = '-'
                x = x.append(last_row)
            return x

        my_df = my_df.dropna(subset=["sequence"])
        my_df = my_df[~my_df["sequence"].str.startswith("missing")]

        key_cols_array = ["srrID", "assembly"]
        t_df = my_df.groupby(key_cols_array, group_keys=False).apply(add_missing_freqs).reset_index(drop=True)

        def get_levels_and_freq_per_group(x):
            global art_combination
            a=[]
            for drg in art_combination:
                if drg['name'] not in list(x['drug']):
                    last_row = x.iloc[-1]
                    last_row['drug'] = drg['name']
                    last_row['class'] = drg['class']
                    last_row['level'] = 0
                    last_row['score'] = 0
                    x = x.append(last_row)
                # else:
                #     a.append(x.loc[x['drug'] == drg['name']])

            return x

        key_cols_array_seq = ["srrID", "assembly", "sequence"]
        t_df = t_df.groupby(key_cols_array_seq, group_keys=False).apply(get_levels_and_freq_per_group)
        t_df = t_df.reset_index(drop=True)

        t_df.to_csv(cached_fl)

    related_drugs = list([x['name'] for x in art_combination])
    t_df = t_df.loc[t_df['drug'].isin(related_drugs)]
    # filter_size = len(related_drugs)
    # t_df = t_df.groupby(["srrID", "assembly"]).filter(lambda x: x['sequence'].count()>filter_size)
    # t_df = t_df.groupby(["srrID", "assembly"]).filter(lambda x: x['Mfmax'].max()>0.1)

    def get_min_levels_per_sequence(x):
        min_level = x['level'].min()
        x = x.loc[x['level'] == min_level]
        my_argmax = x['freq'].fillna(0).idxmax()
        return x.loc[my_argmax:my_argmax]

    key_cols_array_seq = ["srrID", "assembly", "sequence"]
    t_df = t_df.groupby(key_cols_array_seq, group_keys=False).apply(get_min_levels_per_sequence).reset_index(drop=True)

    def get_levels_per_sample(x):
        max_level = x['level'].max()
        x = x.loc[x['level'] == max_level]
        my_argmax = x['freq'].fillna(0).idxmax()
        return x.loc[my_argmax:my_argmax]

    key_cols_array = ["srrID", "assembly"]
    t_df = t_df.groupby(key_cols_array, group_keys=False).apply(get_levels_per_sample).reset_index(drop=True)

    sns.set(style="whitegrid", context="paper", font_scale=1.5)
    g = sns.catplot(data=t_df, y='freq', x='level', kind='violin', cut=0,
                    inner='stick', col='assembly', col_wrap=2, palette='Set2')

    g.set_axis_labels("Severity level", "Abundance (%)")

    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True


def draw_drm_who_severity_per_sample_as_quasispecies(my_df, plot_name):
    global output_path

    cached_fl = "E:\\app\\HexaHedron_scripts\\pan_HIV_analysis\\results\\lanl_runs_who_levels_and_frequencies_per_sample.csv"
    t_df=""
    if os.path.exists(cached_fl):
        t_df = pd.read_csv(cached_fl)
    else:
        def add_missing_freqs(x):
            no_mutations_freq = 100 - x.groupby('sequence')['freq'].first().fillna(0).sum()
            if no_mutations_freq == 100:
                x.iat[0, 11] = 0
                x.iat[0, 89] = 100
            elif no_mutations_freq > 0:
                last_row = x.iloc[-1]
                last_row['freq'] = no_mutations_freq
                last_row['level'] = 0
                last_row['sequence'] = '-'
                x = x.append(last_row)
            return x

        my_df = my_df.dropna(subset=["sequence"])
        my_df = my_df[~my_df["sequence"].str.startswith("missing")]

        key_cols_array = ["srrID", "assembly"]
        t_df = my_df.groupby(key_cols_array, group_keys=False).apply(add_missing_freqs).reset_index(drop=True)

        def get_levels_and_freq_per_group(x):
            a = []
            global art_combination
            for drg in art_combination:
                if drg['name'] not in list(drg['name']):
                    last_row = x.iloc[-1]
                    last_row['drug'] = drg['name']
                    last_row['class'] = drg['class']
                    last_row['level'] = 0
                    last_row['score'] = 0
                    a.append(last_row)
                else:
                    a.append(x.loc[x['drug'] == drg['name']])

            return a

        key_cols_array_seq = ["srrID", "assembly", "sequence"]
        t_df = t_df.groupby(key_cols_array_seq, group_keys=False).apply(get_levels_and_freq_per_group)
        t_df = t_df.reset_index(drop=True)

        t_df.to_csv(cached_fl)


    related_drugs = list([x['name'] for x in art_combination])
    t_df = t_df.loc[t_df['drug'].isin(related_drugs)]
    # filter_size = len(related_drugs)
    # t_df = t_df.groupby(["srrID", "assembly"]).filter(lambda x: x['sequence'].count()>filter_size)
    # t_df = t_df.groupby(["srrID", "assembly"]).filter(lambda x: x['Mfmax'].max() > 0.1)

    def get_max_levels_per_drug(x):
        max_level = x['level'].max()
        x = x.loc[x['level'] == max_level]
        my_argmax = x['freq'].fillna(0).idxmax()
        return x.loc[my_argmax:my_argmax]

    key_cols_array_seq = ["srrID", "assembly", "drug"]
    t_df = t_df.groupby(key_cols_array_seq, group_keys=False).apply(get_max_levels_per_drug).reset_index(drop=True)

    def get_levels_per_sample(x):
        min_level = x['level'].min()
        x = x.loc[x['level'] == min_level]
        my_argmax = x['freq'].fillna(0).idxmax()
        return x.loc[my_argmax:my_argmax]

    key_cols_array = ["srrID", "assembly"]
    t_df = t_df.groupby(key_cols_array, group_keys=False).apply(get_levels_per_sample).reset_index(drop=True)

    sns.set(style="whitegrid", context="paper", font_scale=1.5)
    g = sns.catplot(data=t_df, y='freq', x='level', kind='violin', cut=0,
                    inner='stick', col='assembly', col_wrap=2, palette='Set2')

    g.set_axis_labels("Severity level", "Abundance (%)")

    plt.savefig(os.path.join(output_path, plot_name + ".png"))

    print "Done with " + plot_name
    plt.clf()
    return True


def usage():
    print 'parse_pan_hive_results.py -i <input_file> -o <output_path> -t <HIV phylogenetic tree in newick format>'


def convert_json_column(data):
    return json.loads(data)


def main(argv):

    try:
        opts, args = getopt.getopt(argv, 'i:o:t:h', ['input=', 'output=', 'tree=', 'help'])
    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit(2)
    if len(opts) < 2:
        usage()
        sys.exit(2)

    global output_path
    global tree_path
    my_input = ""
    
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt in ('-i', '--input'):
            my_input = arg
        elif opt in ('-o', '--output'):
            output_path = arg
        elif opt in ('-t', '--tree'):
            tree_path = arg
        else:
            usage()
            sys.exit()
    
    if not os.path.exists(my_input):
        print "Input file doesn't exist"
        sys.exit()
        
    if not os.path.exists(output_path):
        try :
            os.makedirs(output_path)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                print "Missing or corrupted output dir"
                raise

    my_df = pd.read_csv(my_input, escapechar="\\")
    # # Fixed it with new versions of the file
    my_df = fix_consensus_names(my_df, 'Final Reference names')
    # my_df.drop(my_df.columns[36], axis=1, inplace=True)

    # draw_runs_per_technology(my_df)
    # draw_runs_per_subtype(my_df)
    # draw_subtypes_per_technology(my_df)
    # draw_unaligned_per_technology(my_df)
    # draw_nucleotide_diversity_per_technology(my_df)
    # draw_nucleotide_diversity_per_subtype(my_df)
    # draw_div_stats(my_df, 'Platform', 'Platform', 'per_platform')
    # draw_div_stats(my_df, 'Final Reference names', 'Forms', 'per_subtype')
    # draw_subtype_coocurrence(my_df)

    # draw_drm_frequency_bins(my_df, "frequency_histogram")
    # draw_drm_frequency_bins(my_df, "frequency_histogram_up_to_10", my_bins=range(0, 11, 1))
    # draw_drm_frequency_bins(my_df, "frequency_histogram_percentage", percentage=True)
    # draw_drm_frequency_bins(my_df, "frequency_histogram_percentage_up_to_10", my_bins=range(0, 11, 1), percentage=True)
    # draw_drm_frequency_bins_per_class(my_df, "frequency_histogram_per_class")
    # draw_drm_frequency_bins_per_class(my_df, "frequency_histogram_per_class_percentage_", percentage=True)
    # draw_drm_frequency_bins_per_platform(my_df, "frequency_histogram_per_platform")
    # draw_drm_frequency_bins_per_platform(my_df, "frequency_histogram_per_platform_percentage", percentage=True)
    # draw_drm_frequency_abundance_scatter_on_level_and_class(
    #     my_df, 'frequency_abundance_scatter_on_min_lvl_and_assembly', how='min')
    # draw_drm_frequency_abundance_scatter_on_level_and_class(
    #     my_df, 'frequency_abundance_scatter_on_max_lvl_and_assembly', how='max')
    # draw_drm_frequency_abundance_scatter_on_level_and_class(
    #     my_df, 'frequency_abundance_scatter_on_mean_lvl_and_assembly', how='mean')
    # draw_drm_frequency_abundance_scatter_on_level_and_subtype(
    #     my_df, 'frequency_abundance_scatter_on_min_lvl_and_subtype', how='min')
    # draw_drm_frequency_abundance_scatter_on_level_and_subtype(
    #     my_df, 'frequency_abundance_scatter_on_max_lvl_and_subtype', how='max')
    # draw_drm_frequency_abundance_scatter_on_level_and_subtype(
    #     my_df, 'frequency_abundance_scatter_on_mean_lvl_and_subtype', how='mean')

    # my_df.dropna(subset=['mutations'],inplace=True)
    # my_df = remove_scores_from_mutations_df(my_df)
    # my_df['length of mutations'] = my_df['mutations'].apply(lambda x: len(json.loads(x)))
    # draw_drm_severity_and_number_of_mutations_scatter_on_synergistic(
    #     my_df, 'cardinality_and_number_mutations_on_severity_per_sequence', how='min')
    # draw_drm_severity_and_cardinatlity_of_mutations(
    #     my_df, 'cardinality_of_mutations_and_severity', how='min')
    # draw_drm_severity_per_sample(
    #     my_df, 'severity_per_sample', how='min')
    # draw_drm_max_severity_per_sample(
    #     my_df, 'max_severity_per_sample', how='min')
    # draw_drm_severity_and_diversity_per_sample(
    #     my_df, 'max_severity_and_diversity_per_sample')
    # draw_drm_min_severity_and_diversity_per_sample(
    #     my_df, 'min_severity_and_diversity_per_sample')

    global art_combination
    art_combination = [{'class': 'NNRTI', 'name': 'TDF'},
                       {'class': 'NNRTI', 'name': '3TC'},
                       {'class': 'NNRTI', 'name': 'EFV'}]

    # draw_drm_who_severity_per_sample_as_individual_sequences(my_df,'who_resistance_individual_all')
    # draw_drm_who_severity_per_sample_as_quasispecies(my_df,'who_resistance_quasi_all')

    draw_arm_coocurrence(my_df)

    # cur_mut_string = my_df["mutations"]
    # cur_mut_string = replace('\'','"')
    # cur_muts = json.loads(cur_mut)

def draw_arm_coocurrence(my_df):
    global output_path

    func_name = sys._getframe().f_code.co_name
    plot_name = func_name.replace("draw_", "")

    global global_coocurs
    global_coocurs = {}
    global local_coocurs
    local_coocurs = {}
    global sample_local_coocurs
    sample_local_coocurs = {}
    global sample_global_coocurs
    sample_global_coocurs = {}
    global cur_coocurs
    cur_coocurs = {}
    global cur_srr_global_coocurs
    cur_srr_global_coocurs = {}
    global cur_srr_local_coocurs
    cur_srr_local_coocurs = {}
    global prev_srr_id
    prev_srr_id = ""
    global prev_sequence
    prev_sequence = ""
    global prev_assembly
    prev_assembly = ""

    def sum_dictionary(cooccurs, my_dict):
        if len(my_dict) <= 1:
            return
        for c_ref in my_dict:
            for ic_ref in my_dict[c_ref]:
                if c_ref not in cooccurs:
                    cooccurs[c_ref] = {}
                if ic_ref not in cooccurs[c_ref]:
                    cooccurs[c_ref][ic_ref] = [0, 0, 0, 0, 0, 0]
                for nn in range(6):
                    if my_dict[c_ref][ic_ref][nn] > 0:
                        cooccurs[c_ref][ic_ref][nn] += 1

    def flush_dictionary(cooccurs, my_dict, collapse=False):
        if len(my_dict) <= 1:
            return
        for c_ref in my_dict:
            for ic_ref in my_dict:
                if ic_ref == c_ref:
                    continue
                if c_ref not in cooccurs:
                    cooccurs[c_ref] = {}
                if ic_ref not in cooccurs:
                    cooccurs[ic_ref] = {}
                if ic_ref not in cooccurs[c_ref]:
                    cooccurs[c_ref][ic_ref] = [0, 0, 0, 0, 0, 0]
                if c_ref not in cooccurs[ic_ref]:
                    cooccurs[ic_ref][c_ref] = [0, 0, 0, 0, 0, 0]
                for nn in range(6):
                    for inn in range(5, nn-1, -1):
                        if my_dict[c_ref][inn] > 0 and my_dict[ic_ref][inn] > 0:
                            if collapse:
                                cooccurs[c_ref][ic_ref][nn] = 1
                                cooccurs[ic_ref][c_ref][nn] = 1
                            else:
                                cooccurs[c_ref][ic_ref][nn] += 1
                                cooccurs[ic_ref][c_ref][nn] += 1
                            break

    def get_coocurrences(my_row):
        global prev_srr_id
        global prev_sequence
        global prev_assembly
        global cur_coocurs
        global cur_srr_local_coocurs
        global cur_srr_global_coocurs

        my_mut_str = my_row['mutations']
        my_sequence = my_row['sequence']
        my_srr_id = my_row['srrID']
        my_assembly = my_row['assembly']
        if prev_srr_id is "":
            prev_sequence = my_sequence
            prev_srr_id = my_srr_id
            prev_assembly = my_assembly
            cur_coocurs = {}

        if (prev_sequence != my_sequence) or (my_assembly != prev_assembly) or ((my_srr_id != prev_srr_id)):
            if prev_assembly == "global":
                global global_coocurs
                flush_dictionary(global_coocurs, cur_coocurs)
                flush_dictionary(cur_srr_global_coocurs, cur_coocurs, True)
            else:
                global local_coocurs
                flush_dictionary(local_coocurs, cur_coocurs)
                flush_dictionary(cur_srr_local_coocurs, cur_coocurs, True)
                
            if (my_srr_id != prev_srr_id):
                global sample_global_coocurs
                sum_dictionary(sample_global_coocurs, cur_srr_global_coocurs)
                cur_srr_global_coocurs = {}
                global sample_local_coocurs
                sum_dictionary(sample_local_coocurs, cur_srr_local_coocurs)
                cur_srr_local_coocurs = {}
                prev_srr_id = my_srr_id
            if (prev_assembly != my_assembly):
                prev_assembly = my_assembly
            if (prev_sequence != my_sequence):
                prev_sequence = my_sequence
            cur_coocurs = {}

        if pd.isna(my_mut_str) is not True:
            my_level = int(my_row['level'])

            all_muts = json.loads(my_mut_str)
            for muts in all_muts:
                for my_mut in muts["mutations"]:
                    if my_mut not in cur_coocurs:
                        cur_coocurs[my_mut] = [0,0,0,0,0,0]
                    cur_coocurs[my_mut][my_level] = 1

    [get_coocurrences(m_row) for m_index, m_row in my_df.iterrows()]

    def convert_to_array(idx_dict,lvl,assembly):
        m_df = pd.DataFrame(idx_dict)
        for m_key1 in m_df.keys():
            for m_key2 in m_df.keys():
                if m_key1 == m_key2:
                    continue
                if not m_df.isna()[m_key1][m_key2]:
                    m_df[m_key2][m_key1] = np.nan
        m_df = m_df.reset_index()
        m_df = m_df.melt(id_vars=['index']).dropna()
        for i, row in m_df.iterrows():
            m_df.at[i, 'value'] = row['value'][lvl]
        m_df = m_df.sort_values(by='value', ascending=False)
        m_df = m_df[m_df.value != 0]
        global output_path
        m_df.to_csv(os.path.join(output_path, "coocurs_table_"+assembly + "_lvl_"+str(lvl)+".csv"))
        return m_df.iloc[10].value

    def get_output_from_dictionaries(idx_dict, assembly):
        for nn in range(6):
            sbtp_tree_idx = {}
            thrsld = convert_to_array(idx_dict,nn,assembly)
            for key_1 in idx_dict.keys():
                sbtp_tree_idx[key_1] = {}
                for key_2 in idx_dict.keys():
                    if key_2 in idx_dict[key_1] and idx_dict[key_1][key_2][nn] >= thrsld:
                        sbtp_tree_idx[key_1][key_2] = int(idx_dict[key_1][key_2][nn])
                    # else:
                    #     sbtp_tree_idx[key_1][key_2] = 0
            mcl_df = pd.DataFrame(data=sbtp_tree_idx)
            mcl_df = mcl_df.dropna(axis='columns', how='all')
            mcl_df = mcl_df.fillna(int(0))
            # mcl_df = mcl_df.reindex(sbtp_tre e_idx, columns=sbtp_tree_idx)

            mcl_df = mcl_df.sort_index(axis=0).sort_index(axis=1)
            mask = np.zeros_like(mcl_df)
            mask[np.triu_indices_from(mask)] = True
            g = sns.heatmap(mcl_df, annot=True, linewidths=.2, fmt='.0f', mask=mask, cmap="Reds",vmin=0)
            for v in g.get_xticklabels():
                v.set_rotation(60)
            plt.tight_layout()
            plt.savefig(os.path.join(output_path, plot_name + "_" + assembly + "_lvl_"+str(nn)+".png"))
            plt.clf()

    # get_output_from_dictionaries(global_coocurs,"global")
    # get_output_from_dictionaries(local_coocurs,"local")
    get_output_from_dictionaries(sample_global_coocurs,"sample_global")
    get_output_from_dictionaries(sample_local_coocurs,"sample_local")

    print "Done with " + func_name

    return True

def remove_scores_from_cell(my_x):
    j = json.loads(my_x)
    res=[]
    for m in j:
        if 'score' in m:
            del m['score']
        res.append(m)
    return res


def remove_scores_from_mutations_df(my_df):
    my_df['mutations'].apply( lambda x: remove_scores_from_cell(x) )
    return my_df


def draw_div_stats(my_df, my_y_group, my_y_group_title, graph_sub_title):
    #globals
    # draw_diversity_per_subtype(my_df, 'Haplotypes', '# of haplotypes', 'Haplotype_'+graph_sub_title,
    #                            which_assembly="global", y_group=my_y_group, y_group_title=my_y_group_title)
    draw_diversity_per_subtype(my_df, 'Mfmin', 'Mfmin', 'mfmin_'+graph_sub_title, which_assembly="global",
                               y_group = my_y_group, y_group_title = my_y_group_title)
    # draw_diversity_per_subtype(my_df, 'Mfe', 'Mfe', 'mfe_'+graph_sub_title, which_assembly="global",
    #                            y_group=my_y_group, y_group_title=my_y_group_title)
    # draw_diversity_per_subtype(my_df, 'Mfmax', 'Mfmax', 'mfmax_'+graph_sub_title, which_assembly="global",
    #                            y_group=my_y_group, y_group_title=my_y_group_title)
    # draw_diversity_per_subtype(my_df, 'Shannon', 'Shannon\'s entropy', 'shannon_'+graph_sub_title, which_assembly="global",
    #                            y_group=my_y_group, y_group_title=my_y_group_title)
    # draw_diversity_per_subtype(my_df, 'Shannon log(Reads) normalized', 'Shannon\'s entropy (read count normalized)',
    #                            'shannon_read_norm_'+graph_sub_title, which_assembly="global",
    #                            y_group=my_y_group, y_group_title=my_y_group_title)
    # draw_diversity_per_subtype(my_df, 'Shannon log(Haplotypoes) normalized', 'Shannon\'s entropy (haplotype count normalized)',
    #                            'shannon_heplotype_norm_'+graph_sub_title, which_assembly="global",
    #                            y_group=my_y_group, y_group_title=my_y_group_title)
    # draw_diversity_per_subtype(my_df, 'Shannon bias corrected', 'Shannon\'s entropy bias corrected',
    #                            'shannon_bias_corrected_'+graph_sub_title, which_assembly="global",
    #                            y_group=my_y_group, y_group_title=my_y_group_title)
    # draw_diversity_per_subtype(my_df, 'population nucleotide diversity', 'population nucleotide diversity',
    #                            'nucleotide_diversity_'+graph_sub_title, which_assembly="global",
    #                            y_group=my_y_group, y_group_title=my_y_group_title)
    # draw_diversity_per_subtype(my_df, 'FAD', 'FAD', 'fad_per_'+graph_sub_title, which_assembly="global",
    #                            y_group=my_y_group, y_group_title=my_y_group_title)


if __name__ == '__main__':
    main(sys.argv[1:])


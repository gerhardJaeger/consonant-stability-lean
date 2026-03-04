import argparse
import pandas as pd
import itertools
import os
import sys
import numpy as np
from pathlib import Path
from pybtex.database import parse_file
import shutil
from tqdm import tqdm
from copy import deepcopy

CAT = "category"

# airstream mechanism variable names
VOIC = "periodicGlottalSource"
ASPI = "spreadGlottis"
GLOT = "constrictedGlottis"
EJEC = "raisedLarynxEjective"
IMPL = "loweredLarynxImplosive"
#CLIC = "click"
LONG = "long"
FORT = "fortis"
NASA = "nasal" # for +,- as in prenasalized stops
SONO = "sonorant" # for +,- as in prenasalized stops
AIRS = [VOIC, ASPI, GLOT, EJEC, IMPL, SONO, NASA, FORT]#, CLIC, LONG]

# oral cavity variable names
LABI = "labial"
LADE = "labiodental"
CORO = "coronal"
ANTE = "anterior"
DIST = "distributed"
DORS = "dorsal"
HIGH = "high"
LOW  = "low"
FRON = "front"
BACK = "back"
EPIL = "epilaryngealSource"
ROUN = "round"
LAT  = "lateral"
STRI = "strident"
RETR = "retractedTongueRoot"
#DELR = "delayedRelease"
PLOS = [LABI, LADE, CORO, LAT, ANTE, DIST, DORS, HIGH, LOW, FRON, BACK, EPIL, ROUN, RETR]#, CLIC
FRIC = [LABI, LADE, CORO, LAT, STRI, ANTE, DIST, DORS, HIGH, LOW, FRON, BACK, EPIL, ROUN, RETR]#, CLIC, DELR
AFFR = [LABI, LADE, CORO, LAT, STRI, ANTE, DIST, DORS, HIGH, LOW, FRON, BACK, EPIL, ROUN, RETR]#, CLIC, DELR
# note that CLIC nees to be in the oral parameters so that CLICK_TRANSFORMATIONS can be triggered

# special characters
LBOX = chr(827)
LBAR = chr(826)
LCRO = chr(799)
DEVS = [chr(805), chr(778)]
LCOR = chr(841)
LCORS = " " + LCOR
DEVSS = [" " + x for x in DEVS]
DEV_STR = "DEVOICED"
ORALCHR = [" " + LBOX, " " + LBAR, " " + LCRO]
SP_CHRS = ORALCHR + [DEV_STR] + [LCORS]


def main(phoible_loc, transformation_loc, inventories_loc, segments_loc, languages_loc, output_loc, full_headers, suppress_charts):
    global BASIC_TRANSFORMATIONS
    global AIRSTREAM_TRANSFORMATIONS
    BASIC_TRANSFORMATIONS, AIRSTREAM_TRANSFORMATIONS = read_transformation_list(transformation_loc)
        
    data = pd.read_csv(phoible_loc + "/data/phoible.csv", keep_default_na=False, low_memory=False)
    inventories = pd.read_csv(inventories_loc, keep_default_na=False, low_memory=False)
    segments = pd.read_csv(segments_loc, keep_default_na=False, low_memory=False)
    
    # create output_loc if it doesn't exist
    if not suppress_charts:
        Path(output_loc + "charts_full/").mkdir(parents=True, exist_ok=True)
    else:
        Path(output_loc).mkdir(parents=True, exist_ok=True)
    
    # filter data to only include the inventories we are using
    data = data[data.InventoryID.isin(inventories.Inventory)]
    
    # add a functional column saying if something is an oral plosive, affricate, fricative
    data.loc[:, CAT] = data.apply(categorize_phoneme, axis=1)
    segments.loc[:, CAT] = segments.apply(categorize_phoneme, axis=1)
    
    all_plosives = data[(data.category == "plosive")]
    all_fricatives = data[(data.category == "fricative")]
    all_affricates = data[(data.category == "affricate")]
    plosive_segments = segments[(segments.category == "plosive")]
    fricative_segments = segments[(segments.category == "fricative")]
    affricate_segments = segments[(segments.category == "affricate")]
    
    #  remove all glottally articulated stops ʔ, ʔʲ, ʔʷ, etc
    all_plosives = all_plosives[~((all_plosives.coronal == "-") & (all_plosives.constrictedGlottis == "+") & 
                                  (all_plosives.labial != "+") & (all_plosives.delayedRelease == "-") & (all_plosives.dorsal != "+"))]
    plosive_segments = plosive_segments[~((plosive_segments.coronal == "-") & (plosive_segments.constrictedGlottis == "+") & 
                                  (plosive_segments.labial != "+") & (plosive_segments.delayedRelease == "-") & (plosive_segments.dorsal != "+"))]
    
    if languages_loc == "None":
        conv_list = None
        glot_set = set(data.Glottocode) - set(["NA"])
    else:
        conv_list = pd.read_csv(languages_loc)
        glot_set = set(pd.merge(conv_list.glottocode, data.Glottocode, left_on="glottocode", right_on="Glottocode").glottocode)
    
    glot_to_series = dict()
    for glot in tqdm(sorted(glot_set), file=sys.stdout):
        # print(glot, file=sys.stderr)
        # generate plosives chart
        plosives = all_plosives[(all_plosives.InventoryID == 
                        inventories.Inventory[inventories.Glottocode == glot].values[0])]
        plosives_chart, am_plosives = calculate_series(plosives, AIRS, PLOS, plosive_segments, full_headers)
        # generate fricatives chart
        fricatives = all_fricatives[(all_fricatives.InventoryID == 
                        inventories.Inventory[inventories.Glottocode == glot].values[0])]
        fricatives_chart, am_fricatives = calculate_series(fricatives, AIRS, FRIC, fricative_segments, full_headers)
        # generate affricates chart
        affricates = all_affricates[(all_affricates.InventoryID == 
                        inventories.Inventory[inventories.Glottocode == glot].values[0])]
        affricates_chart, am_affricates = calculate_series(affricates, AIRS, FRIC, affricate_segments, full_headers)
        # generate a combined chart
        obstruents_chart = generate_obstruent_chart(plosives, affricates, fricatives,
                                                    am_plosives, am_affricates, am_fricatives, 
                                                    plosive_segments, affricate_segments, fricative_segments,
                                                    full_headers)
        
        # do some error checking
        plosive_num = count_segments(plosives_chart)
        fricative_num = count_segments(fricatives_chart)
        affricate_num = count_segments(affricates_chart)
        if plosive_num != plosives.shape[0]:
            print(glot + ": Number of plosives in chart does not match number in source. Are there indistinguishable segments?", file=sys.stderr)
        if fricative_num != fricatives.shape[0]:
            print(glot + ": Number of fricatives in chart does not match number in source. Are there indistinguishable segments?", file=sys.stderr)
        if affricate_num != affricates.shape[0]:
            print(glot + ": Number of affricates in chart does not match number in source. Are there indistinguishable segments?", file=sys.stderr)
        expected_count = plosive_num + fricative_num + affricate_num
        if count_segments(obstruents_chart) != expected_count:
            if count_segments(obstruents_chart) < expected_count:
                print(glot + ": ERROR: Lost some segments along the way.", file=sys.stderr)
            else:
                print(glot + "ERROR: Duplication of segments detected.", file=sys.stderr)
        
        if not suppress_charts:
            write_chart(plosives_chart, output_loc + "charts_full/" + glot, "plosives", True)
            write_chart(fricatives_chart, output_loc + "charts_full/" + glot, "fricatives", True)
            write_chart(affricates_chart, output_loc + "charts_full/" + glot, "affricates", True)
            write_chart(obstruents_chart, output_loc + "charts_full/" + glot, "obstruents", True)
       
        # write series fullness to glot_to_series dict
        add_to_series_dict(glot_to_series, glot, plosives_chart, fricatives_chart, affricates_chart, obstruents_chart)

    # write the dictionary out to file as a table
    glot_to_series_df = pd.DataFrame.from_dict(glot_to_series, orient="index",
                                               columns=["plosive_series", "plosive_fullness", "plosive_markedness_fullness",
                                                        "fricative_series", "fricative_fullness", "fricative_markedness_fullness",
                                                        "affricate_series", "affricate_fullness", "affricate_markedness_fullness",
                                                        "series", "series_fullness", "series_markedness_fullness",
                                                        "plosive_places", "plosive_places_fullness", "plosive_places_markedness_fullness",
                                                        "fricative_places", "fricative_places_fullness", "fricative_places_markedness_fullness",
                                                        "affricate_places", "affricate_places_fullness", "affricate_places_markedness_fullness",
                                                        "all_places","all_places_fullness", "all_places_markedness_fullness"])
    glot_to_series_df = glot_to_series_df.sort_index()
    glot_to_series_df.to_csv(output_loc + "/series_counts.tsv", index_label="glottocode", sep="\t")


# Function that returns whether a given row is a plosive, fricative, affricate, or other
def categorize_phoneme(row):
    if row.SegmentClass == "consonant" and row.consonantal == "+" and row.sonorant != "+" and \
       row.continuant == "-" and row.approximant == "-" and row.delayedRelease == "-" and row.syllabic == "-" and \
       row.long != "+" and row.long != "-,+" and row.click == "-":
       return "plosive"
    elif row.SegmentClass == "consonant" and row.sonorant != "+" and row.approximant == "-" and \
        (row.continuant == "+" and row.consonantal == "+" and row.syllabic == "-" and \
        row.long != "+" and row.long != "-,+" and row.nasal != "+" and row.click == "-") or (row.continuant == "-,+" and row.consonantal == "+" and row.syllabic == "-" and \
        row.long != "+" and row.long != "-,+" and row.nasal == "+,-" and row.sonorant == "+,-" and row.click == "-"): # nasal fricatives like β̃ I am assuming are primarily nasals, not fricatives
        return "fricative"
    elif row.SegmentClass == "consonant" and row.consonantal == "+" and row.approximant == "-" and row.sonorant != "+" and \
        row.continuant != "+" and row.syllabic == "-" and row.delayedRelease != "-" and \
        row.long != "+" and row.long != "-,+" and row.long != "+,-" and row.click == "-":
        return "affricate"
    else:
        return "other"

# Calculates an entire series (plosives, affricates, fricatives)
# Input: data         the data to be 
#        am_param_set the possible set of parameters for airstream mechanisms
#        op_param_set the possible set of parameters for oral closures
#        segments     the list of possible segments in PHOIBLE
#        full_headers whether to return the full headers or to truncate
# Output: a tuple of
#           chart     the fully articulated chart
#           am        the airstream parameters detected for this data
def calculate_series(data, am_param_set, op_param_set, segments, full_headers):
    am = get_variation_set(am_param_set, data)
    am = add_unicode_airstream_params(data, am)
    op = get_variation_set(op_param_set, data)
    op, op_short = collapse_oral_params(data, am, op, segments)
    am, am_short = collapse_airstream_mechanisms(data, am, op, segments)
    chart = generate_chart(data, am, op, segments)
    if not full_headers:
        chart = replace_headers(chart, am_short, op_short)
    return chart, am


# Adds elements to the glot_to_series dict
# Input: glot_to_series   the dict keeping track of all variables across all glottocodes
#        glot             the current glottocode
#        obstruents_chart an obstruent chart for the current glot
# Output: None, modifies the dict
def add_to_series_dict(glot_to_series, glot, plosives_chart, fricatives_chart, affricates_chart, obstruents_chart):
    plosive_series = len(plosives_chart) - 1
    plosive_fullness = 0
    plosive_markedness_fullness = 0
    for row in plosives_chart[1:-2]:
        plosive_fullness += convert_to_float(row[-2])
        plosive_markedness_fullness += convert_to_float(row[-1])
    plosive_places = len(plosives_chart[0][1:-2])
    plosive_places_fullness = sum([convert_to_float(x) for x in plosives_chart[-2][1:]])
    plosive_places_markedness_fullness = sum([convert_to_float(x) for x in plosives_chart[-1][1:]])
    fricative_series = len(fricatives_chart) - 1 if len(fricatives_chart) > 1 else 0
    fricative_fullness = 0
    fricative_markedness_fullness = 0
    for row in fricatives_chart[1:-2]:
        fricative_fullness += convert_to_float(row[-2])
        fricative_markedness_fullness += convert_to_float(row[-1])
    fricative_places = len(fricatives_chart[0][1:-2]) if len(fricatives_chart) > 1 else 0
    fricative_places_fullness = sum([convert_to_float(x) for x in fricatives_chart[-2][1:]]) if len(fricatives_chart) > 1 else 0
    fricative_places_markedness_fullness = sum([convert_to_float(x) for x in fricatives_chart[-1][1:]]) if len(fricatives_chart) > 1 else 0
    affricate_series = len(affricates_chart) - 1 if len(affricates_chart) > 1 else 0
    affricate_fullness = 0
    affricate_markedness_fullness = 0
    for row in affricates_chart[1:-2]:
        affricate_fullness += convert_to_float(row[-2])
        affricate_markedness_fullness += convert_to_float(row[-1]) 
    affricate_places = len(affricates_chart[0][1:-2]) if len(affricates_chart) > 1 else 0
    affricate_places_fullness = sum([convert_to_float(x) for x in affricates_chart[-2][1:]]) if len(affricates_chart) > 1 else 0
    affricate_places_markedness_fullness = sum([convert_to_float(x) for x in affricates_chart[-1][1:]]) if len(affricates_chart) > 1 else 0
    series = len(obstruents_chart) - 1
    series_fullness = 0
    series_markedness_fullness = 0
    for row in obstruents_chart[1:-2]:
        series_fullness += convert_to_float(row[-2])
        series_markedness_fullness += convert_to_float(row[-1])
    all_places = len(obstruents_chart[0][1:-2])
    all_places_fullness = sum([convert_to_float(x) for x in obstruents_chart[-2][1:]])
    all_places_markedness_fullness = sum([convert_to_float(x) for x in obstruents_chart[-1][1:]])
    glot_to_series[glot] = [plosive_series, plosive_fullness, plosive_markedness_fullness,
                            fricative_series, fricative_fullness, fricative_markedness_fullness,
                            affricate_series, affricate_fullness, affricate_markedness_fullness,
                            series, series_fullness, series_markedness_fullness,
                            plosive_places, plosive_places_fullness, plosive_places_markedness_fullness,
                            fricative_places, fricative_places_fullness, fricative_places_markedness_fullness,
                            affricate_places, affricate_places_fullness, affricate_places_markedness_fullness,
                            all_places, all_places_fullness, all_places_markedness_fullness]


# Reads a transformation list from file
# Input: transformation_loc the location of the transformation list
# Output: a list of dictionary tuples, in the format of:
#         [ ({pre-transformation properties}, {post-transformation properties}),
#           ({pre-transformation properties}, {post-transformation properties}), ... ]
def read_transformation_list(transformation_loc):
    oral_list = []
    airstream_list = []
    with open(transformation_loc) as f:
        for l in f.readlines()[1:]:
            l_spl = l.split("\t")
            pre_t = dict()
            post_t = dict()
            for kv in l_spl[1].split(";"):
                pre_t[kv.split(":")[0].strip()] = kv.split(":")[1].strip()
            for kv in l_spl[2].split(";"):
                post_t[kv.split(":")[0].strip()] = kv.split(":")[1].strip()
            if l_spl[0] == "oral":
                oral_list.append((pre_t, post_t))
            elif l_spl[0] == "airstream":
                airstream_list.append((pre_t, post_t))
    # add LCOR to airstream
    airstream_list.append(({LCORS: "+"}, {LCORS: "-"}))
    airstream_list.append(({DEV_STR: "+"}, {DEV_STR: "-"}))
    return oral_list, airstream_list


# Generates the combined obstruent chart of all plosives, affricates, and fricatives
# Input: plosives        a data frame of all the plosives in the language
#        affricates      a data frame of all the affricates in the language
#        fricatives      a data frame of all the fricatives in the language
#        am_plosives     a list of all the airstream mechanisms of the plosives
#        am_affricates   a list of all the airstream mechanisms of the plosives
#        am_fricatives   a list of all the airstream mechanisms of the plosives
#        all_plosives    a data frame of all plosives
#        all_affricates  a data frame of all affricates
#        all_fricatives  a data frame of all fricatives
# Output: a chart (list of lists) that combines plosives, affricates, and fricatives
def generate_obstruent_chart(plosives, affricates, fricatives, am_plosives, am_affricates, am_fricatives, 
                             plosive_segments, affricate_segments, fricative_segments, 
                             full_headers):
    if affricates.shape[0] == 0 and fricatives.shape[0] == 0: # only plosives
        obstruents = plosives
        ops_plosives = get_variation_set(PLOS, obstruents)
        ops_plosives, ops_plosives_short = collapse_oral_params(obstruents, am_plosives, ops_plosives, plosive_segments)
        for am in am_plosives:
            for am_tup in am:
                am_tup[CAT] = "plosive"
        return generate_chart(obstruents, am_plosives, ops_plosives, plosive_segments)
    if affricates.shape[0] == 0: # only plosives and fricatives
        obstruents = pd.concat([plosives, fricatives])
        ops_plofrics = get_variation_set(FRIC, obstruents)
        for am in am_plosives:
            for am_tup in am:
                am_tup[CAT] = "plosive"
        for am in am_fricatives:
            for am_tup in am:
                am_tup[CAT] = "fricative"
        ops_plofrics, ops_plofrics_short = collapse_oral_params(obstruents, am_plosives + am_fricatives, ops_plofrics, pd.concat([plosive_segments, fricative_segments]))
        return generate_chart(obstruents, am_plosives + am_fricatives, ops_plofrics, pd.concat([plosive_segments, fricative_segments]))
    # treat affricates as plosives
    obstruents = pd.concat([plosives, affricates, fricatives])
    segments = pd.concat([plosive_segments, affricate_segments, fricative_segments])
    ploaffs = pd.concat([plosives, affricates])
    am_plos_aff = []
    for am in add_unicode_airstream_params(ploaffs, get_variation_set(AIRS, ploaffs)):
        new_tup = tuple()
        for am_tup in am:
            for ca in ["plosive", "affricate"]:
                am_plus = deepcopy(am_tup)
                am_plus[CAT] = ca
                new_tup += tuple([am_plus])
        am_plos_aff.append(new_tup)
    op_aff_as_plos = get_variation_set(AFFR, obstruents)
    am_fricatives_aap = deepcopy(am_fricatives)
    for am in am_fricatives_aap: #already has unicodes added
        for am_tup in am:
            am_tup[CAT] = "fricative"
    # add category info to oral params if required (distinguishes e.g. k from kx, c from cç)
    op_aff_as_plos = add_category(obstruents, am_plos_aff + am_fricatives_aap, op_aff_as_plos)
    am_aff_as_plos, am_aff_as_plos_short = collapse_airstream_mechanisms(obstruents, am_plos_aff + am_fricatives_aap, op_aff_as_plos, segments)
    op_coll_aff_as_plos, op_coll_aff_as_plos_short = collapse_oral_params(obstruents, am_aff_as_plos, op_aff_as_plos, segments)
    aff_as_plos = generate_chart(obstruents, am_aff_as_plos, op_coll_aff_as_plos, segments)
    if not full_headers:
        aff_as_plos = replace_headers(aff_as_plos, am_aff_as_plos_short, op_coll_aff_as_plos_short)
    # if affricates are all coronal, assume they are stop-like
    # otherwise, check to see which table is denser: stop-like or series-like
    if ((affricates[CORO] != "+") & (affricates[CORO] != "+,-")).any():
        op_aff_as_series = get_variation_set(AFFR, obstruents)
        to_delete = set()
        for i in range(0, len(op_aff_as_series)):
            for j in range(i+1, len(op_aff_as_series)):
                if op_aff_as_series[i] == op_aff_as_series[j]:
                    to_delete.add(j)
        for i in sorted(to_delete, reverse=True):
            del op_aff_as_series[i]
        am_plosives_aas = deepcopy(am_plosives)
        for am in am_plosives_aas:
            for am_tup in am:
                am_tup[CAT] = "plosive"
        am_affricates_aas = []
        for am in am_affricates:
            new_tup = tuple()
            for am_tup in am:
                am_tup[CAT] = "affricate"
                new_tup += tuple([am_tup])
            am_affricates_aas.append(new_tup)
        am_fricatives_aas = deepcopy(am_fricatives)
        for am in am_fricatives_aas:
            for am_tup in am:
                am_tup[CAT] = "fricative"
        am_aff_as_series, am_aff_as_series_short = collapse_airstream_mechanisms(obstruents, am_plosives_aas + am_affricates_aas + am_fricatives_aas, op_aff_as_series, segments)
        op_coll_aff_as_series, op_coll_aff_as_series_short = collapse_oral_params(obstruents, am_aff_as_series, op_aff_as_series, segments)
        aff_as_series = generate_chart(obstruents, am_aff_as_series, op_coll_aff_as_series, segments)
        if not full_headers:
            aff_as_series = replace_headers(aff_as_series, am_aff_as_series_short, op_coll_aff_as_series_short)
        if count_empties(aff_as_plos) < count_empties(aff_as_series):
            return aff_as_plos
        else:
            print("Treating affricates as a series in language " + list(set(plosives.Glottocode))[0], file=sys.stderr)
            return aff_as_series
    else:
        return aff_as_plos


# Adds category information to oral_params if necessary to distinguish segments
#  If some affricate/plosive pair exists that is otherwise identical except for category
#  (i.e. typically delayedRelease) information, then the input oral parameter that
#  selects both of these segments is split into two different parameters: one of which is
#  category:affricate | category:fricative; and one of which is category:plosive.
#  The affricate | fricative parameter ensures affricates will match fricatives first
#  before plosives do.
# Input: data             a data frame of all the obstruents in the language
#        airstream_params a list of tuples of airstream parameters
#        oral_params      a list of tuples of oral parameters
# Output: new_oral_params either identical to oral_params or augmented with category information
def add_category(data, airstream_params, oral_params):
    new_oral_params = []
    for op in oral_params:
        cat_split = False
        for am in airstream_params:
            op_am_segments = get_segments(data, am, op)
            if not cat_split and op_am_segments.shape[0] > 1 and len(set(op_am_segments.category)) > 1 : # we need to separate affr/fric from plos
                cat_split = True
                op_affr_fric = tuple()
                op_plos = tuple()
                for op_tup in op:
                    for val in ["affricate", "fricative"]:
                        op_tup_copy = deepcopy(op_tup)
                        op_tup_copy[CAT] = val
                        op_affr_fric += tuple([op_tup_copy])
                    op_tup_copy = deepcopy(op_tup)
                    op_tup_copy[CAT] = "plosive"
                    op_plos += tuple([op_tup_copy])
                new_oral_params.append(op_affr_fric)
                new_oral_params.append(op_plos)
        if not cat_split:
            new_oral_params.append(op)
    return new_oral_params


# Manipulates category information to an oral_params transformation if necessary.
#  If category information is added to oral_params, it means that add_category() has
#  been called at some point, and there is some pair of segments distinguishable only
#  according to category information for that particular parameterization. However, 
#  for transformation rules to still apply as designed, the transformation needs to
#  be "blind" to this past. The consequence is that if category:plosive is present, and a 
#  transformation rule applies, then it should be modified to category:fricative. This 
#  maintains the distinct affricate column, while allowing for fricatives to be in either.
#  Fricatives need to be permitted so that the possibility of t -> θ still occurs, e.g. 
#  Col. 1 (affricate|fricative): ts/dz/s/z, Col. 2 (plosive|fricative): t/d/θ/ð
# Input:  op         a list of tuples of airstream parameters
# Output: either identical to op or modified with category:fricative, if op had category:plosive
def tf_op_category(op):
    plosive_in_op = False
    for op_tup in op:
        if CAT in op_tup and op_tup[CAT] == "plosive":
            plosive_in_op = True
    if not plosive_in_op:
        return op
    op_fricative = deepcopy(op)
    for op_tup in op_fricative:
        if CAT in op_tup and op_tup[CAT] == "plosive":
            op_tup[CAT] = "fricative"
    return op_fricative


# counts all "" and "-" in chart
# Input: chart a list of lists that represents the phonological chart
# Output: a count of all ""s and "-"s in the chart
def count_empties(chart):
    return sum([x.count("") for x in chart[1:-2]]) + sum([x.count("-") for x in chart[1:-2]])


# Input 1: data, a pandas DataFrame of phonological variables
#          params, a tuple of dictionaries (interpreted as logical OR ) that lists parameters
# Input 2: data, a pandas Dataframe of phonological variables
#          airstream_params, a tuple of dictionaries (interpreted as logical OR ) that lists parameters
#          oral_params, a tuple of dictionaries (interpreted as logical OR ) that lists parameters
# Output:  segments, a pandas dataframe containing all the segments that filtered for by the parameter inputs
def get_segments(*args):
    if len(args) == 2:
        data = args[0]
        params_no_unicode, params_unicodes = remove_unicode(args[1])
        if len(params_unicodes) == 0:
            return get_segments_helper(data, params_no_unicode)
        else:
            segments_no_unicode = get_segments_helper(data, params_no_unicode)
            ret_segments = pd.DataFrame(columns=data.columns)
            for uc_dict in params_unicodes:
                this_dict_segments = deepcopy(segments_no_unicode)
                for uc in uc_dict:
                    if uc == DEV_STR: # we have to look for two different characters
                        if uc_dict[uc] == "+":
                            this_dict_segments = this_dict_segments[(this_dict_segments.Phoneme.str.contains(DEVS[0])) | (this_dict_segments.Phoneme.str.contains(DEVS[1]))]
                        else: # "-"
                            this_dict_segments = this_dict_segments[~((this_dict_segments.Phoneme.str.contains(DEVS[0])) | (this_dict_segments.Phoneme.str.contains(DEVS[1])))]
                    else: # it's a normal unicode character
                        uc_chr = uc[1]
                        if uc_dict[uc] == "+":
                            this_dict_segments = this_dict_segments[this_dict_segments.Phoneme.str.contains(uc_chr)]
                        else: # "-"
                            this_dict_segments = this_dict_segments[~this_dict_segments.Phoneme.str.contains(uc_chr)]
                ret_segments = pd.concat([ret_segments, this_dict_segments]) if ret_segments.shape[0] != 0 else this_dict_segments
            return ret_segments
    elif len(args) == 3:
        data = args[0]
        airstream_param_no_unicode, airstream_unicodes = remove_unicode(args[1])
        oral_param_no_unicode, oral_unicodes = remove_unicode(args[2])
        param_dicts = []
        for ap in airstream_param_no_unicode:
            for op in oral_param_no_unicode:
                param_dicts.append({**ap, **op})
        segments = pd.DataFrame(columns=data.columns)
        for p in param_dicts:
            new_segments = get_segments_helper(data, p)
            segments = pd.concat([segments, new_segments]) if segments.shape[0] != 0 else new_segments
        if len(airstream_unicodes) == 0 and len(oral_unicodes) == 0: # no unicodes, we're done
            return segments.drop_duplicates()
        else: # both unicode tuples have something in them
            ret_segments = pd.DataFrame(columns=data.columns)
            for uc_dict in airstream_unicodes:
                for uc_dict_2 in oral_unicodes:
                    this_combo_segments = deepcopy(segments)
                    combo_dict = {**uc_dict, **uc_dict_2}
                    for uc in combo_dict:
                        if uc == DEV_STR: # we have to look for two different characters
                            if uc_dict[uc] == "+":
                                this_combo_segments = this_combo_segments[((this_combo_segments.Phoneme.str.contains(DEVS[0])) | this_combo_segments.Phoneme.str.contains(DEVS[1]))]
                            else: # "-"
                                this_combo_segments = this_combo_segments[~((this_combo_segments.Phoneme.str.contains(DEVS[0])) | this_combo_segments.Phoneme.str.contains(DEVS[1]))]
                        else:
                            uc_chr = uc[1]
                            if combo_dict[uc] == "+":
                                this_combo_segments = this_combo_segments[this_combo_segments.Phoneme.str.contains(uc_chr)]
                            else: # "-"
                                this_combo_segments = this_combo_segments[~this_combo_segments.Phoneme.str.contains(uc_chr)]
                    ret_segments = pd.concat([ret_segments, this_combo_segments]) if ret_segments.shape[0] != 0 else this_combo_segments
            return ret_segments
    else:
        raise Exception("get_segments() only accepts 2 or 3 arguments")


# recursive helper function for get_segments_helper()
def get_segments_helper(data, params):
    if type(params) is dict: #easy case
        if any(x in params for x in SP_CHRS):
            non_sp_chrs_params = {k:params[k] for k in params.keys() if k not in SP_CHRS}
            ret = data[(data[list(non_sp_chrs_params)] == pd.Series(non_sp_chrs_params, dtype="string")).all(axis=1)]
            for sp_chr in {k:params[k] for k in params.keys() if k in SP_CHRS}:
                if params[sp_chr] == "+":
                    ret = ret[ret.Phoneme.str.contains(sp_chr[1])]
                else: # "-"
                    ret = ret[~ret.Phoneme.str.contains(sp_chr[1])]
            return ret
        return data[(data[list(params)] == pd.Series(params, dtype="string")).all(axis=1)]
    else: # tuple
        segments = pd.DataFrame(columns=data.columns)
        for p in params:
            if any([x in p for x in SP_CHRS]):
                non_sp_chrs_p = {k:p[k] for k in p.keys() if k not in SP_CHRS}
                p_segments = data[(data[list(non_sp_chrs_p)] == pd.Series(non_sp_chrs_p, dtype="string")).all(axis=1)]
                for sp_chr in {k:p[k] for k in p.keys() if k in SP_CHRS}:
                    if p[sp_chr] == "+":
                        p_segments = p_segments[p_segments.Phoneme.str.contains(sp_chr[1])]
                    else: # "-"
                        p_segments = p_segments[~p_segments.Phoneme.str.contains(sp_chr[1])]
                segments = pd.concat([segments, p_segments]) if segments.shape[0] != 0 else p_segments
            else:
                if segments.shape[0] != 0:
                    segments = pd.concat([segments, data[(data[list(p)] == pd.Series(p, dtype="string")).all(axis=1)]])
                else:
                    segments = data[(data[list(p)] == pd.Series(p, dtype="string")).all(axis=1)]
        return segments.drop_duplicates()


# collapses airstream parameters according to a small number of rules
# Input: data             the data to be examined
#        airstream_params the airstream parameters
# Output: new_airstream_params the airstream_params after collapsing them
def collapse_airstream_mechanisms(data, airstream_params, oral_params, segments):
    new_airstream_params = deepcopy(airstream_params)
    short_airstream_params = deepcopy(airstream_params)
    params_to_segments = dict()
    params_to_theoretical = dict()
    for i in range(0, len(new_airstream_params)):
        params_to_segments[i] = [set(get_segments(data, new_airstream_params[i]).Phoneme)]
        params_to_theoretical[i] = set(get_segments(segments, new_airstream_params[i]).Phoneme)
    for tf in AIRSTREAM_TRANSFORMATIONS:
        to_delete = set()
        for i in range(0, len(new_airstream_params)):
            am_tup = new_airstream_params[i]
            transformed_dicts = transform_dicts(am_tup, tf)
            if len(transformed_dicts) > 0:
                am_phonemes = params_to_segments[i]
                am_tf_phonemes = set(get_segments(data, tuple(transformed_dicts)).Phoneme)
                # if the transformation is even theoretically useless, go to the next rule
                theoretical_tf_segments = get_segments(segments, tuple(transformed_dicts))
                if (len(am_tf_phonemes) == 0 and theoretical_tf_segments.shape[0] == 0):
                    continue
                matching_index = -1
                partial_match = False # a partial match is bad -> leads to double counting
                for j in [x for x in params_to_segments.keys() if x != i]:
                    if any([all(y in am_tf_phonemes for y in x) for x in params_to_segments[j]]):
                        matching_index = j
                    # if there is some set in params_to_segments which contains all of 
                    # am_tf_phonemes, but also something else, we will want to skip this
                    if len(am_tf_phonemes) > 0 and any(all([x in segset for x in am_tf_phonemes]) and not all([x in am_tf_phonemes for x in segset]) for segset in params_to_segments[j]):
                        partial_match = True
                if partial_match:
                    continue
                if matching_index == -1:
                    new_param_set = am_tup + tuple([x for x in transformed_dicts if x not in am_tup])
                else:
                    new_param_set = am_tup + tuple([x for x in new_airstream_params[matching_index] if x not in am_tup]) + tuple([x for x in transformed_dicts if x not in am_tup and x not in new_airstream_params[matching_index]])
                # if the transformation generates a collision, go to next tf rule
                am_and_tf_phonemes = get_segments(data, new_param_set)
                theoretical_segments = set(get_segments(segments, new_param_set).Phoneme)
                collision = False
                for op in oral_params:
                    op_filter = get_segments(am_and_tf_phonemes, op)
                    if op_filter.shape[0] > 1:
                        collision = True
                        break
                if collision:
                    continue
                # if this rule is non-matching, check and see if it generates a theoretical duplicant
                if matching_index == -1:
                    for j in [x for x in params_to_theoretical.keys() if x != i]:
                        if any([x in params_to_theoretical[j] for x in theoretical_segments]):
                            collision = True
                            break
                if collision:
                    continue
                # no collisions, if we get here, this transformation is good
                to_delete.add(i)
                new_airstream_params.append(new_param_set)
                params_to_theoretical[len(new_airstream_params)-1] = theoretical_segments
                if matching_index == -1:
                    params_to_segments[len(new_airstream_params)-1] = am_phonemes
                    short_airstream_params.append(short_airstream_params[i])
                else:
                    to_delete.add(matching_index)
                    params_to_segments[len(new_airstream_params)-1] = am_phonemes + params_to_segments[matching_index]
                    # modify short params
                    am_short_params = short_airstream_params[i]
                    short_transformed_dicts = transform_dicts(am_short_params, tf)
                    if matching_index == -1:
                        new_short_param_set = am_short_params + tuple([x for x in short_transformed_dicts if x not in am_short_params])
                    else:
                        new_short_param_set = am_short_params + tuple([x for x in short_airstream_params[matching_index] if x not in am_short_params]) + tuple([x for x in short_transformed_dicts if x not in am_short_params and x not in short_airstream_params[matching_index]])
                    short_airstream_params.append(new_short_param_set)
        if len(to_delete) > 0:
            for td in sorted(to_delete, reverse=True):
                del new_airstream_params[td]
                del params_to_segments[td]
                del params_to_theoretical[td]
                del short_airstream_params[td]
            pts_keys = sorted([x for x in params_to_segments.keys()])
            for i in range(0, len(new_airstream_params)):
                params_to_segments[i] = params_to_segments.pop(pts_keys[i])
                params_to_theoretical[i] = params_to_theoretical.pop(pts_keys[i])
    return new_airstream_params, short_airstream_params


# collapses oral params according to transformation rules
# transformation rules are performed in the order given
# two oral_params are collapsed with an OR statement, if they match a transformation
#    rule and if there are no collisions
# Input: data             the data to be examined
#        airstream_params the airstream parameters
#        oral_params      the oral parameters, as they exist
#        segments         a Pandas dataframe of all possible relevant segments
# Output: a tuple of two items:
#         new_oral_params the oral_params after collapsing them with transformation rules
#         new_short_oral_params the oral params without all the do-nothing rules
def collapse_oral_params(data, airstream_params, oral_params, segments):
    new_oral_params = add_unicode_oral_params(data, airstream_params, deepcopy(oral_params))
    new_oral_params.sort(key = sort_op)
    short_oral_params = deepcopy(new_oral_params)
    params_to_segments = dict()
    params_to_theoretical = dict()
    for i in range(0, len(new_oral_params)):
        params_to_segments[i] = [set(get_segments(data, new_oral_params[i]).Phoneme)]
        params_to_theoretical[i] = set(get_segments(segments, new_oral_params[i]).Phoneme)
    # set the transformation list based on the presence of clicks
    tfs = BASIC_TRANSFORMATIONS
    # perform transformations
    for tf in tfs:
        to_delete = set()
        for i in range(0, len(new_oral_params)):
            op_tup = new_oral_params[i]
            transformed_dicts = tf_op_category(transform_dicts(op_tup, tf))
            if len(transformed_dicts) > 0: # something matched the transformation rule
                if transformed_dicts in new_oral_params or all([x in op_tup for x in transformed_dicts]): # this is a copy of something that exists
                    continue
                op_phonemes = params_to_segments[i]
                op_tf_phonemes = set(get_segments(data, tuple(transformed_dicts)).Phoneme)
                # if the transformation is even theoretically useless, go to the next rule
                theoretical_tf_segments = get_segments(segments, tuple(transformed_dicts))
                if (len(op_tf_phonemes) == 0 and theoretical_tf_segments.shape[0] == 0):
                    continue
                matching_index = -1
                partial_match = False # a partial match is bad -> leads to double counting
                for j in [x for x in params_to_segments.keys() if x != i]:
                    if any([all(y in op_tf_phonemes for y in x) for x in params_to_segments[j]]):
                        matching_index = j
                    # if there is some set in params_to_segments which contains all of 
                    # op_tf_phonemes, but also something else, we will want to skip this
                    if len(op_tf_phonemes) > 0 and any(all([x in segset for x in op_tf_phonemes]) and not all([x in op_tf_phonemes for x in segset]) for segset in params_to_segments[j]):
                        partial_match = True
                if partial_match:
                    continue
                if matching_index == -1:
                    new_param_set = op_tup + tuple([x for x in transformed_dicts if x not in op_tup])
                else:
                    new_param_set = op_tup + tuple([x for x in new_oral_params[matching_index] if x not in op_tup]) + tuple([x for x in transformed_dicts if x not in op_tup and x not in new_oral_params[matching_index]])
                # if the transformation generates a collision, go to next tf rule
                op_and_tf_phonemes = get_segments(data, new_param_set)
                theoretical_segments = set(get_segments(segments, new_param_set).Phoneme)
                collision = False
                for am in airstream_params:
                    am_filter = get_segments(op_and_tf_phonemes, am)
                    if am_filter.shape[0] > 1:
                        collision = True
                        break
                if collision:
                    continue
                # if this rule is non-matching, check and see if it generates a theoretical duplicant
                if matching_index == -1:
                    for j in [x for x in params_to_theoretical.keys() if x != i]:
                        if any([x in params_to_theoretical[j] for x in theoretical_segments]):
                            collision = True
                            break
                if collision:
                    continue
                # if we get here, this transformation is good
                to_delete.add(i)
                new_oral_params.append(new_param_set)
                params_to_theoretical[len(new_oral_params)-1] = theoretical_segments
                if matching_index == -1:
                    params_to_segments[len(new_oral_params)-1] = op_phonemes
                    short_oral_params.append(short_oral_params[i])
                else:
                    to_delete.add(matching_index)
                    params_to_segments[len(new_oral_params)-1] = op_phonemes + params_to_segments[matching_index]
                    # modify short params
                    op_short_params = short_oral_params[i]
                    short_transformed_dicts = tf_op_category(transform_dicts(op_short_params, tf))
                    if matching_index == -1:
                        new_short_param_set = op_short_params + tuple([x for x in short_transformed_dicts if x not in op_short_params])
                    else:
                        new_short_param_set = op_short_params + tuple([x for x in short_oral_params[matching_index] if x not in op_short_params]) + tuple([x for x in short_transformed_dicts if x not in op_short_params and x not in short_oral_params[matching_index]])
                    short_oral_params.append(new_short_param_set)
        if len(to_delete) > 0:
            for td in sorted(to_delete, reverse=True):
                del new_oral_params[td]
                del params_to_segments[td]
                del params_to_theoretical[td]
                del short_oral_params[td]
            pts_keys = sorted([x for x in params_to_segments.keys()])
            for i in range(0, len(new_oral_params)):
                params_to_segments[i] = params_to_segments.pop(pts_keys[i])
                params_to_theoretical[i] = params_to_theoretical.pop(pts_keys[i])
            # resort all variables
            old_order = [i for i in range(len(new_oral_params))]
            new_order = sorted([i for i in range(len(new_oral_params))], key = lambda x: sort_op(short_oral_params[x]))
            new_oral_params = [new_oral_params[i] for i in new_order]
            short_oral_params = [short_oral_params[i] for i in new_order]
            params_to_segments = {i: params_to_segments[new_order[i]] for i in old_order}
            params_to_theoretical = {i: params_to_theoretical[new_order[i]] for i in old_order}
    return new_oral_params, short_oral_params


# transforms a tuple of dicts according to a transformation rule
# Input: param_tup a tuple of dicts which represents phonological parameters
#        tf        the transformation rule, expressed as a 2-tuple of dicts
# Output: transformed_dicts the input param_tup modified according to tf
def transform_dicts(param_tup, tf):
    transformed_dicts = []
    for op in param_tup:
        convert_useful = False
        for feature in tf[1].keys():
            if feature in op.keys() and op[feature] != tf[1][feature]:
                convert_useful = True
                break
        if not convert_useful:
            return transformed_dicts
        this_converts = True
        for feature in tf[0].keys():
            if feature in op.keys() and op[feature] != tf[0][feature]:
                this_converts = False
                break            
        if this_converts: # generate transformed oral parameter, add to transformed_dicts
            op_tf = deepcopy(op)
            something_new = False
            for feature in tf[1].keys():
                op_tf[feature] = tf[1][feature]
                if feature in op.keys() and op_tf[feature] != op[feature]:
                    something_new = True
            if something_new:
                transformed_dicts.append(op_tf)
    return transformed_dicts


# adds unicode characters to the oral parameterization
# the unicode characters being considered are those in ORALCHR
# Input: data             the data to be examined
#        airstream_params the airstream parameters
#        oral_params      the oral parameters, as they exist
# Output: new_oral_params the oral_params augmented with unicode characters, 
#                         if there are any
def add_unicode_oral_params(data, airstream_params, oral_params):
    sp_chrs_present = set()
    for s in set(data.Phoneme):
        if len(s) > 1:
            for sc in ORALCHR:
                if sc[1] in s:
                    sp_chrs_present.add(sc)
    if len(sp_chrs_present) == 0:
        return oral_params
    sp_chr_matters = []
    sp_chr_combos = []
    for op in oral_params:
        op_filter = get_segments(data, op)
        this_scs = set()
        for am in airstream_params:
            am_filter = get_segments(op_filter, am)
            if am_filter.shape[0] > 1:
                for sc in sp_chrs_present:
                    if any([sc[1] in x for x in am_filter.Phoneme]):
                        this_scs.add(sc)
                        if op not in sp_chr_matters:
                            sp_chr_matters.append(op)
        if len(this_scs) != 0:
            sp_chr_combos.append(this_scs)
    if len(sp_chr_matters) == 0:
        return oral_params
    new_oral_params = [x for x in oral_params if x not in sp_chr_matters]
    for nop in new_oral_params:
        for nop_tup in nop:
            for sp in sp_chrs_present:
                if all([sp[1] in x for x in get_segments(data, nop_tup).Phoneme]):
                    nop_tup[sp] = "+"
                else:
                    nop_tup[sp] = "-"
    for i in range(0, len(sp_chr_matters)):
        op = sp_chr_matters[i]
        sp_combos = list(sp_chr_combos[i])
        for combo in list(itertools.product(['+','-'], repeat=len(sp_combos))):
            new_oral_params.append(deepcopy(op))
            for j in range(0, len(sp_combos)):
                for nop in new_oral_params[-1]:
                    nop[sp_combos[j]] = combo[j]
    return new_oral_params


# adds unicode characters to the airstream parameterization
# the unicode characters being considered are those in DEVS / DEVSS and LCOR / LCORS
# Input: data             the data to be examined
#        airstream_params the airstream_params as they exist
# Output: new_airstream_params the airstream_params augmented with unicode characters,
#                              if there are any
def add_unicode_airstream_params(data, airstream_params):
    dev_present = False
    lcor_present = False
    for s in set(data.Phoneme):
        if len(s) > 1:
            for d in DEVS:
                if d in s:
                    dev_present = True
            if LCOR in s:
                lcor_present = True
    if not dev_present and not lcor_present:
        return airstream_params
    ret_airstream_params = deepcopy(airstream_params)
    if dev_present:
        dev_matters = []
        for am in ret_airstream_params:
            am_filter = get_segments(data, am)
            devo_here = False
            if am_filter.shape[0] > 1 and (any([DEVS[0] in x for x in am_filter.Phoneme]) or any([DEVS[1] in x for x in am_filter.Phoneme])):
                if am not in dev_matters:
                    dev_matters.append(am)
                    devo_here = True
        if len(dev_matters) > 0:
            new_airstream_params = []
            for changed_param in dev_matters:
                ret_airstream_params.remove(changed_param)
            for i in range(0, len(dev_matters)):
                nodev = deepcopy(dev_matters[i])
                for optup in nodev:
                    optup[DEV_STR] = "-"
                if get_segments(data, nodev).shape[0] > 0:
                    new_airstream_params.append(nodev)
                devsym = deepcopy(dev_matters[i])
                for optup in devsym:
                    optup[DEV_STR] = "+"
                new_airstream_params.append(devsym)
            ret_airstream_params += new_airstream_params
    if lcor_present:
        lcor_matters = []
        for am in ret_airstream_params:
            am_filter = get_segments(data, am)
            lcor_here = False
            if am_filter.shape[0] > 1 and any([LCOR in x for x in am_filter.Phoneme]):
                if am not in lcor_matters:
                    lcor_matters.append(am)
        if len(lcor_matters) > 0:
            new_airstream_params = []
            for changed_param in lcor_matters:
                ret_airstream_params.remove(changed_param)
            for i in range(0, len(lcor_matters)):
                nolcor = deepcopy(lcor_matters[i])
                for optup in nolcor:
                    optup[LCORS] = "-"
                if get_segments(data, nolcor).shape[0] > 0:
                    new_airstream_params.append(nolcor)
                lcorp = deepcopy(lcor_matters[i])
                for optup in lcorp:
                    optup[LCORS] = "+"
                new_airstream_params.append(lcorp)
            ret_airstream_params += new_airstream_params
    return ret_airstream_params


# generates a phonological chart, given its paramaters
# Input: data             the data to be arranged into a chart
#        airstream_params a list of tuples of dicts of the airstream parameters
#        oral_params      a list of tuples of dicts of the oral_params
#        segments         a data frame of all relevant segments in all languages
# Output: the data rearranged into a chart (list of lists)
def generate_chart(data, airstream_params, oral_params, segments):
    if data.shape[0] == 0:
        return [""]
    header = [""]
    for oral_param in oral_params:
        header.append(" | ".join([";".join([":".join([k, v]) for k, v in x.items()]) for x in oral_param]))
    header += ["count", "markedness_count"]
    all_series = []
    all_markedness_series = []
    for air_mech in airstream_params:
        series = [" | ".join([";".join([":".join([k, v]) for k, v in tup.items()]) for tup in air_mech])]
        markedness_series = []
        am_filter = get_segments(data, air_mech)
        for op in oral_params:
            op_phonemes = list(get_segments(am_filter, op).Phoneme.values)
            op_filt = tuple([x for x in op])
            poss_segments = get_segments(segments, air_mech, op_filt)
            markedness_series.append(sum(poss_segments.Count))
            if len(op_phonemes) != 0:
                series.append("//".join(op_phonemes))
                if len(op_phonemes) > 1:
                    print(am_filter.Glottocode.values[0] + ": indistinguishable segments " + "/".join(op_phonemes), file=sys.stderr)
            elif poss_segments.shape[0] != 0: # at least one other segment like this exists
                series.append("")
            else:
                series.append("-")
        all_series.append(series)
        all_markedness_series.append(markedness_series)
    # append row counts
    for i in range(len(all_series)):
        segs = all_series[i][1:]
        poss = [i for i in range(len(segs)) if segs[i] != "-"]
        pres = [i for i in range(len(segs)) if segs[i] != "-" and segs[i] != ""]
        markedness_total = sum(all_markedness_series[i])
        markedness_pres = sum(all_markedness_series[i][j] for j in pres)
        all_series[i].append(str(len(pres)) + "/" + str(len(poss)))
        all_series[i].append(str(markedness_pres/markedness_total))
    # append column counts
    col_counts = [""]
    col_markedness_counts = [""]
    for j in range(1,len(all_series[0])-2): # the first column is airstream mechanism, last 2 we just added
        col = [row[j] for row in all_series]
        poss = [i for i in range(len(col)) if col[i] != "-"]
        pres = [i for i in range(len(col)) if col[i] != "-" and col[i] != ""]
        markedness_total = sum(row[j-1] for row in all_markedness_series)
        markedness_pres = sum(all_markedness_series[i][j-1] for i in pres)
        col_counts.append(str(len(pres)) + "/" + str(len(poss)))
        col_markedness_counts.append(str(markedness_pres/markedness_total))
    return [header] + all_series + [col_counts] + [col_markedness_counts]


# replaces the airstream_param and oral_param headers in a chart - used for abbreviations
# note that oral_params must be the exact same lengths as the original chart
# Input: chart            a chart (i.e. outputted by generate_chart)
#        airstream_params the replacement airstream_params (must be in the same order as original)
#        oral_params      the replacement oral_params (must be in the same order as original)
# Output: the chart with the oral_params and airstream_params replaced
def replace_headers(chart, airstream_params, oral_params):
    if chart == [""]:
        return chart
    if len(chart) != len(airstream_params) + 3 or len(chart[0]) != len(oral_params) + 3:
        print(len(chart))
        print(len(airstream_params))
        print(len(chart[0]))
        print(len(oral_params))
        raise Exception("airstream_params and oral_params must match dimensions of chart")
    header = [""]
    for oral_param in oral_params:
        header.append(" | ".join([";".join([":".join([k, v]) for k, v in x.items()]) for x in oral_param]))
    header += ["count", "markedness_count"]
    new_chart = []
    for i in range(0,len(airstream_params)):
        new_chart.append([shorter_airstream_text(airstream_params[i])] + chart[i+1][1:])
    return [header] + new_chart + chart[-2:]


# collapses category: information if it's the only variation, otherwise returns dictionary as a string
# Input: am  a tuple of dicts of airstream_parameters
# Output: a string representation of am, with category: information collapsed
def shorter_airstream_text(am):
    if len(am) == 1:
        return ";".join([":".join([k, v]) for k, v in am[0].items()])
    only_cat_differs = True
    for i in range(0,len(am)):
        for j in range(i+1, len(am)):
            for k in [x for x in am[i].keys() if x != CAT]:
                if am[i][k] != am[j][k]:
                    only_cat_differs = False
                    break
            if not only_cat_differs:
                break
        if not only_cat_differs:
            break
    if only_cat_differs:
        return ";".join([":".join([k, v]) for k, v in am[0].items() if k != CAT] + ["category:" + "|".join(sorted(set([am[i][CAT] for i in range(0, len(am))])))])
    else:
        return " | ".join([";".join([":".join([k, v]) for k, v in tup.items()]) for tup in am])


# removes the unicode values (defined in SP_CHRS) from the param dict
# Input: param a tuple of dicts describing parameters
# Output: tup, sprchs
#         tup    the param input minus all unicode characters
#         sprchs a tuple of dicts of all the SP_CHRS removed from param, and their values,
#                e.g. {(chr1: "+", chr2: "-"), (chr1: "-", chr2: "_+"")}
def remove_unicode(param):
    if type(param) == dict:
        return tuple([{k: param[k] for k in param.keys() - SP_CHRS}]), tuple([{x:param[x] for x in SP_CHRS if x in param}])
    else: # tuple
        ret_tup = tuple()
        ret_spchrs = tuple()
        for p_tup in param:
            add_tup = {k: p_tup[k] for k in p_tup.keys() - SP_CHRS}
            add_spchrs = {k: p_tup[k] for k in SP_CHRS if k in p_tup.keys()}
            if not any(add_tup == x for x in ret_tup):
                ret_tup += tuple([add_tup])
            if not any(add_spchrs == x for x in ret_spchrs):
                ret_spchrs += tuple([add_spchrs])
        return ret_tup, ret_spchrs


# writes the chart to file
# Input: chart      a list of lists giving the phonological chart
#        path_loc   the location where to write the file
#        file_name  where to write out to file
#        headerless boolean whether to write a headerless version of the file
# Output: None. Writes the chart to path_loc/file_name.tsv
def write_chart(chart, path_loc, file_name, headerless = False):
    if not os.path.exists(path_loc):
        os.makedirs(path_loc, exist_ok=True)
    with open(path_loc + "/" + file_name + ".tsv", "w+") as f:
        for line in chart:
            _ = f.writelines("\t".join(line) + "\n")
    if headerless:
        with open(path_loc + "/" + file_name + "_noheaders.tsv", "w+") as f:
            for line in chart[1:]:
                _ = f.writelines("\t".join(line[1:]) + "\n")


# generates the set of variations which are necessary to fully distinguish input data
# Input: feature_set  the set of features the variation of which will be evaluated 
#        data         the data against which variation will be evaluated
# Output: A list of tuples of dicts which give all the variation of features in 
#         the feature_set for this data.  Tuples because they are tuples everywhere else.
#         E.g., [({labial:+}), ({labial:-, coronal:+}), ({labial:-, coronal:-})]
def get_variation_set(feature_set, data):
    features = dict()
    features_with_zero = dict()
    predictable_features = dict()
    for f in feature_set:
        features_with_zero[f] = set(data[f])
        features[f] = set(data[f]) - {"0"}  # remove not applicables ('0')
        predictable_features[f] = list()
    # if there's any points of variation, we will do some stuff
    if any([len(x) > 1 for x in features.values()]):
        # remove invariant features
        for f in feature_set:
            if len(features[f]) == 1:
                del features[f]
                del features_with_zero[f]
        # find predictable features
        for var1 in features_with_zero.keys():
            if len(features_with_zero[var1]) > 1:  # find co-variation
                for var2 in features_with_zero.keys() - {var1}:
                    if len(features_with_zero[var2]) > 1:
                        redundant = True
                        for val in features_with_zero[var1]:
                            if len(set(data[var2][data[var1] == val])) > 1:
                                redundant = False
                        if redundant:
                            predictable_features[var1].append(var2)
        # find mutually predictable features (=redundant)
        # TODO: mutually_predictable can be reintroduced for shortening the output
        # but it leads to stochastic differences in outputs, depending
        # on which feature is deleted, interacting with tranformation rules
        # and must be omitted
        # mutually_predictable = set()
        # for a in feature_set:
        #     if len(predictable_features[a]) > 0:
        #         for b in predictable_features[a]:
        #             if a in predictable_features[b]:
        #                 predictable_features[a].remove(b)
        #                 predictable_features[b].remove(a)
        #                 mutually_predictable.add(frozenset([a, b]))
        # remove redundant features : but never CLIC
        # for mp in mutually_predictable:
        #     redundant = set(mp)
        #     # thou shalt not remove stridents, clicks, implosives:
        #     # stridents: we need to know /s/ may not match with a plosive
        #     # clicks: we need this for click transformations
        #     # implosives: we need this for implosive rules
        #     if STRI in redundant:
        #         redundant.remove(STRI)
        #     if CLIC in redundant:
        #         redundant.remove(CLIC)
        #     while len(redundant) > 1:
        #         f = redundant.pop()
        #         if f in features.keys():
        #             del features[f]
        #         if f in features_with_zero.keys():
        #             del features_with_zero[f]
        #         if f in predictable_features.keys():
        #             del predictable_features[f]
        #         for v in predictable_features.keys():
        #             if f in predictable_features[v]:
        #                 predictable_features[v].remove(f)
        # clean up predictability matrix
        for f in feature_set:
            if f in predictable_features.keys():
                if len(predictable_features[f]) == 0:
                    del predictable_features[f]
    # generate all combinations
    full_list = get_variation_set_helper(list(features_with_zero), dict(), data)
    # return this as a list of tuples
    return [tuple([x]) for x in full_list] # full_list


# recursive helper function for get_variation_set()
def get_variation_set_helper(features, up_to_now_dict, filtered_data):
    if len(features) == 0:
        return [up_to_now_dict]
    ret_set = []
    for val in set(filtered_data[features[0]].values):   
        next_dict = deepcopy(up_to_now_dict)
        next_dict[features[0]] = val
        downstream_dicts = get_variation_set_helper(features[1:], next_dict, filtered_data[filtered_data[features[0]] == val])
        for d in downstream_dicts:
            ret_set.append({**up_to_now_dict, **d})
    return ret_set


# converts a string of type "x/y", where x and y are integers to a float
# Input: frac_str a string of type "x/y"
# Output: frac_str as a float
def convert_to_float(frac_str):
    if "/" not in frac_str:
        return float(frac_str)
    try:
        num, denom = frac_str.split("/")
        return float(num) / float(denom)
    except ZeroDivisionError:
        print("Error: 0/0", file=sys.stderr)
        return 0


# Calculates the number of segments in a chart (for error checking)
# Input: chart, a chart
# Output: the number of segments
def count_segments(chart):
    if len(chart) == 0:
        return 0
    return sum([int(x[-2].split('/')[0]) for x in chart[1:-2]]) # first row: params, last two: place counts


# Used to sort a dictionary of oral_params
# Input: tup an oral parameter tuple of dicts
# Output: a string which will be sorted
def sort_op(op):
    return (collapse_plus_minus(op, DORS, rev=True), 
            collapse_plus_minus(op, HIGH), 
            collapse_plus_minus(op, BACK, rev=True), 
            collapse_plus_minus(op, LABI), 
            collapse_plus_minus(op, LADE, rev=True),
            collapse_plus_minus(op, ANTE), 
            collapse_plus_minus(op, DIST), 
            collapse_plus_minus(op, STRI, rev=True),
            collapse_plus_minus(op, FRON), 
            collapse_plus_minus(op, LOW, rev=True), 
            collapse_plus_minus(op, RETR, rev=True), 
            collapse_plus_minus(op, EPIL, rev=True),
            collapse_plus_minus(op, ROUN, rev=True)
            )


# Helper function for sorting oral_params
# Input: op    an oral parameter tuple of dicts
#        feat  the feature to collapse the minuses and pluses
# Return: a string which represents all the values the parameter can take, sorted with + high and - low (and 0 even lower)
def collapse_plus_minus(op, feat, rev=False):
    vals = set()
    for op_tup in op:
        if feat in op_tup and op_tup[feat] != "0":
            if (feat == LABI or feat == FRON) and op_tup[feat].endswith("-,+"): # I want + to the front, then -, then -,+
                vals.add("z," + op_tup[feat])
            elif rev:
                vals.add(op_tup[feat].replace("+","v").replace("-","+").replace("v","-"))
            else:
                vals.add(op_tup[feat])
    if len(vals) == 0:
        return "0"
    if len(vals) == 1:
        return list(vals)[0]
    else:
        # privilege "clean" parameterization
        return(";".join([x for x in vals if len(x) == 1] + sorted([x for x in vals if len(x) > 1])))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p", "--phoible", help="Path to the PHOIBLE directory. Default is ../phoible/", default=str(Path(__file__).parent / ".." / "phoible"), required=False, type=str
    )
    parser.add_argument(
        "-t", "--transformations", help="Path to a tsv which describes the transformation rules.", default=str(Path(__file__).parent / ".." / "transformations" / "transformation_rules.tsv"), required=False, type=str
    )
    parser.add_argument(
        "-i", "--inventories", help="Path to a csv of the inventories being used, normally generated through process-inventories. Required.", required=True, type=str,
    )
    parser.add_argument(
        "-s", "--segments", help="Path to a csv of the markedness-weighted segment counts, normally generated through process-inventories. Required.", required=True, type=str,
    )
    parser.add_argument(
        "-l", "--languages", help="Path to a csv with a column named 'glottocode' listing the languages to process. If not provided, all of PHOIBLE is processed.", default="None", required=False, type=str,
    )
    parser.add_argument(
        "-o", "--output", help="Folder to write charts and counts to at the end of process.", default="None", required=True, type=str,
    )
    parser.add_argument(
        "-f", "--fullheaders", help="Write the full headers (not the human readable condensed ones) to the charts.", action="store_true"
    )
    parser.add_argument(
        "--suppress-charts", help="Suppresses the output of charts.", action="store_true"
    )
    args = parser.parse_args()
    if args.output[-1] != "/":
        args.output = args.output + "/"
    main(args.phoible, args.transformations, args.inventories, args.segments, args.languages, args.output, args.fullheaders, args.suppress_charts)
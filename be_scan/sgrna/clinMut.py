from Bio import Entrez, SeqIO
import xml.etree.ElementTree as ET
import time
import re
import pandas as pd
from tqdm import tqdm


#Translating three letter amino acid code to one letter aa code
aa_three_to_one = {"Ala": "A", "Arg": "R", "Asn": "N",
                  "Asp": "D", "Cys": "C", "Glu": "E",
                  "Gln": "Q", "Gly": "G", "His": "H",
                  "Ile": "I", "Leu": "L", "Lys": "K",
                  "Met": "M", "Phe": "F", "Pro": "P",
                  "Ser": "S", "Thr": "T", "Trp": "W",
                  "Tyr": "Y", "Val": "V"}

#Find Amino Acid Changes
def find_muts(df, cols):
    """[Summary]
    Given a dataframe and a list of columns containing mutations, 
    find all unique mutations in the dataframe.

    Parameters
    ------------
    df: pd.DataFrame
        Dataframe containing mutations (generated from design_library)
    cols: list of columns to search for
    
    Returns
    ------------
    mut_list: list of mutations
        List of unique mutations in all columns searched

    """
    mutations = set()
    for col in cols:
        lsts = list(map(str, (df[col])))
        mutations.update((set([mut for lst in lsts for mut in lst.split(";")])))
    mut_list = list(mutations)
    #Clean up
    mut_list.remove("nan")
    mut_list.remove("")
    
    return mut_list

def split_into_batches(lst, batch_size):
    """[Summary]
        Given a list, split the list into batches of a given size
    """
    for i in range(0, len(lst), batch_size):
        yield lst[i:i + batch_size]

def process_mutation(mut):
    """[Summary]
        Given a mutation, process the mutation to be used in a clinvar query using
        Boolean searches
    """
    if "AND" in mut:
        return "(" + mut + ")"
    else:
        return mut

def build_query(mut_list, gene):
    """[Summary]
        Builds a clinvar query, splitting query into batches of 1000.
    """
    #if contains "/", replace with AND, changing . to *, adding p. for query purposes
    mut_list_clean = list(map(lambda x: "p." + x.replace(".", "*").replace("/", " AND p."), mut_list))
    query_list = [
        f"{gene} AND ({' OR '.join([process_mutation(m) for m in batch])})"
        for batch in split_into_batches(mut_list_clean, 1000)
    ]
    return query_list
    

#convert 3 letter to 1 letter amino acid code:
def convert_aa(mut_three):
    """[Summary]
    Given a three letter amino acid mutation, convert to one letter amino acid code.

    Parameters
    ------------
    mut_three: str
        Mutation written in 3-letter amino acid code
    Returns
    ------------
    mut_one: str
        Mutation written in 1-letter amino acid code
    """
    mut_from = mut_three[:3]
    mut_loc = mut_three[3:-3]
    mut_to = mut_three[-3:]

    aa_from = aa_three_to_one[mut_from]
    aa_to = aa_three_to_one[mut_to]
    
    return aa_from + mut_loc + aa_to


#Find all the clinical mutations from ClinVar (IDS)
def fetch_clinvar_query(gene, mut_list, email = "christinelee@fas.harvard.edu"):
    """[Summary]
    Given a gene name the amino acid changes, 

    Parameters
    ------------
    gene: str
        Gene name 
    protein_change: str
        Amino acid conversion

    Returns
    ------------
    
    """
    #handle = Entrez.esearch(db = "clinvar", term = gene + ": (" + protein_change + ")")
    #handle = Entrez.esearch(db = "clinvar", term = gene, retmax = 5000)
    
    #Sample:
    #handle = Entrez.esearch(db = "clinvar", term = "CFTR: (p.Arg117His)")
    #search_term = "EZH2 AND (p.R746G OR p.E745K )"

    #Set email
    Entrez.email = email

    ids = []
    #Build boolean query based on list of genes
    search_term = build_query(mut_list, gene)
    for term in tqdm(search_term, desc="Fetching ClinVar UIDs"):
        handle = Entrez.esearch(db = "clinvar", term = term)
        record = Entrez.read(handle)
        ids.append(record["IdList"])
        handle.close()

        time.sleep(0.4) #be nice to the server!
    return [id for id_lst in ids for id in id_lst]


#Given a list of IDS, query the summary, given a batch size of 500.
def fetch_clinvar_summaries(uids, batch_size=500):
    """[Summary]
    Given a list of IDs, query summary of each ID, slicing the dataset into batches of
    500. 

    Parameters
    ------------
    uids: list of IDs
        List of IDs for the Clinvar search
    batch_size: int
        Number of IDs to query, maximum 500 to reduce load on Clinvar

    Returns
    ------------
    results: list of XML files
        List of XML files containing the summary of each ID, batched by 500
    """
    if(batch_size > 500):
        raise ValueError("Batch size cannot exceed 500")
    
    results = []
    # Split UIDs into batches
    for i in tqdm(range(0, len(uids), batch_size), desc="Fetching ClinVar summaries"):
        batch_uids = uids[i:i+batch_size]
        try:
            # Fetch data in bulk using esummary
            handle = Entrez.esummary(db="clinvar", id=",".join(batch_uids), retmode="xml")
            # Read and decode the bytes response to a string
            results.append(handle.read().decode("utf-8"))
            handle.close()
        except Exception as e:
            print(f"Error fetching batch starting at UID {i}: {e}")
        time.sleep(1) #be nice to the server!

    return results

#Parse Clinvar XML
def parse_clinvar_XML(xmls, gene):
    """[Summary]
    Given a list of XML files, parse each XML file, returning a dictionary with
    key: accession number, value: amino acid change, converted to one letter code
    """

    parsed_dict = {}
    xml_string = ''.join(xmls)

    root = ET.fromstring(xml_string)
    for doc in root.findall('.//DocumentSummary'):
        accession = doc.findtext('accession')
        if not accession:
            continue #skip entry if there is no accesssion

        #in case there are multiple gene, see if this is correct, CHECK WITH CALVIN
        gene_symbols = [gene.findtext('symbol') for gene in doc.findall('.//genes/gene')]

        if gene in gene_symbols:
            title = doc.findtext('title')
            if title:
                #Regex for finding amino acid change
                regex = "\(p\.[A-za-z]{3}\d+[A-za-z]{3}\)"
                matches = re.findall(regex, title)
                
                for match in matches:
                    #string parsing based on clinvar specifications
                    match = match.replace("(p.", "")
                    match = match.replace(")", "")
                    conv_mut = convert_aa(match)
                    parsed_dict[conv_mut] = accession         
                
    #print(parsed_dict)
    print("Clinvar XML parsed")
    return(parsed_dict)


def get_clin_muts(df, cols, gene, email):
    #Run method
    mut_list = find_muts(df, cols)
    ids = fetch_clinvar_query(gene, mut_list, email)
    files = fetch_clinvar_summaries(ids)
    return(parse_clinvar_XML(files, gene))


def insert_muts(df, col, mut_dict):
    #given a dictionary of mutations, generated from get_clin_muts, 
    #add a new column to the dataframe containing the dictionary value 
    #given that the mutation is in the column of mutations.
    for key, value in mut_dict.items():
        #print(key)
        df.loc[df[col].str.contains(key, na = False), [f"Clinvar_{col}"]] = key + ": " + value
    #df.to_csv("../../outputs/annotated_guides_clinvar_check.csv")
    return df


#Testing module
df = pd.read_csv("../../outputs/annotated_guides.csv")
cols = ["AtoG_mutations", "CtoT_mutations", "ACtoGT_mutations"]
parsed = get_clin_muts(df = df, cols = cols, gene = "AR", email = "christinelee1@college.harvard.edu")
insert_muts(df, "CtoT_mutations", parsed)
print(parsed)
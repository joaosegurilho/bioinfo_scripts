#%%
"""Explore the existence of an Interpro hit in a filtered ncbi taxonomy.  

Usage:
    python prot_w_dom.py tax_id_to_filter interpro_id output_dir

Output:
    interpro_processed_query.pkl, rendered_tree.pdf

"""
from sys import argv, version_info
if not version_info.major >= 3 and not version_info.minor >= 8:
    raise "Need python version >= 3.8"

from dataclasses import dataclass
import os
import pickle
import requests
from time import sleep

try:
    from ete3 import NCBITaxa, TextFace, SeqMotifFace, TreeStyle, faces, AttrFace, TreeStyle, NodeStyle
except Exception as e:
    raise f"Exception:{e}. Need to install ete3"
else:
    ncbi = NCBITaxa()

try:
    from tqdm import tqdm
except Exception as e:
    print(f"Exception:{e}. Tqdm instalation optional/recommended. Continuing without.")
    TQDM_ON = False
else:
    TQDM_ON = True


# https://www.ebi.ac.uk/interpro/api/taxonomy/uniprot/entry/interpro/

@dataclass
class Protein:
    prot_acc: str
    prot_len: int
    prot_tax: int
    dom_loc: list


def get_query_from_interpro(hit_id):

    base = "https://www.ebi.ac.uk/interpro/api/protein/uniprot/entry/interpro/"
    base_url = "/".join([base,hit_id])
    session = requests.Session()

    def get_pg(url):
        sleep(1.0)
        try:     
            r = session.get(url,timeout=10.0)
        except:
            print(f"The search: {r.json()} didnt work.")
            return (None, None)
        else:
            if r.ok:
                _temp_obj = r.json()
                next_page = _temp_obj['next']
                return _temp_obj, next_page

    obj, next_page = get_pg(base_url)
    total =  round(int(obj['count']) / len(obj['results'])) 
    with tqdm(desc="Query interpro: ", total=total, initial=1, leave=True, position=0, disable=TQDM_ON) as pbar:
        while next_page != None:
            try:
                new_obj, next_page = get_pg(next_page)
            except:
                continue
            else:
                for entry in new_obj['results']:
                    obj['results'].append(entry)
                pbar.update(1)

    return obj['results']

def process_interpro_obj(ip_obj, tax_to_keep):
    processed_obj = []

    for protein in ip_obj:
        prot_acc = protein['metadata']['accession']
        locations = [(x['start'],x['end']) for x in protein['entries'][0]['entry_protein_locations'][0]['fragments']]
        prot_len = protein['metadata']['length']
        tax_acc = int(protein['metadata']['source_organism']['taxId'])

        if tax_acc in tax_to_keep:
            processed_obj.append(
                Protein(prot_acc=prot_acc, prot_len=prot_len, prot_tax=tax_acc, dom_loc=locations)
                )
    
    return processed_obj

def make_tree(ipr_obj, main_tax_id, queried_dom, save_to_dir):
    colors = ['LightCyan','Khaki','LightSalmon','Lavender']
                    #Blastopirellula marina; Gemmata obscuriglobus; Gimesia maris; Planctopirus limnophila
    gene_manipul = ["314230", "214688", "344747", "521674"]
    
    def gen_tree(tree_file):
        t = tree_file
        # Customize the tree
        def layout(node):
            if node.is_root():
                for i, name in enumerate(set(node.get_children())):
                    name.img_style["bgcolor"] = colors[i-1]

            # If node is a leaf, add the nodes scientific name
            if node.is_leaf():
                faces.add_face_to_node(AttrFace("sci_name"), node, column=0)

                if node.name in gene_manipul:
                    style_leaf = NodeStyle()
                    style_leaf["bgcolor"] = "LightGreen"
                    style_leaf["fgcolor"] = "LightGreen"
                    node.set_style(style_leaf)
            

        # print(t.features)
        def qk_seq(length, doms):
            #How to normalize all prot len and motif coord???
            dom_motifs = [
                [dom[0],dom[1],'()', None, 10, "DodgerBlue", 'DodgerBlue', f"arial|6|black|{queried_dom}"] 
                for dom in doms
                ] 
            return {'seq':("-"*length), 'motifs':dom_motifs}
        
        dict_tax_w_prot = {}
        for prot in ipr_obj:
            seqFace = SeqMotifFace(**qk_seq(prot.prot_len,prot.dom_loc))
            (t & str(prot.prot_tax)).add_face(seqFace, 1, "aligned")

            textFace = TextFace(prot.prot_acc,tight_text=True)
            (t & str(prot.prot_tax)).add_face(textFace, 0, "aligned")
            
            if prot.prot_tax in dict_tax_w_prot.keys():
                dict_tax_w_prot[prot.prot_tax].append(prot.prot_acc)
            else:
                dict_tax_w_prot[prot.prot_tax] = [prot.prot_acc]
        
        for prot_tax in dict_tax_w_prot.keys():
            (t & str(prot_tax)).add_feature("num_prots",len(dict_tax_w_prot[prot_tax]))

            if len(dict_tax_w_prot[prot_tax]) == 1:
                (t & str(prot_tax)).add_feature("prot_acc",dict_tax_w_prot[prot_tax][0])
            else:
                (t & str(prot_tax)).add_feature("prot_acc",",".join(dict_tax_w_prot[prot_tax]))
        
        # #prune the tree keeping the node names specified [Cannot make it work.......]
        # leaves_with_prot = [p.prot_tax for p in ipr_obj if p.prot_tax in t.get_leaf_names()]
        # for node in tqdm(t.iter_leaves()):
        #     if node.name not in leaves_with_prot:
        #         # node.delete(prevent_nondicotomic=False,preserve_branch_length=False)
        #         node.up.remove_child(node)

        ts = TreeStyle()
        ts.layout_fn = layout
        ts.scale = 100
        ts.branch_vertical_margin = 10
        return t, ts

    final_tree, ts = gen_tree(
        ncbi.get_descendant_taxa(
            main_tax_id,
            intermediate_nodes=False,
            return_tree=True
            )
        )

    final_tree.render(save_to_dir + f"{queried_dom}_in_{main_tax_id}.pdf", w=1000, dpi=300, tree_style=ts)
    final_tree.show(tree_style=ts, name=queried_dom)

def main(tax_id:int, query_dom_id:str, save_to_dir:str):
    descendents = ncbi.get_descendant_taxa(tax_id)
    
    path_file = f"ip_obj_{tax_id}_{query_dom_id}.pkl"
    if os.path.exists(path_file):
        with open(path_file, 'rb') as f:
            process_itrpro_obj = pickle.load(f)
    else:
        itrpro_obj = get_query_from_interpro(query_dom_id)
        process_itrpro_obj = process_interpro_obj(itrpro_obj, descendents)
        with open(path_file, 'wb') as f:
            pickle.dump(process_itrpro_obj, f, pickle.HIGHEST_PROTOCOL)
    
    make_tree(process_itrpro_obj, tax_id, query_dom_id, save_to_dir)
    print("Finished script.")


if __name__ == "__main__":
    try:
        tax_id = int(argv[1]) 
        query_dom_id = str(argv[2])
        if len(argv) == 3:
            save_to_dir = str(argv[3])
        else:
            save_to_dir = "./"
    except:
        raise "Need a tax_id and a query domain, none provided."
    else:
        main(tax_id,query_dom_id,save_to_dir)

# %%

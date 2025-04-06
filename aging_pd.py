import os
import pandas as pd
import scipy.stats as stats
import itertools
import networkx as nx
import matplotlib.pyplot as plt


DESKTOP_PATH = "/Users/mariannamilano/Desktop/"
file_path = os.path.join(DESKTOP_PATH, "f1.xlsx")
na_folder = os.path.join(DESKTOP_PATH, "NA")


os.makedirs(na_folder, exist_ok=True)


df = pd.read_excel(file_path)


age_groups = ['20-29', '30-39', '40-49', '50-59', '60-69', '70-79']


non_gene_columns = ['Descriptio', 'SEX', 'AGE', 'CTRL']
gene_columns = [col for col in df.columns if col not in non_gene_columns]


selected_genes = []
for gene in gene_columns:
    group_means = [df[df['AGE'] == age][gene].mean() for age in age_groups]
    increasing = all(group_means[i] < group_means[i+1] for i in range(len(group_means)-1))
    decreasing = all(group_means[i] > group_means[i+1] for i in range(len(group_means)-1))
    if increasing != decreasing:
        selected_genes.append(gene)


results_with_pvalues = []
for gene in selected_genes:
    age_numeric = []
    expression_values = []
    for idx, age in enumerate(age_groups):
        values = df[df['AGE'] == age][gene].dropna().tolist()
        expression_values.extend(values)
        age_numeric.extend([idx] * len(values))
    pearson_r, pearson_p = stats.pearsonr(age_numeric, expression_values)
    group_data = [df[df['AGE'] == age][gene].dropna().values for age in age_groups]
    kruskal_stat, kruskal_p = stats.kruskal(*group_data)
    results_with_pvalues.append((gene, pearson_p, kruskal_p))


result_txt = os.path.join(DESKTOP_PATH, "risultato_geni.txt")
with open(result_txt, 'w') as f:
    f.write("Gene\tPearson_p\tKruskal_p\n")
    for gene, p1, p2 in results_with_pvalues:
        f.write(f"{gene}\t{p1:.4e}\t{p2:.4e}\n")



def generate_edgelist(data, gene_list, filename):
    edges = []
    for gene1, gene2 in itertools.combinations(gene_list, 2):
        try:
            stat, pval = stats.wilcoxon(data[gene1], data[gene2])
            if pval > 0.05:
                edges.append((gene1, gene2))
        except Exception:
            continue  
    df_edge = pd.DataFrame(edges, columns=["Source", "Target"])
    df_edge.to_csv(filename, index=False, sep='\t')
    return df_edge


edgelist_path = os.path.join(na_folder, "edgelist_completa.txt")
edgelist_completa = generate_edgelist(df[selected_genes], selected_genes, edgelist_path)


group_edgelists = {}
for i, age in enumerate(age_groups, start=1):
    group_data = df[df['AGE'] == age][selected_genes]
    file_path = os.path.join(na_folder, f"gruppo{i}.txt")
    edgelist = generate_edgelist(group_data, selected_genes, file_path)
    group_edgelists[f"gruppo{i}"] = edgelist


def plot_network(edgelist_df, title, output_file):
    G = nx.Graph()
    G.add_edges_from(edgelist_df.values)
    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(G, seed=42)
    nx.draw(G, pos, with_labels=True, node_size=800, font_size=10, node_color='skyblue', edge_color='gray')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()


plot_network(edgelist_completa, "Rete Completa", os.path.join(na_folder, "rete_completa.png"))


for i, (group_name, edgelist_df) in enumerate(group_edgelists.items(), start=1):
    plot_network(edgelist_df, f"Rete {group_name}", os.path.join(na_folder, f"rete_{group_name}.png"))




import requests
import os


def map_gene_to_string_id(gene_name, species=9606):
    url = "https://string-db.org/api/json/get_string_ids"
    params = {
        "identifiers": gene_name,
        "species": species
    }
    response = requests.get(url, params=params)
    if response.ok:
        data = response.json()
        if data:
            return data[0]['stringId'], data[0]['preferredName']
    return None, gene_name  


def get_go_bp_annotations_full(string_id):
    url = "https://string-db.org/api/json/annotation"
    params = {
        "identifiers": string_id,
        "species": 9606,
        "annotation": "GO:BP"
    }
    response = requests.get(url, params=params)
    if response.ok:
        return [(entry['term'], entry['description']) for entry in response.json()]
    else:
        return []


go_output_path = os.path.join(DESKTOP_PATH, "go_annotations.txt")


with open(go_output_path, 'w') as f:
    f.write("Gene\tMapped_Name\tGO_Biological_Process_Annotations\n")
    for gene in selected_genes:
        string_id, mapped_name = map_gene_to_string_id(gene)
        if string_id:
            annotations = get_go_bp_annotations_full(string_id)
            if annotations:
                annotation_text = "; ".join([f"{go_id}, {desc}" for go_id, desc in annotations])
            else:
                annotation_text = "N/A"
        else:
            annotation_text = "ID non trovato"
        f.write(f"{gene}\t{mapped_name}\t{annotation_text}\n")



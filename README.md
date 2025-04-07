# AGeNA
# Gene Expression Trend & Network Analysis with GO Annotations

This Python script performs a full gene expression analysis pipeline starting from GTEx-like gene expression data. It detects genes with age-related trends, performs statistical tests, builds interaction networks, and retrieves Gene Ontology Biological Process (GO:BP) annotations using STRING.

---

##  What the Script Does

1. **Input**: Reads a gene expression file (Excel format), with:
   - Rows = samples
   - Columns = gene expression values
   - Additional metadata columns: `AGE`, `SEX`, `CTRL`, `Descriptio`

2. **Expression Trend Analysis**:
   - Groups samples by `AGE`: 20–29, 30–39, ..., 70–79
   - For each gene, computes average expression per group
   - Identifies genes with strictly increasing or decreasing trends across age

3. **Statistical Tests**:
   - Pearson correlation test (age vs expression)
   - Kruskal-Wallis test (differences across age groups)

4. **Network Analysis**:
   - Computes pairwise Wilcoxon test among selected genes
   - If p-value > 0.05, an edge is added (non-significant difference)
   - Builds edge list for entire dataset and per age group
   - Saves visual network plots using `networkx` and `matplotlib`

5. **GO Annotation with STRING**:
   - Maps genes to STRING identifiers
   - Retrieves GO:BP terms (e.g., "axon guidance", "neuron differentiation")
   - Outputs gene → GO term associations with both GO ID and description

6. **Repeat for Both Sexes**:
   - Entire pipeline is executed separately for `SEX=1` (males) and `SEX=2` (females)

---

##  Input Format

- **File**: `f1.xlsx`
- **Location**: Desktop (`/Users/mariannamilano/Desktop/f1.xlsx`)
- **Required columns**:
  - `AGE`: must contain age groups as strings (e.g., `20-29`)
  - `SEX`: 1 = Male, 2 = Female
  - Gene columns: expression values (e.g., `LRP8`, `BACE1`, etc.)

---

##  Output Files

All results are saved to your Desktop:

| File/Folder | Description |
|-------------|-------------|
| `risultato_1.txt`, `risultato_2.txt` | Genes with increasing/decreasing trends + p-values (males/females) |
| `NA_1/`, `NA_2/` | Network Analysis output for males/females |
| `edgelist_completa.txt` | Global gene-gene edges based on Wilcoxon test |
| `gruppo1.txt` → `gruppo6.txt` | Edge lists for each age group |
| `rete_completa.png` | Network graph for all samples |
| `rete_gruppoX.png` | Network graphs by age group |
| `go_annotations.txt` | GO:BP annotations from STRING (GO ID + description) |

---

##  Requirements

You can install the required packages with:

```bash
pip install pandas scipy matplotlib networkx openpyxl requests

import polars as pl

# Tools in the row/column order of figure b (this IS the ranking:
# top/left = highest mean, bottom/right = lowest -> gives the red/blue split)
tools = [
    "CPT-1", "BayesDel", "REVEL", "ClinPred", "AlphaMissense",
    "CADD Raw", "PolyPhen2", "ESM1v", "GPN-MSA", "Vertebrate PhyloP",
]

# Symmetric significance, upper-triangle (i<j) read off the image.
# '?' marks cells where the two mirrored halves disagreed -> verify against source.
sig_upper = {
    ("CPT-1", "BayesDel"): "*",
    ("CPT-1", "REVEL"): "*",
    ("CPT-1", "ClinPred"): "*",
    ("CPT-1", "AlphaMissense"): "",
    ("CPT-1", "CADD Raw"): "***",          # ? (light red cell, but ***-starred)
    ("CPT-1", "PolyPhen2"): "***",
    ("CPT-1", "ESM1v"): "***",
    ("CPT-1", "GPN-MSA"): "***",
    ("CPT-1", "Vertebrate PhyloP"): "***",
    ("BayesDel", "REVEL"): "",
    ("BayesDel", "ClinPred"): "",
    ("BayesDel", "AlphaMissense"): "",
    ("BayesDel", "CADD Raw"): "",
    ("BayesDel", "PolyPhen2"): "***",       # ?
    ("BayesDel", "ESM1v"): "**",
    ("BayesDel", "GPN-MSA"): "***",
    ("BayesDel", "Vertebrate PhyloP"): "***",
    ("REVEL", "ClinPred"): "",
    ("REVEL", "AlphaMissense"): "",
    ("REVEL", "CADD Raw"): "*",
    ("REVEL", "PolyPhen2"): "***",
    ("REVEL", "ESM1v"): "**",
    ("REVEL", "GPN-MSA"): "***",
    ("REVEL", "Vertebrate PhyloP"): "***",
    ("ClinPred", "AlphaMissense"): "",
    ("ClinPred", "CADD Raw"): "*",
    ("ClinPred", "PolyPhen2"): "***",
    ("ClinPred", "ESM1v"): "**",
    ("ClinPred", "GPN-MSA"): "***",
    ("ClinPred", "Vertebrate PhyloP"): "***",
    ("AlphaMissense", "CADD Raw"): "*",     # ?
    ("AlphaMissense", "PolyPhen2"): "***",
    ("AlphaMissense", "ESM1v"): "**",       # ?
    ("AlphaMissense", "GPN-MSA"): "***",
    ("AlphaMissense", "Vertebrate PhyloP"): "***",
    ("CADD Raw", "PolyPhen2"): "***",       # ?
    ("CADD Raw", "ESM1v"): "",              # confirmed: no stars
    ("CADD Raw", "GPN-MSA"): "",
    ("CADD Raw", "Vertebrate PhyloP"): "***",
    ("PolyPhen2", "ESM1v"): "",
    ("PolyPhen2", "GPN-MSA"): "*",
    ("PolyPhen2", "Vertebrate PhyloP"): "***",
    ("ESM1v", "GPN-MSA"): "*",              # ?
    ("ESM1v", "Vertebrate PhyloP"): "***",
    ("GPN-MSA", "Vertebrate PhyloP"): "***",
}

def get_heat_df() -> pl.DataFrame:
    """Long-format tool-vs-tool table read off figure b.

    Columns: Tool_X, Tool_Y, mean_diff (+1 red / -1 blue / 0 diagonal), sig.
    """
    rank = {t: i for i, t in enumerate(tools)}
    rows = []
    for x in tools:
        for y in tools:
            if x == y:
                sig, direction = "", 0
            else:
                key = (x, y) if rank[x] < rank[y] else (y, x)
                sig = sig_upper[key]
                direction = -1 if rank[x] < rank[y] else 1   # +1 red (X-Y>0), -1 blue
            rows.append({
                "Tool_X": x, 
                "Tool_Y": y, 
                "full_data_dir": direction, 
                "full_data_sig": sig
            })
    return pl.DataFrame(rows)


if __name__ == "__main__":
    print(get_heat_df())

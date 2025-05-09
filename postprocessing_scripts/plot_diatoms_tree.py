from Bio import Phylo
from io import StringIO
import re
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def remove_node_labels(newick_str, pattern=r'\)N[0-9]+'):
    return re.sub(pattern, ')', newick_str)

def italic_labels_factory(ref_species_set):
    """
    If clade is terminal => italicize the name; for reference species => bold+italic.
    Otherwise => no label (so internal nodes have no text).
    """
    def italic_labels(clade):
        if not clade.name:
            return None
        if clade.is_terminal():
            # Replace underscores with '\ ' so LaTeX math mode shows an actual space
            label_clean = clade.name.replace('_', r'\ ')
            if clade.name in ref_species_set:
                # Bold + italic
                return rf"$\mathbf{{\mathit{{{label_clean}}}}}$"
            else:
                # Normal italic
                return rf"$\mathit{{{label_clean}}}$"
        else:
            return None  # hide internal node labels
    return italic_labels

def make_label_colors(ref_species, outgroup_dict=None):
    """
    Return a dictionary { cladeName : colorString }
    IMPORTANT: cladeName must match the actual clade.name in the tree file.
    """
    label_colors = {}
    if outgroup_dict:
        label_colors.update(outgroup_dict)
    for species in ref_species:
        # references => red
        label_colors[species] = 'red'
    return label_colors

def plot_phylogenetic_tree(
    tree,
    output_file,
    ref_species,
    base_label_colors=None,
    sublineage_annotations=None,
    sublineage_adjustments=None):
    ref_species_set = set(ref_species)
    label_func = italic_labels_factory(ref_species_set)

    fig, ax = plt.subplots(figsize=(24, 22))

    # Sort tree so outgroups (blue labels) appear at the bottom
    if base_label_colors:
        def sort_key(clade):
            if clade.name and clade.name in base_label_colors and base_label_colors[clade.name] == "blue":
                return 0  # Outgroup (blue) comes first (bottom)
            return 1  # Others come after
        tree.ladderize(reverse=True)  # first reverse (default ascending), then custom
        tree.clade.clades.sort(key=sort_key)

    # Draw the tree
    Phylo.draw(
        tree,
        label_func=label_func,
        label_colors=base_label_colors,
        branch_labels=lambda c: c.branch_length,
        do_show=False,
        axes=ax
    )

    # Remove axes/spines
    ax.set_xlabel(None)
    ax.set_ylabel(None)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ["top", "right", "bottom", "left"]:
        ax.spines[spine].set_visible(False)

    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # (Optional) Sublineage rectangles
    if sublineage_annotations:
        sublineage_adjustments = sublineage_adjustments or {}
        short_factor = 0.95
        new_x_max = x_min + short_factor * (x_max - x_min)

        pastel_palette = [
            "#FED9A6",  "#DECBE4", "#CCEBC5", "#FBB4AE","#FED9A6","#aec6cf","#CCEBC5"
               
        ]

        for i, (label, start, end) in enumerate(sublineage_annotations):
            if label in sublineage_adjustments:
                top_adj, bot_adj = sublineage_adjustments[label]
            else:
                top_adj, bot_adj = (0, 0)

            y_low = start + bot_adj
            y_high = end + top_adj
            color = pastel_palette[i % len(pastel_palette)]
            rect = Rectangle(
                (x_min, y_low),
                (new_x_max - x_min),
                (y_high - y_low),
                facecolor=color,
                edgecolor=None,
                alpha=0.6,
                zorder=-1
            )
            ax.add_patch(rect)

        # Sublineage labels
        offset_x = 0.02 * (x_max - x_min)
        label_x = (x_min + short_factor*(x_max - x_min)) + offset_x
        for (label, start, end) in sublineage_annotations:
            if sublineage_adjustments and label in sublineage_adjustments:
                top_adj, bot_adj = sublineage_adjustments[label]
            else:
                top_adj, bot_adj = (0, 0)
            y_low = start + bot_adj
            y_high = end + top_adj

            ax.text(
                label_x,
                0.5 * (y_low + y_high),
                label,
                va='center',
                ha='left',
                fontsize=14,
                color='black'
            )

    # Scale bar (updated to 0.5)
    scale_bar_len = 0.5  # changed from 1.0
    x_margin = 0.10 * (x_max - x_min)
    y_margin = 0.02 * (y_max - y_min)
    bar_x_start = x_min + x_margin
    bar_x_end = bar_x_start + scale_bar_len
    bar_y = y_min + y_margin

    ax.plot([bar_x_start, bar_x_end], [bar_y, bar_y], color='black', lw=1, zorder=2)
    cap_size = 0.003 * (y_max - y_min)
    ax.plot([bar_x_start, bar_x_start], [bar_y - cap_size, bar_y + cap_size], color='black', lw=2, zorder=2)
    ax.plot([bar_x_end, bar_x_end], [bar_y - cap_size, bar_y + cap_size], color='black', lw=2, zorder=2)

    ax.text(
        0.5 * (bar_x_start + bar_x_end),
        bar_y + 2 * cap_size,
        "0.5",  # changed label from "1.0" to "0.5"
        ha='center',
        va='bottom',
        fontsize=12,
        zorder=2
    )

    # Shrink & shift numeric branch-length labels
    numeric_pattern = re.compile(r'^[0-9]+(\.[0-9]+)?$')
    x_shift = -0.02 * (x_max - x_min)
    for txt_obj in ax.texts:
        txt = txt_obj.get_text().strip()
        if numeric_pattern.match(txt):
            txt_obj.set_fontsize(8)
            x_old, y_old = txt_obj.get_position()
            txt_obj.set_position((x_old + x_shift, y_old))
        else:
            # This is the TAXON label, increase its font size for better visibility
            txt_obj.set_fontsize(14)  # Adjust taxon label font size
    fig.savefig(output_file, format='jpeg', dpi=300, bbox_inches='tight')

## Coloring is regrettably not working... we tried.

# 1) Reference species => set them to red
ref_species = [
    "Chaetoceros_tenuissimus",
    "Cylindrotheca_closterium",
    "Fragilaria_crotonensis",
    "Mayamaea_pseudoterrestris",
    "Nitzschia_inconspicua",
    "Phaeodactylum_tricornutum",
    "Pseudo-nitzschia_multistriata",
    "Seminavis_robusta",
    "Thalassiosira_pseudonana"
]

# 2) Outgroup => blue
outgroup = {
    "Phytophthora_cinnamomi": "blue",
    "Bremia_lactucae": "blue",
    "Phytophthora_infestans": "blue",
    "Phytophthora_nicotianae": "blue",
    "Phytophthora_ramorum": "blue",
    "Phytophthora_sojae": "blue"
}

# 3) Sublineage intervals etc.
sublineage_annotations = [
    ("Coscinodiscophyceae", 0, 34.5),
    ("Bacillariophyceae", 34.5, 50.5),
    ("Fragilariophyceae", 50.5, 55.5),
    ("Mediophyceae ", 55.5,56.5),
    ("Coscinodiscophyceae", 56.5, 58.5),
    ("Outgroup: Oomycota",      58.5, 65.5)
]

sublineage_adjustments = {
    "Outgroup: Oomycota":        (0,0),
    "Coscinodiscophyceae":        (0,0),
    "Bacillariophyceae":        (0,0),
    "Fragilariophyceae":        (0,0),
    "Mediophyceae ":        (0,0),
    "Fragilariophyceae":        (0,0),
    "Coscinodiscophyceae":        (0,0)
    }


output_file = "Figure4_tree_with_sublineages.jpeg"

tree = Phylo.read("phytophtora_cleaned.nwk", "newick")

# OPTIONAL: remove node labels from Newick text
# newick_str = remove_node_labels(tree.format("newick"))
# tree = Phylo.read(StringIO(newick_str), "newick")

base_label_colors = make_label_colors(ref_species, outgroup)

plot_phylogenetic_tree(
    tree=tree,
    output_file=output_file,
    ref_species=ref_species,
    base_label_colors=base_label_colors,
    sublineage_annotations=sublineage_annotations,
    sublineage_adjustments=sublineage_adjustments
)

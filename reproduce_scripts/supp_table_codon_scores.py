import sys

def make_heatmap_tex_table(in_fname, out_fname):
    in_file = open(in_fname, "r")
    out_file = open(out_fname, "w")
    dt = []
    for line in in_file:
        line = line.strip().split()
        line = [float(elt) for elt in line]
        line = [round(elt, 4) for elt in line]
        dt.append(line)
    alpha = "ACGT"
    cods = [x + y + z for x in alpha for y in alpha for z in alpha]
    stops = ["TAA", "TAG", "TGA"]
    header ="codon," + ",".join([str(elt) for elt in range(-7,6)]) + "\n"
    out_file.write(header)
    for idx, cod in enumerate(cods):
        if cod not in stops:
            line = ",".join([cod] + [str(elt) for elt in dt[idx]]) + "\n"
            out_file.write(line)
    in_file.close()
    out_file.close()

if __name__ == "__main__":
    in_fname = sys.argv[1]
    out_fname = sys.argv[2]
    make_heatmap_tex_table(in_fname, out_fname)

process CIRCOS_SV {
    tag "$meta.id"
    label 'process_medium'
    
    conda 'bioconda::circos=0.69.9'
    //container "" // No container available
    
    input:
    tuple val(meta), path(vcf)
    path reference_fai
    val caller
    
    output:
    tuple val(meta), path("*.png"), emit: plot
    tuple val(meta), path("*_mqc.html"), emit: multiqc
    
    script:
    def prefix = "${meta.id}.${caller}"
        """
export REFERENCE_FAI="${reference_fai}"
export VCF_PATH="${vcf}"
export PREFIX="${prefix}"
export CALLER="${caller}"

python << 'PYCODE'
import gzip
import os

REFERENCE_FAI = os.environ["REFERENCE_FAI"]
VCF_PATH      = os.environ["VCF_PATH"]
PREFIX        = os.environ["PREFIX"]
CALLER        = os.environ["CALLER"]

# --------------------------------------
# 1) Build karyotype from reference .fai
# --------------------------------------
with open("karyotype.txt", "w") as k:
    with open(REFERENCE_FAI) as fai:
        for line in fai:
            chrom, length, *_ = line.strip().split()
            k.write(f"chr - {chrom} {chrom} 0 {length} chr{chrom}\\n")

# --------------------------------------
# 2) Parse VCF and build links.txt
# --------------------------------------
def smart_open(filename):
    with open(filename, 'rb') as f:
        magic = f.read(2)
        f.seek(0)
        if magic == b'\\x1f\\x8b':  # gzip magic number
            return gzip.open(filename, 'rt')
        else:
            return open(filename, 'r')

links = []
with smart_open(VCF_PATH) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\\t")
        info = dict(
            item.split("=") if "=" in item else (item, True)
            for item in fields[7].split(";")
        )

        svtype = info.get("SVTYPE", "")
        if svtype in ["BND", "TRA", "INV"]:
            chrom1 = fields[0]
            pos1   = fields[1]
            chrom2 = info.get("CHR2", chrom1)
            pos2   = info.get("END", pos1)
            links.append(f"{chrom1} {pos1} {pos1} {chrom2} {pos2} {pos2}\\n")

with open("links.txt", "w") as l:
    l.writelines(links)

# --------------------------------------
# 3) Write circos config
# --------------------------------------
with open("circos.conf", "w") as conf:
    conf.write('''
karyotype = karyotype.txt

<links>
<link>
file          = links.txt
radius        = 0.95r
bezier_radius = 0.1r
thickness     = 2
</link>
</links>

<image>
<<include etc/image.conf>>
</image>

<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
''')
PYCODE

# --------------------------------------
# 4) Run circos (now in conda env, so in PATH)
# --------------------------------------
echo "DEBUG: PATH=\$PATH"
echo "DEBUG: which circos:"
which circos || echo "DEBUG: circos NOT FOUND in PATH"
echo "DEBUG: conda env:"
conda info --envs || true

circos -conf circos.conf -outputfile ${prefix}.circos.png

# --------------------------------------
# 5) MultiQC custom content HTML
# --------------------------------------
cat > ${prefix}_circos_mqc.html <<EOF
<div id="circos_plot">
<h3>Circos Plot - ${caller}</h3>
<img src="${prefix}.circos.png" width="800px">
</div>
EOF
    """
}
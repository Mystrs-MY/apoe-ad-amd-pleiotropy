BEGIN {
    FS = OFS = "\t"
    alpha_count = 0
    pqtl_count = 0
}

NR == 1 { next }

{
    rsfield = $4
    gsub(/;/, ",", rsfield)
    alpha_variant = ""
    if (("," rsfield ",") ~ /,rs429358,/) alpha_variant = "rs429358"
    if (("," rsfield ",") ~ /,rs7412,/) alpha_variant = "rs7412"

    if (alpha_variant != "") {
        print normalization, analysis_role, gene, seqid, uniprot, source_file, alpha_variant, $0 >> alpha_out
        alpha_count++
    }

    if ($8 ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ && ($8 + 0) < p_primary) {
        strict = (($8 + 0) < p_strict) ? "true" : "false"
        print normalization, analysis_role, gene, seqid, uniprot, source_file, strict, $0 >> pqtl_out
        pqtl_count++
    }
}

END {
    print normalization, gene, seqid, alpha_count, pqtl_count
}

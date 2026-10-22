# Genes whose CeNGEN expression values are unreliable because they were
# overexpressed to build the transgenic sorting strains (warning on
# cengen.org). They are excluded from the clustering features in addition to
# the 199-gene stress list unless the clustering runner gets --keep-transgenes.
transgene_ids <- c(
  "lin-15A"  = "WBGene00023498",
  "lin-15B"  = "WBGene00023497",
  "rol-6"    = "WBGene00004397",
  "unc-119"  = "WBGene00006843",
  "pha-1"    = "WBGene00004010",
  "unc-54"   = "WBGene00006789",
  "eat-4"    = "WBGene00001135",
  "dpy-20"   = "WBGene00001079",
  "unc-47"   = "WBGene00006783",
  "cho-1"    = "WBGene00000501",
  "unc-53"   = "WBGene00006788",
  "gcy-35"   = "WBGene00001555",
  "C30A5.16" = "WBGene00206533",
  "saeg-2"   = "WBGene00011964",
  "F38B6.2"  = "WBGene00018172",
  "C30F8.3"  = "WBGene00016259",
  "cex-1"    = "WBGene00023407"
)


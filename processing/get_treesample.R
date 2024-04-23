library(ape)
library(fpeek)

fname = "analysis/m_relaxed_262E8_empirical/BEAST_2_m-papua_m.trees"

n_lines = peek_count_lines(fname)

n_lines_head = 1038
phylo_head = readLines(fname, n_lines_head)

con = file(fname)
open(con)
phylo_tail = scan(con,
                  skip = 50000,
                  what = "char(0)",
                  sep = "\n",
                  quiet = TRUE)

phylo = c(phylo_head, phylo_tail, "End;")

tail(phylo)

tmp = writeLines(phylo, "analysis/m_relaxed_262E8_empirical/tmp.tree")

tree = read.nexus(
  file = "analysis/m_relaxed_262E8_empirical/tmp.tree"
  )

write.nexus(tree, file = "analysis/m_relaxed_262E8_empirical/m_962trees.tree")

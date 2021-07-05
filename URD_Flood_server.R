#!/usr/bin/env Rscript

flood.result <- floodPseudotime(object.merged.urd, root.cells=root.cells, n=10, minimum.cells.flooded=2, verbose=T)

saveRDS(flood.result, file=tempfile(pattern = "flood-", tmpdir = ".", fileext = ".rds"))

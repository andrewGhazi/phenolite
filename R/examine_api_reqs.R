og = fread("~/projects/phenolite/data/adj_2026-09-24T21-35-31-668Z.csv")

x = resp$body |> rawToChar() |> fromJSON() 

dat = x$data |> qDT()

dat$rows[[1]]$fields

og |> sbt(species_2 %like% "Asclepias tuber")

# I think I can see dat, species name, flower_phenophase, and leaf_phenophase.
# That should theoretically? be everything needed.

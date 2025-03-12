# test different color schemes


#Hiroshige
species_colors <- c(
  piab = "#2A5674",
  fasy = "#4682B4",
  quro = "#7EA6CE",
  pisy = "#A7C6DD",
  soau = "#C7DDE5",
  acps = "#F4E8C8",
  potr = "#E8C099",
  abal = "#E28F60",
  besp = "#D0593E",
  lade = "#9C2E2E"
)

# archambeau
species_colors <- c(
  piab = "#1F2041",
  fasy = "#4A2670",
  quro = "#7C3688",
  pisy = "#A64E97",
  soau = "#C7749C",
  acps = "#E39F9E",
  potr = "#F3C8A7",
  abal = "#F5E3B4",
  besp = "#E7F0BC",
  lade = "#C1D982"
)

# vanGogh3
species_colors <- c(
  piab = "#355F3A",
  fasy = "#4A7B41",
  quro = "#6E984B",
  pisy = "#92B45B",
  soau = "#B1CC69",
  acps = "#CCE07B",
  potr = "#E3EE91",
  abal = "#F1F6A4",
  besp = "#F8F9C3",
  lade = "#FAFBE6"
)

library(MetBrewer)
MetBrewer::colorblind_palettes

MetBrewer::colorblind.friendly("Hiroshige")
met.brewer(name="Hiroshige", n=10, type="discrete")
met.brewer(name="VanGogh3", n=10, type="continuous")   # greens
met.brewer(name="VanGogh1", n=10, type="continuous")
my_cols <- met.brewer(name="Derain", n=10, type="continuous")

met.brewer(name="Archambault", n=10, type="continuous")
met.brewer(name="Demuth", n=10, type="discrete")

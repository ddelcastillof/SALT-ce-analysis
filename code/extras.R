library(sf)
library(rnaturalearth)
library(rnaturalearthhires)
library(ggplot2)

# Cargar los estados (departamentos) de Perú
peru_dep <- ne_states(country = "Peru", returnclass = "sf")

# Definir color de relleno según departamento
peru_dep$fill_col <- ifelse(peru_dep$name_en == "Tumbes",
                            "#4b2e83",    # morado para Tumbes
                            "gray90")     # gris pálido para los demás

# Graficar
peru_map <- ggplot(peru_dep) +
  geom_sf(aes(fill = fill_col),
          color  = "#4b2e83",  # color de borde morado
          size   = 0.4) +      # grosor de línea de borde
  scale_fill_identity() +             # usar exactamente los colores definidos
  theme_void() +                      # quitar elementos innecesarios
  theme(
    plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(peru_map,
filename = "Tumbes_Location.jpeg",
width = 25, height = 40, units = "cm",
dpi = 1200,
path = "output/figs")

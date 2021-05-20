library(sf)
library(tidygraph)
library(igraph)
library(dplyr)
library(tibble)
library(ggplot2)
library(units)
library(tmap)
library(osmdata)
library(rgrass7)
library(link2GI)
library(nabor)
library(geojsonio)


a <- geojson_read("data/Kalikiri_RoadNetwork.geojson",  what = "sp")%>% st_as_sf()%>% st_transform(st_crs(b))
# hem_kkiri = st_as_sf(spdf)
# spdf = st_transform(j,st_crs(kalikiri_center))

# kalikiri <- opq(bbox =  c(78.70315, 13.56041, 78.87547, 13.86580)) %>% 
#   add_osm_feature(key = 'highway') %>% 
#   osmdata_sf() %>% 
#   osm_poly2line()
# 
# b <- kalikiri$osm_lines %>% 
#   select(highway)


a = st_cast(a,"LINESTRING")
count = 0
for(i in 1:nrow(a)){
  
  if(is.na(a$road_name[i])){
    a$road_name[i] = count
    count = count + 1
  }
}
x = a
edges <- x %>%
  mutate(edgeID = c(1:n()))

nodes <- edges %>%
  st_coordinates() %>%
  as_tibble() %>%
  rename(edgeID = L1) %>%
  group_by(edgeID) %>%
  slice(c(1, n())) %>%
  ungroup() %>%
  mutate(start_end = rep(c('start', 'end'), times = n() / 2)) %>%
  mutate(xy = paste(.$X, .$Y)) %>%
  mutate(nodeID = group_indices(., factor(xy, levels = unique(xy)))) %>%
  select(-xy)

source_nodes <- nodes %>%
  filter(start_end == 'start') %>%
  pull(nodeID)

target_nodes <- nodes %>%
  filter(start_end == 'end') %>%
  pull(nodeID)

edges = edges %>%
  mutate(from = source_nodes, to = target_nodes)

nodes <- nodes %>%
  distinct(nodeID, .keep_all = TRUE) %>%
  select(-c(edgeID, start_end)) %>%
  st_as_sf(coords = c('X', 'Y')) %>%
  st_set_crs(st_crs(edges))

tbl_graph(nodes = nodes,
          edges = as_tibble(edges),
          directed = directed)



h = tbl_graph(nodes = nodes,
              edges = as_tibble(edges),
              directed = directed) %>%
              activate(edges) %>%
              mutate(length = st_length(geometry))

h = h%>% activate(edges) %>%
  as_tibble() %>%
  st_as_sf() %>%
  group_by(road_name) %>%
  summarise(length = sum(length))


graph = h %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree()) %>%
  mutate(betweenness = centrality_betweenness(weights = length)) %>%
  activate(edges) %>%
  mutate(betweenness = centrality_edge_betweenness(weights = length))






ggplot() +
  geom_sf(data = graph %>% activate(edges) %>% as_tibble() %>% st_as_sf(), col = 'grey50') + 
  geom_sf(data = graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf(), aes(col = (betweenness), size = betweenness)) +
  scale_colour_viridis_c(option = 'inferno') +
  scale_size_continuous(range = c(0,6))















# "+proj=longlat +datum=WGS84"

spdf <- geojson_read("data/Kalikiri_RoadNetwork.geojson",  what = "sp")

hem_kkiri = st_as_sf(spdf)

spdf = st_transform(j,st_crs(kalikiri_center))



if(!"remotes" %in% installed.packages()) {
  install.packages("remotes")
}

cran_pkgs = c(
  "sf",
  "tidygraph",
  "igraph",
  "osmdata",
  # "dplyr",
  "tibble",
  "ggplot2",
  "units",
  "tmap",
  "rgrass7",
  "link2GI",
  "nabor"
)

remotes::install_cran(cran_pkgs)

# x 78.70315 78.87547
# y 13.56041 13.86580

kalikiri <- opq(bbox =  c(78.70315, 13.56041, 78.87547, 13.86580)) %>% 
  add_osm_feature(key = 'highway') %>% 
  osmdata_sf() %>% 
  osm_poly2line()

kalikiri_center <- kalikiri$osm_lines %>% 
  select(highway)

ggplot(data = kalikiri_center) + geom_sf()

ggplot(data = a) + geom_sf()


edges <- kalikiri_center %>%
  mutate(edgeID = c(1:n()))

edges



nodes <- edges %>%
  st_coordinates() %>%
  as_tibble() %>%
  rename(edgeID = L1) %>%
  group_by(edgeID) %>%
  slice(c(1, n())) %>%
  ungroup() %>%
  mutate(start_end = rep(c('start', 'end'), times = n()/2))

nodes


nodes <- nodes %>%
  mutate(xy = paste(.$X, .$Y)) %>%
  mutate(nodeID = group_indices(., factor(xy, levels = unique(xy)))) %>%
  select(-xy)

nodes



source_nodes <- nodes %>%
  filter(start_end == 'start') %>%
  pull(nodeID)

target_nodes <- nodes %>%
  filter(start_end == 'end') %>%
  pull(nodeID)

edges = edges %>%
  mutate(from = source_nodes, to = target_nodes)

edges

nodes <- nodes %>%
  distinct(nodeID, .keep_all = TRUE) %>%
  select(-c(edgeID, start_end)) %>%
  st_as_sf(coords = c('X', 'Y')) %>%
  st_set_crs(st_crs(edges))

nodes


graph = tbl_graph(nodes = nodes, edges = as_tibble(edges), directed = FALSE)

graph



graph <- graph %>%
  activate(edges) %>%
  mutate(length = st_length(geometry))

graph


ggplot() +
  geom_sf(data = graph %>% activate(edges) %>% as_tibble() %>% st_as_sf()) + 
  geom_sf(data = graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf(), size = 1)



library(tmap)
tmap_mode('view')

tm_shape(graph %>% activate(edges) %>% as_tibble() %>% st_as_sf()) +
  tm_lines() +
  tm_shape(graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf()) +
  tm_dots() +
  tmap_options(basemaps = 'OpenStreetMap')



graph <- graph %>%
  activate(nodes) %>%
  mutate(degree = centrality_degree()) %>%
  mutate(betweenness = centrality_betweenness(weights = length)) %>%
  activate(edges) %>%
  mutate(betweenness = centrality_edge_betweenness(weights = length))

graph

ggplot() +
  geom_sf(data = graph %>% activate(edges) %>% as_tibble() %>% st_as_sf(), col = 'grey50') + 
  geom_sf(data = graph %>% activate(nodes) %>% as_tibble() %>% st_as_sf(), aes(col = betweenness, size = betweenness)) +
  scale_colour_viridis_c(option = 'inferno') +
  scale_size_continuous(range = c(0,4))

ggplot() +
  geom_sf(data = graph %>% activate(edges) %>% as_tibble() %>% st_as_sf(), aes(col = betweenness, size = betweenness)) +
  scale_colour_viridis_c(option = 'inferno') +
  scale_size_continuous(range = c(0,4))


library(ggmap)

# For google map, you have to give the center of the window you are looking at.
# Possibility for the map type argument: terrain / satellite / roadmap / hybrid

library(ggmap)
# get the map info
map <- get_googlemap()

# 13.651884320475968, 78.80673318475489
map <-  get_openstreetmap(bbox = c(13.56041,78.70315,13.86580,78.87547))
map <- get_stamenmap( bbox = c(left = 78.70315, bottom = 13.56041, right = 78.87547, top = 13.86580), zoom = 4, maptype = "terrain")

ggmap(map)

ndvi = readRDS("ndvi_kalikiri.rds")

xx = ndvi[[1]]
# xx = raster(xx)
library(broom)
library(raster)
xxdf = as.data.frame(xx)

res_df <- as.data.frame((xx),
                        xy = TRUE, na.rm = TRUE)


res_df <- as.data.frame(ndvi[[1]],
                        xy = TRUE, na.rm = TRUE)



ggplot(data = res_df, aes(x = x, y = y)) +
  geom_raster(aes(fill = X2011.12.19)) +
  scale_fill_viridis_c()+
  geom_sf(data = graph %>% activate(edges) %>% as_tibble() %>% st_as_sf(), aes(col = betweenness, size = betweenness)) +
  scale_colour_viridis_c(option = 'inferno') +
  scale_size_continuous(range = c(0,4))


# ggplot(res_df, aes(x = x, y = y)) +
  ggplot()+ geom_raster(data = res_df,aes(x=x,y=y,fill = X2011.12.19)) +
    scale_fill_viridis_c()+
    geom_sf(data = graph %>% activate(edges) %>% as_tibble() %>% st_as_sf(), aes(col = betweenness, size = betweenness)) +
    scale_colour_viridis_c(option = 'inferno') +
    scale_size_continuous(range = c(0,4))
  

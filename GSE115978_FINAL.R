#Dataset del repositori GEO.
#GSE115978.

#En aquest dataset, ens parla dels Immune checkpoint inhibitors (ICI), que produeixen respostes d'alta duració en pacients de melanoma. En aquest dataset s'ha utilitzat single-cell RNA-seq de 31 tumors de melanoma, i diferents mètodes computacionals per tal de veure estats cel·lulars malignes que promousen evasió immunitària.

#Els resultats extrets, indiquen que existeix un programa de resistència expressat per les cèl·lules malignes associat a l'exclusió de cèl·lules T i l'evasió immunitària. També descobreixen els autors que la inhibició de CDK4/6 podria ser una estratègia terapèutica efectiva per superar aquesta resistència.

#La majoria de les cèl·lules individuals han estat extretes de resseccions tumorals fresques. Aïllades cèl·lules immune i no immune mitjançant FACS (Fluorescence-activated cell sorting) basades en la tinció de CD45. I modelades amb un protocol SMART-Seq2.




#Ens col·loquem al directori on tenim els fitxers.


setwd("~/Github/test")




#Descarreguem els diferents fitxers

#annotations:
#El fitxer anotacions, ens servirà per filtrar, classificar les dades de scRNA-seq, basant-nos en diferents propietats de les diferents cèl·lules, com l'estat de la malaltia, el tipus de cèl·lules o el perfil d'expressió gènica.



library(utils)


#tpm: tenim l'expressió en una matriu.


download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978_cell.annotations.csv.gz", "GSE115978_cell.annotations.csv.gz")
annotations <- read.csv("GSE115978_cell.annotations.csv.gz", header=T)


library(curl)
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978_tpm.csv.gz"
destfile <- "GSE115978_tpm.csv.gz"

handle <- new_handle(timeout = 500)
curl_download(url, destfile, handle = handle)


#counts <- read.csv("GSE115978_counts.csv", header=T, row.names=1)
tpm <- read.csv(destfile, header=T, row.names=1)


#Carreguem les diferents llibreries que seran necessàries 



library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(tidyr)




#Creem l'objecte Seurat.
#counts=li posem la matriu expressió. On hi ha com a columnes lès cèl·lules, i com a noms de les files els diferents gens.
#Ara en aquest pas no filtrem res, perquè ho necessitem per afegir les anotacions.


tres <- CreateSeuratObject(counts=tpm)


#Ens canvia els "_" dels noms d'identificació de les cèl·lules per "-".


#A l'opció @meta.data, és on hi tenim l'informació d'interès.

#nCount: nombre de contejos per cada cèl·lula del dataset. És una mesura de l'expressió total del s gens en cada cèl·lula.
#nFeature: nombre de gens únics que són detectats en cada cèl·lula. És una mesura de la diversitat transcripcional de cada cèl·lula.
#Les variables de la metadata són importants a considerar quan es realitzen anàlisis "downstream"

#View(tres@meta.data)



#Reemplaçem totes les cèl·lules que no estan anotades amb un NA. 


annotations$cell.types <- gsub("\\?", "NA", annotations$cell.types)




#Li podem afegir les diferents anotacions del document "annotations" a la metadata de l'objecte Seurat. Per tal de tenir les diferents cèl·lules identificades.
#Cell.types= tipus de cèl·lules diferents, ja identificades segons el seu grup.
 # Mal: Que són les malignes.
 # I els diferents tipus.
  
#Treatment.group= 
 # no tractament (naive)
  #post tractament (després)
#Samples= 
  #Diferencia els pacients
#Cohort=
  #New: Es refereix a les dades obtingudes de tumors primaris de melanoma i metàstasis de ganglis linfàtics. 
 # Tirosh: Es refereix a les dades obtingudes de tumors primaris de melanoma i gànglis linfàtics. TIROSH

#La informació de les cèl·lules, que ens pot servir d'interès. S'HA DE VEURE.
#exemples:
#cy79_p1_CD45_neg_PDL1_pos_AS_C1_R1_G01_S73_comb
#CY88_3_D02_S614_comb
#Veure si podem filtrar per obtenir alguna altra informació.

#Afegim les diferents dades de les anotacions a la metadata, per poder filtrar segons la diferent informació donada.


tres@meta.data$cell.types <- annotations$cell.types
tres@meta.data$treatment.group <- annotations$treatment.group
tres@meta.data$samples <- annotations$samples
tres@meta.data$cohort <- annotations$Cohort




#Comencem el control de qualitat, per determinar que es comporti el dataset com volem.
#Borrem els NA.


tres <- na.omit(tres)


#Filtrem per les dades que tinguin un nombre mínim de genes detectats en cada cèlula, per tal de treure el soroll de les que tinguin un nFeatures petit.
#A vegades també es filtra per proporció de molècules provinent del ribosoma, però en aquest cas no tenim seqüències d'ARN.

tres <- subset(tres, subset = nFeature_RNA > 200)




#Normalitzem les dades amb la funció que ens aporta Seurat per fer-ho.

tres <- NormalizeData(tres)




#Identifiquem els gens més variables d'entre les cèl·lules del dataset. Aquests ens donen informació per distingir el tipus cel·lular, o els estats.
#Utilitzem el mètode vst (variance stabilizing transformation), per tal d'identificar els gens amb la major variància, aquest mètode és robuts per la majoria de datasets.
#Com a nfeatures= establim el nombre de gens que volem com a variables, per defecte, li posem els 2000, es pot ajustar segons el tamany i la complexitat del dataset. 



tres <- FindVariableFeatures(tres, selection.method = "vst", nfeatures = 2000)



#Escalem les dades per tal de normalitzar l'expressió dels gens en les cèl·lules i esborrar qualsevol variació que no volem, per tal de tenir un anàlisis més acurat.


tres <- ScaleData(tres)







#VEURE SI FUNCIONA PERQUÈ VEIG QUE No
#pbmc <- JackStraw(tres, num.replicate = 100)
#pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
#JackStrawPlot(pbmc)


#Utilitzem la funció de RunPCA per tal d'identificar les fonts de variació més importants en el dataset.


tres <- RunPCA(tres, features = VariableFeatures(object = tres))



#Els gens que trobem com a resultat positiu al realitzar el RunPCA tendeixen a ser regulats en cèl·lules per un alt valor d'aquell mateix PC. upregulated.
#En canvi els de signe negatiu downregulated.
#Aquests gens d'interès ens permetràn identificar els diferents subgrups cèl·lulars, en cas que no estiguin etiquetats.




#Representem l'ElbowPlot, per determinar les dimensions.

ElbowPlot(tres)


#Utilitzem les 20 dimensions.





#Amb la funció FindNeighbors per trobar les cèl·lules veïnes basant-nos en similituds en perfils d'expressió gèniques. Bàsicament, calcula SNN graph (shared nearest neighbor)


tres <- FindNeighbors(tres, dims = 1:20)






#Una vegada utilitzada la funció anterior, utilitzem aquesta per buscar els clusters de les cèl·lules basant-nos en l'estructura gràfica. El cluster correspon a grups de cèl·lules que són més semblants a altres parlant sobre expressió gènica.

#Realitzem el cluster amb una resolució de 0.5.


tres <- FindClusters(tres, resolution = 0.5)



#Realitzem la funció RunUMAP, per tal d'identificar els clusters en una menor dimensió


tres <- RunUMAP(tres, dims=1:20)





#Podem visualitzar els resultats dels clusters.
#Agrupem els diferents clusters (segons les diferents categories que tenim a les anotacions)

#Agrupem per les diferents columnes de la metadata.
#Utilitzem la funció DimPlot, per tal de visualitzar scRNA-seq en una dimensionalitat reduida.


cluster_samples <- DimPlot(tres, reduction="umap", group.by="samples")
cluster_cell.types <- DimPlot(tres, reduction="umap", group.by="cell.types")
cluster_treatment.group <- DimPlot(tres, reduction="umap", group.by="treatment.group")
cluster_cohort <- DimPlot(tres, reduction="umap", group.by="cohort")
cluster_resolució_0.5 <- DimPlot(tres, reduction = "umap", group.by = "RNA_snn_res.0.5")
cluster_orig.ident <- DimPlot(tres, reduction = "umap", group.by = "orig.ident")


t.cell <- subset(tres, cell.types=="T.cell")
DimPlot(t.cell, reduction="umap", group.by="treatment.group")

Mal <- subset(tres, cell.types=="Mal")
malignes.treatment <- DimPlot(Mal, reduction="umap", group.by="treatment.group")
malignes2 <- DimPlot(Mal, reduction="umap", group.by="cell.types")


post <- subset(tres,treatment.group=="post.treatment" )
naive <- subset(tres, treatment.group=="treatment.naive")
post <- DimPlot(post, reduction="umap", group.by="cell.types")
naive <- DimPlot(naive, reduction="umap", group.by="cell.types")


malignant <- subset(tres, samples %in% c("Mel78","Mel79", "Mel88", "Mel71", "Mel81","Mel80", "Mel89", "Mel194","Mel102","Mel110", "Mel103", "Mel106","Mel98", "Mel129pa"))

immune <- subset(tres, cell.types %in% c("T.CD8", "B.cell", "CAF", "Endo", "Macrophage", "NK", "T.CD4", "T.CD8", "T.Cell"))

no <- DimPlot(malignant, reduction="umap", group.by="samples")
immune.treatment <- DimPlot(immune, reduction="umap", group.by="treatment.group")
immune2 <- DimPlot(immune, reduction="umap", group.by="cell.types")


post | naive
View(tres@meta.data)



#El primer gràfic és el post, el segon és el naive. És a dir, post tractament i pre-tractament. Podem veure que les cèl·lules "Mal" que entenc que són les malignes, al post tractament les d'abaix a la dreta no hi són i n'hi ha menys.




malignes2 | immune2






malignes.treatment | immune.treatment



#Diferents cèl·lules malignes amb les immune segons el tipus i segons el tractament.


#He trobat curiós de mirar aquests. Suposu que hi ha d'haver altra informació que podem extreure, o que m'ha passat per alt en aquesta primera ullada.
#Mirar a veure quines més convinacions hi ha. Perquè se'n poden extreure moltes la veritat.



#PODEM FILTRAR PER TOTS ELS DIFERENTS TIPUS DE MELANOMA, I VEURE COM SÓN LES DIFERENTS RESPOSTES, HEM DE DECIDIR QUE ES FILTRA.







#Veure si podem fer algun FeaturePlot
#AMB ELS DIFERENTS GENS I VEURE LES REPRESENTACIONS
#AIXÒ ÉS EL QUE SERIA INTERESSANT NO???
#Veure una mica què i com ho volem.
#FeaturePlot(qqq, features = c("CD3E", "CD4", "CD8A", "FOXP3"), reduction = "umap")

#Representem els diferents clusters


cluster_samples
cluster_cell.types
cluster_treatment.group
cluster_cohort 
cluster_resolució_0.5
cluster_orig.ident



#VEURE COM ELS PODEM FILTRAR PELS DIFERENTS GRUPS DE CLUSTER
#AIXÒ HO PODEM MIRAR
#VEURE COM PODEM FER EL CLUSTER SENSE ANOTACIONS PER TAL D'IDENTIFICAR ELS DIFERENTS GRUPS ETIQUETATS.
#EN AQUEST CAS, PODREM COMPARAR.




#Ja hem acabat de fer l'anàlisi single-cell

#Ara realitzem l'anàlisi amb Monocle3.

#Amb el paquet SeuratWrappers, convertim l'objecte seurat en un nou cell dataset.
#Veure si podem fer algun canvi més o que obtenim etc.

sccc <- as.cell_data_set(tres)





#Amb colData, ens retorna les columnes i les files.

#colData(sccc)
#fData(sccc)

#Per veure els diferents nomes de les files.

#rownames(fData(sccc))[1:10]



#Li posem el nom dels gens a fData (el data frame que conté la metadata)


fData(sccc)$gene_short_name <- rownames(fData(sccc))

#counts(sccc)



#Creem una nova variable, que es diu "recreate.partition", creem una partició de totes les mostres de l'objecte "scc" en un sol grup. Creem un vector del mateix tamany que el número de cèl·lules de la matriu d'expressió.
#L'objectiu és inicialitzar una nova variable per tal d'agregar a l'objecte "sccc"


recreate.partition <- rep(1,length(sccc@colData@rownames))



#S'assignen noms als elements del vector anterior, utilitzant els noms de les files (gens). Per tal de que cada cèl·lula en la matriu d'expressió tingui una etiqueta única en la nova variable.



names(recreate.partition) <- sccc@colData@rownames


#Convertim el vector generat en un factor. Per tal de representar una variable categòrica. 


recreate.partition <- as.factor(recreate.partition)


#En aquest cas, hem assignat totes les cèl·lules a un mateix grup, que el representem amb el número 1.




#Hem fet el codi anterior per tal d'assignar la variable com a partició en el cluster UMAP. L'objecte sccc@clusters és resultat del clustering del conjunt de dades de l'expressió gènica. L'assignació de la variable recreate.partition, ens permet etiquetar els clusters amb una nova variable, ja que les cèl·lules s'han assignat al mateix grup.



sccc@clusters$UMAP$partitions <- recreate.partition



 

#Assignem a la variable "list_cluster" els valors de la columna "active.ident" de la matriu anterior de l'ojbecte Seurat.



list_cluster <- tres@active.ident


#list_cluster <- tres@meta.data$cell.types




#Assignem a la columna clusters, els valors de "list_cluster". És a dir, assignem a cada cèl·lula el valor que ens ha donat el cluster realitzat anteriorment.
#En aquest cas, ens ha donat 23 valors de cluster diferent.




sccc@clusters$UMAP$clusters <- list_cluster

#head(sccc@clusters$UMAP$clusters)






#Aquest codi, ens assigna les coordenades UMAP generades per l'objecte Seurat al nou objecte creat "sccc".

sccc@int_colData@listData$reducedDims$UMAP <- tres@reductions$umap@cell.embeddings




#Realitzem les diferents trajectòries`amb plot_cells, que s'utilitza per generar un scatter, on cada punt representa una cèl·lula individual (en espai PCA o UMAP), el color del punt representa el cluster assignat de la cèl·lula.


cluster.before.trajectory <- plot_cells(sccc, 
                                        color_cells_by = "cluster", #Les cèl·lules són colorejades a partir del cluster assignat.
                                        label_groups_by_cluster = FALSE, #Indiquem si volem etiquetar a partir del cluster o no.
                                        group_label_size = 5) + theme(legend.position="right")  #Ens indica la mida de la lletra de les etiquetes del cluster









#Utiltizem la funció per visualitzar els clusters generats a partir de la partició de les cèl·lules.


cluster.names <- plot_cells(sccc,
                            color_cells_by = "partition",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  
  theme(legend.position="right")



cluster.before.trajectory | cluster.names


#unique(tres@meta.data$redefined_cluster)


#scale_color_manual(values=c("red", "blue", "green", "maroon", "yellow", "grey", "cyon")) 



#Utilitzem la funció learn_graph per tal de realitzar la representació gràfica de les diferents relacions entre les cèl·lules dels conjunts de dades. El gràfic es construeix a partir de les similituds entre les cèl·lues mitjançant els seus perfils d'expressió gènica.

#Li diem que la partició anterior no es tindrà en compte. És a dir, s'estan considerant totes les cèl·lules en l'anàlisis sense distinció de grup o cluster. 


sccc <- learn_graph(sccc, use_partition=FALSE)






plot_cells(sccc,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_branch_points=FALSE,
           label_roots=FALSE,
           label_leaves=FALSE,
           group_label_size=5)
```



#Triem quin grup de cluster és el root.
#Hem d'escollir quin seria. NO TINC CLAR COM HO HEM DE DEFINIR, en aquest cas estem definint el cluster que té com a número assignat el 5.

sccc <- order_cells(sccc, reduction_method="UMAP", root_cells=colnames(sccc[,clusters(sccc)==5]))





#Veure com s'ordenen segons el pseudotime.
#Hem de mirar, si els diferents grups, i el que segueixen.


plot_cells(sccc,
           color_cells_by="pseudotime",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves=FALSE)





#Aquí filtrem les dades segons el tractament, si n'hi ha hagut, o si no n'hi ha hagut. Per tal de comparar les diferències de tractament.



filter <- subset(sccc,sccc@metadata$treatment.group=="post.treatment")

filter2 <- subset(sccc,sccc@metadata$treatment.group=="treatment.naive")


filter <- reduce_dimension(filter)
filter <- cluster_cells(filter, reduction_method = "UMAP")
filter <- learn_graph(filter)

#Posem com a root les malignes?
filter <- order_cells(filter, reduction_method="UMAP", root_cells=colnames(filter[,clusters(filter) %in% c("3","4")]))


tractament <- plot_cells(filter, color_cells_by = "pseudotime", 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves=FALSE)


filter2 <- reduce_dimension(filter2)
filter2 <- cluster_cells(filter2, reduction_method = "UMAP")
filter2 <- learn_graph(filter2)

#Posem com a root les malignes?
filter2 <- order_cells(filter2, reduction_method="UMAP", root_cells=colnames(filter2[,clusters(filter2) %in% c("3","4")]))


no_tractament <- plot_cells(filter2, color_cells_by = "pseudotime", 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves=FALSE)


#FILTRAT TRACTAMENT I NO TRACTAMENT EN CELL.TYPES



tractament | no_tractament




#Podem observar diferències! 

#Ara hem d'establir roots, i anar perfilant prim a veure que més podem fer ??
  
  
  
 # PROVA DE QUE ES PODEN FER SUBSETS AMB L'OBJECTE MONOCLE. SREIA INTERESSANT VEURE QUINES COMPARACIONS PODEM FER, I COM LES PODEM ANAR FENT.

#Com a roots, hem establerts els valors del cluster que sabem quins son.

filter <- subset(sccc,sccc@metadata$treatment.group=="post.treatment" & sccc@metadata$cell.types %in% c("T.CD8", "T.CD4", "T.cell"))



filter <- reduce_dimension(filter)
filter <- cluster_cells(filter, reduction_method = "UMAP")
filter <- learn_graph(filter)

#filter <- order_cells(filter, reduction_method = "UMAP")

filter <- order_cells(filter, reduction_method="UMAP", root_cells=colnames(filter[,clusters(filter) %in% c("0","1")]))


mal_post <- plot_cells(filter, color_cells_by = "pseudotime", 
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves=FALSE)

filter2 <- subset(sccc,sccc@metadata$treatment.group=="treatment.naive" & sccc@metadata$cell.types %in% c("T.CD8", "T.CD4", "T.cell"))




filter2 <- reduce_dimension(filter2)
filter2 <- cluster_cells(filter2, reduction_method = "UMAP")
filter2 <- learn_graph(filter2)

#filter2 <- order_cells(filter2, reduction_method = "UMAP")

#filter2 <- order_cells(filter2, reduction_method="UMAP", root_cells=colnames(filter2[,clusters(filter2) %in% c("3","4")]))
filter2 <- order_cells(filter2, reduction_method="UMAP", root_cells=colnames(filter2[,clusters(filter2) %in% c("0","1")]))



mal_notract <- plot_cells(filter2, color_cells_by = "pseudotime", 
                       label_groups_by_cluster = FALSE,
                       label_branch_points = FALSE,
                       label_roots = FALSE,
                       label_leaves=FALSE)



mal_post | mal_notract



provua <- subset(tres, cell.types %in% c("T.CD8", "T.CD4", "T.cell") & treatment.group=="post.treatment")
provua2 <- subset(tres, cell.types %in% c("T.CD8", "T.CD4", "T.cell") & treatment.group=="treatment.naive")
malignes.treatment <- DimPlot(provua, reduction="umap", group.by="cell.types")
malignes.naive <- DimPlot(provua2, reduction="umap", group.by="cell.types")


malignes.treatment |malignes.naive





#Ja em miraré com utilitzar la llibreria dynverse. Però de moment, veiem que monocle3 funciona.



  
# HSC Interactome<img src="https://github.com/SzegezdiLab/HSC_Interactome/raw/main/hex.png" alt="hex sticker" width="130" align='right'/>



<br>This repository contains the code for a Cell-Cell Interactome shiny app for exploring cell-cell interactions between different cell types.
The app takes four uploads, all in the form of an .rds file (must be lowercase); 
1.	Interactome Data: An rds file containing at least the following columns of data:
*	‘source’: Interacting cell expressing a ligand
*	‘target’: Target cell the source cell is interacting with
* ‘ligand’: Ligand being expressed by the source cell
*	‘receptor’: Receptor expressed by target cell
*	‘timepoints’: Timepoints at which the above are expressed. Covered in this app are ‘Healthy’, ‘Diagnosis’, ‘Relapse’ and ‘Post-treatment’. 
*	‘aggregate_rank_(Timepoint)’ : an aggregate rank score generated using liana_aggregate function. Timepoint should be replaced with Healthy, Diagnosis etc. and is case sensitive
2.	A Ligand-Receptor Matrix: An rds file -- NB: Only required for generation of Violin plots, other plots will load without uploading this data
3.	Meta Data: An rds file containing UMI barcodes, their associated timepoint and their cell type. NB: Only required for generation of Violin plots, other plots will load without uploading this data
4.	A Colour Palette: An rds file that is a list of 2 lists; ‘timepoint’ and ‘celltype’. 
•	‘timepoint’ contains each timepoint e.g ‘Healthy’ and an associated colour. 
•	‘celltype’ contains each celltype e.g ‘Mesenchymal’ and an associated colour.

This app is a continuation of the HSC-Interactome shiny app made by Ennis S et. al as part of a study that can be found here: Cell-cell interactome of the hematopoietic niche and its changes in acute myeloid leukemia. *Ennis S et. al., iScience, 2023. DOI: [10.1016/j.isci.2023.106943](https://doi.org/10.1016/j.isci.2023.106943)*

To play with the original shiny app, [click here](https://sarahennis.shinyapps.io/HSC_Interactome/)

OR

You can download this repository and run the app locally by following the instructions below...


---

### Run docker container to host the app locally

The repository contains a Dockerfile and a compose.yaml file to build and run a RStudio server image to host and run the ShinyApp locally. To do so you need [Docker](https://www.docker.com/) and [Docker compose plugin](https://docs.docker.com/compose/) installed. 

Once the repository has been cloned, go inside the HSC_Interactome folder and run:

```bash
docker compose up
```
It will build and run the container.

With your browser go to: [http://localhost:8888/](http://localhost:8888/)

Login to RStudio:

- Name: rstudio
- Password: pass

Once the page has loaded, open one of the files (`server.R` or `ui.R`) and click on `Run App` on the top right corner to initialize the app.

If you have any questions or problems running the app, feel free to raise an issue [here](https://github.com/SzegezdiLab/HSC_Interactome/issues).

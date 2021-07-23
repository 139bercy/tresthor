## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(tresthor)
library(tidyverse)

## ----show_model , echo = FALSE,comment=NA-------------------------------------
cat(readLines(model_source_example(TRUE)), sep = '\n')

## ----create_model_1, eval = FALSE---------------------------------------------
#  ## Création rapide
#  create_model(model_name = "mon_modele",
#               model_source = "mon_dossier/modele.txt")
#  
#  ## Creation détaillée
#  create_model(model_name = "mon_modele",
#               model_source = "mon_dossier/modele.txt",
#               algo = FALSE, rcpp = TRUE , rcpp_path = "mon_dossier/sousdossier_rcpp")
#  

## ----create_model_2, eval = FALSE---------------------------------------------
#  ## Création manuelle
#  create_model(model_name = "mon_modele_manuel",
#               endogenous = c("x","y","z"),
#               exogenous = c("a","b","c"),
#               coefficients = c("coef_1", "coef_2"),
#               equations = c(
#                 "delta(1,log(x)) = 0.5* delta(1,log(a)) - 0.2 * delta(1,log(z))",
#                 "mon_equation1 :delta(1,log(y)) = coef_1 + coef_2 * delta(1,log(b)) - 0.3 * lag(x,1)",
#                 "calcul_de_z:z = lag(c,1) + x + y")  )
#  

## ----visualisation_eq1, include=FALSE,echo=FALSE------------------------------
## Création manuelle
create_model(model_name = "model_source_example",
             model_source = model_source_example(TRUE)  ) 


## ----visualisation_eq, echo=TRUE,eval=FALSE-----------------------------------
#  ## Création manuelle
#  create_model(model_name = "model_source_example",
#               model_source = model_source_example(TRUE)  )
#  

## ----visualisation_eq_3, echo=TRUE--------------------------------------------
model_source_example@equation_list %>% select(id,name,equation) 


## ----interrogation_mod--------------------------------------------------------
equation_info_model(c("equilibrium","eq_1"),model_source_example)

var_info_model(c("endovar4","cf3"),model_source_example)

## ----model_mod, eval = FALSE--------------------------------------------------
#  
#  model_endo_exo_switch(base_model = model_source_example,
#                        new_model_name = "mse_version2",
#                        new_endo = c("exovar2"),
#                        new_exo= c("endovar2"))
#  
#  model_equations_remove(base_model = model_source_example,
#                         new_model_name = "mse_version3",
#                         equations_to_remove = c("eq_2") ,
#                         endos_to_remove = c('endovar2')) ## Dans ce exemple endovar2 deviendra une variable exogène
#  
#  ##On crée une nouvelle équation
#  create_equation('new_equation',formula = "new_var = cf_4 *exovar2/endovar3",endogenous = "new_var",coefflist = c('cf_4'))
#  
#  ##On l'ajoute au modèle
#  model_equations_add(base_model = model_source_example,
#                      new_model_name = "mse_version4",
#                      thor_equations_add = list(new_equation))

## ----save_model, eval = FALSE-------------------------------------------------
#  save_model(mse_version3,folder_path = "mes_modeles")
#  
#  ##Deux manières de charger un modèle : nom + dossier
#  load_model(model = "mse_version3",folder_path = "mes_modeles")
#  load_model(file = "mes_models/mse_version3.rds")
#  
#  export_model(mse_version4,filename = "mes_modeles/nouveau_mod.txt")
#  

## ----format_date--------------------------------------------------------------
dates_annuelles <- seq.Date(as.Date("2010-01-01"), by="year", length.out = 5)
print(dates_annuelles)
cat(reformat_date(dates_annuelles))

dates_semestrielles <- seq.Date(as.Date("2010-01-01"), by="6 month", length.out = 5)
print(dates_semestrielles)
cat(reformat_date(dates_semestrielles))

dates_trimestrielles <- seq.Date(as.Date("2010-01-01"), by="quarter", length.out = 5)
print(dates_trimestrielles)
cat(reformat_date(dates_trimestrielles))

dates_mensuelles <- seq.Date(as.Date("2010-01-01"), by="month", length.out = 5)
print(dates_mensuelles)
# les mois s'afficheront dans la langue de l'installation de R de l'utilisateur
# puisque l'ordre alphabétique ne sera pas conservé, on ne peut pas utiliser ce format pour les données mensuelles comme indicateur temporel
cat(reformat_date(dates_mensuelles)) 


## ----add_coeff, eval = FALSE--------------------------------------------------
#  coeff_df <- data.frame(name  = c("coeff1","coeff2"),
#                         value = c(0.23,0.78) )
#  
#  donnees_completes <- add_coeffs(listcoeff = coeff_df, database = donnees_base ,
#                                  pos.coeff.name = 1, pos.coeff.value = 2, ##Specifier les colonnes où se trouvent le nom des coefficients et leur valeur.
#                                  overwrite = TRUE) #Pour remplacer les anciennes valeurs des coefficients
#  
#  

## ----thor_solver, eval = FALSE------------------------------------------------
#  donness_prevision <- thor_solver(model_source_example,
#                                  first_period = "2000Q1", last_period = "2020Q4",
#                                  database = mes_donnees, index_time = "dateQ" )

## ----system_expr , include=FALSE----------------------------------------------
sys_eq <-  noquote(paste(expression("a + log(b) = 2\\\\a + c = 4\\\\log(b) + c = 3")))
eq <-  noquote(paste(expression("log(x) + x = 4 + exo_1")))


## ----exemple_resolution_systeme_eval , results='hide'-------------------------
create_model(model_name = "exemple",
             equations = c("a + log(b) = 2",
                           "a + c = 4",
                           "log(b) + c = 3"),
                           endogenous = c("a","b","c") ,algo = FALSE)

donnees <- data.frame(obs = c(0,1), ## obs servira d'indicateur de temps
                        a = c(0,NA),
                        b = c(1,NA), ## b est utilisée en log donc 0 n'est pas une valeur initiale valide 
                        c = c(0,NA) )

solution <- thor_solver(exemple,first_period = 1,last_period = 1,
                          database = donnees , index_time = "obs") %>% 
            filter(obs == 1) %>% select(all_of(exemple@endo_list))



## ----print_solution , echo = FALSE--------------------------------------------
plop_sol <- paste( names(solution),"=", round(solution[1,],4),collapse = "\\\\" )
solution_exp <- noquote(paste(plop_sol))

## ----exemple_resolution_equation_eval , results='hide'------------------------
create_equation("eq_exemple",
                formula = "log(x) + x = 4 + exo_1",
                endogenous = "x")

donnees <- data.frame(obs = c(0,1), ## obs servira d'indicateur de temps
                        x = c(1,NA),
                    exo_1 = c(2,2)) ## bien renseigner la valeur de l'exogène sur les deux observations, celle de obs = 1 est utilisé pour la résolution

solution <- thor_equation_solver(eq_exemple,first_period = 1,last_period = 1,
                          database = donnees , index_time = "obs") %>% 
            filter(obs == 1) %>% select(paste0(eq_exemple@endogenous,".simul"))



## ----print_solution_equation , echo = FALSE-----------------------------------
plop_sol <- paste( eq_exemple@endogenous,"=", round(solution[1,1],4) )
solution_exp <- noquote(paste(plop_sol))

## ----quicksolve---------------------------------------------------------------
quick_solve("x^2 + log(x) - 3 = 0", "x", init = 1)

## ----create_equation----------------------------------------------------------
## création de l'équation manuelle
create_equation(equation_name = "eq_prix",
                formula = "y = b0 +a1*delta(1,x1) + a2*lag(x2,1) + residu",
                endogenous = "y",coefflist =   c("a1","a2","b0"))
## les coefficients sont réordonnés par ordre d'apparition dans la formule:
eq_prix@coefflist


## ----equation_from_model, eval = FALSE----------------------------------------
#  ## création de l'équation à partir du modèle
#  create_model("Opale",
#               model_source = system.file("Opale","opale.txt",package = "tresthor") )
#  
#  ## On récupère l'équation de consommation des ménages d'Opale
#  equation_from_model(equation_name_or_id = "conso_m",model = Opale,
#                      endogenous = "td_p3m_d7_ch",new_name = "conso_opale")

## ----eq_conso_ex--------------------------------------------------------------
create_equation("conso_ex",
                formula = "delta(1,log(conso)) = c0 + c1 * delta(1,log(rdb)) + c2 * delta(1,chom) +
                            m*(log(lag(conso,-1)) - lag(residu,-1)-lt0 - lt1*log(lag(rdb,1)))+
                            delta(1,residu)"  ,
                coefflist = c("c0","c1","c2","m","lt0","lt1"),endogenous = "conso")

## ----data_conso---------------------------------------------------------------
## téléchargement des données et constitution de la base
library(rdbnomics)
ids <- set_names(c("INSEE/CNT-2014-PIB-EQB-RF/T.CNT-EQUILIBRE_PIB.S14.P3.D-CNT.VALEUR_ABSOLUE.FE.L.EUROS.CVS-CJO",
                   "INSEE/CNT-2014-CSI/T.CNT-SOLDES_COMPTABLES_SECTEURS_INST.S14.SO.B6.VALEUR_ABSOLUE.FE.EUROS.CVS-CJO",
                   "INSEE/CHOMAGE-TRIM-NATIONAL/T.CTTXC.TAUX.FR-D976.0.00-.POURCENT.CVS"),
                 c("conso","rdb","chom"))

donnees_conso <- ids %>% imap(~{rdb(ids=.x) %>%  
                                select(period,value) %>%  
                                filter(period >= as.Date("1990-01-01") & period <= as.Date("2020-10-01")) %>% 
                                rename(date = period, !!.y := value) %>% as.data.frame()}) %>% 
                          reduce(full_join,by="date") %>% 
                          mutate(residu = 0, # on ajoute une colonne pour la variable résiduelle
                                 date = reformat_date(date)) #optionnel, pour transformer le format de date 


## ----formula_with_coeffs------------------------------------------------------
## création d'une base données avec des coefficients quelconques
data_coeffs <- data.frame(c0 = 0.5,c1=0.3,c2=-0.01, m = -0.4,lt0=0.01,lt1=0.98765432)

cat(formula_with_coeffs(conso_ex,database = data_coeffs,quiet = TRUE))

## on peut ajouter ces coeffficients à la base de données existante :
donnees_conso_c <- data_coeffs %>% t() %>% as.data.frame() %>% mutate(names = rownames(.)) %>%
                         add_coeffs(donnees_conso,pos.coeff.name = 2 , pos.coeff.value = 1)

# donnees_conso_c <- cbind(donnees_conso, data_coeffs) # fonctionne aussi car les coefficientsne sont pas dans la base et le format de data_coeffs s'y prête.

### et utiliser la base complète pour la fonction :
# cat(formula_with_coeffs(conso_ex,database = data_coeffs,quiet = TRUE))

## ----formula_latex_ex---------------------------------------------------------
cat(formula_latex(conso_ex))

cat(formula_latex(formula_with_coeffs(conso_ex,database = donnees_conso_c,quiet = TRUE)))

## ---- echo=FALSE--------------------------------------------------------------
plop<-formula_latex(formula_with_coeffs(conso_ex,database = donnees_conso_c,quiet = TRUE))
plop2<-formula_latex(conso_ex)

## -----------------------------------------------------------------------------
conso_solved <- thor_equation_solver(conso_ex , first_period = "2019Q1", last_period = "2020Q4",
                                      database = donnees_conso_c,index_time = "date")
tail(conso_solved,10)

## ----quick_estim--------------------------------------------------------------

conso_data <- quick_estim(thor_equation = conso_ex , database = donnees_conso,
                          estim_start = "1991Q1", estim_end = "2005Q4",
                          coeff_lt = c("lt0","lt1"))

nouvelle_formule<-formula_latex(formula_with_coeffs(conso_ex,database=conso_data,round_digits = 4,quiet = TRUE))



## ----quick_estim_all, eval = FALSE--------------------------------------------
#  ### On suppose que le modèle "Mod_Prev" contient les équations nommées suivantes : "invest" , "exports", "imports"
#  
#  ## on crée la liste suivante
#  
#  list_equation <- list(invest = list(endogenous_name = "inv_s",
#                                      estim_start = "2000Q1", estim_end = "2018Q4" ,
#                                      coef_lt = NULL , const = TRUE) ,
#                        exports = list(endogenous_name = "exp_tot" ,
#                                      estim_start =  "2001Q1", estim_end = "2017Q4",
#                                      coef_lt = c("b_lt0","b_lt1"), const = FALSE) ,
#                        imports = list(endogenous_name = "imp_tot" ,
#                                      estim_start =  "2000Q1", estim_end = "2018Q2",
#                                      coef_lt = c("a_lt1","a_lt2","a_lt3"), const = TRUE))
#  
#  data_quick_estim <- quick_estim_all(info_equations = list_equation , thoR.model = Mod_Prev,
#                                      database = mes_donnees,index_time = "date")

## ----sim_obs------------------------------------------------------------------

simul_data <- simulate_equation(conso_ex, database = conso_data, 
                                start_sim = "1991Q2", end_sim = "2020Q4",index_time = "date", 
                                residual_var = "residu")

tail(simul_data,10)

## ----graph_simobs-------------------------------------------------------------
graph_sim_obs(simul_data , start_plot = "2010Q1",type = "lvl")

## ----graph_sim_obs_pipe-------------------------------------------------------

conso_ex %>% simulate_equation( database = conso_data, 
                                start_sim = "2010Q1", end_sim = "2020Q4",index_time = "date", 
                                residual_var = "residu") %>% 
              graph_sim_obs(start_plot = "2009Q1",end_plot = "2019Q4",
                            title = "Equation de consommation des ménages",type = 'g',
                            colours = c("hotpink","olivedrab3"))


## ----graph sim_obs2-----------------------------------------------------------
conso_ex %>% simulate_equation( database = conso_data, 
                                start_sim = "1991Q1", end_sim = "2020Q4",index_time = "date", 
                                residual_var = "residu") %>% graph_sim_obs(start_plot = "2000Q1", type = "lvl")

## ----graph sim_obs3-----------------------------------------------------------
conso_ex %>% simulate_equation( database = conso_data, 
                                start_sim = "1991Q1", end_sim = "2019Q4",index_time = "date", 
                                residual_var = "residu") %>% graph_sim_obs(start_plot = "2000Q1", type = "g")

## ----contribs-----------------------------------------------------------------
essai1 <- dyn_contribs(conso_ex,conso_data,"2000Q1","2020Q4",residual_var = "residu")

essai2 <- dyn_contribs(conso_ex,conso_data,"1992Q1","2020Q4",residual_var = "residu")

## ----contribs data------------------------------------------------------------
tail(essai1,5) 

## ----graph_contrib------------------------------------------------------------
graph_contrib(essai1,start_plot = "2018Q1",end_plot = "2019Q4")

## ----graph_contrib2-----------------------------------------------------------
essai2 %>% graph_contrib(start_plot = "2017Q1",end_plot = "2019Q4",
                         title = "Contributions dynamiques trimestrielles :\n consommation des ménages",
                         name_endogenous = "Consommation des ménages",
                         colour_line = "springgreen4",
                         colours_bar = c("hotpink","hotpink4","lightpink","lightgoldenrod"), ##on note que lightpink ne sera pas utilisé
                         bar_labels = c("Revenu disponible brut des ménages","Taux de chômage","Inexpliqué"),
                         legend_n_row = 3)

## ----index year , eval=TRUE---------------------------------------------------
### exemples de création de vecteur pour l'argument index_year

## Si le colonne de dates trimestrielles est au format date (ie as.Date("2010-07-01")) :

vecteur_dates_trim <- seq.Date(from = as.Date("2010-07-01"),by = "quarter",length.out = 8)
index_annees <- as.numeric(format(vecteur_dates_trim,"%Y"))

vecteur_dates_trim
index_annees

## Si le colonne de dates trimestrielles est au format simplifié (ie "2010Q3"):
vecteur_dates_trim2 <- reformat_date(vecteur_dates_trim)
index_annees2 <- as.numeric(gsub("(.*)Q[1-4]$","\\1",vecteur_dates_trim2))

vecteur_dates_trim2
index_annees2


## ----contrib_ann--------------------------------------------------------------
contrib_ann <- yearly_contrib(essai2 , index_year = as.numeric(gsub("(.*)Q[1-4]$","\\1",essai1$date)))

tail(contrib_ann,5)

## ----graph contribann---------------------------------------------------------
graph_contrib(contrib_ann, start_plot = 2000,end_plot = 2020)

## ----graph contribann2--------------------------------------------------------
conso_ex %>% 
  dyn_contribs(conso_data,"1992Q1","2020Q4",residual_var = "residu") %>% 
    yearly_contrib(index_year = as.numeric(gsub("(.*)Q[1-4]$","\\1",essai1$date))) %>% 
      graph_contrib(start_plot = 2010,end_plot = 2020,
                    title = "Contributions dynamiques annuelles:\nconsommation des ménages",
                         name_endogenous = "Consommation des ménages",
                         colour_line = "red3",
                         colours_bar = c("darkblue","steelblue","lightskyblue","seashell3"), 
                         bar_labels = c("Revenu disponible brut des ménages","Taux de chômage","Inexpliqué"),
                         legend_n_row = 3)
  

## ----quickplots---------------------------------------------------------------

quick_plot(c("conso","rdb"),conso_data )
quick_plot(c("conso","rdb"),conso_data,
           start = "2005Q1", growth_rate = TRUE, 
           colours = c("hotpink","springgreen3","gold3"))

quick_plot(c("dlog_obs","g_obs"), simul_data,start = "2019Q1" , 
           title = "Comparaison des méthodes de calul du taux de croissance", 
           colours = c('cornflowerblue','springgreen3') )



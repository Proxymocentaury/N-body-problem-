include("nbodylib.jl")  
#include("nbodylibMT.jl") # liberia che sfrutta la parallelizzazione (avviare con:  julia -t n main.jl ) 


t_parm= 1e-3

N_dt= 10e12 # numero passi massimi della simulazione
tmax= 2*pi  # tempo massimo della simulazione 

snapshot_freq=500 #N_dt #indica ogni quanti passi vengono archiviati i risultati
check_freq= 1e3 #indica ogni quanti passi viene stampato il tempo e l'energia

errore_soglia= 1e-3


# --------- estrazioni dei dati ------------

mass, coord, vel= prelevo_dati_csv() 

stato_finale=simulazione_nostop(t_parm, tmax, mass, coord, vel, N_dt, check_freq, snapshot_freq)
#println(stato_finale)


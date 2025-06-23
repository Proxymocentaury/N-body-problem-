include("nbodylib.jl")  # liberia che calcola accelerazioni e jerk insieme (velocizza l'esecuzione)
#include("nbodylibMT.jl") # liberia che sfrutta la parallelizzazione (avviare con:  julia -t 2 main.jl ) 
include("anima.jl") 

#-------------------------simulazioni-------------------------------

#simulazione_nostop(t_parm, tmax, mass, coord, vel,  N_dt, check_freq, snapshot_freq)
#simulazione_break(t_parm, tmax, mass, coord, vel, N_dt, errore_soglia)


#- - - - - - - - - - - - - definizione dei parametri per la simulazione- - - - - - - - - - - - - - - - -


t_parm= 1e-3

N_dt= 10e10 # numero passi massimi della simulazione
tmax= 2*pi  # tempo massimo della simulazione 

snapshot_freq=50 #N_dt #indica ogni quanti passi vengono archiviati i risultati
check_freq= 100000 #indica ogni quanti passi viene stampato il tempo e l'energia

errore_soglia= 1e-3


# --------- estrazioni dei dati ------------

mass, coord, vel= prelevo_dati_csv() 
#mass, coord, vel= prelevo_dati_csv_MT( "data.csv" ) # estrae dati in file csv, vale solo quando si usa nbodylibMT
#mass, coord, vel= prelevo_dati_jld2( "corpi.jld2" )  # estrae dati in file jld2

#-----------------------------------------

stato_finale=simulazione_nostop(t_parm, tmax, mass, coord, vel, N_dt, check_freq, snapshot_freq)
#println(stato_finale)

crea_gif_2d("data_output.jld2", "output.gif", 30, 1.1)

#@time A = simulazione_nostop(t_parm, tmax, mass, coord, vel, N_dt, check_freq,snapshot_freq)

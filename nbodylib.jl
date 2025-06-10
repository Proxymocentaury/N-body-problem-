using CSV
using DataFrames
using DelimitedFiles
using LinearAlgebra

###############################
#Forse sarebbe meglio mettere G=1 ?
###############################

G = 6.67430e-11


################################################
# funzioni per la lettura e scrittura dei dati 
################################################

# funzione per prelevare i dati da un file csv
function prelevo_dati_csv()

    # Leggi il file CSV senza intestazioni
    df = CSV.read("data.csv", DataFrame, header=false)

    data_matrix = Matrix(df)  # Converte il DataFrame in una matrice

    mass = data_matrix[:, 1]   
    coord = data_matrix[:, 2:4] 
    vel = data_matrix[:, 5:7] 

    return mass, coord, vel
    
end


# Funzione per scrivere la matrice su file CSV
function snapshot_csv(io, t, A)
    println(io, "tempo: $t")
    
    # Scrivi ogni riga della matrice con elementi separati da virgola e spazio
    for i in axes(A, 1)
        println(io, join(A[i, :], ", "))
    end
    
    println(io) # Riga vuota dopo ogni blocco di dati
    
    flush(io) # Assicurati che i dati vengano effettivamente scritti sul disco
end

# Funzione per prelevare i dati da un file JLD2
function prelevo_dati_jld2(jld2_file::String)
    # Carica i dati dal file JLD2
    @load jld2_file A

    mass = A[:, 1]
    coord = A[:, 2:4]
    vel = A[:, 5:7]

    return mass, coord, vel
end




############################################################
# funzioni per calcolare le grandezze fisiche 
############################################################


#funzione per calcolare l'energia cinetica di singola particella
function calcolo_Ek(mass,vel)
    N=length(mass)
    Ek=0
    for i in 1:N
        Ek += 0.5*mass[i]*sum(vel[i, :] .^2) 
    end
    return Ek
    
end



#funzione per calcolare l'energia potenziale del sistema
function calcolo_Ep(mass, coord)
    N=length(mass)
    p=0

    for i in 1:N-1
        for j in i+1:N
            p+= -G*mass[i]*mass[j]/( sqrt(sum((coord[i, :] - coord[j, :]) .^ 2)) )
        end
        
    end
    return p
end



#funzione per calcolare il tempo di collisione
function calcolo_tcoll(coord, vel, acc)
    N = size(coord,1)

    coll_time =  Inf  #inf  

    for i in 1:N-1
        for j in i+1:N
            r_ij = coord[j, :] - coord[i, :] # vettore ri- rj
            v_ij = vel[j, :] - vel[i, :] # vettore vi- vj
            a_ij = acc[j, :] - acc[i, :] # vettore ai- aj
            
            r2 = sum(r_ij .^ 2)  # Distanza al quadrato
            v2 = sum(v_ij .^ 2)  # Velocità relativa al quadrato
            a2 = sum(a_ij .^ 2)  # Accelerazione relativa al quadrato
            
            if v2 != 0 
                tcoll1 = sqrt( r2 / v2 )
                coll_time = min(coll_time, tcoll1)
            end
            
            if a2 != 0 
                tcoll2 = (r2 / a2)^0.25
                coll_time = min(coll_time, tcoll2)
            end
        end
    end
     
    return coll_time  
end



# calcola l'erore relativo dell'energia
function diagnostica(Ein, E)
    return abs( (Ein - E)/Ein )
end


#################
# funzioni il calcolo delle matrici per l'integrazione
################

# funzione per calcolare la matrice accellerazione
function calcolo_accelerazioni(mass, coord)
    N = length(mass)
    acc = zeros(N, 3)

    for i in 1:N
        for j in 1:N
            if i != j
                r_ij = coord[j, :] - coord[i, :]  # Vettore posizione relativa
                r_norm2 = sum(r_ij .^ 2)  # Distanza al quadrato
                if r_norm2 > 1e-10  
                    acc[i, :] += G * mass[j] * r_ij / r_norm2^1.5
                end
            end
        end
    end
    return acc
end



# funzione per calcolare la matrice jerk 
function calcolo_jerk(mass, coord, vel)
    N = length(mass) 
    jerk = zeros(N, 3) 

    for i in 1:N
        for k in 1:N
            if k != i
                r_ij = coord[k, :] - coord[i, :]  
                v_ij = vel[k, :] - vel[i, :]      
                r2 = sum(r_ij .^ 2)               
                r = sqrt(r2)                      

                if r2 > 0  
                    v_rad = dot(r_ij, v_ij)  

                    jerk[i, :] .+= G * mass[k] * (v_ij / r^3 - 3 * v_rad * r_ij / r^5)
                end
            end
        end
    end

    return jerk
end




###########################
# funzioni integrazione 
###########################


function evoluzione_coordinate(dt, coord, vel, acc, jerk)
    new_coord = copy(coord)  
    N = size(new_coord, 1)   
    for i in 1:N
        new_coord[i, :] .= coord[i, :] .+ vel[i, :] * dt .+ 0.5 * acc[i, :] * dt^2 .+ jerk[i, :] * dt^3 / 6
    end
    return new_coord
end



function evoluzione_velocità(dt, vel, acc, jerk)
    new_vel = copy(vel)  
    N = size(new_vel, 1)  
    for i in 1:N
        new_vel[i, :] .= vel[i, :] .+ acc[i, :] * dt .+ 0.5 * jerk[i, :] * dt^2
    end
    return new_vel
end




# restituise la matrice delle coordinate e delle velocià aggiornate e correte
function evolve_step(dt,mass, coord_old, vel_old, acc_old, jerk_old )

    coord_new=evoluzione_coordinate(dt, coord_old, vel_old, acc_old, jerk_old)
    vel_new=evoluzione_velocità(dt, vel_old, acc_old, jerk_old)

    acc_new=calcolo_accelerazioni(mass, coord_new)
    jerk_new=calcolo_jerk(mass, coord_new, vel_new)

    vel_new = vel_old + 1/2*(acc_new + acc_old)*dt +  1/12*(-jerk_new + jerk_old)*dt^2
    coord_new = coord_old + (vel_old) * dt + 1/4 * (acc_new + acc_old) * dt^2  - 1/12 * (jerk_new) * dt^3 

    return coord_new, vel_new, acc_new, jerk_new
end


#---- simulazioni, salvano i risultati su un file JLD2

using JLD2

#simulazione: porta a termine i calcoli fino al tempo richiesto
function simulazione_nostop(t_parm, tmax, mass, coord, vel,  N_dt)

    acc=calcolo_accelerazioni(mass,coord)
    jerk=calcolo_jerk(mass, coord, vel)

    Ein= calcolo_Ek(mass, vel) + calcolo_Ep(mass, coord) # energia totale del sistema allo stato iniziale
 

    t = 0  
    ti = 0  

    # Inizializza le variabili per evitare UndefVarError
    coord_new, vel_new, acc_new, jerk_new = coord, vel, acc, jerk

    A = hcat(mass, coord, vel)

    jldopen("data_output.jld2", "w") do io
        io["snapshot_0"] = (t=t, A=A)

        while ti < N_dt && t < tmax
            dt = calcolo_tcoll(coord, vel, acc) * t_parm
            coord_new, vel_new, acc_new, jerk_new = evolve_step(dt, mass, coord, vel, acc, jerk)

            coord, vel, acc, jerk = coord_new, vel_new, acc_new, jerk_new
            ti += 1
            t += dt

            if ti % check_freq == 0
                E = calcolo_Ek(mass, vel_new) + calcolo_Ep(mass, coord_new)
                E_err=abs( (E-Ein)/Ein )
                println("tempo: $t, energia: $E, (energia stato iniziale) $Ein,  errore energia: $E_err")
            end

            if ti % snapshot_freq == 0
                A = hcat(mass, coord, vel)
                io["snapshot_$ti"] = (t=t, A=A)
                E = calcolo_Ek(mass, vel_new) + calcolo_Ep(mass, coord_new)
                E_err=abs( (E-Ein)/Ein )
                println("tempo: $t, energia: $E, (energia stato iniziale) $Ein,  errore energia: $E_err")
            end
        end

        println("\nrisultati finali")   
        E = calcolo_Ek(mass, vel_new) + calcolo_Ep(mass, coord_new)
        println("tempo: $t, passo: $ti")
        err = diagnostica(Ein, E)
        println("errore relativo dell'energia: $err")

        A = hcat(mass, coord, vel)
        io["final_snapshot"] = (t=t, A=A)
    end 

    return A
end


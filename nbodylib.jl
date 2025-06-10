using CSV, JLD2, DataFrames
using DelimitedFiles
using LinearAlgebra

###############################
#Forse sarebbe meglio mettere G=1 ?
###############################

G = 1 #6.67430e-11


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


# Funzione per scrivere i dati nel file CSV con intestazioni in una riga
function scrivi_dati_csv(io, t, mass, coord, vel)
    # Scrivi l'intestazione solo una volta
    if t == 0
        header = ["tempo"]
        for i in 1:length(mass)
            append!(header, ["massa$i", "coord_x$i", "coord_y$i", "vel_x$i", "vel_y$i"])
        end
        println(io, join(header, ","))
    end

    # Scrivi i dati
    data = Float64[t]
    N = length(mass)
    for i in 1:N
        append!(data, [mass[i], coord[i, 1], coord[i, 2], vel[i, 1], vel[i, 2]])
    end
    println(io, join(data, ","))
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


# funzione combinata per calcolare sia le accelerazioni che i jerk
function calcolo_acc_jerk(mass, coord, vel)
    N = length(mass)
    acc = zeros(N, 3)
    jerk = zeros(N, 3)

    # Calcola le accelerazioni e i jerk in un'unica iterazione
    for i in 1:N
        for k in 1:N
            if k != i
                r_ij = coord[k, :] - coord[i, :]  # Vettore posizione relativa
                v_ij = vel[k, :] - vel[i, :]      # Vettore velocità relativa
                r2 = sum(r_ij .^ 2)               # Distanza al quadrato
                r = sqrt(r2)                      # Distanza

                if r2 > 0
                    # Calcola l'accelerazione per la particella i dovuta alla particella k
                    acc[i, :] += G * mass[k] * r_ij / r2^1.5

                    # Calcola il termine di jerk per la particella i dovuto alla particella k
                    v_rad = dot(r_ij, v_ij)  # Velocità radiale relativa
                    jerk[i, :] += G * mass[k] * (v_ij / r^3 - 3 * v_rad * r_ij / r^5)
                end
            end
        end
    end

    return acc, jerk
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

    acc_new, jerk_new = calcolo_acc_jerk(mass, coord_new, vel_new)
    

    vel_new = vel_old + 1/2*(acc_new + acc_old)*dt +  1/12*(-jerk_new + jerk_old)*dt^2
    coord_new = coord_old + (vel_old) * dt + 1/4 * (acc_new + acc_old) * dt^2  - 1/12 * (jerk_new) * dt^3 

    return coord_new, vel_new, acc_new, jerk_new
end


######################################################
#---- simulazioni, salvano i risultati su un file JLD2
######################################################



#simulazione: porta a termine i calcoli fino al tempo richiesto
function simulazione_nostop(t_parm, tmax, mass, coord, vel,  N_dt, check_freq, snapshot_freq)

    acc, jerk = calcolo_acc_jerk(mass, coord, vel)

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


#simulazione: porta a termine i calcoli fino al tempo richiesto
function simulazione_break(t_parm, tmax, mass, coord, vel,  N_dt, errore_soglia, check_freq, snapshot_freq)

    acc, jerk = calcolo_acc_jerk(mass, coord, vel)

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

            E=calcolo_Ek(mass, vel) + calcolo_Ep(mass, coord)
            E_err= diagnostica(Ein, E)

            if E_err >= errore_soglia  
                println("Errore sopra soglia! Terminazione anticipata al passo $ti.")
                println("tempo: $t, energia: $E, (energia stato iniziale) $Ein,  errore energia: $E_err")
                A=hcat(mass, coord, vel) 
                io["snapshot_$ti"] = (t=t, A=A)
                return  A 
            end



            if ti % check_freq == 0
                E = calcolo_Ek(mass, vel_new) + calcolo_Ep(mass, coord_new)
                E_err=abs( (E-Ein)/Ein )
                println("tempo: $t, energia: $E, (energia stato iniziale) $Ein,  errore energia: $E_err")
            end

            if ti % snapshot_freq == 0
                A = hcat(mass, coord, vel)
                io["snapshot_$ti"] = (t=t, A=A)
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


# simulazione creata in modo tale che l'errore dell'energia non superi il valore soglia
function simulazione_autocorrect(t_parm, tmax, mass, coord, vel, N_dt, errore_soglia, snapshot_freq)
    acc, jerk = calcolo_acc_jerk(mass, coord, vel)
    Ein = calcolo_Etot(mass, coord, vel)

    t = 0
    ti = 0
    E_start = Ein
    E = Ein

    coord_new, vel_new, acc_new, jerk_new = coord, vel, acc, jerk
    A = hcat(mass, coord, vel)

    jldopen("data_output.jld2", "w") do io
        io["snapshot_0"] = (t=t, A=A)

        while ti < N_dt && t < tmax
            dt = calcolo_tcoll(coord, vel, acc) * t_parm

            # questo ciclo, dimezza il passo temporale finché l'errore dell'energia è inferiore all'errore soglia
            while diagnostica(E_start, E) > errore_soglia
                dt *= 0.5
                coord_new, vel_new, acc_new, jerk_new = evolve_step(dt, mass, coord, vel, acc, jerk)
                E = calcolo_Etot(mass, coord_new, vel_new)
            end

            coord_new, vel_new, acc_new, jerk_new = evolve_step(dt, mass, coord, vel, acc, jerk)
            E = calcolo_Etot(mass, coord_new, vel_new)

            coord, vel, acc, jerk = coord_new, vel_new, acc_new, jerk_new
            ti += 1
            t += dt

            if ti % snapshot_freq == 0
                A = hcat(mass, coord, vel)
                io["snapshot_$ti"] = (t=t, A=A)
                println("tempo: $t, energia: $E, errore energia: $(diagnostica(Ein, E))")
            end
        end

        println("\nrisultati finali")
        println("passo $ti: Errore energia: $(diagnostica(Ein, E)), dt = $dt")
        println("tempo: $t")
        A = hcat(mass, coord, vel)
        io["final_snapshot"] = (t=t, A=A)
    end

    return hcat(mass, coord, vel)
end


function simulazione_con_ricalcolo(t_parm, tmax, mass, coord, vel, N_dt, soglia_err, snapshot_freq)
    acc, jerk = calcolo_acc_jerk(mass, coord, vel)
    Ein = calcolo_Etot(mass, coord, vel)

    t = 0
    ti = 0

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

            E = calcolo_Etot(mass, coord_new, vel_new)
            err = diagnostica(Ein, E)

            # Ricalcola finché l'errore non è accettabile
            while err > soglia_err
                dt /= 2  # Dimezza il tempo dello step
                coord_new, vel_new, acc_new, jerk_new = evolve_step(dt, mass, coord, vel, acc, jerk)
                coord, vel, acc, jerk = coord_new, vel_new, acc_new, jerk_new
                t += dt  # Aggiungi il nuovo dt al tempo totale

                # Ricalcola energia e errore con il nuovo passo
                E = calcolo_Etot(mass, coord_new, vel_new)
                err = diagnostica(Ein, E)
            end

            # Altrimenti, salva lo stato e continua
            if ti % snapshot_freq == 0
                A = hcat(mass, coord, vel)
                io["snapshot_$ti"] = (t=t, A=A)
                println("tempo: $t, energia: $E, errore energia: $err")
            end
        end

        println("\nrisultati finali")
        E = calcolo_Etot(mass, coord_new, vel_new)
        println("tempo: $t, passo: $ti")
        err = diagnostica(Ein, E)
        println("errore relativo dell'energia: $err")

        A = hcat(mass, coord, vel)
        io["final_snapshot"] = (t=t, A=A)
    end

    return A
end


function simulazione_con_ricalcolo_adattivo(t_parm, tmax, mass, coord, vel, N_dt, soglia_err, snapshot_freq, rid=2.0, dt_min=1e-10)
    acc, jerk = calcolo_acc_jerk(mass, coord, vel)
    Ein = calcolo_Etot(mass, coord, vel)

    t = 0
    ti = 0

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

            E = calcolo_Etot(mass, coord_new, vel_new)
            err = diagnostica(Ein, E)

            # Se l'errore è troppo grande, riduci il passo (con un fattore di riduzione)
            while err > soglia_err
                println("\nErrore superiore alla soglia ($err > $soglia_err), ricalcolando con step ridotto.")
                dt *= 1 / rid  # Riduci il passo
                dt = max(dt, dt_min)

                coord_new, vel_new, acc_new, jerk_new = evolve_step(dt, mass, coord, vel, acc, jerk)
                coord, vel, acc, jerk = coord_new, vel_new, acc_new, jerk_new
                t += dt  # Aggiungi il nuovo dt al tempo totale

                # Ricalcola energia e errore con il nuovo passo
                E = calcolo_Etot(mass, coord_new, vel_new)
                err = diagnostica(Ein, E)
            end

            # Salva lo stato e continua
            if ti % snapshot_freq == 0
                A = hcat(mass, coord, vel)
                io["snapshot_$ti"] = (t=t, A=A)
                println("tempo: $t, energia: $E, errore energia: $err")
            end
        end

        println("\nrisultati finali")
        E = calcolo_Etot(mass, coord_new, vel_new)
        println("tempo: $t, passo: $ti")
        err = diagnostica(Ein, E)
        println("errore relativo dell'energia: $err")

        A = hcat(mass, coord, vel)
        io["final_snapshot"] = (t=t, A=A)
    end

    return A
end



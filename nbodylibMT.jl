using CSV
using DataFrames
using DelimitedFiles
using LinearAlgebra
using Base.Threads
using JLD2

G = 1 #6.67430e-11

#######################
#GESTIONE LETTURA E SCRITTURA DATI
#######################

function prelevo_dati_csv_MT(file::String)
    df = CSV.read(file, DataFrame, header=false)
    data_matrix = Matrix(df)

    mass = Vector{Float64}(undef, size(data_matrix, 1))
    coord = Matrix{Float64}(undef, size(data_matrix, 1), 3)
    vel = Matrix{Float64}(undef, size(data_matrix, 1), 3)

    @threads for i in 1:size(data_matrix, 1)
        mass[i] = data_matrix[i, 1]
        coord[i, :] = data_matrix[i, 2:4]
        vel[i, :] = data_matrix[i, 5:7]
    end

    return mass, coord, vel
end



# Funzione per prelevare i dati da un file JLD2 ha come argomento il nome del file
function prelevo_dati_jld2(jld2_file::String)
    # Carica i dati dal file JLD2
    @load jld2_file A

    mass = A[:, 1]
    coord = A[:, 2:4]
    vel = A[:, 5:7]

    return mass, coord, vel
end



##########################
#CALCOLO ENERGIA ED ERRORE
##########################


function calcolo_Ek(mass, vel)
    N = length(mass)
    # Ek = Threads.Atomic{Float64}(0.0)
    # Ek = 0.0
    Ekt = zeros(Threads.nthreads())
    @threads for i in 1:N
        # Threads.atomic_add!(Ek, 0.5 * mass[i] * sum(vel[i, :] .^ 2))
        Ekt[Threads.threadid()] += 0.5 * mass[i] * sum(vel[i, :] .^ 2)
    end

    return sum(Ekt)
end




function calcolo_Ep(mass, coord)
    N = length(mass)
    # Ep = Threads.Atomic{Float64}(0.0)
    # Ep = 0.0
    Ept = zeros(Threads.nthreads())

    #= @threads=# for i in 1:N-1
        for j in i+1:N
            r = sqrt(sum( (coord[i, :] - coord[j, :]) .^ 2))
            # Threads.atomic_add!(Ep, -G * mass[i] * mass[j] / r)
            Ept[Threads.threadid()] += - mass[i] * mass[j] / r
        end
    end

    return sum(Ept)
end




function calcolo_Etot(mass, coord, vel) 
    # Ek, Ep = @sync begin
    #     t1 = Threads.@spawn calcolo_Ek(mass, vel)
    #     t2 = Threads.@spawn calcolo_Ep(mass, coord)
    #     Ek = fetch(t1)
    #     Ep = fetch(t2)
    #     Ek, Ep
    # end
    Ek =  calcolo_Ek(mass, vel)
    Ep =  calcolo_Ep(mass, coord)        
    return Ek + Ep
end




function diagnostica(Ein, E)
    return abs((Ein - E) / Ein)
end


#################
# FUNZIONI PER IL CALCOLO DELLE MATRICI PER L'INTEGRAZIONE, E PER IL TEMPO DI COLLISIONE
################



function calcola_vij(v, i, j)
    return v[j, :] - v[i, :]
end


# funzione per calcolare il tempo di collisione
function calcolo_tcoll(coord, vel, acc)
    N = size(coord, 1)
    coll_time = Threads.Atomic{Float64}(Inf)
    @sync begin
        @spawn for i in 1:N-1
            @inbounds for j in i+1:N
                r_ij = coord[j, :] - coord[i, :]
                v_ij = vel[j, :] - vel[i, :]
                a_ij = acc[j, :] - acc[i, :]

                r2 = sum(r_ij .^ 2)
                v2 = sum(v_ij .^ 2)
                a2 = sum(a_ij .^ 2)

                if v2 != 0
                    tcoll1 = sqrt(r2 / v2)
                    Threads.atomic_min!(coll_time, tcoll1)
                end

                if a2 != 0
                    tcoll2 = (r2 / a2)^0.25
                    Threads.atomic_min!(coll_time, tcoll2)
                end
            end
        end
    end
    return coll_time[]
end


# funzione per calcolare accelletrazioni e jerk 
function calcolo_acc_jerk(mass, coord, vel)
    N = length(mass)
    acc = zeros(N, 3)
    jerk = zeros(N, 3)

    # @sync begin
        @threads for i in 1:N
            @inbounds for k in 1:i-1 
                    r_ij = [coord[k, 1] - coord[i, 1];coord[k, 2] - coord[i, 2]; coord[k, 3] - coord[i, 3]] 
                    r = sqrt(sum(r_ij .^ 2))

                    # if r > 0
                        v_ij = [vel[k, 1] - vel[i, 1];vel[k, 2] - vel[i, 2]; vel[k, 3] - vel[i, 3]] 
                        acc[i, 1] +=  mass[k] * r_ij[1] / r^3
                        acc[i, 2] += mass[k] * r_ij[2] / r^3
                        acc[i, 3] +=  mass[k] * r_ij[3] / r^3
                        v_rad = dot(r_ij, v_ij)
                        jerk[i, 1] +=  mass[k] * (v_ij[1] / r^3 - 3 * v_rad * r_ij[1] / r^5)
                        jerk[i, 2] +=  mass[k] * (v_ij[2] / r^3 - 3 * v_rad * r_ij[2] / r^5)
                        jerk[i, 3] +=  mass[k] * (v_ij[3] / r^3 - 3 * v_rad * r_ij[3] / r^5)
                    # end
            end

            @inbounds for k in i+1:N
                r_ij = [coord[k, 1] - coord[i, 1];coord[k, 2] - coord[i, 2]; coord[k, 3] - coord[i, 3]] 
                r = sqrt(sum(r_ij .^ 2))

                # if r > 0
                    v_ij = [vel[k, 1] - vel[i, 1];vel[k, 2] - vel[i, 2]; vel[k, 3] - vel[i, 3]] 
                    acc[i, 1] +=  mass[k] * r_ij[1] / r^3
                    acc[i, 2] +=  mass[k] * r_ij[2] / r^3
                    acc[i, 3] +=  mass[k] * r_ij[3] / r^3
                    v_rad = dot(r_ij, v_ij)
                    jerk[i, 1] +=  mass[k] * (v_ij[1] / r^3 - 3 * v_rad * r_ij[1] / r^5)
                    jerk[i, 2] +=  mass[k] * (v_ij[2] / r^3 - 3 * v_rad * r_ij[2] / r^5)
                    jerk[i, 3] +=  mass[k] * (v_ij[3] / r^3 - 3 * v_rad * r_ij[3] / r^5)
                # end
            end
        end
    # end

    return acc, jerk
end


# funzione per calcolare accelletrazioni e jerk 
function calcolo_acc_jerk2(mass, coord, vel)
    N = length(mass)
    acc = zeros(N, 3)
    jerk = zeros(N, 3)

    @sync begin
        @threads for i in 1:N
            for k in 1:N
                if k != i
                    r_ij = coord[k, :] - coord[i, :]
                    v_ij = vel[k, :] - vel[i, :]
                    r2 = sum(r_ij .^ 2)
                    r = sqrt(r2)

                    if r2 > 0
                        acc[i, :] += G * mass[k] * r_ij / r2^1.5
                        v_rad = dot(r_ij, v_ij)
                        jerk[i, :] += G * mass[k] * (v_ij / r^3 - 3 * v_rad * r_ij / r^5)
                    end
                end
            end
        end
    end

    return acc, jerk
end




function evoluzione_coordinate_velocità(dt, coord, vel, acc, jerk)
    
    N = size(coord, 1)
    new_coord = Matrix{Float64}(undef, N, 3)
    new_vel = Matrix{Float64}(undef, N, 3)
    @threads for i in 1:N
        new_coord[i, :] .= coord[i, :] .+ vel[i, :] * dt .+ 0.5 * acc[i, :] * dt^2 .+ jerk[i, :] * dt^3 / 6
        new_vel[i, :] .= vel[i, :] .+ acc[i, :] * dt .+ 0.5 * jerk[i, :] * dt^2
    end
    return new_coord, new_vel
    
end



function evolve_step(dt, mass, coord_old, vel_old, acc_old, jerk_old)
    coord_new, vel_new = evoluzione_coordinate_velocità(dt, coord_old, vel_old, acc_old, jerk_old)


    acc_new, jerk_new = calcolo_acc_jerk(mass, coord_new, vel_new)

    @threads for i in 1:length(mass)
        vel_new[i, :] .= vel_old[i, :] .+ 1/2 * (acc_new[i, :] + acc_old[i, :]) * dt .+ 1/12 * (-jerk_new[i, :] + jerk_old[i, :]) * dt^2
        coord_new[i, :] .= coord_old[i, :] .+ vel_old[i, :] * dt .+ 1/4 * (acc_new[i, :] + acc_old[i, :]) * dt^2 .- 1/12 * jerk_new[i, :] * dt^3
    end

    return coord_new, vel_new, acc_new, jerk_new
end



################
#SIMULAZIONI
################



# simululazione che non si interrompe fino al tempo richiesto
function simulazione_nostop(t_parm, tmax, mass, coord, vel, N_dt, check_freq, snapshot_freq)
    
    acc, jerk = calcolo_acc_jerk(mass, coord, vel)
    Ein = calcolo_Etot(mass, coord, vel)
    println("energia stato iniziale: $Ein")

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

            if ti % check_freq == 0
                E = calcolo_Etot(mass, coord, vel)
                E_err = abs((E - Ein) / Ein)
                println("tempo: $t, passo: $ti, energia: $E, errore energia: $E_err")
            end

            if ti % snapshot_freq == 0
                A = hcat(mass, coord, vel)
                io["snapshot_$ti"] = (t=t, A=A)
            end


        end

        println("\nrisultati finali")
        E = calcolo_Etot(mass, coord, vel)
        println("tempo: $t, passo: $ti")
        err = diagnostica(Ein, E)
        println("errore relativo dell'energia: $err")

        A = hcat(mass, coord, vel)
        io["final_snapshot"] = (t=t, A=A)
    end

    return A
end


# Interrompe la simulazione quando l'errore supera la soglia minima tollerabile
function simulazione_break(t_parm, tmax, mass, coord, vel, N_dt, errore_soglia, check_freq, snapshot_freq)

    acc, jerk = calcolo_acc_jerk(mass, coord, vel)
    Ein = calcolo_Etot(mass, coord, vel)
    println("energia stato iniziale: $Ein")

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

            E = calcolo_Etot(mass, coord, vel)
            
            E_err = abs((Ein-E)/Ein)

            if E_err >= errore_soglia
                println("Errore sopra soglia! Terminazione anticipata al passo $ti.")
                println("tempo: $t, energia: $E, errore energia: $E_err")
                A = hcat(mass, coord, vel)
                io["snapshot_$ti"] = (t=t, A=A)
                return A
            end

            if ti % check_freq == 0
                println("tempo: $t, passo: $ti energia: $E, errore energia: $E_err")
            end

            if ti % snapshot_freq == 0
                A = hcat(mass, coord, vel)
                io["snapshot_$ti"] = (t=t, A=A)
            end
        end


        println("\nrisultati finali")
        E = calcolo_Etot(mass, coord, vel)
        println("tempo: $t, passo: $ti")
        err = diagnostica(Ein, E)
        println("errore relativo dell'energia: $err")

        A = hcat(mass, coord, vel)
        io["final_snapshot"] = (t=t, A=A)
    end

    return A

end


# simulazioni creata in modo tale che l'errore dell'energia non superi il valore soglia
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

            # questo ciclo, dimezza il passo temporale finché l'errore dell'energia è superiore all'errore soglia
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
using CSV
using DataFrames
using DelimitedFiles
using LinearAlgebra
using Base.Threads
using JLD2


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
    Ek =  calcolo_Ek(mass, vel)
    Ep =  calcolo_Ep(mass, coord)        
    return Ek + Ep
end



function calcolo_energia_totale(mass, coord, vel)
    N = length(mass)
    Ept = zeros(Threads.nthreads())
    Ekt = zeros(Threads.nthreads())

    @threads for i in 1:N
        Ekt[Threads.threadid()] += 0.5 * mass[i] * sum(vel[i, :] .^ 2)
    end

    pairs = [(i, j) for i in 1:N-1 for j in i+1:N]
    @threads for k in 1:length(pairs)
        i, j = pairs[k]
        r = sqrt(sum((coord[i, :] - coord[j, :]).^2))
        Ept[Threads.threadid()] += - mass[i] * mass[j] / r
    end

    return sum(Ekt) + sum(Ept)
end




function diagnostica(Ein, E)
    return abs((Ein - E) / Ein)
end


#################
# FUNZIONI PER IL CALCOLO DELLE MATRICI PER L'INTEGRAZIONE, E PER IL TEMPO DI COLLISIONE
################


function calcolo_tcoll_acc_jerk_bool(mass, coord, vel, x)
    
    N = length(mass)
    acc = zeros(N, 3)
    jerk = zeros(N, 3)
    coll_time = Threads.Atomic{Float64}(Inf)

    
        @threads for i in 1:N
            
            @inbounds for k in 1:i-1 

                r_ij = [coord[k, 1] - coord[i, 1];coord[k, 2] - coord[i, 2]; coord[k, 3] - coord[i, 3]] 
                r2 = sum(r_ij .^ 2)
                r = sqrt(r2) 

                v_ij = [vel[k, 1] - vel[i, 1];vel[k, 2] - vel[i, 2]; vel[k, 3] - vel[i, 3]] 
                 
                # calcolo tempo collisione (r/v) se richiesto
                if x==true
                    v2 = sum(v_ij .^ 2)

                    if v2 != 0 
                        tcoll1 = sqrt(r2 / v2)
                        Threads.atomic_min!(coll_time, tcoll1)
                    end

                end

                
                
                ## calcolo accelerazione e jerk

                acc[i, 1] +=  mass[k] * r_ij[1] / r^3
                acc[i, 2] += mass[k] * r_ij[2] / r^3                        
                acc[i, 3] +=  mass[k] * r_ij[3] / r^3

                v_rad = dot(r_ij, v_ij)

                jerk[i, 1] +=  mass[k] * (v_ij[1] / r^3 - 3 * v_rad * r_ij[1] / r^5)
                jerk[i, 2] +=  mass[k] * (v_ij[2] / r^3 - 3 * v_rad * r_ij[2] / r^5)
                jerk[i, 3] +=  mass[k] * (v_ij[3] / r^3 - 3 * v_rad * r_ij[3] / r^5)

                if x==true

                    a_ij = acc[k, :] - acc[i, :]
                    a2 = sum(a_ij .^ 2)

                    if a2 != 0
                        tcoll2 = (r2 / a2)^0.25
                        Threads.atomic_min!(coll_time, tcoll2)
                    end

                end
                    
            end



            @inbounds for k in i+1:N 

                r_ij = [coord[k, 1] - coord[i, 1];coord[k, 2] - coord[i, 2]; coord[k, 3] - coord[i, 3]] 
                r2 = sum(r_ij .^ 2)
                r = sqrt(r2) 
                
                v_ij = [vel[k, 1] - vel[i, 1];vel[k, 2] - vel[i, 2]; vel[k, 3] - vel[i, 3]] 

                # calcolo tempo collisione (r/v) se richiesto
                if x==true
                    
                    v2 = sum(v_ij .^ 2)

                    if v2 != 0 
                        tcoll1 = sqrt(r2 / v2)
                        Threads.atomic_min!(coll_time, tcoll1)
                    end

                end
                
                
                ## calcolo accelerazione e jerk

                acc[i, 1] +=  mass[k] * r_ij[1] / r^3
                acc[i, 2] += mass[k] * r_ij[2] / r^3                        
                acc[i, 3] +=  mass[k] * r_ij[3] / r^3

                v_rad = dot(r_ij, v_ij)

                jerk[i, 1] +=  mass[k] * (v_ij[1] / r^3 - 3 * v_rad * r_ij[1] / r^5)
                jerk[i, 2] +=  mass[k] * (v_ij[2] / r^3 - 3 * v_rad * r_ij[2] / r^5)
                jerk[i, 3] +=  mass[k] * (v_ij[3] / r^3 - 3 * v_rad * r_ij[3] / r^5)

                if x==true

                    a_ij = acc[k, :] - acc[i, :]
                    a2 = sum(a_ij .^ 2)

                    if a2 != 0
                        tcoll2 = (r2 / a2)^0.25
                        Threads.atomic_min!(coll_time, tcoll2)
                    end

                end

            end
   
        end

        if x==true
            return coll_time[], acc, jerk 
        else
            return acc, jerk
            
        end
        
end


function predict_stept(t_parm, mass,coord, vel)

    N = size(coord, 1)
    new_coord = Matrix{Float64}(undef, N, 3)
    new_vel = Matrix{Float64}(undef, N, 3)

    tcoll, acc, jerk = calcolo_tcoll_acc_jerk_bool(mass, coord, vel, true)
    
    dt=tcoll*t_parm

    @threads for i in 1:N
        new_coord[i, :] .= coord[i, :] .+ vel[i, :] * dt .+ 0.5 * acc[i, :] * dt^2 .+ jerk[i, :] * dt^3 / 6
        new_vel[i, :] .= vel[i, :] .+ acc[i, :] * dt .+ 0.5 * jerk[i, :] * dt^2
    end



    return dt, new_coord, new_vel, acc, jerk    

end


function correct_stept(dt, mass, coord_old, vel_old, coord_new ,vel_new, acc_old, jerk_old)

    N = size(coord, 1)
    acc_new, jerk_new = calcolo_tcoll_acc_jerk_bool(mass, coord_old, vel_old, false)
    
    @threads for i in 1:N
        vel_new[i, :] .= vel_old[i, :] .+ 1/2 * (acc_new[i, :] + acc_old[i, :]) * dt .+ 1/12 * (-jerk_new[i, :] + jerk_old[i, :]) * dt^2
        coord_new[i, :] .= coord_old[i, :] .+ vel_old[i, :] * dt .+ 1/4 * (acc_new[i, :] + acc_old[i, :]) * dt^2 .- 1/12 * jerk_new[i, :] * dt^3
    end

    return coord_new, vel_new

end





function simulazione_nostop(t_parm, tmax, mass, coord, vel, N_dt, check_freq, snapshot_freq)

    Ein = calcolo_Etot(mass, coord, vel)
    println("Energia iniziale: ", Ein)
    
    t=0
    ti=0
    
    A= hcat(mass, coord, vel)

    jldopen("data_output.jld2", "w") do io
        io["snapshot_0"] = (t=t, A=A)

        while ti < N_dt && t < tmax

            dt, new_coord, new_vel, acc, jerk = predict_stept(t_parm, mass,coord, vel)
            coord, vel = correct_stept(dt, mass, coord, vel, new_coord ,new_vel, acc, jerk)

            ti +=1
            t +=dt


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

        A = hcat(mass, coord, vel)
        io["final_snapshot"] = (t=t, A=A)
    
    end

    println("\nrisultati finali")
    E = calcolo_Etot(mass, coord, vel)
    println("tempo: $t, passo: $ti")
    err = diagnostica(Ein, E)
    println("errore relativo dell'energia: $err")


    return A  
end



function simulazione_break(t_parm, tmax, mass, coord, vel, N_dt, check_freq, snapshot_freq, errore_soglia)

    Ein = calcolo_Etot(mass, coord, vel)
    println("Energia iniziale: ", Ein)
    
    t=0
    ti=0
    

    A= hcat(mass, coord, vel)

    jldopen("data_output.jld2", "w") do io
        io["snapshot_0"] = (t=t, A=A)


        while ti < N_dt && t < tmax

            dt, new_coord, new_vel, acc, jerk = predict_stept(t_parm, mass,coord, vel)
            coord, vel = correct_stept(dt, mass, coord, vel, new_coord ,new_vel, acc, jerk)

            ti +=1
            t +=dt

            E = calcolo_Etot(mass, coord, vel)
            E_err = abs((E - Ein) / Ein)

            if E_err > errore_soglia
                println("Errore sopra soglia! Terminazione anticipata al passo $ti.")
                println("tempo: $t, passo: $ti, energia: $E, errore energia: $E_err")
                A = hcat(mass, coord, vel)
                io["snapshot_$ti"] = (t=t, A=A)
                return A
            end
 


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

        A = hcat(mass, coord, vel)
        io["final_snapshot"] = (t=t, A=A)
    
    end

    println("\nrisultati finali")
    E = calcolo_Etot(mass, coord, vel)
    println("tempo: $t, passo: $ti")
    err = diagnostica(Ein, E)
    println("errore relativo dell'energia: $err")


    return A  
end









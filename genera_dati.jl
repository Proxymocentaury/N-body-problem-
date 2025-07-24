using Random, Plots, JLD2

###distribuzioni


#----generazione dati all'interno di una sfera
function random_points_in_sphere(N, R)

    points = zeros(N, 3)
    
    for i in 1:N
        u1, u2, u3 = rand(), rand(), rand()
        
        r = R * u1^(1/3)  # Distribuzione cubica per r
        theta = acos(1 - 2 * u2)  # Distribuzione corretta per theta
        phi = 2 * π * u3  # Distribuzione uniforme per phi
        
        # Conversione in coordinate cartesiane
        x = r * sin(theta)*cos(phi)
        y = r * sin(theta)*sin(phi)
        z = r * cos(theta)
        
        points[i, :] = [x, y, z]
    end
    
    return points
end


function random_points_in_shell(N, R)

    points = zeros(N, 3)
    
    for i in 1:N
         u2, u3 =  rand(), rand()
        
        theta = acos(1 - 2 * u2)  # Distribuzione corretta per theta
        phi = 2 * π * u3  # Distribuzione uniforme per phi
        
        # Conversione in coordinate cartesiane
        x = R * sin(theta)*cos(phi)
        y = R * sin(theta)*sin(phi)
        z = R * cos(theta)
        
        points[i, :] = [x, y, z]
    end
    
    return points
end



#----generazione dati all'interno di un cilindro
function random_points_in_cylinder(N, R, H)
    points = zeros(N, 3)
    
    for i in 1:N
        # Dichiarazione esplicita delle variabili
        local u1, u2, u3, r, theta, z
        
        u1, u2, u3 = rand(), rand(), rand()
        
        r = R * sqrt(u1)       # Distribuzione uniforme in area per r
        theta = 2 * π * u2     # Distribuzione uniforme per theta
        z = H * (u3 - 0.5)     # Distribuzione uniforme per z (centrato in 0)
        
        # Conversione in coordinate cartesiane
        x = r * cos(theta)
        y = r * sin(theta)
        
        points[i, :] = [x, y, z]
    end
    
    return points
end




#------------creazione dati----------------------------------------



# Funzione per generare un set di dati di n corpi con posizioni e velocità casuali
function genera_dati_corpi(N, mass_min, mass_max, pos_min, pos_max, vel_min, vel_max)
    # Limiti per la massa, posizione e velocità
   
    # Inizializza le variabili
    mass = rand(Float64, N) * (mass_max - mass_min) .+ mass_min
    pos = rand(Float64, N, 3) * (pos_max - pos_min) .+ pos_min
    vel = rand(Float64, N, 3) * (vel_max - vel_min) .+ vel_min

    # Unisce massa, posizione e velocità in un'unica matrice A
    A = hcat(mass, pos, vel)

    # Salva i dati in un file JLD2
    @save "corpi.jld2" A

    # Ritorna i dati generati (opzionale, non necessario per l'uso successivo)
    return A
end



function random_mass_in_sphere(N::Int, R, mass_max, mass_min, vel_min, vel_max)
  
    # Inizializza le variabili
    mass = rand(Float64, N) * (mass_max - mass_min) .+ mass_min
    pos = random_points_in_sphere(N, R)
    vel = rand(Float64, N, 3) * (vel_max - vel_min) .+ vel_min

    # Unisce massa, posizione e velocità in un'unica matrice A
    A = hcat(mass, pos, vel)

    # Salva i dati in un file JLD2
    @save "corpi2.jld2" A

    # Ritorna i dati generati (opzionale, non necessario per l'uso successivo)
    return A

end







function random_mass_in_shell(N::Int, R, mass_max, mass_min, vel_min, vel_max)
  
    # Inizializza le variabili
    mass = rand(Float64, N) * (mass_max - mass_min) .+ mass_min
    pos = random_points_in_shell(N, R)
    vel = rand(Float64, N, 3) * (vel_max - vel_min) .+ vel_min

    # Unisce massa, posizione e velocità in un'unica matrice A
    A = hcat(mass, pos, vel)

    # Salva i dati in un file JLD2
    @save "corpi.jld2" A

    # Ritorna i dati generati (opzionale, non necessario per l'uso successivo)
    return A

end



function random_mass_in_cylinder(N, R, h, mass_max, mass_min, vel_min, vel_max)
  
    # Inizializza le variabili
    mass = rand(Float64, N) * (mass_max - mass_min) .+ mass_min
    pos = random_points_in_cylinder(N, R, h)
    vel = rand(Float64, N, 3) * (vel_max - vel_min) .+ vel_min

    # Unisce massa, posizione e velocità in un'unica matrice A
    A = hcat(mass, pos, vel)

    # Salva i dati in un file JLD2
    @save "corpi.jld2" A

    # Ritorna i dati generati (opzionale, non necessario per l'uso successivo)
    return A

end



#----------------------------------------------------

function due_ammassi_sferici(N1, N2, R1, R2, distanza; 
                             mass_min=1.0, mass_max=1.0, 
                             vel_min=0.0, vel_max=0.0, 
                             fileout="corpi.jld2")

    # Primo ammasso centrato in (-distanza/2, 0, 0)
    pos1 = random_points_in_sphere(N1, R1) .+ [-distanza/2 0 0]
    mass1 = rand(Float64, N1) * (mass_max - mass_min) .+ mass_min
    vel1 = rand(Float64, N1, 3) * (vel_max - vel_min) .+ vel_min

    # Secondo ammasso centrato in (+distanza/2, 0, 0)
    pos2 = random_points_in_sphere(N2, R2) .+ [distanza/2 0 0]
    mass2 = rand(Float64, N2) * (mass_max - mass_min) .+ mass_min
    vel2 = rand(Float64, N2, 3) * (vel_max - vel_min) .+ vel_min

    # Unisci tutto
    mass = vcat(mass1, mass2)
    pos = vcat(pos1, pos2)
    vel = vcat(vel1, vel2)

    A = hcat(mass, pos, vel)
    @save fileout A
    return A
end


#----------------------------------------------------

N = 5000

#genera_dati_corpi(N, 1, 1, -50, 50, 0, 0)

random_mass_in_sphere(N, 1, 1/N, 1/N, 0, 0)


#-----plot-----

function plot_dati_generati(file::String)
    # Carica i dati dal file JLD2
    @load file A

    # Estrai le posizioni e le masse
    pos = A[:, 2:4]

    # Plot delle posizioni in 3D
    scatter(pos[:, 1], pos[:, 2], pos[:, 3], markersize=1, label="", title="Distribuzione dei corpi", xlabel="x", ylabel="y", zlabel="z")
end


# Estrai i dati generati e plottali
plot_dati_generati("corpi2.jld2")

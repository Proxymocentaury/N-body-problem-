
using JLD2, Plots, Colors

function crea_gif_3d(file_jld2, output_gif, nfps, bodysize)
    # Apri il file e ottieni la lista delle chiavi (snapshot disponibili)
    io = jldopen(file_jld2, "r")
    snapshots = filter(x -> startswith(x, "snapshot_"), keys(io)) |> collect
    close(io)  # Chiudi il file dopo aver ottenuto le chiavi

    # Ordina gli snapshot per numero
    sorted_snapshots = sort(snapshots, by=x -> parse(Int, split(x, "_")[end]))

    # Inizializza la lista dei tempi e delle coordinate
    tempi = Float64[]
    all_positions = []

    for snap in sorted_snapshots
        io = jldopen(file_jld2, "r")
        data = io[snap]
        close(io)

        push!(tempi, data.t)  # Salva il tempo
        positions = data.A[:, 2:4]  # Estrai X, Y, Z
        push!(all_positions, positions)
    end

    # Numero di corpi nel sistema
    num_corpi = size(all_positions[1], 1)
    colori = distinguishable_colors(num_corpi)  # Colori per ogni corpo

    # Determina il range massimo per scalare gli assi
    max_range = maximum([maximum(abs.(pos)) for pos in all_positions]) * 1.1

    # Creazione dell'animazione
    anim = @animate for i in 1:length(tempi)
        pos = all_positions[i]

        scatter3d(pos[:, 1], pos[:, 2], pos[:, 3],
                  xlabel="X", ylabel="Y", zlabel="Z",
                  title="Tempo: $(round(tempi[i], digits=2))",
                  markersize=bodysize, legend=false, 
                  color=colori, xlims=(-max_range, max_range),
                  ylims=(-max_range, max_range), zlims=(-max_range, max_range))
    end

    # Salva l'animazione come GIF
    gif(anim, output_gif, fps=nfps)
end



function crea_gif_2d(file_jld2, output_gif, nfps, bodysize)
    # Apri il file e ottieni la lista delle chiavi (snapshot disponibili)
    io = jldopen(file_jld2, "r")
    snapshots = filter(x -> startswith(x, "snapshot_"), keys(io)) |> collect
    close(io)  # Chiudi il file dopo aver ottenuto le chiavi

    # Ordina gli snapshot per numero
    sorted_snapshots = sort(snapshots, by=x -> parse(Int, split(x, "_")[end]))

    # Inizializza la lista dei tempi e delle coordinate
    tempi = Float64[]
    all_positions = []

    for snap in sorted_snapshots
        io = jldopen(file_jld2, "r")
        data = io[snap]
        close(io)

        push!(tempi, data.t)  # Salva il tempo
        positions = data.A[:, 2:3]  # Estrai X, Y (assumiamo che le colonne 2 e 3 siano le coordinate)
        push!(all_positions, positions)
    end

    # Numero di corpi nel sistema
    num_corpi = size(all_positions[1], 1)
    colori = distinguishable_colors(num_corpi)  # Colori per ogni corpo

    # Determina il range massimo per scalare gli assi
    max_range = maximum([maximum(abs.(pos)) for pos in all_positions]) * 1.1

    # Creazione dell'animazione
    anim = @animate for i in 1:length(tempi)
        pos = all_positions[i]

        scatter(pos[:, 1], pos[:, 2],
                xlabel="X", ylabel="Y",
                title="Tempo: $(round(tempi[i], digits=2))",
                markersize=bodysize, legend=false, 
                color=colori, xlims=(-max_range, max_range),
                ylims=(-max_range, max_range))
    end

    # Salva l'animazione come GIF
    gif(anim, output_gif, fps=nfps)
end




function crea_gif_2d_traiet(file_jld2, output_gif, nfps, bodysize)
    # Apri il file e ottieni la lista delle chiavi (snapshot disponibili)
    io = jldopen(file_jld2, "r")
    snapshots = filter(x -> startswith(x, "snapshot_"), keys(io)) |> collect
    close(io)  # Chiudi il file dopo aver ottenuto le chiavi

    # Ordina gli snapshot per numero
    sorted_snapshots = sort(snapshots, by=x -> parse(Int, split(x, "_")[end]))

    # Inizializza la lista dei tempi e delle coordinate
    tempi = Float64[]
    all_positions = []

    for snap in sorted_snapshots
        io = jldopen(file_jld2, "r")
        data = io[snap]
        close(io)

        push!(tempi, data.t)  # Salva il tempo
        positions = data.A[:, 2:3]  # Estrai X, Y (colonne 2 e 3)
        push!(all_positions, positions)
    end

    # Numero di corpi nel sistema
    num_corpi = size(all_positions[1], 1)
    colori = distinguishable_colors(num_corpi)  # Colori per ogni corpo

    # Determina il range massimo per scalare gli assi
    max_range = maximum([maximum(abs.(pos)) for pos in all_positions]) * 1.1

    # Creazione dell'animazione
    anim = @animate for i in 1:length(tempi)
        scatter([], [], xlabel="X", ylabel="Y",
                title="Tempo: $(round(tempi[i], digits=2))",
                legend=false, xlims=(-max_range, max_range),
                ylims=(-max_range, max_range))

        # Disegna le traiettorie fino al frame corrente
        for j in 1:num_corpi
            traj_x = [all_positions[k][j, 1] for k in 1:i]
            traj_y = [all_positions[k][j, 2] for k in 1:i]
            plot!(traj_x, traj_y, linestyle=:dot, color=colori[j], linewidth=1.5)  # Linea tratteggiata o punti
        end

        # Disegna i corpi nel frame corrente
        pos = all_positions[i]
        scatter!(pos[:, 1], pos[:, 2], markersize=bodysize, color=colori)
    end

    # Salva l'animazione come GIF
    gif(anim, output_gif, fps=nfps)
end





function crea_gif_3d_traiet(file_jld2, output_gif, nfps, bodysize)
    # Apri il file e ottieni la lista delle chiavi (snapshot disponibili)
    io = jldopen(file_jld2, "r")
    snapshots = filter(x -> startswith(x, "snapshot_"), keys(io)) |> collect
    close(io)  # Chiudi il file dopo aver ottenuto le chiavi

    # Ordina gli snapshot per numero
    sorted_snapshots = sort(snapshots, by=x -> parse(Int, split(x, "_")[end]))

    # Inizializza la lista dei tempi e delle coordinate
    tempi = Float64[]
    all_positions = []

    for snap in sorted_snapshots
        io = jldopen(file_jld2, "r")
        data = io[snap]
        close(io)

        push!(tempi, data.t)  # Salva il tempo
        positions = data.A[:, 2:4]  # Estrai X, Y, Z (colonne 2, 3 e 4)
        push!(all_positions, positions)
    end

    # Numero di corpi nel sistema
    num_corpi = size(all_positions[1], 1)
    colori = distinguishable_colors(num_corpi)  # Colori per ogni corpo

    # Determina il range massimo per scalare gli assi
    max_range = maximum([maximum(abs.(pos)) for pos in all_positions]) * 1.1

    # Creazione dell'animazione
    anim = @animate for i in 1:length(tempi)
        plot3d([], [], [], xlabel="X", ylabel="Y", zlabel="Z", 
               title="Tempo: $(round(tempi[i], digits=2))",
               legend=false, xlims=(-max_range, max_range),
               ylims=(-max_range, max_range),
               zlims=(-max_range, max_range))

        # Disegna le traiettorie fino al frame corrente
        for j in 1:num_corpi
            traj_x = [all_positions[k][j, 1] for k in 1:i]
            traj_y = [all_positions[k][j, 2] for k in 1:i]
            traj_z = [all_positions[k][j, 3] for k in 1:i]
            plot3d!(traj_x, traj_y, traj_z, linestyle=:dot, color=colori[j], linewidth=1.5)
        end

        # Disegna i corpi nel frame corrente
        pos = all_positions[i]
        scatter3d!(pos[:, 1], pos[:, 2], pos[:, 3], markersize=bodysize, color=colori)
    end

    # Salva l'animazione come GIF
    gif(anim, output_gif, fps=nfps)
end




function crea_gif_3d_reso(file_jld2, output_gif, nfps, bodysize)


    io = jldopen(file_jld2, "r")
    snapshots = filter(x -> startswith(x, "snapshot_"), keys(io)) |> collect
    close(io)

    sorted_snapshots = sort(snapshots, by = x -> parse(Int, split(x, "_")[end]))

    tempi = Float64[]
    all_positions = []

    for snap in sorted_snapshots
        io = jldopen(file_jld2, "r")
        data = io[snap]
        close(io)

        push!(tempi, data.t)
        positions = data.A[:, 2:4]
        push!(all_positions, positions)
    end

    num_corpi = size(all_positions[1], 1)
    colori = distinguishable_colors(num_corpi)
    max_range = maximum([maximum(abs.(pos)) for pos in all_positions]) * 1.1

    anim = @animate for i in 1:length(tempi)
        pos = all_positions[i]
        scatter3d(pos[:, 1], pos[:, 2], pos[:, 3],
                  xlabel="X", ylabel="Y", zlabel="Z",
                  title="Tempo: $(round(tempi[i], digits=2))",
                  markersize=bodysize, legend=false,
                  color=colori,
                  xlims=(-max_range, max_range),
                  ylims=(-max_range, max_range),
                  zlims=(-max_range, max_range),
                  size=(1000, 1000))  # Aumenta risoluzione
    end

    gif(anim, output_gif, fps=nfps)  # rimosso dpi
end


function crea_gif_3d_adapt(file_jld2, output_gif, nfps, bodysize)

    io = jldopen(file_jld2, "r")
    snapshots = filter(x -> startswith(x, "snapshot_"), keys(io)) |> collect
    close(io)

    sorted_snapshots = sort(snapshots, by = x -> parse(Int, split(x, "_")[end]))

    tempi = Float64[]
    all_positions = []

    for snap in sorted_snapshots
        io = jldopen(file_jld2, "r")
        data = io[snap]
        close(io)

        push!(tempi, data.t)
        positions = data.A[:, 2:4]
        push!(all_positions, positions)
    end

    num_corpi = size(all_positions[1], 1)
    colori = distinguishable_colors(num_corpi)

    anim = @animate for i in 1:length(tempi)
        pos = all_positions[i]

        # Calcola limiti dinamici per l'asse
        lim = maximum(abs.(pos)) * 1.1

        scatter3d(pos[:, 1], pos[:, 2], pos[:, 3],
                  xlabel = "X", ylabel = "Y", zlabel = "Z",
                  title = "Tempo: $(round(tempi[i], digits=2))",
                  markersize = bodysize, legend = false,
                  color = colori,
                  xlims = (-lim, lim),
                  ylims = (-lim, lim),
                  zlims = (-lim, lim),
                  size = (1000, 1000))  # alta risoluzione
    end

    gif(anim, output_gif, fps = nfps)
end



#crea_gif_2d_traiet("data_output.jld2", "output_gif", 50, 4)


# Esempio di utilizzo
#crea_gif_2d("data_output.jld2", "animazione_2d.gif", 30, 5)


# Esegui la funzione per generare la GIF
#crea_gif_3d("data_output.jld2", "orbite_animazione.gif", 500, 0.5)

# feelpp_mesh_partitioner bash completion script
_feelpp_mesh_partitioner() {
    local cur prev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    opts="--part --dim --shape --order --by-markers --by-markers-desc --ifile --ofile --odir"

    case "${prev}" in
        --dim)
            COMPREPLY=( $(compgen -W "1 2 3" -- ${cur}) )
            return 0
            ;;
        --shape)
            COMPREPLY=( $(compgen -W "Simplex" -- ${cur}) )
            return 0
            ;;
        --order)
            COMPREPLY=( $(compgen -W "1" -- ${cur}) )
            return 0
            ;;
        --ifile)
            _filedir
            return 0
            ;;
        --ofile)
            _filedir
            return 0
            ;;
        --odir)
            _filedir -d
            return 0
            ;;
        --part)
            # If you have a known list of partitions, you can complete them here
            ;;
        --by-markers|--by-markers-desc)
            # If you have a known list of markers, you can complete them here
            ;;
        *)
            ;;
    esac

    COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
}

complete -F _feelpp_mesh_partitioner feelpp_mesh_partitioner


# feelpp_mesh_exporter bash completion script

_feelpp_mesh_exporter() {
    local cur prev opts
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    opts="--dim --shape --scalar_expr --vectorial_expr --gmsh.hsize --gmsh.filename"

    case "${prev}" in
        --dim)
            COMPREPLY=( $(compgen -W "1 2 3" -- ${cur}) )
            return 0
            ;;
        --shape)
            COMPREPLY=( $(compgen -W "Simplex" -- ${cur}) )  # Assuming 'Simplex' is the only shape, add more if available
            return 0
            ;;
        --scalar_expr|--vectorial_expr)
            # You might want to leave these as free text, as expressions can be quite complex
            ;;
        --gmsh.hsize)
            # Accepts numerical values; no specific completion required unless you have specific suggestions
            ;;
        --gmsh.filename)
            _filedir
            return 0
            ;;
        *)
            ;;
    esac

    COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
}

complete -F _feelpp_mesh_exporter feelpp_mesh_exporter


# Shared options completion for all feelpp_toolbox applications
_feelpp_toolbox_common() {
    local toolbox="$1"
    local cur="$2" opts
    # Add all common options shared among toolbox applications here
    opts="--$toolbox.filename --config-file --case --case.dim --case.discretization --$toolbox.ksp-monitor --$toolbox.ksp-view --$toolbox.ksp-type --$toolbox.pc-view --$toolbox.pc-type --$toolbox.snes-monitor --$toolbox.snes-type"

    case "${prev}" in
        --case.dim)
            COMPREPLY=( $(compgen -W "2 3" -- ${cur}) )
            return 0
            ;;
        --case.discretization)
            COMPREPLY=( $(compgen -W "P1 P2 P3" -- ${cur}) )
            return 0
            ;;
        --$toolbox.ksp-monitor|--$toolbox.ksp-view|--$toolbox.pc-view|--$toolbox.snes-monitor)
            COMPREPLY=( $(compgen -W "0 1" -- ${cur}) )
            return 0
            ;;
        --$toolbox.ksp-type|--$toolbox.pc-type|--$toolbox.snes-type)
            # This completion would ideally offer all valid PETSc options, which may be too numerous to list here.
            # You may want to leave these as free text or source a list from somewhere if possible.
            ;;
        --$toolbox.filename|--config-file|--case)
            _filedir
            return 0
            ;;
        *)
            ;;
    esac

    COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
}

# Completion function for feelpp_toolbox_fluid
_feelpp_toolbox_fluid() {
    local cur prev
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Call common options completion
    _feelpp_toolbox_common fluid "${cur}"

    # Add specific options for feelpp_toolbox_fluid
    case "${prev}" in
        # Specific options for fluid toolbox
        *)
            ;;
    esac
}

_feelpp_toolbox_solid() {
    local cur prev
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Call common options completion
    _feelpp_toolbox_common solid "${cur}"

    # Add specific options for feelpp_toolbox_solid
    case "${prev}" in
        # Specific options for fluid toolbox
        *)
            ;;
    esac
}

_feelpp_toolbox_heat() {
    local cur prev
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Call common options completion
    _feelpp_toolbox_common heat "${cur}"

    # Add specific options for feelpp_toolbox_heat
    case "${prev}" in
        # Specific options for fluid toolbox
        *)
            ;;
    esac
}

_feelpp_toolbox_fsi() {
    local cur prev
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Call common options completion
    _feelpp_toolbox_common fsi "${cur}"

    # Add specific options for feelpp_toolbox_fsi
    case "${prev}" in
        # Specific options for fluid toolbox
        *)
            ;;
    esac
}

_feelpp_toolbox_heatfluid() {
    local cur prev
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Call common options completion
    _feelpp_toolbox_common heatfluid "${cur}"

    # Add specific options for feelpp_toolbox_heatfluid
    case "${prev}" in
        # Specific options for fluid toolbox
        *)
            ;;
    esac
}

_feelpp_toolbox_thermoelectric() {
    local cur prev
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Call common options completion
    _feelpp_toolbox_common thermoelectric "${cur}"

    # Add specific options for feelpp_toolbox_thermoelectric
    case "${prev}" in
        # Specific options for fluid toolbox
        *)
            ;;
    esac
}

_feelpp_toolbox_electric() {
    local cur prev
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Call common options completion
    _feelpp_toolbox_common electric "${cur}"

    # Add specific options for feelpp_toolbox_electric
    case "${prev}" in
        # Specific options for fluid toolbox
        *)
            ;;
    esac
}
# Similar functions for other toolbox applications...
# _feelpp_toolbox_solid(), _feelpp_toolbox_fsi(), etc.

# Register the completion functions
complete -F _feelpp_toolbox_fluid feelpp_toolbox_fluid
complete -F _feelpp_toolbox_solid feelpp_toolbox_solid
complete -F _feelpp_toolbox_fsi feelpp_toolbox_fsi_2d
complete -F _feelpp_toolbox_fsi feelpp_toolbox_fsi_3d
complete -F _feelpp_toolbox_heat feelpp_toolbox_heat
complete -F _feelpp_toolbox_heatfluid feelpp_toolbox_heatfluid
complete -F _feelpp_toolbox_thermoelectric feelpp_toolbox_thermoelectric

#!/usr/bin/env python3

"""
make_etk_par.py

Generate an Einstein Toolkit BHNS parfile from an Elliptica
BHNS_properties.txt file.

The script reads:

  1. BHNS_properties.txt
  2. an ETK parfile template

and fills the template.

The Elliptica properties file is expected to contain entries such as:

    Project = ...
    BH_irreducible_mass = ...
    NS_baryonic_mass = ...
    BHNS_separation = ...

and the final/converged centre-of-mass positions:

    BH_x_CM = ...
    BH_y_CM = ...
    NS_x_CM = ...
    NS_y_CM = ...

If BHNS_properties.txt contains multiple "# iteration = N" blocks,
parse_key_value_file keeps overwriting each key as it scans forward,
so the LAST block in the file naturally wins -- no explicit block
splitting needed, as long as keys repeat identically across blocks.

EOS
---

The script supports:

    polytropic
    piecewise_polytropic

For a piecewise-polytropic EOS with

    NS_EoS_Gamma   = [G0, G1, G2, G3]
    NS_EoS_rho0_th = [0, rho1, rho2, rho3]

the actual transition densities are:

    rho1, rho2, rho3

i.e. rho0_th[1:].

WARNING: the piecewise-polytropic parameter names used below
(IllinoisGRMHD::K_ppoly_tab0, rho_ppoly_tab_in, Gamma_ppoly_tab_in,
Gamma_th; EOS_Omni::n_pieces, hybrid_k0, hybrid_rho, hybrid_gamma,
hybrid_gamma_th) are NOT verified against the actual thorn param.ccl
files. Check them against arrangements/GRHayLET/IllinoisGRMHD/param.ccl
and arrangements/EinsteinEOS/EOS_Omni/param.ccl before trusting a
piecewise-EOS run.

USAGE
-----

    python make_etk_par.py \
        --elliptica-par path/to/BHNS_properties.txt \
        --template ETK_template.par \
        --output ETK_run.par
"""

import argparse
import os
import re
import sys


# ---------------------------------------------------------------------
# 1. Generic parsing utilities
# ---------------------------------------------------------------------

def parse_key_value_file(path):
    """
    Parse a simple Elliptica-style text file containing

        key = value

    into a dictionary of strings.

    Lines beginning with '#' and blank lines are ignored.
    """

    params = {}

    line_re = re.compile(
        r"^\s*([A-Za-z0-9_\/]+)\s*=\s*(.*?)\s*$"
    )

    with open(path, "r") as f:

        for raw_line in f:

            # Remove comments.
            line = raw_line.split("#", 1)[0].rstrip()

            if not line.strip():
                continue

            m = line_re.match(line)

            if m:

                key = m.group(1)
                value = m.group(2)

                params[key] = value

    return params


def parse_scalar(raw):
    """
    Parse a scalar float from a string.

    Examples
    --------
    '80.'       -> 80.0
    '+55.0'     -> 55.0
    '[8.9e-2]'  -> 0.089
    """

    raw = raw.strip()

    try:

        return float(raw)

    except ValueError:

        pass

    # Extract the first floating-point number.
    m = re.search(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?",
        raw
    )

    if m:

        return float(m.group(0))

    raise ValueError(
        f"Could not parse scalar from {raw!r}"
    )


def parse_array(raw):
    """
    Parse an Elliptica array such as

        [1.0, 2.0, 3.0]

    into a list of floats.
    """

    raw = raw.strip()

    # Remove brackets if present.
    if raw.startswith("[") and raw.endswith("]"):

        raw = raw[1:-1]

    if not raw.strip():

        return []

    values = []

    for item in raw.split(","):

        item = item.strip()

        if not item:
            continue

        values.append(parse_scalar(item))

    return values


def get_scalar(params, key, default=None):

    if key not in params:

        if default is not None:
            return default

        raise KeyError(
            f"Missing required Elliptica parameter: {key}"
        )

    return parse_scalar(params[key])


def get_array(params, key, default=None):

    if key not in params:

        if default is not None:
            return default

        raise KeyError(
            f"Missing required Elliptica parameter: {key}"
        )

    return parse_array(params[key])


# ---------------------------------------------------------------------
# 2. Elliptica positions
# ---------------------------------------------------------------------
def read_positions(params):
    """
    Read the individual BH/NS centre-of-mass positions and the
    physical BHNS system centre-of-mass position from
    BHNS_properties.txt.

    Elliptica reports:

        BH_x_CM, BH_y_CM
        NS_x_CM, NS_y_CM

    as the individual object positions, and

        BHNS_x_CM, BHNS_y_CM

    as the centre-of-mass position of the complete binary.

    The ETK initial-data setup is constructed in the BHNS COM frame,
    so the positions supplied to the ETK template are

        x_ETK = x_Elliptica - BHNS_x_CM
        y_ETK = y_Elliptica - BHNS_y_CM

    This is only a coordinate translation. It does not modify any
    physical initial-data quantities such as masses, momenta,
    separation, spins, etc.

    Returns
    -------
    dict
        {
            "BH_POS_X": ...,
            "BH_POS_Y": ...,
            "NS_POS_X": ...,
            "NS_POS_Y": ...,
            "BHNS_COM_X": ...,
            "BHNS_COM_Y": ...
        }
    """

    required = [
        "BH_x_CM",
        "BH_y_CM",
        "NS_x_CM",
        "NS_y_CM",
        "BHNS_x_CM",
        "BHNS_y_CM",
    ]

    missing = [
        key
        for key in required
        if key not in params
    ]

    if missing:
        raise KeyError(
            "BHNS_properties.txt is missing the required "
            f"position entries: {missing}"
        )

    # Individual object positions in the original Elliptica frame.
    bh_x = get_scalar(params, "BH_x_CM")
    bh_y = get_scalar(params, "BH_y_CM")

    ns_x = get_scalar(params, "NS_x_CM")
    ns_y = get_scalar(params, "NS_y_CM")

    # Physical centre of mass of the complete BHNS system.
    system_x_cm = get_scalar(params, "BHNS_x_CM")
    system_y_cm = get_scalar(params, "BHNS_y_CM")

    # Translate the coordinate origin to the BHNS COM.
    bh_x_com = bh_x - system_x_cm
    bh_y_com = bh_y - system_y_cm

    ns_x_com = ns_x - system_x_cm
    ns_y_com = ns_y - system_y_cm

    return {
        "BH_POS_X": bh_x_com,
        "BH_POS_Y": bh_y_com,
        "NS_POS_X": ns_x_com,
        "NS_POS_Y": ns_y_com,
        "BHNS_COM_X": system_x_cm,
        "BHNS_COM_Y": system_y_cm,
    }


# ---------------------------------------------------------------------
# 3. EOS
# ---------------------------------------------------------------------

def make_eos_block(params):
    """
    Construct the complete EOS parameter block from Elliptica
    properties.

    Supported EOS types
    --------------------

    polytropic
    piecewise_polytropic

    The returned string replaces the complete {{EOS_BLOCK}}
    placeholder in the ETK template.
    """

    eos_type = params.get(
        "NS_EoS_type",
        ""
    ).strip().lower()

    if eos_type == "polytropic":

        return make_polytropic_eos_block(params)

    elif eos_type == "piecewise_polytropic":

        return make_piecewise_polytropic_eos_block(params)

    else:

        raise ValueError(
            f"Unsupported NS_EoS_type = "
            f"{params.get('NS_EoS_type')!r}. "
            "Currently supported: 'polytropic' and "
            "'piecewise_polytropic'."
        )


def make_polytropic_eos_block(params):

    K_values = get_array(
        params,
        "NS_EoS_K0"
    )

    Gamma_values = get_array(
        params,
        "NS_EoS_Gamma"
    )

    if len(K_values) != 1:
        raise ValueError(
            "Polytropic EOS should contain exactly one "
            f"NS_EoS_K0 value, got {K_values}"
        )

    if len(Gamma_values) != 1:
        raise ValueError(
            "Polytropic EOS should contain exactly one "
            f"NS_EoS_Gamma value, got {Gamma_values}"
        )

    K = K_values[0]
    Gamma = Gamma_values[0]

    block = f"""
########################################################################
#       EOS: SINGLE POLYTROPE
########################################################################

Eos_omni::poly_gamma = {Gamma:.15e}
Eos_omni::poly_k     = {K:.15e}

illinoisgrmhd::neos     = 1
illinoisgrmhd::gamma_th = {Gamma:.15e}
illinoisgrmhd::k_poly   = {K:.15e}

id_converter_ilgrmhd::gamma_initial = {Gamma:.15e}
id_converter_ilgrmhd::k_initial     = {K:.15e}
"""

    return block.strip()


def make_piecewise_polytropic_eos_block(params):
    """
    Construct the EOS block for a piecewise-polytropic EOS
    using the current GRHayL/GRHayLib interface.

    Elliptica supplies:

        NS_EoS_K0       = [K0]
        NS_EoS_Gamma    = [Gamma0, Gamma1, ..., GammaN-1]
        NS_EoS_rho0_th  = [0, rho1, rho2, ..., rhoN-1]

    GRHayLib expects:

        GRHayLib::neos
        GRHayLib::k_ppoly0
        GRHayLib::Gamma_ppoly_in[i]   for i = 0 ... N-1
        GRHayLib::rho_ppoly_in[i]     for i = 0 ... N-2
        GRHayLib::Gamma_th

    The first Elliptica rho value is the lower boundary rho=0
    and is therefore NOT passed to rho_ppoly_in.

    GRHayLib computes all K_i for i > 0 automatically by
    enforcing pressure continuity.
    """

    K_values = get_array(
        params,
        "NS_EoS_K0"
    )

    Gamma_values = get_array(
        params,
        "NS_EoS_Gamma"
    )

    rho_values = get_array(
        params,
        "NS_EoS_rho0_th"
    )

    # -------------------------------------------------------------
    # Validate K0
    # -------------------------------------------------------------

    if len(K_values) != 1:
        raise ValueError(
            "Piecewise-polytropic EOS should contain exactly one "
            f"NS_EoS_K0 value, got {K_values}"
        )

    # -------------------------------------------------------------
    # Number of pieces
    # -------------------------------------------------------------

    n_pieces = len(Gamma_values)

    if n_pieces < 1:
        raise ValueError(
            "Piecewise-polytropic EOS must contain at least one "
            "Gamma value."
        )

    # -------------------------------------------------------------
    # Validate rho array
    #
    # Elliptica convention:
    #
    #   [0, rho1, rho2, ..., rho_(N-1)]
    #
    # Therefore len(rho_values) == n_pieces.
    # -------------------------------------------------------------

    if len(rho_values) != n_pieces:
        raise ValueError(
            "Inconsistent piecewise-polytropic EOS:\n"
            f"  number of Gamma values = {n_pieces}\n"
            f"  number of rho0_th values = {len(rho_values)}\n"
            f"  Gamma = {Gamma_values}\n"
            f"  rho0_th = {rho_values}"
        )

    if abs(rho_values[0]) > 1e-14:
        raise ValueError(
            "Expected the first NS_EoS_rho0_th entry to be zero "
            "for a piecewise-polytropic EOS, got "
            f"{rho_values[0]}"
        )

    # -------------------------------------------------------------
    # Extract physical parameters
    # -------------------------------------------------------------

    K0 = K_values[0]

    # The first entry is rho = 0.
    # GRHayLib::rho_ppoly_in contains only the N-1
    # transition densities.
    transition_rho = rho_values[1:]

    # Thermal Gamma.
    #
    # Keep the existing choice used by your setup.
    Gamma_th = 2.0

    # -------------------------------------------------------------
    # Construct GRHayLib block
    # -------------------------------------------------------------

    block_lines = []

    block_lines.append(
        "########################################################################"
    )
    block_lines.append(
        "#       EOS: PIECEWISE POLYTROPE"
    )
    block_lines.append(
        "#       GRHayL / GRHayLib interface"
    )
    block_lines.append(
        "########################################################################"
    )
    block_lines.append("")

    block_lines.append(
        "# Number of cold piecewise-polytropic pieces"
    )
    block_lines.append(
        f"GRHayLib::neos = {n_pieces}"
    )
    block_lines.append("")

    block_lines.append(
        "# K for the first, lowest-density piece."
    )
    block_lines.append(
        "# GRHayLib computes all remaining K_i from continuity."
    )
    block_lines.append(
        f"GRHayLib::k_ppoly0 = {K0:.15e}"
    )
    block_lines.append("")

    block_lines.append(
        "# Gamma_i for each polytropic piece"
    )

    for i, gamma in enumerate(Gamma_values):
        block_lines.append(
            f"GRHayLib::Gamma_ppoly_in[{i}] = "
            f"{gamma:.15e}"
        )

    block_lines.append("")

    block_lines.append(
        "# Density boundaries between consecutive pieces"
    )

    for i, rho in enumerate(transition_rho):
        block_lines.append(
            f"GRHayLib::rho_ppoly_in[{i}] = "
            f"{rho:.15e}"
        )

    block_lines.append("")

    block_lines.append(
        "# Thermal Gamma"
    )
    block_lines.append(
        f"GRHayLib::Gamma_th = {Gamma_th:.15e}"
    )

    return "\n".join(block_lines)


# ---------------------------------------------------------------------
# 4. Grid / AMR heuristics
# ---------------------------------------------------------------------

def estimate_grid(separation, m_bh, m_ns, max_refinement_levels=8):
    """
    Simple starting-point grid prescription.

    These are heuristics and should be tuned against the actual
    convergence/cost requirements of the simulation.

    CRITICAL CONSTRAINT: Carpet requires the domain size divided by
    the coarse grid spacing to be an exact integer at every
    refinement level it constructs via convergence_factor=2 -- not
    just at the coarsest level. So the number of coarse-grid
    half-cells (grid_outer / dx) is forced to be a multiple of
    2**max_refinement_levels, and dx is then derived so that
    grid_outer / dx is exact by construction, rather than picking
    grid_outer and dx independently and hoping they divide evenly.
    """

    total_mass = m_bh + m_ns

    target_outer = separation * 8.4
    target_dx = max(total_mass * 0.18, 0.5)

    divisor = 2 ** max_refinement_levels
    n_half_cells = round((target_outer / target_dx) / divisor) * divisor
    if n_half_cells < divisor:
        n_half_cells = divisor

    grid_dx = target_outer / n_half_cells       # exact by construction
    grid_dx = round(grid_dx, 8)
    grid_outer = round(grid_dx * n_half_cells, 8)

    region_radii = [

        round(
            grid_outer * 0.446,
            0
        ),

        round(
            grid_outer * 0.223,
            0
        ),

        round(
            grid_outer * 0.111,
            0
        ),

        round(
            grid_outer * 0.0558,
            1
        ),

        round(
            grid_outer * 0.0297,
            1
        ),

        round(
            grid_outer * 0.023,
            0
        ),

        round(
            grid_outer * 0.015,
            0
        ),
    ]

    ah_init_radius = round(
        m_bh * 0.9,
        3
    )

    ah_max_radius = round(
        m_bh * 1.4,
        3
    )

    return {

        "GRID_OUTER": grid_outer,

        "GRID_DX": grid_dx,

        "REGION_R1": region_radii[0],
        "REGION_R2": region_radii[1],
        "REGION_R3": region_radii[2],
        "REGION_R4": region_radii[3],
        "REGION_R5": region_radii[4],
        "REGION_R6": region_radii[5],
        "REGION_R7": region_radii[6],

        "AH_INIT_RADIUS": ah_init_radius,
        "AH_MAX_RADIUS": ah_max_radius,
    }


# ---------------------------------------------------------------------
# 5. Template filling
# ---------------------------------------------------------------------

def fill_template(template_text, tokens):
    """
    Replace {{TOKEN}} occurrences.

    Raises if a template token has no supplied value.

    Warns about supplied tokens that are not used.
    """

    used = set()

    def replace(match):

        key = match.group(1)

        used.add(key)

        if key not in tokens:

            raise KeyError(
                f"No value supplied for template token "
                f"{{{{{key}}}}}"
            )

        return str(tokens[key])

    result = re.sub(
        r"\{\{([A-Za-z0-9_]+)\}\}",
        replace,
        template_text
    )

    unused = set(tokens.keys()) - used

    if unused:

        print(
            "[warn] unused tokens (not present in template): "
            f"{sorted(unused)}",
            file=sys.stderr
        )

    return result


# ---------------------------------------------------------------------
# 6. Main
# ---------------------------------------------------------------------

def main():

    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    ap.add_argument(
        "--elliptica-par",
        required=True,
        help=(
            "Elliptica BHNS_properties.txt file containing "
            "the input properties and converged positions"
        )
    )

    ap.add_argument(
        "--template",
        default="ETK_template.par",
        help="ETK parfile template"
    )

    ap.add_argument(
        "--output",
        required=True,
        help="Output ETK parfile"
    )

    ap.add_argument(
        "--run-title",
        default=None,
        help=(
            "Override Cactus::cctk_run_title "
            "(defaults to Project from BHNS_properties.txt)"
        )
    )

    args = ap.parse_args()

    # ---------------------------------------------------------------
    # Check input properties file
    # ---------------------------------------------------------------

    properties_path = os.path.abspath(
        args.elliptica_par
    )

    # ---------------------------------------------------------------
    # Checkpoint path for the ETK template
    #
    # The checkpoint itself is NOT opened or parsed. We only construct
    # its expected path from the directory containing
    # BHNS_properties.txt.
    # ---------------------------------------------------------------

    checkpoint_path = os.path.join(
        os.path.dirname(properties_path),
        "checkpoint.dat"
    )

    if not os.path.isfile(properties_path):

        raise FileNotFoundError(
            f"Elliptica properties file not found: "
            f"{properties_path}"
        )

    # ---------------------------------------------------------------
    # Elliptica input properties
    # ---------------------------------------------------------------

    params = parse_key_value_file(
        properties_path
    )

    print(
        "[info] Reading Elliptica properties from:",
        file=sys.stderr
    )

    print(
        f"       {properties_path}",
        file=sys.stderr
    )

    # ---------------------------------------------------------------
    # Basic properties
    # ---------------------------------------------------------------

    project = params.get(
        "Project",
        "elliptica_bhns_run"
    )

    run_title = (
        args.run_title
        or project
    )

    m_bh = get_scalar(
        params,
        "BH_irreducible_mass"
    )

    m_ns = get_scalar(
        params,
        "NS_baryonic_mass"
    )

    separation = get_scalar(
        params,
        "BHNS_separation"
    )

    # ---------------------------------------------------------------
    # Elliptica positions
    # ---------------------------------------------------------------

    positions = read_positions(
        params
    )

    print(
        "[info] Using Elliptica CM positions from "
        "BHNS_properties.txt:",
        file=sys.stderr
    )

    print(
        "[info] Elliptica individual CM positions:",
        file=sys.stderr
    )

    print(
        "       BH = "
        f"({get_scalar(params, 'BH_x_CM'):.15f}, "
        f"{get_scalar(params, 'BH_y_CM'):.15f})",
        file=sys.stderr
    )

    print(
        "       NS = "
        f"({get_scalar(params, 'NS_x_CM'):.15f}, "
        f"{get_scalar(params, 'NS_y_CM'):.15f})",
        file=sys.stderr
    )

    print(
        "[info] Elliptica BHNS system COM:",
        file=sys.stderr
    )

    print(
        "       COM = "
        f"({positions['BHNS_COM_X']:.15f}, "
        f"{positions['BHNS_COM_Y']:.15f})",
        file=sys.stderr
    )

    print(
        "[info] Using BHNS COM frame for ETK:",
        file=sys.stderr
    )

    print(
        "       BH = "
        f"({positions['BH_POS_X']:.15f}, "
        f"{positions['BH_POS_Y']:.15f})",
        file=sys.stderr
    )

    print(
        "       NS = "
        f"({positions['NS_POS_X']:.15f}, "
        f"{positions['NS_POS_Y']:.15f})",
        file=sys.stderr
    )

    # ---------------------------------------------------------------
    # EOS
    # ---------------------------------------------------------------

    eos_block = make_eos_block(
        params
    )

    print(
        "[info] Elliptica EOS block:",
        file=sys.stderr
    )

    print(
        eos_block,
        file=sys.stderr
    )

    # ---------------------------------------------------------------
    # Grid
    # ---------------------------------------------------------------

    grid = estimate_grid(
        separation,
        m_bh,
        m_ns
    )

    # ---------------------------------------------------------------
    # Token dictionary
    # ---------------------------------------------------------------

    tokens = {

        "RUN_TITLE": run_title,

        "ID_TYPE": "BHNS",

        "CHECKPOINT_PATH": checkpoint_path,

        # -----------------------------------------------------------
        # Converged Elliptica CM positions
        # -----------------------------------------------------------

        "BH_POS_X": (
            f"{positions['BH_POS_X']:.15f}"
        ),

        "BH_POS_Y": (
            f"{positions['BH_POS_Y']:.15f}"
        ),

        "NS_POS_X": (
            f"{positions['NS_POS_X']:.15f}"
        ),

        "NS_POS_Y": (
            f"{positions['NS_POS_Y']:.15f}"
        ),

        # -----------------------------------------------------------
        # EOS
        # -----------------------------------------------------------

        "EOS_BLOCK": eos_block,

        # -----------------------------------------------------------
        # Grid
        # -----------------------------------------------------------

        **grid,
    }

    # ---------------------------------------------------------------
    # Fill template
    # ---------------------------------------------------------------

    template_path = os.path.abspath(
        args.template
    )

    if not os.path.isfile(template_path):

        raise FileNotFoundError(
            f"ETK template not found: "
            f"{template_path}"
        )

    with open(template_path, "r") as f:

        template_text = f.read()

    filled = fill_template(
        template_text,
        tokens
    )

    # ---------------------------------------------------------------
    # Write output
    # ---------------------------------------------------------------

    output_path = os.path.abspath(
        args.output
    )

    with open(output_path, "w") as f:

        f.write(filled)

    print(
        f"[info] wrote {output_path}",
        file=sys.stderr
    )


if __name__ == "__main__":
    main()

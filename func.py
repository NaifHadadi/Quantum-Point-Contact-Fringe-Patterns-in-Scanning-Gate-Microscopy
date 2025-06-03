
import kwant
import numpy as np
import matplotlib.pyplot as plt


def make_system(system_length=1):
    """Create a Kwant system with a central potential and tip potential.
    
    Args:
        system_length: Length of the system (default: 1)
    
    Returns:
        A Kwant system with attached leads
    """
    system = kwant.Builder()
    lattice = kwant.lattice.square(norbs=1)

    def central_potential(site, gate_voltage=0):
        """Potential at the central site."""
        return gate_voltage

    def tip_potential(site, tip_voltage=1.0, x_position=1): 
        """Potential at the tip position, zero elsewhere."""
        return tip_voltage if site.pos == (x_position, 0) else 0
    
    # Define the system geometry
    # Left part of the system
    system[(lattice(-1, y) for y in range(-8, 9))] = 0
    # Central site with gate potential
    system[lattice(0, 0)] = central_potential
    # Right part with tip potential
    system[(lattice(x, y) for x in range(1, system_length + 1) 
            for y in range(-8, 9))] = tip_potential
    
    # Set hopping between neighboring sites
    system[lattice.neighbors()] = -1

    # Define and attach leads
    symmetry = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(symmetry)
    lead[(lattice(1, y) for y in range(-8, 9))] = 0
    lead[lattice.neighbors()] = -1

    system.attach_lead(lead)
    system.attach_lead(lead.reversed(), add_cells=0)
    
    return system



def make_qpc_with_tip(w=8, length=30, tip_pos=(10, 0), tip_voltage=1.0):
    """Construct a QPC system with a scanning gate tip at a defined position."""
    lat = kwant.lattice.square(norbs=1)
    sys = kwant.Builder()

    # Central scattering region
    def potential(site, gate_voltage=0):
        return gate_voltage

    def tip(site):
        return tip_voltage if site.pos == tip_pos else 0

    # Left lead sites
    sys[(lat(-1, y) for y in range(-w, w + 1))] = 0

    # Central region
    sys[lat(0, 0)] = potential

    # Right scattering region (where tip moves)
    sys[(lat(x, y) for x in range(1, length + 1)
         for y in range(-w, w + 1))] = tip

    # Hopping terms
    sys[lat.neighbors()] = -1

    # Lead definition
    sym = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(sym)
    lead[(lat(1, y) for y in range(-w, w + 1))] = 0
    lead[lat.neighbors()] = -1

    # Attach leads
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed(), add_cells=0)

    return sys, lat






def study1(w=10, pos=(30, 0)):
    sys = kwant.Builder()
    lat = kwant.lattice.square()

    def Pot(site, Vg): return Vg
    def Hop(site1, site2, tc): return tc
    def Vtip(site, v): return v if site.pos == pos else 0

    sys[lat(0, 0)] = Pot
    sys[(lat(-1, y) for y in range(-w, w + 1))] = 0
    sys[(lat(1, y) for y in range(-w, w + 1))] = 0
    sys[lat.neighbors()] = -1
    sys[lat(0, 0), lat(-1, 0)] = Hop
    sys[lat(0, 0), lat(1, 0)] = Hop
    sys[lat(*pos)] = Vtip

    sym = kwant.TranslationalSymmetry((-1, 0))
    lead = kwant.Builder(sym)
    lead[(lat(-2, y) for y in range(-w, w + 1))] = 0
    lead[lat.neighbors()] = -1
    sys.attach_lead(lead)
    sys.attach_lead(lead.reversed(), add_cells=0)

    return sys

def make_and_run_study(tc_values, Vg_ranges, divisors, energy=-3.8, N=100, w=10, pos=(30, 0)):
    sys = study1(w=w, pos=pos)
    sysf = sys.finalized()
    plt.figure(figsize=(7, 5))

    for tc, Vg_range, divisor in zip(tc_values, Vg_ranges, divisors):
        if N is not None:
            Vg_values = np.linspace(Vg_range[0], Vg_range[1], N)
        else:
            Vg_values = [x / divisor for x in range(int(Vg_range[0]*divisor), int(Vg_range[1]*divisor))]

        T_list = []
        for Vg in Vg_values:
            params = dict(Vg=Vg, tc=tc, v=0)
            smatrix = kwant.smatrix(sysf, energy=energy, params=params)
            T_list.append(smatrix.transmission(0, 1))

        plt.plot(Vg_values, T_list, label=f"tc = {tc}", linewidth=2)

    plt.ylim(0,1.2)
    plt.xlabel("Vg", fontsize=16)
    plt.ylabel("Transmission", fontsize=16)
    plt.title(f"Transmission vs Vg at E = {energy}")
    plt.legend()
    plt.tight_layout()
    plt.show()

import radis
import astropy.units as u 
from  astropy.units import cds 
from src.models.payload import Payload

ALLOWED_UNITS = {
    "m": u.m,
    "cm": u.cm,
    "mm": u.mm,
    "nm": u.nm,
    "um": u.um,
    "1/u.cm": 1 / u.cm,
    # Add other allowed units here
}

def calculate_spectrum(payload: Payload):
    print(">> Payload : ")
    print(payload)
    
    if payload.wavelength_units not in ALLOWED_UNITS:
        raise ValueError(f"Invalid wavelength unit: {payload.wavelength_units}")
    if payload.pressure_units not in ALLOWED_UNITS:
        raise ValueError(f"Invalid pressure unit: {payload.pressure_units}")
    if payload.path_length_units not in ALLOWED_UNITS:
        raise ValueError(f"Invalid path length unit: {payload.path_length_units}")
    
    spectrum = radis.calc_spectrum(
        payload.min_wavenumber_range * ALLOWED_UNITS[payload.wavelength_units],
        payload.max_wavenumber_range * ALLOWED_UNITS[payload.wavelength_units],
        molecule=[species.molecule for species in payload.species],
        mole_fraction={
            species.molecule: species.mole_fraction for species in payload.species
        },
        # TODO: Hard-coding "1,2,3" as the isotopologue for the time-being
        isotope={species.molecule: "1,2,3" for species in payload.species},
        pressure=payload.pressure * ALLOWED_UNITS[payload.pressure_units],
        Tgas=payload.tgas,
        Tvib=payload.tvib,
        Trot=payload.trot,
        path_length=payload.path_length * ALLOWED_UNITS[payload.path_length_units],
        export_lines=False,
        wstep="auto",
        databank=payload.database,
        use_cached=True,
    )
    return spectrum

import numpy as np
from os import path
import pandas as pd


def gaia_exoplanets_cross(gaia_filename, crossmatch_dir, save_gaia_id=False, return_data=False, save_spherical=True):
    """
    Cross match Gaia dataset with NASA exoplanet dataset.

    :param gaia_filename: Name of the file to read Gaia information from
    :param crossmatch_dir: Path to crossmatch directory
    :param save_gaia_id: Bool to save source_id values in a separate file
    :param return_data: Return crossmatched gaia data for further use
    :param save_spherical: Save crossmatched data to a CSV format

    :return: Cross matched dataset
    """

    # Path to downloaded datasets
    datasets_dir = "data/initial_datasets"

    # Read Exoplanets data: [SR] changed skiprows from 28 to 24 bc exoplanet datafile seems to have changed
    exoplanets = pd.read_csv(path.join(datasets_dir, "exoplanets.csv"), skiprows=24,
                             usecols=["pl_name", "hostname", "gaia_id", "pl_orbper", "pl_orbsmax", "pl_bmasse"]) #path.joinpath, *paths) joins directories togeter so that you can create a new directory that is a floder withing a folder etc.
    #exoplanets = pd.read_csv(path.join(datasets_dir, "exoplanets.csv"), skiprows=28,
                             #usecols=["pl_name", "hostname", "gaia_id", "pl_orbper", "pl_orbsmax", "pl_bmasse"])
    # Process Exoplanets data
    exoplanets.dropna(subset=["gaia_id"], inplace=True)
    exoplanets["source_id"] = exoplanets["gaia_id"].str.rsplit(" ", n=1, expand=True)[1].astype("int64") #The source_id is the first part of the gaia_id (split up if it has spaces in it). Splits the source_id on the whitespace. Max one split can occur. Expands the split strings into separate columns (i.e. the source_id gets separated into two columns if it has a space in it.
    exoplanets.drop(["gaia_id"], axis=1, inplace=True) #deletes all entries in the gaia_id column and replaces them with 'None'
    exoplanets["Host"] = exoplanets["hostname"].str.replace(" ", "") #the exoplanets' Host is the hostname without spaces
    exoplanets.drop_duplicates(subset=["Host"], inplace=True)

    # Read Gaia data
    gaia = pd.read_csv(path.join(datasets_dir, gaia_filename))

    # Add Gaia information to Exoplanet hosts
    exoplanets = pd.merge(exoplanets, gaia, on="source_id") #gaia information is added to the list of exoplanets. The gaia data is matched to the exoplanet list by the source_id
    exoplanets.drop(["pl_name", "hostname"], axis=1, inplace=True) #deletes the planet name and hostname from the exoplanets df. The source_id remains

    # Remove exoplanet hosts from Gaia. Why?
    gaia = gaia[~gaia["source_id"].isin(exoplanets["source_id"])]
    # Further reduction of the data
    gaia = gaia[4.5 < gaia["parallax"] / gaia["parallax_error"]] #removes stars with with parallax/parallax error > 4.5
    # like a signal to noise cut

    # Concatenate exoplanet hosts back, however at the top of the dataframe. This way for testing purposes we later
    # iterate only over first 1065 entries that are exoplanet hosts.
    gaia = pd.concat([exoplanets, gaia]) # adding the exoplanet list back into the gaia df, at the top

    # Calculate distance in pc and drop any stars with negative or null distance
    gaia["distance_pc"] = (1. / gaia["parallax"]) * 1000 #closely aligned sources are only occasionally resolved in Gaia, confusion in observation-to-source matching can lead to spurious parallax values which are either very large or have a negative value very far away from zero
    gaia = gaia[gaia["distance_pc"] > 0] #returns all of gaia for which distance_pc > 0 and overwrites the gaia df with it. Gets rid of all entries where distance_pc <= 0. For these entries, the solution returned by gaia is unphysical so we want to ditch it

    # Convert from degrees to pc
    gaia["ra"] = (gaia["ra"] * np.pi) / 180.
    gaia["dec"] = (gaia["dec"] * np.pi) / 180.

    if save_gaia_id: #if we have entered save_gaia_id=True into the function
        gaia[["source_id", "Host"]].to_csv(path.join(crossmatch_dir, f"{gaia_filename.split('.')[0]}_star_labels.csv"), index=False) #create a csv containing the source_id and Host

    # Drop all unnecessary data leaving only 6 coordinates, their errors and distance
    gaia.drop(gaia.columns[:5], axis=1, inplace=True) #keep columns 0-5

    # Save transformed data to a new file
    if save_spherical: #if we have entered save_sperical=True into the function
        gaia.to_csv(path.join(crossmatch_dir, f"{gaia_filename.split('.')[0]}_exoplanet_cross_spherical.csv"),
                    index=False) #save the gaia data to a csv as it is (in spherical coords). join.(...) is used to make the name of the csv
        exoplanets.to_csv(path.join(crossmatch_dir, f"{gaia_filename.split('.')[0]}_exoplanet_hosts.csv"), index=False) #save the exoplanets data to a csv as it is (spherical)

    if return_data: #return_data: Return crossmatched gaia data for further use
        return gaia #return the gaia df


def transform_to_cart(gaia, table_name, crossmatch_dir, setting="6d", predicted_radial_velocity=None):
    """

    :param gaia: Gaia dataset
    :param table_name: Name of the Gaia table used
    :param setting: Option used to convert either from 5D or 6D spherical to cartesian
    :param predicted_radial_velocity: Optional - True when using predicted radial velocity coordinate
    :return: Return Gaia dataset with converted coordinates
    """
    # The transformation you do depends on what information you have available. Whether radial velocity is or isn't available

    # First 3 coordinates remain the same for all options
    gaia["x"], gaia["y"], gaia["z"] = sph2cart(gaia["distance_pc"], gaia["ra"], gaia["dec"])

    if predicted_radial_velocity:
        # 5D coords with predicted radial velocity [hard]
        gaia["vx"], gaia["vy"], gaia["vz"] = vsph2cart(predicted_radial_velocity, gaia["pmra"], gaia["pmdec"],
                                                       gaia["distance_pc"], gaia["ra"], gaia["dec"])
        gaia.drop(gaia.columns[:-6], axis=1, inplace=True)
        setting = "predicted_rv"

    if setting == "6d":
        # Convert from spherical to cartesian coordinates [easy]
        gaia["vx"], gaia["vy"], gaia["vz"] = vsph2cart(gaia["dr2_radial_velocity"], gaia["pmra"], gaia["pmdec"],
                                                       gaia["distance_pc"], gaia["ra"], gaia["dec"])
        gaia.drop(gaia.columns[:-6], axis=1, inplace=True)
    elif setting == "5d_drop_rv":
        # 5D coords [average]
        gaia[["vx", "vy"]]= gaia[["pmra", "pmdec"]]
        gaia.drop(gaia.columns[:-5], axis=1, inplace=True)
    elif setting == "5d_drop_vz":
        gaia["vx"], gaia["vy"], gaia["vz"] = vsph2cart(gaia["dr2_radial_velocity"], gaia["pmra"], gaia["pmdec"],
                                                       gaia["distance_pc"], gaia["ra"], gaia["dec"])
        gaia.drop(gaia.columns[:-6], axis=1, inplace=True)
        gaia.drop(["vz"], axis=1, inplace=True)
    elif setting == "5d_drop_vy":
        # 5D coords [average]
        gaia["vx"], gaia["vy"], gaia["vz"] = vsph2cart(gaia["dr2_radial_velocity"], gaia["pmra"], gaia["pmdec"],
                                                       gaia["distance_pc"], gaia["ra"], gaia["dec"])
        gaia.drop(gaia.columns[:-6], axis=1, inplace=True)
        gaia.drop(["vy"], axis=1, inplace=True)
    elif setting == "5d_drop_vx":
        # 5D coords [average]
        gaia["vx"], gaia["vy"], gaia["vz"] = vsph2cart(gaia["dr2_radial_velocity"], gaia["pmra"], gaia["pmdec"],
                                                       gaia["distance_pc"], gaia["ra"], gaia["dec"])
        gaia.drop(gaia.columns[:-6], axis=1, inplace=True)
        gaia.drop(["vx"], axis=1, inplace=True)

    # Save to CSV
    gaia.to_csv(path.join(crossmatch_dir, f"{table_name}_exoplanet_cross_cartesian_{setting}.csv"), index=False)

    return gaia


# Function to convert positions from spherical to cartesian coordinates
def sph2cart(r, ra, dec):
    x = r * np.cos(ra) * np.cos(dec)
    y = r * np.sin(ra) * np.cos(dec)
    z = r * np.sin(dec)

    return x, y, z


# Function to convert positions from cartesian to spherical coordinates
def cart2sph(x, y, z):
    r = np.sqrt(x * x + y * y + z * z)
    dec = np.arcsin(z / r)
    ra = np.arctan2(y, x)

    return r, ra, dec


# Function to convert velocities from cartesian to spherical coordinates
def vcart2vsph(vx, vy, vz, x, y, z):
    r = np.sqrt(x * x + y * y * z * z)
    R = np.sqrt(x * x + y * y)
    rdot = x * vx + y * vy + z * vz #radial velocity
    rdot /= r 
    radot = vx * y - vy * x #velocity in ra
    radot /= R * R
    decdot = z * (x * vx + y * vy) - R * R * vz #velocity in dec
    decdot /= (R * r * r)

    return rdot, radot, decdot


# Function to convert velocities from spherical to cartesian coordinates
def vsph2cart(rdot, radot, decdot, r, ra, dec):
    xdot = np.cos(ra) * np.cos(dec) * rdot - r * np.sin(ra) * np.cos(dec) * radot - r * np.cos(ra) * np.sin(dec) * decdot #x velocity
    ydot = np.sin(ra) * np.cos(dec) * rdot + r * np.cos(ra) * np.cos(dec) * radot - r * np.sin(ra) * np.sin(dec) * decdot #y velocity
    zdot = np.sin(dec) * rdot + r * np.cos(dec) * decdot #z velocity

    return xdot, ydot, zdot

from astroquery.gaia import Gaia


class GaiaDataset:
    # a class is a set of functions and variables that you will use many times
    """
    Use this class to fetch GAIA dataset from ESA.
    """

    def __init__(self, query=None, table="gaiaedr3.gaia_source", filename="gaiaedr3"):
        # the init function, initialisation step. 
        # Cretes the self, query and table variables (defaults to thsoe inside the brackets)
        self.query = query # creating a variable called query
        self.filename = filename # creting a variable called filename
        self.table = table # etc

    def __login(self):
        """
        Initiate UI frame for user to log in to ESA servers.
        """

        Gaia.login_gui()

    def __default_query(self):
        """
        Default query used only if no query supplied by the user.
        """

        self.query = f"""
        select source_id, ra, ra_error, dec, dec_error, parallax, parallax_error, pmra, pmra_error, pmdec, pmdec_error,
        dr2_radial_velocity, dr2_radial_velocity_error from {self.table} where source_id is not null and ra is 
        not null and ra_error is not null and dec is not null and dec_error is not null and parallax is not null and 
        parallax_error is not null and pmra is not null and pmra_error is not null and pmdec is not null and 
        pmdec_error is not null and dr2_radial_velocity is not null and dr2_radial_velocity_error is not null
        """

    def get_gaia(self, query=None):
        """
        Fetch data from GAIA table
        :param query: SQL query supplied by user
        """

        self.__login() # runs the self._login() function above
        if query is None:
            self.__default_query()
            # use the default

        if self.filename[-4:] != ".csv":
            self.filename = self.filename + ".csv"

        Gaia.launch_job_async(self.query).get_results().to_pandas().to_csv("data/initial_datasets/" + self.filename,
                                                                           index=False)

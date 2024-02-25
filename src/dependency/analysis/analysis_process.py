class AnalysisProcess:

    def __init__(self):
        
        # ================================================================================
        #
        # Update the attributes
        #
        # ================================================================================
        
        # Attributes for the file reading
        self.list_flood_type       = ['inuncoast', 'inunriver']
        self.list_climate_scenario = ['historical', 'rcp4p5', 'rcp8p5']
        self.list_subsidence       = ['nosub', 'wtsub']
        self.list_model            = ['000000000WATCH', '00000NorESM1-M', '0000GFDL-ESM2M', '0000HadGEM2-ES', '00IPSL-CM5A-LR', 'MIROC-ESM-CHEM']
        self.list_year             = ['hist', '2030', '2050', '2080']
        self.list_returnperiod     = ['rp0002', 'rp0005', 'rp0010', 'rp0025', 'rp0050', 'rp0100', 'rp0250', 'rp0500', 'rp1000']
        self.list_projection       = ['0', '0_perc_05', '0_perc_50']

        # Attributes for the analysis
        self.flooding_depth_threshold = 0.1

        return